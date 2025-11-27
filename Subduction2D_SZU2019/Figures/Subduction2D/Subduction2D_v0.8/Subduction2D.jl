# Load script dependencies
using Pkg
Pkg.activate(".")

using GeoParams, GLMakie

const isCUDA = false

@static if isCUDA
    using CUDA
end
using JustRelax, JustRelax.JustRelax2D, JustRelax.DataIO

const backend = @static if isCUDA
    CUDABackend # Options: CPUBackend, CUDABackend, AMDGPUBackend
else
    JustRelax.CPUBackend # Options: CPUBackend, CUDABackend, AMDGPUBackend
end

using ParallelStencil, ParallelStencil.FiniteDifferences2D

@static if isCUDA
    @init_parallel_stencil(CUDA, Float64, 2)
else
    @init_parallel_stencil(Threads, Float64, 2)
end

using JustPIC, JustPIC._2D
# Threads is the default backend,
# to run on a CUDA GPU load CUDA.jl (i.e. "using CUDA") at the beginning of the script,
# and to run on an AMD GPU load AMDGPU.jl (i.e. "using AMDGPU") at the beginning of the script.
const backend_JP = @static if isCUDA
    CUDABackend # Options: CPUBackend, CUDABackend, AMDGPUBackend
else
    JustPIC.CPUBackend # Options: CPUBackend, CUDABackend, AMDGPUBackend
end

# Load file with all the rheology configurations
include("Subduction2D_setup.jl")
include("Subduction2D_rheology.jl")

## SET OF HELPER FUNCTIONS PARTICULAR FOR THIS SCRIPT --------------------------------

import ParallelStencil.INDICES
const idx_k = INDICES[2]
macro all_k(A)
    return esc(:($A[$idx_k]))
end

function copyinn_x!(A, B)
    @parallel function f_x(A, B)
        @all(A) = @inn_x(B)
        return nothing
    end

    return @parallel f_x(A, B)
end

# Initial pressure profile - not accurate
@parallel function init_P!(P, ρg, z)
    @all(P) = abs(@all(ρg) * @all_k(z)) * <(@all_k(z), 0.0)
    return nothing
end
## END OF HELPER FUNCTION ------------------------------------------------------------

## BEGIN OF MAIN SCRIPT --------------------------------------------------------------
function main(li, origin, phases_GMG, igg; nx=16, ny=16, figdir="figs2D", do_vtk=false)
    # Physical domain ------------------------------------
    ni = nx, ny           # number of cells
    di = @. li / ni       # grid steps
    grid = Geometry(ni, li; origin=origin) # JustRelax.jl/src/topology/Topology.jl
    (; xci, xvi) = grid # nodes at the center and vertices of the cells
    # ----------------------------------------------------

    # Physical properties using GeoParams ----------------
    rheology = init_rheology_nonNewtonian() # JustRelax.jl/miniapps/subduction/2D/Subduction2D_rheology.jl
    dt = 10.0e3 * 3600 * 24 * 365 # diffusive CFL timestep limiter
    # ----------------------------------------------------

    # Initialize particles -------------------------------
    nxcell = 40
    max_xcell = 60
    min_xcell = 20
    particles = init_particles( # From JustPIC.jl/src/Particles/particles_utils.jl. Makes Particles struct: src/particles.jl. Don't know why its not linked.
        backend_JP, nxcell, max_xcell, min_xcell, xvi, di, ni
    )
    subgrid_arrays = SubgridDiffusionCellArrays(particles) # JustRelax.jl/src/particles/subgrid_diffusion.jl
    grid_vxi = velocity_grids(xci, xvi, di) # JustRelax.jl/src/topology/Topology.jl
    # material phase & temperature

    pPhases, pT = init_cell_arrays(particles, Val(2)) # JustPIC.jl/src/Particles/particles_utils.jl

    # particle fields for the stress rotation. 
    pτ = StressParticles(particles) # Struct, src/stress_rotation/types.jl
    particle_args = (pT, pPhases, unwrap(pτ)...)
    particle_args_reduced = (pT, unwrap(pτ)...)

    # Assign particles phases anomaly
    phases_device = PTArray(backend)(phases_GMG) # JustRelax.jl/src/JustRelax.jl
    phase_ratios = phase_ratios = PhaseRatios(backend_JP, length(rheology), ni) # JustPIC.jl/src/PhaseRatios/PhaseRatios.jl
    init_phases!(pPhases, phases_device, particles, xvi) # JustRelax.jl/miniapps/subduction/2D/Subduction2D_rheology.jl
    update_phase_ratios!(phase_ratios, particles, xci, xvi, pPhases) # JustPIC.jl/src/PhaseRatios/utils.jl
    # ----------------------------------------------------

    # STOKES ---------------------------------------------
    # Allocate arrays needed for every Stokes problem
    stokes = StokesArrays(backend, ni) # JustRelax.jl/src/types/stokes.jl
    pt_stokes = PTStokesCoeffs(li, di; ϵ_abs=1.0e-4, ϵ_rel=1.0e-4, Re=1.0e0, r=0.7, CFL=0.9 / √2.1) # Re=3π, r=0.7 # JustRelax.jl/src/types/stokes.jl
    # ----------------------------------------------------

    # TEMPERATURE PROFILE --------------------------------
    Ttop = 20 + 273
    Tbot = maximum(T_GMG)
    thermal = ThermalArrays(backend, ni) # JustRelax.jl/src/types/heat_diffusion.jl
    @views thermal.T[2:(end-1), :] .= PTArray(backend)(T_GMG) # JustRelax.jl/src/JustRelax.jl
    thermal_bc = TemperatureBoundaryConditions(; # JustRelax.jl/src/boundaryconditions/types.jl
        no_flux=(left=true, right=true, top=false, bot=false),
    )
    thermal_bcs!(thermal, thermal_bc) # JustRelax.jl/src/boundaryconditions/BoundaryConditions.jl
    @views thermal.T[:, end] .= Ttop
    @views thermal.T[:, 1] .= Tbot
    temperature2center!(thermal) # JustRelax.jl/src/Interpolations.jl
    # ----------------------------------------------------

    # Buoyancy forces
    ρg = ntuple(_ -> @zeros(ni...), Val(2))
    compute_ρg!(ρg[2], phase_ratios, rheology, (T=thermal.Tc, P=stokes.P)) # JustRelax.jl/src/rheology/BuoyancyForces.jl
    stokes.P .= PTArray(backend)(reverse(cumsum(reverse((ρg[2]) .* di[2], dims=2), dims=2), dims=2)) # JustRelax.jl/src/JustRelax.jl

    # Rheology
    args0 = (T=thermal.Tc, P=stokes.P, dt=Inf)
    viscosity_cutoff = (1.0e18, 1.0e23)
    compute_viscosity!(stokes, phase_ratios, args0, rheology, viscosity_cutoff) # JustRelax.jl/src/rheology/Viscosity.jl

    # PT coefficients for thermal diffusion
    pt_thermal = PTThermalCoeffs( # JustRelax.jl/src/thermal_diffusion/DiffusionPT_coefficients.jl
        backend, rheology, phase_ratios, args0, dt, ni, di, li; ϵ=1.0e-8, CFL=0.95 / √2
    )

    # Boundary conditions
    flow_bcs = VelocityBoundaryConditions(; # JustRelax.jl/src/boundaryconditions/types.jl
        free_slip=(left=true, right=true, top=true, bot=true),
        free_surface=false,
    )
    flow_bcs!(stokes, flow_bcs) # apply boundary conditions, JustRelax.jl/src/boundaryconditions/BoundaryConditions.jl
    update_halo!(@velocity(stokes)...) # Update the halo of the given GPU/CPU-array(s). ~/.julia/packages/ImplicitGlobalGrid/dgclG/src/update_halo.jl

    # IO -------------------------------------------------
    # if it does not exist, make folder where figures are stored
    if do_vtk
        vtk_dir = joinpath(figdir, "vtk")
        take(vtk_dir)
    end
    take(figdir)
    # ----------------------------------------------------

    local Vx_v, Vy_v
    if do_vtk
        Vx_v = @zeros(ni .+ 1...)
        Vy_v = @zeros(ni .+ 1...)
    end

    T_buffer = @zeros(ni .+ 1)
    Told_buffer = similar(T_buffer)
    dt₀ = similar(stokes.P)
    for (dst, src) in zip((T_buffer, Told_buffer), (thermal.T, thermal.Told))
        copyinn_x!(dst, src)
    end
    grid2particle!(pT, xvi, T_buffer, particles) # JustPIC.jl/ext/JustPICCUDAExt.jl

    τxx_v = @zeros(ni .+ 1...)
    τyy_v = @zeros(ni .+ 1...)

    # Time loop
    t, it = 0.0, 0

    while it < 500 # run only for 5 Myr

        # interpolate fields from particle to grid vertices
        particle2grid!(T_buffer, pT, xvi, particles) # JustPIC.jl/ext/JustPICAMDGPUExt.jl
        @views T_buffer[:, end] .= Ttop
        @views T_buffer[:, 1] .= Tbot
        @views thermal.T[2:(end-1), :] .= T_buffer
        thermal_bcs!(thermal, thermal_bc) # JustRelax.jl/src/boundaryconditions/BoundaryConditions.jl
        temperature2center!(thermal) # JustRelax.jl/src/Interpolations.jl

        # interpolate stress back to the grid
        stress2grid!(stokes, pτ, xvi, xci, particles) # JustRelax.jl/src/stress_rotation/stress_rotation_particles.jl

        # Stokes solver ----------------
        args = (; T=thermal.Tc, P=stokes.P, dt=Inf)
        t_stokes = @elapsed begin
            out = solve!( # JustRelax.jl/src/stokes/Stokes2D.jl
                stokes,
                pt_stokes,
                di,
                flow_bcs,
                ρg,
                phase_ratios,
                rheology,
                args,
                dt,
                igg;
                kwargs=(
                    iterMax=100.0e3,
                    nout=2.0e3,
                    viscosity_cutoff=viscosity_cutoff,
                    free_surface=false,
                    viscosity_relaxation=1.0e-2,
                )
            )
        end

        # print some stuff
        println("Stokes solver time             ")
        println("   Total time:      $t_stokes s")
        println("   Time/iteration:  $(t_stokes / out.iter) s")

        # rotate stresses
        rotate_stress!(pτ, stokes, particles, xci, xvi, dt) # JustRelax.jl/src/stress_rotation/stress_rotation_particles.jl
        # compute time step
        dt = compute_dt(stokes, di) * 0.8 # JustRelax.jl/src/Utils.jl
        # compute strain rate 2nd invariant - for plotting
        tensor_invariant!(stokes.ε) # JustRelax.jl/src/stokes/StressKernels.jl
        # ------------------------------

        # Thermal solver ---------------
        heatdiffusion_PT!( # JustRelax.jl/src/thermal_diffusion/DiffusionPT_solver.jl
            thermal,
            pt_thermal,
            thermal_bc,
            rheology,
            args,
            dt,
            di;
            kwargs=(
                igg=igg,
                phase=phase_ratios,
                iterMax=50.0e3,
                nout=1.0e2,
                verbose=true,
            )
        )
        subgrid_characteristic_time!( # JustRelax.jl/src/particles/subgrid_diffusion.jl
            subgrid_arrays, particles, dt₀, phase_ratios, rheology, thermal, stokes, xci, di
        )
        centroid2particle!(subgrid_arrays.dt₀, xci, dt₀, particles) # JustPIC.jl/src/Interpolations/centroid_to_particle.jl
        subgrid_diffusion!( # JustPIC.jl/src/Physics/subgrid_diffusion.jl !!!
            pT, thermal.T, thermal.ΔT, subgrid_arrays, particles, xvi, di, dt
        )
        # ------------------------------

        # Advection --------------------
        # advect particles in space
        advection!(particles, RungeKutta2(), @velocity(stokes), grid_vxi, dt) # JustPIC.jl/src/Particles/Advection/advection.jl
        # advect particles in memory
        move_particles!(particles, xvi, particle_args) # JustPIC.jl/src/Particles/move.jl
        # check if we need to inject particles
        # need stresses on the vertices for injection purposes
        center2vertex!(τxx_v, stokes.τ.xx) # JustRelax.jl/src/Interpolations.jl
        center2vertex!(τyy_v, stokes.τ.yy)
        inject_particles_phase!( # JustPIC.jl/src/Particles/injection.jl
            particles,
            pPhases,
            particle_args_reduced,
            (T_buffer, τxx_v, τyy_v, stokes.τ.xy, stokes.ω.xy),
            xvi
        )

        # update phase ratios
        update_phase_ratios!(phase_ratios, particles, xci, xvi, pPhases) # JustPIC.jl/src/PhaseRatios/utils.jl

        @show it += 1
        t += dt

        # Data I/O and plotting ---------------------
        if it == 1 || rem(it, 10) == 0
            # checkpointing(figdir, stokes, thermal.T, η, t)
            (; η_vep, η) = stokes.viscosity
            if do_vtk
                velocity2vertex!(Vx_v, Vy_v, @velocity(stokes)...) # JustRelax.jl/src/Interpolations.jl
                data_v = (;
                    T=Array(T_buffer),
                    τII=Array(stokes.τ.II),
                    εII=Array(stokes.ε.II),
                    Vx=Array(Vx_v),
                    Vy=Array(Vy_v),
                )
                data_c = (;
                    P=Array(stokes.P),
                    η=Array(η_vep),
                )
                velocity_v = (
                    Array(Vx_v),
                    Array(Vy_v),
                )
                save_vtk(
                    joinpath(vtk_dir, "vtk_" * lpad("$it", 6, "0")),
                    xvi,
                    xci,
                    data_v,
                    data_c,
                    velocity_v;
                    t=t
                )
            end

            # Make particles plottable
            p = particles.coords
            ppx, ppy = p
            pxv = ppx.data[:] ./ 1.0e3
            pyv = ppy.data[:] ./ 1.0e3
            clr = pPhases.data[:]
            # clr      = pT.data[:]
            idxv = particles.index.data[:]

            # Make Makie figure
            ar = 3
            fig = Figure(size=(1200, 900), title="t = $t")
            ax1 = Axis(fig[1, 1], aspect=ar, title="T [K]  (t=$(t / (1.0e6 * 3600 * 24 * 365.25)) Myrs)")
            ax2 = Axis(fig[2, 1], aspect=ar, title="Phase")
            ax3 = Axis(fig[1, 3], aspect=ar, title="log10(εII)")
            ax4 = Axis(fig[2, 3], aspect=ar, title="log10(η)")
            # Plot temperature
            h1 = heatmap!(ax1, xvi[1] .* 1.0e-3, xvi[2] .* 1.0e-3, Array(thermal.T[2:(end-1), :]), colormap=:batlow)
            # Plot particles phase
            h2 = scatter!(ax2, Array(pxv[idxv]), Array(pyv[idxv]), color=Array(clr[idxv]), markersize=1)
            # Plot 2nd invariant of strain rate
            # h3  = heatmap!(ax3, xci[1].*1e-3, xci[2].*1e-3, Array(log10.(stokes.ε.II)) , colormap=:batlow)
            h3 = heatmap!(ax3, xci[1] .* 1.0e-3, xci[2] .* 1.0e-3, Array((stokes.τ.II)), colormap=:batlow)
            # Plot effective viscosity
            h4 = heatmap!(ax4, xci[1] .* 1.0e-3, xci[2] .* 1.0e-3, Array(log10.(stokes.viscosity.η_vep)), colormap=:batlow)
            hidexdecorations!(ax1)
            hidexdecorations!(ax2)
            hidexdecorations!(ax3)
            Colorbar(fig[1, 2], h1)
            Colorbar(fig[2, 2], h2)
            Colorbar(fig[1, 4], h3)
            Colorbar(fig[2, 4], h4)
            linkaxes!(ax1, ax2, ax3, ax4)
            fig
            save(joinpath(figdir, "$(it).png"), fig)

        end
        # ------------------------------

    end

    return nothing
end

# ## END OF MAIN SCRIPT ----------------------------------------------------------------
do_vtk = true # set to true to generate VTK files for ParaView
figdir = "Figures/Subduction2D/Subduction2D_v0.8"
n = 32
nx, ny = n * 2, n
li, origin, phases_GMG, T_GMG = GMG_subduction_2D(nx + 1, ny + 1)
igg = if !(JustRelax.MPI.Initialized()) # initialize (or not) MPI grid
    IGG(init_global_grid(nx, ny, 1; init_MPI=true)...)
else
    igg
end

# List of files you want to copy
script_files = [
    "Subduction2D.jl",
    "Subduction2D_rheology.jl",
    "Subduction2D_setup.jl",
]

# Ensure the figdir directory exists
isdir(figdir) || mkpath(figdir)

# Get the directory of the currently-running script
basepath = @__DIR__     # e.g. miniapps/subduction/2D/

# Copy each script into the figdir folder
for f in script_files
    src = joinpath(basepath, f)
    dest = joinpath(figdir, f)
    cp(src, dest; force=true)
end


main(li, origin, phases_GMG, igg; figdir=figdir, nx=nx, ny=ny, do_vtk=do_vtk);

# BvA, November 2025
# This script is based on JustRelax.jl/miniapps/subduction/2D.
# Here, the model is bundled with the _setup and _rheology file and
# is being tuned to the van Dinther (2019) paper on Mega Thrusts.


let
    try
        import Pkg
        pkgdir = joinpath(@__DIR__, "..") # access the Project.toml one level above
        Pkg.activate(pkgdir)
        Pkg.instantiate()
        @info "Successfully activated and instantiated project in $pkgdir"
    catch e
        @warn "Failed to activate/instantiate project: $e"
    end
end
using JustRelax

using GeoParams, CairoMakie

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
include("Subduction2D_setup_SZU2019.jl")
include("Subduction2D_rheology_SZU2019.jl")

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

function make_figure(
        it, t, dt,
        xvi, xci,
        T_buffer, ρg,
        stokes,
        Vx_limited, Vy_limited,
        pxv, pyv, clr, idxv,
        xmin, xmax, zmin, zmax,
        figpath
    )

    # Add isotherms
    isotherms_C = [100, 150, 350, 450, 900, 1300]
    isotherms_K = isotherms_C .+ 273
    
    # Match the T_buffer slicing
    xT = xvi[1][2:end-1] .* 1e-3
    zT = xvi[2] .* 1e-3
    TT = Array(T_buffer[2:end-1, :])

    # Velocity plot limits at 10 cm/yr (unit is in m/s)
    vmin = -0.1 / (365.25 * 24 * 3600)
    vmax = +0.1 / (365.25 * 24 * 3600)

    # Custom color range for viscosity
    visc_min = log10(1.0e18)
    visc_max = log10(1.0e24)

    # Custom color range for strain rate
    ε_min = -16
    ε_max = -12

    # Convert time to Myr and dt to Kyr for the title
    t_Myr = t / (3600 * 24 * 365.25 * 1.0e6)
    t_Myr_round = round(t_Myr; digits = 2)
    t_Kyr = round(dt / (3600 * 24 * 365.25 * 1000); digits = 2)

    fig = Figure(size = (1500, 600))
    Label(
        fig[0, :],
        "Subduction2D - SZU2019 | t = $t_Myr_round Myr | dt = $t_Kyr Kyr",
        fontsize = 28,
        tellwidth = false
    )

    ar = 2.8
    ax1 = Axis(fig[1, 1], aspect = ar, title = "Material Phase")
    ax2 = Axis(fig[2, 1], aspect = ar, title = "Temperature [K]")
    ax3 = Axis(fig[3, 1], aspect = ar, title = "Density [kg/m³]")
    ax4 = Axis(fig[1, 3], aspect = ar, title = "Vx [m/s]")
    ax5 = Axis(fig[2, 3], aspect = ar, title = "Vy [m/s]")
    ax6 = Axis(fig[3, 3], aspect = ar, title = "log10(τII) [Pa]")
    ax7 = Axis(fig[1, 5], aspect = ar, title = "log10(εII)")
    ax8 = Axis(fig[2, 5], aspect = ar, title = "log10(η)")
    ax9 = Axis(fig[3, 5], aspect = ar, title = "log10(η_vep)")

    # Apply zoom limits
    for ax in (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9)
        xlims!(ax, xmin, xmax)
        ylims!(ax, zmin, zmax)
    end

    # ---- Plots ----
    # Material phase
    h1 = scatter!(ax1, pxv[idxv], pyv[idxv], color = clr[idxv], markersize = 1)

    contour!(
        ax1,
        xT,
        zT,
        TT;
        levels = isotherms_K,
        labels=true,
        labelsize = 12,
        color = :white,
        linewidth = 1.5
    )

    # Temperature
    h2 = heatmap!(
        ax2,
        xvi[1] .* 1.0e-3,
        xvi[2] .* 1.0e-3,
        Array(T_buffer[2:(end - 1), :]),
        colormap = :vik
    )

    contour!(
        ax2,
        xT,
        zT,
        TT;
        levels = isotherms_K,
        labels=true,
        labelsize = 12,
        color = :white,
        linewidth = 1.5
    )


    # Density
    h3 = heatmap!(
        ax3,
        xci[1] .* 1.0e-3,
        xci[2] .* 1.0e-3,
        Array(ustrip.(ρg[2] ./ 9.81));
        colormap = :vik,
        colorrange = (2500, 3300)
    )

    contour!(
        ax3,
        xT,
        zT,
        TT;
        levels = isotherms_K,
        labels=true,
        labelsize = 12,
        color = :white,
        linewidth = 1.5
    )

    # Velocity x-direction
    h4 = heatmap!(
        ax4,
        xvi[1] .* 1.0e-3,
        xvi[2] .* 1.0e-3,
        Vx_limited;
        colormap = :vik,
        colorrange = (vmin, vmax)
    )

    # Velocity y-direction
    h5 = heatmap!(
        ax5,
        xvi[1] .* 1.0e-3,
        xvi[2] .* 1.0e-3,
        Vy_limited;
        colormap = :vik,
        colorrange = (vmin, vmax)
    )

    # Second invariant of the stress tensor
    h6 = heatmap!(
        ax6,
        xci[1] .* 1.0e-3,
        xci[2] .* 1.0e-3,
        Array(log10.(stokes.τ.II)),
        colormap = :batlow
    )

    # Second invariant of the strain rate tensor
    h7 = heatmap!(
        ax7,
        xci[1] .* 1.0e-3,
        xci[2] .* 1.0e-3,
        Array(log10.(stokes.ε.II));
        colormap = :vik,
        colorrange = (ε_min, ε_max)
    )

    contour!(
        ax7,
        xT,
        zT,
        TT;
        levels = isotherms_K,
        labels=true,
        labelsize = 12,
        color = :white,
        linewidth = 1.5
    )

    # Viscosity
    h8 = heatmap!(
        ax8,
        xci[1] .* 1.0e-3,
        xci[2] .* 1.0e-3,
        Array(log10.(stokes.viscosity.η));
        colormap = :vik,
        colorrange = (visc_min, visc_max)
    )
    contour!(
        ax8,
        xT,
        zT,
        TT;
        levels = isotherms_K,
        labels=true,
        labelsize = 12,
        color = :white,
        linewidth = 1.5
    )

    # Visco_plasto_elastic viscosity
    h9 = heatmap!(
        ax9,
        xci[1] .* 1.0e-3,
        xci[2] .* 1.0e-3,
        Array(log10.(stokes.viscosity.η_vep));
        colormap = :vik,
        colorrange = (visc_min, visc_max)
    )
    contour!(
        ax9,
        xT,
        zT,
        TT;
        levels = isotherms_K,
        labels=true,
        labelsize = 12,
        color = :white,
        linewidth = 1.5
    )


    # Colorbars
    for (row, col, h) in [
            (1, 2, h1),
            (2, 2, h2),
            (3, 2, h3),
            (1, 4, h4),
            (2, 4, h5),
            (3, 4, h6),
            (1, 6, h7),
            (2, 6, h8),
            (3, 6, h9),
        ]
        gl = fig[row, col] = GridLayout()
        Colorbar(gl[1, 1], h; height = 120)
    end

    linkaxes!(ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9)

    save(figpath, fig)
    return nothing
end
## END OF HELPER FUNCTION ------------------------------------------------------------


## BEGIN OF MAIN SCRIPT --------------------------------------------------------------
function main(li, origin, phases_GMG, igg; nx = 16, ny = 16, figdir = "figs2D", do_vtk = false)

    # Physical domain ------------------------------------
    ni = nx, ny           # number of cells
    di = @. li / ni       # grid steps
    grid = Geometry(ni, li; origin = origin)
    (; xci, xvi) = grid # nodes at the center and vertices of the cells
    # ----------------------------------------------------

    # Physical properties using GeoParams ----------------
    rheology = init_rheologies()
    dt = 10.0e3 * 3600 * 24 * 365 # diffusive CFL timestep limiter, seconds
    # ----------------------------------------------------

    # Initialize particles -------------------------------
    nxcell = 40
    max_xcell = 60
    min_xcell = 20
    particles = init_particles(
        backend_JP, nxcell, max_xcell, min_xcell, xvi...
    )
    subgrid_arrays = SubgridDiffusionCellArrays(particles)
    # velocity grids
    grid_vxi = velocity_grids(xci, xvi, di)
    # material phase & temperature
    pPhases, pT = init_cell_arrays(particles, Val(2))

    # particle fields for the stress rotation.
    pτ = StressParticles(particles)
    particle_args = (pT, pPhases, unwrap(pτ)...)
    particle_args_reduced = (pT, unwrap(pτ)...)

    # Assign particles phases anomaly
    phases_device = PTArray(backend)(phases_GMG)
    phase_ratios = phase_ratios = PhaseRatios(backend_JP, length(rheology), ni)
    init_phases!(pPhases, phases_device, particles, xvi)
    update_phase_ratios!(phase_ratios, particles, xci, xvi, pPhases)
    # ----------------------------------------------------

    # STOKES ---------------------------------------------
    # Allocate arrays needed for every Stokes problem
    stokes = StokesArrays(backend, ni)
    pt_stokes = PTStokesCoeffs(li, di; ϵ_abs = 1.0e-4, ϵ_rel = 1.0e-4, Re = 1.0e0, r = 0.7, CFL = 0.9 / √2.1) # Re=3π, r=0.7
    # ----------------------------------------------------

    # TEMPERATURE PROFILE --------------------------------
    Ttop = 20 + 273
    Tbot = maximum(T_GMG)
    thermal = ThermalArrays(backend, ni)
    @views thermal.T[2:(end - 1), :] .= PTArray(backend)(T_GMG)
    thermal_bc = TemperatureBoundaryConditions(;
        no_flux = (left = true, right = true, top = false, bot = false),
    )
    thermal_bcs!(thermal, thermal_bc)
    @views thermal.T[:, end] .= Ttop
    @views thermal.T[:, 1] .= Tbot
    temperature2center!(thermal)
    # ----------------------------------------------------

    # Buoyancy forces
    ρg = ntuple(_ -> @zeros(ni...), Val(2))
    compute_ρg!(ρg[2], phase_ratios, rheology, (T = thermal.Tc, P = stokes.P))
    stokes.P .= PTArray(backend)(reverse(cumsum(reverse((ρg[2]) .* di[2], dims = 2), dims = 2), dims = 2))

    # Rheology
    args0 = (T = thermal.Tc, P = stokes.P, dt = Inf)
    viscosity_cutoff = (1.0e18, 5.0e23)
    compute_viscosity!(stokes, phase_ratios, args0, rheology, viscosity_cutoff)

    # PT coefficients for thermal diffusion
    pt_thermal = PTThermalCoeffs(
        backend, rheology, phase_ratios, args0, dt, ni, di, li; ϵ = 1.0e-8, CFL = 0.95 / √2
    )
    ###### Work in progress, not ready for GPU yet. ######
    ###### Requires JustRelax custom patch
    # Bert added: Special box of nodes where eastward slip of subducting plate is fixed.
    # Work in progress, since v0.10 and working on CPU in v0.86.
    # Will be replaced with a GPU-compatible version.

    # Vx array and sizes
    Vx = stokes.V.Vx
    nxV, nyV = size(Vx)

    # Coordinates of Vx nodes
    xV = grid_vxi[1][1]   # length nxV
    zV = grid_vxi[1][2]   # length nyV

    # Dirichlet mask
    mask = falses(nxV, nyV)

    for j in 1:nyV
        for i in 1:nxV
            if 170e3 < xV[i] < 190e3 &&
            -68e3 < zV[j] < -24.5e3
                mask[i, j] = true
            end
        end
    end

    # Boundary conditions
    flow_bcs = VelocityBoundaryConditions(;
        free_slip    = (left = true, right = true, top = true, bot = true),
        free_surface = false,
        dirichlet    = (; constant = 7.5e-2 / (3600 * 24 * 365.25), mask = mask),
    )

    flow_bcs!(stokes, flow_bcs) # apply boundary conditions, custom edit
    update_halo!(@velocity(stokes)...)
    ######### End WIP
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
    grid2particle!(pT, xvi, T_buffer, particles)

    τxx_v = @zeros(ni .+ 1...)
    τyy_v = @zeros(ni .+ 1...)

    # Time loop
    t, it = 0.0, 0

    while it < 500 # run only for 5 Myr

        # interpolate fields from particle to grid vertices
        particle2grid!(T_buffer, pT, xvi, particles)
        @views T_buffer[:, end] .= Ttop
        @views T_buffer[:, 1] .= Tbot
        @views thermal.T[2:(end - 1), :] .= T_buffer
        thermal_bcs!(thermal, thermal_bc)
        temperature2center!(thermal)

        # interpolate stress back to the grid
        stress2grid!(stokes, pτ, xvi, xci, particles)

        # Stokes solver ----------------
        args = (; T = thermal.Tc, P = stokes.P, dt = Inf)
        t_stokes = @elapsed begin
            out = solve!(
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
                kwargs = (
                    iterMax = 100.0e3,
                    nout = 2.0e3,
                    viscosity_cutoff = viscosity_cutoff,
                    free_surface = false,
                    viscosity_relaxation = 1.0e-2,
                )
            )
        end

        # print some stuff
        println("Stokes solver time             ")
        println("   Total time:      $t_stokes s")
        println("   Time/iteration:  $(t_stokes / out.iter) s")

        # rotate stresses
        rotate_stress!(pτ, stokes, particles, xci, xvi, dt)
        # compute time step
        dt = compute_dt(stokes, di) * 0.8
        # compute strain rate 2nd invariant - for plotting
        tensor_invariant!(stokes.ε)
        # ------------------------------

        # Thermal solver ---------------
        heatdiffusion_PT!(
            thermal,
            pt_thermal,
            thermal_bc,
            rheology,
            args,
            dt,
            di;
            kwargs = (
                igg = igg,
                phase = phase_ratios,
                iterMax = 50.0e3,
                nout = 1.0e2,
                verbose = true,
            )
        )
        subgrid_characteristic_time!(
            subgrid_arrays, particles, dt₀, phase_ratios, rheology, thermal, stokes, xci, di
        )
        centroid2particle!(subgrid_arrays.dt₀, xci, dt₀, particles)
        subgrid_diffusion!(
            pT, thermal.T, thermal.ΔT, subgrid_arrays, particles, xvi, di, dt
        )
        # ------------------------------

        # Advection --------------------
        # advect particles in space
        advection_MQS!(particles, RungeKutta2(), @velocity(stokes), grid_vxi, dt)        # advect particles in memory
        # advect particles in memory
        move_particles!(particles, xvi, particle_args)
        # check if we need to inject particles
        # need stresses on the vertices for injection purposes
        center2vertex!(τxx_v, stokes.τ.xx)
        center2vertex!(τyy_v, stokes.τ.yy)
        inject_particles_phase!(
            particles,
            pPhases,
            particle_args_reduced,
            (T_buffer, τxx_v, τyy_v, stokes.τ.xy, stokes.ω.xy),
            xvi
        )

        # update phase ratios
        update_phase_ratios!(phase_ratios, particles, xci, xvi, pPhases)

        @show it += 1
        t += dt

        # # Data I/O and plotting ---------------------
        # if it == 1 #|| rem(it, 10) == 0
        # checkpointing(figdir, stokes, thermal.T, η, t)
        (; η_vep, η) = stokes.viscosity
        if do_vtk
            velocity2vertex!(Vx_v, Vy_v, @velocity(stokes)...)
            data_v = (;
                T = Array(T_buffer),
                τII = Array(stokes.τ.II),
                εII = Array(stokes.ε.II),
                Vx = Array(Vx_v),
                Vy = Array(Vy_v),
            )
            data_c = (;
                P = Array(stokes.P),
                η = Array(η_vep),
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
                t = t
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

        # --- New figure: velocity with limited range ---
        vmin = -0.1 / (365.25 * 24 * 3600)     # -10 cm/yr
        vmax = +0.1 / (365.25 * 24 * 3600)     # +10 cm/yr
        Vx_raw = Array(stokes.V.Vx[2:(end - 1), :])
        Vx_limited = clamp.(Vx_raw, vmin, vmax)
        Vy_raw = Array(stokes.V.Vy[2:(end - 1), :])
        Vy_limited = clamp.(Vy_raw, vmin, vmax)


        # --- ZOOM REGION ---
        xmin_zoom, xmax_zoom = 800, 1100
        zmin_zoom, zmax_zoom = -80, 0

        # --- FULL DOMAIN ---
        xmin_full = minimum(xvi[1]) * 1.0e-3
        xmax_full = maximum(xvi[1]) * 1.0e-3
        zmin_full = minimum(xvi[2]) * 1.0e-3
        zmax_full = maximum(xvi[2]) * 1.0e-3

        # ZOOMED FIGURE
        make_figure(
            it, t, dt,
            xvi, xci,
            T_buffer, ρg,
            stokes,
            Vx_limited, Vy_limited,
            pxv, pyv, clr, idxv,
            xmin_zoom, xmax_zoom, zmin_zoom, zmax_zoom,
            joinpath(figdir, "$(lpad(it, 2, "0"))_zoom.png")
        )

        # FULL DOMAIN FIGURE
        make_figure(
            it, t, dt,
            xvi, xci,
            T_buffer, ρg,
            stokes,
            Vx_limited, Vy_limited,
            pxv, pyv, clr, idxv,
            xmin_full, xmax_full, zmin_full, zmax_full,
            joinpath(figdir, "$(lpad(it, 2, "0"))_full.png")
        )
    end

    # ------------------------------

    return nothing
end


# ## END OF MAIN SCRIPT ----------------------------------------------------------------
do_vtk = true # set to true to generate VTK files for ParaView
version = "v0.135"
figdir = "Subduction2D_SZU2019/Figures/Subduction2D_SZU2019/SZU2019_$version"
println(version)
n = 32 
nx, ny = n*2, n
li, origin, phases_GMG, T_GMG = GMG_subduction_2D(nx + 1, ny + 1)
igg = if !(JustRelax.MPI.Initialized()) # initialize (or not) MPI grid
    IGG(init_global_grid(nx, ny, 1; init_MPI = true)...)
else
    igg
end

# List of files you want to copy
script_files = [
    basename(@__FILE__),
    "Subduction2D_setup_SZU2019.jl",
    "Subduction2D_rheology_SZU2019.jl",
]

# Ensure the figdir directory exists
isdir(figdir) || mkpath(figdir)

# Get the directory of the currently-running script
basepath = @__DIR__
prefix = "used_"

# Copy each script into the figdir folder
for f in script_files
    src = joinpath(basepath, f)
    name, ext = splitext(f)
    dest = joinpath(figdir, prefix * name * ext)
    cp(src, dest; force = true)
end


main(li, origin, phases_GMG, igg; figdir = figdir, nx = nx, ny = ny, do_vtk = do_vtk);

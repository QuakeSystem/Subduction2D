# Load script dependencies
using GeoParams, CairoMakie, LinearAlgebra, HDF5
const isCUDA = true

remote = true
if remote
    working_dir = "/scratch/tectonics/bert/Subduction2D"
else
    working_dir = "/Users/5723272/SD/Subduction2D"
end

@static if isCUDA
    using CUDA
end

using JustRelax, JustRelax.JustRelax2D, JustRelax.DataIO

const backend = @static if isCUDA
    CUDABackend # Options: CPUBackend, CUDABackend, AMDGPUBackend
    const backend_JR = CUDABackend
else
    JustRelax.CPUBackend # Options: CPUBackend, CUDABackend, AMDGPUBackend
    const backend_JR = CPUBackend
end

using ParallelStencil, ParallelStencil.FiniteDifferences2D

@static if isCUDA
    @init_parallel_stencil(CUDA, Float64, 2)
else
    @init_parallel_stencil(Threads, Float64, 2)
end

using JustPIC
# Threads is the default backend,
# to run on a CUDA GPU load CUDA.jl (i.e. "using CUDA") at the beginning of the script,
# and to run on an AMD GPU load AMDGPU.jl (i.e. "using AMDGPU") at the beginning of the script.
const backend_JP = @static if isCUDA
    CUDABackend # Options: CPUBackend, CUDABackend, AMDGPUBackend
else
    JustPIC.CPUBackend # Options: CPUBackend, CUDABackend, AMDGPUBackend
end

include(joinpath(working_dir, "LocalModules/visualisation.jl"))
include(joinpath(working_dir, "LocalModules/VelocityBoxes.jl"))
include(joinpath(working_dir, "LocalModules/HelperFunctions.jl"))
include(joinpath(working_dir, "LocalModules/NonuniformGrid.jl"))

# Load file with all the rheology configurations
setup_file = "Subduction2D_setup.jl"
rheology_file = "Subduction2D_rheology.jl"
include(setup_file)
include(rheology_file)

# Velocity box application kernels -- needs @init_parallel_stencil (above)
# to have already run, and needs VelBox2D (from VelocityBoxes.jl, above).
include(joinpath(working_dir,"LocalModules/VelocityBoxKernels.jl"))

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
function main(
    li,
    origin,
    phases_GMG,
    T_GMG,
    igg;
    xvi = nothing,
    xci = nothing,
    nx = 16,
    ny = 16,
    ref_grid = 0,
    version = nothing,
)

    # Physical domain ------------------------------------
    ni = nx, ny           # number of cells

    # non-uniform grid with refinement
    if ref_grid == 1
        grid = Geometry(
            PTArray(backend_JR),
            xvi...,
        )
        di_min =  (min(minimum.(grid.di.center)...),
        min(minimum.(grid.di.vertex)...))
    else
        grid = Geometry(ni, li; origin = origin)
        di_min = @. li / ni       # grid steps
    end
    
    (; xci, xvi) = grid # nodes at the center and vertices of the cells

    # ----------------------------------------------------
    # Set flags and parameters for visualization and output and create folders for output
    vis = prepare_visualisation(ni, version=version)

    # Physical properties using GeoParams ----------------
    rheology = init_rheologies_start()
    dt = 25.0e3 * 3600 * 24 * 365 # diffusive CFL timestep limiter
    dt_max = 25.0e3 * 3600 * 24 * 365 # diffusive CFL timestep limiter
    # ----------------------------------------------------

    # Initialize particles -------------------------------
    nxcell = 40
    max_xcell = 60
    min_xcell = 20
    # JustPIC's `init_particles` expects host coordinate vectors (its parameter is
    # named `xi_vel_cpu`) and uploads to the backend itself. With the refined grid,
    # `grid.xi_vel` already lives on the GPU, so hand it CPU copies to avoid scalar
    # indexing in `add_periodic_ghost_nodes`.
    particles = init_particles(
        backend_JP, nxcell, max_xcell, min_xcell,
        map(t -> Array.(t), grid.xi_vel)...
    )
    subgrid_arrays = SubgridDiffusionCellArrays(particles; loc = :center)
    # grid_vxi = velocity_grids(xci, xvi, grid.di.vertex)
    grid_vxi = velocity_grids_gpu(xci, xvi, grid.di.vertex)
    # material phase & temperature
    pPhases, pT = init_cell_arrays(particles, Val(2))

    # particle fields for the stress rotation
    pτ = StressParticles(particles)
    particle_args = (pT, pPhases, unwrap(pτ)...)
    particle_args_reduced = (pT, unwrap(pτ)...)

    # Assign particles phases anomaly
    phases_device = PTArray(backend)(phases_GMG)
    phase_ratios = phase_ratios = PhaseRatios(backend_JP, length(rheology), ni)
    init_phases!(pPhases, phases_device, particles, xvi)
    update_phase_ratios!(phase_ratios, particles, pPhases)
    # ----------------------------------------------------

    # STOKES ---------------------------------------------
    # Allocate arrays needed for every Stokes problem
    stokes = StokesArrays(backend, ni)
    # ----------------------------------------------------

    # TEMPERATURE PROFILE --------------------------------
    Ttop = 0 + 273
    Tbot = maximum(T_GMG)
    thermal = ThermalArrays(backend, ni)
    vertex2center!(thermal.T, PTArray(backend)(T_GMG); ghost_x = true, ghost_y = true)
    thermal_bc = TemperatureBoundaryConditions(;
        no_flux = (left = true, right = true, top = false, bot = false),
        constant_value = (left = false, right = false, top = Ttop, bot = Tbot),
    )
    thermal_bcs!(thermal, thermal_bc)
    # ----------------------------------------------------

    # Buoyancy forces
    ρg = ntuple(_ -> @zeros(ni...), Val(2))
    compute_ρg!(ρg[2], phase_ratios, rheology, (T = thermal.T, P = stokes.P))
    if ref_grid == 0
        stokes.P .= PTArray(backend)(reverse(cumsum(reverse((ρg[2]) .* di_min[2], dims = 2), dims = 2), dims = 2))
    else
        # Lithostatic pressure integrates vertical body force using local cell dy (vertex spacing).
        stokes.P .= PTArray(backend)(reverse(cumsum(reverse((ρg[2]) .* reshape(grid.di.vertex[2], 1, :), dims = 2), dims = 2), dims = 2))
    end

    # Rheology
    args0 = (T = thermal.T, P = stokes.P, dt = Inf)
    viscosity_cutoff = (1.0e18, 1.0e24)
    compute_viscosity!(stokes, phase_ratios, args0, rheology, viscosity_cutoff)
    center2vertex!(stokes.viscosity.ηv, stokes.viscosity.η)
    # ----------------------------------------------------

    # PT coefficients for thermal diffusion
    pt_thermal = PTThermalCoeffs(
        backend, rheology, phase_ratios, args0, dt, ni, di_min, li; ϵ = 1.0e-8, CFL = 0.95 / √2
    )

    # Boundary conditions (includes the internal velocity-box dirichlet region,
    # rebuilt each timestep below as the boxes ramp up -- see VelocityBoxKernels.jl)
    flow_bcs = velocity_box_flow_bcs(
        stokes, grid, vel_boxes_2D;
        free_slip = (left = true, right = true, top = true, bot = true),
        free_surface = false,
    )
    flow_bcs!(stokes, flow_bcs) # apply boundary conditions
    update_halo!(@velocity(stokes)...)

    # visualization prep moved to LocalModules/visualisation.jl
    T_buffer = thermal.T[2:(end - 1), 2:(end - 1)]
    dt₀ = similar(stokes.P)
    centroid2particle!(pT, T_buffer, particles)

    τxx_v = @zeros(ni .+ 1...)
    τyy_v = @zeros(ni .+ 1...)

    dyrel = DYREL(backend, stokes, rheology, phase_ratios, grid.di, dt; ϵ = 1.0e-3)

    # Time loop
    t, it = 0.0, 0
    while it <= 2000 || t < 4e6 * (3600 * 24 * 365.25)  # run only for 4 Myrs
        if it == 5
            vel_boxes_2D[1] = VelBox2D(vel_boxes_2D[1].cenx, vel_boxes_2D[1].cenz, vel_boxes_2D[1].widthx, vel_boxes_2D[1].widthz, 2.5 * 0.01 / (3600*24*365), vel_boxes_2D[1].vy, true, vel_boxes_2D[1].has_vy)
            rheology = init_rheologies()
        elseif it == 10
            vel_boxes_2D[1] = VelBox2D(vel_boxes_2D[1].cenx, vel_boxes_2D[1].cenz, vel_boxes_2D[1].widthx, vel_boxes_2D[1].widthz, 5.0 * 0.01 / (3600*24*365), vel_boxes_2D[1].vy, true, vel_boxes_2D[1].has_vy)
        elseif it == 15
            vel_boxes_2D[1] = VelBox2D(vel_boxes_2D[1].cenx, vel_boxes_2D[1].cenz, vel_boxes_2D[1].widthx, vel_boxes_2D[1].widthz, 7.5 * 0.01 / (3600*24*365), vel_boxes_2D[1].vy, true, vel_boxes_2D[1].has_vy)
        end
        # Get flat views of the raw data
        phases_flat = pPhases.data[:]   # all particle phase values
        temps_flat  = pT.data[:]        # all particle temperatures
        index_flat  = particles.index.data[:]  # true = active particle
        # Find active air particles
        air_mask = (phases_flat .== 2.0) .& index_flat
        @show sum(air_mask)
        @show extrema(temps_flat[air_mask])
        @show mean(temps_flat[air_mask])   # needs Statistics
        # Set air particle temperatures to 273 K
        pT.data[air_mask] .= 273.0

        # interpolate fields from particle to grid vertices
        particle2centroid!(T_buffer, pT, particles; ghost_1 = false, ghost_2 = false, ghost_3 = false)
        @views thermal.T[2:end-1, 2:end-1] .= T_buffer
        thermal_bcs!(thermal, thermal_bc)
        # interpolate stress back to the grid
        stress2grid!(stokes, pτ, particles)

        # Rebuild flow_bcs' internal velocity-box dirichlet region so it picks up any
        # change to vel_boxes_2D (e.g. the ramped-up prescribed velocity at it==5,10,15);
        # solve_DYREL! itself re-asserts the prescribed value every DR iteration.
        flow_bcs = velocity_box_flow_bcs(
            stokes, grid, vel_boxes_2D;
            free_slip = (left = true, right = true, top = true, bot = true),
            free_surface = false,
        )
        update_halo!(@velocity(stokes)...)

        # Stokes solver ----------------
        args = (; T = thermal.T, P = stokes.P, dt = Inf)
        t_stokes = @elapsed begin
            out = solve_DYREL!(
                stokes,
                ρg,
                dyrel,
                flow_bcs,
                phase_ratios,
                rheology,
                args,
                grid,
                dt,
                igg;
                kwargs = (;
                    verbose_PH = true,
                    verbose_DR = true,
                    iterMax = 50.0e2,
                    rel_drop = 1.0e-2,
                    nout = 400,
                    λ_relaxation_PH = 1,
                    λ_relaxation_DR = 1,
                    viscosity_relaxation = 1.0e-2,
                    viscosity_cutoff = viscosity_cutoff,
                )
            )
        end
        # print some stuff
        println("Stokes solver time             ")
        println("   Total time:      $t_stokes s")
        # println("   Time/iteration:  $(t_stokes / out.iter) s")

        # rotate stresses
        rotate_stress!(pτ, stokes, particles, dt)
        # compute time step
        dt_plot = dt
        dt = compute_dt(stokes, di_min, dt_max) #* 0.8
        # compute strain rate 2nd invartian - for plotting
        tensor_invariant!(stokes.τ)
        tensor_invariant!(stokes.ε)
        tensor_invariant!(stokes.ε_pl)
        # ------------------------------

        # Thermal solver ---------------
        heatdiffusion_PT!(
            thermal,
            pt_thermal,
            thermal_bc,
            rheology,
            args,
            dt,
            grid;
            kwargs = (
                igg = igg,
                phase = phase_ratios,
                iterMax = 50.0e3,
                nout = 1.0e2,
                verbose = true,
            )
        )
        subgrid_characteristic_time!(
            subgrid_arrays, particles, dt₀, phase_ratios, rheology, thermal, stokes
        )
        centroid2particle!(subgrid_arrays.dt₀, dt₀, particles)
        subgrid_diffusion_centroid!(
            pT, T_buffer, thermal.ΔT, subgrid_arrays, particles, dt
        )
        # ------------------------------

        # Advection --------------------
        # advect particles in space
        advection_MQS!(particles, RungeKutta2(), @velocity(stokes), dt)
        # advect particles in memory
        move_particles!(particles, particle_args)
        # check if we need to inject particles
        # need stresses on the vertices for injection purposes
        inject_particles_phase!(
            particles,
            pPhases,
            particle_args_reduced,
            (T_buffer, stokes.τ.xx_v, stokes.τ.yy_v, stokes.τ.xy, stokes.ω.xy)
        )

        # update phase ratios
        update_phase_ratios!(phase_ratios, particles, pPhases)
        ### PARAVIEW PLOTTING
        if it == 1 || rem(it, 25) == 0
            checkpointing_jld2(vis.checkpoint, stokes, thermal, t, dt; it = it)
            checkpointing_particles(vis.checkpoint, particles; phases = pPhases, phase_ratios = phase_ratios, particle_args = particle_args, particle_args_reduced = particle_args_reduced, t = t, dt = dt, it = it)
        end
        (; η_vep, η) = stokes.viscosity
        if vis.do_vtk && (it == 1 || rem(it, vis.vtk_every) == 0)
            velocity2vertex!(vis.Vx_v, vis.Vy_v, @velocity(stokes)...)
            # Reconstruct compact phase "shapes" on the grid from particle phase ratios.
            phase_vertex = [argmax(p) for p in Array(phase_ratios.vertex)]
            Rx_c = zeros(size(stokes.P))
            Ry_c = zeros(size(stokes.P))
            @views Rx_c[axes(stokes.R.Rx, 1), axes(stokes.R.Rx, 2)] .= Array(stokes.R.Rx)
            @views Ry_c[axes(stokes.R.Ry, 1), axes(stokes.R.Ry, 2)] .= Array(stokes.R.Ry)

            data_v = (;
                τII            = Array(stokes.τ.II),
                εII            = Array(stokes.ε.II),
                Vx             = Array(vis.Vx_v),
                Vy             = Array(vis.Vy_v),
                phase_vertex   = phase_vertex,
                ResT           = Array(thermal.ResT),
                log10_absResT  = log10.(abs.(Array(thermal.ResT)) .+ 1e-30),
                dτ_ρ           = Array(pt_thermal.dτ_ρ),
                log10_dτ_ρ     = log10.(abs.(Array(pt_thermal.dτ_ρ)) .+ 1e-30),
                θr_dτ          = Array(pt_thermal.θr_dτ),
            )
            data_c = (;
                T_buffer   = Array(T_buffer),
                thermal_T = Array(thermal.T[2:end-1,2:end-1]),
                P   = Array(stokes.P),
                P0   = Array(stokes.P0),
                density = Array(ustrip.(ρg[2]) ./ 9.81),
                divV   = Array(stokes.∇V),
                η_vep   = Array(η_vep),
                Rx  = Array(Rx_c),
                Ry  = Array(Ry_c),
                Rmag = sqrt.(Rx_c .^ 2 .+ Ry_c .^ 2),
                log10_η     = log10.(Array(stokes.viscosity.η)),
                log10_η_vep = log10.(Array(stokes.viscosity.η_vep)),
                λ         = Array(stokes.λ),
                EII_pl    = Array(stokes.EII_pl),
                EVol_pl   = Array(stokes.EVol_pl),
                ε_vol_pl  = Array(stokes.ε_vol_pl),
                ΔPψ = Array(stokes.ΔPψ),
                τxx = Array(stokes.τ.xx),
                τyy = Array(stokes.τ.yy),
                τxy = Array(stokes.τ.xy_c),
                εxx = Array(stokes.ε.xx),
                εyy = Array(stokes.ε.yy),
                εxy = Array(stokes.ε.xy_c),
                εII = Array(stokes.ε.II),
                ΔT            = Array(thermal.ΔT[2:end-1,2:end-1]),
                adiabatic     = Array(thermal.adiabatic),
                dT_dt         = Array(thermal.dT_dt),
                H             = Array(thermal.H),
                shear_heating = Array(thermal.shear_heating),
            )
            velocity_v = (
                Array(vis.Vx_v),
                Array(vis.Vy_v),
            )
            path_vtk = joinpath(vis.vtk_dir, "vtk_" * lpad("$it", 6, "0"))
            save_vtk(
                path_vtk,
                (Array(xvi[1]), Array(xvi[2])),
                (Array(xci[1]), Array(xci[2])),
                data_v,
                data_c,
                velocity_v;
                t = t,
                pvd=joinpath(vis.vtk_dir, vis.pvd_name)
            )
            # Optional particle point-cloud output (large files).
            if vis.save_particle_points && (it == 1 || rem(it, vis.particle_vtk_every) == 0)
                save_particles(
                    particles,
                    pPhases;
                    fname = joinpath(vis.vtk_dir, "particles_" * lpad("$it", 6, "0")),
                    t = t,
                )
            end
            

        end

        if vis.pictures == true && (it == 1 || rem(it, vis.picture_every) == 0)
            # Make particles plottable
            p = particles.coords
            ppx, ppy = p
            pxv = Array(ppx.data[:] ./ 1.0e3)
            pyv = Array(ppy.data[:] ./ 1.0e3)
            clr = Array(pPhases.data[:])
            # clr      = pT.data[:]
            idxv = Array(particles.index.data[:])

            # --- New figure: velocity with limited range ---
            vmin = -0.1 / (365.25 * 24 * 3600)     # -10 cm/yr
            vmax = +0.1 / (365.25 * 24 * 3600)     # +10 cm/yr
            velocity2vertex!(vis.Vx_v, vis.Vy_v, @velocity(stokes)...)
            Vx_limited = clamp.(Array(vis.Vx_v), vmin, vmax)
            Vy_limited = clamp.(Array(vis.Vy_v), vmin, vmax)


            # --- ZOOM REGION ---
            xmin_zoom, xmax_zoom = 800, 1100
            ymin_zoom, ymax_zoom = -80, 0

            # --- FULL DOMAIN ---
            xmin_full = minimum(xvi[1]) * 1.0e-3
            xmax_full = maximum(xvi[1]) * 1.0e-3
            ymin_full = minimum(xvi[2]) * 1.0e-3
            ymax_full = maximum(xvi[2]) * 1.0e-3

            # FULL DOMAIN FIGURE
            make_figure(
                it, t, dt_plot,
                xvi, xci,
                T_buffer, ρg,
                stokes,
                Vx_limited, Vy_limited,
                pxv, pyv, clr, idxv,
                xmin_full, xmax_full, ymin_full, ymax_full,
                joinpath(vis.figdir, "full", "full_$(lpad(it, 2, "0")).png"), version=version
            )

            # ZOOMED FIGURE
            make_figure(
                it, t, dt_plot,
                xvi, xci,
                T_buffer, ρg,
                stokes,
                Vx_limited, Vy_limited,
                pxv, pyv, clr, idxv,
                xmin_zoom, xmax_zoom, ymin_zoom, ymax_zoom,
                joinpath(vis.figdir, "zoom", "zoom_$(lpad(it, 2, "0")).png"), version=version
            )
        end
        @show it += 1
        t += dt
    end
        # ------------------------------
    return nothing
end

## END OF MAIN SCRIPT ----------------------------------------------------------------
version = get(ENV, "SLURM_JOB_NAME", "unknown_version")
# version = "v0.357_computedensfix_interactivenode"
println("version is $version")

# MODEL SETUP
nx, ny = 900, 220
# Choose grid type: original uniform grid (ref_grid=0) or non-uniform logistic grid (ref_grid=1)
ref_grid = 1 # 0: original uniform grid, 1: non-uniform logistic grid

# GENERATE GRID
li, origin, phases_GMG, T_GMG, xvi, xci = GMG_subduction_2D_with_coords(
    nx + 1,
    ny + 1;
    ref_grid = ref_grid,
)

# Initialize MPI grid (or not)
igg = if !(JustRelax.MPI.Initialized()) # initialize (or not) MPI grid
    IGG(init_global_grid(nx, ny, 1; init_MPI = true)...)
else
    igg
end
function edge_pad(x::AbstractVector, lo_val, hi_val)
    lo = similar(x, 1); lo .= lo_val
    hi = similar(x, 1); hi .= hi_val
    return vcat(lo, x, hi)
end

function velocity_grids_gpu(xci, xvi, di)
    dxW = sum(@view di[1][1:1]);   dyW = sum(@view di[2][1:1])
    dxE = sum(@view di[1][end:end]); dyE = sum(@view di[2][end:end])

    x0 = sum(@view xci[1][1:1]);   xN = sum(@view xci[1][end:end])
    y0 = sum(@view xci[2][1:1]);   yN = sum(@view xci[2][end:end])

    xghost = edge_pad(xci[1], x0 - dxW, xN + dxE)
    yghost = edge_pad(xci[2], y0 - dyW, yN + dyE)

    grid_vx = xvi[1], yghost
    grid_vy = xghost, xvi[2]

    return grid_vx, grid_vy
end
main(
    li,
    origin,
    phases_GMG,
    T_GMG,
    igg;
    xvi = xvi,
    xci = xci,
    nx = nx,
    ny = ny,
    version = version,
    ref_grid = ref_grid,
);

# Load script dependencies
using GeoParams, CairoMakie
include("../utils/visualisation.jl")

const isCUDA = true

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
include("../Subduction2D_SZU2019/Subduction2D_setup_nonuniform.jl")
include("../Subduction2D_SZU2019/Subduction2D_rheology.jl")

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
# PREPARE VISUALIZATION SETTINGS
function prepare_visualisation(ni; version=nothing)
    # SETTINGS FOR VISUALIZATION AND OUTPUT
    do_vtk   = true # set to true to generate VTK files for ParaView
    pictures = true # set to true to generate PNG figures of particles and fields using Makie
    # IF VTK OUTPUT YES
    pvd_name = "Subduction2D"
    figdir   = "Subduction2D_SZU2019/Figures/Subduction2D_nonuniform/$version"
    save_particle_points = false # set to true to save particle point clouds as VTK files (can generate large files)
    vtk_every = 1 # save VTK every N iterations
    particle_vtk_every = 1 # save particle VTK every N iterations


    if do_vtk == true
        vtk_dir = joinpath(figdir, "vtk")
        if isfile(joinpath(vtk_dir, "$pvd_name.pvd"))
            rm(joinpath(vtk_dir, "$pvd_name.pvd"))
        end
        take(vtk_dir)
        checkpoint = joinpath(figdir, "checkpoint")
        take(checkpoint)
    end
    vis=(;do_vtk,vtk_dir,pvd_name ,figdir,save_particle_points,vtk_every,particle_vtk_every,pictures,checkpoint,Vx_v = @zeros(ni .+ 1...), Vy_v = @zeros(ni .+ 1...),)

    return vis
end
# VELOCITY BOXES ROUTINES
@parallel_indices (i, j) function _apply_vel_box_Vx!(
    Vx,
    xvx,
    yvx,
    cenx,
    cenz,
    halfx,
    halfz,
    vx_val,
)
    if i ≤ size(Vx, 1) && j ≤ size(Vx, 2)
        x = xvx[i]
        z = yvx[j]
        if abs(x - cenx) ≤ halfx && abs(z - cenz) ≤ halfz
            @inbounds Vx[i, j] = vx_val
        end
    end
    return nothing
end

@parallel_indices (i, j) function _apply_vel_box_Vy!(
    Vy,
    xvy,
    yvy,
    cenx,
    cenz,
    halfx,
    halfz,
    vy_val,
)
    if i ≤ size(Vy, 1) && j ≤ size(Vy, 2)
        x = xvy[i]
        z = yvy[j]
        if abs(x - cenx) ≤ halfx && abs(z - cenz) ≤ halfz
            @inbounds Vy[i, j] = vy_val
        end
    end
    return nothing
end

@parallel_indices (i, j) function _mark_vbox_mask_Vx!(
    mask_vbox_x,
    xvx,
    yvx,
    cenx,
    cenz,
    halfx,
    halfz,
)
    if i ≤ size(mask_vbox_x, 1) && j ≤ size(mask_vbox_x, 2)
        # mask indices (i,j) correspond to velocity DoFs at (i+1,j+1)
        ii = i + 1
        jj = j + 1
        if ii ≤ length(xvx) && jj ≤ length(yvx)
            x = xvx[ii]
            z = yvx[jj]
            if abs(x - cenx) ≤ halfx && abs(z - cenz) ≤ halfz
                @inbounds mask_vbox_x[i, j] = 1
            end
        end
    end
    return nothing
end

@parallel_indices (i, j) function _mark_vbox_mask_Vy!(
    mask_vbox_y,
    xvy,
    yvy,
    cenx,
    cenz,
    halfx,
    halfz,
)
    if i ≤ size(mask_vbox_y, 1) && j ≤ size(mask_vbox_y, 2)
        # mask indices (i,j) correspond to velocity DoFs at (i+1,j+1)
        ii = i + 1
        jj = j + 1
        if ii ≤ length(xvy) && jj ≤ length(yvy)
            x = xvy[ii]
            z = yvy[jj]
            if abs(x - cenx) ≤ halfx && abs(z - cenz) ≤ halfz
                @inbounds mask_vbox_y[i, j] = 1
            end
        end
    end
    return nothing
end

# Velocity boxes are applied on the same staggered coordinates as the Stokes solver.
# In the new Geometry API these coordinates are stored in `grid.xi_vel`:
# - `grid.xi_vel[1]` are the coordinates for Vx (x-face, z)
# - `grid.xi_vel[2]` are the coordinates for Vy (x, z-face)
# so the box region is applied to the correct velocity DoFs.
function apply_vel_boxes!(
    stokes,
    grid,
    boxes::Vector{VelBox2D},
)
    isempty(boxes) && return nothing

    Vx, Vy = @velocity(stokes)
    grid_vx, grid_vy = grid.xi_vel
    xvx, yvx = grid_vx
    xvy, yvy = grid_vy

    # reset velocity-box masks: 0 ⇒ no box (free)
    stokes.mask_vbox_x.mask .= 0
    stokes.mask_vbox_y.mask .= 0

    for box in boxes
        halfx = box.widthx / 2
        halfz = box.widthz / 2

        if box.has_vx
            nx = length(xvx)
            ny = length(yvx)
            @parallel (@idx (nx, ny)) _apply_vel_box_Vx!(
                Vx, xvx, yvx, box.cenx, box.cenz, halfx, halfz, box.vx
            )
            @parallel (@idx (nx, ny)) _mark_vbox_mask_Vx!(
                stokes.mask_vbox_x.mask,
                xvx,
                yvx,
                box.cenx,
                box.cenz,
                halfx,
                halfz,
            )
        end

        if box.has_vy
            nx = length(xvy)
            ny = length(yvy)
            @parallel (@idx (nx, ny)) _apply_vel_box_Vy!(
                Vy, xvy, yvy, box.cenx, box.cenz, halfx, halfz, box.vy
            )
            @parallel (@idx (nx, ny)) _mark_vbox_mask_Vy!(
                stokes.mask_vbox_y.mask,
                xvy,
                yvy,
                box.cenx,
                box.cenz,
                halfx,
                halfz,
            )
        end
    end

    return nothing
end
## END OF HELPER FUNCTION ------------------------------------------------------------

## BEGIN OF MAIN SCRIPT --------------------------------------------------------------
function main(
    li,
    origin,
    phases_GMG,
    igg;
    xvi,
    xci,
    nx = 16,
    ny = 16,
    ref_grid = 0,
    version = nothing,
)

    # Physical domain ------------------------------------
    ni = nx, ny           # number of cells
    di = @. li / ni       # grid steps
    grid = Geometry(ni, li; origin = origin)
    # non-uniform grid with refinement
    if ref_grid == 1
        grid = Geometry(
            PTArray(backend_JR),
            xvi...,
        )
    end
    (; xci, xvi) = grid # nodes at the center and vertices of the cells
    di_min = min(
        min(minimum.(grid.di.center)...),
        min(minimum.(grid.di.vertex)...),
        )
    di1 = grid.di 
    # ----------------------------------------------------
    # Set flags and parameters for visualization and output and create folders for output
    vis = prepare_visualisation(ni, version=version)
    # Physical properties using GeoParams ----------------
    rheology = init_rheologies()
    dt = 25.0e3 * 3600 * 24 * 365 # diffusive CFL timestep limiter
    dt_max = 25.0e3 * 3600 * 24 * 365 # diffusive CFL timestep limiter
    # ----------------------------------------------------

    # Initialize particles -------------------------------
    nxcell = 40
    max_xcell = 60
    min_xcell = 20
    particles = init_particles(
        backend_JP, nxcell, max_xcell, min_xcell, grid.xi_vel...
    )
    subgrid_arrays = SubgridDiffusionCellArrays(particles)
    # grid_vxi = velocity_grids(xci, xvi, di)
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
    if ref_grid == 0
        stokes.P .= PTArray(backend)(reverse(cumsum(reverse((ρg[2]) .* di[2], dims = 2), dims = 2), dims = 2))
    else
    # Lithostatic pressure integrates vertical body force using local cell dy (vertex spacing).
    stokes.P .= PTArray(backend)(reverse(cumsum(reverse((ρg[2]) .* reshape(di1.vertex[2], 1, :), dims = 2), dims = 2), dims = 2))
    end
    # Rheology
    args0 = (T = thermal.Tc, P = stokes.P, dt = Inf)
    viscosity_cutoff = (1.0e18, 1.0e23)
    compute_viscosity!(stokes, phase_ratios, args0, rheology, viscosity_cutoff)
    center2vertex!(stokes.viscosity.ηv, stokes.viscosity.η)
    # ----------------------------------------------------

    # PT coefficients for thermal diffusion
    pt_thermal = PTThermalCoeffs(
        backend, rheology, phase_ratios, args0, dt, ni, di1.vertex, li; ϵ = 1.0e-8, CFL = 0.95 / √2
    )

    # Boundary conditions
    flow_bcs = VelocityBoundaryConditions(;
        free_slip = (left = true, right = true, top = true, bot = true),
        free_surface = false,
    )
    flow_bcs!(stokes, flow_bcs) # apply boundary conditions
    update_halo!(@velocity(stokes)...)

    T_buffer = @zeros(ni .+ 1)
    Told_buffer = similar(T_buffer)
    dt₀ = similar(stokes.P)
    for (dst, src) in zip((T_buffer, Told_buffer), (thermal.T, thermal.Told))
        copyinn_x!(dst, src)
    end
    grid2particle!(pT, T_buffer, particles)

    τxx_v = @zeros(ni .+ 1...)
    τyy_v = @zeros(ni .+ 1...)

    dyrel = DYREL(backend, stokes, rheology, phase_ratios, di1, dt; ϵ = 1.0e-3)

    # Time loop
    t, it = 0.0, 0
    while it < 1000 # run only for 5 Myrs

        # interpolate fields from particle to grid vertices
        particle2grid!(T_buffer, pT, particles)
        @views T_buffer[:, end] .= Ttop
        @views T_buffer[:, 1] .= Tbot
        @views thermal.T[2:(end - 1), :] .= T_buffer
        thermal_bcs!(thermal, thermal_bc)
        temperature2center!(thermal)

        # interpolate stress back to the grid
        stress2grid!(stokes, pτ, particles)

        # Prescribe velocity boxes before solve so solver finds a solution consistent with them
        apply_vel_boxes!(stokes, grid, vel_boxes_2D)
        update_halo!(@velocity(stokes)...)

        # Stokes solver ----------------
        args = (; T = thermal.Tc, P = stokes.P, dt = Inf)
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
                    apply_velocity_box = stokes -> apply_vel_boxes!(stokes, grid, vel_boxes_2D),
                    viscosity_cutoff = (1.0e18, 1.0e23),
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
        subgrid_diffusion!(
            pT, thermal.T, thermal.ΔT, subgrid_arrays, particles, dt
        )
        # ------------------------------

        # Advection --------------------
        # advect particles in space
        advection_MQS!(particles, RungeKutta2(), @velocity(stokes), dt)
        # advect particles in memory
        move_particles!(particles, particle_args)
        # check if we need to inject particles
        # need stresses on the vertices for injection purposes
        # center2vertex!(τxx_v, stokes.τ.xx)
        # center2vertex!(τyy_v, stokes.τ.yy)
        inject_particles_phase!(
            particles,
            pPhases,
            particle_args_reduced,
            (T_buffer, stokes.τ.xx_v, stokes.τ.yy_v, stokes.τ.xy, stokes.ω.xy)
        )

        # update phase ratios
        update_phase_ratios!(phase_ratios, particles, pPhases)

        @show it += 1
        t += dt

        ### PARAVIEW PLOTTING
        if it >= 0 #it == 1 || rem(it, 5) == 0
            # checkpointing_jld2(vis.checkpoint, stokes, thermal, t, dt; it = it)
            # checkpointing_particles(vis.checkpoint, particles; phases = pPhases, phase_ratios = phase_ratios, particle_args = particle_args, particle_args_reduced = particle_args_reduced, t = t, dt = dt, it = it)
            (; η_vep, η) = stokes.viscosity
            if vis.do_vtk && (it == 1 || rem(it, vis.vtk_every) == 0)
                Vx_v = vis.Vx_v
                Vy_v = vis.Vy_v
                velocity2vertex!(Vx_v, Vy_v, @velocity(stokes)...)
                # Reconstruct compact phase "shapes" on the grid from particle phase ratios.
                phase_vertex = [argmax(p) for p in Array(phase_ratios.vertex)]
                Rx_c = zeros(size(stokes.P))
                Ry_c = zeros(size(stokes.P))
                @views Rx_c[axes(stokes.R.Rx, 1), axes(stokes.R.Rx, 2)] .= Array(stokes.R.Rx)
                @views Ry_c[axes(stokes.R.Ry, 1), axes(stokes.R.Ry, 2)] .= Array(stokes.R.Ry)

                data_v = (;
                    T = Array(T_buffer),
                    τII = Array(stokes.τ.II),
                    εII = Array(stokes.ε.II),
                    Vx = Array(Vx_v),
                    Vy = Array(Vy_v),
                    phase_vertex = phase_vertex,
                )
                data_c = (;
                    P   = Array(stokes.P),
                    η   = Array(η_vep),
                    Rx  = Array(Rx_c),
                    Ry  = Array(Ry_c),
                    Rmag = sqrt.(Rx_c .^ 2 .+ Ry_c .^ 2),
                )
                velocity_v = (
                    Array(Vx_v),
                    Array(Vy_v),
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

            if vis.pictures == true
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
                    joinpath(vis.figdir, "full_$(lpad(it, 2, "0")).png"), version=version
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
                    joinpath(vis.figdir, "zoom_$(lpad(it, 2, "0")).png"), version=version
                )
            end
        end
        # ------------------------------

    end

    return nothing
end

## END OF MAIN SCRIPT ----------------------------------------------------------------
version = "v0.274_newresolution_also_v5cmyr"
println("version is $version")
# MODEL SETUP
# n = 256
# nx, ny = n * 4, n
n = 32
nx, ny = n * 10, round(Int, n * 1.5 * 1.5) # increased vertical size by 50%
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

main(
    li,
    origin,
    phases_GMG,
    igg;
    xvi,
    xci,
    nx = nx,
    ny = ny,
    version = version,
    ref_grid = ref_grid,
);

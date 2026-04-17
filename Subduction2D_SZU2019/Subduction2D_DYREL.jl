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

const isCUDA = true

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
        xmin, xmax, ymin, ymax,
        figpath; version = nothing
    )

    # Add isotherms
    isotherms_C = [10, 50, 100, 150, 350, 450, 900, 1300]
    isotherms_K = isotherms_C .+ 273
    
    # Vertex grid (for T, V, etc.)
    xv = xvi[1][2:end-1] .* 1e-3
    yv = xvi[2] .* 1e-3
    Tv = Array(T_buffer[2:end-1, :])

    xi_mask = xmin .≤ xv .≤ xmax
    yi_mask = ymin .≤ yv .≤ ymax
    xv_zoom = xv[xi_mask]
    yv_zoom = yv[yi_mask]
    Tv_zoom = Tv[xi_mask, yi_mask]

    # Cell-centered grid (for density, stress, viscosity)
    xc = xci[1] .* 1e-3
    yc = xci[2] .* 1e-3
    ρ  = Array(ustrip.(ρg[2])) ./ 9.81
    isograds_density = [2700, 2800, 2900, 3100, 3200, 3250, 3300, 3350]

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

    fig[0, :] = GridLayout()

    Label(
        fig[0, 1],
        "Subduction2D - SZU2019 | version: $version | t = $t_Myr_round Myr | dt = $t_Kyr Kyr | it = $it",
        fontsize = 28,
        halign = :center,
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
    ax8 = Axis(fig[2, 5], aspect = ar, title = "log10(εII_pl)")
    ax9 = Axis(fig[3, 5], aspect = ar, title = "log10(η_vep)")

    # Apply zoom limits
    for ax in (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9)
        xlims!(ax, xmin, xmax)
        ylims!(ax, ymin, ymax)
    end

    # ---- Plots ----
    # Material phase
    h1 = scatter!(ax1, pxv[idxv], pyv[idxv], color = clr[idxv], markersize = 1)

    contour!(
        ax1,
        xv_zoom,
        yv_zoom,
        Tv_zoom;
        levels = isotherms_K,
        labels=true,
        labelsize = 12,
        color = :white,
        linewidth = 1.5
    )

    # Temperature
    h2 = heatmap!(
        ax2,
        xc,
        yc,
        Array(T_buffer[2:(end - 1), :]),
        colorrange = (253, 1718),
        colormap = :vik
    )

    contour!(
        ax2,
        xv_zoom,
        yv_zoom,
        Tv_zoom;
        levels = isotherms_K,
        labels=true,
        labelsize = 12,
        color = :white,
        linewidth = 1.5
    )


    # Density
    h3 = heatmap!(
        ax3,
        xc,
        yc,
        ρ;
        colormap = :vik,
        colorrange = (2500, 3300)
    )

    contour!(
        ax3,
        xc,
        yc,
        ρ;
        levels = isograds_density,
        labels=true,
        labelsize = 12,
        color = :white,
        linewidth = 1.5
    )

    # Velocity x-direction
    h4 = heatmap!(
        ax4,
        xv,
        yv,
        Vx_limited;
        colormap = :vik,
        colorrange = (vmin, vmax)
    )

    # Velocity y-direction
    h5 = heatmap!(
        ax5,
        xv,
        yv,
        Vy_limited;
        colormap = :vik,
        colorrange = (vmin, vmax)
    )

    # Second invariant of the stress tensor
    h6 = heatmap!(
        ax6,
        xc,
        yc,
        Array(log10.(stokes.τ.II)),
        colormap = :batlow
    )

    # Second invariant of the strain rate tensor
    h7 = heatmap!(
        ax7,
        xc,
        yc,
        Array(log10.(stokes.ε.II));
        colormap = :vik,
        colorrange = (ε_min, ε_max)
    )

    contour!(
        ax7,
        xv_zoom,
        yv_zoom,
        Tv_zoom;
        levels = isotherms_K,
        labels=true,
        labelsize = 12,
        color = :white,
        linewidth = 1.5
    )

    # plastic strain rate
    h8 = heatmap!(
        ax8,
        xc,
        yc,
        Array(log10.(stokes.EII_pl));
        colormap = :vik,
        colorrange = (ε_min, ε_max)
    )
    contour!(
        ax8,
        xv_zoom,
        yv_zoom,
        Tv_zoom;
        levels = isotherms_K,
        labels=true,
        labelsize = 12,
        color = :white,
        linewidth = 1.5
    )

    # Visco_plasto_elastic viscosity
    h9 = heatmap!(
        ax9,
        xc,
        yc,
        Array(log10.(stokes.viscosity.η_vep));
        colormap = :vik,
        colorrange = (visc_min, visc_max)
    )
    contour!(
        ax9,
        xv,
        yv,
        Tv;
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

# Velocity boxes are applied on the same staggered grid as the Stokes solver:
# - Vx lives at (x face, z) with size (ni[1]+1, ni[2]+2); grid_vx = (xvi[1], yVx).
# - Vy lives at (x, z face) with size (ni[1]+2, ni[2]+1); grid_vy = (xVy, xvi[2]).
# There are no separate "cell-center" or "vertex" velocity arrays; Vx and Vy are the
# only velocity DoFs, and apply_vel_boxes! correctly uses grid.grid_v so the box
# region is applied to the right indices.
function apply_vel_boxes!(
    stokes,
    grid::Geometry{2, T},
    boxes::Vector{VelBox2D},
) where {T}
    isempty(boxes) && return nothing

    Vx, Vy = @velocity(stokes)
    grid_vx, grid_vy = grid.grid_v
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
function main(li, origin, phases_GMG, igg; nx = 16, ny = 16, figdir = "figs2D", do_vtk = false, version = nothing)

    # Physical domain ------------------------------------
    ni = nx, ny           # number of cells
    di = @. li / ni       # grid steps
    grid = Geometry(ni, li; origin = origin)
    (; xci, xvi) = grid # nodes at the center and vertices of the cells
    # ----------------------------------------------------

    # Physical properties using GeoParams ----------------
    rheology = init_rheologies()
    dt = 0.5e3 * 3600 * 24 * 365.25
    dt_max = 10e3 * 3600 * 24 * 365.25 # diffusive CFL timestep limiter
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

    # particle fields for the stress rotation
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
    # ----------------------------------------------------

    # TEMPERATURE PROFILE --------------------------------
    Ttop = -10 + 273
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
    viscosity_cutoff = (5.0e19, 5.0e23)
    compute_viscosity!(stokes, phase_ratios, args0, rheology, viscosity_cutoff)
    center2vertex!(stokes.viscosity.ηv, stokes.viscosity.η)
    # ----------------------------------------------------

    # PT coefficients for thermal diffusion
    pt_thermal = PTThermalCoeffs(
        backend, rheology, phase_ratios, args0, dt, ni, di, li; ϵ = 1.0e-8, CFL = 0.95 / √2
    )

    # Boundary conditions
    flow_bcs = VelocityBoundaryConditions(;
        free_slip = (left = false, right = false, top = false, bot = false),
        no_slip = (left = true, right = true, top = true, bot = true),
        free_surface = false,
    )
    flow_bcs!(stokes, flow_bcs) # apply boundary conditions
    update_halo!(@velocity(stokes)...)

    # IO -------------------------------------------------
    # if it does not exist, make folder where figures are stored
    if do_vtk
        vtk_dir = joinpath(figdir, "vtk")
        take(vtk_dir)
        checkpoint = joinpath(figdir, "checkpoint")
        take(checkpoint)
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

    dyrel = DYREL(backend, stokes, rheology, phase_ratios, di, dt; ϵ = 1.0e-3)

    # Time loop
    t, it = 0.0, 0
    while it < 1000

        # interpolate fields from particle to grid vertices
        particle2grid!(T_buffer, pT, xvi, particles)
        @views T_buffer[:, end] .= Ttop
        @views T_buffer[:, 1] .= Tbot
        @views thermal.T[2:(end - 1), :] .= T_buffer
        thermal_bcs!(thermal, thermal_bc)
        temperature2center!(thermal)

        # interpolate stress back to the grid
        stress2grid!(stokes, pτ, xvi, xci, particles)

        # Prescribe velocity boxes before solve so solver finds a solution consistent with them
        apply_vel_boxes!(stokes, grid, vel_boxes_2D)
        update_halo!(@velocity(stokes)...)

        # Stokes solver: re-apply velocity boxes after every V update so the solver
        # keeps the prescribed velocities in the box (otherwise each iteration overwrites them).

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
                di,
                dt,
                igg;
                kwargs = (;
                    verbose_PH = true,
                    verbose_DR = false,
                    iterMax = 50.0e3,
                    rel_drop = 1.0e-2,
                    nout = 400,
                    λ_relaxation_PH = 1,
                    λ_relaxation_DR = 1,
                    viscosity_relaxation = 2.0e-2,
                    viscosity_cutoff = viscosity_cutoff,
                    apply_velocity_box = stokes -> apply_vel_boxes!(stokes, grid, vel_boxes_2D),
                )
            )
        end

        # Enforce velocity boxes again so the final velocity field (used for advection) matches the prescription
        apply_vel_boxes!(stokes, grid, vel_boxes_2D)
        update_halo!(@velocity(stokes)...)
        # print some stuff
        println("Stokes solver time             ")
        println("   Total time:      $t_stokes s")
        t_Kyr = round(dt / (3600 * 24 * 365.25 * 1000); digits = 2)
        println("Timestep:           $t_Kyr")
        # rotate stresses
        rotate_stress!(pτ, stokes, particles, xci, xvi, dt)
        # compute time step
        dt_plot = dt
        dt = compute_dt(stokes, di, dt_max) * 0.8
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
            di;
            kwargs = (
                igg = igg,
                phase = phase_ratios,
                iterMax = 50.0e3,
                nout = 2.0e2,
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
        advection_MQS!(particles, RungeKutta2(), @velocity(stokes), grid_vxi, dt)
        # advect particles in memory
        move_particles!(particles, xvi, particle_args)
        # check if we need to inject particles
        # need stresses on the vertices for injection purposes
        # center2vertex!(τxx_v, stokes.τ.xx)
        # center2vertex!(τyy_v, stokes.τ.yy)
        inject_particles_phase!(
            particles,
            pPhases,
            particle_args_reduced,
            (T_buffer, stokes.τ.xx_v, stokes.τ.yy_v, stokes.τ.xy, stokes.ω.xy),
            xvi
        )

        # update phase ratios
        update_phase_ratios!(phase_ratios, particles, xci, xvi, pPhases)

        @show it += 1
        t += dt

        # Data I/O and plotting ---------------------
        if it == 1 || it == 2 || rem(it, 5) == 0
            # checkpointing_jld2(checkpoint, stokes, thermal, t, dt; it = it)
            # checkpointing_particles(checkpoint, particles; phases = pPhases, phase_ratios = phase_ratios, particle_args = particle_args, particle_args_reduced = particle_args_reduced, t = t, dt = dt, it = it)
            (; η_vep, η) = stokes.viscosity
            if do_vtk
                velocity2vertex!(Vx_v, Vy_v, @velocity(stokes)...)
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
            pxv = Array(ppx.data[:] ./ 1.0e3)
            pyv = Array(ppy.data[:] ./ 1.0e3)
            clr = Array(pPhases.data[:])
            # clr      = pT.data[:]
            idxv = Array(particles.index.data[:])

            # --- New figure: velocity with limited range ---
            vmin = -0.1 / (365.25 * 24 * 3600)     # -10 cm/yr
            vmax = +0.1 / (365.25 * 24 * 3600)     # +10 cm/yr
            Vx_raw = Array(stokes.V.Vx[2:(end - 1), :])
            Vx_limited = clamp.(Vx_raw, vmin, vmax)
            Vy_raw = Array(stokes.V.Vy[2:(end - 1), :])
            Vy_limited = clamp.(Vy_raw, vmin, vmax)


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
                joinpath(figdir, "full_$(lpad(it, 2, "0")).png"), version=version
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
                joinpath(figdir, "zoom_$(lpad(it, 2, "0")).png"), version=version
            )


            # --- FULL DOMAIN ---
            xmin_full = minimum(xvi[1]) * 1.0e-3
            xmax_full = maximum(xvi[1]) * 1.0e-3
            ymin_full = minimum(xvi[2]) * 1.0e-3
            ymax_full = maximum(xvi[2]) * 1.0e-3


    end

        # ------------------------------
    end 
    return nothing
end

## END OF MAIN SCRIPT ----------------------------------------------------------------
do_vtk = true # set to true to generate VTK files for ParaView
version = "v0.261_noslip_andhigherdt_cutoff_5e19_resx2"
figdir = "Subduction2D_SZU2019/Figures/Subduction2D_DYREL/dyrel_$version"
println(version)
# n=144
# n = 80
# n = 32
nx, ny = round(Int, 2*1466) , round(Int, 2*270)
# nx, ny = n * 10, round(Int, n * 1.5)

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


main(li, origin, phases_GMG, igg; figdir = figdir, nx = nx, ny = ny, do_vtk = do_vtk, version = version);

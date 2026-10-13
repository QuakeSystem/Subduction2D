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
const backend_JP = @static if isCUDA
    CUDA.CUDABackend # Options: JustPIC.CPU, CUDA.CUDABackend, AMDGPU.ROCBackend
else
    JustPIC.CPU # Options: JustPIC.CPU, CUDA.CUDABackend, AMDGPU.ROCBackend
end
# GPU-safe copy for CellArrays (works around JustPIC 0.7.14 resample! scalar indexing)
Base.copy(x::JustPIC.CellArrays.CellArray) = typeof(x)(copy(x.data), size(x))

# Workaround JustPIC 0.7.14 / GridGeometryUtils 0.3.1: recentering loses ~1e-10 m on
# refined grids, so chain endpoints miss the cell edges. Every caller passes a chord
# spanning the full cell width, so snap its x-ends onto the side edges exactly.
const GGU = JustPIC.GridGeometryUtils
@inline function JustPIC.recenter_on_cell(s::GGU.Segment, r::GGU.Rectangle{T}) where {T}
    oy = r.center[2]
    hl = r.l / 2
    s_local = GGU.Segment(GGU.Point(-hl, s.p1[2] - oy), GGU.Point(hl, s.p2[2] - oy))
    return s_local, GGU.Rectangle((zero(T), zero(T)), r.l, r.h; θ = zero(T))
end
# JustRelax's markerchain phase correction (update_phases_given_markerchain!) assumes a
# uniform grid twice: it picks rows with a uniform dy (misses the surface) and finds the
# chain column with a uniform-dx `cell_index` (compares particles against the wrong part
# of the chain, empties cells -> NaN). Patch both. The GPU path runs its own copy of these
# functions inside the CUDA extension module, so patch that module as well.
const JR2D_modules = @static if isCUDA
    (JustRelax.JustRelax2D, Base.get_extension(JustRelax, :JustRelaxCUDAExt).JustRelax2D)
else
    (JustRelax.JustRelax2D,)
end
for M in JR2D_modules
    @eval begin
        # search all rows (clamped inside the kernel)
        $M.find_minmax_cell_indices(chain_yi, origin_y, di) = (11, 10^9)
        # same function, but with JustPIC's refined-grid (bisection) column lookup
        @inline function $M.is_above_chain(xq, yq, coords, cell_vertices)
            I = JustPIC.parent_cell_index(xq, cell_vertices, JustPIC.midpoint_seed(cell_vertices))
            x_cell, y_cell = coords[1][I], coords[2][I]
            ychain = if 1 < I[1] < length(cell_vertices) - 1
                JustPIC.interp1D_inner(xq, x_cell, y_cell, coords, I)
            else
                JustPIC.interp1D_extremas(xq, x_cell, y_cell)
            end
            return yq > ychain
        end
    end
end

# Load script dependencies
using GeoParams, CairoMakie, LinearAlgebra, HDF5, PoissonGrids
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

"""
    topo_y(x_km)

Initial air/rock interface matching the "Air temperature" polygon in
`GMG_subduction_2D_with_coords`: flat at -12.5 km for x ≤ 854 km,
flat at -8 km for x ≥ 954 km, linearly ramped between them.
Input and output are both in km.
"""
function topo_y(x_km::Real)
    x_left, x_right = 854.0, 954.0
    z_left, z_right = -12.5, -8.0
    if x_km <= x_left
        return z_left
    elseif x_km >= x_right
        return z_right
    else
        t = (x_km - x_left) / (x_right - x_left)
        return z_left + t * (z_right - z_left)
    end
end

function particle_report(pT, particles, thermal, xvi; nshow = 15, label = "")
    xv, yv = Array(xvi[1]), Array(xvi[2])
    nx, ny = length(xv) - 1, length(yv) - 1
    xc = 0.5 .* (xv[1:(end - 1)] .+ xv[2:end])
    yc = 0.5 .* (yv[1:(end - 1)] .+ yv[2:end])
    # Ghost centroid *positions* only (mirrored across the boundary vertex, one
    # extra half-cell out), matching the layout of the padded (nx+2, ny+2)
    # `thermal.T` array. This placement is just geometry and holds for every
    # boundary alike -- it says nothing about what *value* ends up at the ghost.
    # That value comes from `T` below (the real, already-BC'd field): top/bot
    # are `constant_value` (Tghost = 2*face_value - Tint), left/right are
    # `no_flux` (Tghost = Tint, a plain copy) -- see thermal_bcs!/free_slip!.
    xcg = [2xv[1] - xc[1]; xc; 2xv[end] - xc[end]]
    ycg = [2yv[1] - yc[1]; yc; 2yv[end] - yc[end]]

    T = Array(thermal.T)                      # (nx+2, ny+2), matches xcg × ycg
    Tmin, Tmax = extrema(@view T[2:(end - 1), 2:(end - 1)])

    println("\n=== particle report: $label ===")
    println("particle fields: ", propertynames(particles))
    hasproperty(particles, :di)  && println("particles.di type: ", typeof(particles.di))
    hasproperty(particles, :xci) && println("particles.xci[1][1:3] = ", Array(particles.xci[1])[1:3],
                                            "   my ghosted xc[1:3] = ", xcg[1:3])
    same = size(particles.index) == (nx, ny)
    println("particle cell grid ", size(particles.index), "  vs grid cells ", (nx, ny))

        # copy particle data to host once (works for CPU and GPU)
    idxd = Array(particles.index.data)
    pxd  = Array(particles.coords[1].data)
    pyd  = Array(particles.coords[2].data)
    Td   = Array(pT.data)
    ncell = prod(size(particles.index))
    @assert size(idxd, 1) == ncell "unexpected CellArray layout: $(size(idxd))"
    CI = CartesianIndices(size(particles.index))

    # JustPIC pads its cell grid by one cell per side: particle cell (i, j) = grid cell (i-1, j-1)
    same = size(particles.index) == (nx + 2, ny + 2)

    rows = NamedTuple[]
    nact = nabove = nbelow = nwrong = nout = 0
    maxerr = 0.0
    for c in 1:ncell, ip in axes(idxd, 2)
        idxd[c, ip, 1] || continue
        nact += 1
        i, j = Tuple(CI[c])
        x, y, Tp = pxd[c, ip, 1], pyd[c, ip, 1], Td[c, ip, 1]
        interior = same && 2 <= i <= nx + 1 && 2 <= j <= ny + 1

        (xv[1] <= x <= xv[end] && yv[1] <= y <= yv[end]) || (nout += 1)
        interior && !(xv[i - 1] <= x <= xv[i] && yv[j - 1] <= y <= yv[j]) && (nwrong += 1)

        # reference bilinear interpolation on the ghosted centroid grid
        ic = clamp(searchsortedlast(xcg, x), 1, length(xcg) - 1)
        jc = clamp(searchsortedlast(ycg, y), 1, length(ycg) - 1)
        tx = (x - xcg[ic]) / (xcg[ic + 1] - xcg[ic])
        ty = (y - ycg[jc]) / (ycg[jc + 1] - ycg[jc])
        Tref = (1 - tx) * (1 - ty) * T[ic, jc]     + tx * (1 - ty) * T[ic + 1, jc] +
               (1 - tx) * ty       * T[ic, jc + 1] + tx * ty       * T[ic + 1, jc + 1]
        err = Tp - Tref
        maxerr = max(maxerr, abs(err))
        Tp > Tmax + 1e-6 && (nabove += 1)
        Tp < Tmin - 1e-6 && (nbelow += 1)

        if abs(err) > 1.0 || Tp > Tmax + 1e-6 || Tp < Tmin - 1e-6
            # Which face's ghost (if any) this interpolation bracket touches.
            # A corner interpolation can touch two faces at once, hence a Vector.
            touches = Symbol[]
            ic == 1               && push!(touches, :left)
            ic == length(xcg) - 1 && push!(touches, :right)
            jc == 1               && push!(touches, :bot)
            jc == length(ycg) - 1 && push!(touches, :top)
            push!(rows, (; x_km = x / 1e3, y_km = y / 1e3, cell = (i, j),
                dx_km = interior ? (xv[i] - xv[i - 1]) / 1e3 : NaN,
                dy_km = interior ? (yv[j] - yv[j - 1]) / 1e3 : NaN,
                Tp, Tref, err, touches))
        end
    end

    println("active particles: $nact | outside domain: $nout | in wrong cell: $nwrong")
    println("grid interior T range: ", (Tmin, Tmax))
    println("particles above grid max: $nabove | below grid min: $nbelow")
    nghost = count(r -> !isempty(r.touches), rows)
    println("max |pT - Tref| = $maxerr | flagged: ", length(rows),
            " (", round(100 * length(rows) / max(nact, 1), digits = 1), "% of active) | touching a ghost: ", nghost)
    if !isempty(rows)
        for face in (:left, :right, :top, :bot)
            n = count(r -> face in r.touches, rows)
            n > 0 && println("  touching $face ghost: $n")
        end
        println("flagged x range [km]: ", extrema(r.x_km for r in rows),
                "  y range [km]: ", extrema(r.y_km for r in rows))
        sort!(rows; by = r -> -abs(r.err))
        println("worst $nshow:")
        for r in rows[1:min(nshow, end)]
            tag = isempty(r.touches) ? "" : "  [$(join(r.touches, ","))]"
            println("  x=", round(r.x_km, digits = 2), " y=", round(r.y_km, digits = 2),
                    " cell=", r.cell, " dx=", round(r.dx_km, digits = 2), " dy=", round(r.dy_km, digits = 2),
                    " | pT=", round(r.Tp, digits = 1), " Tref=", round(r.Tref, digits = 1),
                    " err=", round(r.err, digits = 1), tag)
        end
    end
    return rows
end
## END OF HELPER FUNCTION ------------------------------------------------------------

## BEGIN OF MAIN SCRIPT --------------------------------------------------------------
function main(
    li,
    origin,
    phases_GMG,
    T_GMG,
    igg;
    vertices = nothing,
    nx = 16,
    ny = 16,
    ref_grid = 0,
    version = nothing,
)

    # Physical domain ------------------------------------
    ni = nx, ny            # number of cells
    di = @. li / ni        # grid step in x- and -y
    grid0 = Geometry(ni, li; origin = origin)
    di_min = @. li / ni       # grid steps
    # non-uniform grid with refinement
    grid = Geometry(
        PTArray(backend_JR),
        vertices...
        )
    di_min = min(
        min(minimum.(grid.di.center)...),
        min(minimum.(grid.di.vertex)...),
    )
    (; xci, xvi) = grid # nodes at the center and vertices of the cells
    # ----------------------------------------------------
    # visualisation subroutine
    vis = prepare_visualisation(ni, version=version)
    let
        dx, dz = Array(grid.di.vertex[1]), Array(grid.di.vertex[2])
        ar = [dx[i] / dz[j] for i in eachindex(dx), j in eachindex(dz)]
        println("max cell aspect ratio dx/dz = ", maximum(ar), " at ", Tuple(argmax(ar)))
    end
    # Physical properties using GeoParams ----------------
    rheology = init_rheologies_start()
    dt = 25.0e3 * 3600 * 24 * 365 # diffusive CFL timestep limiter
    dt_max = 25e3 * 3600 * 24 * 365 # diffusive CFL timestep limiter
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
        map(t -> Array.(t), grid.xi_vel)...;
        periodic = (false, false)
    )
    subgrid_arrays = SubgridDiffusionCellArrays(particles; loc = :center)
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
    # ### == Markerchain
    initial_elevation = -8.0e3
    n_chain = round(Int, nx+1)  # or whatever resolution you want for the interface
    chain = init_markerchain(backend_JP, nxcell, min_xcell, max_xcell, PTArray(backend)(xvi[1]), initial_elevation)

    # # chain = init_markerchain(backend_JP, nxcell, min_xcell, max_xcell, xvi[1], initial_elevation)
    air_phase = 2
    topo_y_vec = topo_y.(xvi[1] ./ 1.0e3) .* 1.0e3
    fill_chain_from_vertices!(chain, PTArray(backend)(topo_y_vec))
    @show extrema(Array(chain.h_vertices))
    @show size(Array(chain.h_vertices))
    @show extrema(topo_y_vec)
    # update_phases_given_markerchain!(pPhases, chain, particles, origin, grid.di, air_phase)
    # # end

    update_phase_ratios!(phase_ratios, particles, pPhases)
    # ----------------------------------------------------
    ϕ = RockRatio(backend, ni)
    # update_rock_ratio!(ϕ, phase_ratios, air_phase)
    compute_rock_fraction!(ϕ, chain, xvi, grid.di.vertex)

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
        # no_flux = (left = true, right = true, top = true, bot = true),
        no_flux = (left = true, right = true, top = false, bot = false),
        constant_value = (left = false, right = false, top = Ttop, bot = Tbot),
        # constant_value = (left = false, right = false, top = false, bot = false),
    )
    thermal_bcs!(thermal, thermal_bc)
    # ----------------------------------------------------

    # Buoyancy forces
    ρg = ntuple(_ -> @zeros(ni...), Val(2))
    compute_ρg!(ρg[2], phase_ratios, rheology, (T = thermal.T, P = stokes.P))
    # Lithostatic pressure: integrate ρg down each column into the cell-centered P.
    # `compute_lithostatic_pressure!` adds the -w/2 half-cell correction and tolerates
    # a scalar or a per-cell dz vector, so it works for both grid modes. A single
    # non-finite ρg cell no longer smears a NaN over the whole column the way the old
    # hand-rolled `reverse(cumsum(reverse(...)))` did.
    if ref_grid == 0
        compute_lithostatic_pressure!(stokes.P, ρg[2], di_min[2], igg)
    else
        # Non-uniform vertical spacing: pass the per-cell heights (length ny).
        compute_lithostatic_pressure!(stokes.P, ρg[2], grid.di.vertex[2], igg)
    end

    # Rheology
    args0 = (T = thermal.T, P = stokes.P, dt = Inf)
    viscosity_cutoff = (1.0e18, 1.0e24)
    compute_viscosity!(stokes, phase_ratios, args0, rheology, viscosity_cutoff)
    center2vertex!(stokes.viscosity.ηv, stokes.viscosity.η)
    # ----------------------------------------------------

    # PT coefficients for thermal diffusion
    pt_thermal = PTThermalCoeffs(
        backend, rheology, phase_ratios, args0, dt, ni, minimum.(grid.di.vertex), li; ϵ = 1.0e-12, CFL = 0.7 / √2
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
    centroid2particle!(pT, thermal.T, particles)
    # dt₀ = similar(thermal.T)
    # AI assisted implement. verify.
    dt₀_center = similar(stokes.P)   # plain (nx, ny), matches what subgrid_characteristic_time! expects
    dt₀         = similar(thermal.T) # ghost-padded (nx+2, ny+2), for centroid2particle! below

    # Time loop
    t, it = 0.0, 0

    dyrel = DYREL(backend, stokes, rheology, phase_ratios, grid.di, dt; ϵ = 1.0e-3)

    while t < 4.4e6 * (3600 * 24 * 365.25)  # run only for 4 Myrs
        if it == 4
            rheology = init_rheologies()
        end
        # Get flat views of the raw data
        phases_flat = pPhases.data[:]   # all particle phase values
        temps_flat  = pT.data[:]        # all particle temperatures
        index_flat  = particles.index.data[:]  # true = active particle
        air_mask = (phases_flat .== 2.0) .& index_flat
        pT.data[air_mask] .= 273.0
        pT_active = Array(pT.data)[Array(particles.index.data)]

        particle2centroid!(thermal.T, pT, particles)
        thermal_bcs!(thermal, thermal_bc)
        thermal.T .= clamp.(thermal.T, Ttop - 50.0, Tbot + 50.0)

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
                    verbose_DR = false,
                    iterMax = 50.0e2,
                    rel_drop = 1.0e-2,
                    nout = 400,
                    λ_relaxation_PH = 1,
                    λ_relaxation_DR = 1,
                    viscosity_relaxation = 1.0e-2,
                    viscosity_cutoff = viscosity_cutoff,
                    free_surface = true,
                    debug=false
                )
            )
        end

        println("Stokes solver time             ")
        println("   Total time:      $t_stokes s")
        # println("   Time/iteration:  $(t_stokes / out.iter) s")

        # rotate stresses
        rotate_stress!(pτ, stokes, particles, dt)
        # compute time step
        dt_plot = dt
        dt = compute_dt(stokes, di_min, dt_max)
        println("dt kyr= ", dt/365.25/24/3600/1000, "  finite = ", isfinite(dt), "  positive = ", dt > 0)

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
                iterMax = 50.0e5,
                nout = 1.0e2,
                verbose = false,
            )
        )
        thermal.T .= clamp.(thermal.T, Ttop - 50.0, Tbot + 50.0)


        subgrid_characteristic_time!(
            subgrid_arrays, particles, dt₀_center, phase_ratios, rheology, thermal, stokes
        )
        @views dt₀[2:end-1, 2:end-1] .= dt₀_center   # copy into the padded array at the correct offset
        @views dt₀[1, :]   .= dt₀[2, :]
        @views dt₀[end, :] .= dt₀[end - 1, :]
        @views dt₀[:, 1]   .= dt₀[:, 2]
        @views dt₀[:, end] .= dt₀[:, end - 1]
        centroid2particle!(subgrid_arrays.dt₀, dt₀, particles)

        subgrid_diffusion_centroid!(
            pT, thermal.T, thermal.ΔT, subgrid_arrays, particles, dt
        )
        pT.data   .= clamp.(pT.data, Ttop - 50.0, Tbot + 50.0)

        # ------------------------------
        # Advection --------------------
        # advect particles in space
        advection_MQS!(particles, RungeKutta2(), @velocity(stokes), dt)
        # advect particles in memory
        move_particles!(particles, particle_args)
        # Inject phase labels first, then initialize every newly injected particle field
        # through the regular centroid/vertex interpolation paths.
        inject_particles_phase!(particles, pPhases, (pT, ), (thermal.T,))

        centroid2particle!(pτ.τ_normal[1], stokes.τ.xx, particles; ghosted = false)
        centroid2particle!(pτ.τ_normal[2], stokes.τ.yy, particles; ghosted = false)
        grid2particle!(pτ.τ_shear[1], stokes.τ.xy, particles; ghost_1 = false, ghost_2 = false)
        grid2particle!(pτ.ω[1], stokes.ω.xy, particles; ghost_1 = false, ghost_2 = false)

        # # advect marker chain
        # semi-Lagrangian update with slope limiter: stops the chain forming a near-vertical
        # cliff at the plate boundary (x ≈ 854 km); mean height free to change
        semilagrangian_advection_markerchain!(chain, RungeKutta2(), @velocity(stokes), grid.xi_vel, xvi, dt;
            max_slope_angle = 30.0, conserve_mean = false)
        # advect_markerchain!(chain, RungeKutta2(), @velocity(stokes), grid.xi_vel, dt)
        update_phases_given_markerchain!(pPhases, chain, particles, origin, di, air_phase)
        update_phase_ratios!(phase_ratios, particles, pPhases)
        # update_phase_ratios!(phase_ratios, particles, pPhases)
        update_rock_ratio!(ϕ, phase_ratios, air_phase)

        compute_rock_fraction!(ϕ, chain, xvi, grid.di.vertex)        # if it == 1
        #     @show extrema(Array(chain.h_vertices))
        #     @show extrema(Array(xvi[1]))
        #     @show Array(chain.h_vertices)
        # end

        # update phase ratios
        particle_report(pT, particles, thermal, xvi; label = "end_of_time_step")
        let h = Array(chain.h_vertices), x = Array(chain.cell_vertices)
            s = atand.(abs.(diff(h) ./ diff(x)))
            println("chain max slope = ", maximum(s), "° at x = ", x[argmax(s)] / 1e3, " km")
        end


        @show it += 1
        t += dt
        ### PARAVIEW PLOTTING
        if it == 1 || rem(it, 25) == 0
            checkpointing_hdf5(vis.checkpoint, stokes, thermal.T, t, dt)
            # checkpointing_particles(vis.checkpoint, particles; phases = pPhases, phase_ratios = phase_ratios, particle_args = particle_args, particle_args_reduced = particle_args_reduced, t = t, dt = dt, it = it)

            # ppx, ppy = particles.coords
            # fname = joinpath(vis.figdir, "marker_$(lpad(it,6,'0')).h5")
            # HDF5.h5open(fname, "w") do file
            #     mg = HDF5.create_group(file, "MarkerGroup")
            #     mg["markx"] = Array(ppx.data[:])
            #     mg["marky"] = Array(ppy.data[:])
            #     mg["markt"] = Array(pPhases.data[:])

            #     modelg = HDF5.create_group(file, "ModelGroup")
            #     modelg["Model"] = [t]
            # end
        end
        (; η_vep, η) = stokes.viscosity
        save_marker_chain(
            joinpath(vis.vtk_dir, "chain_" * lpad("$it", 6, "0")),
            Array(chain.cell_vertices),
            Array(chain.h_vertices);
            pvd = joinpath(vis.vtk_dir, "MarkerChain.pvd"),
            t = t,
        )
        # save_marker_chain_hdf5(joinpath(vis.vtk_dir, "chain_$(lpad(it,3,'0'))"), chain; t = t)
        if vis.do_vtk && (it == 1 || rem(it, vis.vtk_every) == 0)
            velocity2vertex!(vis.Vx_v, vis.Vy_v, @velocity(stokes)...)
            # Reconstruct compact phase "shapes" on the grid from particle phase ratios.
            pv = Array(phase_ratios.vertex.data)        # (cells, nphases, 1)
            phase_vertex = reshape([argmax(@view pv[c, :, 1]) for c in axes(pv, 1)], size(phase_ratios.vertex))
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
                # T_buffer   = Array(T_buffer),
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
                thermal.T[2:end-1,2:end-1], ρg,
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
                thermal.T[2:end-1,2:end-1], ρg,
                stokes,
                Vx_limited, Vy_limited,
                pxv, pyv, clr, idxv,
                xmin_zoom, xmax_zoom, ymin_zoom, ymax_zoom,
                joinpath(vis.figdir, "zoom", "zoom_$(lpad(it, 2, "0")).png"), version=version
            )
        end
    end
        # ------------------------------
    return nothing
end

## END OF MAIN SCRIPT ----------------------------------------------------------------
version = get(ENV, "SLURM_JOB_NAME", "unknown_version")
# version = "v0.389_justpic0.7.14_interactive"
println("version is $version")

# MODEL SETUP
nx, ny = 1049, 269
# nx, ny = 3000, 520
# nx, ny = 300, 52

# Choose grid type: original uniform grid (ref_grid=0) or non-uniform logistic grid (ref_grid=1)
ref_grid = 1 # 0: original uniform grid, 1: non-uniform logistic grid

# GENERATE GRID
li, origin, phases_GMG, T_GMG, (xv, yv) = GMG_subduction_2D_with_coords(
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
    T_GMG,
    igg;
    vertices = (xv, yv),
    nx = nx,
    ny = ny,
    version = version,
    ref_grid = ref_grid,
);


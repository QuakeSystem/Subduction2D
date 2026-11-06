# PREPARE VISUALIZATION SETTINGS
function prepare_visualisation(ni; do_vtk=true, pictures = true,  pvd_name = "Subduction2D", data_folder = "Subduction2D_SZU2019/data", JRversion = "Subduction2D_pre_JRv0.8.1", version=nothing, save_particle_points = true, vtk_every = 5,  particle_vtk_every = 5, picture_every = 5)
    figdir   = joinpath(data_folder, JRversion, version)
    if do_vtk == true
        vtk_dir = joinpath(figdir, "vtk")
        if isfile(joinpath(vtk_dir, "$pvd_name.pvd"))
            rm(joinpath(vtk_dir, "$pvd_name.pvd"))
        end
        take(vtk_dir)
        checkpoint = joinpath(figdir, "checkpoint")
        take(checkpoint)
    end
    vis=(;do_vtk,vtk_dir,pvd_name ,figdir,save_particle_points,vtk_every,particle_vtk_every,pictures,picture_every,checkpoint,Vx_v = @zeros(ni .+ 1...), Vy_v = @zeros(ni .+ 1...),)

    return vis
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

"""
    open_surface_track(path)

Open the per-step surface-track file and write its header. The rows follow the SZU2019
marker-track column order (`eachdt_markertrack_lphr0.txt`), so the same loaders read both;
`#` header lines are skipped by `np.loadtxt`.
"""
function open_surface_track(path)
    io = open(path, "w")
    println(io, "# time_s x_m depth_m exx_1_s exy_1_s phase dt_s vertex_id vx_m_s vy_m_s")
    println(io, "# columns 1-8 as SZU2019 (exx, exy = strain-rate components); columns 9-10 added: velocity")
    println(io, "# marker-chain surface at fixed grid-vertex x (vertex_id = vertex index); depth = -h, positive down")
    println(io, "# exx, exy, vx, vy sampled half a cell below the surface (rock side); vx, vy in m/s, vy > 0 is up")
    println(io, "# phase = dominant non-air phase of the cell at that sample point; dt_s = step that moves this state")
    return io
end

# Bilinear interpolation of A on the rectilinear (possibly refined) grid xs × ys.
function bilinear_rect(A, xs, ys, x, y)
    i = clamp(searchsortedlast(xs, x), 1, length(xs) - 1)
    j = clamp(searchsortedlast(ys, y), 1, length(ys) - 1)
    tx = (x - xs[i]) / (xs[i + 1] - xs[i])
    ty = (y - ys[j]) / (ys[j + 1] - ys[j])
    return (1 - tx) * (1 - ty) * A[i, j] + tx * (1 - ty) * A[i + 1, j] +
           (1 - tx) * ty * A[i, j + 1] + tx * ty * A[i + 1, j + 1]
end

"""
    write_surface_track!(io, t, dt, chain, stokes, phase_ratios, grid, air_phase; xlim)

Append one row per marker-chain vertex with `xlim[1] ≤ x ≤ xlim[2]`: the surface height
`h(t)` with the strain rate, velocity and `dt` that move it during this step. Fields are
interpolated half a cell below the surface (strain rates from cell centres, velocities from
the staggered Vx/Vy grids), so they are rock values rather than a mix with the air across the
interface. Columns 1-8 follow the SZU2019 marker track; vx, vy are appended.
"""
function write_surface_track!(io, t, dt, chain, stokes, phase_ratios, grid, air_phase; xlim = (600.0e3, 1500.0e3))
    xv, yv = Array(grid.xvi[1]), Array(grid.xvi[2])
    nx, ny = length(xv) - 1, length(yv) - 1
    h = Array(chain.h_vertices)
    Vx, Vy = Array(stokes.V.Vx), Array(stokes.V.Vy)
    (xVx, yVx), (xVy, yVy) = map(g -> Array.(g), grid.xi_vel)
    @assert size(Vx) == (length(xVx), length(yVx)) && size(Vy) == (length(xVy), length(yVy))
    xc, yc = Array(grid.xci[1]), Array(grid.xci[2])
    exx, exy = Array(stokes.ε.xx), Array(stokes.ε.xy_c) # both at cell centres
    @assert size(exx) == size(exy) == (length(xc), length(yc))
    pc = Array(phase_ratios.center.data) # (cells, nphases, 1), cells in column-major (i, j) order

    for i in eachindex(xv)
        x = xv[i]
        xlim[1] ≤ x ≤ xlim[2] || continue
        # rock-side sample point: half the local cell height below the surface
        jh = clamp(searchsortedlast(yv, h[i]), 1, ny)
        y = h[i] - 0.5 * (yv[jh + 1] - yv[jh])
        vx = bilinear_rect(Vx, xVx, yVx, x, y)
        vy = bilinear_rect(Vy, xVy, yVy, x, y)
        εxx = bilinear_rect(exx, xc, yc, x, y)
        εxy = bilinear_rect(exy, xc, yc, x, y)
        # dominant non-air phase of the cell containing the sample point
        c = clamp(searchsortedlast(xv, x), 1, nx) + (clamp(searchsortedlast(yv, y), 1, ny) - 1) * nx
        phase, best = air_phase, 0.0
        for p in axes(pc, 2)
            p == air_phase && continue
            pc[c, p, 1] > best && ((phase, best) = (p, pc[c, p, 1]))
        end
        println(io, t, ' ', x, ' ', -h[i], ' ', εxx, ' ', εxy, ' ', phase, ' ', dt, ' ', i, ' ', vx, ' ', vy)
    end
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
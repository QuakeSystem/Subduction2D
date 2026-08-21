"""
NonuniformGrid.jl

Generic 1D logistic grid-refinement machinery: given a domain [x0, x1], build a
non-uniform set of vertices that is finer (down to a target cell size) inside
a window around `x_center` and coarser outside it, using a logistic (sigmoid)
stretching function.
"""

using Statistics

"""
    logistic_vertices(n_cells, L, x0; x_center, w_ref, refine_factor, k=4.0)

Build `n_cells` non-uniform cell widths over a 1D domain of length `L`
starting at `x0`, refined around `x_center`. Cells within roughly `w_ref` of
`x_center` are near their base (fine) width; cells far from it are stretched
by up to `refine_factor`. `k` controls how sharp the logistic transition
between the two regimes is (higher = sharper).
"""
function logistic_vertices(
        n_cells::Int,
        L::Float64,
        x0::Float64;
        x_center::Float64,
        w_ref::Float64,
        refine_factor::Float64,
        k::Float64 = 4.0,
    )
    # Compute cell widths with a logistic stretching weight, then integrate.
    widths_lin = fill(L / n_cells, n_cells)
    x_centers_u = [x0 + (i - 0.5) * (L / n_cells) for i in 1:n_cells]

    local_width_weight(x) = begin
        d = abs(x - x_center) / w_ref
        s = 1 / (1 + exp(-k * (d - 1.0)))
        1 + (refine_factor - 1) * s
    end

    widths = similar(widths_lin)
    @inbounds for i in 1:n_cells
        widths[i] = widths_lin[i] * local_width_weight(x_centers_u[i])
    end
    # Renormalize to preserve the total domain length exactly.
    widths .*= L / sum(widths)

    vertices = zeros(n_cells + 1)
    vertices[1] = x0
    @inbounds for i in 1:n_cells
        vertices[i + 1] = vertices[i] + widths[i]
    end

    return vertices, widths, x_centers_u
end

"""
    nonuniform_coords_1d(n_vertices, x0, x1; ref_grid, refine_factor, w_ref_ratio,
                          x_center_frac, k=4.0, verbose=0)

Build `n_vertices` non-uniform grid vertices over `[x0, x1]`, refined around
`x0 + x_center_frac * (x1 - x0)`. `w_ref_ratio` and `x_center_frac` are given
as fractions of the domain length so this scales naturally with domain size.

If `ref_grid == 0` or `refine_factor == 1.0`, returns a plain `LinRange` --
this preserves the type that `Geometry` (and downstream JustPIC advection,
which currently dispatches on it) expects for a uniform grid.

Set `verbose = 1` to print the resulting min/max/average cell size as sanity check
"""
function nonuniform_coords_1d(
        n_vertices::Int,
        x0::Float64,
        x1::Float64;
        ref_grid::Int,
        refine_factor::Float64,
        w_ref_ratio::Float64, # width of refined region as a fraction of total domain length
        x_center_frac::Float64,
        k::Float64 = 4.0, # logistic stretching parameter (higher k = sharper transition)
        verbose::Int = 0,
    )
    if ref_grid == 0 || refine_factor == 1.0
        # Preserve the original type returned by `Geometry` (LinRange/StepRangeLen),
        # which downstream JustPIC advection currently dispatches on.
        return LinRange(x0, x1, n_vertices)
    end

    n_cells = n_vertices - 1
    L = x1 - x0
    x_center = x0 + x_center_frac * L
    w_ref = w_ref_ratio * abs(L)

    vertices, widths, _ = logistic_vertices(
        n_cells,
        abs(L),
        x0;
        x_center = x_center,
        w_ref = w_ref,
        refine_factor = refine_factor,
        k = k,
    )

    if verbose == 1
        println("[grid verbose] non-uniform grid from $(x0) to $(x1):")
        println("  max cell size = $(maximum(widths))")
        println("  min cell size = $(minimum(widths))")
        println("  average cell size = $(mean(widths))")
        println("  refine_factor = $(refine_factor)")
        println("  w_ref_ratio = $(w_ref_ratio) -> w_ref = $(w_ref)")
        println("  center = $(x_center)")
    end

    return vertices
end
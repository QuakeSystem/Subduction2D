"""
NonuniformGrid.jl

1D grid-refinement machinery for the 2D subduction model, built on
PoissonGrids.jl's `solve_grid`: given a domain `[x0, x1]`, a fixed number of
cells, and a scalar monitor function `M(x)`, `solve_grid` iteratively
equidistributes vertices so that cells are narrow where `M` is large and wide
where `M` is small.

`poisson_coords_1d` is the only entry point used downstream -- it wraps
`solve_grid` with the uniform-grid fallback (`ref_grid == 0`) and verbose
reporting that `GMG_subduction_2D_with_coords` (in Subduction2D_setup.jl)
relies on. The monitor itself (e.g. `window_monitor`/`tanh_monitor` from
PoissonGrids) is built by the caller, since it differs between the x- and
z-axes.
"""

using PoissonGrids
using Statistics

"""
    poisson_coords_1d(n_vertices, x0, x1; ref_grid, monitor, tol=5e-5, maxiter=20_000, verbose=0, label="axis")

Build `n_vertices` vertices over `[x0, x1]` by equidistributing `monitor` over
`n_vertices - 1` cells with PoissonGrids' `solve_grid`.

If `ref_grid == 0`, returns a plain `LinRange` instead -- this preserves the
type that `Geometry` (and downstream JustPIC advection, which currently
dispatches on it) expects for a uniform grid.

`tol` defaults to `5e-5` rather than `solve_grid`'s own `1e-7` default: for a
sharp monitor (small `κ`, large contrast `α`) over ~900 cells, the residual
plateaus somewhere in the `1e-6`-`1e-5` range and stops improving regardless
of `maxiter`, without the resulting grid changing at all once it gets there
-- `5e-5` converges in under two seconds with comfortable margin above that
floor, while `1e-5` and tighter intermittently throw `ConvergenceError` for
some monitor/zone combinations even at `maxiter = 60_000`.

Set `verbose = 1` to print the resulting min/max/average cell size as a sanity
check.
"""
function poisson_coords_1d(
        n_vertices::Int,
        x0::Float64,
        x1::Float64;
        ref_grid::Int,
        monitor,
        tol::Float64 = 5.0e-5,
        maxiter::Int = 20_000,
        verbose::Int = 0,
        label::String = "axis",
    )
    if ref_grid == 0
        # Preserve the original type returned by `Geometry` (LinRange/StepRangeLen),
        # which downstream JustPIC advection currently dispatches on.
        return LinRange(x0, x1, n_vertices)
    end

    nc = n_vertices - 1
    vertices = solve_grid(x0, x1, monitor, nc; verbose = verbose == 1, tol = tol, maxiter = maxiter)

    if verbose == 1
        widths = diff(vertices)
        println("[grid verbose] non-uniform grid ($label) from $(x0) to $(x1):")
        println("  n cells = $nc")
        println("  max cell size = $(maximum(widths))")
        println("  min cell size = $(minimum(widths))")
        println("  average cell size = $(mean(widths))")
    end

    return vertices
end

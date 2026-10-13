#=
nonuniform_grid.jl

Nonuniform 1D vertex coordinates for the 2D subduction model, built with
PoissonGrids.jl. Specify the resolution you want and where; the cell count
follows from that.
=#
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using PoissonGrids
using Statistics
using CairoMakie
# --- optional check plot -----------------------------------------------------
function plot_cellsizes(x, z; file = joinpath(@__DIR__, "poisson_grid_cellsizes.png"))
    @eval using CairoMakie
    mid(v) = 0.5 .* (v[1:end-1] .+ v[2:end])
    fig = Figure(size = (900, 350))
    ax1 = Axis(fig[1, 1]; xlabel = "x [km]", ylabel = "Δx [km]", yscale = log10)
    ax2 = Axis(fig[1, 2]; xlabel = "z [km]", ylabel = "Δz [km]", yscale = log10)
    scatterlines!(ax1, mid(x), diff(x); markersize = 3)
    scatterlines!(ax2, mid(z), diff(z); markersize = 3)
    save(file, fig)
    return fig
end

"""
    refined_axis(xmin, xmax, M, zone, h_fine; nq, verbose, label)

Vertices on [xmin, xmax] equidistributing monitor `M`, with the number of cells
chosen so every cell inside `zone = (lo, hi)` is at most `h_fine` wide.
Solves M dx/dξ = const directly by inverting the cumulative integral of M.
"""
function refined_axis(xmin, xmax, M, zone, h_fine;
        nq = 200_001, verbose = true, label = "axis")
    xs = range(xmin, xmax; length = nq)
    dx = step(xs)
    Ms = M.(xs)

    # cumulative integral of M (trapezoid)
    F = zeros(nq)
    for i in 2:nq
        F[i] = F[i-1] + 0.5 * (Ms[i] + Ms[i-1]) * dx
    end
    ∫M = F[end]

    # h(x) = ∫M / (nc M(x)), so the widest cell in the zone is where M is smallest
    M_zone = minimum(M, range(zone[1], zone[2]; length = 2001))
    nc = ceil(Int, ∫M / (M_zone * h_fine))

    # vertex k sits where F = (k / nc) ∫M
    v = Vector{Float64}(undef, nc + 1)
    v[1], v[end] = xmin, xmax
    j = 1
    for k in 1:nc-1
        t = k / nc * ∫M
        while F[j+1] < t
            j += 1
        end
        v[k+1] = xs[j] + (t - F[j]) / (F[j+1] - F[j]) * dx
    end

    if verbose
        h = diff(v)
        mid = 0.5 .* (v[1:end-1] .+ v[2:end])
        inzone = zone[1] .<= mid .<= zone[2]
        ratio = maximum(max.(h[2:end] ./ h[1:end-1], h[1:end-1] ./ h[2:end]))
        println("[$label] $nc cells ($(nc + 1) vertices)")
        println("  finest / coarsest cell   = $(minimum(h)) / $(maximum(h))")
        println("  widest cell in $zone = $(maximum(h[inzone]))  (target $h_fine)")
        println("  max neighbour size ratio = $(round(ratio; digits = 4))")
    end
    return v
end

# --- resolution targets [km] -------------------------------------------------
h_fine   = 0.5                  # inside the refined zones
h_coarse = 5.0                  # background, away from the zones
α = h_coarse / h_fine - 1       # monitor contrast (cell width ∝ 1/M)

# Transition sharpness [1/km]. Larger κ = shorter transition but stronger
# cell-to-cell stretching; check the printed neighbour ratio (aim ≲ 1.05).
κ = 0.05
δ = 1.5 / κ   # push each tanh edge outward so the zone sits on the plateau (M ≥ ~0.95 max)

# --- x: refine 800–1100 km ---------------------------------------------------
x_zone = (800.0, 1100.0)
Mx = window_monitor(α, κ, (x_zone[2] - x_zone[1]) / 2 + δ, sum(x_zone) / 2)
x = refined_axis(0.0, 1500.0, Mx, x_zone, h_fine; label = "x")

# --- z: refine top 80 km -----------------------------------------------------
z_zone = (-80.0, 0.0)
Mz = tanh_monitor(α, κ, z_zone[1] - δ)   # 1 → 1+α for z > c, so the surface stays fine
z = refined_axis(-260.0, 0.0, Mz, z_zone, h_fine; label = "z")
plot_cellsizes(x, z)


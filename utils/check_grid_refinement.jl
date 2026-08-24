#=
check_grid_refinement.jl
=#
using Statistics
using CairoMakie

function logistic_vertices(
        n_cells::Int,
        L::Float64,
        x0::Float64;
        x_center::Float64,
        w_ref::Float64,
        refine_factor::Float64,
        k::Float64 = 4.0,
    )
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
    widths .*= L / sum(widths)
    vertices = zeros(n_cells + 1)
    vertices[1] = x0
    @inbounds for i in 1:n_cells
        vertices[i + 1] = vertices[i] + widths[i]
    end
    return vertices, widths, x_centers_u
end
function nonuniform_coords_1d(
        n_vertices::Int,
        x0::Float64,
        x1::Float64;
        ref_grid::Int,
        refine_factor::Float64,
        w_ref_ratio::Float64,
        x_center_frac::Float64,
        k::Float64 = 4.0,
        verbose::Int = 0,
    )
    if ref_grid == 0 || refine_factor == 1.0
        return LinRange(x0, x1, n_vertices)
    end
    n_cells = n_vertices - 1
    L = x1 - x0
    x_center = x0 + x_center_frac * L
    w_ref = w_ref_ratio * abs(L)
    vertices, widths, _ = logistic_vertices(
        n_cells, abs(L), x0;
        x_center = x_center, w_ref = w_ref, refine_factor = refine_factor, k = k,
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
function zone_stats(centers, widths_m, lo, hi)
    mask = (centers .>= lo) .& (centers .<= hi)
    zw = widths_m[mask]
    return (min = minimum(zw), max = maximum(zw), mean = mean(zw), n = length(zw))
end
# ------------------------------------------------------------
# Domain + target zoom window (edit these to match your setup)
# ------------------------------------------------------------
x0_km, x1_km = 0.0, 1500.0
z0_km, z1_km = -260.0, 0.0
xlo, xhi = 800.0, 1100.0   # target refined zone in x
zlo, zhi = -80.0, 0.0      # target refined zone in z
target_cell_m = 500.0
# ------------------------------------------------------------
# CURRENT parameters, exactly as in Subduction2D_setup.jl (unchanged):
# ------------------------------------------------------------
nx_vertices = 900 + 1
refine_factor_x = 20.0
w_ref_ratio_x = 660.0 / 1500.0
x_center_frac = 785.0 / 1500.0
k_x = 25.0
ny_vertices = 220 + 1
refine_factor_y = 30.0
w_ref_ratio_y = 170.0 / 260.0
y_center_frac = (-80.0 - (-260.0)) / 260.0
k_y = 20.0
ref_grid = 1
# ------------------------------------------------------------
# Build grid -- same call your real GMG_subduction_2D_with_coords makes
# ------------------------------------------------------------
x = nonuniform_coords_1d(nx_vertices, x0_km, x1_km;
    ref_grid = ref_grid, refine_factor = refine_factor_x, w_ref_ratio = w_ref_ratio_x,
    x_center_frac = x_center_frac, k = k_x, verbose = 1)
z = nonuniform_coords_1d(ny_vertices, z0_km, z1_km;
    ref_grid = ref_grid, refine_factor = refine_factor_y, w_ref_ratio = w_ref_ratio_y,
    x_center_frac = y_center_frac, k = k_y, verbose = 1)
# TRUE physical cell centers and widths -- exactly what xci / Δx actually are downstream
xc = 0.5 .* (x[1:(end - 1)] .+ x[2:end])
xw_m = diff(x) .* 1.0e3
zc = 0.5 .* (z[1:(end - 1)] .+ z[2:end])
zw_m = diff(z) .* 1.0e3
# ------------------------------------------------------------
# Report
# ------------------------------------------------------------
xz = zone_stats(xc, xw_m, xlo, xhi)
zz = zone_stats(zc, zw_m, zlo, zhi)
println("X target zone [$(xlo), $(xhi)] km -> cells: min=$(round(xz.min,digits=1)) m, max=$(round(xz.max,digits=1)) m, mean=$(round(xz.mean,digits=1)) m, n=$(xz.n)")
println("Z target zone [$(zlo), $(zhi)] km -> cells: min=$(round(zz.min,digits=1)) m, max=$(round(zz.max,digits=1)) m, mean=$(round(zz.mean,digits=1)) m, n=$(zz.n)")
println("Global X: min=$(round(minimum(xw_m),digits=1)) m, max=$(round(maximum(xw_m),digits=1)) m, total nx=$(nx_vertices-1)")
println("Global Z: min=$(round(minimum(zw_m),digits=1)) m, max=$(round(maximum(zw_m),digits=1)) m, total ny=$(ny_vertices-1)")
println("Logistic center (x) = $(x0_km + x_center_frac*(x1_km-x0_km)) km  (target zone midpoint = $((xlo+xhi)/2) km)")
println("Logistic center (z) = $(z0_km + y_center_frac*(z1_km-z0_km)) km  (target zone midpoint = $((zlo+zhi)/2) km)")
# ------------------------------------------------------------
# Plot -- against TRUE physical centers, not the uniform pre-stretch grid
# ------------------------------------------------------------
fig = Figure(size = (900, 800))
ax1 = Axis(fig[1, 1], xlabel = "x [km]", ylabel = "Δx [m]",
    title = "X cell size (nx_vertices=$nx_vertices, refine_factor=$refine_factor_x, w_ref_ratio=$(round(w_ref_ratio_x,digits=4)), k=$k_x)")
vspan!(ax1, xlo, xhi, color = (:orange, 0.15))
hlines!(ax1, [target_cell_m], color = :red, linestyle = :dash)
lines!(ax1, xc, xw_m, color = :steelblue, linewidth = 2)
ylims!(ax1, 0, maximum(xw_m) * 1.15)
ax2 = Axis(fig[2, 1], xlabel = "z [km]", ylabel = "Δz [m]",
    title = "Z cell size (ny_vertices=$ny_vertices, refine_factor=$refine_factor_y, w_ref_ratio=$(round(w_ref_ratio_y,digits=4)), k=$k_y)")
vspan!(ax2, zlo, zhi, color = (:orange, 0.15))
hlines!(ax2, [target_cell_m], color = :red, linestyle = :dash)
lines!(ax2, zc, zw_m, color = :seagreen, linewidth = 2)
ylims!(ax2, 0, maximum(zw_m) * 1.15)
save("grid_refinement_260km.png", fig)
println("Saved grid_refinement_260km.png")
 

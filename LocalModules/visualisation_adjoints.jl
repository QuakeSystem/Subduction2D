const CMAP_SENS = Reverse(:RdBu)      # matplotlib RdBu_r
const COLOR_OBS = :magenta
const COLOR_BOX = :black
const MS2CMYR   = 100 * 365.25 * 24 * 3600

## HELPERS ---------------------------------------------------------------------------

# symmetric colour range (-a, a), a the 99th percentile of |A| inside the zoom window
function symmetric_range_in_view(A, xmask, ymask; q = 0.99)
    v = sort!(filter(isfinite, abs.(vec(A[xmask, ymask]))))
    isempty(v) && return (-1.0, 1.0)
    a = v[clamp(ceil(Int, q * length(v)), 1, length(v))]
    return iszero(a) ? (-1.0, 1.0) : (-a, a)
end

# closed outline of a box from its centre and half widths
box_outline(cx, cz, hx, hz) = (
    [cx - hx, cx + hx, cx + hx, cx - hx, cx - hx],
    [cz - hz, cz - hz, cz + hz, cz + hz, cz - hz],
)

# observation box (magenta) and velocity boxes (black), all in km
function draw_boxes!(ax, observation, boxes)
    c, h = observation.center, observation.half_width
    lines!(ax, box_outline(c[1] / 1e3, c[2] / 1e3, h[1] / 1e3, h[2] / 1e3)...; color = COLOR_OBS, linewidth = 1.5)
    for b in boxes
        lines!(ax, box_outline(b.cenx / 1e3, b.cenz / 1e3, b.widthx / 2e3, b.widthz / 2e3)...; color = COLOR_BOX, linewidth = 1.5)
    end
    return ax
end

# observed field without its ghost nodes, with its node coordinates (m)
function observed_field(field, stokes, xvi, xci)
    if field === :Vx
        return Array(stokes.V.Vx)[:, 2:(end - 1)], Array(xvi[1]), Array(xci[2])
    elseif field === :Vy
        return Array(stokes.V.Vy)[2:(end - 1), :], Array(xci[1]), Array(xvi[2])
    elseif field === :P
        return Array(stokes.P), Array(xci[1]), Array(xci[2])
    end
    throw(ArgumentError("observation field must be :Vx, :Vy or :P"))
end

"""
    adjoint_objective(observation, stokes, xvi, xci) -> (; J, n, Jmean)

Objective of the box observation `(; field, center, half_width)` as `solve_DYREL!` defines it:
`J = Σ field` over the nodes within `half_width` of `center` (the same node selection as
JustRelax's `observation_mask`), the number of nodes `n`, and the mean `Jmean = J / n`.
"""
function adjoint_objective(observation, stokes, xvi, xci)
    (; field, center, half_width) = observation
    F, x, z = observed_field(field, stokes, xvi, xci)
    i = findall(xi -> abs(xi - center[1]) ≤ half_width[1], x)
    j = findall(zj -> abs(zj - center[2]) ≤ half_width[2], z)
    n = length(i) * length(j)
    J = sum(F[i, j])
    return (; J, n, Jmean = J / n)
end

# cell areas in km² (non-uniform grid)
cell_area_km2(xvi) = diff(Array(xvi[1]) .* 1e-3) .* diff(Array(xvi[2]) .* 1e-3)'

# air cells (ρ < 1000), as in make_figure
air_mask(ρg) = abs.(Array(ustrip.(ρg[2]))) ./ 9.81 .< 1000.0

"""
    adjoint_sensitivity_fields(xvi, ρg, stokes, stokes_ad, n) -> (; ρ_dJdρ, η_dJdη)

Relative sensitivities of the mean objective J̄ = J / n on the cell centres, per unit area
[cm/yr per km²], air masked: ρ·∂J̄/∂ρ and η·∂J̄/∂η (η the converged viscosity).
"""
function adjoint_sensitivity_fields(xvi, ρg, stokes, stokes_ad, n)
    area = cell_area_km2(xvi)
    air  = air_mask(ρg)
    mask_air(A) = ifelse.(air, NaN, A)
    ρ    = abs.(Array(ustrip.(ρg[2]))) ./ 9.81
    ρ_dJdρ = ρ .* Array(stokes_ad.ρ) .* (MS2CMYR / n) ./ area
    η_dJdη = Array(stokes.viscosity.η) .* Array(stokes_ad.viscosity.η) .* (MS2CMYR / n) ./ area
    return (; ρ_dJdρ = mask_air(ρ_dJdρ), η_dJdη = mask_air(η_dJdη))
end

# Value of the first field called `name` in a phase's material parameters (searched
# recursively through the laws), as Float64; NaN if the phase has no such parameter.
function material_parameter_value(x, name::Symbol, depth = 0)
    depth > 8 && return NaN
    x isa Union{Number, AbstractString, Symbol, AbstractArray, Module, Function} && return NaN
    T = typeof(x)
    (x isa Tuple || isstructtype(T)) || return NaN
    for f in fieldnames(T)
        isdefined(x, f) || continue
        v = getfield(x, f)
        if f === name
            val = hasfield(typeof(v), :val) ? v.val : v        # GeoUnit → value
            return val isa Real ? Float64(val) : NaN
        end
        r = material_parameter_value(v, name, depth + 1)
        isnan(r) || return r
    end
    return NaN
end

phase_label(m, p) = try
    "$p: " * first(unsafe_string(m.Name), 20)
catch
    "phase $p"
end

"""
    adjoint_parameter_fields(xvi, ρg, gradients, rheology, n; air_phase = 2)

For each parameter `p` in `gradients` (from `material_controls`), with J̄ = J / n:

  - `map`: Σ_phases p·∂J̄/∂p per unit area [cm/yr per km²] on the cell centres, air masked;
  - `phase_total`: p·∂J̄/∂p summed over the domain for each phase [cm/yr], i.e. the change of
    J̄ when that phase's value of `p` is raised by 100% everywhere;
  - `values`: the value of `p` in each phase (NaN where the phase does not use it);
  - `use`: the phases that use `p` (the air phase is left out).
"""
function adjoint_parameter_fields(xvi, ρg, gradients, rheology, n; air_phase = 2)
    area    = cell_area_km2(xvi)
    air     = air_mask(ρg)
    nphases = length(rheology)
    fields  = map(keys(gradients)) do name
        g      = Array(gradients[name].center)                   # (nphases, nx, ny), per cell
        values = [material_parameter_value(rheology[p], name) for p in 1:nphases]
        use    = [isfinite(values[p]) && p != air_phase for p in 1:nphases]
        rel    = zeros(size(g, 2), size(g, 3))
        phase_total = zeros(nphases)
        for p in 1:nphases
            use[p] || continue
            layer = values[p] .* g[p, :, :] .* (MS2CMYR / n)     # cm/yr per cell
            rel .+= layer
            phase_total[p] = sum(layer)
        end
        (; map = ifelse.(air, NaN, rel ./ area), phase_total, values, use)
    end
    return NamedTuple{keys(gradients)}(fields)
end

## FIGURES ---------------------------------------------------------------------------

"""
    make_adjoint_figure(it, t, xvi, xci, ρg, stokes, stokes_ad, observation, boxes,
                        xmin, xmax, ymin, ymax, figpath; version = nothing)

Observed field with the observation box (magenta) and velocity boxes (black), and the relative
sensitivities ρ·∂J̄/∂ρ and η·∂J̄/∂η per unit area. Limits in km.
"""
function make_adjoint_figure(
        it, t,
        xvi, xci,
        ρg,
        stokes, stokes_ad,
        observation, boxes,
        xmin, xmax, ymin, ymax,
        figpath; version = nothing
    )

    xc = Array(xci[1]) .* 1e-3
    yc = Array(xci[2]) .* 1e-3
    xi_mask = xmin .≤ xc .≤ xmax
    yi_mask = ymin .≤ yc .≤ ymax

    obj = adjoint_objective(observation, stokes, xvi, xci)
    S   = adjoint_sensitivity_fields(xvi, ρg, stokes, stokes_ad, obj.n)

    # observed field, in cm/yr for a velocity and MPa for the pressure
    field = observation.field
    F, xF, yF = observed_field(field, stokes, xvi, xci)
    F_scale, F_unit = field === :P ? (1.0e-6, "MPa") : (MS2CMYR, "cm/yr")
    F_plot  = F .* F_scale
    xF, yF  = xF .* 1e-3, yF .* 1e-3
    Jmean   = round(obj.Jmean * F_scale; sigdigits = 4)
    F_range = symmetric_range_in_view(F_plot, xmin .≤ xF .≤ xmax, ymin .≤ yF .≤ ymax; q = 1.0)

    t_Myr_round = round(t / (3600 * 24 * 365.25 * 1.0e6); digits = 3)

    fig = Figure(size = (900, 750), backgroundcolor = :white)
    fig[0, :] = GridLayout()
    Label(
        fig[0, 1:2],
        "Subduction2D - SZU2019 adjoint | version: $version | t = $t_Myr_round Myr | it = $it\n" *
            "⟨J⟩ = mean $field in box = $Jmean $F_unit ($(obj.n) nodes)",
        fontsize = 16, halign = :center, tellwidth = false
    )

    ar = 2.8
    tick_opts = (
        xticksmirrored = true, yticksmirrored = true,
        xtickalign = 1, ytickalign = 1,
        xgridvisible = false, ygridvisible = false,
    )
    ax1 = Axis(fig[1, 1], aspect = ar, title = "$field [$F_unit]"; tick_opts...)
    ax2 = Axis(fig[2, 1], aspect = ar, title = "ρ ∂⟨J⟩/∂ρ [cm/yr per km²]"; tick_opts...)
    ax3 = Axis(fig[3, 1], aspect = ar, title = "η ∂⟨J⟩/∂η [cm/yr per km²]"; tick_opts...)

    # ---- Plots ----
    h1 = heatmap!(ax1, xF, yF, F_plot; colormap = CMAP_SENS, colorrange = F_range)
    h2 = heatmap!(ax2, xc, yc, S.ρ_dJdρ; colormap = CMAP_SENS, colorrange = symmetric_range_in_view(S.ρ_dJdρ, xi_mask, yi_mask))
    h3 = heatmap!(ax3, xc, yc, S.η_dJdη; colormap = CMAP_SENS, colorrange = symmetric_range_in_view(S.η_dJdη, xi_mask, yi_mask))

    for ax in (ax1, ax2, ax3)
        draw_boxes!(ax, observation, boxes)
        xlims!(ax, xmin, xmax)
        ylims!(ax, ymin, ymax)
    end

    # Colorbars
    for (row, h) in [(1, h1), (2, h2), (3, h3)]
        Colorbar(fig[row, 2], h; height = 120, ticklabelsize = 10)
    end

    colsize!(fig.layout, 1, Aspect(1, ar))   # column as wide as the panels: colourbars next to them
    linkaxes!(ax1, ax2, ax3)
    resize_to_layout!(fig)
    mkpath(dirname(figpath))
    save(figpath, fig)
    return nothing
end

"""
    make_adjoint_parameter_figure(it, t, xvi, xci, ρg, stokes, gradients, rheology,
                                  observation, boxes, xmin, xmax, ymin, ymax, figpath;
                                  air_phase = 2, version = nothing)

One row per parameter in `gradients`: left the map of Σ_phases p·∂J̄/∂p per km², right a bar
per phase of p·∂J̄/∂p summed over the domain (the change of J̄ for +100% of that phase's `p`).
"""
function make_adjoint_parameter_figure(
        it, t,
        xvi, xci,
        ρg,
        stokes, gradients, rheology,
        observation, boxes,
        xmin, xmax, ymin, ymax,
        figpath; air_phase = 2, version = nothing
    )
    isempty(gradients) && return nothing

    xc = Array(xci[1]) .* 1e-3
    yc = Array(xci[2]) .* 1e-3
    xi_mask = xmin .≤ xc .≤ xmax
    yi_mask = ymin .≤ yc .≤ ymax

    obj     = adjoint_objective(observation, stokes, xvi, xci)
    P       = adjoint_parameter_fields(xvi, ρg, gradients, rheology, obj.n; air_phase)
    names   = keys(P)
    nphases = length(rheology)
    labels  = [phase_label(rheology[p], p) for p in 1:nphases]

    field = observation.field
    F_scale, F_unit = field === :P ? (1.0e-6, "MPa") : (MS2CMYR, "cm/yr")
    Jmean = round(obj.Jmean * F_scale; sigdigits = 4)
    t_Myr_round = round(t / (3600 * 24 * 365.25 * 1.0e6); digits = 3)

    fig = Figure(size = (1300, 230 * length(names) + 80), backgroundcolor = :white)
    fig[0, :] = GridLayout()
    Label(
        fig[0, 1:3],
        "Subduction2D - SZU2019 adjoint | version: $version | t = $t_Myr_round Myr | it = $it\n" *
            "⟨J⟩ = mean $field in box = $Jmean $F_unit | relative sensitivities p ∂⟨J⟩/∂p",
        fontsize = 16, halign = :center, tellwidth = false
    )

    ar = 2.8
    tick_opts = (
        xticksmirrored = true, yticksmirrored = true,
        xtickalign = 1, ytickalign = 1,
        xgridvisible = false, ygridvisible = false,
    )

    # ---- Plots ----
    axs = Axis[]
    for (row, name) in enumerate(names)
        S = P[name]
        ax = Axis(fig[row, 1], aspect = ar, title = "$name ∂⟨J⟩/∂$name [cm/yr per km²]"; tick_opts...)
        h = heatmap!(ax, xc, yc, S.map; colormap = CMAP_SENS, colorrange = symmetric_range_in_view(S.map, xi_mask, yi_mask))
        draw_boxes!(ax, observation, boxes)
        xlims!(ax, xmin, xmax)
        ylims!(ax, ymin, ymax)
        Colorbar(fig[row, 2], h; height = 120, ticklabelsize = 10)
        push!(axs, ax)

        # per-phase totals (whole domain, not only the view)
        phases = findall(S.use)
        axb = Axis(
            fig[row, 3]; title = "per phase, +100% $name [cm/yr]",
            yticks = (1:length(phases), labels[phases]), yreversed = true,
            yticklabelsize = 9, xgridvisible = false, ygridvisible = false,
        )
        if !isempty(phases)
            vals = S.phase_total[phases]
            barplot!(
                axb, 1:length(phases), vals; direction = :x,
                color = ifelse.(vals .>= 0, :firebrick, :steelblue),
            )
            vlines!(axb, 0; color = :black, linewidth = 1)
        end
    end
    colsize!(fig.layout, 1, Aspect(1, ar))   # column as wide as the panels: colourbars next to them
    colsize!(fig.layout, 3, Fixed(380))

    length(axs) > 1 && linkaxes!(axs...)
    resize_to_layout!(fig)
    mkpath(dirname(figpath))
    save(figpath, fig)
    return nothing
end

## OUTPUT ----------------------------------------------------------------------------

"""
    save_adjoint_sensitivities(path, it, t, xvi, xci, ρg, stokes, stokes_ad, observation,
                               gradients = (;); rheology = nothing, air_phase = 2)

Write the sensitivities of one adjoint solve to HDF5: coordinates in km, J and J̄ in the
unit of the observed field (cm/yr or MPa), the raw ∂J/∂ρ and ∂J/∂η per cell (SI, of the
sum J), the relative maps of J̄ per km² (cm/yr), and for each parameter `p` in `gradients`
a group `parameters/<p>` with the raw phase-resolved `dJdp` (nphases × nx × ny, SI, of J)
and, when `rheology` is given, `phase_values`, `rel_area` and `rel_phase_total`.
"""
function save_adjoint_sensitivities(
        path, it, t, xvi, xci, ρg, stokes, stokes_ad, observation, gradients = (;);
        rheology = nothing, air_phase = 2,
    )
    obj = adjoint_objective(observation, stokes, xvi, xci)
    S   = adjoint_sensitivity_fields(xvi, ρg, stokes, stokes_ad, obj.n)
    F_scale = observation.field === :P ? 1.0e-6 : MS2CMYR
    mkpath(dirname(path))
    HDF5.h5open(path, "w") do f
        f["it"]   = it
        f["t_s"]  = t
        f["J"]    = obj.J * F_scale
        f["Jmean"] = obj.Jmean * F_scale
        f["n_nodes"] = obj.n
        f["observation_field"]      = string(observation.field)
        f["observation_center"]     = collect(Float64, observation.center)
        f["observation_half_width"] = collect(Float64, observation.half_width)
        f["xc_km"] = Array(xci[1]) .* 1e-3
        f["zc_km"] = Array(xci[2]) .* 1e-3
        f["xv_km"] = Array(xvi[1]) .* 1e-3
        f["zv_km"] = Array(xvi[2]) .* 1e-3
        f["dJdrho"] = Array(stokes_ad.ρ)
        f["dJdeta"] = Array(stokes_ad.viscosity.η)
        f["rho_dJdrho_area"] = S.ρ_dJdρ
        f["eta_dJdeta_area"] = S.η_dJdη
        if !isempty(gradients)
            params = HDF5.create_group(f, "parameters")
            P = isnothing(rheology) ? nothing :
                adjoint_parameter_fields(xvi, ρg, gradients, rheology, obj.n; air_phase)
            for name in keys(gradients)
                g = HDF5.create_group(params, string(name))
                g["dJdp"] = Array(gradients[name].center)
                isnothing(P) && continue
                g["phase_values"]    = P[name].values
                g["rel_area"]        = P[name].map
                g["rel_phase_total"] = P[name].phase_total
            end
        end
    end
    return path
end

# In-run figures. Colour schemes follow the Cargese26 post-processing scripts
# (utils/plot_{lithology,stress,viscosity}_from_vtk.py):
#   lithology  - phase palette grouped by domain (mantle / oceanic / continental)
#   τII        - viridis, linear in MPa, 0 .. 99th percentile in view
#   η, η_vep   - reversed batlow (weak = light/pink, strong = dark blue), log10 18..25
#   velocities - red/blue diverging (RdBu_r), cm/yr
#   isotherms  - 10, 100, 150, 350, 450, 1300 °C, labelled; air (ρ < 1000) masked

# phase => (label, colour), as in plot_lithology_from_vtk.py / Subduction2D_setup.jl
const PHASE_STYLE = [
    1 => ("Mantle",              "#c9b8dd"),
    2 => ("Air",                 "#ffffff"),
    3 => ("Mantle weak zone",    "#7b5ea7"),
    4 => ("Ocean. crust/interf.", "#2a9d8f"),
    5 => ("Ocean. gabbro",       "#1d5c4f"),
    6 => ("Felsic crust",        "#e9a46a"),
    7 => ("Sediments",           "#f4d35e"),
]
const PHASE_CMAP   = cgrad([Makie.parse(Makie.Colorant, c) for (_, (_, c)) in PHASE_STYLE], length(PHASE_STYLE); categorical = true)
const ISOTHERMS_C  = [10, 100, 150, 350, 450, 1300]
const CMAP_DIVERGE = Reverse(:RdBu)     # matplotlib RdBu_r
const CMAP_TAUII   = :viridis
const CMAP_VISC    = Reverse(:batlow)

# 99th percentile of the finite values of A inside the zoom window (fallback 1.0)
function pct99_in_view(A, xmask, ymask)
    v = filter(isfinite, vec(A[xmask, ymask]))
    isempty(v) && return 1.0
    sort!(v)
    return max(v[clamp(ceil(Int, 0.99 * length(v)), 1, length(v))], eps())
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

    Tc = Array(T_buffer)
    xv = Array(xvi[1]) .* 1e-3
    yv = Array(xvi[2]) .* 1e-3
    xc = Array(xci[1]) .* 1e-3
    yc = Array(xci[2]) .* 1e-3

    xi_mask = xmin .≤ xc .≤ xmax
    yi_mask = ymin .≤ yc .≤ ymax
    xc_zoom = xc[xi_mask]
    yc_zoom = yc[yi_mask]
    Tc_zoom = Tc[xi_mask, yi_mask] .- 273.15   # °C

    ρ   = abs.(Array(ustrip.(ρg[2]))) ./ 9.81
    air = ρ .< 1000.0
    mask_air(A) = ifelse.(air, NaN, A)
    isograds_density = [2700, 2800, 2900, 3100, 3200, 3250, 3300, 3350]

    # cell-centre fields, air masked
    τII_MPa  = mask_air(Array(stokes.τ.II) ./ 1.0e6)
    logεII   = mask_air(log10.(Array(stokes.ε.II)))
    logη     = mask_air(log10.(Array(stokes.viscosity.η)))
    logη_vep = mask_air(log10.(Array(stokes.viscosity.η_vep)))
    τII_max  = pct99_in_view(τII_MPa, xi_mask, yi_mask)

    # velocities in cm/yr (input already clamped to ±10 cm/yr, in m/s)
    ms2cmyr = 100 * 365.25 * 24 * 3600
    Vx_cmyr = Array(Vx_limited) .* ms2cmyr
    Vy_cmyr = Array(Vy_limited) .* ms2cmyr
    vlim    = 10.0

    visc_range = (18.0, 25.0)
    ε_range    = (-15.0, -12.0)

    t_Myr_round = round(t / (3600 * 24 * 365.25 * 1.0e6); digits = 3)
    t_Kyr       = round(dt / (3600 * 24 * 365.25 * 1000); digits = 2)

    fig = Figure(size = (1500, 600), backgroundcolor = :white)
    fig[0, :] = GridLayout()
    Label(
        fig[0, 1:6],
        "Subduction2D - SZU2019 | version: $version | t = $t_Myr_round Myr | dt = $t_Kyr Kyr | it = $it",
        fontsize = 20, halign = :center, tellwidth = false
    )

    ar = 2.8
    tick_opts = (
        xticksmirrored = true, yticksmirrored = true,
        xtickalign = 1, ytickalign = 1,
        xgridvisible = false, ygridvisible = false,
    )
    ax1 = Axis(fig[1, 1], aspect = ar, title = "Lithology (markers)"; tick_opts...)
    ax2 = Axis(fig[2, 1], aspect = ar, title = "Temperature [°C]"; tick_opts...)
    ax3 = Axis(fig[3, 1], aspect = ar, title = "Density [kg/m³]"; tick_opts...)
    ax4 = Axis(fig[1, 3], aspect = ar, title = "Vx [cm/yr]"; tick_opts...)
    ax5 = Axis(fig[2, 3], aspect = ar, title = "Vy [cm/yr]"; tick_opts...)
    ax6 = Axis(fig[3, 3], aspect = ar, title = "τII [MPa]"; tick_opts...)
    ax7 = Axis(fig[1, 5], aspect = ar, title = "log10 εII [1/s]"; tick_opts...)
    ax8 = Axis(fig[2, 5], aspect = ar, title = "log10 η [Pa s]"; tick_opts...)
    ax9 = Axis(fig[3, 5], aspect = ar, title = "log10 η_vep [Pa s]"; tick_opts...)

    for ax in (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9)
        xlims!(ax, xmin, xmax)
        ylims!(ax, ymin, ymax)
    end

    # labelled isotherms; black on light maps, white on dark ones
    iso!(ax, color) = contour!(
        ax, xc_zoom, yc_zoom, Tc_zoom;
        levels = ISOTHERMS_C, color = color, linewidth = 0.8,
        labels = true, labelsize = 8, labelcolor = color,
        labelformatter = l -> "$(round(Int, l)) °C",
    )

    # ---- Plots ----
    # Lithology
    h1 = scatter!(
        ax1, pxv[idxv], pyv[idxv];
        color = clr[idxv], colormap = PHASE_CMAP,
        colorrange = (0.5, length(PHASE_STYLE) + 0.5), markersize = 1,
    )
    iso!(ax1, :black)

    # Temperature
    h2 = heatmap!(ax2, xc, yc, Tc .- 273.15; colormap = :lajolla, colorrange = (0, 1450))
    iso!(ax2, :black)

    # Density
    h3 = heatmap!(ax3, xc, yc, ρ; colormap = :vik, colorrange = (2500, 3300))
    contour!(ax3, xc, yc, ρ; levels = isograds_density, color = :white, linewidth = 0.8)

    # Velocities
    h4 = heatmap!(ax4, xv, yv, Vx_cmyr; colormap = CMAP_DIVERGE, colorrange = (-vlim, vlim))
    h5 = heatmap!(ax5, xv, yv, Vy_cmyr; colormap = CMAP_DIVERGE, colorrange = (-vlim, vlim))

    # Second invariant of the deviatoric stress
    h6 = heatmap!(ax6, xc, yc, τII_MPa; colormap = CMAP_TAUII, colorrange = (0.0, τII_max), highclip = :yellow)
    iso!(ax6, :black)

    # Second invariant of the strain rate
    h7 = heatmap!(ax7, xc, yc, logεII; colormap = :batlow, colorrange = ε_range)
    iso!(ax7, :white)

    # Viscosities
    h8 = heatmap!(ax8, xc, yc, logη; colormap = CMAP_VISC, colorrange = visc_range)
    iso!(ax8, :white)
    h9 = heatmap!(ax9, xc, yc, logη_vep; colormap = CMAP_VISC, colorrange = visc_range)
    iso!(ax9, :white)

    # Colorbars
    for (row, col, h) in [
            (2, 2, h2), (3, 2, h3),
            (1, 4, h4), (2, 4, h5), (3, 4, h6),
            (1, 6, h7), (2, 6, h8), (3, 6, h9),
        ]
        Colorbar(fig[row, col], h; height = 120, ticklabelsize = 10)
    end
    # lithology legend as a labelled categorical colourbar
    Colorbar(
        fig[1, 2], h1; height = 120, ticklabelsize = 8,
        ticks = (1:length(PHASE_STYLE), ["$p: $(lbl)" for (p, (lbl, _)) in PHASE_STYLE]),
    )

    linkaxes!(ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9)
    mkpath(dirname(figpath))
    save(figpath, fig)
    return nothing
end

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
    isotherms_C = [100, 150, 350, 450, 900, 1300]
    isotherms_K = isotherms_C .+ 273
    Tv = Array(T_buffer)          # (321, 73) 
    xv = Array(xvi[1]) .* 1e-3   # length 321 (nx+1)
    yv = Array(xvi[2]) .* 1e-3   # length 73  (ny+1)
    xc = Array(xci[1]) .* 1e-3   # length 320 (nx)
    yc = Array(xci[2]) .* 1e-3   # length 72  (ny)

    xi_mask = xmin .≤ xv .≤ xmax
    yi_mask = ymin .≤ yv .≤ ymax
    xv_zoom = xv[xi_mask]
    yv_zoom = yv[yi_mask]
    Tv_zoom = Tv[xi_mask, yi_mask]


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
        labels=false,
        labelsize = 12,
        color = :white,
        linewidth = 1.0
    )

    # Temperature
    h2 = heatmap!(
        ax2,
        xv,
        yv,
        Tv,
        colorrange = (253, 1718),
        colormap = :vik
    )

    contour!(
        ax2,
        xv_zoom,
        yv_zoom,
        Tv_zoom;
        levels = isotherms_K,
        labels=false,
        labelsize = 12,
        color = :white,
        linewidth = 1.0
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
        linewidth = 1.0
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
        labels=false,
        labelsize = 12,
        color = :white,
        linewidth = 1.0
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
        labels=false,
        labelsize = 12,
        color = :white,
        linewidth = 1.0

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
        labels=false,
        labelsize = 12,
        color = :white,
        linewidth = 1.0
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
    mkpath(dirname(figpath))

    save(figpath, fig)
    return nothing
end
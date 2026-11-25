#!/usr/bin/env python3
"""
plot_surface_spacetime.py

Space-time plot of surface vertical velocity (x vs time, colour = vy in mm/yr, up > 0), in
the style of the "Surface velocities" figure of stm-rsf-vis/2019_GPS_surface_displacement.py.
One vy panel per run, side by side.

  --jr-dhdt  JustRelax marker chain as topography_timeseries.csv (from extract_markerchain.py),
         or its run/vtk directory. For runs without a surface track. Plots the rate of surface
         height change at fixed x, dh/dt [mm/yr] (uplift > 0), computed as in
         plot_topography_comparison.py --quantity dhdt (centred difference over
         --dhdt-window-kyr) and interpolated to a 1 km x grid. dh/dt is not vy: they differ by
         vx * dh/dx, mostly near the trench. These panels come first (leftmost).

  --jr   JustRelax surface track (vtk/surface_track.txt written by Subduction2D.jl, or its
         run/vtk directory). Velocities are the solver's, sampled just below the surface at
         fixed grid x (column 10).
  --szu  SZU2019 marker track (eachdt_markertrack_lphr0.txt, or its run directory).
         Velocities from marker displacement between steps, as in the GPS script, keeping
         The track is first subsampled to steps at least --szu-min-dt-kyr apart (default 1 kyr,
         the geodynamic-phase step), so the shorter steps of the transition phase after ~3.7 Myr
         are included at the same sampling instead of being dropped.

Usage:
    python3 plot_surface_spacetime.py [--jr-dhdt <run>/vtk] --jr <run> [--szu <SZU run>] \\
        --xlim 800 1250 --vlim 5 [--titles "JustRelax" "SZU2019"] -o spacetime.png
"""

import argparse
import os
import textwrap

import numpy as np
import matplotlib.pyplot as plt

from plot_topography_comparison import read_columns, SECONDS_PER_YR, load_markerchain_csv, height_rate

try:
    from cmcrameri import cm
    CMAP = cm.vik
except ImportError:
    CMAP = "RdBu_r"
    print("cmcrameri not found, using RdBu_r instead of vik")

MS_TO_MMYR = 1.0e3 * SECONDS_PER_YR


def find(path, names):
    if os.path.isfile(path):
        return path
    for n in names:
        p = os.path.join(path, n)
        if os.path.isfile(p):
            return p
    raise SystemExit(f"no {names[0]} found for {path}")


def blocks(data):
    """Reshape rows (one block of equal size per timestep) to (nt, npts, ncols); drops an incomplete last block."""
    npts = np.argmax(data[:, 0] != data[0, 0])
    nt = len(data) // npts
    return data[:nt * npts].reshape(nt, npts, data.shape[1])


def load_jr(path):
    d = blocks(read_columns(path, (0, 1, 9)))
    return d[0, :, 1] / 1e3, d[:, 0, 0] / SECONDS_PER_YR, d[:, :, 2] * MS_TO_MMYR


def load_szu(path, min_dt_kyr=1.0):
    d = blocks(read_columns(path, (0, 1, 2)))
    t, x, depth = d[:, 0, 0], d[:, :, 1], d[:, :, 2]
    # keep steps at least ~min_dt apart (small tolerance for rounding of the 1 kyr steps)
    min_dt = 0.99 * min_dt_kyr * 1e3 * SECONDS_PER_YR
    keep, last = [0], t[0]
    for i in range(1, len(t)):
        if t[i] - last >= min_dt:
            keep.append(i)
            last = t[i]
    t, x, depth = t[keep], x[keep], depth[keep]
    dt = np.diff(t)[:, None]
    vy = -np.diff(depth, axis=0) / dt * MS_TO_MMYR  # depth positive down -> vy positive up
    return x[:-1] / 1e3, t[:-1] / SECONDS_PER_YR, vy


def load_jr_dhdt(path, window_kyr, dx_km=1.0):
    rates = height_rate(load_markerchain_csv(path), window_kyr * 1e3 * SECONDS_PER_YR)
    if not rates:
        raise SystemExit(f"{path}: run shorter than the {window_kyr} kyr dh/dt window")
    t = np.array(sorted(rates))
    x0, x1 = (min(np.nanmin(rates[k][0]) for k in t), max(np.nanmax(rates[k][0]) for k in t))
    x = np.arange(np.ceil(x0), np.floor(x1) + dx_km / 2, dx_km)
    r = np.array([np.interp(x, *rates[k], left=np.nan, right=np.nan) for k in t])
    return x, t / SECONDS_PER_YR, r


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--jr", nargs="*", default=[], help="JustRelax surface_track.txt or run/vtk dir")
    ap.add_argument("--jr-dhdt", nargs="*", default=[],
                    help="JustRelax marker-chain topography_timeseries.csv or run/vtk dir (plots dh/dt)")
    ap.add_argument("--dhdt-window-kyr", type=float, default=20.0,
                    help="time window of the dh/dt centred difference [kyr] (default 20)")
    ap.add_argument("--szu", nargs="*", default=[], help="SZU marker track or run dir")
    ap.add_argument("--xlim", type=float, nargs=2, default=(800, 1250), metavar=("XMIN", "XMAX"))
    ap.add_argument("--vlim", type=float, default=5.0, help="vy colour limit [mm/yr]")
    ap.add_argument("--titles", nargs="*",
                    help="panel titles, one per run in plotting order (all --jr-dhdt, then --jr, then --szu runs); "
                         "'\\n' in a title starts a new line; default: run directory names")
    ap.add_argument("--szu-min-dt-kyr", type=float, default=1.0,
                    help="SZU velocity sampling interval [kyr] (default 1, the geodynamic step)")
    ap.add_argument("--time-unit", choices=("yr", "kyr", "Myr"), default="yr", help="time axis unit (default yr)")
    ap.add_argument("--tmax-myr", type=float, help="fixed upper time limit [Myr] (default: end of the longest run)")
    ap.add_argument("--fontsize", type=float, default=16, help="base font size: axes, ticks, colour bars (default 16)")
    ap.add_argument("--title-fontsize", type=float, help="panel title size (default: --fontsize); figure title is 4 pt larger")
    ap.add_argument("-o", "--out", default="surface_spacetime.png")
    args = ap.parse_args()
    plt.rcParams.update({"font.size": args.fontsize})
    title_fs = args.title_fontsize or args.fontsize
    tscale = {"yr": 1.0, "kyr": 1e3, "Myr": 1e6}[args.time_unit]  # loaders return time in years

    csv_names = ("topography_timeseries.csv", "topography_csv/topography_timeseries.csv",
                 "vtk/topography_csv/topography_timeseries.csv")
    runs = [(os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(f))))),
             load_jr_dhdt(f, args.dhdt_window_kyr), "dh/dt")
            for f in (find(p, csv_names) for p in args.jr_dhdt)]
    runs += [(os.path.basename(os.path.dirname(os.path.dirname(os.path.abspath(f)))), load_jr(f), "vy")
            for f in (find(p, ("surface_track.txt", "vtk/surface_track.txt")) for p in args.jr)]
    runs += [(os.path.basename(os.path.dirname(os.path.abspath(f))), load_szu(f, args.szu_min_dt_kyr), "vy")
             for f in (find(p, ("eachdt_markertrack_lphr0.txt",)) for p in args.szu)]
    if not runs:
        raise SystemExit("pass --jr-dhdt, --jr and/or --szu")
    if args.titles:
        if len(args.titles) != len(runs):
            raise SystemExit(f"--titles: got {len(args.titles)} titles for {len(runs)} panels")
        runs = [(title, data, q) for title, (_, data, q) in zip(args.titles, runs)]

    # common time axis for all panels; a run that ends earlier stays white above its last step
    fig, axs = plt.subplots(1, len(runs), figsize=(7 * len(runs), 14), squeeze=False, sharey=True)
    for (name, (x, t, vy), quantity), ax in zip(runs, axs[0]):
        X = np.broadcast_to(x, vy.shape) if x.ndim == 1 else x
        T = np.broadcast_to(t[:, None] / tscale, vy.shape)
        pc = ax.pcolormesh(X, T, vy, cmap=CMAP, vmin=-args.vlim, vmax=args.vlim, shading="nearest")
        fig.colorbar(pc, ax=ax).set_label("$v_y$ (mm/yr)" if quantity == "vy" else "$dh/dt$ (mm/yr)")
        ax.set_xlim(args.xlim)
        ax.set_xlabel("X (km)")
        # explicit '\n' line breaks are kept; only lines that are still too long are wrapped
        lines = name.replace("\\n", "\n").split("\n")
        ax.set_title("\n".join(textwrap.fill(l, 34) for l in lines), fontsize=title_fs)
    axs[0, 0].set_ylabel(f"Time ({args.time_unit})")
    tmax = args.tmax_myr * 1e6 if args.tmax_myr else max(t.max() for _, (_, t, _), _ in runs)
    axs[0, 0].set_ylim(min(t.min() for _, (_, t, _), _ in runs) / tscale, tmax / tscale)
    fig.suptitle("Surface vertical velocity" if all(q == "vy" for *_, q in runs)
                 else "Surface uplift and subsidence", fontsize=title_fs + 4, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.98))  # room for the suptitle
    fig.savefig(args.out, dpi=200, bbox_inches="tight")
    print(f"Saved {args.out}")


if __name__ == "__main__":
    main()

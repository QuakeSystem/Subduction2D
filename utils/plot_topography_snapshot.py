#!/usr/bin/env python3
"""
plot_topography_snapshot.py

Overlay the surface topography of several runs at one moment in a single panel.

  --szu      SZU2019 run directories. The profile is the marker-track step nearest to
             --time-myr, or, with --szu-step, nearest to the time of that output step
             (ModelGroup/Model[0] in lphr<step>.gzip.h5, needs h5py).
  --jr       JustRelax inputs as for plot_topography_comparison.py (run/vtk dir,
             topography_timeseries.csv or surface_track.txt); profile nearest --time-myr.
             Density-surface CSVs are drawn dashed.

Usage:
    python3 plot_topography_snapshot.py --szu /scratch/tectonics/bert/SZU2019-gps* \\
        --jr <run>/vtk <run>/vtk/topography_density_csv/topography_timeseries.csv \\
        --time-myr 4 --szu-step 217 -o topography_4Myr.png   # or just --time-myr 3.5
"""

import argparse
import os

import numpy as np
import matplotlib.pyplot as plt

from plot_topography_comparison import (
    SECONDS_PER_MYR, JR_CSV_CANDIDATES, JR_TRACK_CANDIDATES,
    resolve_jr, resolve_szu, run_label, load_markertrack, load_markerchain_csv,
)


def szu_time_myr(run_dir, step, default):
    """Model time [Myr] of SZU output step `step`, or `default` if it cannot be read."""
    path = os.path.join(run_dir, f"lphr{step}.gzip.h5")
    try:
        import h5py
        with h5py.File(path, "r") as f:
            return f["ModelGroup/Model"][0] / 1.0e6  # stored in years
    except (ImportError, OSError, KeyError):
        print(f"  could not read time of step {step} from {path}, using {default} Myr")
        return default


def nearest(profiles, t_myr):
    times = np.array(sorted(profiles))
    t = times[np.argmin(np.abs(times - t_myr * SECONDS_PER_MYR))]
    return t / SECONDS_PER_MYR, profiles[t]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--szu", nargs="*", default=[])
    ap.add_argument("--szu-file", default="eachdt_markertrack_lphr0.txt")
    ap.add_argument("--szu-step", type=int, help="match SZU runs to this output step (e.g. 217 = ~4 Myr)")
    ap.add_argument("--jr", nargs="*", default=[])
    ap.add_argument("--time-myr", type=float, default=4.0)
    ap.add_argument("--xlim", type=float, nargs=2, default=(700, 1300), metavar=("XMIN", "XMAX"))
    ap.add_argument("--ylim", type=float, nargs=2, metavar=("YMIN", "YMAX"))
    ap.add_argument("--labels", nargs="*", help="legend names, in plotting order (all --szu runs, then all --jr runs)")
    ap.add_argument("--fontsize", type=float, default=10)
    ap.add_argument("--legend-loc", default="outside",
                    help="'outside' (right of the plot, default) or a matplotlib location inside it, e.g. 'lower right'")
    ap.add_argument("-o", "--out", default="topography_snapshot.png")
    args = ap.parse_args()
    plt.rcParams.update({"font.size": args.fontsize})
    names = iter(args.labels or [])

    fig, ax = plt.subplots(figsize=(16, 5))
    colors = plt.get_cmap("tab10").colors
    n = 0
    for p in args.szu:
        f = resolve_szu(p, args.szu_file)
        if f is None:
            print(f"skipping {p}: no {args.szu_file}")
            continue
        if args.szu_step is None:
            target, tag = args.time_myr, ""
        else:
            target, tag = szu_time_myr(os.path.dirname(f), args.szu_step, args.time_myr), f"step {args.szu_step}, "
        t, (x, h) = nearest(load_markertrack(f), target)
        ax.plot(x, h, color=colors[n % 10], lw=1.5, label=f"{next(names, run_label(f, p))} ({tag}{t:.3f} Myr)")
        n += 1
    for p in args.jr:
        f = resolve_jr(p, JR_CSV_CANDIDATES + JR_TRACK_CANDIDATES)
        if f is None:
            print(f"skipping {p}: no topography CSV or surface track")
            continue
        loader = load_markertrack if f.endswith(".txt") else load_markerchain_csv
        t, (x, h) = nearest(loader(f), args.time_myr)
        density = "topography_density_csv" in f
        ax.plot(x, h, color="k", lw=1.0 if density else 2.0, ls="--" if density else "-",
                label=f"{next(names, run_label(f, p))} ({'density surface' if density else 'marker chain'}, {t:.3f} Myr)")

    ax.set_xlim(args.xlim)
    if args.ylim:
        ax.set_ylim(args.ylim)
    ax.set_xlabel("x [km]")
    ax.set_ylabel("topography height [km]")
    ax.set_title(f"Surface topography at ~{args.time_myr:g} Myr")
    if args.legend_loc == "outside":  # keeps the profiles clear
        ax.legend(fontsize=0.8 * args.fontsize, loc="upper left", bbox_to_anchor=(1.01, 1))
    else:
        ax.legend(fontsize=0.7 * args.fontsize, loc=args.legend_loc, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(args.out, dpi=200, bbox_inches="tight")
    print(f"Saved {args.out}")


if __name__ == "__main__":
    main()

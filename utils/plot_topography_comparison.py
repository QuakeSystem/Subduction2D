#!/usr/bin/env python3
"""
plot_topography_comparison.py

Compare evolving surface topography between any number of model runs, one
panel per run:
  1. SZU2019-gps reference-code runs (--szu): marker-track output
     (eachdt_markertrack_lphr0.txt), fixed-identity markers tracked through
     time. Columns (no header, space-delimited):
       time_s, x_m, depth_m, vx_m_s, vy_m_s, phase, dt_s, marker_id
     `depth_m` is positive-down; converted to height = -depth_m here to match
     the usual up-is-positive topography convention.
  2. This project's JustRelax runs (--jr): marker-chain output, already
     extracted by extract_markerchain.py into topography_timeseries.csv
     (time_s, time_kyr, x_km, h_km) -- one full x-profile per saved timestep.

Each argument may be a file or a run directory:
  --szu DIR  -> DIR/<--szu-file>  (default eachdt_markertrack_lphr0.txt)
  --jr  DIR  -> DIR/topography_csv/topography_timeseries.csv, or
                DIR/vtk/topography_csv/topography_timeseries.csv
Missing inputs are skipped with a warning, so a shell glob over all runs works.

For each run, plots a subsample of topography profiles (x vs height) as
overlaid lines, colored by time. All panels share x/y limits.

Usage:
    python3 plot_topography_comparison.py \\
        --szu /scratch/tectonics/bert/SZU2019-gps* \\
        --jr  Subduction2D_SZU2019/data/<run>/vtk \\
        -o topography_comparison.png
"""

import argparse
import csv
import math
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

SECONDS_PER_YR = 3600 * 24 * 365.25
SECONDS_PER_KYR = SECONDS_PER_YR * 1000
SECONDS_PER_MYR = SECONDS_PER_YR * 1.0e6

JR_CSV_CANDIDATES = (
    "topography_timeseries.csv",
    os.path.join("topography_csv", "topography_timeseries.csv"),
    os.path.join("vtk", "topography_csv", "topography_timeseries.csv"),
)


def resolve_szu(path, fname):
    if os.path.isdir(path):
        path = os.path.join(path, fname)
    return path if os.path.isfile(path) else None


def resolve_jr(path):
    if os.path.isfile(path):
        return path
    for cand in JR_CSV_CANDIDATES:
        p = os.path.join(path, cand)
        if os.path.isfile(p):
            return p
    return None


def run_label(path, input_path):
    """Short panel label: the run directory name (skipping vtk/topography_csv levels)."""
    d = input_path if os.path.isdir(input_path) else os.path.dirname(path)
    d = os.path.normpath(os.path.abspath(d))
    while os.path.basename(d) in ("vtk", "topography_csv"):
        d = os.path.dirname(d)
    return os.path.basename(d)


def load_markertrack(path):
    """Return dict: time_s -> (x_km, height_km), sorted by x, one entry per saved timestep."""
    data = np.loadtxt(path, usecols=(0, 1, 2))
    t, x_m, depth_m = data[:, 0], data[:, 1], data[:, 2]
    height_km = -depth_m / 1.0e3
    x_km = x_m / 1.0e3

    # Group by time with one sort instead of a mask per timestep (~10^4 steps).
    order = np.lexsort((x_km, t))
    t, x_km, height_km = t[order], x_km[order], height_km[order]
    tu, start = np.unique(t, return_index=True)
    bounds = np.append(start, len(t))
    return {tu[i]: (x_km[bounds[i]:bounds[i + 1]], height_km[bounds[i]:bounds[i + 1]])
            for i in range(len(tu))}


def load_markerchain_csv(path):
    """Return dict: time_s -> (x_km, h_km), sorted by x, one entry per saved timestep."""
    by_time = {}
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            t = float(row["time_s"])
            by_time.setdefault(t, []).append((float(row["x_km"]), float(row["h_km"])))

    profiles = {}
    for t, pts in by_time.items():
        pts.sort()
        xs, hs = zip(*pts)
        profiles[t] = (np.array(xs), np.array(hs))
    return profiles


def plot_profiles(ax, profiles, n_lines, cmap_name):
    times = np.array(sorted(profiles))
    idx = np.unique(np.linspace(0, len(times) - 1, min(n_lines, len(times))).astype(int))
    sel_times = times[idx]

    if sel_times[-1] >= SECONDS_PER_MYR:
        unit_s, unit_label = SECONDS_PER_MYR, "Myr"
    else:
        unit_s, unit_label = SECONDS_PER_KYR, "kyr"

    cmap = plt.get_cmap(cmap_name)
    norm = Normalize(vmin=sel_times[0] / unit_s, vmax=sel_times[-1] / unit_s)

    for t in sel_times:
        x, h = profiles[t]
        ax.plot(x, h, color=cmap(norm(t / unit_s)), lw=1.0)

    ax.set_xlabel("x [km]")
    ax.set_ylabel("topography height [km]")

    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cb = plt.colorbar(sm, ax=ax)
    cb.set_label(f"time [{unit_label}]")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--szu", nargs="*", default=[],
                        help="SZU2019 run directories or marker-track .txt files")
    parser.add_argument("--szu-file", default="eachdt_markertrack_lphr0.txt",
                        help="Marker-track filename inside each --szu directory")
    parser.add_argument("--jr", nargs="*", default=[],
                        help="JustRelax run/vtk directories or topography_timeseries.csv files")
    parser.add_argument("-o", "--out", default="topography_comparison.png")
    parser.add_argument("-n", "--n-lines", type=int, default=20, help="Number of profile lines per panel")
    parser.add_argument("--ncols", type=int, default=3, help="Panels per row")
    parser.add_argument("--cmap", default="viridis")
    args = parser.parse_args()

    # (label, kind, loader, path) for every input that actually exists.
    runs = []
    for p in args.szu:
        f = resolve_szu(p, args.szu_file)
        if f is None:
            print(f"WARNING: no {args.szu_file} for {p}, skipping")
            continue
        runs.append((run_label(f, p), "SZU2019 marker track", load_markertrack, f))
    for p in args.jr:
        f = resolve_jr(p)
        if f is None:
            print(f"WARNING: no topography_timeseries.csv for {p} "
                  f"(run extract_markerchain.py on its vtk/ dir first), skipping")
            continue
        runs.append((run_label(f, p), "JustRelax marker chain", load_markerchain_csv, f))

    if not runs:
        sys.exit("No inputs found; pass --szu and/or --jr")

    loaded = []
    for label, kind, loader, f in runs:
        print(f"Loading {f} ...")
        profiles = loader(f)
        print(f"  {len(profiles)} timesteps")
        loaded.append((label, kind, profiles))

    ncols = min(args.ncols, len(loaded))
    nrows = math.ceil(len(loaded) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.5 * ncols, 4.5 * nrows),
                             sharex=True, sharey=True, squeeze=False)
    axes = axes.ravel()

    for ax, (label, kind, profiles) in zip(axes, loaded):
        plot_profiles(ax, profiles, args.n_lines, args.cmap)
        ax.set_title(f"{label}\n({kind})", fontsize=10)
        ax.tick_params(labelleft=True, labelbottom=True)
    for ax in axes[len(loaded):]:
        ax.set_visible(False)

    fig.suptitle("Topography profile evolution (early = dark, late = yellow)")
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"Saved {args.out}")


if __name__ == "__main__":
    main()

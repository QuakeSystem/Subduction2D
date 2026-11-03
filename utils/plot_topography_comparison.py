#!/usr/bin/env python3
"""
plot_topography_comparison.py

Compare evolving surface topography between any number of model runs, one
panel per run:
  1. SZU2019-gps reference-code runs (--szu): marker-track output
     (eachdt_markertrack_lphr0.txt), fixed-identity markers tracked through
     time. Columns (no header, space-delimited), as written by loadevps_omp.c:
       time_s, x_m, depth_m, exx_1_s, exy_1_s, phase, dt_s, marker_id
     (columns 4-5 are marker strain-rate components, not velocities).
     `depth_m` is positive-down; converted to height = -depth_m here to match
     the usual up-is-positive topography convention.
  2. This project's JustRelax runs (--jr), either
     - marker-chain output extracted by extract_markerchain.py (or the density
       surface from extract_surface_from_vtk.py) into topography_timeseries.csv
       (time_s, time_kyr, x_km, h_km) -- one full x-profile per saved timestep, or
     - the per-step surface track vtk/surface_track.txt written by Subduction2D.jl:
       the SZU column layout above (same meaning for columns 1-8, marker_id =
       grid-vertex index) plus vx_m_s, vy_m_s (vy > 0 up) in columns 9-10.

Each argument may be a file or a run directory:
  --szu DIR  -> DIR/<--szu-file>  (default eachdt_markertrack_lphr0.txt)
  --jr  DIR  -> DIR/topography_csv/topography_timeseries.csv, or
                DIR/vtk/topography_csv/topography_timeseries.csv
                or, if neither exists (and always for --quantity vx/vy),
                DIR/surface_track.txt or DIR/vtk/surface_track.txt
Missing inputs are skipped with a warning, so a shell glob over all runs works.

For each run, plots a subsample of topography profiles (x vs height) as
overlaid lines, colored by time. All panels share x/y limits.

With --quantity dhdt, plots the rate of surface height change at fixed x,
dh/dt [mm/yr] (uplift > 0), instead of the height. It is a centred difference
over a time window (--window-kyr, default 100 kyr): per-step differences of the
chain are dominated by noise. Note this is the Eulerian rate dh/dt, not the
material vertical velocity vy (they differ by vx * dh/dx). It is computed the
same way for SZU and JustRelax runs, so the panels compare like with like.

With --quantity vx or vy, plots the material surface velocity [mm/yr] from the
JustRelax surface track (columns 9-10). SZU marker tracks have no velocity
columns, so SZU inputs are skipped with a warning for these quantities.

Usage:
    python3 plot_topography_comparison.py \\
        --szu /scratch/tectonics/bert/SZU2019-gps* \\
        --jr  Subduction2D_SZU2019/data/<run>/vtk \\
        -o topography_comparison.png
"""

import argparse
import csv
import io
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
JR_TRACK_CANDIDATES = ("surface_track.txt", os.path.join("vtk", "surface_track.txt"))

# 0-based velocity columns of the surface track (SZU layout + vx, vy)
VELOCITY_COLUMNS = {"vx": 8, "vy": 9}
MS_TO_MMYR = 1.0e3 * SECONDS_PER_YR


def resolve_szu(path, fname):
    if os.path.isdir(path):
        path = os.path.join(path, fname)
    return path if os.path.isfile(path) else None


def resolve_jr(path, candidates=JR_CSV_CANDIDATES):
    if os.path.isfile(path):
        return path
    for cand in candidates:
        p = os.path.join(path, cand)
        if os.path.isfile(p):
            return p
    return None


def run_label(path, input_path):
    """Short panel label: the run directory name (skipping vtk/topography_csv levels)."""
    d = input_path if os.path.isdir(input_path) else os.path.dirname(path)
    d = os.path.normpath(os.path.abspath(d))
    while os.path.basename(d) == "vtk" or os.path.basename(d).startswith("topography_"):
        d = os.path.dirname(d)
    return os.path.basename(d)


def count_columns(path):
    """Number of columns in the first data (non-'#') line of a whitespace-delimited file."""
    with open(path) as f:
        for line in f:
            if line.strip() and not line.lstrip().startswith("#"):
                return len(line.split())
    return 0


def read_columns(path, usecols):
    """np.loadtxt that tolerates a file still being written (partial last line) and rows
    with too few fields (glitched lines), which are dropped with a note."""
    with open(path) as f:
        text = f.read()
    if not text.endswith("\n"):
        text = text[:text.rfind("\n") + 1]  # drop the partially written last line
    try:
        return np.loadtxt(io.StringIO(text), usecols=usecols)
    except ValueError:
        need = max(usecols) + 1
        lines = [l for l in text.splitlines()
                 if l.strip() and not l.lstrip().startswith("#") and len(l.split()) >= need]
        print(f"  note: dropped malformed rows from {path}")
        return np.loadtxt(io.StringIO("\n".join(lines)), usecols=usecols)


def load_markertrack(path, column=None, scale=1.0):
    """Return dict: time_s -> (x_km, value), sorted by x, one entry per saved timestep.

    value is the height [km] (from depth_m) when column is None, else column `column`
    (0-based) times `scale`.
    """
    if column is None:
        data = read_columns(path, (0, 1, 2))
        value = -data[:, 2] / 1.0e3
    else:
        data = read_columns(path, (0, 1, column))
        value = data[:, 2] * scale
    t, x_km = data[:, 0], data[:, 1] / 1.0e3

    # Group by time with one sort instead of a mask per timestep (~10^4 steps).
    order = np.lexsort((x_km, t))
    t, x_km, value = t[order], x_km[order], value[order]
    tu, start = np.unique(t, return_index=True)
    bounds = np.append(start, len(t))
    return {tu[i]: (x_km[bounds[i]:bounds[i + 1]], value[bounds[i]:bounds[i + 1]])
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


def height_rate(profiles, window_s):
    """Return dict: t_centre -> (x_km, dh/dt [mm/yr]) from profiles differenced over ~window_s."""
    times = np.array(sorted(profiles))
    rates = {}
    for t in times:
        t0, t1 = t - window_s / 2, t + window_s / 2
        if t0 < times[0] or t1 > times[-1]:
            continue
        i0 = times[np.argmin(np.abs(times - t0))]
        i1 = times[np.argmin(np.abs(times - t1))]
        if i1 <= i0:
            continue
        x1, h1 = profiles[i1]
        x0, h0 = profiles[i0]
        h0_on_x1 = np.interp(x1, x0, h0, left=np.nan, right=np.nan)
        # km / s -> mm / yr
        rates[0.5 * (i0 + i1)] = (x1, (h1 - h0_on_x1) / (i1 - i0) * 1.0e6 * SECONDS_PER_YR)
    return rates


def plot_profiles(ax, profiles, n_lines, cmap_name, ylabel="topography height [km]"):
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
    ax.set_ylabel(ylabel)

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
    parser.add_argument("--xlim", type=float, nargs=2, metavar=("XMIN", "XMAX"), help="x range [km]")
    parser.add_argument("--ylim", type=float, nargs=2, metavar=("YMIN", "YMAX"),
                        help="y range [km for h, mm/yr otherwise]")
    parser.add_argument("--quantity", choices=("h", "dhdt", "vx", "vy"), default="h",
                        help="h: surface height (default); dhdt: rate of height change [mm/yr]; "
                             "vx, vy: surface velocity from the JustRelax surface track [mm/yr]")
    parser.add_argument("--window-kyr", type=float, default=100.0,
                        help="time window for --quantity dhdt [kyr] (default 100)")
    args = parser.parse_args()

    # (label, kind, loader, path) for every input that actually exists.
    velocity = args.quantity in VELOCITY_COLUMNS
    runs = []
    for p in args.szu:
        f = resolve_szu(p, args.szu_file)
        if f is None:
            print(f"WARNING: no {args.szu_file} for {p}, skipping")
            continue
        if velocity and count_columns(f) <= VELOCITY_COLUMNS[args.quantity]:
            print(f"WARNING: {f} has no {args.quantity} column (SZU marker tracks store strain rates), skipping")
            continue
        runs.append((run_label(f, p), "SZU2019 marker track", load_markertrack, f))
    for p in args.jr:
        f = resolve_jr(p, JR_TRACK_CANDIDATES if velocity else JR_CSV_CANDIDATES + JR_TRACK_CANDIDATES)
        if f is None:
            what = "surface_track.txt" if velocity else "topography_timeseries.csv (run extract_markerchain.py on its vtk/ dir first)"
            print(f"WARNING: no {what} for {p}, skipping")
            continue
        if f.endswith(".txt"):
            if velocity and count_columns(f) <= VELOCITY_COLUMNS[args.quantity]:
                print(f"WARNING: {f} has no {args.quantity} column, skipping")
                continue
            runs.append((run_label(f, p), "JustRelax surface track", load_markertrack, f))
        else:
            kind = "JustRelax density surface" if "topography_density_csv" in f else "JustRelax marker chain"
            runs.append((run_label(f, p), kind, load_markerchain_csv, f))

    if not runs:
        sys.exit("No inputs found; pass --szu and/or --jr")

    loaded = []
    for label, kind, loader, f in runs:
        print(f"Loading {f} ...")
        if velocity:
            profiles = loader(f, column=VELOCITY_COLUMNS[args.quantity], scale=MS_TO_MMYR)
        else:
            profiles = loader(f)
        print(f"  {len(profiles)} timesteps")
        if args.quantity == "dhdt":
            profiles = height_rate(profiles, args.window_kyr * SECONDS_PER_KYR)
            if not profiles:
                print(f"  WARNING: run shorter than the {args.window_kyr} kyr window, skipping")
                continue
        loaded.append((label, kind, profiles))

    if not loaded:
        sys.exit("Nothing to plot for --quantity " + args.quantity)
    ncols = min(args.ncols, len(loaded))
    nrows = math.ceil(len(loaded) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.5 * ncols, 4.5 * nrows),
                             sharex=True, sharey=True, squeeze=False)
    axes = axes.ravel()

    for ax, (label, kind, profiles) in zip(axes, loaded):
        ylabel = {"h": "topography height [km]", "dhdt": "dh/dt [mm/yr]",
                  "vx": "surface vx [mm/yr]", "vy": "surface vy [mm/yr]"}[args.quantity]
        plot_profiles(ax, profiles, args.n_lines, args.cmap, ylabel)
        if args.quantity != "h":
            ax.axhline(0, color="0.5", lw=0.5)
        ax.set_title(f"{label}\n({kind})", fontsize=10)
        ax.tick_params(labelleft=True, labelbottom=True)
    for ax in axes[len(loaded):]:
        ax.set_visible(False)
    if args.xlim:
        axes[0].set_xlim(args.xlim)
    if args.ylim:
        axes[0].set_ylim(args.ylim)

    title = {"h": "Topography profile evolution",
             "dhdt": f"Rate of surface height change dh/dt, {args.window_kyr:g} kyr window",
             "vx": "Surface horizontal velocity vx (rock side)",
             "vy": "Surface vertical velocity vy (rock side, up > 0)"}[args.quantity]
    fig.suptitle(f"{title} (early = dark, late = yellow)")
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"Saved {args.out}")


if __name__ == "__main__":
    main()

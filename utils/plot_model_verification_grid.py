#!/usr/bin/env python3
"""
plot_model_verification_grid.py

The model-verification figures of plot_model_verification.py combined into one grid: one
row per run, columns lithology | effective viscosity | deviatoric tau_xx | v_y, with shared
x and z axes, one colour bar per column above the top row and the first run's time top left.

Each run is given as  --run LABEL {jr,szu} PATH STEP  (STEP may be 'latest' for jr).

Usage:
    python3 plot_model_verification_grid.py \\
        --run Published szu /scratch/tectonics/bert/SZU2019-gps 217 \\
        --run Simplified szu /scratch/tectonics/bert/SZU2019-gps-noerossedserpfricyn 217 \\
        --run "Current model" jr <run>/vtk 1775 \\
        --vlim-vy 1 -o "z_output/Model verification/model_verification_grid.png"
"""

import argparse
import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

from plot_model_verification import (CATEGORIES, CMAP_DIV, CMAP_ETA, TEMP_LEVELS,
                                     latest_step, load_jr, load_szu)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", nargs=4, action="append", required=True, metavar=("LABEL", "TYPE", "PATH", "STEP"))
    ap.add_argument("--xlim", type=float, nargs=2, default=(800, 1100))
    ap.add_argument("--ylim", type=float, nargs=2, default=(-80, 0))
    ap.add_argument("--clim-eta", type=float, nargs=2, default=(18, 25))
    ap.add_argument("--vlim-normal", type=float, default=50)
    ap.add_argument("--vlim-vy", type=float, default=1)
    ap.add_argument("--res-km", type=float, default=0.2)
    ap.add_argument("--zticks", type=float, default=20, help="depth tick spacing [km] (default 20)")
    ap.add_argument("--fontsize", type=float, default=20)
    ap.add_argument("-o", "--out", default="model_verification_grid.png")
    args = ap.parse_args()
    plt.rcParams.update({"font.size": args.fontsize})

    data = []
    for label, kind, path, step in args.run:
        if kind == "jr":
            d = load_jr(path, latest_step(path) if step == "latest" else int(step), args.xlim, args.ylim, args.res_km)
        else:
            d = load_szu(path, int(step), args.xlim, args.ylim, args.res_km)
        print(f"{label}: t = {d['t_myr']:.3f} Myr")
        data.append((label, d))

    ids = sorted(CATEGORIES)
    lith_cmap = ListedColormap([CATEGORIES[i][1] for i in ids])
    lith_norm = BoundaryNorm(np.append(np.array(ids) - 0.5, ids[-1] + 0.5), lith_cmap.N)
    columns = [  # (field, title, pcolormesh kwargs, colour bar label, isotherm colour)
        ("lith", "Lithology", dict(cmap=lith_cmap, norm=lith_norm), None, "k"),
        ("eta", "Effective viscosity", dict(cmap=CMAP_ETA, vmin=args.clim_eta[0], vmax=args.clim_eta[1]),
         "log$_{10}$ $\\eta_{eff}$ (Pa s)", "w"),
        ("txx", "Deviatoric $\\tau_{xx}$", dict(cmap=CMAP_DIV, vmin=-args.vlim_normal, vmax=args.vlim_normal),
         "$\\tau_{xx}$ (MPa)", "k"),
        ("vy", "Vertical velocity $v_y$", dict(cmap=CMAP_DIV, vmin=-args.vlim_vy, vmax=args.vlim_vy),
         "$v_y$ (cm/yr)", "k"),
    ]

    nrow, ncol = len(data), len(columns)
    panel_w = 9.0
    panel_h = panel_w * (args.ylim[1] - args.ylim[0]) / (args.xlim[1] - args.xlim[0])
    fig, axs = plt.subplots(nrow, ncol, figsize=(ncol * panel_w, nrow * (panel_h + 0.5) + 1.4),
                            sharex=True, sharey=True, squeeze=False, layout="constrained")
    fig.get_layout_engine().set(h_pad=0.02, hspace=0.0)  # tight rows
    for r, (label, d) in enumerate(data):
        for c, (field, title, kw, _, iso) in enumerate(columns):
            ax = axs[r, c]
            x, y, v = d[field]
            im = ax.pcolormesh(x, y, v, shading="nearest", **kw)
            tx, ty, tv = d["T"]
            cs = ax.contour(tx, ty, tv, levels=TEMP_LEVELS, colors=iso, linewidths=0.8)
            ax.clabel(cs, inline=True, fontsize=0.6 * args.fontsize, fmt="%d C")
            ax.set_title(f"{label}: {title}")
            ax.set_xlim(args.xlim)
            ax.set_ylim(args.ylim)
            ax.set_aspect("equal")
            ax.set_yticks(np.arange(args.ylim[0], args.ylim[1] + 0.1, args.zticks))
            if r == nrow - 1:
                ax.set_xlabel("x (km)")
            if c == 0:
                ax.set_ylabel("z (km)")
            if r == nrow - 1:  # one colour bar per column, above the top row
                cb = fig.colorbar(im, ax=axs[:, c], location="top", shrink=0.8, pad=0.02)
                if columns[c][3]:
                    cb.set_label(columns[c][3])
                else:  # lithology: no colour bar; its slot shows the time instead
                    cb.ax.clear()
                    cb.ax.axis("off")
                    cb.ax.text(0.5, 0.5, f"Time: {data[0][1]['t_myr']:.3f} Myr",
                               ha="center", va="center", transform=cb.ax.transAxes,
                               fontsize=args.fontsize * 1.1, fontweight="bold")

    for label, d in data:  # the figure shows the first run's time; the others go to the terminal
        print(f"  {label}: {d['t_myr']:.3f} Myr")
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"Saved {args.out}")


if __name__ == "__main__":
    main()

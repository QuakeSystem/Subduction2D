#!/usr/bin/env python3
"""
plot_viscosity_from_vtk.py

Effective visco-elasto-plastic viscosity log10(eta_vep) of one JustRelax output step
(vtk_<step>.vtr), with temperature contours, in the same layout as the stress and lithology
plots. Air (density < 1000 kg/m^3) is masked. Colormap: Crameri's batlow (perceptually
uniform and colour-blind friendly; falls back to viridis without cmcrameri).

Usage:
    python3 plot_viscosity_from_vtk.py <run>/vtk --step 1175 [--xlim 800 1100 --ylim -80 0]
        [--clim 18 25] [--cmap batlow] [--fontsize 20]
"""

import argparse
import os

import numpy as np
import matplotlib.pyplot as plt

from extract_surface_from_vtk import read_vtr_arrays
from plot_stress_from_vtk import vtk_time_myr, TEMP_LEVELS


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("vtk_dir")
    ap.add_argument("--step", type=int, required=True)
    ap.add_argument("--xlim", type=float, nargs=2, default=(800, 1100))
    ap.add_argument("--ylim", type=float, nargs=2, default=(-80, 0))
    ap.add_argument("--clim", type=float, nargs=2, default=(18, 25), help="log10(eta_vep) colour range (default 18 25)")
    ap.add_argument("--cmap", default="batlow", help="cmcrameri colormap name (default batlow)")
    ap.add_argument("--fontsize", type=float, default=12)
    ap.add_argument("-o", "--outdir", default=".")
    args = ap.parse_args()
    plt.rcParams.update({"font.size": args.fontsize})

    try:
        from cmcrameri import cm
        cmap = getattr(cm, args.cmap)
    except ImportError:
        print("cmcrameri not found, using viridis")
        cmap = "viridis"

    fname = f"vtk_{args.step:06d}.vtr"
    xv, yv, a = read_vtr_arrays(os.path.join(args.vtk_dir, fname), ["η_vep", "thermal_T", "density"])
    xc, yc = 0.5 * (xv[:-1] + xv[1:]) / 1e3, 0.5 * (yv[:-1] + yv[1:]) / 1e3
    eta = np.where(a["density"] < 1000.0, np.nan, np.log10(a["η_vep"]))

    fig, ax = plt.subplots(figsize=(18, 6), constrained_layout=True)
    im = ax.pcolormesh(xc, yc, eta, cmap=cmap, vmin=args.clim[0], vmax=args.clim[1], shading="nearest")
    fig.colorbar(im, ax=ax, pad=0.0).set_label("log$_{10}$ $\\eta_{vep}$ (Pa s)")
    cs = ax.contour(xc, yc, a["thermal_T"] - 273.15, levels=TEMP_LEVELS, colors="w", linewidths=0.8)
    ax.clabel(cs, inline=True, fontsize=0.7 * args.fontsize, fmt="%d C")
    ax.set(xlabel="x (km)", ylabel="z (km)", xlim=args.xlim, ylim=args.ylim,
           title="Effective viscosity $\\eta_{vep}$")
    ax.set_aspect("equal")
    fig.suptitle(f"Time: {vtk_time_myr(args.vtk_dir, fname):.3f} Myr")
    os.makedirs(args.outdir, exist_ok=True)
    out = os.path.join(args.outdir, f"viscosity_vep_{args.step}.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"Saved {out}")


if __name__ == "__main__":
    main()

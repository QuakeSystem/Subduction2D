#!/usr/bin/env python3
"""
plot_lithology_from_vtk.py

Marker lithology of one JustRelax output step: the particles of particles_<step>.vtu are
rasterised onto a pixel grid (each pixel takes the phase of a particle inside it), coloured
by phase grouped into domains as defined in Subduction2D_setup.jl, with temperature
contours from vtk_<step>.vtr.

  Air          2
  Mantle       1 asthenosphere + oceanic and continental lithospheric mantle, 3 weak zone
  Oceanic      4 hydrated top crust / thrust interface, 5 gabbro
  Continental  6 felsic crust, 7 passive-margin sediments

With --outline, the units are drawn as outlines only (no fill, no legend) with dashed
isotherms, e.g. as a background for other figures; add --transparent for a transparent PNG.

Usage:
    python3 plot_lithology_from_vtk.py <run>/vtk --step 1175 --xlim 800 1100 --ylim -80 0
    python3 plot_lithology_from_vtk.py <run>/vtk --step 1175 --outline --transparent
    python3 plot_lithology_from_vtk.py <run>/vtk --step 1175 --xlim 0 1500 --ylim -260 0 \\
        --res-km 1.5 --figsize 26 7          # whole domain
"""

import argparse
import os
import re

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch

from extract_surface_from_vtk import _decode_block, read_vtr_arrays
from plot_stress_from_vtk import vtk_time_myr, TEMP_LEVELS

# phase -> (label, colour); colours grouped by domain
PHASES = {
    2: ("Air", "#ffffff"),
    1: ("Mantle", "#c9b8dd"),
    3: ("Mantle weak zone", "#7b5ea7"),
    4: ("Oceanic: hydrated top crust / interface", "#2a9d8f"),
    5: ("Oceanic: gabbro", "#1d5c4f"),
    6: ("Continental: felsic crust", "#e9a46a"),
    7: ("Continental: sediments", "#f4d35e"),
}


def read_vtu_points_phase(path):
    """Particle coordinates (km) and phase from a WriteVTK particles .vtu file."""
    with open(path, "rb") as f:
        data = f.read()
    start = data.find(b"<AppendedData")
    raw = data[data.find(b"_", start) + 1:]
    header = data[:start].decode("utf-8", "replace")
    off = {m.group(1): int(m.group(2))
           for m in re.finditer(r'<DataArray [^>]*Name="([^"]+)"[^>]*offset="(\d+)"', header)}
    pts = _decode_block(raw, off["Points"], "<f4").reshape(-1, 3)
    phase = _decode_block(raw, off["phase"], "<f4")
    x, y = pts[:, 0].astype(float), pts[:, 1].astype(float)
    if np.nanmax(np.abs(x)) > 1.0e4:  # written in metres rather than km
        x, y = x / 1e3, y / 1e3
    return x, y, np.rint(phase).astype(int)


def rasterise(x, y, phase, xlim, ylim, res):
    """Pixel image of particle phases on [xlim] x [ylim] with pixel size res (all km)."""
    nx = int(np.ceil((xlim[1] - xlim[0]) / res))
    ny = int(np.ceil((ylim[1] - ylim[0]) / res))
    keep = (x >= xlim[0]) & (x < xlim[1]) & (y >= ylim[0]) & (y < ylim[1])
    i = ((x[keep] - xlim[0]) / res).astype(int)
    j = ((y[keep] - ylim[0]) / res).astype(int)
    img = np.full((ny, nx), np.nan)
    img[j, i] = phase[keep]  # one particle's phase per pixel
    xp = xlim[0] + (np.arange(nx) + 0.5) * res
    yp = ylim[0] + (np.arange(ny) + 0.5) * res
    return xp, yp, img


def fill_empty(img):
    """Give pixels without a particle the phase of the nearest pixel that has one."""
    from scipy import ndimage
    idx = ndimage.distance_transform_edt(np.isnan(img), return_distances=False, return_indices=True)
    return img[tuple(idx)]


def clean_phases(img, size=5):
    """Fill empty pixels and remove single-pixel speckle, so unit outlines follow the unit
    boundaries rather than particle noise."""
    from scipy import ndimage
    return ndimage.median_filter(fill_empty(img), size=size)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("vtk_dir")
    ap.add_argument("--step", type=int, required=True)
    ap.add_argument("--xlim", type=float, nargs=2, default=(800, 1100))
    ap.add_argument("--ylim", type=float, nargs=2, default=(-80, 0))
    ap.add_argument("--res-km", type=float, default=0.2, help="pixel size [km] (default 0.2)")
    ap.add_argument("--fontsize", type=float, default=12, help="base font size (default 12)")
    ap.add_argument("--figsize", type=float, nargs=2, default=(18, 6), metavar=("W", "H"), help="inches")
    ap.add_argument("--outline", action="store_true", help="unit outlines only, no fill (background figure)")
    ap.add_argument("--transparent", action="store_true", help="save with a transparent background")
    ap.add_argument("-o", "--outdir", default=".")
    args = ap.parse_args()
    plt.rcParams.update({"font.size": args.fontsize})

    x, y, phase = read_vtu_points_phase(os.path.join(args.vtk_dir, f"particles_{args.step:06d}.vtu"))
    xp, yp, img = rasterise(x, y, phase, args.xlim, args.ylim, args.res_km)

    fname = f"vtk_{args.step:06d}.vtr"
    xv, yv, a = read_vtr_arrays(os.path.join(args.vtk_dir, fname), ["thermal_T"])
    xc, yc = 0.5 * (xv[:-1] + xv[1:]) / 1e3, 0.5 * (yv[:-1] + yv[1:]) / 1e3

    ids = sorted(PHASES)
    cmap = ListedColormap([PHASES[p][1] for p in ids])
    norm = BoundaryNorm(np.append(np.array(ids) - 0.5, ids[-1] + 0.5), cmap.N)

    fig, ax = plt.subplots(figsize=args.figsize, constrained_layout=True)
    if args.outline:
        clean = clean_phases(img)
        for p in ids:  # each unit's boundary; shared boundaries are drawn twice on the same line
            if np.any(clean == p):
                ax.contour(xp, yp, (clean == p).astype(float), levels=[0.5], colors="k", linewidths=1.2)
        iso_style = dict(colors="0.4", linewidths=1.0, linestyles="--")
    else:
        ax.pcolormesh(xp, yp, fill_empty(img), cmap=cmap, norm=norm, shading="nearest")
        iso_style = dict(colors="k", linewidths=0.8)
    cs = ax.contour(xc, yc, a["thermal_T"] - 273.15, levels=TEMP_LEVELS, **iso_style)
    ax.clabel(cs, inline=True, fontsize=0.7 * args.fontsize, fmt="%d C")
    ax.set(xlabel="x (km)", ylabel="z (km)", xlim=args.xlim, ylim=args.ylim,
           title="Lithology (outlines)" if args.outline else "Lithology (markers)")
    ax.set_aspect("equal")
    if not args.outline:
        present = [p for p in ids if np.any(img == p)]
        ax.legend(handles=[Patch(facecolor=PHASES[p][1], edgecolor="0.5", label=f"{p}: {PHASES[p][0]}") for p in present],
                  loc="upper left", bbox_to_anchor=(1.01, 1), fontsize=0.85 * args.fontsize)
    fig.suptitle(f"Time: {vtk_time_myr(args.vtk_dir, fname):.3f} Myr")
    os.makedirs(args.outdir, exist_ok=True)
    out = os.path.join(args.outdir, f"lithology_{'outline_' if args.outline else ''}{args.step}.png")
    fig.savefig(out, dpi=200, bbox_inches="tight", transparent=args.transparent)
    print(f"Saved {out}")


if __name__ == "__main__":
    main()

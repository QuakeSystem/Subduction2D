#!/usr/bin/env python3
"""
plot_stress_from_vtk.py

Plot the deviatoric normal (tau_xx) and shear (tau_xy) stress and the second invariant of
the deviatoric stress (tau_II) of one output step as three separate figures, with temperature contours, in the style of the marker stress plots of the
reference visualisation routines. Air (density < 1000 kg/m^3) is masked.

The input is either
  - a JustRelax vtk/ directory: vtk_<step>.vtr (cell-centred fields), or
  - an SZU2019 run directory: lphr<step>.gzip.h5 (nodal sxx, sxy, tk, ro; needs h5py).
    tau_II is computed at cell centres as in loadevps_omp.c (postproc):
    sqrt(sxx^2 + (mean of the 4 surrounding sxy)^2).
    SZU depth y (positive down from the top of the sticky air) is plotted as z = -y,
    the same frame as JustRelax (z = 0 at the top of the domain); sxy is negated accordingly.

Usage:
    python3 plot_stress_from_vtk.py <run>/vtk --step 1175 \\
        --xlim 800 1100 --ylim -80 0 [--vlim-normal 100 --vlim-shear 50] [-o outdir]
    python3 plot_stress_from_vtk.py /scratch/tectonics/bert/SZU2019-gps-noerossedserpfricyn --step 105
"""

import argparse
import os
import xml.etree.ElementTree as ET

import numpy as np
import matplotlib.pyplot as plt

from extract_surface_from_vtk import read_vtr_arrays

try:
    from cmcrameri import cm
    CMAP, CMAP_II = cm.vik, cm.batlow
except ImportError:
    CMAP, CMAP_II = "RdBu_r", "viridis"

SECONDS_PER_MYR = 3600 * 24 * 365.25 * 1e6
TEMP_LEVELS = [10, 100, 150, 350, 450, 1300]  # deg C


def vtk_time_myr(vtk_dir, fname):
    for ds in ET.parse(os.path.join(vtk_dir, "Subduction2D.pvd")).getroot().iter("DataSet"):
        if ds.get("file") == fname:
            return float(ds.get("timestep")) / SECONDS_PER_MYR
    return float("nan")


def load_jr(vtk_dir, step):
    fname = f"vtk_{step:06d}.vtr"
    xv, yv, a = read_vtr_arrays(os.path.join(vtk_dir, fname), ["τxx", "τxy", "τII", "thermal_T", "density"])
    xc, yc = 0.5 * (xv[:-1] + xv[1:]) / 1e3, 0.5 * (yv[:-1] + yv[1:]) / 1e3
    air = a["density"] < 1000.0
    tII = (xc, yc, np.where(air, np.nan, a["τII"]))  # cell centres, like the other fields
    return xc, yc, a["τxx"], a["τxy"], tII, a["thermal_T"], a["density"], vtk_time_myr(vtk_dir, fname)


def load_szu(run_dir, step):
    import h5py
    with h5py.File(os.path.join(run_dir, f"lphr{step}.gzip.h5"), "r") as f:
        t_myr = f["ModelGroup/Model"][0] / 1e6  # years
        gx, gy = f["ModelGroup/gx"][:], f["ModelGroup/gy"][:]
        # nodal arrays are stored x-major (index = i*ny + j); transpose to (ny, nx) like the VTK fields
        g = lambda k: f["NodeGroup/" + k][:].reshape(len(gx), len(gy)).T
        sxx, sxy, ro, tk = g("sxx"), g("sxy"), g("ro"), g("tk")
    # tau_II at cell centres, as SZU computes sii: sxx of the cell, sxy averaged from its 4 nodes
    sxy4 = 0.25 * (sxy[:-1, :-1] + sxy[1:, :-1] + sxy[:-1, 1:] + sxy[1:, 1:])
    sii = np.sqrt(sxx[:-1, :-1] ** 2 + sxy4 ** 2)
    air4 = np.minimum.reduce([ro[:-1, :-1], ro[1:, :-1], ro[:-1, 1:], ro[1:, 1:]]) < 1000.0
    tII = (0.5 * (gx[:-1] + gx[1:]) / 1e3, -0.5 * (gy[:-1] + gy[1:]) / 1e3, np.where(air4, np.nan, sii))
    # SZU y points down; flipping to z = -y also flips the sign of the xy component
    return gx / 1e3, -gy / 1e3, sxx, -sxy, tII, tk, ro, t_myr


def plot_field(xc, yc, field, T_c, title, label, vlim, xlim, ylim, t_myr, path,
               field_xy=None, cmap=CMAP, vmin=None):
    """field_xy: (x, y) of `field` if it lives on other points than the temperature (xc, yc)."""
    fx, fy = field_xy if field_xy is not None else (xc, yc)
    fig, ax = plt.subplots(figsize=(18, 6), constrained_layout=True)
    im = ax.pcolormesh(fx, fy, field, cmap=cmap, vmin=-vlim if vmin is None else vmin, vmax=vlim, shading="nearest")
    fig.colorbar(im, ax=ax, pad=0.0).set_label(label)
    cs = ax.contour(xc, yc, T_c, levels=TEMP_LEVELS, colors="k", linewidths=0.8)
    ax.clabel(cs, inline=True, fontsize=0.7 * plt.rcParams["font.size"], fmt="%d C")
    ax.set(xlabel="x (km)", ylabel="z (km)", xlim=xlim, ylim=ylim, title=title)
    ax.set_aspect("equal")
    fig.suptitle(f"Time: {t_myr:.3f} Myr")
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("vtk_dir", help="JustRelax vtk/ directory or SZU2019 run directory")
    ap.add_argument("--step", type=int, required=True)
    ap.add_argument("--xlim", type=float, nargs=2, default=(800, 1100))
    ap.add_argument("--ylim", type=float, nargs=2, default=(-80, 0))
    ap.add_argument("--vlim-normal", type=float, help="colour limit tau_xx [MPa] (default: 99th pct in view)")
    ap.add_argument("--vlim-shear", type=float, help="colour limit tau_xy [MPa] (default: 99th pct in view)")
    ap.add_argument("--vlim-ii", type=float, help="upper colour limit tau_II [MPa] (default: 99th pct in view)")
    ap.add_argument("--fields", nargs="+", choices=("normal", "shear", "ii"), default=("normal", "shear", "ii"),
                    help="which figures to make (default: all three)")
    ap.add_argument("--fontsize", type=float, default=10, help="base font size (default 10)")
    ap.add_argument("-o", "--outdir", default=".")
    args = ap.parse_args()
    plt.rcParams.update({"font.size": args.fontsize})

    szu = os.path.isfile(os.path.join(args.vtk_dir, f"lphr{args.step}.gzip.h5"))
    xc, yc, sxx, sxy, (xi, yi, sii), T_K, rho, t_myr = (load_szu if szu else load_jr)(args.vtk_dir, args.step)
    air = rho < 1000.0
    txx = np.where(air, np.nan, sxx / 1e6)
    txy = np.where(air, np.nan, sxy / 1e6)
    T_c = T_K - 273.15
    tag = f"{os.path.basename(os.path.normpath(args.vtk_dir))}_" if szu else ""

    view = (xc[None, :] >= args.xlim[0]) & (xc[None, :] <= args.xlim[1]) & \
           (yc[:, None] >= args.ylim[0]) & (yc[:, None] <= args.ylim[1])
    vn = args.vlim_normal or np.nanpercentile(np.abs(txx[view]), 99)
    vs = args.vlim_shear or np.nanpercentile(np.abs(txy[view]), 99)
    tii = sii / 1e6
    view_ii = (xi[None, :] >= args.xlim[0]) & (xi[None, :] <= args.xlim[1]) & \
              (yi[:, None] >= args.ylim[0]) & (yi[:, None] <= args.ylim[1])
    vii = args.vlim_ii or np.nanpercentile(tii[view_ii], 99)

    os.makedirs(args.outdir, exist_ok=True)
    if "normal" in args.fields:
        plot_field(xc, yc, txx, T_c, "Normal stress $\\tau_{xx}$ (deviatoric)", "$\\tau_{xx}$ (MPa)", vn,
               args.xlim, args.ylim, t_myr, os.path.join(args.outdir, f"stress_normal_{tag}{args.step}.png"))
    if "shear" in args.fields:
        plot_field(xc, yc, txy, T_c, "Shear stress $\\tau_{xy}$", "$\\tau_{xy}$ (MPa)", vs,
               args.xlim, args.ylim, t_myr, os.path.join(args.outdir, f"stress_shear_{tag}{args.step}.png"))
    if "ii" in args.fields:
        plot_field(xc, yc, tii, T_c, "Second invariant of the deviatoric stress $\\tau_{II}$", "$\\tau_{II}$ (MPa)", vii,
               args.xlim, args.ylim, t_myr, os.path.join(args.outdir, f"stress_II_{tag}{args.step}.png"),
               field_xy=(xi, yi), cmap=CMAP_II, vmin=0.0)


if __name__ == "__main__":
    main()

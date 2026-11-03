#!/usr/bin/env python3
"""
plot_model_verification.py

Model-verification figures for one run, JustRelax or SZU2019, all on the same domain and
colour scales so the runs can be compared side by side:

  lithology        marker lithology grouped into the setup's domains
  viscosity        log10 effective viscosity (JustRelax eta_vep, SZU nu), batlow_r
  normal_stress    deviatoric normal stress tau_xx (MPa), vik
  vy               vertical velocity (cm/yr, up > 0), vik

each with temperature contours, air masked. One PNG per field: <field>_<label>.png.

Inputs
  --jr  <run>/vtk  --step N|latest   vtk_<step>.vtr + particles_<step>.vtu
  --szu <SZU run>  --step N          lphr<step>.gzip.h5 (needs h5py). SZU depth (y down)
        is plotted as z = -y; tau_xx keeps its sign, vy is negated (y up).

Usage:
    python3 plot_model_verification.py --szu /scratch/tectonics/bert/SZU2019-gps --step 217 \\
        --label Published -o "z_output/Model verification"
    python3 plot_model_verification.py --jr <run>/vtk --step latest --label "Current model" -o ...
"""

import argparse
import glob
import os
import re

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch

from extract_surface_from_vtk import read_vtr_arrays
from plot_lithology_from_vtk import read_vtu_points_phase, rasterise, fill_empty
from plot_stress_from_vtk import vtk_time_myr, TEMP_LEVELS

try:
    from cmcrameri import cm
    CMAP_DIV, CMAP_ETA = cm.vik, cm.batlow_r
except ImportError:
    print("cmcrameri not found, using RdBu_r / viridis_r")
    CMAP_DIV, CMAP_ETA = "RdBu_r", "viridis_r"

MS_TO_CMYR = 100.0 * 3600 * 24 * 365.25

# lithology categories shared by both codes: id -> (label, colour)
CATEGORIES = {
    0: ("Air", "#ffffff"),
    1: ("Continental: sediments", "#f4d35e"),
    2: ("Continental: felsic crust", "#e9a46a"),
    3: ("Continental: hydrated felsic crust", "#b5652a"),
    4: ("Oceanic: hydrated top crust / interface", "#2a9d8f"),
    5: ("Oceanic: gabbro", "#1d5c4f"),
    6: ("Mantle", "#c9b8dd"),
    7: ("Mantle weak / shear zone", "#7b5ea7"),
    8: ("Hydrated / serpentinised mantle", "#b07aa1"),
    9: ("Molten rock", "#d1495b"),
}
# JustRelax phases (Subduction2D_setup.jl) -> category
JR_CATEGORY = {1: 6, 2: 0, 3: 7, 4: 4, 5: 5, 6: 2, 7: 1}
# SZU2019 rock types (init.t3c rock list) -> category; 18 and up are molten rocks
SZU_CATEGORY = {0: 0, 1: 0, 2: 1, 3: 1, 4: 1, 5: 2, 6: 2, 7: 4, 8: 5, 9: 6, 10: 6, 11: 8,
                12: 7, 13: 8, 14: 6, 15: 5, 16: 3, 17: 3}


def to_category(phase_img, mapping, default):
    out = np.full(phase_img.shape, np.nan)
    finite = np.isfinite(phase_img)
    vals = phase_img[finite].astype(int)
    out[finite] = [mapping.get(v, default) for v in vals]
    return out


# ---------------------------------------------------------------------------------------
# loaders: return a dict with km coordinates and fields on their own points

def latest_step(vtk_dir):
    steps = [int(re.search(r"(\d+)", os.path.basename(p)).group(1))
             for p in glob.glob(os.path.join(vtk_dir, "particles_*.vtu"))]
    steps = [s for s in steps if os.path.isfile(os.path.join(vtk_dir, f"vtk_{s:06d}.vtr"))]
    return max(steps)


def load_jr(vtk_dir, step, xlim, ylim, res_km):
    fname = f"vtk_{step:06d}.vtr"
    xv, yv, a = read_vtr_arrays(os.path.join(vtk_dir, fname),
                                ["thermal_T", "density", "η_vep", "τxx", "Vy"])
    xv, yv = xv / 1e3, yv / 1e3
    xc, yc = 0.5 * (xv[:-1] + xv[1:]), 0.5 * (yv[:-1] + yv[1:])
    air_c = a["density"] < 1000.0
    # a vertex is air if all its neighbouring cells are
    pad = np.pad(air_c, 1, mode="edge")
    air_v = pad[:-1, :-1] & pad[1:, :-1] & pad[:-1, 1:] & pad[1:, 1:]
    x, y, phase = read_vtu_points_phase(os.path.join(vtk_dir, f"particles_{step:06d}.vtu"))
    xp, yp, img = rasterise(x, y, phase, xlim, ylim, res_km)
    return dict(
        t_myr=vtk_time_myr(vtk_dir, fname),
        T=(xc, yc, a["thermal_T"] - 273.15),
        lith=(xp, yp, to_category(fill_empty(img), JR_CATEGORY, 6)),
        eta=(xc, yc, np.where(air_c, np.nan, np.log10(a["η_vep"]))),
        txx=(xc, yc, np.where(air_c, np.nan, a["τxx"] / 1e6)),
        vy=(xv, yv, np.where(air_v, np.nan, a["Vy"] * MS_TO_CMYR)),
    )


def decode_szu_composition(mtype, nx, ny):
    """Run-length decoding of VisMarkerGroup/Mtype as in SZU's vis_hdf5.m (column-major)."""
    out, k = [], 0
    while k < len(mtype):
        v = mtype[k]
        if v == -2:
            out.extend([mtype[k + 2]] * mtype[k + 1])
            k += 3
        elif v == -1:
            out.append(0)
            k += 1
        else:
            out.append(v)
            k += 1
    return np.array(out[:nx * ny], dtype=float).reshape((ny, nx), order="F")


def load_szu(run_dir, step, xlim, ylim, res_km):
    import h5py
    with h5py.File(os.path.join(run_dir, f"lphr{step}.gzip.h5"), "r") as f:
        model = f["ModelGroup/Model"][:]
        gx, gy = f["ModelGroup/gx"][:] / 1e3, f["ModelGroup/gy"][:] / 1e3
        g = lambda k: f["NodeGroup/" + k][:].reshape(len(gx), len(gy)).T  # x-major -> (ny, nx)
        tk, ro, nu, sxx, vy = g("tk"), g("ro"), g("nu"), g("sxx"), g("vy")
        mtype = f["VisMarkerGroup/Mtype"][:].astype(int)
    air = ro < 1000.0
    z = -gy
    # composition grid: resolution = finest grid step; pixel (r, c) sits at ((c+1) res, (r+1) res)
    res = np.min(np.diff(gy))
    nxv, nyv = int(model[1] / 1e3 / res + 1), int(model[2] / 1e3 / res + 1)
    comp = decode_szu_composition(mtype, nxv, nyv)
    cx, cz = (np.arange(nxv) + 1) * res, -(np.arange(nyv) + 1) * res
    return dict(
        t_myr=model[0] / 1e6,
        T=(gx, z, tk - 273.15),
        lith=(cx, cz, to_category(comp, SZU_CATEGORY, 9)),
        eta=(gx, z, np.where(air, np.nan, np.log10(nu))),
        txx=(gx, z, np.where(air, np.nan, sxx / 1e6)),
        vy=(gx, z, np.where(air, np.nan, -vy * MS_TO_CMYR)),
    )


# ---------------------------------------------------------------------------------------

def figure(d, field_xyz, title, args, cbar_label=None, iso_colour="k", **pcolor_kw):
    fig, ax = plt.subplots(figsize=(18, 6), constrained_layout=True)
    x, y, v = field_xyz
    im = ax.pcolormesh(x, y, v, shading="nearest", **pcolor_kw)
    if cbar_label:
        fig.colorbar(im, ax=ax, pad=0.0).set_label(cbar_label)
    tx, ty, tv = d["T"]
    cs = ax.contour(tx, ty, tv, levels=TEMP_LEVELS, colors=iso_colour, linewidths=0.8)
    ax.clabel(cs, inline=True, fontsize=0.7 * args.fontsize, fmt="%d C")
    ax.set(xlabel="x (km)", ylabel="z (km)", xlim=args.xlim, ylim=args.ylim, title=f"{args.label}: {title}")
    ax.set_aspect("equal")
    fig.suptitle(f"Time: {d['t_myr']:.3f} Myr")
    return fig, ax


def save(fig, args, field):
    slug = re.sub(r"[^A-Za-z0-9]+", "_", args.label).strip("_")
    out = os.path.join(args.outdir, f"{field}_{slug}.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--jr", help="JustRelax vtk/ directory")
    src.add_argument("--szu", help="SZU2019 run directory")
    ap.add_argument("--step", required=True, help="output step, or 'latest' (JustRelax)")
    ap.add_argument("--label", required=True, help="run name used in titles and file names")
    ap.add_argument("--xlim", type=float, nargs=2, default=(800, 1100))
    ap.add_argument("--ylim", type=float, nargs=2, default=(-80, 0))
    ap.add_argument("--clim-eta", type=float, nargs=2, default=(18, 25), help="log10 viscosity range")
    ap.add_argument("--vlim-normal", type=float, default=50, help="tau_xx colour limit [MPa]")
    ap.add_argument("--vlim-vy", type=float, default=5, help="vy colour limit [cm/yr]")
    ap.add_argument("--res-km", type=float, default=0.2, help="JustRelax lithology pixel size [km]")
    ap.add_argument("--fontsize", type=float, default=20)
    ap.add_argument("-o", "--outdir", default=".")
    args = ap.parse_args()
    plt.rcParams.update({"font.size": args.fontsize})
    os.makedirs(args.outdir, exist_ok=True)

    if args.jr:
        step = latest_step(args.jr) if args.step == "latest" else int(args.step)
        d = load_jr(args.jr, step, args.xlim, args.ylim, args.res_km)
    else:
        d = load_szu(args.szu, int(args.step), args.xlim, args.ylim, args.res_km)
    print(f"{args.label}: t = {d['t_myr']:.3f} Myr")

    # lithology
    ids = sorted(CATEGORIES)
    cmap = ListedColormap([CATEGORIES[i][1] for i in ids])
    norm = BoundaryNorm(np.append(np.array(ids) - 0.5, ids[-1] + 0.5), cmap.N)
    fig, ax = figure(d, d["lith"], "Lithology", args, cmap=cmap, norm=norm)
    lx, ly, lv = d["lith"]
    inview = (lx[None, :] >= args.xlim[0]) & (lx[None, :] <= args.xlim[1]) & \
             (ly[:, None] >= args.ylim[0]) & (ly[:, None] <= args.ylim[1])
    present = [i for i in ids if np.any(lv[inview] == i)]
    ax.legend(handles=[Patch(facecolor=CATEGORIES[i][1], edgecolor="0.5", label=CATEGORIES[i][0]) for i in present],
              loc="upper left", bbox_to_anchor=(1.01, 1), fontsize=0.85 * args.fontsize)
    save(fig, args, "lithology")

    fig, _ = figure(d, d["eta"], "Effective viscosity", args, "log$_{10}$ $\\eta_{eff}$ (Pa s)",
                    iso_colour="w", cmap=CMAP_ETA, vmin=args.clim_eta[0], vmax=args.clim_eta[1])
    save(fig, args, "viscosity")

    fig, _ = figure(d, d["txx"], "Deviatoric normal stress $\\tau_{xx}$", args, "$\\tau_{xx}$ (MPa)",
                    cmap=CMAP_DIV, vmin=-args.vlim_normal, vmax=args.vlim_normal)
    save(fig, args, "normal_stress")

    fig, _ = figure(d, d["vy"], "Vertical velocity $v_y$", args, "$v_y$ (cm/yr)",
                    cmap=CMAP_DIV, vmin=-args.vlim_vy, vmax=args.vlim_vy)
    save(fig, args, "vy")


if __name__ == "__main__":
    main()

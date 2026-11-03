#!/usr/bin/env python3
"""
extract_surface_from_vtk.py

Extract the model's air/rock surface from the gridded JustRelax VTK output
(vtk_NNNNNN.vtr + Subduction2D.pvd), as an alternative to the marker chain.

For every grid column, the surface is where the cell-centred `density` field,
scanned from the top down, first rises through a threshold (default
1000 kg/m^3, i.e. well between air and any rock phase). The height is linearly
interpolated between the two bracketing cell centres, so it resolves the
surface below one cell height: cell densities are phase-ratio mixtures, so a
partly-air cell has an intermediate density.

Writes the same long-format CSV as extract_markerchain.py, so
plot_topography_comparison.py can plot it directly (pass the CSV file to --jr):
    <vtk dir>/<outdir>/topography_timeseries.csv
    columns: time_s, time_kyr, x_km, h_km

The .vtr files are parsed directly (raw appended, zlib-compressed, UInt64
headers, as written by WriteVTK.jl). The `vtk` Python reader segfaults on them
(duplicate array names), so it is not used.

Usage:
    python3 extract_surface_from_vtk.py /path/to/run/vtk
    python3 extract_surface_from_vtk.py /path/to/run/vtk --threshold 1500 --smooth-km 2
"""

import argparse
import csv
import os
import re
import sys
import xml.etree.ElementTree as ET
import zlib

import numpy as np

SECONDS_PER_KYR = 3600 * 24 * 365.25 * 1000
VTK_DTYPES = {"Float32": "<f4", "Float64": "<f8", "Int32": "<i4", "Int64": "<i8", "UInt8": "<u1"}


def read_pvd(pvd_path):
    """Return an ordered list of (time_seconds, filename) pairs, duplicates removed."""
    entries, seen = [], set()
    for ds in ET.parse(pvd_path).getroot().iter("DataSet"):
        key = (ds.get("file"), float(ds.get("timestep")))
        if key not in seen:
            seen.add(key)
            entries.append((key[1], key[0]))
    return entries


def _decode_block(raw, offset, dtype):
    """Decode one zlib-compressed appended DataArray starting at `offset` (UInt64 header)."""
    hdr = np.frombuffer(raw, dtype="<u8", count=3, offset=offset)
    nblocks = int(hdr[0])
    csizes = np.frombuffer(raw, dtype="<u8", count=nblocks, offset=offset + 24)
    pos = offset + 24 + 8 * nblocks
    chunks = []
    for cs in csizes:
        chunks.append(zlib.decompress(raw[pos:pos + int(cs)]))
        pos += int(cs)
    return np.frombuffer(b"".join(chunks), dtype=dtype)


def read_vtr_arrays(path, wanted):
    """Read vertex coordinates (x, y) and the named arrays from a .vtr file: CellData arrays
    come back as (ny, nx) cell arrays, PointData arrays (if not in CellData) as (ny+1, nx+1)."""
    with open(path, "rb") as f:
        data = f.read()
    start = data.find(b"<AppendedData")
    if start < 0:
        raise ValueError(f"{path}: no <AppendedData> section")
    raw = data[data.find(b"_", start) + 1:]  # appended block starts after the '_' marker
    header = data[:start].decode("utf-8", "replace")

    def find(section, name):
        m = re.search(rf"<{section}>(.*?)</{section}>", header, re.S)
        if m is None:
            return None
        for a in re.finditer(r"<DataArray ([^>]*)/>", m.group(1)):
            attrs = dict(re.findall(r'(\w+)="([^"]*)"', a.group(1)))
            if attrs.get("Name") == name:  # first match wins (names can repeat)
                return _decode_block(raw, int(attrs["offset"]), VTK_DTYPES[attrs["type"]])
        return None

    x, y = find("Coordinates", "x"), find("Coordinates", "y")
    out = {}
    for name in wanted:
        arr = find("CellData", name)
        if arr is not None:
            out[name] = arr.reshape(len(y) - 1, len(x) - 1)  # VTK order: x fastest
            continue
        arr = find("PointData", name)
        if arr is None:
            raise KeyError(f"{path}: array '{name}' not found in CellData or PointData")
        out[name] = arr.reshape(len(y), len(x))
    return x, y, out


def surface_from_density(xv, yv, rho, threshold):
    """Height of the first top-down crossing of `threshold` per column (NaN if none)."""
    yc = 0.5 * (yv[:-1] + yv[1:])
    order = np.argsort(yc)[::-1]  # top -> bottom
    yc, rho = yc[order], rho[order, :]
    h = np.full(rho.shape[1], np.nan)
    for i in range(rho.shape[1]):
        col = rho[:, i]
        hit = np.flatnonzero(col >= threshold)
        if hit.size == 0:
            continue
        j = hit[0]
        if j == 0:  # rock up to the top cell
            h[i] = yc[0]
            continue
        r_lo, r_hi = col[j], col[j - 1]  # r_hi is the (air-side) cell above
        f = (threshold - r_hi) / (r_lo - r_hi) if r_lo != r_hi else 0.0
        h[i] = yc[j - 1] + f * (yc[j] - yc[j - 1])
    return 0.5 * (xv[:-1] + xv[1:]), h


def gaussian_smooth(x, h, sigma):
    """Gaussian smoothing on a non-uniform x grid (same units as x); NaNs ignored."""
    ok = np.isfinite(h)
    out = np.full_like(h, np.nan)
    for i in np.flatnonzero(ok):
        w = np.exp(-0.5 * ((x[ok] - x[i]) / sigma) ** 2)
        out[i] = np.sum(w * h[ok]) / np.sum(w)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("vtk_dir", help="run vtk/ directory containing Subduction2D.pvd")
    ap.add_argument("--pvd", default="Subduction2D.pvd")
    ap.add_argument("--threshold", type=float, default=1000.0, help="density threshold [kg/m^3] (default 1000)")
    ap.add_argument("--smooth-km", type=float, default=0.0, help="optional Gaussian smoothing sigma [km] (default off)")
    ap.add_argument("--outdir", default="topography_density_csv")
    args = ap.parse_args()

    pvd = os.path.join(args.vtk_dir, args.pvd)
    if not os.path.isfile(pvd):
        sys.exit(f"not found: {pvd}")
    outdir = os.path.join(args.vtk_dir, args.outdir)
    os.makedirs(outdir, exist_ok=True)
    combined = os.path.join(outdir, "topography_timeseries.csv")

    entries = read_pvd(pvd)
    nwritten = 0
    with open(combined, "w", newline="") as fc:
        wc = csv.writer(fc)
        wc.writerow(["time_s", "time_kyr", "x_km", "h_km"])
        for t, fname in entries:
            path = os.path.join(args.vtk_dir, fname)
            if not os.path.isfile(path):
                print(f"  skip (missing): {fname}")
                continue
            xv, yv, arrs = read_vtr_arrays(path, ["density"])
            x, h = surface_from_density(xv, yv, arrs["density"], args.threshold)
            if args.smooth_km > 0:
                h = gaussian_smooth(x, h, args.smooth_km * 1e3)
            x_km, h_km = x / 1e3, h / 1e3
            stem = os.path.splitext(os.path.basename(fname))[0]
            with open(os.path.join(outdir, stem + ".csv"), "w", newline="") as fs:
                ws = csv.writer(fs)
                ws.writerow(["x_km", "h_km"])
                ws.writerows(zip(x_km, h_km))
            for xk, hk in zip(x_km, h_km):
                wc.writerow([t, t / SECONDS_PER_KYR, xk, hk])
            nwritten += 1
            print(f"  {fname}: t = {t / SECONDS_PER_KYR:8.1f} kyr, h range [{np.nanmin(h_km):.2f}, {np.nanmax(h_km):.2f}] km")

    print(f"wrote {nwritten} timesteps -> {combined}")


if __name__ == "__main__":
    main()

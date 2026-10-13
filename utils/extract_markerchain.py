#!/usr/bin/env python3
"""
extract_markerchain.py

Extract topography time-series data from a JustRelax marker-chain VTK output
(a directory of chain_NNNNNN.vtp files plus MarkerChain.pvd) into plain CSV,
for plotting outside ParaView.

Reads MarkerChain.pvd for the authoritative (simulation time -> file) mapping,
then each chain_*.vtp for its (x, h) points. Writes, into a subfolder inside
the given vtk/ directory:
  - one CSV per timestep: <stem>.csv with columns x_km, h_km
  - one combined long-format CSV: topography_timeseries.csv with columns
    time_s, time_kyr, x_km, h_km -- the file you want for an h(x,t) plot.

Usage:
    python3 extract_markerchain.py /path/to/run/vtk
    python3 extract_markerchain.py /path/to/run/vtk --outdir topo_csv

Requires the `vtk` Python package (not `pyvista`, which may be broken in
some environments on this system). No Julia dependencies touched.
"""

import argparse
import csv
import os
import sys
import xml.etree.ElementTree as ET

import vtk
from vtk.util.numpy_support import vtk_to_numpy

SECONDS_PER_KYR = 3600 * 24 * 365.25 * 1000


def read_pvd(pvd_path):
    """Return an ordered list of (time_seconds, filename) pairs, duplicates removed."""
    tree = ET.parse(pvd_path)
    entries = []
    seen = set()
    for dataset in tree.getroot().iter("DataSet"):
        fname = dataset.get("file")
        t = float(dataset.get("timestep"))
        key = (fname, t)
        if key in seen:
            continue
        seen.add(key)
        entries.append((t, fname))
    return entries


def read_chain_points(vtp_path):
    """Return the (x, h) points of a chain_*.vtp file as an Nx2 array, in km."""
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(vtp_path)
    reader.Update()
    pts = vtk_to_numpy(reader.GetOutput().GetPoints().GetData())
    return pts[:, 0] / 1.0e3, pts[:, 1] / 1.0e3  # x_km, h_km


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("vtk_dir", help="Directory containing MarkerChain.pvd and chain_*.vtp")
    parser.add_argument(
        "--pvd", default="MarkerChain.pvd",
        help="Name of the marker-chain .pvd file inside vtk_dir (default: MarkerChain.pvd)",
    )
    parser.add_argument(
        "--outdir", default="topography_csv",
        help="Subfolder name to create inside vtk_dir for the extracted CSVs (default: topography_csv)",
    )
    args = parser.parse_args()

    pvd_path = os.path.join(args.vtk_dir, args.pvd)
    if not os.path.isfile(pvd_path):
        sys.exit(f"Not found: {pvd_path}")

    outdir = os.path.join(args.vtk_dir, args.outdir)
    os.makedirs(outdir, exist_ok=True)

    entries = read_pvd(pvd_path)
    print(f"{len(entries)} timesteps found in {pvd_path}")

    timeseries_path = os.path.join(outdir, "topography_timeseries.csv")
    with open(timeseries_path, "w", newline="") as f_series:
        series_writer = csv.writer(f_series)
        series_writer.writerow(["time_s", "time_kyr", "x_km", "h_km"])

        for t, fname in entries:
            vtp_path = os.path.join(args.vtk_dir, fname)
            if not os.path.isfile(vtp_path):
                print(f"  skipping {fname}: file not found")
                continue

            x_km, h_km = read_chain_points(vtp_path)
            time_kyr = t / SECONDS_PER_KYR

            # Per-timestep CSV, for plotting a single profile or an animation.
            stem = os.path.splitext(fname)[0]
            with open(os.path.join(outdir, f"{stem}.csv"), "w", newline="") as f_step:
                step_writer = csv.writer(f_step)
                step_writer.writerow(["x_km", "h_km"])
                step_writer.writerows(zip(x_km, h_km))

            # Append to the combined long-format time series.
            for x, h in zip(x_km, h_km):
                series_writer.writerow([t, time_kyr, x, h])

    print(f"Wrote {len(entries)} per-timestep CSVs and {timeseries_path}")


if __name__ == "__main__":
    main()

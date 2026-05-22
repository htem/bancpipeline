#!/usr/bin/env python3
"""
Export BANC and maleCNS neuropil surfaces as watertight STL files for 3D printing.

Produces four files (by default in ~/Downloads):
  malecns_brain_neuropil.stl
  malecns_vnc_neuropil.stl
  banc_brain_neuropil.stl
  banc_vnc_neuropil.stl

Surfaces are loaded from the malecns and bancr R packages, voxelized,
hole-filled via morphological closing, reconstructed with marching cubes
(guaranteeing watertight output), Taubin-smoothed, and scaled so the
short axis = 6 inches (152.4 mm). Units in mm (STL convention).

Requirements:
  R packages: malecns, bancr, Rvcg, rgl
  Python:     trimesh, numpy, scipy, scikit-image

Usage:
  python banc/utilities/malecns-3dprint-stl.py
  python banc/utilities/malecns-3dprint-stl.py --target-inches 4
  python banc/utilities/malecns-3dprint-stl.py --voxels 150 --smooth 80
  python banc/utilities/malecns-3dprint-stl.py --dataset malecns   # malecns only
  python banc/utilities/malecns-3dprint-stl.py --dataset banc      # banc only
"""

import argparse
import os
import subprocess
import sys
import tempfile

import numpy as np
import trimesh
from scipy import ndimage
from skimage import measure


# malecns surfaces are in 8nm voxel coords; dividing by 8 gives nm-ish units
# bancr surfaces are in raw BANC coordinates (4nm voxel units)
SURFACES = {
    "malecns": {
        "r_package": "malecns",
        "surfaces": {
            "brain": "malecns.surf",
            "vnc": "malecnsvnc.surf",
        },
        "scale_divisor": "c(8, 8, 8)",
        "prefix": "malecns",
    },
    "banc": {
        "r_package": "bancr",
        "surfaces": {
            "brain": "banc_brain_neuropil.surf",
            "vnc": "banc_vnc_neuropil.surf",
        },
        "scale_divisor": "c(4, 4, 4)",
        "prefix": "banc",
    },
}

R_EXPORT_TEMPLATE = (
    'library({r_package})\n'
    'library(Rvcg)\n'
    'library(rgl)\n'
    'surfaces <- list({surface_list})\n'
    'for (name in names(surfaces)) {{\n'
    '  m <- as.mesh3d(surfaces[[name]] / {scale_divisor})\n'
    '  m$material <- NULL\n'
    '  outpath <- file.path("{tmpdir}", paste0(name, ".ply"))\n'
    '  Rvcg::vcgPlyWrite(m, filename = tools::file_path_sans_ext(outpath), binary = TRUE)\n'
    '  message(sprintf("Exported %s: %d faces -> %s", name, ncol(m$it), outpath))\n'
    '}}\n'
)


def repair_mesh(ply_path, name, target_voxels=100, smooth_iters=50):
    """Voxelize -> fill -> marching cubes -> smooth. Guarantees watertight."""
    m = trimesh.load(ply_path)
    bodies = m.split(only_watertight=False)
    bodies.sort(key=lambda b: b.faces.shape[0], reverse=True)
    m = bodies[0]
    print(f"  {name}: {m.faces.shape[0]} input faces")

    pitch = max(m.bounding_box.extents) / target_voxels
    vox = m.voxelized(pitch)
    matrix = ndimage.binary_fill_holes(vox.matrix.copy())
    matrix = ndimage.binary_closing(matrix, iterations=2)
    matrix = ndimage.binary_fill_holes(matrix)

    origin = vox.transform[:3, 3]
    verts, faces, normals, _ = measure.marching_cubes(
        matrix.astype(float), level=0.5, spacing=(pitch, pitch, pitch)
    )
    verts += origin

    watertight = trimesh.Trimesh(vertices=verts, faces=faces, process=True)
    parts = watertight.split(only_watertight=False)
    parts.sort(key=lambda p: p.faces.shape[0], reverse=True)
    watertight = parts[0]

    trimesh.smoothing.filter_taubin(watertight, iterations=smooth_iters)
    trimesh.repair.fix_normals(watertight)

    print(f"  {name}: {watertight.faces.shape[0]} faces, watertight={watertight.is_watertight}")
    return watertight


def scale_to_inches(mesh, target_inches):
    """Scale mesh so short axis = target_inches. Assumes input in raw units / 1000 = mm."""
    mesh.vertices /= 1000.0
    extents = mesh.bounding_box.extents
    mesh.vertices *= (target_inches * 25.4 / min(extents))
    return mesh


def export_dataset(dataset_key, args, tmpdir):
    """Export all surfaces for one dataset (malecns or banc)."""
    ds = SURFACES[dataset_key]
    surface_list = ", ".join(
        f"{k} = {v}" for k, v in ds["surfaces"].items()
    )
    r_script = R_EXPORT_TEMPLATE.format(
        r_package=ds["r_package"],
        surface_list=surface_list,
        scale_divisor=ds["scale_divisor"],
        tmpdir=tmpdir,
    )

    print(f"\nExporting {dataset_key} surfaces from R...")
    result = subprocess.run(
        ["Rscript", "-e", r_script],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        sys.exit(1)
    print(result.stderr.strip())

    for part in ds["surfaces"]:
        ply_path = os.path.join(tmpdir, f"{part}.ply")
        label = f"{dataset_key}_{part}"
        mesh = repair_mesh(
            ply_path, label,
            target_voxels=args.voxels,
            smooth_iters=args.smooth,
        )
        mesh = scale_to_inches(mesh, args.target_inches)

        ext = mesh.bounding_box.extents
        print(f"  {label}: {ext[0]:.1f} x {ext[1]:.1f} x {ext[2]:.1f} mm "
              f"({ext[0]/25.4:.1f} x {ext[1]/25.4:.1f} x {ext[2]/25.4:.1f} in)")

        out_name = f"{ds['prefix']}_{part}_neuropil.stl"
        out_path = os.path.join(args.outdir, out_name)
        mesh.export(out_path)
        sz = os.path.getsize(out_path) / 1024 / 1024
        print(f"  Saved: {out_path} ({mesh.faces.shape[0]} faces, {sz:.1f} MB)")

        raw = trimesh.load(ply_path)
        raw_parts = raw.split(only_watertight=False)
        raw_parts.sort(key=lambda b: b.faces.shape[0], reverse=True)
        raw = raw_parts[0]
        raw = scale_to_inches(raw, args.target_inches)
        raw_out = os.path.join(args.outdir, f"{ds['prefix']}_{part}_neuropil_original.stl")
        raw.export(raw_out)
        rsz = os.path.getsize(raw_out) / 1024 / 1024
        print(f"  Saved: {raw_out} ({raw.faces.shape[0]} faces, {rsz:.1f} MB)")


def main():
    parser = argparse.ArgumentParser(
        description="Export BANC and maleCNS neuropil surfaces as 3D-print STLs"
    )
    parser.add_argument("--target-inches", type=float, default=6.0,
                        help="Short-axis size in inches (default: 6)")
    parser.add_argument("--voxels", type=int, default=100,
                        help="Voxel grid resolution on long axis (default: 100)")
    parser.add_argument("--smooth", type=int, default=50,
                        help="Taubin smoothing iterations (default: 50)")
    parser.add_argument("--outdir", type=str, default=os.path.expanduser("~/Downloads"),
                        help="Output directory (default: ~/Downloads)")
    parser.add_argument("--dataset", type=str, default="all",
                        choices=["all", "malecns", "banc"],
                        help="Which dataset to export (default: all)")
    args = parser.parse_args()

    datasets = list(SURFACES.keys()) if args.dataset == "all" else [args.dataset]

    with tempfile.TemporaryDirectory() as tmpdir:
        for ds in datasets:
            export_dataset(ds, args, tmpdir)

    print("\nDone.")


if __name__ == "__main__":
    main()

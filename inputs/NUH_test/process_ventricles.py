"""
Split ventricles_open.vtu into LV_endo, LV_epi, and RV PLY files,
downsampled to UKBB resolution (~4-5K points per surface).

Pipeline:
  1. Extract surface and split into LV and RV by connected components
  2. Split LV into endo and epi using normal-based classification
  3. Downsample each surface to target vertex counts
  4. Save as PLY

Usage:
    python process_ventricles.py
    python process_ventricles.py --input ventricles_open.vtu --output-dir ./output
"""

import argparse
from pathlib import Path
import numpy as np
import pyvista as pv


def split_by_components(mesh: pv.PolyData):
    """Split surface into LV (larger) and RV (smaller) by connected components."""
    conn = mesh.connectivity()
    labels = conn.point_data['RegionId']
    unique, counts = np.unique(labels, return_counts=True)

    if len(unique) != 2:
        raise ValueError(f"Expected 2 connected components, got {len(unique)}")

    # Larger component is LV, smaller is RV
    lv_label = unique[np.argmax(counts)]
    rv_label = unique[np.argmin(counts)]

    lv_ids = np.where(conn.cell_data['RegionId'] == lv_label)[0]
    rv_ids = np.where(conn.cell_data['RegionId'] == rv_label)[0]

    lv = conn.extract_cells(lv_ids).extract_surface(algorithm=None).triangulate().clean()
    rv = conn.extract_cells(rv_ids).extract_surface(algorithm=None).triangulate().clean()

    print(f"LV: {lv.n_points} points, {lv.n_cells} cells")
    print(f"RV: {rv.n_points} points, {rv.n_cells} cells")
    return lv, rv


def split_lv_endo_epi(lv: pv.PolyData):
    """Split closed LV surface into epi and endo using normal dot product."""
    surf = lv.compute_normals(cell_normals=False, point_normals=True, consistent_normals=True)
    normals = surf.point_data['Normals']
    centroid = surf.points.mean(axis=0)

    # Classify points: outward normal = epi, inward = endo
    vecs = surf.points - centroid
    dots = np.sum(normals * vecs, axis=1)
    point_is_epi = dots > 0

    # Classify faces
    faces = surf.faces.reshape(-1, 4)[:, 1:4]
    v_epi = point_is_epi[faces]
    epi_count = v_epi.sum(axis=1)

    epi_ids = np.where(epi_count == 3)[0]
    endo_ids = np.where(epi_count == 0)[0]
    base_count = ((epi_count > 0) & (epi_count < 3)).sum()

    print(f"LV split: epi={len(epi_ids)}, endo={len(endo_ids)}, base={base_count} (discarded)")

    epi = surf.extract_cells(epi_ids).extract_surface(algorithm=None).triangulate()
    endo = surf.extract_cells(endo_ids).extract_surface(algorithm=None).triangulate()

    epi = epi.connectivity(extraction_mode='largest').clean()
    endo = endo.connectivity(extraction_mode='largest').clean()

    return epi, endo


def downsample(mesh: pv.PolyData, target_points: int):
    """Downsample mesh to approximately target_points vertices."""
    if mesh.n_points <= target_points:
        print(f"  Already at {mesh.n_points} points (target {target_points}), skipping")
        return mesh

    reduction = 1.0 - (target_points / mesh.n_points)
    result = mesh.decimate(reduction)
    result = result.clean()
    print(f"  Decimated: {mesh.n_points} -> {result.n_points} points (target {target_points})")
    return result


def process_ventricles(input_path: str, output_dir: str = None,
                       lv_endo_target: int = 4243,
                       lv_epi_target: int = 4315,
                       rv_target: int = 5230):
    input_path = Path(input_path)
    if output_dir is None:
        output_dir = input_path.parent
    else:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

    # Load and extract surface
    mesh = pv.read(str(input_path))
    if isinstance(mesh, pv.UnstructuredGrid):
        surf = mesh.extract_surface(algorithm=None).triangulate()
    else:
        surf = mesh.triangulate()

    print(f"Input: {surf.n_points} points, {surf.n_cells} cells")

    # Step 1: Split LV and RV
    lv, rv = split_by_components(surf)

    # Step 2: Split LV into endo and epi
    lv_epi, lv_endo = split_lv_endo_epi(lv)
    print(f"LV epi:  {lv_epi.n_points} points, {lv_epi.n_cells} cells")
    print(f"LV endo: {lv_endo.n_points} points, {lv_endo.n_cells} cells")

    # Step 3: Downsample to UKBB resolution
    print("\nDownsampling...")
    lv_endo = downsample(lv_endo, lv_endo_target)
    lv_epi = downsample(lv_epi, lv_epi_target)
    rv = downsample(rv, rv_target)

    # Step 4: Save
    paths = {
        'LV_endo_ex_1.ply': lv_endo,
        'LV_epi_ex_1.ply': lv_epi,
        'RV_ex_1.ply': rv,
    }
    for name, m in paths.items():
        out = output_dir / name
        m.save(str(out))
        print(f"Saved {name}: {m.n_points} points, {m.n_cells} cells")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split ventricles VTU into LV_endo, LV_epi, RV PLY files")
    parser.add_argument("--input", "-i",
                        default="C:/Users/CANG_DONGCHENG/InSilicoHeartGen/inputs/NUH_test/ventricles_open.vtu")
    parser.add_argument("--output-dir", "-o", default=None)
    parser.add_argument("--lv-endo-target", type=int, default=4243)
    parser.add_argument("--lv-epi-target", type=int, default=4315)
    parser.add_argument("--rv-target", type=int, default=5230)
    args = parser.parse_args()

    process_ventricles(args.input, args.output_dir,
                       args.lv_endo_target, args.lv_epi_target, args.rv_target)

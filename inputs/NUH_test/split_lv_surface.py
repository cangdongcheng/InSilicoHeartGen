"""
Split a closed LV surface mesh into epicardial and endocardial surfaces.

Usage:
    python split_lv_surface.py LV.vtk
    python split_lv_surface.py LV.vtk --output-dir ./output

Approach:
  1. Classify each POINT as epi or endo using the dot product of its
     normal with the vector from mesh centroid to the point.
  2. Classify each FACE by its vertices:
     - All 3 vertices epi   → epi face
     - All 3 vertices endo  → endo face
     - Mixed                → base transition face (discarded)
  3. This produces open surfaces with clean boundary loops at the base,
     which is what the downstream pipeline (add_lid_LV) expects.
"""

import argparse
from pathlib import Path
import numpy as np
import pyvista as pv


def split_lv_surface(input_path: str, output_dir: str = None):
    mesh = pv.read(input_path)
    input_path = Path(input_path)

    if output_dir is None:
        output_dir = input_path.parent
    else:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

    stem = input_path.stem

    # Ensure we have a PolyData surface
    if isinstance(mesh, pv.UnstructuredGrid):
        surf = mesh.extract_surface(algorithm=None)
    else:
        surf = mesh

    # Triangulate to ensure consistent face handling
    surf = surf.triangulate()

    # Compute point normals
    surf = surf.compute_normals(cell_normals=False, point_normals=True, consistent_normals=True)
    point_normals = surf.point_data['Normals']

    # Centroid of the mesh
    centroid = surf.points.mean(axis=0)

    # Classify each POINT as epi or endo
    vecs = surf.points - centroid
    dots = np.sum(point_normals * vecs, axis=1)
    point_is_epi = dots > 0  # True = epi, False = endo

    # Get face connectivity (triangles only)
    faces = surf.faces.reshape(-1, 4)[:, 1:4]  # each row: [v0, v1, v2]

    # Classify each face by its vertices
    v_epi = point_is_epi[faces]  # (n_cells, 3) boolean
    epi_count = v_epi.sum(axis=1)  # 0, 1, 2, or 3

    epi_face_mask = epi_count == 3   # all vertices are epi
    endo_face_mask = epi_count == 0  # all vertices are endo
    # mixed faces (1 or 2) are base transition → discarded

    print(f"Face classification: epi={epi_face_mask.sum()}, endo={endo_face_mask.sum()}, "
          f"base={((epi_count > 0) & (epi_count < 3)).sum()} (discarded)")

    # Extract surfaces
    epi_ids = np.where(epi_face_mask)[0]
    endo_ids = np.where(endo_face_mask)[0]

    epi_surf = surf.extract_cells(epi_ids).extract_surface(algorithm=None).triangulate()
    endo_surf = surf.extract_cells(endo_ids).extract_surface(algorithm=None).triangulate()

    # Keep only the largest connected component
    epi_surf = epi_surf.connectivity(extraction_mode='largest')
    endo_surf = endo_surf.connectivity(extraction_mode='largest')

    # Clean up: remove unused points
    epi_surf = epi_surf.clean()
    endo_surf = endo_surf.clean()

    # Save as PLY
    epi_path = output_dir / f"{stem}_epi.ply"
    endo_path = output_dir / f"{stem}_endo.ply"

    epi_surf.save(str(epi_path))
    endo_surf.save(str(endo_path))

    print(f"Epicardium:  {epi_surf.n_points} points, {epi_surf.n_cells} cells → {epi_path}")
    print(f"Endocardium: {endo_surf.n_points} points, {endo_surf.n_cells} cells → {endo_path}")

    return epi_surf, endo_surf


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split closed LV surface into epi and endo PLY files")
    parser.add_argument("input", help="Path to the input VTK surface mesh")
    parser.add_argument("--output-dir", "-o", default=None, help="Output directory (default: same as input)")
    args = parser.parse_args()

    split_lv_surface(args.input, args.output_dir)

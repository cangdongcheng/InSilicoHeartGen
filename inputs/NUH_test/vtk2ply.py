"""
Convert a VTK surface mesh to PLY format.

Usage:
    python vtk2ply.py RV_ex_1.vtk
    python vtk2ply.py RV_ex_1.vtk -o output.ply
"""

import argparse
from pathlib import Path
import pyvista as pv


def vtk2ply(input_path: str, output_path: str = None):
    mesh = pv.read(input_path)

    if isinstance(mesh, pv.UnstructuredGrid):
        surf = mesh.extract_surface(algorithm=None).triangulate().clean()
    else:
        surf = mesh.triangulate().clean()

    if output_path is None:
        output_path = str(Path(input_path).with_suffix('.ply'))

    surf.save(output_path)
    print(f"Saved: {surf.n_points} points, {surf.n_cells} cells → {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert VTK mesh to PLY")
    parser.add_argument("input", help="Input VTK file")
    parser.add_argument("-o", "--output", default=None, help="Output PLY path (default: same name with .ply)")
    args = parser.parse_args()

    vtk2ply(args.input, args.output)

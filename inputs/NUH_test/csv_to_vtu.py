"""
Compile InSilicoHeartGen CSV output files into a single VTU file.

Usage:
    python csv_to_vtu.py C:/Users/CANG_DONGCHENG/InSilicoHeartGen/outputs/NUH_test/1/ensi1/CSVFiles
    python csv_to_vtu.py <csv_dir> --output result.vtu
"""

import argparse
from pathlib import Path
import numpy as np
import pyvista as pv


def csv_to_vtu(csv_dir: str, output: str = None):
    csv_dir = Path(csv_dir)

    # Find the prefix (e.g. "1_coarse")
    xyz_files = list(csv_dir.glob("*_xyz.csv"))
    if not xyz_files:
        raise FileNotFoundError("No *_xyz.csv found in " + str(csv_dir))
    prefix = xyz_files[0].name.rsplit("_xyz.csv", 1)[0]
    print(f"Prefix: {prefix}")

    # Load mesh geometry
    points = np.loadtxt(csv_dir / f"{prefix}_xyz.csv", delimiter=",")
    cells = np.loadtxt(csv_dir / f"{prefix}_tetra.csv", delimiter=",", dtype=np.int64)
    # Convert from 1-based (MATLAB) to 0-based (VTK) indexing
    cells = cells - 1
    print(f"Points: {points.shape}, Cells: {cells.shape}")

    # Build unstructured grid (VTK tetra = cell type 10)
    n_cells = cells.shape[0]
    cell_array = np.hstack([np.full((n_cells, 1), 4, dtype=np.int64), cells]).ravel()
    cell_types = np.full(n_cells, 10, dtype=np.uint8)  # VTK_TETRA = 10
    grid = pv.UnstructuredGrid(cell_array, cell_types, points)

    # Cobiveco scalar node fields only
    for name in ["cobiveco-ab", "cobiveco-rt", "cobiveco-tm", "cobiveco-tv"]:
        fpath = csv_dir / f"{prefix}_nodefield_{name}.csv"
        if fpath.exists():
            data = np.loadtxt(fpath, delimiter=",")
            grid.point_data[name] = data
            print(f"  Added point field: {name} {data.shape}")

    # Output
    if output is None:
        output = str(csv_dir.parent / f"{prefix}_fields.vtu")

    grid.save(output)
    print(f"\nSaved: {output}")
    print(f"  {grid.n_points} points, {grid.n_cells} cells")
    print(f"  Point fields: {list(grid.point_data.keys())}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compile CSV fields into VTU")
    parser.add_argument("csv_dir", help="Path to CSVFiles directory")
    parser.add_argument("--output", "-o", default=None, help="Output VTU path")
    args = parser.parse_args()

    csv_to_vtu(args.csv_dir, args.output)

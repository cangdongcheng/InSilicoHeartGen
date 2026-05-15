"""
Generate glyph balls VTU at specific cobiveco landmark positions.
Reads the fields VTU, uses KDTree on cobiveco coordinates to find
the nearest mesh node for each target, then writes sphere glyphs.

Usage:
    python generate_glyph_balls.py <path_to_fields.vtu>
"""

import sys
import numpy as np
from scipy.spatial import KDTree
import pyvista as pv

# ---- Landmark targets ----
# Columns: ab, rt, tm, tv
targets = np.array([
    [0.62, 0.075, 1.0, 0],  # LV Mid-posterior superior
    [0.55, 0.175, 1.0, 0],  # LV Mid-posterior inferior
    [0.85, 0.570, 1.0, 0],  # LV Basal anterior paraseptal
    [0.40, 0.825, 1.0, 0],  # LV Mid-septal
    [0.45, 0.825, 1.0, 1],  # RV Mid-septal
    [0.80, 0.440, 1.0, 1],  # RV Free wall anterior
    [0.85, 0.250, 1.0, 1],  # RV Free wall posterior
])

labels = [
    "LV Mid-posterior superior",
    "LV Mid-posterior inferior",
    "LV Basal anterior paraseptal",
    "LV Mid-septal",
    "RV Mid-septal",
    "RV Free wall anterior",
    "RV Free wall posterior",
]

# Sphere radius (adjust to mesh scale, mm)
SPHERE_RADIUS = 2.0

def main():
    if len(sys.argv) < 2:
        print("Usage: python generate_glyph_balls.py <path_to_fields.vtu>")
        sys.exit(1)

    vtu_path = sys.argv[1]
    mesh = pv.read(vtu_path)
    print(f"Mesh: {mesh.n_points} nodes, {mesh.n_cells} cells")

    # Extract cobiveco array (columns: ab, rt, tm, tv)
    cobiveco = mesh.point_data["Cobiveco"]
    assert cobiveco.shape[1] == 4, f"Expected 4 columns, got {cobiveco.shape[1]}"
    print(f"Cobiveco range: {cobiveco.min(axis=0)} to {cobiveco.max(axis=0)}")

    # Build KDTree on cobiveco coordinates
    tree = KDTree(cobiveco)

    # Query nearest mesh node for each target
    distances, indices = tree.query(targets)

    # Get cartesian coordinates of matched nodes
    xyz = mesh.points[indices]

    print("\nMatched landmarks:")
    for i, (label, dist, idx) in enumerate(zip(labels, distances, indices)):
        cobi_actual = cobiveco[idx]
        print(f"  {label}")
        print(f"    target:  ab={targets[i,0]:.3f} rt={targets[i,1]:.3f} tm={targets[i,2]:.3f} tv={targets[i,3]:.3f}")
        print(f"    matched: ab={cobi_actual[0]:.3f} rt={cobi_actual[1]:.3f} tm={cobi_actual[2]:.3f} tv={cobi_actual[3]:.3f}")
        print(f"    xyz: [{xyz[i,0]:.2f}, {xyz[i,1]:.2f}, {xyz[i,2]:.2f}]  dist={dist:.4f}")

    # Build glyph spheres
    centers = pv.PolyData(xyz)
    centers["Label_ID"] = np.arange(len(labels))
    centers["Label"] = labels
    # TV=0 → LV, TV=1 → RV
    centers["Chamber"] = np.array(["LV" if t[3] == 0 else "RV" for t in targets])

    glyphs = centers.glyph(geom=pv.Sphere(radius=SPHERE_RADIUS), orient=False, scale=False)

    # Write output next to input
    out_path = vtu_path.replace(".vtu", "_glyphs.vtp")
    glyphs.save(out_path)
    print(f"\nWritten: {out_path}")

if __name__ == "__main__":
    main()

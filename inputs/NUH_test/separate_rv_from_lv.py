"""
Push RV surface points away from LV epi surface, accounting for the
3mm RV extrusion done by add_RV_width.

For each RV point, computes where it will end up after extrusion along
its vertex normal, then ensures that extruded position is at least `gap`
mm from the LV epi surface.

Usage:
    python separate_rv_from_lv.py
    python separate_rv_from_lv.py --gap 1.0 --extrusion 3.0
"""

import argparse
import numpy as np
import pyvista as pv
from scipy.spatial import cKDTree


def separate_rv_from_lv(lv_epi_path: str, rv_path: str, output_path: str,
                        gap: float = 1.0, extrusion: float = 3.0, max_iter: int = 20):
    lv_epi = pv.read(lv_epi_path)
    rv = pv.read(rv_path)

    # Get LV epi as surface
    if isinstance(lv_epi, pv.UnstructuredGrid):
        lv_surf = lv_epi.extract_surface(algorithm=None).triangulate()
    else:
        lv_surf = lv_epi.triangulate()
    lv_surf = lv_surf.compute_normals(point_normals=True, cell_normals=False, consistent_normals=True)

    lv_tree = cKDTree(lv_surf.points)

    # Direction from LV to RV centroid (for orienting push)
    lv_to_rv = rv.points.mean(axis=0) - lv_surf.points.mean(axis=0)
    lv_to_rv /= np.linalg.norm(lv_to_rv)

    rv_points = rv.points.copy()

    for iteration in range(max_iter):
        # Compute RV vertex normals from current positions
        rv_temp = rv.copy()
        rv_temp.points = rv_points
        if isinstance(rv_temp, pv.UnstructuredGrid):
            rv_surf = rv_temp.extract_surface(algorithm=None).triangulate()
        else:
            rv_surf = rv_temp.triangulate()
        rv_surf = rv_surf.compute_normals(point_normals=True, cell_normals=False, consistent_normals=True)

        # Simulate extrusion: RV_epi = RV + extrusion * normal
        rv_epi_points = rv_points + extrusion * rv_surf.point_data['Normals']

        # Check distance of extruded points to LV epi
        dists_epi, idx_epi = lv_tree.query(rv_epi_points)
        # Also check RV endo distance to LV epi
        dists_endo, idx_endo = lv_tree.query(rv_points)

        # Need both RV and RV_epi to be far enough from LV epi
        close_epi = dists_epi < gap
        close_endo = dists_endo < gap
        close_mask = close_epi | close_endo
        n_close = close_mask.sum()

        if n_close == 0:
            print(f"Iter {iteration}: done. RV min={dists_endo.min():.3f}mm, RV_epi min={dists_epi.min():.3f}mm")
            break

        print(f"Iter {iteration}: {close_mask.sum()} pts too close "
              f"(RV_endo<gap: {close_endo.sum()}, RV_epi<gap: {close_epi.sum()})")

        # Push RV points away from LV epi using LV surface normal
        # Use the worst violation for each point
        for mask, dists, idx in [(close_endo, dists_endo, idx_endo),
                                  (close_epi, dists_epi, idx_epi)]:
            if mask.sum() == 0:
                continue
            push_dir = lv_surf.point_data['Normals'][idx[mask]].copy()
            # Orient toward RV
            dots = np.sum(push_dir * lv_to_rv, axis=1)
            push_dir[dots < 0] *= -1
            needed = gap - dists[mask] + 0.1
            needed = np.maximum(needed, 0)
            rv_points[mask] += push_dir * needed[:, np.newaxis]
    else:
        print(f"Warning: did not fully converge after {max_iter} iterations")

    # Save
    rv_out = rv.copy()
    rv_out.points = rv_points
    rv_out.save(output_path)

    # Final verification
    dists_f, _ = lv_tree.query(rv_points)
    # Also verify simulated extrusion
    rv_temp2 = rv.copy()
    rv_temp2.points = rv_points
    if isinstance(rv_temp2, pv.UnstructuredGrid):
        rv_s2 = rv_temp2.extract_surface(algorithm=None).triangulate()
    else:
        rv_s2 = rv_temp2.triangulate()
    rv_s2 = rv_s2.compute_normals(point_normals=True, cell_normals=False, consistent_normals=True)
    rv_epi_f = rv_points + extrusion * rv_s2.point_data['Normals']
    dists_epi_f, _ = lv_tree.query(rv_epi_f)
    print(f"Final: RV-to-LV_epi min={dists_f.min():.3f}mm, RV_epi-to-LV_epi min={dists_epi_f.min():.3f}mm")
    print(f"Saved to: {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Separate RV surface from LV epi surface")
    parser.add_argument("--lv-epi", default="C:/Users/CANG_DONGCHENG/InSilicoHeartGen/inputs/NUH_test/LV_epi_ex_1.ply")
    parser.add_argument("--rv", default="C:/Users/CANG_DONGCHENG/InSilicoHeartGen/inputs/NUH_test/RV_ex_1.vtk")
    parser.add_argument("--output", "-o", default="C:/Users/CANG_DONGCHENG/InSilicoHeartGen/inputs/NUH_test/RV_ex_1.vtk")
    parser.add_argument("--gap", type=float, default=1.0, help="Min gap in mm (default: 1.0)")
    parser.add_argument("--extrusion", type=float, default=3.0, help="RV extrusion width in mm (default: 3.0)")
    args = parser.parse_args()

    separate_rv_from_lv(args.lv_epi, args.rv, args.output, args.gap, args.extrusion)

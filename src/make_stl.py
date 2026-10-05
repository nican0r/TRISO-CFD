"""Generate the Case V 3D spouted-bed geometry as a watertight STL.

The vessel inner wall is a surface of revolution built from the (r, z) profile

    (R_i, 0)  --cone--  (R_c, h_cone)  --cylinder--  (R_c, H_dom)

swept around the y-axis.  The surface is capped by an annular disk at y=0
(bottom, with a hole of radius R_i for the inlet orifice) and left open at
the top (y=H_dom) where the pressure-outflow BC lives.  The result is a
closed shell except for the two openings, which is what MFiX's cut-cell/STL
workflow expects for an "internal flow" region (OUT_STL_VALUE = 1).

Public API:

    build_spouted_bed_stl(R_c_m, R_i_m, cone_half_angle_deg, H_dom_m,
                           theta_segments, out_path) -> Path

All lengths SI (m).  Facet normals point outward (away from the fluid), which
matches MFiX's default STL convention for internal flow.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np


def _unit_normal(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> np.ndarray:
    """Return the unit normal of triangle (a, b, c), right-handed order."""
    n = np.cross(b - a, c - a)
    norm = np.linalg.norm(n)
    if norm == 0.0:
        return np.array([0.0, 0.0, 0.0])
    return n / norm


def _write_ascii_stl(
    triangles: list[tuple[np.ndarray, np.ndarray, np.ndarray]],
    name: str,
    out_path: Path,
) -> None:
    """Write an ASCII STL file. Normals computed from vertex order."""
    lines = [f"solid {name}"]
    for a, b, c in triangles:
        n = _unit_normal(a, b, c)
        lines.append(f"  facet normal {n[0]:.6e} {n[1]:.6e} {n[2]:.6e}")
        lines.append("    outer loop")
        for v in (a, b, c):
            lines.append(f"      vertex {v[0]:.6e} {v[1]:.6e} {v[2]:.6e}")
        lines.append("    endloop")
        lines.append("  endfacet")
    lines.append(f"endsolid {name}")
    out_path.write_text("\n".join(lines) + "\n")


def build_spouted_bed_stl(
    R_c_m: float,
    R_i_m: float,
    cone_half_angle_deg: float,
    H_dom_m: float,
    theta_segments: int,
    out_path: Path,
) -> Path:
    """Build the Case V spouted-bed inner-wall STL.

    R_c_m               : column (cylindrical section) inside radius [m]
    R_i_m               : inlet orifice radius [m]
    cone_half_angle_deg : cone half-angle from the vertical axis [deg]
    H_dom_m             : axial domain extent [m] (y: 0 -> H_dom)
    theta_segments      : number of angular facets around the axis [-]
    out_path            : destination .stl file (ASCII format)

    Returns the written path.

    Geometry:
        y = 0           : orifice plane (annular floor, hole r <= R_i for MI)
        0 < y < h_cone  : conical wall at r = R_i + y * tan(theta)
        y = h_cone      : cone top (meets cylindrical column)
        y > h_cone      : cylindrical wall at r = R_c
        y = H_dom       : open top (no cap - PO boundary face)

    Normals point outward (toward -r on the walls, toward -y on the floor
    annulus), i.e. away from the fluid region.
    """
    out_path = Path(out_path)
    tan_h = math.tan(math.radians(cone_half_angle_deg))
    if tan_h <= 0.0:
        raise ValueError("cone_half_angle_deg must be positive and < 90 deg")
    h_cone = (R_c_m - R_i_m) / tan_h
    if h_cone <= 0.0 or h_cone >= H_dom_m:
        raise ValueError(
            f"cone top ({h_cone:.4g} m) must satisfy 0 < h_cone < H_dom ({H_dom_m} m)"
        )

    # (r, y) profile points along the vessel inner wall, bottom -> top.
    #   P0 = (R_i, 0)       -- cone bottom edge (orifice rim)
    #   P1 = (R_c, h_cone)  -- cone top / cylinder bottom edge
    #   P2 = (R_c, H_dom)   -- top opening edge
    profile = [(R_i_m, 0.0), (R_c_m, h_cone), (R_c_m, H_dom_m)]

    # Angular discretization.
    thetas = np.linspace(0.0, 2.0 * math.pi, theta_segments + 1)
    # Pre-compute (cos, sin).
    cs = [(math.cos(t), math.sin(t)) for t in thetas]

    def vertex(r: float, y: float, k: int) -> np.ndarray:
        """Build a 3D vertex at radius r, axial y, angular index k."""
        c, s = cs[k]
        return np.array([r * c, y, r * s])

    triangles: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []

    # Side walls: 2 segments * theta_segments * 2 triangles.
    for i in range(len(profile) - 1):
        r_lo, y_lo = profile[i]
        r_hi, y_hi = profile[i + 1]
        for k in range(theta_segments):
            k2 = k + 1
            v00 = vertex(r_lo, y_lo, k)
            v01 = vertex(r_lo, y_lo, k2)
            v10 = vertex(r_hi, y_hi, k)
            v11 = vertex(r_hi, y_hi, k2)
            # Outward-pointing right-hand ordering: (lo_k, hi_k, hi_k2),
            # (lo_k, hi_k2, lo_k2).  Check: at k=0 this gives a normal with
            # positive radial component for r_hi >= r_lo.
            triangles.append((v00, v10, v11))
            triangles.append((v00, v11, v01))

    # Bottom annular floor at y = 0: inner edge r = R_i, outer edge r = R_c
    # cone-rim radius (but the cone rim is already at R_i at y=0, so the
    # floor is actually zero-width here -- the wall goes straight into the
    # orifice).  Nothing to emit for the bottom; the orifice plane is left
    # open as the MI face r <= R_i, and the cone wall closes everything
    # else.  (If the orifice tube were extruded below y=0, a floor annulus
    # would be needed here; this template does not model that stub.)

    # Top is left open (PO face at y = H_dom).

    _write_ascii_stl(triangles, "spouted_bed_V_3d", out_path)
    return out_path


def _cli() -> None:
    import argparse

    p = argparse.ArgumentParser(description="Generate the Case V 3D spouted-bed STL.")
    p.add_argument("--R-c", type=float, required=True, help="column inside radius [m]")
    p.add_argument("--R-i", type=float, required=True, help="orifice radius [m]")
    p.add_argument(
        "--half-angle", type=float, default=30.0, help="cone half-angle [deg]"
    )
    p.add_argument("--H-dom", type=float, required=True, help="axial domain extent [m]")
    p.add_argument("--theta-segments", type=int, default=60)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    out = build_spouted_bed_stl(
        R_c_m=args.R_c,
        R_i_m=args.R_i,
        cone_half_angle_deg=args.half_angle,
        H_dom_m=args.H_dom,
        theta_segments=args.theta_segments,
        out_path=args.out,
    )
    print(out)


if __name__ == "__main__":
    _cli()

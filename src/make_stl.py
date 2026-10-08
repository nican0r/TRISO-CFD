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

All lengths SI (m).  Facet normals point INWARD (toward the fluid region),
per MFiX's OUT_STL_VALUE = 1.0 internal-flow convention: the F_AT classifier
in cartesian_grid/intersect.f assigns F_AT < 0 (fluid) to nodes where
dot(vec_to_node, NORM_FACE) > 0, which requires NORM_FACE to point into the
fluid interior.  NETL's shipped tutorials/pic/spouted_bed_3d/geometry_0001.stl
has all sidewall normals with nr < 0 (inward toward the axis), confirming
this convention.
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
    close_bottom: bool = False,
    close_top: bool = False,
    L_stub_m: float = 0.0,
) -> Path:
    """Build the Case V spouted-bed inner-wall STL.

    R_c_m               : column (cylindrical section) inside radius [m]
    R_i_m               : inlet orifice radius [m]
    cone_half_angle_deg : cone half-angle from the vertical axis [deg]
    H_dom_m             : axial domain extent [m] (y: 0 -> H_dom)
    theta_segments      : number of angular facets around the axis [-]
    out_path            : destination .stl file (ASCII format)
    close_bottom        : if True, emit an annular floor at y=0 from
                          r=R_i to r=R_c (keeps the orifice open at r<R_i).
                          Normals point in +y (inward, toward the fluid above
                          the floor), per MFiX OUT_STL_VALUE=1.0 convention.
    close_top           : if True, emit a disk cap at y=H_dom from r=0
                          to r=R_c. Normals point in -y (inward, toward the
                          fluid below the cap).
    L_stub_m            : if > 0, prepend an inlet stub tube of radius R_i
                          from y = -L_stub to y = 0 to the profile.  Used
                          by the Day-10 fix to replicate NETL's PIC
                          spouted_bed_3d tutorial topology so the MI BC
                          can sit at the box bottom (y = -L_stub) and the
                          near-axis cells have STL sidewall facets within
                          1-2 cells (seeds F_AT propagation along J).

    Returns the written path.

    Geometry (side walls always; floor / cap only when requested):
        y = 0           : orifice plane
        0 < y < h_cone  : conical wall at r = R_i + y * tan(theta)
        y = h_cone      : cone top (meets cylindrical column)
        y > h_cone      : cylindrical wall at r = R_c
        y = H_dom       : top

    Normals point INWARD (toward +r on the walls -- toward the axis; toward
    +y on the floor annulus -- into the vessel; toward -y on the top cap --
    into the vessel), i.e. toward the fluid region, per MFiX OUT_STL_VALUE=1.0
    convention.

    Day-10 note.  The Day-8 run on the planar quadric slab used an open
    (side-walls-only) STL and only tested STL writing; MFiX never consumed
    it.  On Day 10 the STL IS consumed by MFiX cut-cell preprocessing,
    which runs an inside/outside classifier (``F_AT`` in
    cartesian_grid/intersect.f) that leaves cells with no local facet
    neighbour as UNDEFINED and marks them BLOCKED.  An open-top /
    open-floor shell therefore mis-classifies the near-axis cells inside
    the vessel.  ``close_bottom=True`` + ``close_top=True`` make the STL
    watertight (except for the orifice hole), restoring the correct
    fluid-cell count.
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
    # When L_stub_m > 0, prepend an inlet stub-tube profile point at
    # (R_i, -L_stub) so the stub cylinder sidewalls are swept together
    # with the cone + column.  The stub is open at y = -L_stub (MI BC
    # sits at the Cartesian box bottom there; MFiX handles the "closure"
    # from a BC perspective).
    if L_stub_m > 0.0:
        profile = [
            (R_i_m, -L_stub_m),
            (R_i_m, 0.0),
            (R_c_m, h_cone),
            (R_c_m, H_dom_m),
        ]
    else:
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
            # Inward-pointing right-hand ordering: (lo_k, hi_k2, hi_k),
            # (lo_k, lo_k2, hi_k2).  Reversing the second and third vertices
            # flips the cross-product, giving nr < 0 (toward the axis) on the
            # cylinder and nr < 0 (toward the interior) on the cone.
            triangles.append((v00, v11, v10))   # inward normal
            triangles.append((v00, v01, v11))   # inward normal

    # Bottom annular floor at y = 0: a ring from r = R_i (orifice rim) to
    # r = R_c (vessel wall).  The wall at y = 0 has r = R_i (cone bottom
    # edge) so the floor annulus spans r in [R_i, R_c].  Facet normals
    # point in +y (inward, toward the fluid above the floor), per MFiX
    # OUT_STL_VALUE=1.0 convention.
    # Only emitted when close_bottom=True (Day 10 run with closed STL).
    if close_bottom:
        for k in range(theta_segments):
            k2 = k + 1
            # Inner ring edge (r=R_i), outer ring edge (r=R_c), at y=0.
            vi0 = vertex(R_i_m, 0.0, k)
            vi1 = vertex(R_i_m, 0.0, k2)
            vo0 = vertex(R_c_m, 0.0, k)
            vo1 = vertex(R_c_m, 0.0, k2)
            # Right-handed ordering for +y normal (counter-clockwise when
            # viewed from above, i.e. inward toward the fluid above the
            # floor): (inner_k2, outer_k, inner_k), (inner_k2, outer_k2,
            # outer_k).  Verified by cross product: n_y > 0.
            triangles.append((vi1, vo0, vi0))   # inward (+y) normal
            triangles.append((vi1, vo1, vo0))   # inward (+y) normal

    # Top is left open (PO face at y = H_dom covers the full cross-section;
    # a disk cap would collide with the PO BC region).  close_top is
    # accepted for API symmetry but kept False by callers.
    if close_top:
        v_axis = np.array([0.0, H_dom_m, 0.0])
        for k in range(theta_segments):
            k2 = k + 1
            vo0 = vertex(R_c_m, H_dom_m, k)
            vo1 = vertex(R_c_m, H_dom_m, k2)
            # Right-handed ordering for -y normal (inward, toward the fluid
            # below the cap): (axis, outer_k, outer_k2).  Cross product
            # gives n_y = -(R_c^2 * sin(dtheta)) < 0 (inward).
            triangles.append((v_axis, vo0, vo1))   # inward (-y) normal

    _write_ascii_stl(triangles, "spouted_bed_V_3d", out_path)
    return out_path


def build_spouted_bed_stl_implicit(
    R_c_m: float,
    R_i_m: float,
    cone_half_angle_deg: float,
    H_dom_m: float,
    L_stub_m: float,
    dx_sample_m: float,
    out_path: Path,
    top_extend_m: float = 0.0,
) -> Path:
    """Build the Case V spouted-bed STL via implicit-boolean sampling.

    Samples a signed-distance function (SDF) on a dense Cartesian grid,
    extracts the zero isosurface with marching cubes (via PyVista), then
    rewrites the triangles through ``_write_ascii_stl`` so the output has
    MFiX-standard ASCII STL formatting and outward-pointing normals.

    All inputs SI:

    R_c_m               : column (cylindrical section) inside radius [m]
    R_i_m               : inlet orifice radius [m]
    cone_half_angle_deg : cone half-angle from the vertical axis [deg]
    H_dom_m             : axial domain extent above the orifice [m]
                          (vessel runs y: 0 -> H_dom)
    L_stub_m            : inlet stub-tube length below the orifice [m]
                          (stub runs y: -L_stub -> 0 at radius R_i).  Must
                          be > 0 for Day-10 fix 3 (closes the SDF below the
                          orifice so marching cubes produces a watertight
                          solid; MFiX cuts the MI through the bottom cap).
    dx_sample_m         : SDF sampling grid cell size [m].  Should be
                          <= 0.5 * (smallest solver cell) so every solver
                          cell edge intersects at least one STL facet
                          (notes/day10_summary.md fix 3).
    out_path            : destination .stl file (ASCII format)
    top_extend_m        : if > 0, extends the freeboard cylinder's upper
                          cap above y = H_dom by this amount so the PO
                          plane at y = H_dom cuts through STL side walls
                          (not through the top cap facets, which MFiX
                          discards as "flush with box face").  5 mm is
                          sufficient at typical mesh resolutions.
                          Default 0.0 keeps backward compatibility.

    Returns the written path.

    SDF primitives (all axisymmetric around the y-axis):
      - inlet stub tube : cylinder r = R_i, y in [-L_stub, 0]
      - cone frustum    : r = R_i at y = 0 up to r = R_c at y = h_cone
                          where h_cone = (R_c - R_i) / tan(theta)
      - freeboard       : cylinder r = R_c, y in [h_cone, H_dom]
    Each primitive's SDF is the max of (r - R_bound, y_bot - y, y - y_top),
    i.e. negative inside the solid, zero on the boundary, positive outside.
    The vessel SDF is the min of the three (union).  Marching cubes at the
    zero isosurface returns the watertight closed boundary, including a
    bottom-disk cap at y = -L_stub and a top-disk cap at y = H_dom.  MFiX
    cuts the MI/PO BCs through these caps via its BC-region logic.
    """
    out_path = Path(out_path)
    tan_h = math.tan(math.radians(cone_half_angle_deg))
    if tan_h <= 0.0:
        raise ValueError("cone_half_angle_deg must be positive and < 90 deg")
    h_cone_m = (R_c_m - R_i_m) / tan_h
    if h_cone_m <= 0.0 or h_cone_m >= H_dom_m:
        raise ValueError(
            f"cone top ({h_cone_m:.4g} m) must satisfy 0 < h_cone < H_dom "
            f"({H_dom_m} m)"
        )
    if L_stub_m <= 0.0:
        raise ValueError("L_stub_m must be > 0 for the implicit SDF path")
    if dx_sample_m <= 0.0:
        raise ValueError("dx_sample_m must be positive")

    # Lazy import so pytest collection doesn't require pyvista at import time
    # of the module (keeps CLI/Day-8 users decoupled from the marching-cubes
    # dependency).  pyvista is already a project dependency used by
    # src/post/read_vtk.py so this is a free call.
    import pyvista as pv

    # Sampling grid.  Margin of 2 * dx on all sides so the surface closes
    # on the box faces instead of leaving dangling edges where the solid
    # meets the sampling-grid boundary.  Day-10 R1: optionally extend the
    # freeboard upper bound by top_extend_m so the top cap sits above the
    # Cartesian box face (y_max_m = H_dom_m) and the PO plane cuts through
    # side-wall facets.
    margin = 2.0 * dx_sample_m
    y_top_m = H_dom_m + top_extend_m
    x0 = -R_c_m - margin
    x1 = R_c_m + margin
    y0 = -L_stub_m - margin
    y1 = y_top_m + margin
    z0 = -R_c_m - margin
    z1 = R_c_m + margin

    # Point counts (not cell counts): ImageData.dimensions is number of
    # grid points along each axis.  Use ceil so no face is skipped.
    nx = int(math.ceil((x1 - x0) / dx_sample_m)) + 1
    ny = int(math.ceil((y1 - y0) / dx_sample_m)) + 1
    nz = int(math.ceil((z1 - z0) / dx_sample_m)) + 1

    # Build sample point coordinates per axis.
    xs = x0 + dx_sample_m * np.arange(nx)
    ys = y0 + dx_sample_m * np.arange(ny)
    zs = z0 + dx_sample_m * np.arange(nz)

    # Mesh them: use ij indexing so X[i,j,k] etc.
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    R = np.sqrt(X * X + Z * Z)

    # Primitive SDFs (negative inside, positive outside).
    # Stub tube: cylinder radius R_i, y in [-L_stub, 0].
    sdf_stub = np.maximum.reduce(
        [R - R_i_m, -L_stub_m - Y, Y - 0.0]
    )

    # Cone frustum: r_bound = R_i + y * tan(theta); y in [0, h_cone].
    r_cone_bound = R_i_m + Y * tan_h
    sdf_cone = np.maximum.reduce(
        [R - r_cone_bound, 0.0 - Y, Y - h_cone_m]
    )

    # Freeboard cylinder: radius R_c, y in [h_cone, H_dom + top_extend].
    sdf_free = np.maximum.reduce(
        [R - R_c_m, h_cone_m - Y, Y - y_top_m]
    )

    # Union: min of all primitive SDFs.
    sdf = np.minimum.reduce([sdf_stub, sdf_cone, sdf_free]).astype(np.float32)

    # PyVista ImageData flattening order is Fortran-order (first axis varies
    # fastest when written to point_data).
    grid = pv.ImageData(
        dimensions=(nx, ny, nz),
        spacing=(dx_sample_m, dx_sample_m, dx_sample_m),
        origin=(x0, y0, z0),
    )
    grid.point_data["sdf"] = sdf.flatten(order="F")

    surface = grid.contour(isosurfaces=[0.0], scalars="sdf")
    surface = surface.triangulate()

    # Extract triangles and rewrite via the project's ASCII STL helper so
    # the file matches MFiX's expected header/footer and so normals are
    # recomputed from vertex order.  PyVista returns a PolyData whose
    # faces have mixed polygon sizes; after triangulate() they're all 3-gons.
    faces = surface.regular_faces  # (n_tri, 3) int ndarray
    pts = np.asarray(surface.points)
    if faces is None or len(faces) == 0:
        raise RuntimeError(
            "marching cubes produced no triangles; check SDF sign convention"
        )

    # Determine which winding gives outward-pointing normals (toward +sdf,
    # i.e. away from the vessel interior).  Sample the SDF at (centroid +
    # eps * computed_normal) for one triangle; if that sample is negative,
    # the normal points inward and we flip all triangles.
    centroids = (pts[faces[:, 0]] + pts[faces[:, 1]] + pts[faces[:, 2]]) / 3.0
    e1 = pts[faces[:, 1]] - pts[faces[:, 0]]
    e2 = pts[faces[:, 2]] - pts[faces[:, 0]]
    tri_normals = np.cross(e1, e2)
    norms = np.linalg.norm(tri_normals, axis=1)
    norms[norms == 0.0] = 1.0
    tri_normals /= norms[:, None]

    probe_eps = 0.25 * dx_sample_m
    probe = centroids + probe_eps * tri_normals
    # Evaluate the analytic union-SDF at the probe points.
    r_probe = np.sqrt(probe[:, 0] ** 2 + probe[:, 2] ** 2)
    yp = probe[:, 1]
    def _sdf_point(rp, yp):
        s_stub = np.maximum.reduce([rp - R_i_m, -L_stub_m - yp, yp - 0.0])
        rc = R_i_m + yp * tan_h
        s_cone = np.maximum.reduce([rp - rc, 0.0 - yp, yp - h_cone_m])
        s_free = np.maximum.reduce([rp - R_c_m, h_cone_m - yp, yp - y_top_m])
        return np.minimum.reduce([s_stub, s_cone, s_free])
    sdf_at_probe = _sdf_point(r_probe, yp)
    # Majority vote for INWARD normals per MFiX OUT_STL_VALUE=1.0 convention.
    # A probe at (centroid + eps*normal) landing in sdf < 0 (inside the fluid
    # cavity) means the current normal points INTO the fluid -- which is what
    # MFiX requires (verified in cartesian_grid/intersect.f F_AT sign test,
    # and against NETL tutorials/pic/spouted_bed_3d/geometry_0001.stl which
    # has all sidewall nr < 0).  Flip only when the majority points OUTWARD
    # (probes land in sdf > 0, inside_frac < 0.5).
    inside_frac = float(np.mean(sdf_at_probe < 0.0))
    flip = inside_frac < 0.5

    triangles: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
    for f in faces:
        a, b, c = pts[f[0]], pts[f[1]], pts[f[2]]
        if flip:
            triangles.append((a, c, b))
        else:
            triangles.append((a, b, c))

    _write_ascii_stl(triangles, "spouted_bed_V_3d_implicit", out_path)
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

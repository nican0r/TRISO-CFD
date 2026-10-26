"""Generate the Case V and Case C26 3D spouted-bed geometries as STLs.

The vessel inner wall is a surface of revolution built from the (r, z) profile

    (R_i, 0)  --cone--  (R_c, h_cone)  --cylinder--  (R_c, H_dom)

swept around the y-axis.  The surface is capped by an annular disk at y=0
(bottom, with a hole of radius R_i for the inlet orifice) and left open at
the top (y=H_dom) where the pressure-outflow BC lives.  The result is a
closed shell except for the two openings, which is what MFiX's cut-cell/STL
workflow expects for an "internal flow" region (OUT_STL_VALUE = 1).

Public API:

    build_spouted_bed_stl_implicit(R_c_m, R_i_m, cone_half_angle_deg, H_dom_m,
                                   L_stub_m, dx_sample_m, out_path, ...) -> Path
        Case V: stub + cone + column via SDF sampling and marching cubes.
    build_cone_column_stl_implicit(R_c_m, R_t_m, cone_half_angle_deg, H_dom_m,
                                   dx_sample_m, out_path, top_extend_m) -> Path
        Case C26: closed cone + column for a gas point source (no inlet stub).

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


# Minimum ||cross(V2-V1, V3-V1)|| magnitude required for a triangle to be
# considered non-degenerate.  MFiX's IS_POINT_INSIDE_FACET
# (get_stl_data.f:2386) computes Inv_denom = 1 / (|V0|^2 * |V1|^2 - (V0.V1)^2)
# = 1 / ||cross||^2 with no guard.  Marching cubes on an SDF whose zero-level
# touches grid points emits "sliver" facets whose two edges are coincident or
# collinear.  Those facets pass MFiX's STL filter when the written normal has
# unit length (|cross| normalized away), then trigger SIGFPE in the barycentric
# check.  Writing [0,0,0] for crosses below this threshold lets MFiX's
# NORM > TOL_STL = 1e-6 filter drop them cleanly.
_CROSS_TOL_M2 = 1.0e-12


def _unit_normal(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> np.ndarray:
    """Return the unit normal of triangle (a, b, c), right-handed order.

    Returns ``[0, 0, 0]`` for degenerate triangles (|cross| < _CROSS_TOL_M2)
    so downstream STL consumers (MFiX) can filter them via the zero-normal check.
    """
    n = np.cross(b - a, c - a)
    norm = np.linalg.norm(n)
    if norm < _CROSS_TOL_M2:
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


def build_spouted_bed_stl_implicit(
    R_c_m: float,
    R_i_m: float,
    cone_half_angle_deg: float,
    H_dom_m: float,
    L_stub_m: float,
    dx_sample_m: float,
    out_path: Path,
    top_extend_m: float = 0.0,
    R_stub_m: float | None = None,
    taper_length_m: float = 0.0,
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
    R_stub_m            : inlet-stub cylinder radius [m].  If None (the
                          default), R_stub = R_i (stub tube narrows to
                          match the orifice).  Day-10 R2 Option 4 passes
                          R_stub > R_i so the stub-tube cross-section has
                          multiple fluid cells on the MI face; this
                          escapes MFiX 26.1.2's F_AT cut-cell classifier
                          failure mode when R_i / dx <= 1.  The cone
                          frustum still narrows from R_c (at y = h_cone)
                          down to the physical R_i (at y = 0), so the
                          orifice remains the paper's dimension; a
                          horizontal shelf (ring annulus) connects the
                          stub top (r = R_stub, y = 0-) to the cone apex
                          (r = R_i, y = 0+).  Callers should
                          correspondingly scale bc_v_g by (R_i/R_stub)**2
                          to preserve the orifice jet velocity.
    taper_length_m      : if > 0, replaces the hard shelf at y = 0 with a
                          linear taper of this axial length [m].  The
                          vessel inner-wall radius narrows continuously
                          from R_stub (at y = -taper_length_m) down to
                          R_i (at y = 0).  Day-10 R2 Option 2 (2026-10-09):
                          the Option-4 hard shelf triggered a shock on
                          AWS (Vg >> 100 m/s in a stub-interior cell) that
                          stalled the solver; a 1-cell taper smooths the
                          area contraction and eliminates the shock.  Must
                          satisfy taper_length_m <= L_stub_m.  Default 0.0
                          preserves bit-identical behaviour with the
                          Day-10 R1/R2-Option-4 hard-shelf build.

    Returns the written path.

    SDF primitives (all axisymmetric around the y-axis):
      - inlet stub tube : cylinder r = R_stub, y in [-L_stub, 0]
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
    if R_stub_m is None:
        R_stub_m = R_i_m
    if R_stub_m < R_i_m:
        raise ValueError(
            f"R_stub_m ({R_stub_m}) must be >= R_i_m ({R_i_m}); a stub tube "
            "narrower than the orifice is non-physical."
        )
    if taper_length_m < 0.0:
        raise ValueError(
            f"taper_length_m ({taper_length_m}) must be >= 0."
        )
    if taper_length_m > L_stub_m:
        raise ValueError(
            f"taper_length_m ({taper_length_m}) must be <= L_stub_m "
            f"({L_stub_m}); the taper must fit inside the stub region."
        )

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
    # Note R_c_m is always the widest radius (freeboard > cone top > stub);
    # the sampling box only needs to clear it in x,z.
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
    # Stub tube: cylinder radius R_stub, y in [-L_stub, 0].  R_stub >= R_i
    # (Day-10 R2 Option 4); when R_stub == R_i the stub tapers to the
    # orifice and the behaviour is identical to the Day-10 R1 build.
    # Day-10 R2 Option 2: when taper_length_m > 0, the stub inner wall
    # linearly narrows from R_stub (at y = -taper_length_m) down to R_i
    # (at y = 0) so the SDF matches a converging-diverging nozzle instead
    # of a hard shelf.  taper_length_m == 0 recovers the hard-shelf build.
    if taper_length_m > 0.0:
        # r_bound(y) = R_stub for y <= -taper, then linearly down to R_i at y=0.
        r_bound_stub = np.where(
            Y >= -taper_length_m,
            R_stub_m + (R_i_m - R_stub_m) * (Y + taper_length_m) / taper_length_m,
            R_stub_m,
        )
    else:
        r_bound_stub = R_stub_m
    sdf_stub = np.maximum.reduce(
        [R - r_bound_stub, -L_stub_m - Y, Y - 0.0]
    )

    # Cone frustum: r_bound = R_i + y * tan(theta), clamped to R_c.
    # Day-10 R5: the cone's caps OVERLAP its neighbours by seam_m instead of
    # sharing the planes y = 0 and y = h_cone.  With shared caps, the union
    # SDF is exactly 0 across the whole orifice disc whenever a sampling-grid
    # plane lands on y = 0 (it does: ys = -L_stub - 2*dx + k*dx), and marching
    # cubes emitted a membrane there that MFiX treated as a wall (R4 STL: 144
    # facets across r = 0.24-1.84 mm; see notes/day10_r5_orifice_clog.md).
    # Below y = 0 the cone radius is < R_i <= R_stub and above h_cone it is
    # clamped to R_c, so the overlap moves no wall.
    seam_m = dx_sample_m
    r_cone_bound = np.minimum(R_i_m + Y * tan_h, R_c_m)
    sdf_cone = np.maximum.reduce(
        [R - r_cone_bound, -seam_m - Y, Y - (h_cone_m + seam_m)]
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

    # Drop degenerate facets (coincident or near-collinear vertices) that
    # marching cubes emits whenever the SDF has values of exactly 0 on a
    # grid point.  MFiX's STL reader (get_stl_data.f:973-992) only ignores
    # facets whose written normal has magnitude below TOL_STL=1e-6; a sliver
    # with |cross|=1e-21 still gets normalized to a unit vector and passes
    # that filter, then triggers divide-by-zero in IS_POINT_INSIDE_FACET
    # (get_stl_data.f:2386, Inv_denom = 1 / ||cross||^2).  Filtering here
    # removes the sliver entirely so the written STL is safe for any
    # downstream STL consumer, not just MFiX.
    e1 = pts[faces[:, 1]] - pts[faces[:, 0]]
    e2 = pts[faces[:, 2]] - pts[faces[:, 0]]
    tri_normals = np.cross(e1, e2)
    norms = np.linalg.norm(tri_normals, axis=1)
    valid_mask = norms >= _CROSS_TOL_M2
    faces = faces[valid_mask]
    tri_normals = tri_normals[valid_mask]
    norms = norms[valid_mask]
    if faces.shape[0] == 0:
        raise RuntimeError(
            "all marching-cubes faces were dropped as degenerate; check SDF"
        )

    # Determine which winding gives outward-pointing normals (toward +sdf,
    # i.e. away from the vessel interior).  Sample the SDF at (centroid +
    # eps * computed_normal) for one triangle; if that sample is negative,
    # the normal points inward and we flip all triangles.
    centroids = (pts[faces[:, 0]] + pts[faces[:, 1]] + pts[faces[:, 2]]) / 3.0
    tri_normals /= norms[:, None]

    probe_eps = 0.25 * dx_sample_m
    probe = centroids + probe_eps * tri_normals
    # Evaluate the analytic union-SDF at the probe points.
    r_probe = np.sqrt(probe[:, 0] ** 2 + probe[:, 2] ** 2)
    yp = probe[:, 1]
    def _sdf_point(rp, yp):
        # Day-10 R2 Option 2 taper: match the gridded SDF above so the
        # inward/outward probe sees the same geometry as marching cubes.
        if taper_length_m > 0.0:
            rb_stub = np.where(
                yp >= -taper_length_m,
                R_stub_m + (R_i_m - R_stub_m) * (yp + taper_length_m) / taper_length_m,
                R_stub_m,
            )
        else:
            rb_stub = R_stub_m
        s_stub = np.maximum.reduce([rp - rb_stub, -L_stub_m - yp, yp - 0.0])
        rc = np.minimum(R_i_m + yp * tan_h, R_c_m)
        s_cone = np.maximum.reduce([rp - rc, -seam_m - yp, yp - (h_cone_m + seam_m)])
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


def build_cone_column_stl_implicit(
    R_c_m: float,
    R_t_m: float,
    cone_half_angle_deg: float,
    H_dom_m: float,
    dx_sample_m: float,
    out_path: Path,
    top_extend_m: float = 0.0,
) -> Path:
    """Closed-bottom cone + column STL for gas-point-source decks (Case C26).

    Same implicit-SDF + marching-cubes method as
    :func:`build_spouted_bed_stl_implicit`, but with **no inlet stub**.  The
    cone narrows from R_c down to a throat of radius R_t at y = 0, where a
    flat cap closes the vessel.  Gas enters through an MFiX point source
    (``ps_*`` keywords) inside the vessel, not through a boundary, so there is
    no opening to cut a mass inflow through.  See notes/c26_conry2026_case.md.

    All inputs SI:

    R_c_m               : column (cylindrical section) inside radius [m]
    R_t_m               : cone throat radius at the closed bottom (y = 0) [m]
    cone_half_angle_deg : cone half-angle from the vertical axis [deg]
    H_dom_m             : vessel top / pressure-outflow plane [m]
    dx_sample_m         : SDF sampling step [m]
    out_path            : destination ASCII .stl
    top_extend_m        : push the top cap this far above y = H_dom so the PO
                          plane cuts side-wall facets (spouted-bed recipe
                          Rule 19).

    Primitives (axisymmetric about the y-axis, negative inside):
      - cone   : r <= R_t + y tan(theta), y in [0, h_cone + seam]
      - column : r <= R_c,                y in [h_cone, H_dom + top_extend]
    with h_cone = (R_c - R_t) / tan(theta).  The cone's TOP cap overlaps the
    column by ``seam = dx_sample`` (the Day-10 R5 fix): with a shared plane the
    union SDF is exactly zero across the disc whenever a sampling plane lands
    on it, and marching cubes emits a spurious internal membrane.  The cone's
    BOTTOM cap at y = 0 is a real wall and is kept.

    Normals point INWARD (MFiX OUT_STL_VALUE = 1.0 internal-flow convention),
    chosen by the same majority-vote probe as the V_3d builder.
    """
    out_path = Path(out_path)
    tan_h = math.tan(math.radians(cone_half_angle_deg))
    if not (0.0 < cone_half_angle_deg < 90.0):
        raise ValueError("cone_half_angle_deg must be in (0, 90) deg")
    if R_t_m <= 0.0:
        raise ValueError(f"R_t_m ({R_t_m}) must be > 0")
    if R_t_m >= R_c_m:
        raise ValueError(f"R_t_m ({R_t_m}) must be < R_c_m ({R_c_m})")
    h_cone_m = (R_c_m - R_t_m) / tan_h
    if h_cone_m >= H_dom_m:
        raise ValueError(
            f"cone top ({h_cone_m:.4g} m) must lie below H_dom ({H_dom_m} m)"
        )
    if dx_sample_m <= 0.0:
        raise ValueError("dx_sample_m must be positive")
    if dx_sample_m >= R_t_m:
        raise ValueError(
            f"dx_sample_m ({dx_sample_m}) must be < R_t_m ({R_t_m}) or the "
            "throat is not resolved by the sampling grid"
        )
    if top_extend_m < 0.0:
        raise ValueError("top_extend_m must be >= 0")

    import pyvista as pv

    seam_m = dx_sample_m
    y_top_m = H_dom_m + top_extend_m

    def sdf_fn(rr: np.ndarray, yy: np.ndarray) -> np.ndarray:
        r_cone = np.minimum(R_t_m + yy * tan_h, R_c_m)
        s_cone = np.maximum.reduce(
            [rr - r_cone, 0.0 - yy, yy - (h_cone_m + seam_m)]
        )
        s_col = np.maximum.reduce([rr - R_c_m, h_cone_m - yy, yy - y_top_m])
        return np.minimum(s_cone, s_col)

    margin = 2.0 * dx_sample_m
    x0, x1 = -R_c_m - margin, R_c_m + margin
    y0, y1 = -margin, y_top_m + margin
    nx = int(math.ceil((x1 - x0) / dx_sample_m)) + 1
    ny = int(math.ceil((y1 - y0) / dx_sample_m)) + 1
    xs = x0 + dx_sample_m * np.arange(nx)
    ys = y0 + dx_sample_m * np.arange(ny)
    zs = xs.copy()
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    sdf = sdf_fn(np.sqrt(X * X + Z * Z), Y).astype(np.float32)

    grid = pv.ImageData(
        dimensions=(nx, ny, nx),
        spacing=(dx_sample_m, dx_sample_m, dx_sample_m),
        origin=(x0, y0, x0),
    )
    grid.point_data["sdf"] = sdf.flatten(order="F")
    surface = grid.contour(isosurfaces=[0.0], scalars="sdf").triangulate()
    faces = surface.regular_faces
    pts = np.asarray(surface.points)
    if faces is None or len(faces) == 0:
        raise RuntimeError("marching cubes produced no triangles; check SDF")

    # Drop degenerate slivers (see build_spouted_bed_stl_implicit).
    e1 = pts[faces[:, 1]] - pts[faces[:, 0]]
    e2 = pts[faces[:, 2]] - pts[faces[:, 0]]
    tri_n = np.cross(e1, e2)
    norms = np.linalg.norm(tri_n, axis=1)
    keep = norms >= _CROSS_TOL_M2
    faces, tri_n, norms = faces[keep], tri_n[keep], norms[keep]
    if faces.shape[0] == 0:
        raise RuntimeError("all marching-cubes faces were degenerate")

    # Majority vote: flip winding unless most normals already point INTO the
    # fluid (probe at centroid + eps * n lands where sdf < 0).
    centroids = (pts[faces[:, 0]] + pts[faces[:, 1]] + pts[faces[:, 2]]) / 3.0
    probe = centroids + 0.25 * dx_sample_m * (tri_n / norms[:, None])
    inside_frac = float(
        np.mean(sdf_fn(np.hypot(probe[:, 0], probe[:, 2]), probe[:, 1]) < 0.0)
    )
    flip = inside_frac < 0.5

    triangles: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
    for f in faces:
        a, b, c = pts[f[0]], pts[f[1]], pts[f[2]]
        triangles.append((a, c, b) if flip else (a, b, c))
    _write_ascii_stl(triangles, "cone_column_C26_implicit", out_path)
    return out_path


def _cli() -> None:
    import argparse

    p = argparse.ArgumentParser(
        description="Generate the Case V 3D spouted-bed STL (implicit SDF path)."
    )
    p.add_argument("--R-c", type=float, required=True, help="column inside radius [m]")
    p.add_argument("--R-i", type=float, required=True, help="orifice radius [m]")
    p.add_argument(
        "--half-angle", type=float, default=30.0, help="cone half-angle [deg]"
    )
    p.add_argument("--H-dom", type=float, required=True, help="axial domain extent [m]")
    p.add_argument("--L-stub", type=float, required=True, help="inlet stub length [m]")
    p.add_argument(
        "--dx-sample", type=float, required=True, help="SDF sampling step [m]"
    )
    p.add_argument(
        "--top-extend", type=float, default=0.0, help="top-cap extension above H_dom [m]"
    )
    p.add_argument(
        "--R-stub", type=float, default=None, help="stub radius [m] (default: R_i)"
    )
    p.add_argument(
        "--taper-length", type=float, default=0.0, help="stub taper length [m]"
    )
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    out = build_spouted_bed_stl_implicit(
        R_c_m=args.R_c,
        R_i_m=args.R_i,
        cone_half_angle_deg=args.half_angle,
        H_dom_m=args.H_dom,
        L_stub_m=args.L_stub,
        dx_sample_m=args.dx_sample,
        out_path=args.out,
        top_extend_m=args.top_extend,
        R_stub_m=args.R_stub,
        taper_length_m=args.taper_length,
    )
    print(out)


if __name__ == "__main__":
    _cli()

"""Parametric C26 vessel with a one-cell full-face nozzle (2026-10-27).

The C26 authors defined the gas inlet as "mass flow with velocity vector and
point.  Velocity mag is calculated based on the local mesh characteristics"
(personal communication, 2026-10-27; notes/c26_conry2026_case.md "C26 authors
on the point source; full-face nozzle").  Their 2.2.1 sizing argument
(dx = 2.1 mm, the closest value to the 2.38 mm nozzle) then makes ONE cell
face the nozzle: the gas leaves the source cell through its whole top face,
dx^2 = 4.41 mm^2 against the nozzle's 4.45 mm^2, at U = Q / dx^2.

The vessel is the closed-apex cone of src/make_stl_c26_apex.py with a square
bore one cell wide cut into its tip, on a grid whose axis is a cell centre
(odd cell count, 25 across):

    fluid = apex cone  U  bore
    bore  = |x|, |z| < a = dx/2 + delta,  y > y_floor = y_s - delta

The 30 deg wall runs down unchanged until it meets the bore walls, along the
curve where the cone radius equals the bore's radius at that azimuth,
y = (a / max(|cos t|, |sin t|) - R_t) / tan(theta): -0.23 mm at the middle of
each side, +0.54 mm at the corners for C26.  There is no horizontal surface
except the bore floor under the source cell (an earlier version joined a duct
to the cone with a flat ledge at y = 0.754 mm; replaced, notes).  Above the
first apex ring (y = 1.90 mm for C26) the facets are the apex builder's.

* Why the bore sits delta OUTSIDE the cell.  MFiX classifies cut cells by
  their corner nodes (a fluid pocket whose 8 corners are solid is blocked;
  memory "MFiX cut cells need a fluid corner node").  With the walls at
  dx/2 + delta every corner node of the source cell is fluid and none of its
  12 edges is cut, so it is a STANDARD cell with its full top face open.
  Each neighbour in the source row, and the cell below, keeps a delta-thick
  sliver (delta/dx of its volume) plus, for the four face neighbours, the
  cone's protrusion past the bore's side just below y = 0 (R_t > a there);
  below TOL_SMALL_CELL = 1 % (cartesian_grid_init_namelist.f:704) MFiX flags
  them BLOCKED (get_cut_cell_volume_area.f:616) and set_flags.f leaves the
  shared faces walls, so the gas leaves the source cell only upward.  delta
  must exceed TOL_STL = 1e-6 m, below which MFiX drops an edge crossing as a
  corner hit (get_stl_data.f:1799).
* Watertight, inward normals, analytic vertices on the true cone and bore
  planes.  The cone/bore curve carries a vertex at each of the apex
  builder's N azimuths plus the four bore corners; the bore walls are
  vertical strips under it, the floor a fan from its centre, and a collar
  joins the curve to the first apex ring.

Public API:

    fullface_nozzle_geometry(R_t_m, cone_half_angle_deg, cell_half_m,
                             nozzle_y_s_m, wall_offset_m)
    build_cone_column_stl_fullface(R_c_m, R_t_m, cone_half_angle_deg, H_dom_m,
                                   out_path, top_extend_m, sagitta_max_m,
                                   band_height_m, cell_half_m, nozzle_y_s_m,
                                   wall_offset_m)
    inside_fullface_solid(x, y, z, R_c_m, R_t_m, cone_half_angle_deg, y_top_m,
                          half_width_m, y_floor_m)

All lengths SI (m).
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from src.make_stl import _CROSS_TOL_M2, _write_ascii_stl
from src.make_stl_c26_apex import (
    _PROBE_EPS_M,
    _band_edges,
    inside_apex_solid,
    n_theta_for_sagitta,
    sagitta_m,
)

STL_SOLID_NAME = "cone_column_C26_fullface"

# MFiX 23.2 cut-cell limits the nozzle is checked against
# (cartesian_grid/cartesian_grid_init_namelist.f defaults; not set in the deck).
MFIX_TOL_SMALL_CELL = 0.01      # fraction of a cell's volume (line 704)
MFIX_TOL_STL_M = 1.0e-6         # edge-crossing tolerance [m] (line 748)


def _segment_area(R: float, d: float) -> float:
    """Area of the circle r < R beyond the chord at distance d from its centre."""
    if R <= d:
        return 0.0
    return R * R * math.acos(d / R) - d * math.sqrt(R * R - d * d)


def fullface_nozzle_geometry(
    R_t_m: float,
    cone_half_angle_deg: float,
    cell_half_m: float,
    nozzle_y_s_m: float,
    wall_offset_m: float,
) -> dict:
    """Bore of the full-face nozzle and where it meets the apex cone.

    cell_half_m   : half the source cell's width, dx / 2 [m]
    nozzle_y_s_m  : the source cell's bottom face [m] (its top is
                    nozzle_y_s_m + 2 cell_half_m)
    wall_offset_m : delta, bore walls and floor this far outside the cell [m]

    Returns half_width_m (a), y_floor_m, y_nozzle_top_m, y_curve_min_m (mid-
    side) / y_curve_max_m (corners), nozzle_face_area_m2 (the cell's top
    face), bore_area_m2, sliver_fraction (delta / dx), protrusion_fraction
    (cone outside the bore inside one face neighbour of the source cell, / dx^3)
    and neighbour_fluid_fraction (their sum).
    """
    if not (0.0 < cone_half_angle_deg < 90.0):
        raise ValueError("cone_half_angle_deg must be in (0, 90) deg")
    if R_t_m <= 0.0 or cell_half_m <= 0.0:
        raise ValueError("R_t_m and cell_half_m must be > 0")
    if wall_offset_m <= MFIX_TOL_STL_M:
        raise ValueError(f"wall_offset_m ({wall_offset_m}) must exceed MFiX TOL_STL "
                         f"({MFIX_TOL_STL_M} m)")
    tan_h = math.tan(math.radians(cone_half_angle_deg))
    dx = 2.0 * cell_half_m
    a = cell_half_m + wall_offset_m
    y_top = nozzle_y_s_m + dx
    y_floor = nozzle_y_s_m - wall_offset_m
    if -R_t_m / tan_h <= y_floor:
        raise ValueError("the cone's apex lies below the bore floor")
    y_cmin = (a - R_t_m) / tan_h
    y_cmax = (math.sqrt(2.0) * a - R_t_m) / tan_h
    if y_cmax <= y_top:
        raise ValueError("the cone is wider than the bore's corners at the source cell's top "
                         "face: the source row would open sideways")
    # Cone outside the bore, inside the +x face neighbour of the source cell
    # (the row y in [y_s, y_top]): the circle segment beyond x = a, from
    # y_curve_min up to the top face (Simpson's rule on a smooth integrand).
    lo, hi = max(y_cmin, nozzle_y_s_m), y_top
    prot = 0.0
    if hi > lo:
        n = 200
        ys = np.linspace(lo, hi, 2 * n + 1)
        f = np.array([_segment_area(R_t_m + y * tan_h, a) for y in ys])
        prot = (hi - lo) / (6 * n) * (f[0] + f[-1] + 4 * f[1:-1:2].sum() + 2 * f[2:-1:2].sum())
    sliver = wall_offset_m / dx
    neighbour = sliver + prot / dx ** 3
    if neighbour >= MFIX_TOL_SMALL_CELL:
        raise ValueError(f"the source cell's face neighbours keep {neighbour:.3%} fluid; it must "
                         f"stay below TOL_SMALL_CELL ({MFIX_TOL_SMALL_CELL:.0%}) so MFiX blocks them")
    return {
        "half_width_m": a,
        "y_floor_m": y_floor,
        "y_nozzle_top_m": y_top,
        "y_curve_min_m": y_cmin,
        "y_curve_max_m": y_cmax,
        "nozzle_face_area_m2": dx ** 2,
        "bore_area_m2": (2.0 * a) ** 2,
        "sliver_fraction": sliver,
        "protrusion_fraction": prot / dx ** 3,
        "neighbour_fluid_fraction": neighbour,
    }


def inside_fullface_solid(
    x: np.ndarray, y: np.ndarray, z: np.ndarray,
    R_c_m: float, R_t_m: float, cone_half_angle_deg: float, y_top_m: float,
    half_width_m: float, y_floor_m: float,
) -> np.ndarray:
    """True where (x, y, z) lies strictly inside the fluid: the closed-apex
    cone / column, or the square bore below the top cap."""
    x, y, z = np.asarray(x), np.asarray(y), np.asarray(z)
    bore = ((y > y_floor_m) & (y < y_top_m)
            & (np.abs(x) < half_width_m) & (np.abs(z) < half_width_m))
    return bore | inside_apex_solid(x, y, z, R_c_m, R_t_m, cone_half_angle_deg, y_top_m)


def _orient(tri: tuple[int, int, int], verts: np.ndarray, want: np.ndarray) -> tuple[int, int, int]:
    """Return tri, reversed if its normal points away from ``want``."""
    a, b, c = verts[list(tri)]
    return tri if np.dot(np.cross(b - a, c - a), want) > 0.0 else (tri[0], tri[2], tri[1])


def build_cone_column_stl_fullface(
    R_c_m: float,
    R_t_m: float,
    cone_half_angle_deg: float,
    H_dom_m: float,
    out_path: Path,
    top_extend_m: float,
    sagitta_max_m: float,
    band_height_m: float,
    cell_half_m: float,
    nozzle_y_s_m: float,
    wall_offset_m: float,
) -> tuple[Path, dict]:
    """Write the closed-apex cone + column STL with the square bore; return
    (path, facts).

    R_c_m, R_t_m, cone_half_angle_deg, H_dom_m, top_extend_m, sagitta_max_m,
    band_height_m : as src.make_stl_c26_apex.build_cone_column_stl_apex
    cell_half_m, nozzle_y_s_m, wall_offset_m : as :func:`fullface_nozzle_geometry`

    facts: the nozzle geometry, n_theta, sagitta_Rc_m, n_bands_cone,
    n_bands_column, band heights, y_apex_m, h_cone_m, y_top_m, y_collar_top_m
    (the first apex ring), facet_count, vertex_count, ring_y_m.
    """
    out_path = Path(out_path)
    geo = fullface_nozzle_geometry(R_t_m, cone_half_angle_deg, cell_half_m, nozzle_y_s_m,
                                   wall_offset_m)
    if R_t_m >= R_c_m:
        raise ValueError(f"R_t_m ({R_t_m}) must be < R_c_m ({R_c_m})")
    tan_h = math.tan(math.radians(cone_half_angle_deg))
    h_cone_m = (R_c_m - R_t_m) / tan_h
    if h_cone_m >= H_dom_m:
        raise ValueError(f"cone top ({h_cone_m:.4g} m) must lie below H_dom ({H_dom_m} m)")
    if top_extend_m < 0.0:
        raise ValueError("top_extend_m must be >= 0")
    if sagitta_max_m <= 0.0 or band_height_m <= 0.0:
        raise ValueError("sagitta_max_m and band_height_m must be > 0")
    a, y_floor = geo["half_width_m"], geo["y_floor_m"]
    y_apex_m = -R_t_m / tan_h
    y_top_m = H_dom_m + top_extend_m

    # Rings exactly as the apex builder: the apex fan's ring and everything
    # above it are shared; only the fan below that ring is replaced.
    n_theta = n_theta_for_sagitta(R_c_m, sagitta_max_m)
    if n_theta % 8:
        # theta_i = (i + 1/2) 360 / N lands on 45 deg when N / 8 is not whole.
        raise ValueError(f"n_theta = {n_theta}: a ring vertex would sit on a bore corner's "
                         "azimuth; the builder needs a multiple of 8")
    theta = (np.arange(n_theta) + 0.5) * (2.0 * math.pi / n_theta)
    y_cone = _band_edges(y_apex_m, h_cone_m, band_height_m)
    y_col = _band_edges(h_cone_m, y_top_m, band_height_m)
    ring_y = np.concatenate([y_cone[1:], y_col[1:]])
    ring_r = np.minimum(R_t_m + ring_y * tan_h, R_c_m)
    ring_r[len(y_cone) - 1:] = R_c_m
    n_rings = len(ring_y)
    if geo["y_curve_max_m"] >= ring_y[0]:
        raise ValueError("the bore reaches the first cone ring; lower band_height_m")

    verts: list[tuple[float, float, float]] = []

    def add(p) -> int:
        verts.append(tuple(float(v) for v in p))
        return len(verts) - 1

    ring0 = len(verts)
    for k in range(n_rings):
        for t in theta:
            add((ring_r[k] * math.cos(t), ring_y[k], ring_r[k] * math.sin(t)))
    i_cap = add((0.0, y_top_m, 0.0))

    def ring(k: int, i: int) -> int:
        return ring0 + k * n_theta + (i % n_theta)

    # The cone/bore curve: one vertex per apex azimuth (on that generator) and
    # one per bore corner (azimuth 45 + 90 k deg), sorted by azimuth.  The
    # bore's radius at azimuth t is a / max(|cos t|, |sin t|).
    corner_t = math.pi / 4.0 + np.arange(4) * (math.pi / 2.0)
    curve = []          # (azimuth, vertex index, kind, payload)
    for i, t in enumerate(theta):
        rho = a / max(abs(math.cos(t)), abs(math.sin(t)))
        curve.append((t, add((rho * math.cos(t), (rho - R_t_m) / tan_h, rho * math.sin(t))), "gen", i))
    for k, t in enumerate(corner_t):
        rho = a * math.sqrt(2.0)
        curve.append((t, add((rho * math.cos(t), (rho - R_t_m) / tan_h, rho * math.sin(t))), "corner", k))
    curve.sort(key=lambda c: c[0])
    corner_vid = {c[3]: c[1] for c in curve if c[2] == "corner"}
    gen_vid = {c[3]: c[1] for c in curve if c[2] == "gen"}

    tris: list[tuple[int, int, int]] = []

    def radial_in(*ids: int) -> np.ndarray:
        c = np.mean([verts[j] for j in ids], axis=0)
        return np.array([-c[0], 0.0, -c[2]])

    V = lambda: np.asarray(verts)   # noqa: E731  (vertex table grows as we go)

    # Collar: the cone between the curve and the first ring, one strip per
    # pair of neighbouring generators (three triangles where a corner sits
    # between them).
    corner_between = {}
    for k, t in enumerate(corner_t):
        i = int(math.floor(t / (2.0 * math.pi / n_theta) - 0.5)) % n_theta
        corner_between[i] = k          # corner k lies between generators i and i+1
    for i in range(n_theta):
        j = (i + 1) % n_theta
        ci, cj, ri, rj = gen_vid[i], gen_vid[j], ring(0, i), ring(0, j)
        if i in corner_between:
            ck = corner_vid[corner_between[i]]
            new = [(ci, ck, ri), (ck, rj, ri), (ck, cj, rj)]
        else:
            new = [(ci, cj, rj), (ci, rj, ri)]
        verts_np = V()
        tris += [_orient(t, verts_np, radial_in(*t)) for t in new]

    # Bore walls: side k runs from corner k to corner k+1 (azimuths 45 + 90 k
    # to 135 + 90 k deg) on the plane z = +a, x = -a, z = -a, x = +a; one
    # vertical strip under each curve segment, down to the floor.
    fc = [add((verts[corner_vid[k]][0], y_floor, verts[corner_vid[k]][2])) for k in range(4)]
    boundary: list[int] = []       # floor boundary, in azimuth order
    for k in range(4):
        t0 = corner_t[k]
        gens = sorted((i for i in range(n_theta)
                       if 0.0 < (theta[i] - t0) % (2.0 * math.pi) < math.pi / 2.0),
                      key=lambda i: (theta[i] - t0) % (2.0 * math.pi))
        top = [corner_vid[k]] + [gen_vid[i] for i in gens] + [corner_vid[(k + 1) % 4]]
        bottom = ([fc[k]]
                  + [add((verts[gen_vid[i]][0], y_floor, verts[gen_vid[i]][2])) for i in gens]
                  + [fc[(k + 1) % 4]])
        verts_np = V()
        want = radial_in(*top)
        for n in range(len(top) - 1):
            p, q, s, u = bottom[n], bottom[n + 1], top[n], top[n + 1]
            tris += [_orient((p, q, u), verts_np, want), _orient((p, u, s), verts_np, want)]
        boundary += bottom[:-1]
    # Floor: fan from its centre (fluid above).
    centre = add((0.0, y_floor, 0.0))
    verts_np = V()
    up = np.array([0.0, 1.0, 0.0])
    nb = len(boundary)
    tris += [_orient((centre, boundary[n], boundary[(n + 1) % nb]), verts_np, up)
             for n in range(nb)]
    # Cone and column bands, top-cap fan (winding as the apex builder).
    for k in range(n_rings - 1):
        for i in range(n_theta):
            p, q = ring(k, i), ring(k, i + 1)
            r, s = ring(k + 1, i + 1), ring(k + 1, i)
            tris += [(p, q, r), (p, r, s)]
    kt = n_rings - 1
    tris += [(i_cap, ring(kt, i), ring(kt, i + 1)) for i in range(n_theta)]

    verts_np = V()
    faces = np.asarray(tris, dtype=np.int64)

    # --- self checks (raise rather than write a bad shell) ----------------
    edges = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
    edges.sort(axis=1)
    _, counts = np.unique(edges, axis=0, return_counts=True)
    if not np.all(counts == 2):
        raise RuntimeError("full-face STL is not watertight (edge not shared by exactly 2 facets)")
    A, B, C = verts_np[faces[:, 0]], verts_np[faces[:, 1]], verts_np[faces[:, 2]]
    cross = np.cross(B - A, C - A)
    norm = np.linalg.norm(cross, axis=1)
    if np.any(norm < _CROSS_TOL_M2):
        raise RuntimeError("full-face STL has a degenerate facet")
    n_hat = cross / norm[:, None]
    cen = (A + B + C) / 3.0
    args = (R_c_m, R_t_m, cone_half_angle_deg, y_top_m, a, y_floor)
    p_in = cen + _PROBE_EPS_M * n_hat
    p_out = cen - _PROBE_EPS_M * n_hat
    ok_in = inside_fullface_solid(p_in[:, 0], p_in[:, 1], p_in[:, 2], *args)
    ok_out = ~inside_fullface_solid(p_out[:, 0], p_out[:, 1], p_out[:, 2], *args)
    if not (ok_in.all() and ok_out.all()):
        raise RuntimeError(
            f"full-face STL normals: {int((~ok_in).sum())} not inward, "
            f"{int((~ok_out).sum())} fail the outward probe")
    flat = np.abs(n_hat[:, 1]) > 1.0 - 1e-9
    if not np.all(np.isclose(cen[flat, 1], y_floor) | np.isclose(cen[flat, 1], y_top_m)):
        raise RuntimeError("a horizontal facet other than the bore floor and the top cap")

    triangles = [(verts_np[f[0]], verts_np[f[1]], verts_np[f[2]]) for f in faces]
    _write_ascii_stl(triangles, STL_SOLID_NAME, out_path)

    facts = dict(geo)
    facts.update({
        "n_theta": int(n_theta),
        "sagitta_max_m": float(sagitta_max_m),
        "sagitta_Rc_m": sagitta_m(R_c_m, n_theta),
        "band_height_max_m": float(band_height_m),
        "n_bands_cone": int(len(y_cone) - 1),
        "n_bands_column": int(len(y_col) - 1),
        "band_height_cone_m": float(y_cone[1] - y_cone[0]),
        "band_height_column_m": float(y_col[1] - y_col[0]),
        "y_apex_m": float(y_apex_m),
        "h_cone_m": float(h_cone_m),
        "y_top_m": float(y_top_m),
        "y_collar_top_m": float(ring_y[0]),
        "facet_count": int(len(faces)),
        "vertex_count": int(len(verts_np)),
        "ring_y_m": [float(v) for v in ring_y],
    })
    return out_path, facts

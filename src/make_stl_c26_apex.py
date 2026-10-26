"""Parametric closed-apex cone + column STL for Case C26 (2026-10-24).

A surface of revolution about the y-axis with ANALYTIC vertices (no signed
distance field, no marching cubes), written through
:func:`src.make_stl._write_ascii_stl`.  Profile, bottom to top::

    apex (0, y_apex)  --cone--  (R_c, h_cone)  --column--  (R_c, y_top)

    y_apex = -R_t / tan(theta)          virtual apex of the cone
    r(y)   = R_t + y tan(theta)         same cone equation as the implicit
                                        builder for y >= 0
    h_cone = (R_c - R_t) / tan(theta)   cone / column junction
    y_top  = H_dom + top_extend         top cap, above the PO plane

Design (see notes/c26_conry2026_case.md "Closed apex, parametric STL"):

* Why the apex.  C26 reports neither the throat nor the nozzle geometry; the
  implicit builder closes the cone with a flat cap at y = 0 and puts the point
  source in the first row above it.  Continuing the 30 deg wall to its own
  apex instead makes the open cross-section at y = 0 exactly pi R_t^2, the
  nozzle area, so a mass-only point source placed in the tip gives the nozzle
  velocity at y = 0 by continuity.  For C26 (R_t = 1.1906 mm) the apex is at
  y_apex = -2.0622 mm, 0.038 mm above the grid face at y = -2.1 mm, so the tip
  occupies the single grid row y in [-2.1, 0] mm.  The tip volume is
  (1/3) pi R_t^2 |y_apex| = 3.06 mm^3, shared by the four cells around the
  axis (8 % of a 9.26 mm^3 cell each, above MFiX's TOL_SMALL_CELL = 1 %,
  cartesian_grid_init_namelist.f:704).  The apex vertex is placed exactly at
  (0, y_apex, 0).
* Watertight.  Apex fan + cone bands + column bands + top-cap fan.  Every
  undirected edge is shared by exactly two triangles (asserted before
  writing).  No opening: the gas enters through the point source, not a
  boundary.  Fans at the apex and the cap centre avoid zero-area triangles;
  no quad is ever collapsed.  Day 10 (notes/day10_summary.md "Why PIC +
  NETL-clone geometry") showed MFiX 26.1.2's F_AT classifier marking the
  near-axis column of an OPEN surface-of-revolution STL as blocked; a closed
  shell is the fix.
* Normals INWARD (MFiX OUT_STL_VALUE = 1 internal-flow convention; header of
  src/make_stl.py).  Checked per facet by a centroid probe against the
  analytic solid (centroid + eps n inside, centroid - eps n outside), not by
  a majority vote.
* Angular resolution from a sagitta criterion.  An N-gon inscribed in a
  circle of radius R has chord error (sagitta) s = R (1 - cos(pi/N)).  The
  caller passes the largest acceptable sagitta; for C26 that is 0.1 x the
  particle radius (20.6 um; geometry.yaml stl_apex_sagitta_max_m), i.e.
  facet ridges a tenth of a particle radius, 1/20 of a parcel radius, with
  the contact normal alternating by +/- pi/N about the true one.  N is the
  smallest multiple of 4 satisfying the criterion at R_c (the 4-fold
  rounding keeps the shell symmetric under the grid's x -> -x and z -> -z
  reflections, so the four tip cells get identical volumes).  C26:
  N = 80, s(R_c) = 19.4 um, s(R_t) = 0.92 um, polygon area deficit 0.10 %.
* Vertices sit at theta_i = (i + 1/2) 2 pi / N, i.e. no vertex on the
  planes x = 0 or z = 0 (grid planes on the axis-on-face layout) and the
  wall's largest |x| is R_c cos(pi/N) < R_c = x_max, so no facet edge lies
  on a cell face or a box face.
* Axial bands.  MFiX bins facets to cells by facet BOUNDING BOX: the cut-cell
  intersection loop (cartesian_grid/intersect.f ~800-900) marks every cell
  in the box as a potential cut cell and tests its 12 edges; the DEM wall
  search (des/stl_preproc_des_mod.f BIN_FACETS_TO_DG) sweeps the box +/- 1
  DES cell before the exact triangle/box test.  One facet spanning the
  43.6 mm cone would sweep ~20 rows x most of the cross-section.  The cone
  and column are ruled surfaces, exact along their generators, so axial
  subdivision adds no accuracy; it only localises facets.  Bands of
  <= band_height_m (C26: 2 fluid cells, 4.2 mm) keep every facet's box
  within 3 cell rows; 1-cell bands would double the facet count for no
  gain in accuracy, taller bands grow the bounding-box sweep linearly.
  C26: 11 cone bands (the first is the apex fan) + 42 column bands,
  8,480 facets (implicit: 77,088).

Public API:

    build_cone_column_stl_apex(R_c_m, R_t_m, cone_half_angle_deg, H_dom_m,
                               out_path, top_extend_m, sagitta_max_m,
                               band_height_m) -> (Path, dict)
    n_theta_for_sagitta(R_m, sagitta_max_m) -> int
    inside_apex_solid(x, y, z, R_c_m, R_t_m, cone_half_angle_deg, y_top_m)

All lengths SI (m).
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from src.make_stl import _CROSS_TOL_M2, _write_ascii_stl

STL_SOLID_NAME = "cone_column_C26_apex"

# Probe distance for the per-facet normal check [m].  Must exceed the sagitta
# at R_c (the centroid of a chord facet lies inside the true circle by the
# sagitta), and stay well below every feature size (R_t, band height).
_PROBE_EPS_M = 1.0e-4


def n_theta_for_sagitta(R_m: float, sagitta_max_m: float, multiple_of: int = 4) -> int:
    """Smallest N (a multiple of ``multiple_of``) with R (1 - cos(pi/N)) <= sagitta."""
    if R_m <= 0.0 or sagitta_max_m <= 0.0:
        raise ValueError("R_m and sagitta_max_m must be > 0")
    if sagitta_max_m >= R_m:
        n_min = 3
    else:
        n_min = math.ceil(math.pi / math.acos(1.0 - sagitta_max_m / R_m))
    n_min = max(n_min, 3)
    return int(math.ceil(n_min / multiple_of) * multiple_of)


def sagitta_m(R_m: float, n_theta: int) -> float:
    """Chord error of a regular n_theta-gon inscribed in a circle of radius R."""
    return R_m * (1.0 - math.cos(math.pi / n_theta))


def inside_apex_solid(
    x: np.ndarray, y: np.ndarray, z: np.ndarray,
    R_c_m: float, R_t_m: float, cone_half_angle_deg: float, y_top_m: float,
) -> np.ndarray:
    """True where (x, y, z) lies strictly inside the fluid region of the
    apex-closed vessel: y_apex < y < y_top and r < min(R_t + y tan(theta), R_c)."""
    tan_h = math.tan(math.radians(cone_half_angle_deg))
    y_apex = -R_t_m / tan_h
    r = np.hypot(np.asarray(x), np.asarray(z))
    y = np.asarray(y)
    r_wall = np.minimum(R_t_m + y * tan_h, R_c_m)
    return (y > y_apex) & (y < y_top_m) & (r < r_wall)


def _band_edges(y0: float, y1: float, band_height_m: float) -> np.ndarray:
    """Uniform band boundaries from y0 to y1 with spacing <= band_height_m
    (end points included)."""
    n = max(1, int(math.ceil((y1 - y0) / band_height_m - 1.0e-9)))
    return np.linspace(y0, y1, n + 1)


def build_cone_column_stl_apex(
    R_c_m: float,
    R_t_m: float,
    cone_half_angle_deg: float,
    H_dom_m: float,
    out_path: Path,
    top_extend_m: float,
    sagitta_max_m: float,
    band_height_m: float,
) -> tuple[Path, dict]:
    """Write the closed-apex cone + column STL; return (path, facts).

    R_c_m               : column inside radius [m]
    R_t_m               : cone radius at y = 0 (the open area pi R_t^2) [m]
    cone_half_angle_deg : cone half-angle from the vertical axis [deg]
    H_dom_m             : pressure-outlet plane [m]
    out_path            : destination ASCII .stl
    top_extend_m        : top cap this far above H_dom (>= 0) [m]
    sagitta_max_m       : largest chord error of the faceted wall at R_c [m]
    band_height_m       : largest axial extent of one facet band [m]

    facts: n_theta, sagitta_Rc_m, sagitta_Rt_m, n_bands_cone, n_bands_column,
    band_height_cone_m, band_height_column_m, y_apex_m, h_cone_m, y_top_m,
    open_area_y0_m2 (the N-gon at y = 0), open_area_y0_exact_m2 (pi R_t^2),
    tip_volume_m3 (cone tip, apex to y = 0), tip_volume_polygon_m3,
    facet_count, vertex_count, ring_y_m (every ring height).
    """
    out_path = Path(out_path)
    if not (0.0 < cone_half_angle_deg < 90.0):
        raise ValueError("cone_half_angle_deg must be in (0, 90) deg")
    if R_t_m <= 0.0:
        raise ValueError(f"R_t_m ({R_t_m}) must be > 0")
    if R_t_m >= R_c_m:
        raise ValueError(f"R_t_m ({R_t_m}) must be < R_c_m ({R_c_m})")
    tan_h = math.tan(math.radians(cone_half_angle_deg))
    h_cone_m = (R_c_m - R_t_m) / tan_h
    if h_cone_m >= H_dom_m:
        raise ValueError(f"cone top ({h_cone_m:.4g} m) must lie below H_dom ({H_dom_m} m)")
    if top_extend_m < 0.0:
        raise ValueError("top_extend_m must be >= 0")
    if sagitta_max_m <= 0.0:
        raise ValueError("sagitta_max_m must be > 0")
    if band_height_m <= 0.0:
        raise ValueError("band_height_m must be > 0")

    y_apex_m = -R_t_m / tan_h
    y_top_m = H_dom_m + top_extend_m
    n_theta = n_theta_for_sagitta(R_c_m, sagitta_max_m)
    theta = (np.arange(n_theta) + 0.5) * (2.0 * math.pi / n_theta)
    cos_t, sin_t = np.cos(theta), np.sin(theta)

    # Ring heights: apex (radius 0, a single vertex), cone bands up to h_cone,
    # column bands up to y_top.  The ring at h_cone is shared.
    y_cone = _band_edges(y_apex_m, h_cone_m, band_height_m)      # includes apex
    y_col = _band_edges(h_cone_m, y_top_m, band_height_m)         # includes h_cone
    ring_y = np.concatenate([y_cone[1:], y_col[1:]])               # rings with vertices
    ring_r = np.minimum(R_t_m + ring_y * tan_h, R_c_m)
    ring_r[len(y_cone) - 1:] = R_c_m                               # column rings exactly R_c
    n_rings = len(ring_y)

    # Vertex table: apex, rings (n_rings x n_theta), cap centre.
    verts = np.empty((2 + n_rings * n_theta, 3))
    verts[0] = (0.0, y_apex_m, 0.0)
    for k in range(n_rings):
        base = 1 + k * n_theta
        verts[base:base + n_theta, 0] = ring_r[k] * cos_t
        verts[base:base + n_theta, 1] = ring_y[k]
        verts[base:base + n_theta, 2] = ring_r[k] * sin_t
    i_cap = 1 + n_rings * n_theta
    verts[i_cap] = (0.0, y_top_m, 0.0)

    def ring(k: int, i: int) -> int:
        return 1 + k * n_theta + (i % n_theta)

    # Winding: (x, y, z) right-handed with y up; theta runs x -> z, which is
    # clockwise seen from +y.  With the lower ring first and theta increasing
    # along the first edge, (v_{k,i}, v_{k,i+1}, v_{k+1,i+1}) has an inward
    # normal; the fans follow the same rule with the lower (apex) or upper
    # (cap centre) ring collapsed to a point.  Verified below, every facet.
    tris: list[tuple[int, int, int]] = []
    for i in range(n_theta):                                       # apex fan
        tris.append((0, ring(0, i + 1), ring(0, i)))
    for k in range(n_rings - 1):                                   # wall bands
        for i in range(n_theta):
            a, b = ring(k, i), ring(k, i + 1)
            c, d = ring(k + 1, i + 1), ring(k + 1, i)
            tris.append((a, b, c))
            tris.append((a, c, d))
    kt = n_rings - 1
    for i in range(n_theta):                                       # top-cap fan
        tris.append((i_cap, ring(kt, i), ring(kt, i + 1)))
    faces = np.asarray(tris, dtype=np.int64)

    # --- self checks (raise rather than write a bad shell) ----------------
    # 1. closed 2-manifold: each undirected edge in exactly two triangles.
    edges = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
    edges.sort(axis=1)
    _, counts = np.unique(edges, axis=0, return_counts=True)
    if not np.all(counts == 2):
        raise RuntimeError("apex STL is not watertight (edge not shared by exactly 2 facets)")
    # 2. no degenerate facet.
    A, B, C = verts[faces[:, 0]], verts[faces[:, 1]], verts[faces[:, 2]]
    cross = np.cross(B - A, C - A)
    norm = np.linalg.norm(cross, axis=1)
    if np.any(norm < _CROSS_TOL_M2):
        raise RuntimeError("apex STL has a degenerate facet")
    # 3. every normal points into the fluid (centroid probe both ways).
    n_hat = cross / norm[:, None]
    cen = (A + B + C) / 3.0
    p_in = cen + _PROBE_EPS_M * n_hat
    p_out = cen - _PROBE_EPS_M * n_hat
    ok_in = inside_apex_solid(p_in[:, 0], p_in[:, 1], p_in[:, 2],
                              R_c_m, R_t_m, cone_half_angle_deg, y_top_m)
    ok_out = ~inside_apex_solid(p_out[:, 0], p_out[:, 1], p_out[:, 2],
                                R_c_m, R_t_m, cone_half_angle_deg, y_top_m)
    if not (ok_in.all() and ok_out.all()):
        raise RuntimeError(
            f"apex STL normals: {int((~ok_in).sum())} not inward, "
            f"{int((~ok_out).sum())} fail the outward probe")

    triangles = [(verts[f[0]], verts[f[1]], verts[f[2]]) for f in faces]
    _write_ascii_stl(triangles, STL_SOLID_NAME, out_path)

    area_poly_factor = (n_theta / (2.0 * math.pi)) * math.sin(2.0 * math.pi / n_theta)
    open_exact = math.pi * R_t_m ** 2
    tip_exact = open_exact * (-y_apex_m) / 3.0
    facts = {
        "n_theta": int(n_theta),
        "sagitta_max_m": float(sagitta_max_m),
        "sagitta_Rc_m": sagitta_m(R_c_m, n_theta),
        "sagitta_Rt_m": sagitta_m(R_t_m, n_theta),
        "band_height_max_m": float(band_height_m),
        "n_bands_cone": int(len(y_cone) - 1),
        "n_bands_column": int(len(y_col) - 1),
        "band_height_cone_m": float(y_cone[1] - y_cone[0]),
        "band_height_column_m": float(y_col[1] - y_col[0]),
        "y_apex_m": float(y_apex_m),
        "h_cone_m": float(h_cone_m),
        "y_top_m": float(y_top_m),
        "open_area_y0_exact_m2": open_exact,
        "open_area_y0_m2": open_exact * area_poly_factor,
        "tip_volume_m3": tip_exact,
        "tip_volume_polygon_m3": tip_exact * area_poly_factor,
        "facet_count": int(len(faces)),
        "vertex_count": int(len(verts)),
        "ring_y_m": [float(v) for v in ring_y],
    }
    return out_path, facts


def _cli() -> None:
    import argparse
    import json

    p = argparse.ArgumentParser(description="Closed-apex C26 cone + column STL (parametric).")
    p.add_argument("--R-c", type=float, required=True, help="column inside radius [m]")
    p.add_argument("--R-t", type=float, required=True, help="cone radius at y = 0 [m]")
    p.add_argument("--half-angle", type=float, default=30.0, help="cone half-angle [deg]")
    p.add_argument("--H-dom", type=float, required=True, help="outlet plane [m]")
    p.add_argument("--top-extend", type=float, default=0.0, help="cap above H_dom [m]")
    p.add_argument("--sagitta", type=float, required=True, help="max chord error at R_c [m]")
    p.add_argument("--band-height", type=float, required=True, help="max facet band height [m]")
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    out, facts = build_cone_column_stl_apex(
        R_c_m=a.R_c, R_t_m=a.R_t, cone_half_angle_deg=a.half_angle, H_dom_m=a.H_dom,
        out_path=a.out, top_extend_m=a.top_extend, sagitta_max_m=a.sagitta,
        band_height_m=a.band_height,
    )
    facts.pop("ring_y_m")
    print(out)
    print(json.dumps(facts, indent=2))


if __name__ == "__main__":
    _cli()

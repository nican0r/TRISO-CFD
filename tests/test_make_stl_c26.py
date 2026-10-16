"""Tests for src.make_stl.build_cone_column_stl_implicit (Case C26 vessel)."""

from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np
import pytest

from src.make_stl import build_cone_column_stl_implicit

# C26 geometry (params/geometry.yaml bed_C26), shortened column for speed.
R_C = 0.0252
R_T = 0.5 * 2.38125e-3
HALF = 30.0
H_DOM = 0.060
DX_S = 0.25 * 2.1e-3
TOP = 5.0e-3
H_CONE = (R_C - R_T) / math.tan(math.radians(HALF))


def _read(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Return (normals[n,3], vertices[n,3,3]) from an ASCII STL."""
    txt = path.read_text()
    normals = np.array(
        [[float(v) for v in m.groups()] for m in
         re.finditer(r"facet normal (\S+) (\S+) (\S+)", txt)])
    verts = np.array(
        [[float(v) for v in m.groups()] for m in
         re.finditer(r"vertex (\S+) (\S+) (\S+)", txt)]).reshape(-1, 3, 3)
    return normals, verts


@pytest.fixture(scope="module")
def stl(tmp_path_factory):
    out = tmp_path_factory.mktemp("c26") / "geometry_0001.stl"
    build_cone_column_stl_implicit(
        R_c_m=R_C, R_t_m=R_T, cone_half_angle_deg=HALF, H_dom_m=H_DOM,
        dx_sample_m=DX_S, out_path=out, top_extend_m=TOP,
    )
    return _read(out)


def test_header_and_footer(tmp_path):
    out = tmp_path / "g.stl"
    build_cone_column_stl_implicit(R_C, R_T, HALF, H_DOM, DX_S, out, TOP)
    lines = out.read_text().strip().splitlines()
    assert lines[0].startswith("solid ") and lines[-1].startswith("endsolid ")


def test_bounding_box(stl):
    _, v = stl
    pts = v.reshape(-1, 3)
    r = np.hypot(pts[:, 0], pts[:, 2])
    assert pts[:, 1].min() == pytest.approx(0.0, abs=DX_S)
    assert pts[:, 1].max() == pytest.approx(H_DOM + TOP, abs=DX_S)
    assert r.max() == pytest.approx(R_C, abs=DX_S)


def test_no_degenerate_or_nan_facets(stl):
    n, v = stl
    assert np.isfinite(v).all() and np.isfinite(n).all()
    cross = np.cross(v[:, 1] - v[:, 0], v[:, 2] - v[:, 0])
    assert (np.linalg.norm(cross, axis=1) >= 1e-12).all()


def test_closed_cap_at_throat(stl):
    # Bottom cap: facets at y ~ 0 facing +y (into the fluid), r <= R_t.
    n, v = stl
    c = v.mean(axis=1)
    r = np.hypot(c[:, 0], c[:, 2])
    cap = (np.abs(c[:, 1]) < 0.5 * DX_S) & (np.abs(n[:, 1]) > 0.7)
    assert cap.sum() > 0
    assert r[cap].max() <= R_T + DX_S
    assert (n[cap, 1] > 0).all()           # inward = +y at the bottom


def test_cone_slope(stl):
    _, v = stl
    c = v.mean(axis=1)
    band = (c[:, 1] > 3e-3) & (c[:, 1] < H_CONE - 3e-3)
    r = np.hypot(c[band, 0], c[band, 2])
    r_exp = R_T + c[band, 1] * math.tan(math.radians(HALF))
    assert np.abs(r - r_exp).max() < DX_S


def test_column_radius(stl):
    _, v = stl
    c = v.mean(axis=1)
    band = (c[:, 1] > H_CONE + 3e-3) & (c[:, 1] < H_DOM)
    r = np.hypot(c[band, 0], c[band, 2])
    assert np.abs(r - R_C).max() < DX_S


def test_sidewall_normals_point_inward(stl):
    n, v = stl
    c = v.mean(axis=1)
    side = (c[:, 1] > 2e-3) & (c[:, 1] < H_DOM) & (np.abs(n[:, 1]) < 0.9)
    r = np.hypot(c[side, 0], c[side, 2])
    radial = (n[side, 0] * c[side, 0] + n[side, 2] * c[side, 2]) / r
    assert float(np.mean(radial < 0)) > 0.95


def test_no_membrane_at_cone_column_seam(stl):
    # R5 lesson: no horizontal facets inside the vessel at y = h_cone.
    n, v = stl
    c = v.mean(axis=1)
    r = np.hypot(c[:, 0], c[:, 2])
    seam = ((np.abs(c[:, 1] - H_CONE) < 2 * DX_S) & (np.abs(n[:, 1]) > 0.7)
            & (r < R_C - 2 * DX_S))
    assert seam.sum() == 0


def test_no_internal_horizontal_facets(stl):
    # Away from the caps, any near-horizontal facet must sit on the vessel
    # wall (marching-cubes staircase), never out in the fluid.
    n, v = stl
    c = v.mean(axis=1)
    r = np.hypot(c[:, 0], c[:, 2])
    r_wall = np.minimum(R_T + c[:, 1] * math.tan(math.radians(HALF)), R_C)
    horiz = np.abs(n[:, 1]) > 0.9
    interior = (c[:, 1] > DX_S) & (c[:, 1] < H_DOM + TOP - DX_S)
    off_wall = r < r_wall - 2 * DX_S
    assert (horiz & interior & off_wall).sum() == 0


@pytest.mark.parametrize("kw", [
    {"R_t_m": 0.0}, {"R_t_m": R_C}, {"cone_half_angle_deg": 0.0},
    {"dx_sample_m": 2 * R_T}, {"H_dom_m": 0.01}, {"top_extend_m": -1e-3},
])
def test_invalid_inputs_rejected(tmp_path, kw):
    args = dict(R_c_m=R_C, R_t_m=R_T, cone_half_angle_deg=HALF, H_dom_m=H_DOM,
                dx_sample_m=DX_S, out_path=tmp_path / "x.stl", top_extend_m=TOP)
    args.update(kw)
    with pytest.raises(ValueError):
        build_cone_column_stl_implicit(**args)

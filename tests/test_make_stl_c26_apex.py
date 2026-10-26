"""Tests for src.make_stl_c26_apex (closed-apex parametric C26 vessel) and its
wiring into the C26 case generator (stl_builder = parametric_apex)."""

from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np
import pytest

from src.make_case import REPO_ROOT, _load_yaml, render_case
from src.make_case_c26 import build_C26_3d_params, load_sweep
from src.make_stl_c26_apex import (
    build_cone_column_stl_apex,
    inside_apex_solid,
    n_theta_for_sagitta,
    sagitta_m,
)

C26_TEMPLATE = REPO_ROOT / "cases" / "C26_3d" / "template_dem.mfx.j2"

# C26 geometry (params/geometry.yaml bed_C26 / spouted_bed_C26_3d), full size.
R_C = 0.0252
R_T = 0.5 * 2.38125e-3
HALF = 30.0
H_DOM = 0.21
TOP = 5.0e-3
DX = 2.1e-3
SAGITTA = 2.06e-5           # 0.1 x particle radius (d_p 411.6 um)
BAND = 2 * DX
TAN = math.tan(math.radians(HALF))
H_CONE = (R_C - R_T) / TAN
Y_APEX = -R_T / TAN
Y_TOP = H_DOM + TOP


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
def built(tmp_path_factory):
    out = tmp_path_factory.mktemp("c26apex") / "geometry_0001.stl"
    path, facts = build_cone_column_stl_apex(
        R_c_m=R_C, R_t_m=R_T, cone_half_angle_deg=HALF, H_dom_m=H_DOM,
        out_path=out, top_extend_m=TOP, sagitta_max_m=SAGITTA, band_height_m=BAND,
    )
    n, v = _read(path)
    return path, facts, n, v


# ---------------- resolution choice ----------------

def test_sagitta_hand_values():
    # s = R (1 - cos(pi/N)); 80-gon at 25.2 mm -> 19.4 um; at R_t -> 0.92 um
    assert sagitta_m(R_C, 80) == pytest.approx(1.9428e-5, rel=1e-3)
    assert sagitta_m(R_T, 80) == pytest.approx(9.18e-7, rel=2e-3)


def test_n_theta_from_sagitta_criterion():
    # 0.1 r_p at R_c needs N >= 78 (pi / acos(1 - 2.06e-5 / 0.0252) = 77.7);
    # rounded up to a multiple of 4 -> 80, and 76 would violate it.
    assert n_theta_for_sagitta(R_C, SAGITTA) == 80
    assert sagitta_m(R_C, 80) <= SAGITTA < sagitta_m(R_C, 76)
    assert n_theta_for_sagitta(R_C, 1.0e-3) % 4 == 0
    with pytest.raises(ValueError):
        n_theta_for_sagitta(R_C, 0.0)


# ---------------- the shell ----------------

def test_header_footer_and_facet_count(built):
    path, facts, n, v = built
    lines = path.read_text().strip().splitlines()
    assert lines[0] == "solid cone_column_C26_apex" and lines[-1].startswith("endsolid ")
    # 80 (apex fan) + 52 wall bands x 160 + 80 (cap fan) = 8480, far below the
    # implicit builder's 77,088.
    n_rings = facts["n_bands_cone"] + facts["n_bands_column"]
    expected = 2 * facts["n_theta"] * n_rings
    assert facts["facet_count"] == len(n) == expected == 8480
    assert facts["facet_count"] < 20000
    assert (facts["n_theta"], facts["n_bands_cone"], facts["n_bands_column"]) == (80, 11, 42)


def test_watertight_every_edge_shared_by_two_facets(built):
    _, _, _, v = built
    pts = v.reshape(-1, 3)
    uniq, inv = np.unique(pts, axis=0, return_inverse=True)
    tri = inv.reshape(-1, 3)
    edges = np.concatenate([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]])
    edges.sort(axis=1)
    _, counts = np.unique(edges, axis=0, return_counts=True)
    assert counts.min() == 2 and counts.max() == 2
    # Euler characteristic of a closed genus-0 surface: V - E + F = 2.
    assert len(uniq) - len(counts) + len(tri) == 2


def test_no_degenerate_or_nan_facets(built):
    _, _, n, v = built
    assert np.isfinite(v).all() and np.isfinite(n).all()
    cross = np.cross(v[:, 1] - v[:, 0], v[:, 2] - v[:, 0])
    area2 = np.linalg.norm(cross, axis=1)
    assert area2.min() >= 1e-12
    # written normals are unit and agree with the vertex order
    assert np.allclose(np.linalg.norm(n, axis=1), 1.0, atol=1e-5)
    assert np.allclose(cross / area2[:, None], n, atol=1e-5)


def test_all_normals_point_into_the_fluid(built):
    # Centroid probe: a step along the normal lands inside the analytic fluid
    # region, a step against it lands outside -- for EVERY facet.
    _, _, n, v = built
    c = v.mean(axis=1)
    eps = 1.0e-4
    p_in, p_out = c + eps * n, c - eps * n
    assert inside_apex_solid(p_in[:, 0], p_in[:, 1], p_in[:, 2], R_C, R_T, HALF, Y_TOP).all()
    assert not inside_apex_solid(p_out[:, 0], p_out[:, 1], p_out[:, 2], R_C, R_T, HALF, Y_TOP).any()
    # and every wall normal (cap facets excluded: n = -y exactly) points
    # toward the axis
    r = np.hypot(c[:, 0], c[:, 2])
    side = (r > 1e-6) & (np.abs(n[:, 1]) < 0.99)
    radial = (n[side, 0] * c[side, 0] + n[side, 2] * c[side, 2]) / r[side]
    assert (radial < 0).all() and side.sum() == len(n) - 80


def test_apex_vertex_location(built):
    _, facts, _, v = built
    pts = v.reshape(-1, 3)
    i = np.argmin(pts[:, 1])
    assert facts["y_apex_m"] == pytest.approx(Y_APEX, rel=1e-12)
    assert Y_APEX == pytest.approx(-2.0622e-3, rel=1e-4)
    assert pts[i] == pytest.approx([0.0, Y_APEX, 0.0], abs=1e-9)
    # the apex sits 0.038 mm above the cell face at y = -dx, inside the tip row
    assert -DX < Y_APEX < 0.0
    assert Y_APEX - (-DX) == pytest.approx(3.78e-5, rel=2e-2)
    # the apex fan: exactly n_theta facets touch the apex
    touches = np.isclose(v[:, :, 1], Y_APEX, atol=1e-9).any(axis=1)
    assert touches.sum() == facts["n_theta"]


def _section_polygon_area(v: np.ndarray, y0: float) -> float:
    """Area of the shell's cross-section at the plane y = y0 (points where
    facet edges cross the plane, sorted by azimuth, shoelace)."""
    pts = []
    for tri in v:
        for a, b in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
            if (a[1] - y0) * (b[1] - y0) < 0.0:
                t = (y0 - a[1]) / (b[1] - a[1])
                pts.append(a + t * (b - a))
    pts = np.unique(np.round(np.asarray(pts), 12), axis=0)
    ang = np.arctan2(pts[:, 2], pts[:, 0])
    order = np.argsort(ang)
    x, z = pts[order, 0], pts[order, 2]
    return 0.5 * abs(np.dot(x, np.roll(z, -1)) - np.dot(z, np.roll(x, -1)))


def test_open_area_at_y0_is_the_nozzle_area(built):
    _, facts, _, v = built
    A_exact = math.pi * R_T ** 2
    A_poly = _section_polygon_area(v, 0.0)
    # inscribed 80-gon: deficit 1 - (N/2pi) sin(2pi/N) = 0.10 %, within the
    # sagitta error 2 s(R_t) / R_t = 0.15 %
    assert abs(A_poly - A_exact) / A_exact <= 2.0 * sagitta_m(R_T, 80) / R_T
    assert A_poly == pytest.approx(facts["open_area_y0_m2"], rel=1e-6)
    assert facts["open_area_y0_exact_m2"] == pytest.approx(A_exact)
    assert A_exact == pytest.approx(4.4535e-6, rel=1e-4)
    # tip volume (apex to y = 0): (1/3) pi R_t^2 |y_apex| = 3.06 mm^3
    assert facts["tip_volume_m3"] == pytest.approx(A_exact * (-Y_APEX) / 3.0)
    assert facts["tip_volume_m3"] == pytest.approx(3.061e-9, rel=1e-3)


def test_cone_and_column_radii(built):
    _, _, _, v = built
    pts = v.reshape(-1, 3)
    r = np.hypot(pts[:, 0], pts[:, 2])
    cone = (pts[:, 1] > Y_APEX + 1e-9) & (pts[:, 1] < H_CONE - 1e-9)
    # 1e-8 m: the ASCII STL carries 7 significant digits (%.6e)
    assert np.allclose(r[cone], R_T + pts[cone, 1] * TAN, atol=1e-8)
    col = pts[:, 1] > H_CONE + 1e-9
    assert np.allclose(r[col][r[col] > 1e-9], R_C, atol=1e-8)   # cap centre excluded
    assert r.max() == pytest.approx(R_C, abs=1e-8)
    # no vertex on the grid planes x = 0 / z = 0, so the wall's extreme |x|
    # is R_c cos(pi/N) < R_c = x_max (nothing flush with the box face)
    assert np.abs(pts[:, 0]).max() == pytest.approx(R_C * math.cos(math.pi / 80), abs=1e-8)
    assert np.abs(pts[:, 0]).max() < R_C


def test_top_cap(built):
    _, facts, n, v = built
    pts = v.reshape(-1, 3)
    assert pts[:, 1].max() == pytest.approx(Y_TOP, abs=1e-9)
    assert facts["y_top_m"] == pytest.approx(Y_TOP)
    cap = np.isclose(v[:, :, 1], Y_TOP, atol=1e-9).all(axis=1)
    assert cap.sum() == facts["n_theta"]
    assert np.allclose(n[cap], [0.0, -1.0, 0.0], atol=1e-6)      # inward = -y


def test_rings_avoid_grid_planes(built):
    # No horizontal facet edge lies on a fluid-cell face (y = j dx): keeps
    # MFiX's edge/facet intersections away from exact coincidences.
    _, facts, _, _ = built
    ring_y = np.asarray(facts["ring_y_m"])
    dist = np.abs(ring_y / DX - np.round(ring_y / DX)) * DX
    assert dist.min() > 1.0e-5
    assert facts["band_height_cone_m"] <= BAND and facts["band_height_column_m"] <= BAND
    assert np.isclose(ring_y, H_CONE, atol=1e-12).sum() == 1    # shared junction ring


def test_builder_reports_sagittas(built):
    _, facts, _, _ = built
    assert facts["sagitta_Rc_m"] == pytest.approx(1.943e-5, rel=1e-3) and facts["sagitta_Rc_m"] <= SAGITTA
    assert facts["sagitta_Rt_m"] == pytest.approx(9.18e-7, rel=2e-3)


def test_deterministic_output(built, tmp_path):
    path, _, _, _ = built
    other = tmp_path / "again.stl"
    build_cone_column_stl_apex(R_C, R_T, HALF, H_DOM, other, TOP, SAGITTA, BAND)
    assert other.read_bytes() == path.read_bytes()


@pytest.mark.parametrize("kw", [
    {"R_t_m": 0.0}, {"R_t_m": R_C}, {"cone_half_angle_deg": 0.0}, {"cone_half_angle_deg": 90.0},
    {"H_dom_m": 0.01}, {"top_extend_m": -1e-3}, {"sagitta_max_m": 0.0}, {"band_height_m": 0.0},
])
def test_invalid_inputs_rejected(tmp_path, kw):
    args = dict(R_c_m=R_C, R_t_m=R_T, cone_half_angle_deg=HALF, H_dom_m=H_DOM,
                out_path=tmp_path / "x.stl", top_extend_m=TOP,
                sagitta_max_m=SAGITTA, band_height_m=BAND)
    args.update(kw)
    with pytest.raises(ValueError):
        build_cone_column_stl_apex(**args)


# ---------------- generator wiring ----------------

def _case(case_id: str) -> dict:
    return next(c for c in load_sweep()["cases"] if c["id"] == case_id)


def _render(tmp_path, **kw):
    params = build_C26_3d_params(_case("c02"), **kw)
    out = tmp_path / f"{params['run_name']}.mfx"
    render_case(C26_TEMPLATE, params, out)
    return out.read_text(), params


def _ps_rows(text: str) -> tuple[float, float]:
    v = {k: float(x) for k, x in re.findall(
        r"^\s*(ps_y_[sn])\(1\)\s*=\s*(\S+)", text, flags=re.M)}
    return v["ps_y_s"], v["ps_y_n"]


def test_yaml_default_is_implicit_and_keys_are_documented():
    mesh = _load_yaml(REPO_ROOT / "params" / "geometry.yaml")["spouted_bed_C26_3d"]
    assert mesh["stl_builder"] == "implicit"
    assert mesh["ps_tip_row_offset_cells"] == -1
    assert float(mesh["stl_apex_sagitta_max_m"]) == SAGITTA
    assert int(mesh["stl_apex_band_height_cells"]) == 2


def test_apex_option_moves_the_point_source_to_the_tip_row():
    p = build_C26_3d_params(_case("c02"), stl_builder="parametric_apex")
    dx = p["dx_m"]
    assert p["stl_builder"] == "parametric_apex"
    assert (p["ps_cell_y_s_m"], p["ps_cell_y_n_m"]) == pytest.approx((-dx, 0.0))
    assert p["ps_row_offset_cells"] == -1
    assert p["ps_cell_y_s_m"] <= p["y_apex_m"] < p["ps_cell_y_n_m"]
    assert p["y_apex_m"] == pytest.approx(Y_APEX)
    # same 2 x 2 footprint, same gas, no momentum (mass-only default)
    assert (p["ps_cell_x_w_m"], p["ps_cell_x_e_m"]) == pytest.approx((-dx, dx))
    assert p["ps_v_g_m_s"] == 0.0 and p["ps_momentum"] is False
    # open top face of the tip row = nozzle area -> exit velocity = U_nozzle
    assert p["r_cone_at_ps_top_m"] == pytest.approx(p["R_t_m"])
    assert p["ps_open_area_m2"] == pytest.approx(math.pi * R_T ** 2)
    assert p["U_ps_exit_m_s"] == pytest.approx(p["U_nozzle_m_s"])
    assert p["stl_apex_sagitta_max_m"] == SAGITTA and p["stl_apex_band_height_m"] == pytest.approx(BAND)


def test_implicit_default_is_unchanged():
    p = build_C26_3d_params(_case("c02"))
    dx = p["dx_m"]
    assert p["stl_builder"] == "implicit" and p["y_apex_m"] is None
    assert (p["ps_cell_y_s_m"], p["ps_cell_y_n_m"]) == pytest.approx((0.0, dx))
    assert p["ps_row_offset_cells"] == 0
    with pytest.raises(ValueError):
        build_C26_3d_params(_case("c02"), stl_builder="marching_cubes")


def test_apex_deck_has_the_tip_row_source_and_nothing_else_moves(tmp_path):
    text_apex, _ = _render(tmp_path, stl_builder="parametric_apex", run_name="apex")
    text_impl, _ = _render(tmp_path, run_name="apex")
    assert _ps_rows(text_apex) == pytest.approx((-DX, 0.0))
    assert _ps_rows(text_impl) == pytest.approx((0.0, DX))
    assert "CLOSED-APEX STL" in text_apex and "CLOSED-APEX STL" not in text_impl
    # every other keyword line is identical (the inlet-pressure monitor
    # follows the source cells by design)
    def kw_lines(t):
        return [ln for ln in t.splitlines() if "=" in ln and not ln.strip().startswith("#")]
    diff = [(a, b) for a, b in zip(kw_lines(text_impl), kw_lines(text_apex)) if a != b]
    names = sorted(a.split("=")[0].strip() for a, _ in diff)
    assert names == ["monitor_y_n(1)", "monitor_y_s(1)", "ps_y_n(1)", "ps_y_s(1)"]
    assert len(kw_lines(text_impl)) == len(kw_lines(text_apex))


def test_run_c26_build_stl_apex_facts(tmp_path):
    from src.run_C26_3d import _build_stl
    p = build_C26_3d_params(_case("c02"), stl_builder="parametric_apex")
    facts = _build_stl(p, tmp_path / "geometry_0001.stl")
    assert facts["facet_count"] == 8480 and facts["stl_n_theta"] == 80
    assert re.fullmatch(r"[0-9a-f]{32}", facts["stl_md5"])
    assert facts["stl_open_area_y0_m2"] == pytest.approx(4.4489e-6, rel=1e-4)
    assert facts["stl_tip_volume_m3"] == pytest.approx(3.061e-9, rel=1e-3)

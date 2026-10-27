"""Tests for src.make_stl_c26_fullface (C26 closed-apex vessel with a one-cell
square bore cut into its tip) and its wiring into the C26 case generator
(stl_builder = parametric_fullface)."""

from __future__ import annotations

import json
import math
import re
import sys

import numpy as np
import pytest

from src.make_case import REPO_ROOT, _load_yaml, render_case
from src.make_case_c26 import build_C26_3d_params, load_sweep
from src.make_stl_c26_apex import build_cone_column_stl_apex
from src.make_stl_c26_fullface import (
    _segment_area,
    build_cone_column_stl_fullface,
    fullface_nozzle_geometry,
    inside_fullface_solid,
)
from tests.test_make_stl_c26_apex import _read, _section_polygon_area

C26_TEMPLATE = REPO_ROOT / "cases" / "C26_3d" / "template_dem.mfx.j2"

# C26 geometry (params/geometry.yaml bed_C26 / spouted_bed_C26_3d), full size.
R_C = 0.0252
R_T = 0.5 * 2.38125e-3
HALF = 30.0
H_DOM = 0.21
TOP = 5.0e-3
DX = 2.1e-3
SAGITTA = 2.06e-5
BAND = 2 * DX
DELTA = 0.005 * DX          # stl_fullface_wall_offset_frac
TAN = math.tan(math.radians(HALF))
H_CONE = (R_C - R_T) / TAN
Y_TOP = H_DOM + TOP
A = 0.5 * DX + DELTA
Y_FLOOR = -DX - DELTA
N = 80


def _inside(x, y, z):
    return inside_fullface_solid(x, y, z, R_C, R_T, HALF, Y_TOP, A, Y_FLOOR)


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    out = tmp_path_factory.mktemp("c26fullface") / "geometry_0001.stl"
    path, facts = build_cone_column_stl_fullface(
        R_c_m=R_C, R_t_m=R_T, cone_half_angle_deg=HALF, H_dom_m=H_DOM, out_path=out,
        top_extend_m=TOP, sagitta_max_m=SAGITTA, band_height_m=BAND, cell_half_m=0.5 * DX,
        nozzle_y_s_m=-DX, wall_offset_m=DELTA,
    )
    n, v = _read(path)
    return path, facts, n, v


# ---------------- nozzle geometry ----------------

def test_nozzle_geometry_hand_values():
    g = fullface_nozzle_geometry(R_T, HALF, 0.5 * DX, -DX, DELTA)
    assert g["half_width_m"] == pytest.approx(1.0605e-3)
    assert g["y_floor_m"] == pytest.approx(-2.1105e-3)
    # the cone meets the bore where R_t + y tan30 = a (mid-side) or a sqrt 2
    # (corners): (1.0605 - 1.1906) / 0.5774 = -0.2254 mm, (1.4998 - 1.1906) /
    # 0.5774 = +0.5355 mm
    assert g["y_curve_min_m"] == pytest.approx(-0.2254e-3, rel=1e-3)
    assert g["y_curve_max_m"] == pytest.approx(0.5355e-3, rel=1e-3)
    assert g["nozzle_face_area_m2"] == pytest.approx(4.41e-6)
    assert g["nozzle_face_area_m2"] / (math.pi * R_T ** 2) == pytest.approx(0.9902, rel=1e-3)
    assert g["sliver_fraction"] == pytest.approx(0.005)
    # cone past the bore's side in one face neighbour: ~0.009 mm^3 = 0.09 %
    assert g["protrusion_fraction"] == pytest.approx(9.14e-4, rel=1e-2)
    assert g["neighbour_fluid_fraction"] == pytest.approx(0.00591, rel=1e-2)
    assert g["neighbour_fluid_fraction"] < 0.01


def test_segment_area_hand_values():
    assert _segment_area(1.0, 0.0) == pytest.approx(math.pi / 2)
    assert _segment_area(1.0, 1.0) == 0.0 and _segment_area(1.0, 2.0) == 0.0
    # R_t past a at y = 0: 1.41759 acos(0.89071) - 1.0605 sqrt(0.29293)
    # = 0.66903 - 0.57397 = 0.0950 mm^2
    assert _segment_area(R_T, A) == pytest.approx(9.50e-8, rel=2e-3)


@pytest.mark.parametrize("kw", [
    {"wall_offset_m": 1.0e-6},              # at TOL_STL: crossings dropped
    {"wall_offset_m": 0.0095 * DX},         # neighbours reach TOL_SMALL_CELL: not blocked
    {"R_t_m": 1.6e-3},                      # cone wider than the bore corners at y = 0
    {"R_t_m": 1.25e-3},                     # apex below the bore floor
    {"cone_half_angle_deg": 90.0},
])
def test_nozzle_geometry_rejects_bad_inputs(kw):
    args = dict(R_t_m=R_T, cone_half_angle_deg=HALF, cell_half_m=0.5 * DX, nozzle_y_s_m=-DX,
                wall_offset_m=DELTA)
    args.update(kw)
    with pytest.raises(ValueError):
        fullface_nozzle_geometry(**args)


# ---------------- the shell ----------------

def test_header_and_facet_count(built):
    path, facts, n, _ = built
    lines = path.read_text().strip().splitlines()
    assert lines[0] == "solid cone_column_C26_fullface" and lines[-1].startswith("endsolid ")
    n_bands = facts["n_bands_cone"] + facts["n_bands_column"] - 1   # bands above ring 1
    # collar 2N + 4, bore walls 4 x 2 (N/4 + 1), floor N + 4, bands 2N each, cap N
    expected = (2 * N + 4) + 8 * (N // 4 + 1) + (N + 4) + 2 * N * n_bands + N
    assert facts["facet_count"] == len(n) == expected == 8816
    assert (facts["n_theta"], facts["n_bands_cone"], facts["n_bands_column"]) == (N, 11, 42)


def test_watertight_genus_zero(built):
    _, _, _, v = built
    uniq, inv = np.unique(v.reshape(-1, 3), axis=0, return_inverse=True)
    tri = inv.reshape(-1, 3)
    edges = np.concatenate([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]])
    edges.sort(axis=1)
    _, counts = np.unique(edges, axis=0, return_counts=True)
    assert counts.min() == 2 and counts.max() == 2
    assert len(uniq) - len(counts) + len(tri) == 2


def test_normals_unit_consistent_and_inward(built):
    _, _, n, v = built
    cross = np.cross(v[:, 1] - v[:, 0], v[:, 2] - v[:, 0])
    area2 = np.linalg.norm(cross, axis=1)
    assert area2.min() >= 1e-12
    assert np.allclose(np.linalg.norm(n, axis=1), 1.0, atol=1e-5)
    assert np.allclose(cross / area2[:, None], n, atol=1e-5)
    c = v.mean(axis=1)
    eps = 1.0e-4
    p_in, p_out = c + eps * n, c - eps * n
    assert _inside(p_in[:, 0], p_in[:, 1], p_in[:, 2]).all()
    assert not _inside(p_out[:, 0], p_out[:, 1], p_out[:, 2]).any()


def test_no_ledge_only_floor_and_cap_are_horizontal(built):
    _, _, n, v = built
    flat = np.abs(n[:, 1]) > 1.0 - 1e-6
    y = v[flat][:, :, 1]
    floor = np.isclose(y, Y_FLOOR, atol=1e-9).all(axis=1)
    cap = np.isclose(y, Y_TOP, atol=1e-9).all(axis=1)
    assert (floor | cap).all()
    assert floor.sum() == N + 4 and cap.sum() == N
    assert np.allclose(n[flat][floor], [0.0, 1.0, 0.0], atol=1e-6)


def test_identical_to_the_apex_stl_above_the_collar(built, tmp_path):
    _, facts, _, v = built
    _, apex_facts = build_cone_column_stl_apex(R_C, R_T, HALF, H_DOM, tmp_path / "apex.stl",
                                               TOP, SAGITTA, BAND)
    _, va = _read(tmp_path / "apex.stl")
    y0 = facts["y_collar_top_m"]
    assert y0 == pytest.approx(apex_facts["ring_y_m"][0]) and y0 == pytest.approx(1.906e-3, rel=1e-3)

    def above(vv):
        keep = (vv[:, :, 1] >= y0 - 1e-12).all(axis=1)
        return sorted(tuple(np.round(t, 9).ravel()) for t in vv[keep])
    assert above(v) == above(va)
    assert len(above(v)) == 8816 - (2 * N + 4) - 8 * (N // 4 + 1) - (N + 4)


def test_curve_vertices_lie_on_the_cone_and_the_bore(built):
    _, facts, _, v = built
    pts = np.unique(v.reshape(-1, 3), axis=0)
    mid = pts[(pts[:, 1] > Y_FLOOR + 1e-9) & (pts[:, 1] < facts["y_collar_top_m"] - 1e-9)]
    assert len(mid) == N + 4
    assert np.allclose(np.maximum(np.abs(mid[:, 0]), np.abs(mid[:, 2])), A, atol=1e-9)
    assert np.allclose(np.hypot(mid[:, 0], mid[:, 2]), R_T + mid[:, 1] * TAN, atol=1e-8)
    assert mid[:, 1].min() == pytest.approx(facts["y_curve_min_m"], abs=2e-6)
    assert mid[:, 1].max() == pytest.approx(facts["y_curve_max_m"], abs=1e-9)
    # no curve vertex within 10 um of the grid plane y = 0
    assert np.abs(mid[:, 1]).min() > 1.0e-5


def test_source_cell_is_a_whole_fluid_cell():
    # All 8 corner nodes and every point of the 12 edges lie in the fluid, so
    # MFiX keeps it a standard (uncut) cell: PS volume = dx^3.
    s = np.linspace(0.0, 1.0, 201)
    h = 0.5 * DX
    nodes = [(sx * h, y, sz * h) for sx in (-1, 1) for sz in (-1, 1) for y in (-DX, 0.0)]
    x, y, z = np.array(nodes).T
    assert _inside(x, y, z).all()
    for sx in (-1, 1):
        for sz in (-1, 1):
            assert _inside(sx * h + 0 * s, -DX + DX * s, sz * h + 0 * s).all()
            assert _inside(-h + DX * s, -DX + 0 * s, sz * h + 0 * s).all()
            assert _inside(-h + DX * s, 0.0 * s, sz * h + 0 * s).all()
            assert _inside(sx * h + 0 * s, -DX + 0 * s, -h + DX * s).all()
            assert _inside(sx * h + 0 * s, 0.0 * s, -h + DX * s).all()


def _fluid_fraction_x(x0, y0, z0, n=400):
    """Fluid share of the cell [x0, x0+dx] x [y0, y0+dx] x [z0, z0+dx] with
    x0 >= dx/2: exact fluid length along x on an n x n grid of (y, z) lines
    (the bore and the cone are both intervals in x on such a line)."""
    yy = y0 + DX * (np.arange(n) + 0.5) / n
    zz = z0 + DX * (np.arange(n) + 0.5) / n
    Y, Z = np.meshgrid(yy, zz, indexing="ij")
    hi = np.zeros_like(Y)
    in_bore = (Y > Y_FLOOR) & (np.abs(Z) < A)
    hi = np.where(in_bore, A, hi)
    r = R_T + Y * TAN
    cone_x = np.sqrt(np.clip(r * r - Z * Z, 0.0, None))
    hi = np.where((Y > -R_T / TAN) & (cone_x > hi), cone_x, hi)
    length = np.clip(np.minimum(hi, x0 + DX) - x0, 0.0, None)
    return float(length.mean() / DX)


def test_source_row_neighbours_and_cells_below_stay_under_tol_small_cell():
    # Face neighbour (+x): the 0.5 % bore sliver plus the cone's 0.09 %
    # protrusion; diagonal neighbour: a delta^2 column; cells below: the
    # floor sliver.  All under TOL_SMALL_CELL = 1 %, so MFiX blocks them.
    h = 0.5 * DX
    g = fullface_nozzle_geometry(R_T, HALF, h, -DX, DELTA)
    face = _fluid_fraction_x(h, -DX, -h)
    assert face == pytest.approx(g["neighbour_fluid_fraction"], rel=0.02)
    assert face < 0.01
    assert _fluid_fraction_x(h, -DX, h) < 1e-4                 # diagonal
    assert _fluid_fraction_x(h, -2 * DX, -h) < 1e-4             # below the face neighbour
    # below the source cell: the floor sliver, delta / dx
    yy = np.linspace(-2 * DX, -DX, 20001)
    assert np.mean(_inside(0.0 * yy, yy, 0.0 * yy)) == pytest.approx(DELTA / DX, rel=0.01)


def test_row_above_source_centre_cell_is_whole_fluid():
    h = 0.5 * DX
    gg = np.linspace(-h, h, 41)
    X, Y, Z = np.meshgrid(gg, np.linspace(1e-12, DX, 41), gg, indexing="ij")
    assert _inside(X.ravel(), Y.ravel(), Z.ravel()).all()


def test_sections(built):
    _, _, _, v = built
    # source row: the bore only; y = 0: bore plus the cone's four segments
    assert _section_polygon_area(v, -0.5 * DX) == pytest.approx((2 * A) ** 2, rel=1e-6)
    assert _section_polygon_area(v, 0.0) == pytest.approx(
        (2 * A) ** 2 + 4 * _segment_area(R_T, A), rel=5e-3)
    # above the corners the section is the cone's polygon
    y = 1.0e-3
    assert _section_polygon_area(v, y) == pytest.approx(math.pi * (R_T + y * TAN) ** 2, rel=3e-3)


def test_deterministic_output(built, tmp_path):
    path, _, _, _ = built
    other = tmp_path / "again.stl"
    build_cone_column_stl_fullface(R_C, R_T, HALF, H_DOM, other, TOP, SAGITTA, BAND,
                                   0.5 * DX, -DX, DELTA)
    assert other.read_bytes() == path.read_bytes()


# ---------------- generator wiring ----------------

def _case(case_id: str) -> dict:
    return next(c for c in load_sweep()["cases"] if c["id"] == case_id)


def test_yaml_keys():
    mesh = _load_yaml(REPO_ROOT / "params" / "geometry.yaml")["spouted_bed_C26_3d"]
    assert mesh["stl_builder"] == "implicit"           # default unchanged
    assert float(mesh["stl_fullface_wall_offset_frac"]) == 0.005
    assert "stl_fullface_ledge_clearance_frac" not in mesh
    assert mesh["ps_momentum"] is False                 # production default unchanged


def test_fullface_params():
    p = build_C26_3d_params(_case("c02"), stl_builder="parametric_fullface")
    dx = p["dx_m"]
    assert p["axis_on_cell_face"] is False and p["imax"] == p["kmax"] == 25
    assert (p["x_min_m"], p["x_max_m"]) == pytest.approx((-12.5 * dx, 12.5 * dx))
    assert (p["ps_cell_x_w_m"], p["ps_cell_x_e_m"]) == pytest.approx((-0.5 * dx, 0.5 * dx))
    assert (p["ps_cell_y_s_m"], p["ps_cell_y_n_m"]) == pytest.approx((-dx, 0.0), abs=1e-15)
    assert p["ps_region_x_w_m"] == pytest.approx(-0.5 * dx) and p["ps_point"] is False
    assert p["ps_open_area_m2"] == pytest.approx(dx ** 2)
    assert p["U_ps_exit_m_s"] == pytest.approx(p["Q_actual_m3_s"] / dx ** 2)
    assert p["U_ps_exit_m_s"] == pytest.approx(45.35, rel=1e-3)
    assert p["ps_v_g_m_s"] == 0.0                        # mass-only default
    assert p["nozzle_half_width_m"] == pytest.approx(A)
    assert p["y_bore_curve_min_m"] == pytest.approx(-0.2254e-3, rel=1e-3)
    assert p["nozzle_neighbour_fluid_fraction"] < 0.01
    assert p["r_cone_at_ps_top_m"] is None
    pm = build_C26_3d_params(_case("c02"), stl_builder="parametric_fullface", ps_momentum=True)
    assert pm["ps_v_g_m_s"] == pytest.approx(pm["U_ps_exit_m_s"])
    pa = build_C26_3d_params(_case("c02"), stl_builder="parametric_apex", axis_on_face=True,
                             ps_momentum=True)
    assert pa["ps_v_g_m_s"] == pytest.approx(pa["U_nozzle_m_s"])
    for k in ("dt_s", "dt_max_s", "ps_massflow_g_kg_s", "H_ic_m", "n_parcels_est"):
        assert pm[k] == pytest.approx(pa[k])


def test_fullface_rejects_axis_on_face_and_point():
    with pytest.raises(ValueError):
        build_C26_3d_params(_case("c02"), stl_builder="parametric_fullface", axis_on_face=True)
    with pytest.raises(ValueError):
        build_C26_3d_params(_case("c02"), stl_builder="parametric_fullface", ps_point=True)


def _kw(text: str) -> dict:
    out = {}
    for ln in text.splitlines():
        s = ln.strip()
        if "=" in s and not s.startswith("#") and not re.match(r"d[xz]\(\d+\)", s):
            k, v = s.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def test_fullface_deck_differs_from_the_apex_deck_only_in_x_z_and_ps_v_g(tmp_path):
    def render(**kw):
        p = build_C26_3d_params(_case("c02"), run_name="X", ps_momentum=True, **kw)
        out = tmp_path / "x.mfx"
        render_case(C26_TEMPLATE, p, out)
        return out.read_text()
    full = render(stl_builder="parametric_fullface")
    apex = render(stl_builder="parametric_apex", axis_on_face=True)
    assert "FULL-FACE NOZZLE STL" in full and "FULL-FACE NOZZLE STL" not in apex
    assert "(no ledge)" in full
    a, b = _kw(apex), _kw(full)
    assert set(a) == set(b)
    changed = {k for k in a if a[k] != b[k]}
    allowed = re.compile(r"^(imax|kmax|[xz]_(min|max)|(ic|bc|vtk|ps|monitor)_[xz]_[webt]\(\d+\)"
                         r"|ps_v_g\(1\))$")
    assert changed and all(allowed.match(k) for k in changed), sorted(changed)
    assert "ps_v_g(1)" in changed and float(b["ps_v_g(1)"]) == pytest.approx(45.35, rel=1e-3)
    for k in ("ps_y_s(1)", "ps_y_n(1)", "ps_massflow_g(1)", "dt", "tstop"):
        assert a[k] == b[k]
    assert len(re.findall(r"^\s*dx\(\d+\)", full, flags=re.M)) == 25


def test_run_c26_build_stl_fullface_facts(tmp_path):
    from src.run_C26_3d import _build_stl
    p = build_C26_3d_params(_case("c02"), stl_builder="parametric_fullface")
    facts = _build_stl(p, tmp_path / "geometry_0001.stl")
    assert facts["facet_count"] == 8816 and facts["stl_n_theta"] == N
    assert facts["stl_nozzle_face_area_m2"] == pytest.approx(DX ** 2)
    assert facts["stl_neighbour_fluid_fraction"] == pytest.approx(0.00591, rel=1e-2)
    assert re.fullmatch(r"[0-9a-f]{32}", facts["stl_md5"])


def test_run_c26_cli_fullface_odd_box_and_gridmap(monkeypatch, capsys, tmp_path):
    from src import run_C26_3d
    monkeypatch.setattr(sys, "argv", ["run_C26_3d", "--case", "c02", "--stl-builder",
                                      "parametric_fullface", "--ps-momentum", "--dry-run",
                                      "--out-root", str(tmp_path)])
    run_C26_3d.main()
    m = json.loads(capsys.readouterr().out)
    assert m["imax"] == 25 and m["axis_on_cell_face"] is False
    assert m["gridmap"] == {"i": [12, 13], "j": [104], "k": [12, 13]}
    monkeypatch.setattr(sys, "argv", ["run_C26_3d", "--case", "c02", "--stl-builder",
                                      "parametric_fullface", "--axis-on-face", "--dry-run"])
    with pytest.raises(SystemExit):
        run_C26_3d.main()

"""Tests for the PIV rig model: the C26 layout with the measured values of
data/piv_rig (params *_piv_rig blocks, rim pocket STL, --rig piv_rig)."""

from __future__ import annotations

import csv
import json
import math
import re
import sys

import numpy as np
import pytest

from src.make_case import REPO_ROOT, _load_yaml, render_case
from src.make_case_c26 import build_C26_3d_params, default_gridmap, load_sweep
from src.make_stl_c26_fullface import build_cone_column_stl_fullface, inside_fullface_solid
from tests.test_make_case_c26 import C26_NEW_KEYWORDS, V3D_DEM_TEMPLATE, _keywords
from tests.test_make_stl_c26_apex import _read

TEMPLATE = REPO_ROOT / "cases" / "C26_3d" / "template_dem.mfx.j2"
APPARATUS = REPO_ROOT / "data" / "piv_rig" / "spouted-bed-piv-experiment-log-apparatus.csv"
PARTICLES = REPO_ROOT / "data" / "piv_rig" / "spouted-bed-piv-experiment-log-particles.csv"

# Ideal-gas N2 at 25 C, 1 atm: rho = p M / (R T).
RHO_N2_25C = 101325.0 * 28.0134e-3 / (8.314462618 * 298.15)    # = 1.1451 kg/m^3
TAN30 = math.tan(math.radians(30.0))


def _case(case_id: str) -> dict:
    return next(c for c in load_sweep("piv_rig")["cases"] if c["id"] == case_id)


def _build(case_id: str = "r01", **kw) -> dict:
    return build_C26_3d_params(_case(case_id), rig="piv_rig", **kw)


@pytest.fixture(scope="module")
def r01() -> dict:
    return _build()


@pytest.fixture(scope="module")
def apparatus() -> dict:
    """Apparatus sheet: measurement ID -> value used [mm]."""
    with APPARATUS.open(newline="") as f:
        rows = list(csv.reader(f))
    return {r[0]: float(r[6]) for r in rows if re.fullmatch(r"[FTN]\d+", r[0])}


@pytest.fixture(scope="module")
def particles_sheet() -> dict:
    """Particles sheet: first column -> second column, and the raw samples."""
    with PARTICLES.open(newline="") as f:
        rows = list(csv.reader(f))
    return {
        "derived": {r[0]: r[1] for r in rows if r[0] and len(r) > 1},
        "samples": [(float(r[1]), int(r[2])) for r in rows if r[0].startswith("Sample ")],
        "balls": [float(r[5]) for r in rows if r[4].startswith("Ball ")],
    }


# ---------------- params against the measurement sheets ----------------

def test_geometry_yaml_matches_apparatus_sheet(apparatus):
    bed = _load_yaml(REPO_ROOT / "params" / "geometry.yaml")["bed_piv_rig"]
    assert bed["D_c_m"] * 1e3 == pytest.approx(apparatus["T1"])
    assert bed["D_nozzle_m"] * 1e3 == pytest.approx(apparatus["N1"])
    assert bed["D_throat_m"] * 1e3 == pytest.approx(apparatus["N1"])
    assert bed["D_rim_m"] * 1e3 == pytest.approx(apparatus["F8"])
    assert bed["y_rim_top_m"] * 1e3 == pytest.approx(apparatus["F11"])
    assert bed["L_tube_m"] * 1e3 == pytest.approx(apparatus["T3"])
    assert bed["cone_angle_deg"] == 60        # user decision: 60 deg exact


def test_particles_yaml_matches_particles_sheet(particles_sheet):
    feed = _load_yaml(REPO_ROOT / "params" / "particles.yaml")["case_piv_rig_ss"]
    derived = particles_sheet["derived"]
    assert feed["rho_p_kg_per_m3"] == pytest.approx(float(derived["Density used"]))
    assert feed["d_p_mean_m"] * 1e3 == pytest.approx(float(derived["Diameter used"]))
    # the sheet's derived numbers follow from its raw samples
    d = np.array(particles_sheet["balls"])
    assert d.mean() == pytest.approx(0.822, abs=1e-6)
    assert math.sqrt(feed["d_p_var_m2"]) * 1e6 == pytest.approx(d.std(ddof=1) * 1e3, abs=0.05)
    m_g, n = map(sum, zip(*particles_sheet["samples"]))
    mg_per_ball = 1e3 * m_g / n
    assert mg_per_ball == pytest.approx(2.2059, abs=1e-4)
    # 27,200 balls in the charge make the 60 g
    n_charge = float(derived["Balls in charge, N"].replace(",", ""))
    assert n_charge * mg_per_ball * 1e-3 == pytest.approx(feed["inventory_g"], rel=1e-3)


def test_sweep_flows_and_weight():
    sweep = load_sweep("piv_rig")
    assert [(c["id"], c["flow_slpm"]) for c in sweep["cases"]] == [
        ("r01", 24), ("r02", 26.5), ("r03", 30)]
    assert {c["inventory_g"] for c in sweep["cases"]} == {60}
    assert sweep["base"]["stat_wt"] == 1
    assert sweep["base"]["gas_case"] == "n2_25C"


# ---------------- hand values ----------------

def test_gas_is_nitrogen_at_25c(r01):
    assert r01["gas_short"] == "N2"
    assert r01["rho_g_kg_m3"] == pytest.approx(RHO_N2_25C, rel=2e-3)
    assert r01["rho_ref_kg_m3"] == pytest.approx(RHO_N2_25C, rel=2e-3)
    # N2 at 298 K: 1.78e-5 Pa s (Poling et al. tables) within 2 %
    assert r01["mu_g_pa_s"] == pytest.approx(1.78e-5, rel=0.02)


@pytest.mark.parametrize("cid,slpm,u_nozzle,u_face", [
    # Q = slpm / 60000; orifice pi/4 (2.7 mm)^2 = 5.7256 mm^2; dx^2 = 5.8980 mm^2
    ("r01", 24.0, 69.86, 67.82),
    ("r02", 26.5, 77.14, 74.88),
    ("r03", 30.0, 87.33, 84.78),
])
def test_jet_velocities_hand_values(cid, slpm, u_nozzle, u_face):
    p = _build(cid)
    assert p["U_nozzle_m_s"] == pytest.approx(u_nozzle, abs=0.02)
    assert p["ps_v_g_m_s"] == pytest.approx(u_face, abs=0.02)
    assert p["ps_massflow_g_kg_s"] == pytest.approx(slpm / 60000.0 * p["rho_ref_kg_m3"], rel=1e-9)
    # the paper's step stays below recipe Rule 9 at every flow
    assert p["dt_max_s"] == pytest.approx(1.0e-5)


def test_mesh_hand_values(r01, apparatus):
    dx = 51.0e-3 / 21
    assert r01["dx_m"] == pytest.approx(dx)
    assert (r01["imax"], r01["jmax"], r01["kmax"]) == (23, 150, 23)
    assert r01["x_max_m"] == pytest.approx(11.5 * dx)
    assert r01["x_max_m"] > r01["R_rim_m"]                      # the rim fits the box
    # the source face is the orifice to within 3 %, as C26's is to within 1 %
    assert r01["ps_open_area_m2"] / (math.pi / 4 * 2.7e-3 ** 2) == pytest.approx(1.030, abs=1e-3)
    assert 2.0 < r01["dx_m"] / r01["d_p_m"] < 3.0               # C26 2.2.1 rule
    # the outlet is the whole-row height closest to the tube top
    tube_top = (apparatus["F11"] + apparatus["T3"]) * 1e-3
    assert abs(r01["H_dom_m"] - tube_top) < 0.5 * dx
    assert r01["H_dom_m"] == pytest.approx(146 * dx)


def test_cone_and_rim_heights(r01):
    assert r01["R_c_m"] == pytest.approx(25.5e-3)
    assert r01["R_t_m"] == pytest.approx(1.35e-3)
    assert r01["cone_half_angle_deg"] == 30.0
    assert r01["h_cone_m"] == pytest.approx((27.5e-3 - 1.35e-3) / TAN30)
    assert r01["h_cone_m"] == pytest.approx(45.29e-3, abs=0.01e-3)
    assert r01["y_rim_top_m"] == 0.050


def test_bed_stays_in_the_cone(r01):
    # 60 g / (7585 * 0.60) = 13,184 mm^3 (Particles sheet); its 31.4 mm uses the
    # 29.8 deg half-angle, so 60 deg exact is slightly shorter.
    assert r01["H_static_m"] == pytest.approx(31.2e-3, abs=0.1e-3)
    assert r01["H_static_m"] < r01["H_ic_m"] < r01["h_cone_m"]


def test_ic_region_keeps_the_charge_on_a_grid_row(r01):
    # MFiX snaps ic_y_n to the nearest row and seeds ic_ep_s * V / V_p parcels.
    # Hand values: 0.30 gives 39.94 mm = 16.44 rows -> 16 rows = 38.857 mm;
    # r_top = 1.35 + 38.857 tan 30 = 23.784 mm; V = 24,399 mm^3 ->
    # ic_ep_s = (60 g / 7585 kg/m^3) / V = 0.3242.
    assert r01["H_ic_m"] / r01["dx_m"] == pytest.approx(16.0)
    assert r01["H_ic_m"] == pytest.approx(38.857e-3, abs=1e-6)
    assert r01["eps_s_ic"] == pytest.approx(0.3242, abs=1e-4)
    assert r01["ep_g_ic"] == pytest.approx(1.0 - r01["eps_s_ic"])


def test_single_cell_inlet_is_the_default(r01):
    assert r01["stl_builder"] == "parametric_fullface"
    assert r01["axis_on_cell_face"] is False
    assert r01["ps_footprint_cells"] == 1 and r01["ps_momentum"] is True
    dx = r01["dx_m"]
    assert (r01["ps_cell_x_w_m"], r01["ps_cell_x_e_m"]) == pytest.approx((-dx / 2, dx / 2))
    assert (r01["ps_cell_y_s_m"], r01["ps_cell_y_n_m"]) == pytest.approx((-dx, 0.0), abs=1e-15)
    assert r01["ps_open_area_m2"] == pytest.approx(dx ** 2)
    assert r01["nozzle_neighbour_fluid_fraction"] < 0.01         # MFiX blocks them


def test_parcels_at_weight_one(r01):
    assert r01["cgp_stat_wt"] == 1.0
    assert r01["d_parcel_m"] == pytest.approx(822.0e-6)
    # 60 g / (7585 * pi/6 * E[d^3]) = the sheet's 27,200 balls
    assert r01["n_parcels_est"] == pytest.approx(27200, rel=1e-3)
    # the same parcel mass as a C26 YSZ W = 10 parcel, so the same contact time
    c26 = build_C26_3d_params(next(c for c in load_sweep()["cases"] if c["id"] == "c01"))
    assert r01["dt_rayleigh_s"] == pytest.approx(c26["dt_rayleigh_s"], rel=0.01)


def test_weight_floor_is_per_rig():
    assert build_C26_3d_params({**_case("r01"), "stat_wt": 1}, rig="piv_rig")["cgp_stat_wt"] == 1.0
    c01 = next(c for c in load_sweep()["cases"] if c["id"] == "c01")
    with pytest.raises(ValueError, match="stat_wt must be >= 10"):
        build_C26_3d_params({**c01, "stat_wt": 1})


def test_default_gridmap_for_the_rig_box():
    assert default_gridmap(rig="piv_rig") == ([11, 12], [], [11, 12])


def test_unknown_rig_rejected():
    with pytest.raises(ValueError, match="rig must be one of"):
        build_C26_3d_params(_case("r01"), rig="nope")


def test_rim_needs_the_fullface_builder():
    with pytest.raises(ValueError, match="rim pocket"):
        _build(stl_builder="parametric_apex")


def test_c26_keeps_its_blocks():
    c02 = next(c for c in load_sweep()["cases"] if c["id"] == "c02")
    p = build_C26_3d_params(c02)
    assert p["rig"] == "C26" and p["R_rim_m"] is None and p["gas_short"] == "Ar"
    assert p["run_name"] == "C26_c02"
    assert p["description"].endswith("12 slpm Ar (Conry et al. 2026 model)")


# ---------------- rendered deck ----------------

def _render(tmp_path, **kw) -> tuple[str, dict]:
    p = _build(**kw)
    out = tmp_path / f"{p['run_name']}.mfx"
    render_case(TEMPLATE, p, out)
    return out.read_text(), p


def test_deck_carries_the_rig_values(tmp_path):
    text, p = _render(tmp_path)
    for kw, val in (("ro_s0(1)", 7585.0), ("d_p0(1)", 822.0e-6), ("cgp_stat_wt(1)", 1.0),
                    ("ro_g0", p["rho_g_kg_m3"]), ("mu_g0", p["mu_g_pa_s"]),
                    ("ps_massflow_g(1)", p["ps_massflow_g_kg_s"]),
                    ("ps_v_g(1)", p["ps_v_g_m_s"]),
                    ("ic_psd_mean_dp(2,1)", 822.0e-6), ("ic_psd_stdev(2,1)", 10.3e-6),
                    ("imax", 23), ("jmax", 150), ("kmax", 23),
                    ("y_max", p["H_dom_m"])):
        m = re.search(rf"^\s*{re.escape(kw)}\s*=\s*(\S+)", text, re.M)
        assert m and float(m.group(1)) == pytest.approx(val, rel=1e-5), kw
    assert re.search(r"solids_model\(1\)\s*=\s*'CGP'", text)
    assert re.search(r"drag_type\s*=\s*'TGS'", text)
    assert "24.0 slpm N2" in text and "slpm Ar" not in text
    assert "RIM POCKET" in text and "sweep_C26" not in text


def test_deck_keywords_are_the_verified_c26_set(tmp_path):
    """CLAUDE.md rule 1: the rig deck introduces no keyword."""
    text, _ = _render(tmp_path)
    allowed = _keywords(V3D_DEM_TEMPLATE.read_text()) | C26_NEW_KEYWORDS
    assert not (_keywords(text) - allowed)
    c02 = next(c for c in load_sweep()["cases"] if c["id"] == "c02")
    out = tmp_path / "c26.mfx"
    render_case(TEMPLATE, build_C26_3d_params(c02, stl_builder="parametric_fullface",
                                              ps_momentum=True), out)
    assert _keywords(text) == _keywords(out.read_text())


# ---------------- STL with the rim pocket ----------------

@pytest.fixture(scope="module")
def stl(tmp_path_factory, r01):
    path = tmp_path_factory.mktemp("rig") / "rig.stl"
    _, facts = build_cone_column_stl_fullface(
        R_c_m=r01["R_c_m"], R_t_m=r01["R_t_m"], cone_half_angle_deg=30.0,
        H_dom_m=r01["H_dom_m"], out_path=path, top_extend_m=r01["top_extend_m"],
        sagitta_max_m=r01["stl_apex_sagitta_max_m"], band_height_m=r01["stl_apex_band_height_m"],
        cell_half_m=0.5 * r01["dx_m"], nozzle_y_s_m=r01["ps_cell_y_s_m"],
        wall_offset_m=r01["stl_fullface_wall_offset_m"],
        R_rim_m=r01["R_rim_m"], y_rim_top_m=r01["y_rim_top_m"],
    )
    return _read(path)[1], facts


def _radii_at(tris: np.ndarray, y: float) -> np.ndarray:
    """Radii of the facet vertices at height y."""
    v = tris.reshape(-1, 3)
    return np.hypot(v[:, 0], v[:, 2])[np.isclose(v[:, 1], y)]


def test_rim_stl_is_watertight(stl):
    tris = stl[0]
    v = np.round(tris.reshape(-1, 3) * 1e9).astype(np.int64)
    _, idx = np.unique(v, axis=0, return_inverse=True)
    idx = idx.reshape(-1, 3)
    e = np.sort(np.concatenate([idx[:, [0, 1]], idx[:, [1, 2]], idx[:, [2, 0]]]), axis=1)
    _, counts = np.unique(e, axis=0, return_counts=True)
    assert np.all(counts == 2)


def test_rim_stl_profile(stl, r01):
    tris, facts = stl
    assert facts["n_theta"] == 80 and facts["sagitta_Rrim_m"] <= r01["stl_apex_sagitta_max_m"]
    assert facts["h_cone_m"] == pytest.approx(r01["h_cone_m"])
    # the rim band's bottom ring is at R_rim, the step has R_rim and R_c rings
    assert np.allclose(_radii_at(tris, facts["h_cone_m"]), r01["R_rim_m"])
    at_step = _radii_at(tris, 0.050)
    assert set(np.round(at_step, 7)) == {round(r01["R_rim_m"], 7), round(r01["R_c_m"], 7)}
    # every vertex lies on the true cone, band or column
    v = tris.reshape(-1, 3)
    y, r = v[:, 1], np.hypot(v[:, 0], v[:, 2])
    cone = (y > 1.9e-3) & (y < facts["h_cone_m"] - 1e-9)
    assert np.allclose(r[cone], r01["R_t_m"] + y[cone] * TAN30, atol=1e-9)
    column = (y > 0.050 + 1e-9) & (r > 0.0)        # r = 0: the top cap's fan centre
    assert np.allclose(r[column], r01["R_c_m"], atol=1e-9)


def test_rim_step_is_the_only_new_horizontal_surface_and_faces_down(stl, r01):
    tris = stl[0]
    n = np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0])
    n /= np.linalg.norm(n, axis=1, keepdims=True)
    flat = np.abs(n[:, 1]) > 1.0 - 1e-9
    y = tris[flat, 0, 1]
    assert set(np.round(y, 7)) == {round(r01["y_nozzle_floor_m"], 7), 0.05,
                                   round(r01["H_dom_m"] + r01["top_extend_m"], 7)}
    assert np.all(n[flat][np.isclose(y, 0.05), 1] < 0.0)         # parcels cannot rest on it


def test_inside_solid_with_the_rim(r01):
    def inside(r, y):
        return bool(inside_fullface_solid(
            np.array([r]), np.array([y]), np.array([0.0]), r01["R_c_m"], r01["R_t_m"], 30.0,
            0.3596, r01["nozzle_half_width_m"], r01["y_nozzle_floor_m"],
            r01["R_rim_m"], r01["y_rim_top_m"])[0])
    assert inside(26.0e-3, 47.0e-3)            # in the pocket, wider than the tube
    assert not inside(26.0e-3, 51.0e-3)        # above the step, inside the tube wall
    assert inside(25.0e-3, 51.0e-3)            # in the tube
    assert not inside(27.6e-3, 47.0e-3)        # outside the rim
    assert not inside(20.0e-3, 30.0e-3)        # outside the cone (r_wall = 18.67 mm)


@pytest.mark.parametrize("kw,match", [
    ({"R_rim_m": 27.5e-3}, "go together"),
    ({"R_rim_m": 25.0e-3, "y_rim_top_m": 0.05}, "must be >"),
    ({"R_rim_m": 27.5e-3, "y_rim_top_m": 0.04}, "must lie between"),
])
def test_rim_arguments_rejected(tmp_path, r01, kw, match):
    with pytest.raises(ValueError, match=match):
        build_cone_column_stl_fullface(
            R_c_m=r01["R_c_m"], R_t_m=r01["R_t_m"], cone_half_angle_deg=30.0,
            H_dom_m=r01["H_dom_m"], out_path=tmp_path / "x.stl", top_extend_m=5.0e-3,
            sagitta_max_m=2.2e-5, band_height_m=r01["stl_apex_band_height_m"],
            cell_half_m=0.5 * r01["dx_m"], nozzle_y_s_m=-r01["dx_m"],
            wall_offset_m=r01["stl_fullface_wall_offset_m"], **kw)


# ---------------- CLI ----------------

def test_cli_writes_the_rig_cases(monkeypatch, capsys, tmp_path):
    from src import run_C26_3d
    monkeypatch.setattr(sys, "argv", [
        "run_C26_3d", "--rig", "piv_rig", "--case", "r02", "--out-root", str(tmp_path)])
    run_C26_3d.main()
    run_dir = tmp_path / "r02"
    m = json.loads((run_dir / "manifest.json").read_text())
    assert (run_dir / "RIG_r02.mfx").is_file() and (run_dir / "geometry_0001.stl").is_file()
    assert m["rig"] == "piv_rig" and m["flow_slpm"] == 26.5 and m["is_ofat_baseline"]
    assert m["stl_R_rim_m"] == pytest.approx(27.5e-3) and m["stl_n_bands_rim"] == 1
    assert m["case"].startswith("PIV rig") and m["nranks"] == 4
    assert (run_dir / "gridmap.dat").read_text().splitlines()[0] == "2 1 2"

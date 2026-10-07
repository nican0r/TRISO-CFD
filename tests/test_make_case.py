"""Tests for src.make_case render_case + builders."""

from __future__ import annotations

import math
import re

import pytest

from src.make_case import (
    REPO_ROOT,
    build_fb_sweep_params,
    build_V_2d_params,
    render_case,
)

TEMPLATE = REPO_ROOT / "cases" / "fb_sweep" / "template.mfx.j2"
V2D_TEMPLATE = REPO_ROOT / "cases" / "V_2d" / "template.mfx.j2"


# ---------- fb_sweep (Day 6-7) ----------

def _render(tmp_path, U, imax, jmax, tstop=8.0):
    params = build_fb_sweep_params(
        U_in_m_s=U, imax=imax, jmax=jmax, tstop_s=tstop, run_name="test_case"
    )
    out = tmp_path / "test_case.mfx"
    render_case(TEMPLATE, params, out)
    return out.read_text(), params


def test_render_sets_bc_v_g(tmp_path):
    text, _ = _render(tmp_path, U=0.5, imax=38, jmax=150)
    m = re.search(r"bc_v_g\(1\)\s*=\s*([0-9.eE+\-]+)", text)
    assert m is not None
    assert float(m.group(1)) == pytest.approx(0.5)


def test_render_sets_mesh(tmp_path):
    text, _ = _render(tmp_path, U=0.3, imax=40, jmax=160)
    assert re.search(r"imax\s*=\s*40\b", text)
    assert re.search(r"jmax\s*=\s*160\b", text)


def test_render_sets_case_D_particle(tmp_path):
    text, params = _render(tmp_path, U=0.1, imax=38, jmax=150)
    # Case D (post-retargeting): d_p = 500 um, rho_p = 6050 kg/m3 (S2).
    assert params["d_p_m"] == pytest.approx(5.0e-4)
    assert params["rho_p_kg_m3"] == pytest.approx(6050.0)
    assert re.search(r"d_p0\(1\)\s*=\s*0\.0005\b|d_p0\(1\)\s*=\s*5\.0*e-0?4\b",
                     text, flags=re.IGNORECASE)
    assert re.search(r"ro_s0\(1\)\s*=\s*6050(\.0)?\b", text)


def test_render_sets_air_gas_props(tmp_path):
    text, params = _render(tmp_path, U=0.1, imax=38, jmax=150)
    # Cantera gri30 at 293.15 K, 101325 Pa, N2/O2 = 0.79/0.21 gives
    # rho ~ 1.199 kg/m3, mu ~ 1.83e-5 Pa*s.
    assert 1.15 < params["rho_g_kg_m3"] < 1.25
    assert 1.7e-5 < params["mu_g_pa_s"] < 1.9e-5
    assert "ro_g0" in text and "mu_g0" in text


def test_render_bed_ic_matches_H_static(tmp_path):
    text, params = _render(tmp_path, U=0.1, imax=38, jmax=150)
    H = params["H_static_m"]
    m = re.search(r"ic_y_n\(2\)\s*=\s*([0-9.eE+\-]+)", text)
    assert m is not None and float(m.group(1)) == pytest.approx(H)
    m2 = re.search(r"ic_ep_s\(2,1\)\s*=\s*([0-9.eE+\-]+)", text)
    assert m2 is not None and float(m2.group(1)) == pytest.approx(params["eps_s_bed"])


def test_strict_undefined_catches_typo(tmp_path):
    params = build_fb_sweep_params(U_in_m_s=0.3, imax=38, jmax=150,
                                   tstop_s=8.0, run_name="t")
    del params["tstop_s"]
    with pytest.raises(Exception):
        render_case(TEMPLATE, params, tmp_path / "x.mfx")


def test_fb_sweep_ep_star_matches_case_V(tmp_path):
    """fb_sweep must render the same ep_star as Case V (same particles)."""
    text, params = _render(tmp_path, U=0.3, imax=38, jmax=150)
    assert params["ep_star"] == pytest.approx(0.42)
    assert (1.0 - params["ep_star"]) > params["eps_s_bed"]
    assert re.search(r"\bep_star\s*=\s*0\.42\b", text)


# ---------- Day 8: Case V UTK/ORNL 2D template ----------

def _render_V(tmp_path, U=30.0, tstop=3.0, name="V_2d"):
    params = build_V_2d_params(U_in_m_s=U, tstop_s=tstop, run_name=name)
    out = tmp_path / f"{name}.mfx"
    render_case(V2D_TEMPLATE, params, out)
    return out.read_text(), params


def test_V_2d_coordinates_and_no_k(tmp_path):
    text, _ = _render_V(tmp_path)
    assert re.search(r"coordinates\s*=\s*'CARTESIAN'", text)
    assert re.search(r"cartesian_grid\s*=\s*\.True\.", text)
    assert re.search(r"use_stl\s*=\s*\.False\.", text)
    assert re.search(r"^\s*no_k\s*=\s*\.True\.", text, re.MULTILINE)
    assert re.search(r"n_quadric\s*=\s*1\b", text)
    assert re.search(r"quadric_form\(1\)\s*=\s*'Y_CONE'", text)
    assert re.search(r"half_angle\(1\)\s*=\s*30\.0\b", text)


def test_V_2d_grid_and_orifice_cells(tmp_path):
    text, params = _render_V(tmp_path)
    assert re.search(r"imax\s*=\s*50\b", text)
    assert re.search(r"jmax\s*=\s*200\b", text)
    assert re.search(r"kmax\s*=\s*1\b", text)
    # Orifice must span >= 4 cells (day-8 acceptance).
    dx = (params["x_max_m"] - params["x_min_m"]) / params["imax"]
    orifice_cells = (2.0 * params["R_i_m"]) / dx
    assert orifice_cells >= 4.0
    # Cell size <= 10 * d_p = 5 mm in the particle zone.
    assert dx <= 10.0 * params["d_p_m"]


def test_V_2d_bc_mi_and_po_present(tmp_path):
    text, params = _render_V(tmp_path)
    assert re.search(r"bc_type\(1\)\s*=\s*'MI'", text)
    m = re.search(r"bc_v_g\(1\)\s*=\s*([-+eE.0-9]+)", text)
    assert m is not None and float(m.group(1)) == pytest.approx(params["U_in_m_s"])
    m_xw = re.search(r"bc_x_w\(1\)\s*=\s*([-+eE.0-9]+)", text)
    m_xe = re.search(r"bc_x_e\(1\)\s*=\s*([-+eE.0-9]+)", text)
    assert float(m_xw.group(1)) == pytest.approx(-params["R_i_m"])
    assert float(m_xe.group(1)) == pytest.approx(params["R_i_m"])
    assert re.search(r"bc_type\(2\)\s*=\s*'PO'", text)
    m2 = re.search(r"bc_y_s\(2\)\s*=\s*([-+eE.0-9]+)", text)
    assert m2 is not None and float(m2.group(1)) == pytest.approx(params["H_dom_m"])
    assert re.search(r"bc_type\(3\)\s*=\s*'CG_NSW'", text)
    assert re.search(r"bc_id_q\(1\)\s*=\s*3\b", text)


def test_V_2d_monitors_present(tmp_path):
    text, _ = _render_V(tmp_path)
    assert "V_2d_P_inlet" in text
    assert "V_2d_P_top_bed" in text
    assert "V_2d_solids_inventory" in text
    assert "V_2d_centerline_vg" in text
    # 1 kHz sampling on the pressure monitors (matches S1 and S2's 1000 Hz).
    assert re.search(r"monitor_dt\(1\)\s*=\s*0\.001\b", text)
    assert re.search(r"monitor_dt\(2\)\s*=\s*0\.001\b", text)
    # Inventory monitor is a volume integral (monitor_type = 11).
    assert re.search(r"monitor_type\(3\)\s*=\s*11\b", text)
    assert re.search(r"monitor_ep_s\(3,1\)\s*=\s*\.True\.", text)
    # Centerline v_g probe (monitor_type = 4 arith avg over thin strip at H_mid).
    assert re.search(r"monitor_type\(4\)\s*=\s*4\b", text)
    assert re.search(r"monitor_v_g\(4\)\s*=\s*\.True\.", text)


def test_V_2d_bed_ic_from_inventory(tmp_path):
    """Hand calc: inventory 54.5 g, rho_p = 6050 kg/m^3, eps_s_bed = 0.57
    in a 60 deg cone with 4 mm orifice gives H_static ~ 0.0322 m, which
    puts H0/D_c ~ 0.644 (inside S2's reported 0.50-0.65 range).
    """
    text, params = _render_V(tmp_path)
    # eps_s_bed set to 0.57 (below 1 - ep_star = 0.58 so Schaeffer inactive).
    assert params["eps_s_bed"] == pytest.approx(0.57)
    # Volume balance:
    mass_kg = 0.0545
    rho_p = 6050.0
    V_bed = mass_kg / (rho_p * params["eps_s_bed"])
    R_i = params["R_i_m"]
    tan_t = math.tan(math.radians(params["cone_half_angle_deg"]))
    h = params["H_static_m"]
    r_h = R_i + h * tan_t
    V_check = (math.pi * h / 3.0) * (R_i * R_i + R_i * r_h + r_h * r_h)
    assert V_check == pytest.approx(V_bed, rel=1e-3)
    # H0/D_c inside S2's reported 0.50-0.65 range.
    H0_over_Dc = h / (2.0 * params["R_c_m"])
    assert 0.50 <= H0_over_Dc <= 0.65
    # Rendered IC2 matches.
    m = re.search(r"ic_y_n\(2\)\s*=\s*([-+eE.0-9]+)", text)
    assert float(m.group(1)) == pytest.approx(params["H_static_m"])
    m2 = re.search(r"ic_ep_s\(2,1\)\s*=\s*([-+eE.0-9]+)", text)
    assert float(m2.group(1)) == pytest.approx(0.57)


def test_V_2d_inlet_volumetric_flow(tmp_path):
    """Hand calc: 30 m/s through a 4 mm circular inlet gives Q ~ 22.6 L/min.
    S2 reports ~12 L/min minimum spouting for 500 um ZrO2, so U/U_ms ~ 1.9.
    """
    _, params = _render_V(tmp_path, U=30.0)
    D_i = 2.0 * params["R_i_m"]
    A_inlet = math.pi * (D_i / 2.0) ** 2
    Q_m3_s = params["U_in_m_s"] * A_inlet
    Q_L_min = Q_m3_s * 1000.0 * 60.0
    assert Q_L_min == pytest.approx(22.6, abs=0.5)
    # U/U_ms vs S2's 12 L/min.
    U_over_Ums = Q_L_min / 12.0
    assert 1.5 <= U_over_Ums <= 2.3


def test_V_2d_closure_keywords(tmp_path):
    """S1 Table 1 closures: SYAM_OBRIEN drag, SCHAEFFER friction with
    SIGM_BLEND blending ("modified sigmoidal"), LUN_1984 PDE granular
    energy, phi = 15 deg, ep_star = 0.42.
    """
    text, params = _render_V(tmp_path)
    assert re.search(r"drag_type\s*=\s*'SYAM_OBRIEN'", text)
    assert re.search(r"friction_model\s*=\s*'SCHAEFFER'", text)
    # 'GIDASPOW_PCF' is MFiX 26.1.2's accepted token for scaled-sigmoidal
    # blending (= S1 Table 1's "modified sigmoidal"); the namelist header
    # lists 'SIGM_BLEND' but check_blending_function only accepts
    # 'GIDASPOW_PCF' to activate the SIGM_BLEND flag.
    assert re.search(r"blending_function\s*=\s*'GIDASPOW_PCF'", text)
    assert re.search(r"kt_type\s*=\s*'LUN_1984'", text)
    assert re.search(r"^\s*phi\s*=\s*15(\.0)?\b", text, re.MULTILINE)
    assert re.search(r"^\s*ep_star\s*=\s*0\.42\b", text, re.MULTILINE)
    assert re.search(r"^\s*c_e\s*=\s*0\.9\b", text, re.MULTILINE)
    assert params["phi_deg"] == pytest.approx(15.0)
    assert params["ep_star"] == pytest.approx(0.42)
    assert params["c_e"] == pytest.approx(0.9)


def test_V_2d_particles_match_case_V(tmp_path):
    """Case V particles: 500 um ZrO2 at 6050 kg/m^3 (S2).  Air at ~298 K."""
    _, params = _render_V(tmp_path)
    assert params["d_p_m"] == pytest.approx(5.0e-4)
    assert params["rho_p_kg_m3"] == pytest.approx(6050.0)
    assert 1.15 < params["rho_g_kg_m3"] < 1.25
    assert 1.7e-5 < params["mu_g_pa_s"] < 1.9e-5


def test_V_2d_ep_star_above_eps_s_bed(tmp_path):
    """1 - ep_star must sit above eps_s_bed so Schaeffer is inactive at t=0."""
    _, params = _render_V(tmp_path)
    assert (1.0 - params["ep_star"]) > params["eps_s_bed"]

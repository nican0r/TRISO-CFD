"""Tests for src.make_case render_case + build_fb_sweep_params."""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from src.make_case import (
    REPO_ROOT,
    build_fb_sweep_params,
    build_V_2d_day9_params,
    build_V_2d_params,
    build_V_3d_params,
    render_case,
)

TEMPLATE = REPO_ROOT / "cases" / "fb_sweep" / "template.mfx.j2"
V2D_TEMPLATE = REPO_ROOT / "cases" / "V_2d" / "template.mfx.j2"
V3D_TEMPLATE = REPO_ROOT / "cases" / "V_3d" / "template.mfx.j2"


def _render(tmp_path, U, imax, jmax, tstop=8.0):
    params = build_fb_sweep_params(
        U_in_m_s=U, imax=imax, jmax=jmax, tstop_s=tstop, run_name="test_case"
    )
    out = tmp_path / "test_case.mfx"
    render_case(TEMPLATE, params, out)
    return out.read_text(), params


def test_render_sets_bc_v_g(tmp_path):
    text, _ = _render(tmp_path, U=0.5, imax=38, jmax=150)
    # bc_v_g(1) is the inlet superficial velocity in the generated .mfx.
    m = re.search(r"bc_v_g\(1\)\s*=\s*([0-9.eE+\-]+)", text)
    assert m is not None
    assert float(m.group(1)) == pytest.approx(0.5)


def test_render_sets_mesh(tmp_path):
    text, _ = _render(tmp_path, U=0.3, imax=40, jmax=160)
    assert re.search(r"imax\s*=\s*40\b", text)
    assert re.search(r"jmax\s*=\s*160\b", text)


def test_render_sets_case_D_particle(tmp_path):
    text, params = _render(tmp_path, U=0.1, imax=38, jmax=150)
    # Case D: d_p = 500 um = 5.0e-4 m, rho_p = 6000 kg/m3 (params/particles.yaml).
    assert params["d_p_m"] == pytest.approx(5.0e-4)
    assert params["rho_p_kg_m3"] == pytest.approx(6000.0)
    assert re.search(r"d_p0\(1\)\s*=\s*0\.0005\b|d_p0\(1\)\s*=\s*5\.0*e-0?4\b",
                     text, flags=re.IGNORECASE)
    assert re.search(r"ro_s0\(1\)\s*=\s*6000(\.0)?\b", text)


def test_render_sets_air_gas_props(tmp_path):
    text, params = _render(tmp_path, U=0.1, imax=38, jmax=150)
    # Cantera gri30 at 293.15 K, 101325 Pa, N2/O2 = 0.79/0.21 gives
    # rho ~ 1.199 kg/m3, mu ~ 1.83e-5 Pa*s.
    assert 1.15 < params["rho_g_kg_m3"] < 1.25
    assert 1.7e-5 < params["mu_g_pa_s"] < 1.9e-5
    assert "ro_g0" in text and "mu_g0" in text


def test_render_bed_ic_matches_H_static(tmp_path):
    text, params = _render(tmp_path, U=0.1, imax=38, jmax=150)
    # IC2 is the packed bed; its north edge should be H_static.
    H = params["H_static_m"]
    m = re.search(r"ic_y_n\(2\)\s*=\s*([0-9.eE+\-]+)", text)
    assert m is not None and float(m.group(1)) == pytest.approx(H)
    m2 = re.search(r"ic_ep_s\(2,1\)\s*=\s*([0-9.eE+\-]+)", text)
    assert m2 is not None and float(m2.group(1)) == pytest.approx(params["eps_s_bed"])


def test_strict_undefined_catches_typo(tmp_path):
    # Dropping a key the template needs must raise, not render silently.
    params = build_fb_sweep_params(U_in_m_s=0.3, imax=38, jmax=150,
                                   tstop_s=8.0, run_name="t")
    del params["tstop_s"]
    with pytest.raises(Exception):
        render_case(TEMPLATE, params, tmp_path / "x.mfx")


# ---------- Day 8 : Case V templates ----------

def _render_V(tmp_path, template, builder, U=20.0, tstop=5.0, name="t"):
    params = builder(U_in_m_s=U, tstop_s=tstop, run_name=name)
    out = tmp_path / f"{name}.mfx"
    render_case(template, params, out)
    return out.read_text(), params


def test_V_2d_coordinates_and_no_k(tmp_path):
    text, _ = _render_V(tmp_path, V2D_TEMPLATE, build_V_2d_params, name="V_2d")
    # Day-8 rev 4: V_2d is Cartesian + cartesian_grid cut-cells + Y_CONE
    # quadric + no_k=.True. (true 2D).  Quadrics are compatible with no_k,
    # unlike STL (get_stl_data.f:813), so V_2d runs with kmax=1.
    assert re.search(r"coordinates\s*=\s*'CARTESIAN'", text)
    assert re.search(r"cartesian_grid\s*=\s*\.True\.", text)
    assert re.search(r"use_stl\s*=\s*\.False\.", text)
    assert re.search(r"^\s*no_k\s*=\s*\.True\.", text, re.MULTILINE)
    # Y_CONE quadric must be declared.
    assert re.search(r"n_quadric\s*=\s*1\b", text)
    assert re.search(r"quadric_form\(1\)\s*=\s*'Y_CONE'", text)
    assert re.search(r"half_angle\(1\)\s*=\s*30\.0\b", text)


def test_V_2d_grid_and_orifice_cells(tmp_path):
    text, params = _render_V(tmp_path, V2D_TEMPLATE, build_V_2d_params, name="V_2d")
    assert re.search(r"imax\s*=\s*64\b", text)
    assert re.search(r"jmax\s*=\s*336\b", text)
    assert re.search(r"kmax\s*=\s*1\b", text)
    # Orifice must span >= 4 cells (day-8 acceptance).  Full-width domain:
    # dx = 2*R_c / imax, orifice width 2*R_i.
    dx = (params["x_max_m"] - params["x_min_m"]) / params["imax"]
    orifice_cells = (2.0 * params["R_i_m"]) / dx
    assert orifice_cells >= 4.0


def test_V_2d_bc_mi_and_po_present(tmp_path):
    text, params = _render_V(tmp_path, V2D_TEMPLATE, build_V_2d_params, name="V_2d")
    # BC 1 : MI at y=0, -R_i <= x <= +R_i, bc_v_g = U_in
    assert re.search(r"bc_type\(1\)\s*=\s*'MI'", text)
    m = re.search(r"bc_v_g\(1\)\s*=\s*([-+eE.0-9]+)", text)
    assert m is not None and float(m.group(1)) == pytest.approx(params["U_in_m_s"])
    m_xw = re.search(r"bc_x_w\(1\)\s*=\s*([-+eE.0-9]+)", text)
    m_xe = re.search(r"bc_x_e\(1\)\s*=\s*([-+eE.0-9]+)", text)
    assert float(m_xw.group(1)) == pytest.approx(-params["R_i_m"])
    assert float(m_xe.group(1)) == pytest.approx(params["R_i_m"])
    # BC 2 : PO at the top (y = H_dom)
    assert re.search(r"bc_type\(2\)\s*=\s*'PO'", text)
    m2 = re.search(r"bc_y_s\(2\)\s*=\s*([-+eE.0-9]+)", text)
    assert m2 is not None and float(m2.group(1)) == pytest.approx(params["H_dom_m"])
    # BC 3 : cut-cell NSW tied to the Y_CONE quadric via bc_id_q(1)=3
    assert re.search(r"bc_type\(3\)\s*=\s*'CG_NSW'", text)
    assert re.search(r"bc_id_q\(1\)\s*=\s*3\b", text)


def test_V_2d_monitors_present(tmp_path):
    text, params = _render_V(tmp_path, V2D_TEMPLATE, build_V_2d_params, name="V_2d")
    assert "V_2d_P_inlet" in text
    assert "V_2d_P_top_bed" in text
    assert "V_2d_spout_mid" in text
    assert re.search(r"monitor_p_g\(1\)\s*=\s*\.True\.", text)
    assert re.search(r"monitor_p_g\(2\)\s*=\s*\.True\.", text)
    assert re.search(r"monitor_ep_s\(3,1\)\s*=\s*\.True\.", text)
    assert re.search(r"monitor_v_s\(3,1\)\s*=\s*\.True\.", text)
    # Spout mid-bed probe sits at H_static / 2 and spans |x| <= R_i
    m = re.search(r"monitor_y_s\(3\)\s*=\s*([-+eE.0-9]+)", text)
    assert float(m.group(1)) == pytest.approx(params["H_static_m"] / 2.0)
    m_xw = re.search(r"monitor_x_w\(3\)\s*=\s*([-+eE.0-9]+)", text)
    m_xe = re.search(r"monitor_x_e\(3\)\s*=\s*([-+eE.0-9]+)", text)
    assert float(m_xw.group(1)) == pytest.approx(-params["R_i_m"])
    assert float(m_xe.group(1)) == pytest.approx(params["R_i_m"])


def test_V_2d_bed_ic(tmp_path):
    text, params = _render_V(tmp_path, V2D_TEMPLATE, build_V_2d_params, name="V_2d")
    # Packed-bed IC fills x in [-R_c, +R_c] and y in [0, H_static] at
    # eps_s = 0.59 (= 1 - eps_mf, Missouri S&T RPT minimum fluidization).
    m = re.search(r"ic_y_n\(2\)\s*=\s*([-+eE.0-9]+)", text)
    assert float(m.group(1)) == pytest.approx(params["H_static_m"])
    m = re.search(r"ic_ep_s\(2,1\)\s*=\s*([-+eE.0-9]+)", text)
    assert float(m.group(1)) == pytest.approx(0.59)
    m_xw = re.search(r"ic_x_w\(2\)\s*=\s*([-+eE.0-9]+)", text)
    m_xe = re.search(r"ic_x_e\(2\)\s*=\s*([-+eE.0-9]+)", text)
    assert float(m_xw.group(1)) == pytest.approx(-params["R_c_m"])
    assert float(m_xe.group(1)) == pytest.approx(params["R_c_m"])


def test_V_3d_stl_and_grid(tmp_path):
    text, params = _render_V(tmp_path, V3D_TEMPLATE, build_V_3d_params, name="V_3d")
    assert re.search(r"cartesian_grid\s*=\s*\.True\.", text)
    assert re.search(r"use_stl\s*=\s*\.True\.", text)
    # Grid 31 x 160 x 31
    assert re.search(r"imax\s*=\s*31\b", text)
    assert re.search(r"jmax\s*=\s*160\b", text)
    assert re.search(r"kmax\s*=\s*31\b", text)
    # Cell size ~ 2.45 mm (within 2-3 mm per day-8.md)
    dx = (params["x_max_m"] - params["x_min_m"]) / params["imax"]
    dy = params["H_dom_m"] / params["jmax"]
    dz = (params["z_max_m"] - params["z_min_m"]) / params["kmax"]
    for d in (dx, dy, dz):
        assert 0.002 <= d <= 0.003


def test_V_3d_monitors_present(tmp_path):
    text, _ = _render_V(tmp_path, V3D_TEMPLATE, build_V_3d_params, name="V_3d")
    assert "V_3d_P_inlet" in text
    assert "V_3d_P_top_bed" in text
    assert "V_3d_spout_mid" in text
    assert re.search(r"monitor_ep_s\(3,1\)\s*=\s*\.True\.", text)
    assert re.search(r"monitor_v_s\(3,1\)\s*=\s*\.True\.", text)


def test_V_params_use_case_V_rpt_particle(tmp_path):
    """V_2d / V_3d must carry the Missouri S&T RPT validation inventory:
    2.18 mm soda-lime glass beads, rho_p = 2400 kg/m3, in dry air at 298 K.
    Source: R13 (OSTI 1533774) / R14 (OSTI 1538701) experimental conditions.
    """
    _, p2 = _render_V(tmp_path, V2D_TEMPLATE, build_V_2d_params, name="V_2d")
    _, p3 = _render_V(tmp_path, V3D_TEMPLATE, build_V_3d_params, name="V_3d")
    for p in (p2, p3):
        assert p["d_p_m"] == pytest.approx(2.18e-3)
        assert p["rho_p_kg_m3"] == pytest.approx(2400.0)
        # air @ 25 C / 101 kPa from Cantera: paper reports rho_f = 1.21 kg/m3,
        # mu = 1.81e-5 Pa.s.
        assert 1.15 < p["rho_g_kg_m3"] < 1.25
        assert 1.7e-5 < p["mu_g_pa_s"] < 1.9e-5
        # Internal friction angle 26 deg piped through to MFiX `phi`.
        assert p["phi_deg"] == pytest.approx(26.0)


def test_V_ep_star_above_eps_s_bed(tmp_path):
    """Schaeffer activation threshold 1-ep_star must sit above the IC bed value.

    ep_star is a physical constant (RCP of monodisperse smooth spheres,
    eps_s,RCP ~ 0.636 => ep_star = 0.37). eps_s_bed = 0.59 = 1 - eps_mf
    where eps_mf = 0.41 is the Missouri S&T RPT measured minimum-fluidization
    voidage for the 2.18 mm glass beads. If the inequality flips, the
    frictional-stress term fires at t=0 and the solver collapses dt immediately
    (Day 6-7 blow-up pattern).
    """
    for template, builder in (
        (V2D_TEMPLATE, build_V_2d_params),
        (V3D_TEMPLATE, build_V_3d_params),
    ):
        text, params = _render_V(tmp_path, template, builder, name="t")
        assert params["ep_star"] == pytest.approx(0.37)
        assert params["eps_s_bed"] == pytest.approx(0.59)
        assert (1.0 - params["ep_star"]) > params["eps_s_bed"]
        assert re.search(r"\bep_star\s*=\s*0\.37\b", text)
        # phi = 26 deg (Case V_rpt inertial friction angle).
        assert re.search(r"\bphi\s*=\s*26(\.0)?\b", text)


# ---------- Day 9 : V_2d spouting run (1.2 x U_ms, SYAM_OBRIEN, e = 0.9) ----------

def test_V_2d_day9_U_in_matches_mathur_gishler(tmp_path):
    """U_in at the orifice = 1.2 x U_ms,MG x (D_c/D_i)  [2D slab area ratio].

    Hand calc with Case V_rpt (d_p = 2.18 mm, rho_p = 2400 kg/m^3), air @
    25 C (rho_g ~ 1.18 kg/m^3), 0.076 m / 0.0095 m / 0.10 m bed:
        U_ms,col = (2.18e-3/0.076) * (0.0095/0.076)^(1/3) * sqrt(2*9.80665*0.10
                                                            *(2400-1.18)/1.18)
                 ~ 0.906 m/s
        U_in     = 1.2 * 0.906 * (0.076/0.0095) = 1.2 * 0.906 * 8 ~ 8.70 m/s
    V_2d is a planar slab, so the orifice/column area ratio is D_c/D_i
    (strip widths), not (D_c/D_i)^2 (which is the 3D disc ratio).  See
    notes/day09_debug.md for the switch rationale.
    """
    params = build_V_2d_day9_params(run_name="V_2d_day9_test")
    # Column-superficial U_ms,MG hand value ~ 0.906 m/s (loose tol allows for
    # Cantera gri30 vs. ideal-gas rho_g differences at 298.15 K).
    assert params["U_ms_column_m_s"] == pytest.approx(0.906, rel=1e-2)
    # Orifice-face U_in ~ 8.70 m/s (1.2 x U_ms x D_c/D_i = x8).
    assert params["U_in_m_s"] == pytest.approx(8.70, rel=1e-2)
    assert params["u_over_ums"] == pytest.approx(1.2)


def test_V_2d_day9_renders_drag_and_restitution(tmp_path):
    params = build_V_2d_day9_params(run_name="V_2d_day9_test")
    out = tmp_path / "V_2d_day9_test.mfx"
    render_case(V2D_TEMPLATE, params, out)
    text = out.read_text()
    # drag_type keyword (Day-9 records GIDASPOW after Run 4 SYAM_OBRIEN crash;
    # see notes/day09_debug.md).
    assert re.search(r"drag_type\s*=\s*'GIDASPOW'", text)
    # Restitution coefficient dropped from 0.95 (Day 8) to 0.9 (Day 9).
    assert re.search(r"^\s*c_e\s*=\s*0\.9\b", text, re.MULTILINE)
    # Inlet velocity in the rendered .mfx matches the computed U_in.
    m = re.search(r"bc_v_g\(1\)\s*=\s*([-+eE.0-9]+)", text)
    assert m is not None
    assert float(m.group(1)) == pytest.approx(params["U_in_m_s"], rel=1e-6)
    # Tstop is the day-9 value from params/particles.yaml V_2d_day9.tstop_s
    # (dropped from the step-text 5 s to 2.5 s on the coarsened Day-9 mesh;
    # notes/day09_debug.md Run 2).
    m2 = re.search(r"tstop\s*=\s*([-+eE.0-9]+)", text)
    assert m2 is not None and float(m2.group(1)) == pytest.approx(2.5)
    # Day-9 mesh override: 32 x 168 (orifice spans exactly 4 cells, Day-8 floor).
    assert re.search(r"imax\s*=\s*32\b", text)
    assert re.search(r"jmax\s*=\s*168\b", text)
    dx = (params["x_max_m"] - params["x_min_m"]) / params["imax"]
    orifice_cells = (2.0 * params["R_i_m"]) / dx
    assert orifice_cells == pytest.approx(4.0, abs=1e-9)


def test_fb_sweep_ep_star_matches_case_V(tmp_path):
    """fb_sweep must render the same ep_star as the Case V templates.

    ep_star is a particle property (packing limit of 500 um monodisperse
    spheres); the two cases differ in geometry/gas but share Case D particles,
    so they must share the same ep_star value and threshold.
    """
    text, params = _render(tmp_path, U=0.3, imax=38, jmax=150)
    assert params["ep_star"] == pytest.approx(0.37)
    assert (1.0 - params["ep_star"]) > params["eps_s_bed"]
    assert re.search(r"\bep_star\s*=\s*0\.37\b", text)

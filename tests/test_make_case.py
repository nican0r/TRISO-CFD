"""Tests for src.make_case render_case + builders."""

from __future__ import annotations

import math
import re

import pytest

from src.make_case import (
    REPO_ROOT,
    build_fb_sweep_params,
    build_V_2d_params,
    build_V_3d_params,
    render_case,
)

TEMPLATE = REPO_ROOT / "cases" / "fb_sweep" / "template.mfx.j2"
V2D_TEMPLATE = REPO_ROOT / "cases" / "V_2d" / "template.mfx.j2"
V3D_TEMPLATE = REPO_ROOT / "cases" / "V_3d" / "template.mfx.j2"


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


# ---------- Day 10: Case V UTK/ORNL 3D STL cut-cell template ----------

def _render_V3(tmp_path, U=30.0, tstop=5.5, name="V_3d"):
    params = build_V_3d_params(U_in_m_s=U, tstop_s=tstop, run_name=name)
    out = tmp_path / f"{name}.mfx"
    render_case(V3D_TEMPLATE, params, out)
    return out.read_text(), params


def test_V_3d_cartesian_grid_and_stl(tmp_path):
    """Day-10 template uses STL cut-cell (surface of revolution via
    src/make_stl.py, written to geometry_0001.stl).  The STL is
    side-walls only (no top cap, no bottom floor) because flush closures
    are discarded by MFiX and offsetting the box breaks the MI/PO BC
    flow-plane requirement -- see notes/day10_summary.md for the full
    trade-space.
    """
    text, params = _render_V3(tmp_path)
    assert re.search(r"coordinates\s*=\s*'CARTESIAN'", text)
    assert re.search(r"cartesian_grid\s*=\s*\.True\.", text)
    assert re.search(r"use_stl\s*=\s*\.True\.", text)
    assert re.search(r"out_stl_value\s*=\s*1\.0\b", text)
    # no_k must NOT be .True. (3D case).
    assert not re.search(r"^\s*no_k\s*=\s*\.True\.", text, re.MULTILINE)
    # kmax > 1 (true 3D).
    assert params["kmax"] > 1


def test_V_3d_uniform_mesh_sum(tmp_path):
    """Day-10 R1: UNIFORM 2 mm Cartesian mesh.  x,z box spans ±(R_c + 5 mm
    margin) = ±30 mm so sum_x = sum_z = 60 mm > D_c = 50 mm.  y box spans
    [-L_stub, H_dom] so sum_y = L_stub + H_dom = 220 mm.
    """
    _, params = _render_V3(tmp_path)
    # Symmetric x,z box, encloses vessel.
    sum_x = sum(params["dx_m_list"])
    sum_z = sum(params["dz_m_list"])
    assert sum_x == pytest.approx(sum_z, abs=1e-9)
    assert sum_x >= 2.0 * params["R_c_m"]            # encloses D_c
    assert sum_x == pytest.approx(0.060, abs=1e-9)   # R1 numbers: 60 mm
    # y box = L_stub + H_dom.
    assert sum(params["dy_m_list"]) == pytest.approx(
        params["H_dom_m"] + params["L_stub_m"], abs=1e-9
    )
    assert params["imax"] == len(params["dx_m_list"])
    assert params["jmax"] == len(params["dy_m_list"])
    assert params["kmax"] == len(params["dz_m_list"])
    # Uniform mesh: every cell is dx_uniform_m (R1 = 2 mm).
    dx_u = params["dx_uniform_m"]
    assert dx_u == pytest.approx(2.0e-3)
    assert all(d == pytest.approx(dx_u) for d in params["dx_m_list"])
    assert all(d == pytest.approx(dx_u) for d in params["dy_m_list"])
    assert all(d == pytest.approx(dx_u) for d in params["dz_m_list"])
    # Day-10 R1 cell counts.
    assert params["imax"] == 30
    assert params["jmax"] == 110
    assert params["kmax"] == 30


def test_V_3d_orifice_and_particle_zone_cells(tmp_path):
    """Day-10 R1 (NETL-clone trade-off):
      - cell size <= 10 * d_p = 5 mm in the particle zone (met at 2 mm);
      - "4 cells across the 4 mm inlet" is impractical on a uniform mesh at
        feasible cell counts.  NETL's shipped pic/spouted_bed_3d tutorial
        uses ~4 cells across a 15 mm spout (a 2x cell-to-radius ratio);
        R1 accepts 2 cells across the 4 mm orifice on a uniform 2 mm mesh
        (same cell-to-radius ratio).  This mirrors the recipe that avoids
        the F_AT cut-cell classifier bug (see notes/day10_summary.md R1).
    """
    _, params = _render_V3(tmp_path)
    d_p = params["d_p_m"]
    # All cell sizes <= 10 * d_p = 5 mm in x, y, z.
    assert max(params["dy_m_list"]) <= 10.0 * d_p + 1e-12
    assert max(params["dx_m_list"]) <= 10.0 * d_p + 1e-12
    assert max(params["dz_m_list"]) <= 10.0 * d_p + 1e-12
    # Count orifice cells across x: cells with centre in [-R_i, R_i].  On
    # the uniform 2 mm mesh with R_i = 2 mm, 2 cells span the diameter.
    R_i = params["R_i_m"]
    xs = []
    x = params["x_min_m"]
    for d in params["dx_m_list"]:
        xs.append(x + 0.5 * d)
        x += d
    n_orifice_x = sum(1 for xc in xs if abs(xc) <= R_i)
    assert n_orifice_x >= 2        # R1 NETL-pattern lower bound
    # Same across z.
    zs = []
    z = params["z_min_m"]
    for d in params["dz_m_list"]:
        zs.append(z + 0.5 * d)
        z += d
    n_orifice_z = sum(1 for zc in zs if abs(zc) <= R_i)
    assert n_orifice_z >= 2
    # R1 near-axis cell size is exactly dx_uniform_m = 2 mm.
    near_axis_dx = min(params["dx_m_list"])
    near_axis_dz = min(params["dz_m_list"])
    assert near_axis_dx == pytest.approx(2.0e-3)
    assert near_axis_dz == pytest.approx(2.0e-3)


def test_V_3d_bc_mi_and_po_and_cgnsw(tmp_path):
    """Day-10 R1 (NETL-clone pattern): MI spans the FULL box bottom face;
    the STL stub cylinder (radius R_i) carves out the circular orifice-
    equivalent inflow footprint.  Compare NETL tutorials/pic/spouted_bed_3d
    which uses the same full-face MI pattern.
    """
    text, params = _render_V3(tmp_path)
    # BC_1 CG_NSW (STL walls, no region bounds -> applied to all cut-cells).
    assert re.search(r"bc_type\(1\)\s*=\s*'CG_NSW'", text)
    # BC_2 PO at y = H_dom, full x,z extent.
    assert re.search(r"bc_type\(2\)\s*=\s*'PO'", text)
    m = re.search(r"bc_y_s\(2\)\s*=\s*([-+eE.0-9]+)", text)
    assert m is not None and float(m.group(1)) == pytest.approx(params["H_dom_m"])
    # BC_3 MI at the box bottom (y = -L_stub), FULL BOX BOTTOM FACE.
    assert re.search(r"bc_type\(3\)\s*=\s*'MI'", text)
    m2 = re.search(r"bc_v_g\(3\)\s*=\s*([-+eE.0-9]+)", text)
    assert m2 is not None and float(m2.group(1)) == pytest.approx(params["U_in_m_s"])
    m_xw = re.search(r"bc_x_w\(3\)\s*=\s*([-+eE.0-9]+)", text)
    m_xe = re.search(r"bc_x_e\(3\)\s*=\s*([-+eE.0-9]+)", text)
    m_zb = re.search(r"bc_z_b\(3\)\s*=\s*([-+eE.0-9]+)", text)
    m_zt = re.search(r"bc_z_t\(3\)\s*=\s*([-+eE.0-9]+)", text)
    # Full box face bounds.
    assert float(m_xw.group(1)) == pytest.approx(params["x_min_m"])
    assert float(m_xe.group(1)) == pytest.approx(params["x_max_m"])
    assert float(m_zb.group(1)) == pytest.approx(params["z_min_m"])
    assert float(m_zt.group(1)) == pytest.approx(params["z_max_m"])
    m_ys = re.search(r"bc_y_s\(3\)\s*=\s*([-+eE.0-9]+)", text)
    m_yn = re.search(r"bc_y_n\(3\)\s*=\s*([-+eE.0-9]+)", text)
    assert float(m_ys.group(1)) == pytest.approx(-params["L_stub_m"])
    assert float(m_yn.group(1)) == pytest.approx(-params["L_stub_m"])
    # PIC-specific: bc_pic_mi_const_statwt must be set on the MI BC even
    # though bc_ep_s(3,1)=0.0 (parser requirement; see NETL tutorial).
    m_sw = re.search(r"bc_pic_mi_const_statwt\(3,1\)\s*=\s*([-+eE.0-9]+)", text)
    assert m_sw is not None
    assert float(m_sw.group(1)) == pytest.approx(params["bc_pic_mi_const_statwt"])
    m_eps = re.search(r"bc_ep_s\(3,1\)\s*=\s*([-+eE.0-9]+)", text)
    assert m_eps is not None and float(m_eps.group(1)) == pytest.approx(0.0)


def test_V_3d_monitors_and_vtk(tmp_path):
    text, _ = _render_V3(tmp_path)
    assert "V_3d_P_inlet" in text
    assert "V_3d_P_top_bed" in text
    assert "V_3d_solids_inventory" in text
    assert "V_3d_centerline_vg" in text
    assert re.search(r"monitor_dt\(1\)\s*=\s*0\.001\b", text)
    assert re.search(r"monitor_dt\(2\)\s*=\s*0\.001\b", text)
    assert re.search(r"monitor_type\(3\)\s*=\s*11\b", text)
    # PIC: monitor_ep_s is valid because des_interp_mean_fields projects
    # parcel mass to the Eulerian ep_s mean field.
    assert re.search(r"monitor_ep_s\(3,1\)\s*=\s*\.True\.", text)
    # VTK 1: Eulerian cell-centered background (gas + mean fields).
    assert re.search(r"vtk_filebase\(1\)\s*=\s*'Background'", text)
    assert re.search(r"vtk_vel_g\(1\)\s*=\s*\.True\.", text)
    # VTK 2: PIC parcel point-cloud dump (data='P') -- used by Day 11-12
    # CFD-DEM comparisons.
    assert re.search(r"vtk_filebase\(2\)\s*=\s*'Particles'", text)
    assert re.search(r"vtk_data\(2\)\s*=\s*'P'", text)
    assert re.search(r"vtk_part_diameter\(2\)\s*=\s*\.True\.", text)
    assert re.search(r"vtk_part_vel\(2\)\s*=\s*\.True\.", text)


def test_V_3d_closure_keywords(tmp_path):
    """Day 10 PIC closures.  SYAM_OBRIEN drag is kept from S1 Table 1
    (gas-side closure is physics-agnostic to the solids model).  TFM-only
    keywords (blending_function, kt_type, phi, c_e) must NOT appear.  PIC
    coefficients come from params/geometry.yaml spouted_bed_V_3d_pic
    (NETL tutorials/pic/spouted_bed_3d defaults).
    """
    text, params = _render_V3(tmp_path)
    # Gas-side closure (shared with the TFM template).
    assert re.search(r"drag_type\s*=\s*'SYAM_OBRIEN'", text)
    # Solids model is PIC, not TFM.
    assert re.search(r"solids_model\(1\)\s*=\s*'PIC'", text)
    # SCHAEFFER label is applied to the PIC solids-stress branch.
    assert re.search(r"friction_model\s*=\s*'SCHAEFFER'", text)
    # TFM-only keywords must be absent from the rendered deck.
    assert not re.search(r"blending_function", text)
    assert not re.search(r"\bkt_type\b", text)
    assert not re.search(r"^\s*phi\s*=", text, re.MULTILINE)
    assert not re.search(r"^\s*c_e\s*=", text, re.MULTILINE)
    # PIC closure keywords (verified against installed MFiX 26.1.2
    # reference share/mfix/doc/html/reference/particle_in_cell.html).
    assert re.search(r"fric_exp_pic\s*=\s*3(\.0)?\b", text)
    assert re.search(r"fric_non_sing_fac\s*=\s*1(\.0+)?e-0?7\b", text,
                     re.IGNORECASE)
    assert re.search(r"psfac_fric_pic\s*=\s*100(\.0)?\b", text)
    assert re.search(r"mppic_coeff_en1\s*=\s*0\.85\b", text)
    assert re.search(r"mppic_coeff_en_wall\s*=\s*0\.85\b", text)
    assert re.search(r"mppic_coeff_et_wall\s*=\s*1(\.0)?\b", text)
    assert re.search(r"mppic_velfac_coeff\s*=\s*1(\.0)?\b", text)
    # Shared DES / PIC interpolation keys.
    assert re.search(r"des_interp_mean_fields\s*=\s*\.True\.", text)
    assert re.search(r"des_interp_on\s*=\s*\.True\.", text)
    assert re.search(r"des_interp_scheme\s*=\s*'LINEAR_HAT'", text)
    assert re.search(r"gener_part_config\s*=\s*\.True\.", text)
    # des_epg_clip = 0.42 (NETL tutorial default; equals Case V ep_star,
    # a physics coincidence rather than a derived relationship -- both
    # describe the minimum-fluidization gas volume fraction).
    assert re.search(r"des_epg_clip\s*=\s*0\.42\b", text)
    # Params dict should expose PIC-specific numeric values.
    assert params["fric_exp_pic"] == pytest.approx(3.0)
    assert params["fric_non_sing_fac"] == pytest.approx(1.0e-7)
    assert params["psfac_fric_pic"] == pytest.approx(100.0)
    assert params["mppic_coeff_en1"] == pytest.approx(0.85)
    assert params["mppic_coeff_en_wall"] == pytest.approx(0.85)
    assert params["mppic_coeff_et_wall"] == pytest.approx(1.0)
    # ep_star kept for the IC bed-clearance comparison.  For Case V,
    # ep_star = 0.42 happens to equal the NETL tutorial's des_epg_clip
    # default; this is a physics coincidence (both are min-fluidization
    # gas volume fraction), not a derived identity.
    assert params["ep_star"] == pytest.approx(0.42)
    assert params["des_epg_clip"] == pytest.approx(0.42)


def test_V_3d_bed_ic_matches_case_V(tmp_path):
    """IC2 bed in 3D uses the same H_static as 2D (derived from inventory).
    In 3D, cut-cells intersect the IC box with the STL cone frustum.
    IC4 (Day-10 fix 2 stub tube) sits at y_s = -L_stub, below IC2.
    For PIC each IC that overlaps the solids region must set
    ic_pic_const_statwt so parcels get seeded at startup.
    """
    text, params = _render_V3(tmp_path)
    assert params["eps_s_bed"] == pytest.approx(0.57)
    # Hand calc from inventory 54.5 g, rho_p 6050 kg/m^3, eps_s = 0.57 gives
    # frustum height ~ 32.2 mm (same as V_2d).
    assert params["H_static_m"] == pytest.approx(0.0322, abs=5e-4)
    m = re.search(r"ic_y_n\(2\)\s*=\s*([-+eE.0-9]+)", text)
    assert float(m.group(1)) == pytest.approx(params["H_static_m"])
    # IC2 still starts at y = 0 (bed sits above the stub).
    m_ys2 = re.search(r"ic_y_s\(2\)\s*=\s*([-+eE.0-9]+)", text)
    assert float(m_ys2.group(1)) == pytest.approx(0.0, abs=1e-12)
    # IC2 solids volume fraction renders.
    m_eps2 = re.search(r"ic_ep_s\(2,1\)\s*=\s*([-+eE.0-9]+)", text)
    assert float(m_eps2.group(1)) == pytest.approx(0.57)
    # PIC: every IC that touches the solids region must carry
    # ic_pic_const_statwt(IC, 1); for the bed IC (IC_2) check it matches
    # the YAML-provided value (NETL tutorial default 4.2).
    m_sw = re.search(r"ic_pic_const_statwt\(2,1\)\s*=\s*([-+eE.0-9]+)", text)
    assert m_sw is not None
    assert float(m_sw.group(1)) == pytest.approx(params["ic_pic_const_statwt"])
    # IC_3 covers the stub tube region below the orifice (full box x,z;
    # STL cut-cells block r > R_i so only the stub-interior cells get fluid).
    # The TFM-era pre-opened spout IC was removed for PIC (parcel-vs-gas
    # drag spike at t=0 crashed the first launch); NETL doesn't use one.
    # Stub tube renumbered IC_4 -> IC_3; IC_4 should not render.
    m_ys3 = re.search(r"ic_y_s\(3\)\s*=\s*([-+eE.0-9]+)", text)
    assert m_ys3 is not None
    assert float(m_ys3.group(1)) == pytest.approx(-params["L_stub_m"])
    m_yn3 = re.search(r"ic_y_n\(3\)\s*=\s*([-+eE.0-9]+)", text)
    assert float(m_yn3.group(1)) == pytest.approx(0.0, abs=1e-12)
    assert re.search(r"ic_x_w\(4\)", text) is None


def test_V_3d_dx_rendered_in_template(tmp_path):
    """DX, DY, DZ arrays must be rendered one line per index (dx(0)..dx(IMAX-1))."""
    text, params = _render_V3(tmp_path)
    for i, d in enumerate(params["dx_m_list"]):
        pattern = rf"\bdx\({i}\)\s*=\s*([-+eE.0-9]+)"
        m = re.search(pattern, text)
        assert m is not None, f"dx({i}) missing in template"
        assert float(m.group(1)) == pytest.approx(d, rel=1e-5)
    for j, d in enumerate(params["dy_m_list"]):
        pattern = rf"\bdy\({j}\)\s*=\s*([-+eE.0-9]+)"
        m = re.search(pattern, text)
        assert m is not None, f"dy({j}) missing in template"
        assert float(m.group(1)) == pytest.approx(d, rel=1e-5)


def test_V_3d_stl_stub_circular_cross_section(tmp_path):
    """Day-10 R1 (NETL-clone pattern): the MI spans the FULL box bottom
    face, so the earlier "square MI vs circular orifice" ratio check
    (ratio = 4/pi) no longer applies.  What matters now is that the STL
    stub cylinder, which carves the MI down to a circular inflow
    footprint, has cross-section pi * R_i^2 at y = y_min.  Check that
    orifice-aligned params expose this radius and the circular area is
    the expected 12.57 mm^2.
    """
    _, params = _render_V3(tmp_path)
    R_i = params["R_i_m"]
    # The stub cylinder cross-section at y = y_min is pi R_i^2 (R_i = 2 mm
    # => A_circle = 12.566 mm^2).  The STL cuts the box-bottom MI down to
    # this circle; cells outside r = R_i are BLOCKED by the STL and see
    # no inflow.
    A_circle = math.pi * R_i ** 2
    assert A_circle == pytest.approx(1.2566e-5, rel=5e-3)
    # Orifice-aligned bounds are now distinct from the MI bounds and sit
    # at |x|,|z| <= R_i (used by Monitors 1 and 4 — inlet-plane gas
    # pressure and centerline v_g).  The TFM-era pre-opened spout IC
    # that previously consumed these bounds has been removed for PIC
    # (NETL pattern: inlet jet carves the spout dynamically).
    assert params["orifice_x_w_m"] == pytest.approx(-R_i)
    assert params["orifice_x_e_m"] == pytest.approx(R_i)
    assert params["orifice_z_b_m"] == pytest.approx(-R_i)
    assert params["orifice_z_t_m"] == pytest.approx(R_i)
    # MI bounds now span the full box bottom face (NETL pattern).
    assert params["bc_mi_x_w_m"] == pytest.approx(params["x_min_m"])
    assert params["bc_mi_x_e_m"] == pytest.approx(params["x_max_m"])

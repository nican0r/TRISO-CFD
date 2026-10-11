"""Tests for src.make_case render_case + builders."""

from __future__ import annotations

import math
import re

import pytest

from src.make_case import (
    REPO_ROOT,
    build_fb_sweep_params,
    build_V_2d_params,
    build_V_3d_dem_params,
    build_V_3d_params,
    render_case,
)

TEMPLATE = REPO_ROOT / "cases" / "fb_sweep" / "template.mfx.j2"
V2D_TEMPLATE = REPO_ROOT / "cases" / "V_2d" / "template.mfx.j2"
V3D_TEMPLATE = REPO_ROOT / "cases" / "V_3d" / "template.mfx.j2"
V3D_DEM_TEMPLATE = REPO_ROOT / "cases" / "V_3d" / "template_dem.mfx.j2"


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
    """Day-10 R4: UNIFORM 1 mm Cartesian mesh (was 2 mm under Day-10 R1).  x,z
    box spans ±(R_c + 5 mm margin) = ±30 mm so sum_x = sum_z = 60 mm > D_c
    = 50 mm.  y box spans [-L_stub, H_dom] so sum_y = L_stub + H_dom =
    220 mm.
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
    # Uniform mesh: every cell is dx_uniform_m (Day-10 R4 = 1 mm).
    dx_u = params["dx_uniform_m"]
    assert dx_u == pytest.approx(1.0e-3)
    assert all(d == pytest.approx(dx_u) for d in params["dx_m_list"])
    assert all(d == pytest.approx(dx_u) for d in params["dy_m_list"])
    assert all(d == pytest.approx(dx_u) for d in params["dz_m_list"])
    # Day-10 R4 cell counts (box 60 x 220 x 60 mm at dx = 1 mm).
    assert params["imax"] == 60
    assert params["jmax"] == 220
    assert params["kmax"] == 60


def test_V_3d_orifice_and_particle_zone_cells(tmp_path):
    """Day-10 R4 (2026-10-11):
      - cell size <= 10 * d_p = 5 mm in the particle zone (met at 1 mm);
      - 4 cells span the 4 mm orifice diameter on the uniform 1 mm mesh,
        matching NETL's shipped pic/spouted_bed_3d tutorial resolution ratio
        (D_i / dx = 15 / 4 = 3.75 cells).  The Day-10 R1 2 mm mesh (= 2 cells
        across D_i) could not represent the 30 m/s jet in finite-volume
        cut-cells; smoke validation showed orifice v_y mean = 12.2 m/s
        against the 30 m/s target (notes/day10_summary.md R4).
    """
    _, params = _render_V3(tmp_path)
    d_p = params["d_p_m"]
    # All cell sizes <= 10 * d_p = 5 mm in x, y, z.
    assert max(params["dy_m_list"]) <= 10.0 * d_p + 1e-12
    assert max(params["dx_m_list"]) <= 10.0 * d_p + 1e-12
    assert max(params["dz_m_list"]) <= 10.0 * d_p + 1e-12
    # Count orifice cells across x: cells with centre in [-R_i, R_i].  On
    # the uniform 1 mm mesh with R_i = 2 mm, 4 cells span the diameter.
    R_i = params["R_i_m"]
    xs = []
    x = params["x_min_m"]
    for d in params["dx_m_list"]:
        xs.append(x + 0.5 * d)
        x += d
    n_orifice_x = sum(1 for xc in xs if abs(xc) <= R_i)
    assert n_orifice_x >= 4        # Day-10 R4 NETL-ratio lower bound
    # Same across z.
    zs = []
    z = params["z_min_m"]
    for d in params["dz_m_list"]:
        zs.append(z + 0.5 * d)
        z += d
    n_orifice_z = sum(1 for zc in zs if abs(zc) <= R_i)
    assert n_orifice_z >= 4
    # Day-10 R4 near-axis cell size is exactly dx_uniform_m = 1 mm.
    near_axis_dx = min(params["dx_m_list"])
    near_axis_dz = min(params["dz_m_list"])
    assert near_axis_dx == pytest.approx(1.0e-3)
    assert near_axis_dz == pytest.approx(1.0e-3)


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
    # BC_3 MI at the box bottom (y = -L_stub), FULL BOX BOTTOM FACE.  Day-10 R4
    # (straight-stub refinement): R_stub = R_i so the MI face velocity
    # bc_v_g_mi equals the paper's physical orifice velocity U_in_m_s
    # directly -- no mass-flux scaling.
    assert re.search(r"bc_type\(3\)\s*=\s*'MI'", text)
    m2 = re.search(r"bc_v_g\(3\)\s*=\s*([-+eE.0-9]+)", text)
    assert m2 is not None and float(m2.group(1)) == pytest.approx(
        params["bc_v_g_mi_m_s"]
    )
    # Mass flux check: bc_v_g * (R_stub)**2 == U_in_m_s * R_i**2.
    assert params["bc_v_g_mi_m_s"] * params["R_stub_m"] ** 2 == pytest.approx(
        params["U_in_m_s"] * params["R_i_m"] ** 2
    )
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
    """Day-10 R3 Fix A: packed bed is split into NETL's 4-quadrant strip
    pattern (IC_2 +z, IC_3 -z, IC_4 -x center, IC_5 +x center) that leaves
    a central gas-only column of half-width R_bed_gas_col untouched.  This
    replaces Day-10 R2's single full-box IC_2, which seeded parcels
    directly over the orifice and triggered the LINEAR_HAT projection
    pathology documented in notes/day10_pathology_review.md.  The
    inlet-stub IC_3 of R2 is dropped (redundant with IC_1).
    """
    text, params = _render_V3(tmp_path)
    assert params["eps_s_bed"] == pytest.approx(0.57)
    assert params["H_static_m"] == pytest.approx(0.0322, abs=5e-4)
    R_col = params["R_bed_gas_col_m"]
    # R_bed_gas_col must at least enclose the stub tube and provide a
    # one-cell margin around the orifice.
    assert R_col >= params["R_stub_m"] - 1e-12
    assert R_col >= params["R_i_m"] + params["dx_uniform_m"] - 1e-12
    # Each bed strip IC_2..IC_5 covers 0 <= y <= H_static with ep_s=0.57
    # and the PIC statwt.
    for ic in (2, 3, 4, 5):
        m_ys = re.search(rf"ic_y_s\({ic}\)\s*=\s*([-+eE.0-9]+)", text)
        m_yn = re.search(rf"ic_y_n\({ic}\)\s*=\s*([-+eE.0-9]+)", text)
        m_eps = re.search(rf"ic_ep_s\({ic},1\)\s*=\s*([-+eE.0-9]+)", text)
        m_sw = re.search(rf"ic_pic_const_statwt\({ic},1\)\s*=\s*([-+eE.0-9]+)", text)
        assert m_ys is not None, f"IC_{ic} y_s missing"
        assert float(m_ys.group(1)) == pytest.approx(0.0, abs=1e-12)
        assert float(m_yn.group(1)) == pytest.approx(params["H_static_m"])
        assert float(m_eps.group(1)) == pytest.approx(0.57)
        assert float(m_sw.group(1)) == pytest.approx(params["ic_pic_const_statwt"])
    # IC_2 is the +z strip (z_b = +R_col, full x extent).
    m_zb2 = re.search(r"ic_z_b\(2\)\s*=\s*([-+eE.0-9]+)", text)
    m_xw2 = re.search(r"ic_x_w\(2\)\s*=\s*([-+eE.0-9]+)", text)
    m_xe2 = re.search(r"ic_x_e\(2\)\s*=\s*([-+eE.0-9]+)", text)
    assert float(m_zb2.group(1)) == pytest.approx(R_col)
    assert float(m_xw2.group(1)) == pytest.approx(params["x_min_m"])
    assert float(m_xe2.group(1)) == pytest.approx(params["x_max_m"])
    # IC_3 is the -z strip (z_t = -R_col, full x extent).
    m_zt3 = re.search(r"ic_z_t\(3\)\s*=\s*([-+eE.0-9]+)", text)
    assert float(m_zt3.group(1)) == pytest.approx(-R_col)
    # IC_4 is the -x center strip (x_e = -R_col, |z| <= R_col).
    m_xe4 = re.search(r"ic_x_e\(4\)\s*=\s*([-+eE.0-9]+)", text)
    m_zb4 = re.search(r"ic_z_b\(4\)\s*=\s*([-+eE.0-9]+)", text)
    m_zt4 = re.search(r"ic_z_t\(4\)\s*=\s*([-+eE.0-9]+)", text)
    assert float(m_xe4.group(1)) == pytest.approx(-R_col)
    assert float(m_zb4.group(1)) == pytest.approx(-R_col)
    assert float(m_zt4.group(1)) == pytest.approx(R_col)
    # IC_5 is the +x center strip (x_w = +R_col, |z| <= R_col).
    m_xw5 = re.search(r"ic_x_w\(5\)\s*=\s*([-+eE.0-9]+)", text)
    assert float(m_xw5.group(1)) == pytest.approx(R_col)
    # IC_6 (and beyond) should not render: the inlet stub region is
    # covered by IC_1, not an explicit IC per NETL's pattern.
    assert re.search(r"ic_x_w\(6\)", text) is None
    # Central gas column (|x| <= R_col AND |z| <= R_col) must not be
    # covered by any bed-strip IC: assert no bed IC has both x and z
    # bounds inside the column.
    for ic in (2, 3, 4, 5):
        xw = float(re.search(rf"ic_x_w\({ic}\)\s*=\s*([-+eE.0-9]+)", text).group(1))
        xe = float(re.search(rf"ic_x_e\({ic}\)\s*=\s*([-+eE.0-9]+)", text).group(1))
        zb = float(re.search(rf"ic_z_b\({ic}\)\s*=\s*([-+eE.0-9]+)", text).group(1))
        zt = float(re.search(rf"ic_z_t\({ic}\)\s*=\s*([-+eE.0-9]+)", text).group(1))
        x_inside = (xw >= -R_col - 1e-9) and (xe <= R_col + 1e-9)
        z_inside = (zb >= -R_col - 1e-9) and (zt <= R_col + 1e-9)
        assert not (x_inside and z_inside), (
            f"IC_{ic} seeds parcels inside the central gas column"
        )


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


def test_V_3d_straight_stub_params(tmp_path):
    """Day-10 R4 (2026-10-11): straight-stub geometry.

    R_stub == R_i (no stub widening) and taper_length_m == 0 (no taper
    needed with equal radii).  bc_v_g_mi equals the paper's physical jet
    velocity U_in_m_s directly, without the R2 Option-4 (R_i/R_stub)^2
    scaling that depended on 4x acceleration across a 1-cell taper the
    SIMPLE pressure correction could not deliver.  Enabled by the dx = 1 mm
    mesh refinement, which resolves R_i = 2 mm with 2 cells across the
    radius -- enough for F_AT to seed the stub as fluid.  See
    notes/day10_summary.md R4.
    """
    _, params = _render_V3(tmp_path)
    assert params["taper_length_m"] == pytest.approx(0.0)
    assert params["R_stub_m"] == pytest.approx(params["R_i_m"])
    assert params["bc_v_g_mi_m_s"] == pytest.approx(params["U_in_m_s"])
    assert params["stub_to_orifice_area_ratio"] == pytest.approx(1.0)
    # Mass-flux identity still holds trivially (R_stub == R_i).
    assert params["bc_v_g_mi_m_s"] * params["R_stub_m"] ** 2 == pytest.approx(
        params["U_in_m_s"] * params["R_i_m"] ** 2
    )


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


# ---------- Day 11-12: Case V UTK/ORNL 3D CFD-DEM template ----------

def _render_V3_dem(tmp_path, U=30.0, tstop=5.5, name="V_3d_dem",
                   cgp_factor=10.0):
    params = build_V_3d_dem_params(
        U_in_m_s=U, tstop_s=tstop, run_name=name, cgp_factor=cgp_factor,
    )
    out = tmp_path / f"{name}.mfx"
    render_case(V3D_DEM_TEMPLATE, params, out)
    return out.read_text(), params


def test_V_3d_dem_shares_geometry_with_pic(tmp_path):
    """DEM variant must reuse the Day-10 R1 PIC mesh/geometry verbatim.

    Same uniform 1 mm Cartesian box (Day-10 R4), same per-cell dx/dy/dz lists, same
    H_static (inventory-derived), same L_stub, same orifice.
    """
    _, dem = _render_V3_dem(tmp_path)
    pic = build_V_3d_params()
    for key in ("imax", "jmax", "kmax",
                "x_min_m", "x_max_m", "y_min_m", "y_max_m",
                "z_min_m", "z_max_m",
                "R_c_m", "R_i_m", "L_stub_m",
                "H_dom_m", "H_static_m",
                "eps_s_bed", "ep_star",
                "dx_uniform_m"):
        assert dem[key] == pytest.approx(pic[key]), (
            f"DEM variant diverged from Day-10 PIC on {key}: "
            f"{dem[key]} vs {pic[key]}"
        )
    assert dem["dx_m_list"] == pic["dx_m_list"]
    assert dem["dy_m_list"] == pic["dy_m_list"]
    assert dem["dz_m_list"] == pic["dz_m_list"]


# ---------- Day-10 R4: DMP (MPI) decomposition parameterization ----------

def test_V_3d_nodes_default_serial(tmp_path):
    """Default nodes=(1,1,1) means serial MFiX; template renders those values."""
    text, params = _render_V3(tmp_path)
    assert params["nodes_i"] == 1
    assert params["nodes_j"] == 1
    assert params["nodes_k"] == 1
    assert params["nranks"] == 1
    assert re.search(r"nodesi\s*=\s*1\b", text)
    assert re.search(r"nodesj\s*=\s*1\b", text)
    assert re.search(r"nodesk\s*=\s*1\b", text)


def test_V_3d_nodes_axial_five(tmp_path):
    """Day-10 R4 local-smoke default: 5 ranks on the long y-axis.

    jmax = 220 at dx = 1 mm; 220 % 5 == 0 so (1, 5, 1) divides cleanly
    (44 cells per rank in y).  5 ranks on 10 cores leaves 5 cores for
    interactive use (satisfies the leave-cores-free rule).  6 ranks is
    not axial-decomposable at this mesh (220 % 6 != 0); the alternative
    (3, 2, 1) adds halo surface for marginal gain.
    """
    params = build_V_3d_params(nodes=(1, 5, 1))
    out = tmp_path / "V_3d.mfx"
    render_case(V3D_TEMPLATE, params, out)
    text = out.read_text()
    assert params["nodes_i"] == 1
    assert params["nodes_j"] == 5
    assert params["nodes_k"] == 1
    assert params["nranks"] == 5
    assert re.search(r"nodesj\s*=\s*5\b", text)


def test_V_3d_nodes_three_two_one(tmp_path):
    """Alternative 6-rank decomposition when axial-only is required.

    (3, 2, 1): imax/3 = 20, jmax/2 = 110, kmax = 60.  More halo surface
    than (1, 5, 1) but gives the 6th rank.
    """
    params = build_V_3d_params(nodes=(3, 2, 1))
    assert params["nranks"] == 6


def test_V_3d_nodes_rejects_non_divisor():
    """Day-10 R4 safety: build fails fast if nodes do not divide the mesh."""
    # jmax = 220; 220 % 7 != 0 so (1,7,1) must raise.
    with pytest.raises(ValueError, match="does not divide mesh"):
        build_V_3d_params(nodes=(1, 7, 1))
    # Non-positive component must raise too.
    with pytest.raises(ValueError, match="positive integers"):
        build_V_3d_params(nodes=(0, 1, 1))


def test_V_3d_dem_nodes_pass_through(tmp_path):
    """DEM variant inherits the nodes kwarg from its PIC base."""
    params = build_V_3d_dem_params(nodes=(1, 4, 1), cgp_factor=10.0)
    out = tmp_path / "V_3d_dem.mfx"
    render_case(V3D_DEM_TEMPLATE, params, out)
    text = out.read_text()
    assert params["nodes_j"] == 4
    assert params["nranks"] == 4
    assert re.search(r"nodesj\s*=\s*4\b", text)


def test_V_3d_dem_default_cgp_solids_model(tmp_path):
    """Default (cgp_factor=10) selects CGP with cgp_stat_wt(1)=10."""
    text, params = _render_V3_dem(tmp_path, cgp_factor=10.0)
    assert params["solids_model"] == "CGP"
    assert params["cgp_stat_wt"] == pytest.approx(10.0)
    assert re.search(r"solids_model\(1\)\s*=\s*'CGP'", text)
    assert re.search(r"cgp_stat_wt\(1\)\s*=\s*10(\.0)?\b", text)


def test_V_3d_dem_pure_dem_option(tmp_path):
    """cgp_factor=1 selects pure DEM (no cgp_stat_wt rendered)."""
    text, params = _render_V3_dem(tmp_path, cgp_factor=1.0)
    assert params["solids_model"] == "DEM"
    assert params["cgp_stat_wt"] == pytest.approx(1.0)
    assert re.search(r"solids_model\(1\)\s*=\s*'DEM'", text)
    # Pure DEM must not emit the CGP stat-weight keyword as an active
    # input line (header comment is allowed; strip comments for the check).
    deck_only = "\n".join(
        ln for ln in text.splitlines()
        if not ln.lstrip().startswith("#!")
    )
    assert not re.search(r"^\s*cgp_stat_wt\(1\)\s*=", deck_only, re.MULTILINE)


def test_V_3d_dem_contact_keywords(tmp_path):
    """DEM contact block (kn, kt_fac, mew, des_en_*) must be present with
    values from params/geometry.yaml spouted_bed_V_3d_dem (sourced to
    NETL dem/spouted_bed_pseudo_2d).
    """
    text, params = _render_V3_dem(tmp_path)
    # Normal spring stiffness
    assert re.search(r"^\s*kn\s*=\s*1000(\.0)?\b", text, re.MULTILINE)
    assert re.search(r"^\s*kn_w\s*=\s*1000(\.0)?\b", text, re.MULTILINE)
    # 2/7 ratio for tangential/normal spring (Cundall & Strack 1979)
    m_kt = re.search(r"^\s*kt_fac\s*=\s*([0-9.eE+\-]+)", text, re.MULTILINE)
    assert m_kt is not None
    assert float(m_kt.group(1)) == pytest.approx(2.0 / 7.0, rel=1e-5)
    m_ktw = re.search(r"^\s*kt_w_fac\s*=\s*([0-9.eE+\-]+)", text, re.MULTILINE)
    assert float(m_ktw.group(1)) == pytest.approx(2.0 / 7.0, rel=1e-5)
    # Coulomb friction
    assert re.search(r"^\s*mew\s*=\s*0\.3\b", text, re.MULTILINE)
    assert re.search(r"^\s*mew_w\s*=\s*0\.3\b", text, re.MULTILINE)
    # Normal restitution (particle-particle and particle-wall) = 0.9
    # which also matches case_V.c_e in params/particles.yaml.
    assert re.search(r"des_en_input\(1\)\s*=\s*0\.9\b", text)
    assert re.search(r"des_en_wall_input\(1\)\s*=\s*0\.9\b", text)
    # Damping ratios
    assert re.search(r"des_etat_fac\s*=\s*0\.5\b", text)
    assert re.search(r"des_etat_w_fac\s*=\s*0\.5\b", text)
    # Dict mirrors
    assert params["kn"] == pytest.approx(1000.0)
    assert params["kt_fac"] == pytest.approx(2.0 / 7.0, rel=1e-5)
    assert params["mew"] == pytest.approx(0.3)
    assert params["des_en_input_pp"] == pytest.approx(0.9)
    assert params["des_en_input_pw"] == pytest.approx(0.9)
    assert params["des_etat_fac"] == pytest.approx(0.5)


def test_V_3d_dem_shared_des_block(tmp_path):
    """DEM template must share the gas->parcel interpolation block
    (des_interp_mean_fields / des_interp_on / des_interp_scheme /
    gener_part_config / des_epg_clip) with Day-10 PIC so the drag law
    sees the same interpolated gas velocity.
    """
    text, _ = _render_V3_dem(tmp_path)
    assert re.search(r"des_interp_mean_fields\s*=\s*\.True\.", text)
    assert re.search(r"des_interp_on\s*=\s*\.True\.", text)
    assert re.search(r"des_interp_scheme\s*=\s*'LINEAR_HAT'", text)
    assert re.search(r"gener_part_config\s*=\s*\.True\.", text)
    assert re.search(r"des_oneway_coupled\s*=\s*\.False\.", text)
    assert re.search(r"des_epg_clip\s*=\s*0\.42\b", text)


def test_V_3d_dem_no_pic_or_tfm_keywords(tmp_path):
    """PIC-only and TFM-only keywords must NOT appear in the DEM deck."""
    text, _ = _render_V3_dem(tmp_path)
    # PIC-only (Harris-Crighton + MP-PIC damping + PIC stat-weights)
    assert not re.search(r"\bfric_exp_pic\b", text)
    assert not re.search(r"\bfric_non_sing_fac\b", text)
    assert not re.search(r"\bpsfac_fric_pic\b", text)
    assert not re.search(r"\bmppic_coeff_en1\b", text)
    assert not re.search(r"\bmppic_coeff_en_wall\b", text)
    assert not re.search(r"\bmppic_coeff_et_wall\b", text)
    assert not re.search(r"\bmppic_velfac_coeff\b", text)
    assert not re.search(r"\bic_pic_const_statwt\b", text)
    assert not re.search(r"\bbc_pic_mi_const_statwt\b", text)
    # TFM-only
    assert not re.search(r"\bblending_function\b", text)
    assert not re.search(r"\bkt_type\b", text)
    assert not re.search(r"^\s*phi\s*=", text, re.MULTILINE)
    assert not re.search(r"^\s*c_e\s*=", text, re.MULTILINE)
    # The DEM template does not use the SCHAEFFER friction label; it is a
    # TFM/PIC-only closure flag.  Verify absence.
    assert not re.search(r"\bfriction_model\b", text)


def test_V_3d_dem_particle_count(tmp_path):
    """Hand calc: 54.5 g / (6050 * (pi/6) * (500e-6)^3) = 1.376e5 real
    particles.  With default cgp_factor=10 the parcel count is ~1.4e4.
    """
    import math
    _, params = _render_V3_dem(tmp_path, cgp_factor=10.0)
    expected_real = 0.0545 / (6050.0 * (math.pi / 6.0) * (5.0e-4) ** 3)
    assert params["n_real_particles"] == pytest.approx(expected_real, rel=1e-3)
    # 1.376e5 particles / 10 = ~1.376e4 parcels at the default coarsening.
    assert params["n_parcels_est"] == pytest.approx(expected_real / 10.0, rel=1e-3)
    assert 1.3e5 < params["n_real_particles"] < 1.5e5


def test_V_3d_dem_vtk_particle_id(tmp_path):
    """Residence-time tracking requires vtk_part_id=.True. on the
    Particles VTK dump (otherwise ids are not written and histograms
    cannot follow parcels across frames).
    """
    text, _ = _render_V3_dem(tmp_path)
    assert re.search(r"vtk_filebase\(2\)\s*=\s*'Particles'", text)
    assert re.search(r"vtk_data\(2\)\s*=\s*'P'", text)
    assert re.search(r"vtk_part_id\(2\)\s*=\s*\.True\.", text)
    assert re.search(r"vtk_part_vel\(2\)\s*=\s*\.True\.", text)
    assert re.search(r"vtk_part_diameter\(2\)\s*=\s*\.True\.", text)


def test_V_3d_dem_cgp_factor_validation(tmp_path):
    """cgp_factor < 1 should raise (prevents accidental fractional
    coarsening, which MFiX does not accept)."""
    with pytest.raises(ValueError):
        build_V_3d_dem_params(cgp_factor=0.5)


def test_V_3d_dem_shares_bc_and_ic_with_pic(tmp_path):
    """DEM template must reuse the Day-10 R3 PIC BC/IC layout: full-face
    MI at box bottom, PO at top, five ICs (whole-domain gas + 4 packed
    bed strips that leave a central gas column uncovered -- NETL's
    tutorials/pic/spouted_bed_3d quadrant pattern).  See
    notes/day10_pathology_review.md for why Day-10 R2's single-IC_2
    bed + explicit stub IC_3 was replaced.
    """
    text, params = _render_V3_dem(tmp_path)
    # BC_3 MI across the full box bottom face.
    m_xw = re.search(r"bc_x_w\(3\)\s*=\s*([-+eE.0-9]+)", text)
    m_xe = re.search(r"bc_x_e\(3\)\s*=\s*([-+eE.0-9]+)", text)
    assert float(m_xw.group(1)) == pytest.approx(params["x_min_m"])
    assert float(m_xe.group(1)) == pytest.approx(params["x_max_m"])
    # IC_2..IC_5 all cover y in [0, H_static] with eps_s = 0.57, HEXA lattice.
    R_col = params["R_bed_gas_col_m"]
    for ic in (2, 3, 4, 5):
        m_ys = re.search(rf"ic_y_s\({ic}\)\s*=\s*([-+eE.0-9]+)", text)
        m_yn = re.search(rf"ic_y_n\({ic}\)\s*=\s*([-+eE.0-9]+)", text)
        m_eps = re.search(rf"ic_ep_s\({ic},1\)\s*=\s*([-+eE.0-9]+)", text)
        m_lat = re.search(rf"ic_des_lattice\({ic},1\)\s*=\s*'HEXA'", text)
        assert m_ys is not None
        assert float(m_ys.group(1)) == pytest.approx(0.0, abs=1e-12)
        assert float(m_yn.group(1)) == pytest.approx(params["H_static_m"])
        assert float(m_eps.group(1)) == pytest.approx(0.57)
        assert m_lat is not None, f"IC_{ic} missing ic_des_lattice=HEXA"
    # The +z and -z strips span the full x extent; the center strips are
    # bounded by |z| <= R_col.
    assert float(re.search(r"ic_z_b\(2\)\s*=\s*([-+eE.0-9]+)", text).group(1)) == pytest.approx(R_col)
    assert float(re.search(r"ic_z_t\(3\)\s*=\s*([-+eE.0-9]+)", text).group(1)) == pytest.approx(-R_col)
    assert float(re.search(r"ic_x_e\(4\)\s*=\s*([-+eE.0-9]+)", text).group(1)) == pytest.approx(-R_col)
    assert float(re.search(r"ic_x_w\(5\)\s*=\s*([-+eE.0-9]+)", text).group(1)) == pytest.approx(R_col)
    # No IC_6: the inlet stub region is covered by IC_1, not an explicit IC.
    assert re.search(r"ic_x_w\(6\)", text) is None

"""Tests for Case C26 (Conry et al. 2026): builder, sweep matrix, rendered deck."""

from __future__ import annotations

import math
import re

import pytest

from src.make_case import REPO_ROOT, render_case
from src.make_case_c26 import (
    build_C26_3d_params,
    load_sweep,
    slpm_to_mass_flow,
    write_gridmap,
)

C26_TEMPLATE = REPO_ROOT / "cases" / "C26_3d" / "template_dem.mfx.j2"
V3D_DEM_TEMPLATE = REPO_ROOT / "cases" / "V_3d" / "template_dem.mfx.j2"

# Keywords introduced by the C26 template, each checked against the
# installed MFiX 26.1.2 source on 2026-10-16 (CLAUDE.md rule 1; file/line
# references in notes/c26_conry2026_case.md "Keyword verification").
C26_NEW_KEYWORDS = {
    "ps_x_w", "ps_x_e", "ps_y_s", "ps_y_n", "ps_z_b", "ps_z_t",
    "ps_massflow_g", "ps_u_g", "ps_v_g", "ps_w_g",
    "ic_psd_type", "ic_psd_mean_dp", "ic_psd_stdev",
    "ic_psd_min_dp", "ic_psd_max_dp",
    "vtk_t_start", "vtk_t_end",
    # 2026-10-17: particle inventory monitor (init_namelist.f 4715/4722;
    # check_output_control.f 888-892; tutorials/dem/hopper_3d).
    "monitor_pmass", "monitor_pcount",
}

# Ideal-gas Ar at 25 C, 1 atm: rho = p M / (R T).
_R = 8.314462618
_M_AR = 39.948e-3
RHO_AR_25C = 101325.0 * _M_AR / (_R * 298.15)          # = 1.6329 kg/m^3


def _case(case_id: str) -> dict:
    return next(c for c in load_sweep()["cases"] if c["id"] == case_id)


@pytest.fixture(scope="module")
def p01() -> dict:
    return build_C26_3d_params(_case("c01"))


def _keywords(text: str) -> set[str]:
    """Lower-case MFiX keyword names (indices stripped) in a deck."""
    out = set()
    for line in text.splitlines():
        s = line.strip()
        if not s or s.startswith("#") or s.startswith("!") or "=" not in s:
            continue
        name = s.split("=", 1)[0].strip()
        name = re.sub(r"\(.*\)$", "", name).strip().lower()
        if re.fullmatch(r"[a-z_][a-z0-9_]*", name):
            out.add(name)
    return out


# ---------------- unit conversions / hand values ----------------

def test_slpm_to_mass_flow_hand_value():
    # 10 L/min at 1.6329 kg/m^3 -> 1.6667e-4 m^3/s * 1.6329 = 2.7215e-4 kg/s
    assert slpm_to_mass_flow(10.0, RHO_AR_25C) == pytest.approx(2.7215e-4, rel=1e-4)


def test_slpm_to_mass_flow_rejects_nonpositive():
    with pytest.raises(ValueError):
        slpm_to_mass_flow(0.0, 1.6)
    with pytest.raises(ValueError):
        slpm_to_mass_flow(10.0, -1.0)


def test_argon_density_matches_ideal_gas(p01):
    assert p01["rho_g_kg_m3"] == pytest.approx(RHO_AR_25C, rel=2e-3)


def test_nozzle_velocity_hand_value(p01):
    # Q = 10 slpm = 1.6667e-4 m^3/s (reference state == gas state, 298 K, 1 atm)
    # A = pi (2.38125e-3 / 2)^2 = 4.4535e-6 m^2  ->  U = 37.42 m/s
    assert p01["U_nozzle_m_s"] == pytest.approx(37.42, rel=2e-3)
    assert p01["ps_v_g_m_s"] == pytest.approx(p01["U_nozzle_m_s"])
    assert p01["ps_massflow_g_kg_s"] == pytest.approx(
        p01["U_nozzle_m_s"] * math.pi * p01["R_nozzle_m"] ** 2 * p01["rho_g_kg_m3"],
        rel=1e-9)


@pytest.mark.parametrize("cid,u_expected", [
    ("c01", 37.42), ("c02", 44.91), ("c03", 52.39),
    ("c10", 82.33), ("c11", 97.30), ("c12", 112.27),
])
def test_nozzle_velocity_scales_with_flow(cid, u_expected):
    assert build_C26_3d_params(_case(cid))["U_nozzle_m_s"] == pytest.approx(
        u_expected, rel=2e-3)


# ---------------- mesh / geometry ----------------

def test_mesh_matches_c26(p01):
    # C26 2.2.1: 2.1 mm cells, 24 across the reactor diameter, 100 vertical.
    assert p01["dx_m"] == pytest.approx(2.1e-3)
    assert p01["D_c_m"] / p01["dx_m"] == pytest.approx(24.0)
    assert p01["H_dom_m"] == pytest.approx(0.210)
    assert (p01["imax"], p01["jmax"], p01["kmax"]) == (27, 104, 27)
    assert len(p01["dx_m_list"]) == 27 and len(p01["dy_m_list"]) == 104


def test_axis_is_a_cell_centre(p01):
    # Odd cell count: the vessel axis x = 0 is the centre of cell imax // 2.
    i_axis = p01["imax"] // 2
    x_centre = p01["x_min_m"] + (i_axis + 0.5) * p01["dx_m"]
    assert x_centre == pytest.approx(0.0, abs=1e-12)


def test_point_source_cell_is_inside_cone(p01):
    assert (p01["ps_x_m"], p01["ps_z_m"]) == (0.0, 0.0)
    assert p01["ps_y_m"] == pytest.approx(0.5 * p01["dx_m"])
    # the cone radius at the top of the PS row exceeds the half cell width
    assert p01["r_cone_at_ps_top_m"] > 0.5 * p01["dx_m"]


def test_grid_line_on_cone_cap(p01):
    # y = 0 (the cap) must be a grid line: y_min is an integer number of dx.
    n = -p01["y_min_m"] / p01["dx_m"]
    assert n == pytest.approx(round(n))


def test_cone_height(p01):
    expected = (p01["R_c_m"] - p01["R_t_m"]) / math.tan(math.radians(30.0))
    assert p01["h_cone_m"] == pytest.approx(expected)
    assert p01["h_cone_m"] == pytest.approx(0.04159, rel=1e-3)


def test_bed_heights_and_paper_annotation(p01):
    # 60 g YSZ packed at 0.57 in the 60-deg cone ~ 35 mm; a 90-deg cone would
    # give ~ 24 mm, closer to the C25/C26 Fig. 1 "~25 mm".  Logged, not tuned.
    assert 0.030 < p01["H_static_m"] < 0.040
    assert 0.020 < p01["H_static_if_90deg_m"] < 0.028
    assert p01["H_ic_m"] > p01["H_static_m"]      # loose seeding is taller
    assert p01["H_ic_m"] < p01["H_dom_m"]


def test_nodes_must_divide_mesh():
    with pytest.raises(ValueError):
        build_C26_3d_params(_case("c01"), nodes=(2, 1, 1))   # 27 is odd
    p = build_C26_3d_params(_case("c01"), nodes=(3, 8, 3))
    assert p["nranks"] == 72


# ---------------- time step ----------------

def test_stl_sample_fraction_override_changes_only_the_stl(p01):
    assert p01["dx_sample_m"] == pytest.approx(0.5 * p01["dx_m"])   # geometry.yaml
    p = build_C26_3d_params(_case("c01"), dx_sample_fraction=0.25)
    assert p["dx_sample_m"] == pytest.approx(0.25 * p01["dx_m"])
    for k in ("imax", "jmax", "kmax", "dx_m", "dt_max_s", "ps_cell_y_n_m"):
        assert p[k] == p01[k]


def test_dt_max_is_paper_step_capped_by_rule9(p01):
    # C26 p.9: 0.01 ms; recipe Rule 9: 0.5 dx / U_nozzle (2.8e-5 s at c01).
    assert p01["dt_max_paper_s"] == pytest.approx(1.0e-5)
    assert p01["dt_max_s"] == pytest.approx(1.0e-5)
    assert p01["dt_s"] <= p01["dt_max_s"]
    # SS at 30 slpm: Rule 9 (9.4e-6 s) is tighter than the paper step.
    p12 = build_C26_3d_params(_case("c12"))
    assert p12["dt_max_s"] == pytest.approx(0.5 * p12["dx_m"] / p12["U_nozzle_m_s"])
    assert p12["dt_max_s"] < p01["dt_max_s"]


# ---------------- particles / sweep matrix ----------------

def test_psd_from_measured_variance(p01):
    sigma = math.sqrt(149.1e-12)
    assert p01["psd_mean_dp_m"] == pytest.approx(411.6e-6)
    assert p01["psd_stdev_m"] == pytest.approx(sigma)
    assert p01["psd_min_dp_m"] == pytest.approx(411.6e-6 - 3 * sigma)
    assert p01["psd_max_dp_m"] == pytest.approx(411.6e-6 + 3 * sigma)


def test_sweep_matches_table2():
    cases = load_sweep()["cases"]
    assert [c["id"] for c in cases] == [f"c{i:02d}" for i in range(1, 13)]
    rows = [(c["feedstock"], c["flow_slpm"], c["stat_wt"], c["mew_w"]) for c in cases]
    ysz, ss = "case_C26_ysz", "case_C26_ss"
    assert rows == [
        (ysz, 10, 10, 0.1), (ysz, 12, 10, 0.1), (ysz, 14, 10, 0.1),
        (ysz, 10, 10, 0.3), (ysz, 10, 10, 0.5), (ysz, 10, 10, 0.7),
        (ysz, 10, 5, 0.1), (ysz, 10, 3, 0.1), (ysz, 10, 1, 0.1),
        (ss, 22, 1, 0.1), (ss, 26, 1, 0.1), (ss, 30, 1, 0.1),
    ]


@pytest.mark.parametrize("cid,model", [("c01", "CGP"), ("c08", "CGP"),
                                       ("c09", "DEM"), ("c10", "DEM")])
def test_solids_model_from_stat_wt(cid, model):
    assert build_C26_3d_params(_case(cid))["solids_model"] == model


def test_table1_closures(p01):
    assert p01["drag_type"] == "TGS"
    assert (p01["kn"], p01["kn_w"]) == (1000.0, 1000.0)
    assert (p01["kt_fac"], p01["des_etat_fac"]) == (0.285, 0.5)
    assert (p01["mew"], p01["mew_w"]) == (0.1, 0.1)
    assert (p01["des_en_input_pp"], p01["des_en_input_pw"]) == (0.5, 0.5)
    assert build_C26_3d_params(_case("c06"))["mew_w"] == 0.7


def test_parcel_count_estimate(p01):
    # 60 g / (6000 * pi/6 * E[d^3]) ~ 2.7e5 real particles; W = 10.
    assert p01["n_real_particles"] == pytest.approx(2.7e5, rel=0.03)
    assert p01["n_parcels_est"] == pytest.approx(p01["n_real_particles"] / 10.0)


def test_piv_window(p01):
    assert p01["piv_t_end_s"] == pytest.approx(p01["tstop_s"])
    assert p01["piv_t_end_s"] - p01["piv_t_start_s"] == pytest.approx(300 * 0.002)
    short = build_C26_3d_params(_case("c01"), tstop_s=0.3)
    assert short["piv_t_start_s"] == 0.0


# ---------------- rendered deck ----------------

def _render(tmp_path, cid="c01", **kw):
    params = build_C26_3d_params(_case(cid), **kw)
    out = tmp_path / f"{params['run_name']}.mfx"
    render_case(C26_TEMPLATE, params, out)
    return out.read_text(), params


def test_deck_has_point_source_and_no_inflow(tmp_path):
    text, p = _render(tmp_path)
    m = re.search(r"ps_massflow_g\(1\)\s*=\s*([0-9.eE+\-]+)", text)
    assert m and float(m.group(1)) == pytest.approx(p["ps_massflow_g_kg_s"], rel=1e-5)
    m = re.search(r"ps_v_g\(1\)\s*=\s*([0-9.eE+\-]+)", text)
    assert m and float(m.group(1)) == pytest.approx(p["U_nozzle_m_s"], rel=1e-5)
    assert "'MI'" not in text
    assert not re.search(r"^\s*bc_type\(3\)", text, flags=re.M)


def test_deck_drag_and_psd(tmp_path):
    text, _ = _render(tmp_path)
    assert re.search(r"drag_type\s*=\s*'TGS'", text)
    assert "drag_c1" not in text and "drag_d1" not in text
    assert re.search(r"ic_psd_type\(2,1\)\s*=\s*'NORMAL'", text)
    assert re.search(r"vtk_t_start\(3\)\s*=", text)


def test_deck_cgp_keyword_only_when_coarse_grained(tmp_path):
    text, _ = _render(tmp_path, "c01")
    assert re.search(r"solids_model\(1\)\s*=\s*'CGP'", text)
    assert re.search(r"cgp_stat_wt\(1\)\s*=\s*10", text)
    text9, _ = _render(tmp_path, "c09")
    assert re.search(r"solids_model\(1\)\s*=\s*'DEM'", text9)
    assert not re.search(r"^\s*cgp_stat_wt", text9, flags=re.M)


def test_deck_mesh_lines(tmp_path):
    text, p = _render(tmp_path)
    assert len(re.findall(r"^\s*dx\(\d+\)", text, flags=re.M)) == p["imax"]
    assert len(re.findall(r"^\s*dy\(\d+\)", text, flags=re.M)) == p["jmax"]
    assert len(re.findall(r"^\s*dz\(\d+\)", text, flags=re.M)) == p["kmax"]


def test_deck_point_source_is_one_axis_cell_box(tmp_path, p01):
    """A zero-width PS at a cell centre rounds to a face and lands in a
    corner cut cell (c01 smoke, 2026-10-16).  The region must be the axis
    cell's faces: width dx in x and z, centred on 0, one row in y."""
    text, _ = _render(tmp_path, "c01")
    v = {k: float(x) for k, x in re.findall(
        r"^\s*(ps_[xyz]_[wesnbt])\(1\)\s*=\s*(\S+)", text, flags=re.M)}
    dx = p01["dx_m"]
    assert v["ps_x_e"] - v["ps_x_w"] == pytest.approx(dx)
    assert v["ps_z_t"] - v["ps_z_b"] == pytest.approx(dx)
    assert v["ps_y_n"] - v["ps_y_s"] == pytest.approx(dx)
    assert v["ps_x_w"] + v["ps_x_e"] == pytest.approx(0.0, abs=1e-12)
    assert v["ps_z_b"] + v["ps_z_t"] == pytest.approx(0.0, abs=1e-12)
    # every bound is a grid line
    for key, origin in (("ps_x_w", p01["x_min_m"]), ("ps_y_s", p01["y_min_m"]),
                        ("ps_z_b", p01["z_min_m"])):
        n = (v[key] - origin) / dx
        assert n == pytest.approx(round(n), abs=1e-6)


def test_deck_keywords_are_verified(tmp_path):
    """CLAUDE.md rule 1: every keyword is either already in the verified
    V_3d DEM template or in the 2026-10-16 C26 verified list."""
    allowed = _keywords(V3D_DEM_TEMPLATE.read_text()) | C26_NEW_KEYWORDS
    for cid in ("c01", "c09", "c12"):
        text, _ = _render(tmp_path, cid)
        unknown = _keywords(text) - allowed
        assert not unknown, f"{cid}: unverified keywords {sorted(unknown)}"


# ---------------- non-uniform DMP partition (gridmap.dat) ----------------

GRIDMAP_5x1x5 = ([9, 3, 3, 3, 9], [], [9, 3, 3, 3, 9])


def test_gridmap_sets_nodes_and_skips_divisor_check():
    p = build_C26_3d_params(_case("c01"), gridmap=GRIDMAP_5x1x5)
    # 5 does not divide 27: allowed only because gridmap gives exact sizes.
    assert (p["nodes_i"], p["nodes_j"], p["nodes_k"], p["nranks"]) == (5, 1, 5, 25)
    assert p["gridmap"] == {"i": [9, 3, 3, 3, 9], "j": [104], "k": [9, 3, 3, 3, 9]}
    assert build_C26_3d_params(_case("c01"))["gridmap"] is None


def test_gridmap_validation():
    with pytest.raises(ValueError, match="sum to"):
        build_C26_3d_params(_case("c01"), gridmap=([9, 3, 3, 3, 8], [], [27]))
    with pytest.raises(ValueError, match=">= 3 cells"):
        build_C26_3d_params(_case("c01"), gridmap=([10, 2, 3, 3, 9], [], [27]))
    with pytest.raises(ValueError, match="conflicts"):
        build_C26_3d_params(_case("c01"), nodes=(3, 1, 3), gridmap=GRIDMAP_5x1x5)


def test_gridmap_central_slab_holds_the_axis(p01):
    # The narrow slabs must straddle the vessel axis (cell imax // 2).
    sizes = GRIDMAP_5x1x5[0]
    i_axis = p01["imax"] // 2
    edges = [sum(sizes[:n]) for n in range(len(sizes) + 1)]
    owner = next(n for n in range(len(sizes)) if edges[n] <= i_axis < edges[n + 1])
    assert sizes[owner] == 3 and 0 < owner < len(sizes) - 1


def test_write_gridmap_format(tmp_path):
    """dmp_modules/gridmap_mod.f PARTITION: 'NI NJ NK', then '<rank> <cells>'
    per rank along i, j, k (list-directed)."""
    p = build_C26_3d_params(_case("c01"), gridmap=GRIDMAP_5x1x5)
    lines = write_gridmap(tmp_path, p["gridmap"]).read_text().split("\n")
    assert lines[0] == "5 1 5"
    assert lines[1:6] == ["0 9", "1 3", "2 3", "3 3", "4 9"]
    assert lines[6] == "0 104"
    assert lines[7:12] == ["0 9", "1 3", "2 3", "3 3", "4 9"]

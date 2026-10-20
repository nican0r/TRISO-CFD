"""Parameter builder for Case C26 (Conry et al. 2026 ORNL benchtop rig).

Case C26 reproduces the MFiX CFD-CGDEM model of

    Conry, Dal Forno Chuahy, Oyedeji, Finney, Heldt, Lopez-Honorato (2026),
    Nucl. Eng. Des. 456: 115009 ("C26"),

on the project's MFiX 26.1.2.  It is a separate case from Case V (different
rig, gas, particles, inlet and validation data); see
notes/c26_conry2026_case.md for the design and for every deliberate
deviation from .claude/skills/mfix-spouted-bed-recipe.

Public API:

    load_sweep() -> dict
    slpm_to_mass_flow(flow_slpm, rho_ref_kg_m3) -> float
    resolve_particle(case, particles) -> dict
    build_C26_3d_params(case, tstop_s=None, run_name=None, ...) -> dict
    write_gridmap(run_dir, gridmap) -> Path

All values SI (CLAUDE.md rule 3).  Every physical number comes from
params/*.yaml (CLAUDE.md rule 2).
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Mapping, Sequence

from src.correlations import coated_particle_stages
from src.make_case import (
    PARAMS_DIR,
    _bed_height_from_inventory,
    _load_yaml,
    _uniform_axis,
)

SLPM_TO_M3_S = 1.0e-3 / 60.0   # unit conversion only: 1 L/min in m^3/s


def load_sweep() -> dict:
    """Return params/sweep_C26.yaml (flow / coating-stage / charge matrix)."""
    return _load_yaml(PARAMS_DIR / "sweep_C26.yaml")


def slpm_to_mass_flow(flow_slpm: float, rho_ref_kg_m3: float) -> float:
    """Mass flow [kg/s] for a standard volumetric flow [slpm].

    rho_ref_kg_m3 is the gas density at the MFC's standard reference state
    (params/gases.yaml slpm_reference_C26).
    """
    if flow_slpm <= 0.0:
        raise ValueError(f"flow_slpm must be > 0 (got {flow_slpm})")
    if rho_ref_kg_m3 <= 0.0:
        raise ValueError(f"rho_ref_kg_m3 must be > 0 (got {rho_ref_kg_m3})")
    return float(flow_slpm) * SLPM_TO_M3_S * float(rho_ref_kg_m3)


def _mean_cube_normal(mu: float, sigma: float) -> float:
    """E[d^3] for a normal size distribution: mu^3 + 3 mu sigma^2."""
    return mu ** 3 + 3.0 * mu * sigma ** 2


# Smallest per-rank slab in a gridmap.dat partition [cells].  A uniform
# 9-way split of the 27-cell box already gives 3-cell slabs; MFiX keeps 2
# ghost layers per side, so thinner slabs would be mostly ghost cells.
GRIDMAP_MIN_CELLS = 3


def write_gridmap(run_dir: Path, gridmap: Mapping[str, Sequence[int]]) -> Path:
    """Write MFiX's gridmap.dat (dmp_modules/gridmap_mod.f, PARTITION).

    Format, read list-directed on PE_IO: ``NODESI NODESJ NODESK``, then one
    ``<rank index> <cell count>`` line per rank along i, then j, then k.
    MFiX sizes its partition arrays from the deck's nodesi/j/k, so the deck
    must carry the same rank counts (the builder sets them from the gridmap).
    """
    sizes = [list(gridmap[a]) for a in "ijk"]
    lines = [" ".join(str(len(s)) for s in sizes)]
    for s in sizes:
        lines += [f"{n} {c}" for n, c in enumerate(s)]
    path = Path(run_dir) / "gridmap.dat"
    path.write_text("\n".join(lines) + "\n")
    return path


def resolve_particle(case: Mapping[str, object], particles: Mapping[str, object]) -> dict:
    """Particle properties for one sweep row.

    The row's ``feedstock`` names a params/particles.yaml block (C26 YSZ).
    An optional ``coating_stage`` names the last TRISO layer applied to it
    (particles.yaml ``triso_coating_layers``, deposition order); every layer
    up to and including that one is added by
    :func:`src.correlations.coated_particle_stages`.  The PSD standard
    deviation of the feedstock is kept (uniform layers, assumed).

    Returns label, coating_stage (None when bare), rho_p_kg_m3, d_p_mean_m,
    d_p_sigma_m, psd_type, psd_n_sigma and the feedstock's inventory_g.
    """
    feed_key = str(case["feedstock"])
    if feed_key not in particles:
        raise KeyError(f"feedstock {feed_key!r} not in params/particles.yaml")
    feed = particles[feed_key]
    rho_p = float(feed["rho_p_kg_per_m3"])
    d_mean = float(feed["d_p_mean_m"])
    label = str(feed.get("label", feed_key))
    stage = case.get("coating_stage")
    if stage is not None:
        stage = str(stage)
        layers = particles["triso_coating_layers"]
        names = [str(layer["name"]) for layer in layers]
        if stage not in names:
            raise KeyError(f"coating_stage {stage!r} not in triso_coating_layers {names}")
        applied = layers[: names.index(stage) + 1]
        stages = coated_particle_stages(
            d_mean, rho_p,
            [(layer["name"], float(layer["t_m"]), float(layer["rho_kg_per_m3"]))
             for layer in applied],
        )
        _, d_mean, rho_p = stages[-1]
        label = f"{label} + " + " + ".join(str(layer["name"]) for layer in applied)
    return {
        "feedstock": feed_key,
        "label": label,
        "coating_stage": stage,
        "rho_p_kg_m3": rho_p,
        "d_p_mean_m": d_mean,
        "d_p_sigma_m": math.sqrt(float(feed["d_p_var_m2"])),
        "psd_type": str(feed["psd_type"]),
        "psd_n_sigma": float(feed["psd_n_sigma"]),
        "inventory_g": float(feed["inventory_g"]),
    }


def _gas_props(T_K: float, p_Pa: float, composition: Mapping[str, float]):
    from src.gasprops import gas_props

    return gas_props(T_K, p_Pa, composition)


def build_C26_3d_params(
    case: Mapping[str, object],
    tstop_s: float | None = None,
    run_name: str | None = None,
    description: str = "",
    stl_filename: str = "geometry_0001.stl",
    nodes: tuple[int, int, int] = (1, 1, 1),
    gridmap: tuple[Sequence[int], Sequence[int], Sequence[int]] | None = None,
    dx_sample_fraction: float | None = None,
    ps_momentum: bool | None = None,
    epg_clip: float | None = None,
    interp_scheme: str | None = None,
) -> dict:
    """Assemble the template params for one row of params/sweep_C26.yaml.

    case     : one entry of sweep_C26.yaml ``cases`` -- keys ``id``,
               ``feedstock`` (particles.yaml block), ``flow_slpm``,
               ``inventory_g`` (bed charge [g]; default the feedstock's),
               optional ``coating_stage`` (last TRISO layer on the
               feedstock, see :func:`resolve_particle`) and optional
               ``stat_wt`` (CGP statistical weight W; default and lower
               bound sweep ``base.stat_wt`` = 10).  The
               particle-wall friction factor is geometry.yaml ``mew_w``.
    tstop_s  : simulated end time [s]; default sweep ``base.tstop_s``.
    run_name : MFiX run_name; default ``C26_<id>``.
    nodes    : DMP decomposition (NI, NJ, NK); each must divide its axis.
    gridmap  : optional non-uniform DMP partition, as cell counts per rank
               along i, j and k (e.g. ([9,3,3,3,9], [104], [9,3,3,3,9])).
               Overrides ``nodes`` (which must then match its lengths or be
               left at the default). Written to gridmap.dat by
               :func:`write_gridmap`.
    dx_sample_fraction : STL SDF sampling step as a fraction of dx; default
               geometry.yaml ``dx_sample_fraction``.  Changes only the STL
               facet density (and so the cut-cell intersections), never the
               fluid mesh.
    ps_momentum : True gives the point source ps_v_g = U_nozzle (momentum
               flux m_dot * U_nozzle); False gives ps_v_g = 0, a mass-only
               source.  Default geometry.yaml ``ps_momentum``.
    epg_clip : DES_EPG_CLIP override (gas-fraction floor, 0-1); default
               geometry.yaml ``des_epg_clip``.  Omitting the keyword (MFiX's
               default, no floor) is not offered: on this mesh it fails with
               Error 1100 in sliver cut cells under both LINEAR_HAT and
               SQUARE_DPVM (run_log 2026-10-19, c01_smoke_r5_*).
    interp_scheme : DES_INTERP_SCHEME override ('LINEAR_HAT' or
               'SQUARE_DPVM'); default geometry.yaml ``des_interp_scheme``.
               SQUARE_DPVM also writes DES_INTERP_WIDTH =
               ``des_interp_width_cells`` * dx.

    Geometry (bed_C26, spouted_bed_C26_3d):
      - vessel: cone from a closed throat (r = R_t, y = 0) to the column
        (r = R_c), column up to the pressure outlet at y = H_dom;
      - uniform mesh dx = D_c / n_cells_across_D (2.1 mm); an ODD number of
        cells across the box so the vessel axis is a cell centre;
      - gas enters through ONE point-source cell centred on the axis in the
        first fluid row above the cap (C26 2.2.1).  ps_v_g is 0 for the
        default mass-only source, or the nozzle exit velocity (momentum flux
        m_dot * U_nozzle) with ``ps_momentum``.
    """
    geometry = _load_yaml(PARAMS_DIR / "geometry.yaml")
    particles = _load_yaml(PARAMS_DIR / "particles.yaml")
    gases = _load_yaml(PARAMS_DIR / "gases.yaml")
    numerics = _load_yaml(PARAMS_DIR / "numerics.yaml")["spouted_bed_C26_3d_numerics"]
    sweep = load_sweep()

    bed = geometry["bed_C26"]
    mesh = geometry["spouted_bed_C26_3d"]
    dem = geometry["spouted_bed_C26_3d_dem"]
    scheme = str(interp_scheme or dem["des_interp_scheme"])
    if epg_clip is not None and not 0.0 < float(epg_clip) < 1.0:
        raise ValueError(f"epg_clip must be in (0, 1), got {epg_clip}")
    if scheme not in ("LINEAR_HAT", "SQUARE_DPVM"):
        raise ValueError(f"interp_scheme must be LINEAR_HAT or SQUARE_DPVM, got {scheme!r}")
    base = sweep["base"]

    case_id = str(case["id"])
    particle = resolve_particle(case, particles)
    feed_key = particle["feedstock"]
    flow_slpm = float(case["flow_slpm"])
    # base.stat_wt is both the default W and its lower bound: larger
    # particles get larger parcels, never a smaller W.
    stat_wt_min = float(base["stat_wt"])
    stat_wt = float(case.get("stat_wt", stat_wt_min))
    mew_w = float(dem["mew_w"])
    if stat_wt < stat_wt_min:
        raise ValueError(
            f"stat_wt must be >= {stat_wt_min:g} (got {stat_wt}); larger particles "
            "use larger parcels, not a smaller W")
    if not (0.0 <= mew_w <= 1.0):
        raise ValueError(f"mew_w must be in [0, 1] (got {mew_w})")
    tstop_s = float(base["tstop_s"] if tstop_s is None else tstop_s)
    if tstop_s <= 0.0:
        raise ValueError(f"tstop_s must be > 0 (got {tstop_s})")
    run_name = run_name or f"C26_{case_id}"

    # ---- geometry --------------------------------------------------------
    D_c_m = float(bed["D_c_m"])
    R_c_m = 0.5 * D_c_m
    R_n_m = 0.5 * float(bed["D_nozzle_m"])
    R_t_m = 0.5 * float(bed["D_throat_m"])
    cone_half_angle_deg = 0.5 * float(bed["cone_angle_deg"])
    tan_t = math.tan(math.radians(cone_half_angle_deg))
    h_cone_m = (R_c_m - R_t_m) / tan_t

    n_across = int(mesh["n_cells_across_D"])
    n_vert = int(mesh["n_cells_vertical"])
    n_margin_xz = int(mesh["n_margin_cells_xz"])
    n_below = int(mesh["n_margin_cells_below"])
    if n_across < 2 or n_vert < 2 or n_margin_xz < 0 or n_below < 1:
        raise ValueError("invalid spouted_bed_C26_3d cell counts")
    dx_m = D_c_m / n_across
    H_dom_m = n_vert * dx_m
    if h_cone_m >= H_dom_m:
        raise ValueError(f"cone height {h_cone_m} m exceeds H_dom {H_dom_m} m")

    # Odd cell count across the box -> axis at a cell centre.
    n_xz = n_across + 2 * n_margin_xz + (1 if n_across % 2 == 0 else 0)
    x_half_m = 0.5 * n_xz * dx_m
    dx_m_list = _uniform_axis(dx_m, 2.0 * x_half_m)
    dz_m_list = _uniform_axis(dx_m, 2.0 * x_half_m)
    y_min_m = -n_below * dx_m
    dy_m_list = _uniform_axis(dx_m, H_dom_m - y_min_m)
    imax, jmax, kmax = len(dx_m_list), len(dy_m_list), len(dz_m_list)
    if imax != n_xz or kmax != n_xz or jmax != n_vert + n_below:
        raise AssertionError("mesh construction mismatch")

    if gridmap is not None:
        if len(gridmap) != 3:
            raise ValueError("gridmap needs cell counts for i, j and k")
        # An empty axis means one rank spanning it.
        gridmap = tuple([int(n) for n in sizes] if sizes else [n_cells]
                        for sizes, n_cells in zip(gridmap, (imax, jmax, kmax)))
        lengths = tuple(len(sizes) for sizes in gridmap)
        if tuple(nodes) not in ((1, 1, 1), lengths):
            raise ValueError(f"nodes={tuple(nodes)} conflicts with gridmap ranks {lengths}")
        for axis, sizes, n_cells in zip("ijk", gridmap, (imax, jmax, kmax)):
            if sum(sizes) != n_cells:
                raise ValueError(
                    f"gridmap {axis} sizes {sizes} sum to {sum(sizes)}, not {n_cells}")
            if min(sizes) < GRIDMAP_MIN_CELLS:
                raise ValueError(
                    f"gridmap {axis} sizes {sizes}: each rank needs >= "
                    f"{GRIDMAP_MIN_CELLS} cells")
        nodes = lengths
    ni, nj, nk = (int(nodes[0]), int(nodes[1]), int(nodes[2]))
    if min(ni, nj, nk) < 1:
        raise ValueError(f"nodes must be positive integers; got {nodes}")
    if gridmap is None and (imax % ni or jmax % nj or kmax % nk):
        raise ValueError(
            f"nodes=({ni},{nj},{nk}) does not divide mesh "
            f"(imax={imax}, jmax={jmax}, kmax={kmax})"
        )

    # Point-source cell (one cell, centred on the axis).
    ps_row = int(mesh["ps_row_offset_cells"])
    if ps_row < 0:
        raise ValueError("ps_row_offset_cells must be >= 0")
    ps_y_s_m = ps_row * dx_m
    ps_y_n_m = (ps_row + 1) * dx_m
    ps_y_m = 0.5 * (ps_y_s_m + ps_y_n_m)
    r_cone_at_ps_top = R_t_m + ps_y_n_m * tan_t

    # ---- gas and operating point ----------------------------------------
    gas = gases["cases"][str(base["gas_case"])]
    ref = gases["slpm_reference_C26"]
    rho_g, mu_g = _gas_props(gas["T_K"], gases["pressure_Pa"], gas["composition"])
    rho_ref, _ = _gas_props(ref["T_K"], ref["p_Pa"], ref["composition"])
    m_dot = slpm_to_mass_flow(flow_slpm, rho_ref)
    Q_actual = m_dot / rho_g
    A_nozzle = math.pi * R_n_m ** 2
    U_nozzle = Q_actual / A_nozzle
    U_sup = Q_actual / (math.pi * R_c_m ** 2)
    U_ps_cell = Q_actual / (dx_m * dx_m)   # mean velocity over the PS cell face

    # ---- time step: the paper's step, capped by recipe Rule 9 per case ----
    dt_max_paper_s = float(numerics["dt_max_paper_s"])
    dt_max_s = min(dt_max_paper_s, float(numerics["cfl_target"]) * dx_m / U_nozzle)
    dt_s = float(numerics["dt_init_fraction"]) * dt_max_s

    # ---- particles -------------------------------------------------------
    rho_p = particle["rho_p_kg_m3"]
    d_mean = particle["d_p_mean_m"]
    sigma = particle["d_p_sigma_m"]
    n_sig = particle["psd_n_sigma"]
    d_min = d_mean - n_sig * sigma
    d_max = d_mean + n_sig * sigma
    if d_min <= 0.0:
        raise ValueError(f"{feed_key}: mean - n_sigma*sigma <= 0")
    m_bed = float(case.get("inventory_g", particle["inventory_g"])) * 1.0e-3
    if m_bed <= 0.0:
        raise ValueError(f"inventory_g must be > 0 (got {m_bed * 1e3} g)")
    eps_s_ic = float(dem["eps_s_ic"])
    eps_s_packed = float(dem["eps_s_packed"])

    def _h(eps: float, half_deg: float) -> float:
        return _bed_height_from_inventory(
            mass_kg=m_bed, rho_p_kg_m3=rho_p, eps_s_bed=eps, R_i_m=R_t_m,
            cone_half_angle_deg=half_deg, R_c_m=R_c_m,
        )

    H_ic_m = _h(eps_s_ic, cone_half_angle_deg)
    H_static_m = _h(eps_s_packed, cone_half_angle_deg)
    # Cross-check of the C25/C26 Fig. 1 "~25 mm" annotation (not a model input).
    H_static_if_90deg_m = _h(eps_s_packed, 45.0)
    if H_ic_m >= H_dom_m:
        raise ValueError(f"IC region top H_ic {H_ic_m} m exceeds H_dom {H_dom_m} m")

    n_real = m_bed / (rho_p * math.pi / 6.0 * _mean_cube_normal(d_mean, sigma))
    solids_model = "CGP"   # W >= base.stat_wt (10): always coarse-grained
    m_parcel = stat_wt * rho_p * math.pi / 6.0 * d_mean ** 3
    d_parcel = d_mean * stat_wt ** (1.0 / 3.0)   # CGP parcel diameter [m]
    dt_rayleigh_s = math.pi * math.sqrt(m_parcel / float(dem["kn"]))

    # ---- output windows --------------------------------------------------
    piv_dt = float(numerics["piv_frame_dt_s"])
    piv_n = int(numerics["piv_n_frames"])
    piv_t_end = tstop_s
    piv_t_start = max(0.0, tstop_s - piv_n * piv_dt)

    half = 0.5 * dx_m
    stage_txt = particle["coating_stage"] or "bare"
    return {
        "case_id": case_id,
        "feedstock": feed_key,
        "feedstock_label": particle["label"],
        "coating_stage": particle["coating_stage"],
        "description": description or (
            f"C26 {case_id}: {feed_key} ({stage_txt}), {m_bed * 1e3:g} g, "
            f"{flow_slpm:g} slpm Ar (Conry et al. 2026 model)"
        ),
        "run_name": run_name,
        "run_type": "new",
        "tstop_s": tstop_s,
        "dt_s": dt_s,
        "dt_max_s": dt_max_s,
        "dt_max_paper_s": dt_max_paper_s,
        "dt_min_s": float(numerics["dt_min_s"]),
        "max_nit": int(numerics["max_nit"]),
        "nodes_i": ni, "nodes_j": nj, "nodes_k": nk, "nranks": ni * nj * nk,
        "gridmap": None if gridmap is None else {
            "i": gridmap[0], "j": gridmap[1], "k": gridmap[2]},
        # mesh
        "imax": imax, "jmax": jmax, "kmax": kmax,
        "n_cells_total": imax * jmax * kmax,
        "dx_m": dx_m,
        "dx_m_list": dx_m_list, "dy_m_list": dy_m_list, "dz_m_list": dz_m_list,
        "x_min_m": -x_half_m, "x_max_m": x_half_m,
        "z_min_m": -x_half_m, "z_max_m": x_half_m,
        "y_min_m": y_min_m, "H_dom_m": H_dom_m,
        # vessel / STL
        "stl_filename": stl_filename,
        "D_c_m": D_c_m, "R_c_m": R_c_m, "R_t_m": R_t_m, "R_nozzle_m": R_n_m,
        "cone_half_angle_deg": cone_half_angle_deg,
        "h_cone_m": h_cone_m,
        "dx_sample_m": float(mesh["dx_sample_fraction"] if dx_sample_fraction is None
                             else dx_sample_fraction) * dx_m,
        "top_extend_m": float(mesh["top_extend_m"]),
        # point source
        "ps_x_m": 0.0, "ps_y_m": ps_y_m, "ps_z_m": 0.0,
        "ps_cell_x_w_m": -half, "ps_cell_x_e_m": half,
        "ps_cell_y_s_m": ps_y_s_m, "ps_cell_y_n_m": ps_y_n_m,
        "ps_cell_z_b_m": -half, "ps_cell_z_t_m": half,
        "ps_row_offset_cells": ps_row,
        "r_cone_at_ps_top_m": r_cone_at_ps_top,
        "ps_massflow_g_kg_s": m_dot,
        "ps_momentum": bool(mesh["ps_momentum"] if ps_momentum is None else ps_momentum),
        "ps_v_g_m_s": U_nozzle if (mesh["ps_momentum"] if ps_momentum is None else ps_momentum)
                      else 0.0,
        # operating point (provenance)
        "flow_slpm": flow_slpm,
        "Q_actual_m3_s": Q_actual,
        "U_nozzle_m_s": U_nozzle,
        "U_sup_m_s": U_sup,
        "U_ps_cell_m_s": U_ps_cell,
        "rho_g_kg_m3": float(rho_g), "mu_g_pa_s": float(mu_g),
        "rho_ref_kg_m3": float(rho_ref),
        "gas_case": str(base["gas_case"]),
        "T_g_K": float(gas["T_K"]),
        "p_g_Pa": float(gases["pressure_Pa"]),
        # solids
        "solids_model": solids_model,
        "cgp_stat_wt": stat_wt,
        "rho_p_kg_m3": rho_p,
        "d_p_m": d_mean,
        "psd_type": particle["psd_type"],
        "psd_mean_dp_m": d_mean, "psd_stdev_m": sigma,
        "psd_min_dp_m": d_min, "psd_max_dp_m": d_max,
        "inventory_kg": m_bed,
        "eps_s_ic": eps_s_ic, "ep_g_ic": 1.0 - eps_s_ic,
        "H_ic_m": H_ic_m,
        "eps_s_packed": eps_s_packed,
        "H_static_m": H_static_m,
        "H_static_if_90deg_m": H_static_if_90deg_m,
        "H_mid_m": 0.5 * H_static_m,
        "n_real_particles": n_real,
        "n_parcels_est": n_real / stat_wt,
        "d_parcel_m": d_parcel,
        "dx_over_d_parcel": dx_m / d_parcel,
        "dt_rayleigh_s": dt_rayleigh_s,
        # DEM closure (C26 Table 1; mew_w fixed at the Table 1 baseline)
        "drag_type": str(dem["drag_type"]),
        "kn": float(dem["kn"]), "kn_w": float(dem["kn_w"]),
        "kt_fac": float(dem["kt_fac"]), "kt_w_fac": float(dem["kt_w_fac"]),
        "mew": float(dem["mew"]), "mew_w": mew_w,
        "des_en_input_pp": float(dem["des_en_input_pp"]),
        "des_en_input_pw": float(dem["des_en_input_pw"]),
        "des_etat_fac": float(dem["des_etat_fac"]),
        "des_etat_w_fac": float(dem["des_etat_w_fac"]),
        "des_interp_scheme": scheme,
        "des_interp_width_m": (float(dem["des_interp_width_cells"]) * dx_m
                               if scheme == "SQUARE_DPVM" else None),
        "des_epg_clip": float(dem["des_epg_clip"] if epg_clip is None else epg_clip),
        # output
        "monitor_dt_s": float(numerics["monitor_dt_s"]),
        "vtk_dt_background_s": float(numerics["vtk_dt_background_s"]),
        "vtk_dt_particles_s": float(numerics["vtk_dt_particles_s"]),
        "piv_frame_dt_s": piv_dt,
        "piv_n_frames": piv_n,
        "piv_t_start_s": piv_t_start,
        "piv_t_end_s": piv_t_end,
    }

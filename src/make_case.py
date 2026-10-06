"""Render an MFiX .mfx case file from a Jinja2 template plus a parameter dict.

The template is expected to use only MFiX 26.1.2 keywords that already appear in
``notes/mfx_anatomy.md`` (CLAUDE.md rule 1 - no invented keywords).  All physical
quantities in the ``params`` dict are SI; the template simply substitutes them.

Public API:

    render_case(template_path, params, out_mfx_path) -> Path

    build_fb_sweep_params(U_in, imax, jmax, tstop_s, run_name, description="") -> dict
"""

from __future__ import annotations

import argparse
import copy
from pathlib import Path
from typing import Mapping

import yaml
from jinja2 import Environment, FileSystemLoader, StrictUndefined

REPO_ROOT = Path(__file__).resolve().parent.parent
PARAMS_DIR = REPO_ROOT / "params"


def _load_yaml(path: Path) -> dict:
    with path.open("r") as f:
        return yaml.safe_load(f)


def render_case(
    template_path: Path, params: Mapping[str, object], out_mfx_path: Path
) -> Path:
    """Render ``template_path`` with ``params`` and write to ``out_mfx_path``.

    Uses StrictUndefined so a missing template variable is a hard error, not a
    silent empty string in a .mfx file.  Returns the output path.
    """
    template_path = Path(template_path)
    out_mfx_path = Path(out_mfx_path)
    env = Environment(
        loader=FileSystemLoader(str(template_path.parent)),
        undefined=StrictUndefined,
        keep_trailing_newline=True,
    )
    tmpl = env.get_template(template_path.name)
    out_mfx_path.parent.mkdir(parents=True, exist_ok=True)
    out_mfx_path.write_text(tmpl.render(**dict(params)))
    return out_mfx_path


def build_fb_sweep_params(
    U_in_m_s: float,
    imax: int,
    jmax: int,
    tstop_s: float,
    run_name: str,
    description: str = "",
) -> dict:
    """Assemble a template parameter dict for a single fb_sweep run.

    Reads particle, gas and geometry inputs from the YAML files under
    ``params/`` (CLAUDE.md rule 2 - no magic numbers in code).  Gas density
    and viscosity for the fb_sweep gas case are computed at runtime via
    ``src.gasprops.gas_props`` (Cantera), so switching gas condition does not
    require editing any constants here.

    U_in_m_s   : inlet superficial gas velocity [m/s]
    imax, jmax : grid cells in x and y [-]
    tstop_s    : run end time [s]
    run_name   : MFiX run_name (also used for VTK filebase/output naming)
    description: free-form description string for the .mfx header
    """
    from src.gasprops import gas_props  # local import to keep CLI light

    particles = _load_yaml(PARAMS_DIR / "particles.yaml")
    geometry = _load_yaml(PARAMS_DIR / "geometry.yaml")
    gases = _load_yaml(PARAMS_DIR / "gases.yaml")

    case_D = particles["case_D"]
    fb = particles["fb_sweep"]
    rig = geometry["fb_sweep_2d"]
    gas_case = gases["cases"][fb["gas_case"]]
    rho_g, mu_g = gas_props(
        gas_case["T_K"], gases["pressure_Pa"], gas_case["composition"]
    )

    ep_g_bed = 1.0 - float(fb["eps_s_bed"])

    return {
        "description": description or f"fb_sweep U={U_in_m_s:.4f} m/s",
        "run_name": run_name,
        "tstop_s": float(tstop_s),
        "dt_s": float(fb["dt_init_s"]),
        "imax": int(imax),
        "jmax": int(jmax),
        "width_m": float(rig["width_m"]),
        "height_m": float(rig["height_m"]),
        "H_static_m": float(fb["H_static_m"]),
        "ep_g_bed": ep_g_bed,
        "eps_s_bed": float(fb["eps_s_bed"]),
        "ep_star": float(case_D["ep_star"]),
        "rho_g_kg_m3": float(rho_g),
        "mu_g_pa_s": float(mu_g),
        "d_p_m": float(case_D["d_p_m"]),
        "rho_p_kg_m3": float(case_D["rho_p_kg_per_m3"]),
        "U_in_m_s": float(U_in_m_s),
        "vtk_dt_s": 0.05,
    }


def _common_fluid_solid_params(
    gas_case: str = "air_25C", particle_case: str = "case_V_rpt"
) -> dict:
    """Shared particle + gas properties for the Case V (V_2d/V_3d) templates.

    Defaults to the Missouri S&T RPT validation inventory: 2.18 mm glass beads
    (case_V_rpt) in dry air at 298 K / 101 kPa (air_25C). Returns a dict with
    d_p_m, rho_p_kg_m3, rho_g_kg_m3, mu_g_pa_s suitable for direct substitution
    into the V_2d/V_3d templates.
    """
    from src.gasprops import gas_props

    particles = _load_yaml(PARAMS_DIR / "particles.yaml")
    gases = _load_yaml(PARAMS_DIR / "gases.yaml")
    part = particles[particle_case]
    gas = gases["cases"][gas_case]
    rho_g, mu_g = gas_props(gas["T_K"], gases["pressure_Pa"], gas["composition"])
    return {
        "d_p_m": float(part["d_p_m"]),
        "rho_p_kg_m3": float(part["rho_p_kg_per_m3"]),
        "rho_g_kg_m3": float(rho_g),
        "mu_g_pa_s": float(mu_g),
    }


def build_V_2d_params(
    U_in_m_s: float,
    tstop_s: float,
    run_name: str,
    description: str = "",
    drag_type: str = "SYAM_OBRIEN",
    c_e: float = 0.95,
) -> dict:
    """Assemble a template parameter dict for the Case V 2D quadric spouted bed.

    All physical values come from params/geometry.yaml (spouted_bed_V_2d and
    bed_0p076m blocks), params/particles.yaml (case_D, spouted_bed_V_initial),
    and params/gases.yaml (air_20C via Cantera).

    V_2d uses Cartesian coordinates + cartesian_grid cut-cells + a single
    Y_CONE quadric for the cone wall + no_k=.True. (true 2D).  No STL file
    is needed - the vessel geometry is defined entirely by the quadric plus
    the natural x = +/- R_c domain walls above the cone section.

    U_in_m_s   : superficial gas velocity at the orifice [m/s]
    tstop_s    : run end time [s]
    run_name   : MFiX run_name (also used for VTK/monitor file base names)
    description: free-form description string for the .mfx header
    drag_type  : MFiX drag_type keyword (default SYAM_OBRIEN as per Day-8)
    c_e        : particle-particle restitution coefficient [-] (default 0.95)
    """
    geometry = _load_yaml(PARAMS_DIR / "geometry.yaml")
    particles = _load_yaml(PARAMS_DIR / "particles.yaml")

    bed = geometry["bed_0p076m"]
    v2d = geometry["spouted_bed_V_2d"]
    init = geometry["spouted_bed_V_initial"]
    case_V = particles["case_V_rpt"]
    fb = particles["fb_sweep"]  # reuse dt_init_s convention

    R_c_m = 0.5 * float(bed["D_c_m"])
    R_i_m = 0.5 * float(bed["D_i_m"])
    H_static_m = float(v2d["H_static_m"])
    H_dom_m = float(v2d["H_dom_m"])
    eps_s_bed = float(init["eps_s_bed"])
    # ep_star is a particle property (RCP of monodisperse smooth spheres); for
    # Case V it lives on case_V_rpt and is identical to case_D's value.
    ep_star = float(case_V["ep_star"])
    phi_deg = float(case_V["friction_angle_deg"])

    # Full-width Cartesian domain (same convention as V_3d), kmax=1 slab.
    x_min_m = -R_c_m
    x_max_m = R_c_m
    z_thickness = float(v2d["z_thickness_m"])
    z_min_m = -0.5 * z_thickness
    z_max_m = 0.5 * z_thickness

    params = {
        "description": description or f"V_2d U_in={U_in_m_s:.4f} m/s",
        "run_name": run_name,
        "tstop_s": float(tstop_s),
        "dt_s": float(fb["dt_init_s"]),
        "imax": int(v2d["imax"]),
        "jmax": int(v2d["jmax"]),
        "kmax": int(v2d["kmax"]),
        "x_min_m": x_min_m,
        "x_max_m": x_max_m,
        "z_min_m": z_min_m,
        "z_max_m": z_max_m,
        "zlength_m": z_thickness,
        "R_c_m": R_c_m,
        "R_i_m": R_i_m,
        # Y_CONE quadric parameters.
        "cone_apex_y_m": float(v2d["cone_apex_y_m"]),
        "cone_top_y_m": float(v2d["cone_top_y_m"]),
        "cone_half_angle_deg": float(v2d["cone_half_angle_deg"]),
        # Inlet / spout-core footprint: strip |x| <= R_i across the full slab.
        "bc_mi_x_w_m": -R_i_m,
        "bc_mi_x_e_m": R_i_m,
        "H_static_m": H_static_m,
        "H_dom_m": H_dom_m,
        "H_mid_m": 0.5 * H_static_m,
        "ep_g_bed": 1.0 - eps_s_bed,
        "eps_s_bed": eps_s_bed,
        "ep_star": ep_star,
        "phi_deg": phi_deg,
        "U_in_m_s": float(U_in_m_s),
        "vtk_dt_s": 0.05,
        "drag_type": str(drag_type),
        "c_e": float(c_e),
        # Numerical knobs with safe Day-8 defaults; overridable per case.
        "dt_max_s": 0.01,
        "max_nit": 50,
    }
    params.update(_common_fluid_solid_params("air_25C"))
    return params


def build_V_2d_day9_params(
    run_name: str = "V_2d_day9",
    description: str = "",
) -> dict:
    """Day-9 Case V_2d spouting run parameters.

    Reads every input from params/*.yaml:

      - params/geometry.yaml bed_0p076m: D_c, D_i, (cone geometry)
      - params/particles.yaml case_D: d_p, rho_p
      - params/particles.yaml spouted_bed_V_initial / fb_sweep: H_static, ep_s_bed
      - params/particles.yaml V_2d_day9: u_over_ums, drag_type, c_e, tstop_s
      - params/gases.yaml air_20C via Cantera: rho_g, mu_g

    The orifice velocity imposed in MFiX is the local orifice-face velocity
    that gives the SAME column-superficial velocity in the 2D slab as a 3D
    column operating at 1.2 * U_ms,MG:

        U_ms,col = (d_p/D_c) * (D_i/D_c)^(1/3) * sqrt(2 g H (rho_p - rho)/rho)
        U_in     = u_over_ums * U_ms,col * (D_c / D_i)       [2D slab scaling]

    V_2d is a planar slab, so the orifice is a strip of width 2*R_i and the
    column is a strip of width 2*R_c -- the slab area ratio is D_c/D_i, not
    (D_c/D_i)^2 (that would be the 3D disc ratio).  Using the 3D ratio drives
    bc_v_g up by an extra factor of D_c/D_i = 8, which with Case V_rpt (2.18 mm
    beads, 2400 kg/m^3) gives bc_v_g ~ 70 m/s and causes an immediate solver
    blow-up at the orifice shoulder (|Ug| ~ 7000 m/s, DT < DT_MIN) because the
    gas-particle slip exceeds what Syamlal-O'Brien drag can resolve on a mm-
    scale mesh.  The slab scaling matches the column superficial per unit
    depth, which is the physically meaningful quantity to carry between 2D
    and 3D vessels.  See notes/day09_debug.md Runs 3-4 for the collapsed
    attempt and the switch rationale.
    """
    from src.correlations import ums_mathur_gishler
    from src.gasprops import gas_props

    geometry = _load_yaml(PARAMS_DIR / "geometry.yaml")
    particles = _load_yaml(PARAMS_DIR / "particles.yaml")
    gases = _load_yaml(PARAMS_DIR / "gases.yaml")

    bed = geometry["bed_0p076m"]
    case_V = particles["case_V_rpt"]
    v2d_init = geometry["spouted_bed_V_initial"]
    day9 = particles["V_2d_day9"]
    gas_case = gases["cases"]["air_25C"]
    rho_g, _mu_g = gas_props(
        gas_case["T_K"], gases["pressure_Pa"], gas_case["composition"]
    )

    v2d_geom = geometry["spouted_bed_V_2d"]
    H_static_m = float(v2d_geom["H_static_m"])  # Mathur-Gishler H = static bed

    D_c_m = float(bed["D_c_m"])
    D_i_m = float(bed["D_i_m"])
    U_ms_col_m_s = ums_mathur_gishler(
        d_p=float(case_V["d_p_m"]),
        rho_p=float(case_V["rho_p_kg_per_m3"]),
        rho=float(rho_g),
        D_c=D_c_m,
        D_i=D_i_m,
        H=H_static_m,
    )
    # 2D-slab area ratio (strip widths), not the 3D disc ratio (areas).
    area_ratio = D_c_m / D_i_m
    U_in_m_s = float(day9["u_over_ums"]) * U_ms_col_m_s * area_ratio

    params = build_V_2d_params(
        U_in_m_s=U_in_m_s,
        tstop_s=float(day9["tstop_s"]),
        run_name=run_name,
        description=description
        or (
            f"Day-9 V_2d spouting: U_in={U_in_m_s:.3f} m/s "
            f"({day9['u_over_ums']:.2f}x U_ms,MG={U_ms_col_m_s:.4f} m/s via (Dc/Di)^2)"
        ),
        drag_type=str(day9["drag_type"]),
        c_e=float(day9["c_e"]),
    )
    # Day-9 mesh overrides (notes/day09_debug.md Run 2): coarsen imax/jmax so
    # the orifice stays at the Day-8 >= 4-cell floor while dt_CFL doubles.
    if "imax" in day9:
        params["imax"] = int(day9["imax"])
    if "jmax" in day9:
        params["jmax"] = int(day9["jmax"])
    # Day-9 solver stabilisation (notes/day09_debug.md Run 6):
    #  dt_max_s cap stops MFiX from re-growing dt past the local CFL floor
    #  after Newton recovery; max_nit lifts the inner-iteration ceiling so
    #  stiff startup steps have more head-room before dt_fac shrinks dt.
    if "dt_max_s" in day9:
        params["dt_max_s"] = float(day9["dt_max_s"])
    if "max_nit" in day9:
        params["max_nit"] = int(day9["max_nit"])
    # Day-9 IC override: drop ep_s_bed from the Day-8 settled value to a
    # pre-fluidised value so the TFM solver can start from a dilute state
    # (see notes/day09_debug.md Run 8).  Day-8 ep_s_bed stays untouched.
    if "eps_s_bed" in day9:
        new_eps_s = float(day9["eps_s_bed"])
        params["eps_s_bed"] = new_eps_s
        params["ep_g_bed"] = 1.0 - new_eps_s
    # Record the derived scalars so downstream tooling / tests can introspect
    # them without re-running Cantera.
    params["U_ms_column_m_s"] = float(U_ms_col_m_s)
    params["u_over_ums"] = float(day9["u_over_ums"])
    return params


def build_V_3d_params(
    U_in_m_s: float,
    tstop_s: float,
    run_name: str,
    description: str = "",
) -> dict:
    """Assemble a template parameter dict for the Case V 3D cut-cell spouted bed.

    Pre-condition: cases/V_3d/geometry.stl must exist before running MFiX on
    this case; generate it with src.make_stl.build_spouted_bed_stl using
    bed_0p076m D_c/D_i, cone half-angle 30 deg, H_dom_m from spouted_bed_V_3d.

    U_in_m_s   : superficial gas velocity at the orifice [m/s]
    tstop_s    : run end time [s]
    run_name   : MFiX run_name (also used for VTK/monitor file base names)
    description: free-form description string for the .mfx header
    """
    geometry = _load_yaml(PARAMS_DIR / "geometry.yaml")
    particles = _load_yaml(PARAMS_DIR / "particles.yaml")

    bed = geometry["bed_0p076m"]
    v3d = geometry["spouted_bed_V_3d"]
    init = geometry["spouted_bed_V_initial"]
    case_V = particles["case_V_rpt"]
    fb = particles["fb_sweep"]

    R_c_m = 0.5 * float(bed["D_c_m"])
    R_i_m = 0.5 * float(bed["D_i_m"])
    H_static_m = float(v3d["H_static_m"])
    H_dom_m = float(v3d["H_dom_m"])
    eps_s_bed = float(init["eps_s_bed"])
    # ep_star / phi come from the Case V_rpt particle inventory (R13/R14).
    ep_star = float(case_V["ep_star"])
    phi_deg = float(case_V["friction_angle_deg"])

    # Domain is a cube of side 2*R_c in x and z; the STL carves the vessel.
    x_min_m = -R_c_m
    x_max_m = R_c_m
    z_min_m = -R_c_m
    z_max_m = R_c_m

    params = {
        "description": description or f"V_3d U_in={U_in_m_s:.4f} m/s",
        "run_name": run_name,
        "tstop_s": float(tstop_s),
        "dt_s": float(fb["dt_init_s"]),
        "imax": int(v3d["imax"]),
        "jmax": int(v3d["jmax"]),
        "kmax": int(v3d["kmax"]),
        "x_min_m": x_min_m,
        "x_max_m": x_max_m,
        "z_min_m": z_min_m,
        "z_max_m": z_max_m,
        "H_static_m": H_static_m,
        "H_dom_m": H_dom_m,
        "H_mid_m": 0.5 * H_static_m,
        # Inlet / spout-core footprint: square circumscribing the orifice disc.
        "bc_mi_x_w_m": -R_i_m,
        "bc_mi_x_e_m": R_i_m,
        "bc_mi_z_b_m": -R_i_m,
        "bc_mi_z_t_m": R_i_m,
        "ep_g_bed": 1.0 - eps_s_bed,
        "eps_s_bed": eps_s_bed,
        "ep_star": ep_star,
        "phi_deg": phi_deg,
        "U_in_m_s": float(U_in_m_s),
        "vtk_dt_s": 0.05,
    }
    params.update(_common_fluid_solid_params("air_25C"))
    return params


def _cli() -> None:
    p = argparse.ArgumentParser(description="Render an MFiX .mfx case from a Jinja2 template.")
    p.add_argument(
        "--kind",
        choices=["fb_sweep", "V_2d", "V_3d"],
        default="fb_sweep",
        help="which parameter builder to use (default: fb_sweep)",
    )
    p.add_argument("--template", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--U", type=float, required=True, help="inlet velocity [m/s]")
    p.add_argument("--imax", type=int, default=None, help="(fb_sweep only)")
    p.add_argument("--jmax", type=int, default=None, help="(fb_sweep only)")
    p.add_argument("--tstop", type=float, required=True, help="run end time [s]")
    p.add_argument("--run-name", type=str, required=True)
    p.add_argument("--description", type=str, default="")
    args = p.parse_args()

    if args.kind == "fb_sweep":
        if args.imax is None or args.jmax is None:
            p.error("--imax and --jmax are required for --kind fb_sweep")
        params = build_fb_sweep_params(
            U_in_m_s=args.U,
            imax=args.imax,
            jmax=args.jmax,
            tstop_s=args.tstop,
            run_name=args.run_name,
            description=args.description,
        )
    elif args.kind == "V_2d":
        params = build_V_2d_params(
            U_in_m_s=args.U,
            tstop_s=args.tstop,
            run_name=args.run_name,
            description=args.description,
        )
    elif args.kind == "V_3d":
        params = build_V_3d_params(
            U_in_m_s=args.U,
            tstop_s=args.tstop,
            run_name=args.run_name,
            description=args.description,
        )
    else:
        raise AssertionError(args.kind)

    out = render_case(args.template, params, args.out)
    print(out)


if __name__ == "__main__":
    _cli()

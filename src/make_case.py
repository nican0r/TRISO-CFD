"""Render an MFiX .mfx case file from a Jinja2 template plus a parameter dict.

The template is expected to use only MFiX 26.1.2 keywords that have been
confirmed against the installed namelist reference (CLAUDE.md rule 1 -- no
invented keywords).  All physical quantities in the ``params`` dict are SI;
the template simply substitutes them.

Public API:

    render_case(template_path, params, out_mfx_path) -> Path

    build_fb_sweep_params(U_in, imax, jmax, tstop_s, run_name, description="") -> dict

    build_V_2d_params(U_in_m_s=30.0, tstop_s=3.0, run_name="V_2d", description="") -> dict
"""

from __future__ import annotations

import argparse
import math
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

    Uses StrictUndefined so a missing template variable is a hard error, not
    a silent empty string in a .mfx file.  Returns the output path.
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
    """Assemble a template parameter dict for a single fb_sweep run."""
    from src.gasprops import gas_props

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


def _bed_height_from_inventory(
    mass_kg: float,
    rho_p_kg_m3: float,
    eps_s_bed: float,
    R_i_m: float,
    cone_half_angle_deg: float,
) -> float:
    """Static bed height [m] for a conical bed, derived from inventory.

    The vessel inner profile (above the orifice plane at y=0) is a cone
    r(y) = R_i + y*tan(theta) up to the cone top; this helper assumes the
    bed sits entirely in the cone (confirmed by hand-calc in day-8.md).

    Bed volume in a frustum from y=0 (r=R_i) to y=h (r=r_h):
        V_frustum(h) = (pi*h/3) * (R_i^2 + R_i*r_h + r_h^2)
    with r_h = R_i + h*tan(theta).

    Setting V_frustum(h) = m_p / (rho_p * eps_s_bed) and solving for h
    numerically (bisection on 0 <= h <= 1 m).
    """
    tan_t = math.tan(math.radians(cone_half_angle_deg))
    v_target = mass_kg / (rho_p_kg_m3 * eps_s_bed)

    def v_at(h: float) -> float:
        r_h = R_i_m + h * tan_t
        return (math.pi * h / 3.0) * (R_i_m * R_i_m + R_i_m * r_h + r_h * r_h)

    lo, hi = 0.0, 1.0
    if v_at(hi) < v_target:
        raise ValueError(
            "bed does not fit in a 1 m cone; check inventory / density / eps_s"
        )
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if v_at(mid) < v_target:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-9:
            break
    return 0.5 * (lo + hi)


def _common_fluid_solid_params(gas_case: str = "air_25C") -> dict:
    """Shared particle + gas properties for the Case V template.

    Reads Case V particles (case_V from params/particles.yaml) and the
    chosen gas case from params/gases.yaml, computing rho_g and mu_g via
    Cantera.
    """
    from src.gasprops import gas_props

    particles = _load_yaml(PARAMS_DIR / "particles.yaml")
    gases = _load_yaml(PARAMS_DIR / "gases.yaml")
    part = particles["case_V"]
    gas = gases["cases"][gas_case]
    rho_g, mu_g = gas_props(gas["T_K"], gases["pressure_Pa"], gas["composition"])
    return {
        "d_p_m": float(part["d_p_m"]),
        "rho_p_kg_m3": float(part["rho_p_kg_per_m3"]),
        "rho_g_kg_m3": float(rho_g),
        "mu_g_pa_s": float(mu_g),
    }


def build_V_2d_params(
    U_in_m_s: float = 30.0,
    tstop_s: float = 3.0,
    run_name: str = "V_2d",
    description: str = "",
    drag_type: str = "SYAM_OBRIEN",
) -> dict:
    """Assemble the Case V 2D (planar quadric cut-cell slab) template params.

    All physical values come from params/geometry.yaml (bed_V and
    spouted_bed_V_2d), params/particles.yaml (case_V), and
    params/gases.yaml (air at ambient, S2 operating gas, via Cantera).

    Closures per S1 Table 1:
      - drag_type       = 'SYAM_OBRIEN'   (isothermal Case V default)
      - friction_model  = 'SCHAEFFER'
      - blending_function = 'SIGM_BLEND'  ("modified sigmoidal" blending)
      - kt_type         = 'LUN_1984'      (granular energy as a PDE)
      - phi             = 15 deg          (internal friction angle)
      - ep_star         = 0.42            (void fraction at min fluidization)
      - c_e             = 0.9             (restitution; assumed, S1/S2 silent)

    Geometry: planar Cartesian slab (coordinates='CARTESIAN'), cut-cell
    quadric for the cone wall (n_quadric=1, Y_CONE), no_k=.True. (true 2D).

    U_in_m_s   : inlet jet velocity through the orifice [m/s] (default 30)
    tstop_s    : run end time [s] (default 3.0, Day-9 target)
    run_name   : MFiX run_name (also used for VTK/monitor file base names)
    description: free-form description string for the .mfx header
    drag_type  : MFiX drag_type keyword (default SYAM_OBRIEN per S1 Table 1)
    """
    geometry = _load_yaml(PARAMS_DIR / "geometry.yaml")
    particles = _load_yaml(PARAMS_DIR / "particles.yaml")

    bed = geometry["bed_V"]
    v2d = geometry["spouted_bed_V_2d"]
    case_V = particles["case_V"]

    D_c_m = float(bed["D_c_m"])
    D_i_m = float(bed["D_i_m"])
    R_c_m = 0.5 * D_c_m
    R_i_m = 0.5 * D_i_m
    cone_half_angle_deg = float(v2d["cone_half_angle_deg"])
    H_dom_m = float(v2d["H_dom_m"])

    rho_p_kg_m3 = float(case_V["rho_p_kg_per_m3"])
    inventory_kg = float(case_V["inventory_g"]) * 1.0e-3
    eps_s_bed = float(case_V["eps_s_bed"])
    ep_star = float(case_V["ep_star"])
    phi_deg = float(case_V["friction_angle_deg"])
    c_e = float(case_V["c_e"])

    H_static_m = _bed_height_from_inventory(
        mass_kg=inventory_kg,
        rho_p_kg_m3=rho_p_kg_m3,
        eps_s_bed=eps_s_bed,
        R_i_m=R_i_m,
        cone_half_angle_deg=cone_half_angle_deg,
    )

    # Cone geometry derived scalars.
    tan_t = math.tan(math.radians(cone_half_angle_deg))
    cone_apex_y_m = -R_i_m / tan_t   # y of (virtual) apex, r=0
    cone_top_y_m = (R_c_m - R_i_m) / tan_t  # y where cone meets cylinder

    # Full-width Cartesian domain, kmax=1 slab.
    x_min_m = -R_c_m
    x_max_m = R_c_m
    z_thickness = float(v2d["z_thickness_m"])
    z_min_m = -0.5 * z_thickness
    z_max_m = 0.5 * z_thickness

    params = {
        "description": description or f"V_2d UTK/ORNL U_in={U_in_m_s:.2f} m/s",
        "run_name": run_name,
        "tstop_s": float(tstop_s),
        "dt_s": 1.0e-4,
        "dt_max_s": 1.0e-3,
        "max_nit": 50,
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
        "cone_apex_y_m": cone_apex_y_m,
        "cone_top_y_m": cone_top_y_m,
        "cone_half_angle_deg": cone_half_angle_deg,
        # Inlet strip |x| <= R_i across the full slab.
        "bc_mi_x_w_m": -R_i_m,
        "bc_mi_x_e_m": R_i_m,
        "H_static_m": H_static_m,
        "H_dom_m": H_dom_m,
        "H_mid_m": 0.5 * H_static_m,
        "H0_over_Dc": H_static_m / D_c_m,
        "ep_g_bed": 1.0 - eps_s_bed,
        "eps_s_bed": eps_s_bed,
        "ep_star": ep_star,
        "phi_deg": phi_deg,
        "c_e": c_e,
        "U_in_m_s": float(U_in_m_s),
        "vtk_dt_s": 0.05,
        "monitor_dt_s": 1.0e-3,     # 1 kHz sampling (matches S1 and S2).
        "drag_type": str(drag_type),
        "friction_model": "SCHAEFFER",
        # MFiX 26.1.2's accepted token for the scaled-sigmoidal blending
        # (= S1 Table 1's "modified sigmoidal") is 'GIDASPOW_PCF' -- the
        # namelist header lists 'SIGM_BLEND' as a valid value but
        # check_blending_function in check_solids_continuum.f rejects it
        # and only accepts 'GIDASPOW_PCF' to activate the SIGM_BLEND flag.
        "blending_function": "GIDASPOW_PCF",
        "kt_type": "LUN_1984",               # PDE granular energy per S1 Tbl 1
    }
    params.update(_common_fluid_solid_params("air_25C"))
    return params


def _build_graded_axis(segments: list[dict]) -> list[float]:
    """Flatten a list of {n_cells, d*_m} segments into a per-cell size list.

    segments : [{n_cells: int, dx_m: float}, ...] (or dy_m / dz_m key)
    Returns the per-cell size list [m] of length sum(n_cells).
    """
    out: list[float] = []
    for seg in segments:
        n = int(seg["n_cells"])
        # Accept dx_m, dy_m or dz_m.
        key = next(k for k in seg if k != "n_cells")
        d = float(seg[key])
        out.extend([d] * n)
    return out


def build_V_3d_params(
    U_in_m_s: float = 30.0,
    tstop_s: float = 5.5,
    run_name: str = "V_3d",
    description: str = "",
    drag_type: str = "SYAM_OBRIEN",
    stl_filename: str = "geometry_0001.stl",
) -> dict:
    """Assemble the Case V 3D (STL cut-cell, MP-PIC) template params.

    Day 10 PIC build.  Uses MFiX's Particle-in-cell solids_model with the
    Harris-Crighton continuum stress closure from NETL's canonical
    tutorials/pic/spouted_bed_3d template.  Gas-side closure stays
    SYAM_OBRIEN drag (same as Case V_2d / S1 Table 1); PIC-specific
    coefficients (fric_exp_pic, fric_non_sing_fac, mppic_coeff_en1,
    mppic_coeff_en_wall, mppic_coeff_et_wall, mppic_velfac_coeff,
    psfac_fric_pic) come from params/geometry.yaml spouted_bed_V_3d_pic.
    des_epg_clip is set to 1 - ep_star (= eps_s_max for Case V) so the
    PIC solids-stress clip matches the TFM minimum-fluidization value.

    Why PIC on this geometry.  MFiX 26.1.2's STL cut-cell F_AT classifier
    (cartesian_grid/intersect.f) mis-classifies the near-axis column of
    an open-top surface-of-revolution STL -- a bug that cannot be fixed
    from the input deck.  PIC parcels track in continuous space so the
    solids-phase advection is not dependent on F_AT; PIC also avoids the
    LUN_1984 PDE stiffness that killed the quadric-TFM attempt (see
    notes/day10_summary.md for the full diagnosis).

    Geometry comes from params/geometry.yaml bed_V (R_c = 25 mm,
    R_i = 2 mm, cone half-angle 30 deg) and the Day-10 spouted_bed_V_3d
    block (graded x,z and y meshes, L_stub = 10 mm inlet stub tube).

    U_in_m_s     : physical inlet jet velocity through the 4 mm orifice [m/s]
                   (bc_v_g on the MI slab box |x|<=R_i, |z|<=R_i at
                   y = -L_stub).
    tstop_s      : run end time [s]
    run_name     : MFiX run_name (also used for VTK and monitor file names)
    description  : free-form description string for the .mfx header
    drag_type    : MFiX drag_type keyword (default SYAM_OBRIEN per S1 Table 1)
    stl_filename : expected STL file name (MFiX reads 'geometry.stl' by
                   default; the Day-10 driver writes the STL to the run dir
                   under this name)

    Returns a dict with every template variable, including the per-cell
    DX/DY/DZ lists (``dx_m_list``, ``dy_m_list``, ``dz_m_list``).  SI units
    everywhere (CLAUDE.md rule 3).
    """
    geometry = _load_yaml(PARAMS_DIR / "geometry.yaml")
    particles = _load_yaml(PARAMS_DIR / "particles.yaml")

    bed = geometry["bed_V"]
    v3d = geometry["spouted_bed_V_3d"]
    pic = geometry["spouted_bed_V_3d_pic"]
    case_V = particles["case_V"]

    D_c_m = float(bed["D_c_m"])
    D_i_m = float(bed["D_i_m"])
    cone_half_angle_deg = 0.5 * float(bed["cone_angle_deg"])
    R_c_m = 0.5 * D_c_m
    R_i_m = 0.5 * D_i_m
    H_dom_m = float(v3d["H_dom_m"])
    theta_segments = int(v3d["theta_segments"])
    L_stub_m = float(v3d.get("L_stub_m", 0.0))
    dx_sample_m = float(v3d.get("dx_sample_m", 5.0e-4))
    H_total_m = H_dom_m + L_stub_m

    dx_m_list = _build_graded_axis(v3d["xz_segments"])
    dy_m_list = _build_graded_axis(v3d["y_segments"])
    dz_m_list = _build_graded_axis(v3d["xz_segments"])
    dy_min = min(dy_m_list)
    imax = len(dx_m_list)
    jmax = len(dy_m_list)
    kmax = len(dz_m_list)

    # Sanity: xz-mesh centred on 0, sums must equal domain extents.
    # y-mesh now sums to H_total = L_stub + H_dom (stub tube below y=0).
    sum_x = sum(dx_m_list)
    sum_z = sum(dz_m_list)
    sum_y = sum(dy_m_list)
    if abs(sum_x - D_c_m) > 1e-9 or abs(sum_z - D_c_m) > 1e-9:
        raise ValueError(
            f"xz_segments sum ({sum_x} m) must equal D_c_m ({D_c_m} m)."
        )
    if abs(sum_y - H_total_m) > 1e-9:
        raise ValueError(
            f"y_segments sum ({sum_y} m) must equal H_total "
            f"(L_stub + H_dom = {H_total_m} m)."
        )

    rho_p_kg_m3 = float(case_V["rho_p_kg_per_m3"])
    inventory_kg = float(case_V["inventory_g"]) * 1.0e-3
    eps_s_bed = float(case_V["eps_s_bed"])
    ep_star = float(case_V["ep_star"])
    # PIC does not use phi (internal friction angle) or c_e (binary
    # restitution) directly; the Harris-Crighton continuum stress uses
    # fric_exp_pic / fric_non_sing_fac / psfac_fric_pic for the solid
    # pressure branch and mppic_coeff_en1 / mppic_coeff_en_wall /
    # mppic_coeff_et_wall for the parcel damping branch (see the PIC
    # block below).

    H_static_m = _bed_height_from_inventory(
        mass_kg=inventory_kg,
        rho_p_kg_m3=rho_p_kg_m3,
        eps_s_bed=eps_s_bed,
        R_i_m=R_i_m,
        cone_half_angle_deg=cone_half_angle_deg,
    )

    x_min_m = -R_c_m
    x_max_m = R_c_m
    z_min_m = -R_c_m
    z_max_m = R_c_m
    # Day-10 fix 2: Cartesian box extends from y = -L_stub to y = H_dom so
    # the MI BC sits at the physical box bottom (MFiX error 1100 avoided)
    # and the near-axis column of cells inside the stub has STL sidewall
    # facets within 1-2 cells (seeds F_AT propagation along J).
    y_min_m = -L_stub_m
    y_max_m = H_dom_m
    bc_mi_y_s_m = -L_stub_m
    bc_mi_y_n_m = -L_stub_m

    # Quadric cone geometry: virtual apex at (0, -R_i/tan(half_angle), 0)
    # so the cone has r = R_i at y = 0 (orifice rim) and r = R_c at
    # y = cone_top (where the cylinder takes over).
    tan_t = math.tan(math.radians(cone_half_angle_deg))
    cone_apex_y_m = -R_i_m / tan_t
    cone_top_y_m = (R_c_m - R_i_m) / tan_t

    params = {
        "description": description or f"V_3d UTK/ORNL U_in={U_in_m_s:.2f} m/s",
        "run_name": run_name,
        "tstop_s": float(tstop_s),
        "dt_s": 1.0e-4,
        "dt_max_s": 1.0e-3,
        "max_nit": 50,
        "imax": int(imax),
        "jmax": int(jmax),
        "kmax": int(kmax),
        "x_min_m": x_min_m,
        "x_max_m": x_max_m,
        "y_min_m": y_min_m,
        "y_max_m": y_max_m,
        "z_min_m": z_min_m,
        "z_max_m": z_max_m,
        "R_c_m": R_c_m,
        "R_i_m": R_i_m,
        "L_stub_m": L_stub_m,
        "cone_half_angle_deg": cone_half_angle_deg,
        "cone_apex_y_m": cone_apex_y_m,
        "cone_top_y_m": cone_top_y_m,
        "theta_segments": theta_segments,
        "dx_sample_m": dx_sample_m,
        "bc_mi_x_w_m": -R_i_m,
        "bc_mi_x_e_m": R_i_m,
        "bc_mi_z_b_m": -R_i_m,
        "bc_mi_z_t_m": R_i_m,
        "bc_mi_y_s_m": bc_mi_y_s_m,
        "bc_mi_y_n_m": bc_mi_y_n_m,
        "H_static_m": H_static_m,
        "H_dom_m": H_dom_m,
        "H_mid_m": 0.5 * H_static_m,
        "H0_over_Dc": H_static_m / D_c_m,
        "dy_min_m": float(dy_min),
        "ep_g_bed": 1.0 - eps_s_bed,
        "eps_s_bed": eps_s_bed,
        "ep_star": ep_star,
        "U_in_m_s": float(U_in_m_s),
        "vtk_dt_s": 0.05,
        "monitor_dt_s": 1.0e-3,
        "drag_type": str(drag_type),
        # PIC closures: Harris-Crighton continuum solids stress plus the
        # shared SCHAEFFER friction label used by MFiX's PIC branch.  All
        # numeric values pulled from params/geometry.yaml
        # (spouted_bed_V_3d_pic block), which cites the NETL
        # tutorials/pic/spouted_bed_3d template.
        "friction_model": str(pic.get("friction_model", "SCHAEFFER")),
        "fric_exp_pic": float(pic["fric_exp_pic"]),
        "fric_non_sing_fac": float(pic["fric_non_sing_fac"]),
        "psfac_fric_pic": float(pic["psfac_fric_pic"]),
        "mppic_coeff_en1": float(pic["mppic_coeff_en1"]),
        "mppic_coeff_en_wall": float(pic["mppic_coeff_en_wall"]),
        "mppic_coeff_et_wall": float(pic["mppic_coeff_et_wall"]),
        "mppic_velfac_coeff": float(pic["mppic_velfac_coeff"]),
        "des_interp_scheme": str(pic["des_interp_scheme"]),
        "des_epg_clip": float(pic["des_epg_clip"]),
        "ic_pic_const_statwt": float(pic["ic_pic_const_statwt"]),
        "bc_pic_mi_const_statwt": float(pic["bc_pic_mi_const_statwt"]),
        # Per-cell size arrays rendered into dx(0)..dx(IMAX-1).
        "dx_m_list": dx_m_list,
        "dy_m_list": dy_m_list,
        "dz_m_list": dz_m_list,
        "stl_filename": str(stl_filename),
    }
    params.update(_common_fluid_solid_params("air_25C"))
    return params


def _cli() -> None:
    p = argparse.ArgumentParser(description="Render an MFiX .mfx case from a Jinja2 template.")
    p.add_argument(
        "--kind",
        choices=["fb_sweep", "V_2d", "V_3d"],
        default="fb_sweep",
    )
    p.add_argument("--template", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--U", type=float, default=30.0, help="inlet velocity [m/s]")
    p.add_argument("--imax", type=int, default=None, help="(fb_sweep only)")
    p.add_argument("--jmax", type=int, default=None, help="(fb_sweep only)")
    p.add_argument("--tstop", type=float, default=3.0, help="run end time [s]")
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

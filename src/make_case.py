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
    R_c_m: float | None = None,
) -> float:
    """Static bed height [m] for a conical bed, derived from inventory.

    The vessel inner profile (above the orifice plane at y=0) is a cone
    r(y) = R_i + y*tan(theta) up to the cone top; with R_c_m=None this helper
    assumes the bed sits entirely in the cone (confirmed by hand-calc in
    day-8.md).  When R_c_m is given, the profile continues as a cylinder of
    radius R_c above the cone top (y = (R_c - R_i)/tan(theta)); used for the
    loose PIC IC bed that extends above the cone (Day-10 R5).

    Bed volume in a frustum from y=0 (r=R_i) to y=h (r=r_h):
        V_frustum(h) = (pi*h/3) * (R_i^2 + R_i*r_h + r_h^2)
    with r_h = R_i + h*tan(theta).

    Setting V_frustum(h) = m_p / (rho_p * eps_s_bed) and solving for h
    numerically (bisection on 0 <= h <= 1 m).
    """
    tan_t = math.tan(math.radians(cone_half_angle_deg))
    v_target = mass_kg / (rho_p_kg_m3 * eps_s_bed)

    h_cone = None if R_c_m is None else (R_c_m - R_i_m) / tan_t

    def v_frustum(h: float) -> float:
        r_h = R_i_m + h * tan_t
        return (math.pi * h / 3.0) * (R_i_m * R_i_m + R_i_m * r_h + r_h * r_h)

    def v_at(h: float) -> float:
        if h_cone is None or h <= h_cone:
            return v_frustum(h)
        return v_frustum(h_cone) + math.pi * R_c_m * R_c_m * (h - h_cone)

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


def _uniform_axis(dx_m: float, length_m: float) -> list[float]:
    """Return a flat per-cell size list [m] for a uniform-mesh axis.

    dx_m     : uniform cell size [m]
    length_m : total axis extent [m]

    The number of cells is round(length_m / dx_m); the function raises
    ValueError if the implied extent differs from length_m by more than
    one tenth of a cell (= 10 % of dx_m), which would indicate a
    bad user input (dx_m that does not divide length_m).
    """
    if dx_m <= 0.0 or length_m <= 0.0:
        raise ValueError("dx_m and length_m must both be positive")
    n = int(round(length_m / dx_m))
    if n < 1:
        raise ValueError(f"length_m {length_m} / dx_m {dx_m} rounds to 0 cells")
    implied = n * dx_m
    if abs(implied - length_m) > 0.1 * dx_m:
        raise ValueError(
            f"dx_m {dx_m} does not divide length_m {length_m} (implied "
            f"{implied} vs requested {length_m})"
        )
    return [dx_m] * n


def build_V_3d_params(
    U_in_m_s: float = 30.0,
    tstop_s: float = 5.5,
    run_name: str = "V_3d",
    description: str = "",
    drag_type: str = "SYAM_OBRIEN",
    stl_filename: str = "geometry_0001.stl",
    gas_case: str = "air_25C",
    nodes: tuple[int, int, int] = (1, 1, 1),
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
    R_i = 2 mm, cone half-angle 30 deg) and the Day-10 R1 spouted_bed_V_3d
    block (uniform 2 mm Cartesian mesh, L_stub = 20 mm stub cylinder).

    Day-10 R1 design (NETL-clone, see notes/day10_summary.md "R1 — NETL-
    clone redesign"):
      - uniform 2 mm cells in x, y, z (no grading);
      - Cartesian box in x,z: [-30 mm, +30 mm] (5 mm margin outside R_c);
      - Cartesian box in y: [-L_stub, H_dom] = [-20 mm, 200 mm];
      - MI applied across the FULL box bottom face (y = -L_stub); the
        STL stub-cylinder cross-section carves out the orifice-equivalent
        inflow footprint (circle of radius R_i).

    U_in_m_s     : physical inlet jet velocity through the 4 mm orifice [m/s]
                   (bc_v_g on the full box bottom face at y = -L_stub).
    tstop_s      : run end time [s]
    run_name     : MFiX run_name (also used for VTK and monitor file names)
    description  : free-form description string for the .mfx header
    drag_type    : MFiX drag_type keyword (default SYAM_OBRIEN per S1 Table 1)
    stl_filename : expected STL file name (MFiX reads 'geometry.stl' by
                   default; the Day-10 driver writes the STL to the run dir
                   under this name)
    gas_case     : key in params/gases.yaml selecting the gas composition
                   and temperature.  Default 'air_25C' preserves the Day-10
                   baseline; the Day-13 sweep varies this across air_25C,
                   ar_20C, ar_1300C and ar_h2_1300C.

    Returns a dict with every template variable, including the per-cell
    DX/DY/DZ lists (``dx_m_list``, ``dy_m_list``, ``dz_m_list``).  SI units
    everywhere (CLAUDE.md rule 3).
    """
    geometry = _load_yaml(PARAMS_DIR / "geometry.yaml")
    particles = _load_yaml(PARAMS_DIR / "particles.yaml")
    numerics = _load_yaml(PARAMS_DIR / "numerics.yaml")["spouted_bed_V_3d_numerics"]

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
    L_stub_m = float(v3d.get("L_stub_m", 0.0))
    dx_sample_m = float(v3d.get("dx_sample_m", 1.0e-3))
    dx_uniform_m = float(v3d["dx_uniform_m"])
    # Day-10 R2 Option 4: stub-tube radius decoupled from the physical
    # orifice radius so the MI face has multiple fluid cells.  Fall back to
    # R_i for backward compatibility with R1 decks.
    R_stub_m = float(v3d.get("R_stub_m", R_i_m))
    if R_stub_m < R_i_m:
        raise ValueError(
            f"spouted_bed_V_3d.R_stub_m ({R_stub_m}) must be >= R_i "
            f"({R_i_m}); a stub narrower than the orifice is non-physical."
        )
    # Day-10 R2 Option 2: linear-taper length from R_stub down to R_i.
    # 0 recovers the hard-shelf build; the taper eliminates the shelf-
    # induced velocity shock observed on AWS (see geometry.yaml comment).
    taper_length_m = float(v3d.get("taper_length_m", 0.0))
    if taper_length_m < 0.0:
        raise ValueError(
            f"spouted_bed_V_3d.taper_length_m ({taper_length_m}) must be >= 0."
        )
    if taper_length_m > L_stub_m:
        raise ValueError(
            f"spouted_bed_V_3d.taper_length_m ({taper_length_m}) must be "
            f"<= L_stub_m ({L_stub_m}); the taper must fit inside the stub."
        )
    H_total_m = H_dom_m + L_stub_m

    # Day-10 R1: uniform mesh in x, y, z.  Box in x,z spans ±(R_c + 5 mm
    # margin) = ±30 mm so there is clearance around the vessel wall for
    # the cut-cell preprocessor; this is 6 mm larger than D_c (2 x 30 vs
    # 2 x 25 mm).  Box in y spans [-L_stub, H_dom] = [-20, 200] mm.
    x_half_m = R_c_m + 5.0e-3           # 25 mm + 5 mm margin = 30 mm half-width
    sum_xz_m = 2.0 * x_half_m           # 60 mm box width in x and z
    dx_m_list = _uniform_axis(dx_uniform_m, sum_xz_m)
    dy_m_list = _uniform_axis(dx_uniform_m, H_total_m)
    dz_m_list = _uniform_axis(dx_uniform_m, sum_xz_m)
    dx_min_m = dx_uniform_m
    imax = len(dx_m_list)
    jmax = len(dy_m_list)
    kmax = len(dz_m_list)

    # Sanity: xz-mesh symmetric, sum >= D_c; y-mesh sums to L_stub + H_dom.
    sum_x = sum(dx_m_list)
    sum_z = sum(dz_m_list)
    sum_y = sum(dy_m_list)
    if abs(sum_x - sum_z) > 1e-9:
        raise ValueError(
            f"dx sum ({sum_x} m) must equal dz sum ({sum_z} m)."
        )
    if sum_x < D_c_m:
        raise ValueError(
            f"dx sum ({sum_x} m) must enclose the vessel D_c = {D_c_m} m."
        )
    if abs(sum_y - H_total_m) > 1e-9:
        raise ValueError(
            f"dy sum ({sum_y} m) must equal H_total "
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
    # Day-10 R5: PIC bed quadrants are seeded LOOSE (eps_s_ic_pic, NETL
    # tutorial value) up to the inventory-conserving height H_ic_m, then
    # settle onto the jet.  H_static_m (packed, eps_s_bed) is kept for the
    # monitors and derived scalars.
    eps_s_ic = float(pic["eps_s_ic_pic"])
    H_ic_m = _bed_height_from_inventory(
        mass_kg=inventory_kg,
        rho_p_kg_m3=rho_p_kg_m3,
        eps_s_bed=eps_s_ic,
        R_i_m=R_i_m,
        cone_half_angle_deg=cone_half_angle_deg,
        R_c_m=R_c_m,
    )

    x_min_m = -x_half_m
    x_max_m = x_half_m
    z_min_m = -x_half_m
    z_max_m = x_half_m
    # Day-10 R1 NETL pattern: Cartesian box extends from y = -L_stub to
    # y = H_dom; the MI BC sits at the box bottom (y = -L_stub) across
    # the FULL box face.  The STL (stub cylinder of radius R_i) carves
    # the MI down to the circular orifice-equivalent footprint; the gas
    # cannot enter any cell outside the stub cylinder because the STL
    # marks them as blocked.
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

    # DMP decomposition validation (Day-10 R4 MPI support).  nodes[i]*nodes[j]
    # *nodes[k] = total MPI rank count that mpirun must receive; each
    # component must divide its cell-count axis cleanly so MFiX's domain
    # decomposition does not fault during preprocessing.  AWS Batch overrides
    # these at container launch via docker/entrypoint.sh (idempotent sed),
    # so this is only enforced for the local-render path.
    ni, nj, nk = (int(nodes[0]), int(nodes[1]), int(nodes[2]))
    if ni < 1 or nj < 1 or nk < 1:
        raise ValueError(
            f"nodes must be positive integers; got ({ni}, {nj}, {nk})"
        )
    if imax % ni != 0 or jmax % nj != 0 or kmax % nk != 0:
        raise ValueError(
            f"nodes=({ni},{nj},{nk}) does not divide mesh "
            f"(imax={imax}, jmax={jmax}, kmax={kmax}); each axis must divide "
            f"cleanly or MFiX's DMP decomposition will fault at preprocessing."
        )

    params = {
        "description": description or f"V_3d UTK/ORNL U_in={U_in_m_s:.2f} m/s",
        "run_name": run_name,
        "tstop_s": float(tstop_s),
        "dt_s": float(numerics["dt_s"]),
        "dt_max_s": float(numerics["dt_max_s"]),
        "max_nit": 50,
        "imax": int(imax),
        "jmax": int(jmax),
        "kmax": int(kmax),
        "nodes_i": ni,
        "nodes_j": nj,
        "nodes_k": nk,
        "nranks": ni * nj * nk,
        "x_min_m": x_min_m,
        "x_max_m": x_max_m,
        "y_min_m": y_min_m,
        "y_max_m": y_max_m,
        "z_min_m": z_min_m,
        "z_max_m": z_max_m,
        "R_c_m": R_c_m,
        "R_i_m": R_i_m,
        "R_stub_m": R_stub_m,
        "L_stub_m": L_stub_m,
        "taper_length_m": taper_length_m,
        "cone_half_angle_deg": cone_half_angle_deg,
        "cone_apex_y_m": cone_apex_y_m,
        "cone_top_y_m": cone_top_y_m,
        "dx_sample_m": dx_sample_m,
        "dx_uniform_m": dx_uniform_m,
        # Day-10 R1: MI spans the FULL box bottom face (x_min..x_max,
        # z_min..z_max).  The STL stub cylinder restricts the actual gas
        # inflow to the circular orifice-equivalent footprint.
        "bc_mi_x_w_m": x_min_m,
        "bc_mi_x_e_m": x_max_m,
        "bc_mi_z_b_m": z_min_m,
        "bc_mi_z_t_m": z_max_m,
        "bc_mi_y_s_m": bc_mi_y_s_m,
        "bc_mi_y_n_m": bc_mi_y_n_m,
        # Orifice-aligned square bounds (|x|,|z| <= R_i).  Used by:
        #  - Monitor 1 (inlet-plane P_g probe).
        #  - Monitor 4 (centerline V_g probe at mid-bed).
        # These bounds are NOT applied to the MI BC in Day-10 R1; the MI
        # covers the full box bottom face (see bc_mi_x_w_m etc.).
        "orifice_x_w_m": -R_i_m,
        "orifice_x_e_m": R_i_m,
        "orifice_z_b_m": -R_i_m,
        "orifice_z_t_m": R_i_m,
        # Day-10 R3 Fix A: central gas column half-width used to carve the
        # packed-bed IC into NETL's 4-strip quadrant pattern.  Set to
        # max(R_stub_m, R_i_m + dx_uniform_m) so the gas column at least
        # encloses the stub-tube cross-section AND leaves a 1-cell margin
        # around the orifice (one cell wider than LINEAR_HAT's smearing
        # support, so parcel mass is never projected into orifice-adjacent
        # cells).  Matches the ratio-1 construction NETL uses
        # (column_half_width = R_spout) in tutorials/pic/spouted_bed_3d.
        # R4 experiment (R_i + 2*dx = 6 mm) was reverted because the 25%
        # parcel-count shortfall vs Case V inventory (54.5 g) was a worse
        # paper-fidelity hit than the residual 2-3 stub cells it tried
        # to fix.  See notes/day10_pathology_review.md.
        "R_bed_gas_col_m": max(R_stub_m, R_i_m + dx_uniform_m),
        "H_static_m": H_static_m,
        "H_dom_m": H_dom_m,
        "H_mid_m": 0.5 * H_static_m,
        "H0_over_Dc": H_static_m / D_c_m,
        "dy_min_m": float(dx_min_m),
        "ep_g_bed": 1.0 - eps_s_bed,
        "eps_s_bed": eps_s_bed,
        "H_ic_m": H_ic_m,
        "eps_s_ic": eps_s_ic,
        "ep_g_ic": 1.0 - eps_s_ic,
        "ep_star": ep_star,
        "U_in_m_s": float(U_in_m_s),
        # Day-10 R4 (2026-10-11): straight-stub geometry (R_stub = R_i), so the
        # MI face and the orifice have the same cross-section and bc_v_g_mi
        # equals the paper's physical jet velocity U_in_m_s directly -- no
        # mass-flux rescaling.  The (R_i/R_stub)^2 scaling that lived here
        # under Day-10 R2 Option 4 relied on 4x acceleration through a 1-cell
        # taper that the SIMPLE pressure correction could not deliver on a
        # 2 mm cut-cell mesh (notes/day10_summary.md R4).  The refined 1 mm
        # mesh (geometry.yaml spouted_bed_V_3d.dx_uniform_m) resolves R_i = 2
        # mm with 2 cells across the radius -- enough for F_AT to class the
        # stub as fluid without the Option-4 widening.
        "bc_v_g_mi_m_s": float(U_in_m_s) * (R_i_m / R_stub_m) ** 2,
        "stub_to_orifice_area_ratio": (R_stub_m / R_i_m) ** 2,
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
        # Day-10 R5c: IC layout switch (see geometry.yaml
        # spouted_bed_V_3d_pic.ic_central_gas_column).
        "ic_central_gas_column": bool(pic["ic_central_gas_column"]),
    }
    params.update(_common_fluid_solid_params(gas_case))

    # Day-10 R5c settle-then-ramp MI start-up (geometry.yaml startup_V_3d).
    # v0 is tied to the single-particle terminal velocity of the case's
    # particles in the case's gas so parcels cannot fall through the
    # gas-only stub while the bed settles over the orifice.
    from src.correlations import terminal_velocity

    startup = geometry["startup_V_3d"]
    u_t = terminal_velocity(
        params["d_p_m"], params["rho_p_kg_m3"],
        params["rho_g_kg_m3"], params["mu_g_pa_s"],
    )
    t0 = float(startup["settle_time_s"])
    params.update({
        "ramp_enabled": bool(startup["ramp_enabled"]),
        "u_t_m_s": float(u_t),
        "bc_v_g_ramp_t0_s": t0,
        "bc_v_g_ramp_t1_s": t0 + float(startup["ramp_duration_s"]),
        "bc_v_g_ramp_v0_m_s": float(startup["v0_over_ut"]) * float(u_t),
        "bc_v_g_ramp_v1_m_s": params["bc_v_g_mi_m_s"],
    })
    return params


def build_V_3d_dem_params(
    U_in_m_s: float = 30.0,
    tstop_s: float = 5.5,
    run_name: str = "V_3d_dem",
    description: str = "",
    drag_type: str = "SYAM_OBRIEN",
    stl_filename: str = "geometry_0001.stl",
    cgp_factor: float = 10.0,
    nodes: tuple[int, int, int] = (1, 1, 1),
) -> dict:
    """Assemble the Case V 3D (STL cut-cell, CFD-DEM) template params.

    Day 11-12 build.  Reuses the Day-10 R1 PIC geometry, mesh, STL and
    boundary layout verbatim (per implementation-steps/day-11-12.md: "Reuse
    the Day-10 PIC STL vessel geometry and driver scaffolding ... Switch
    solids_model(1) = 'PIC' -> 'DEM'").  Swaps the PIC Harris-Crighton
    continuum stress for a spring-dashpot soft-sphere DEM contact model
    with values from params/geometry.yaml spouted_bed_V_3d_dem.

    ``cgp_factor`` controls the coarsening:
      - ``cgp_factor == 1``: pure DEM, solids_model(1) = 'DEM', one parcel
        per real particle.  At Case V's 1.376e5 real particles and 500 um
        d_p this is expensive (Rayleigh dt_c ~ 6 us); use for short runs
        or on a cluster.
      - ``cgp_factor > 1``: coarse-grained DEM (CGP), solids_model(1) =
        'CGP', cgp_stat_wt(1) = cgp_factor.  Each parcel represents
        ``cgp_factor`` real particles; the parcel diameter and mass scale
        with cgp_factor^(1/3) and cgp_factor respectively (MFiX internal).
        Day-11-12 default of 10 keeps parcel count ~1.4e4 (same order as
        Day-10 PIC) and raises the Rayleigh dt by sqrt(cgp_factor); see
        notes/day11_summary.md for the full justification.

    U_in_m_s     : physical inlet jet velocity through the 4 mm orifice [m/s]
    tstop_s      : run end time [s]
    run_name     : MFiX run_name (also used for VTK and monitor file names)
    description  : free-form description string for the .mfx header
    drag_type    : MFiX drag_type keyword (default SYAM_OBRIEN per S1 Table 1)
    stl_filename : expected STL file name (same convention as Day-10)
    cgp_factor   : coarse-graining statistical weight [-]; 1 = pure DEM.

    Returns a dict with every template variable, including the per-cell
    DX/DY/DZ lists and all DEM/CGP closure values.  SI units everywhere
    (CLAUDE.md rule 3).
    """
    # Start from the Day-10 PIC params so geometry/mesh/BC/IC stay bit-
    # identical.  Then strip PIC-only keys and add DEM/CGP keys.
    params = build_V_3d_params(
        U_in_m_s=U_in_m_s,
        tstop_s=tstop_s,
        run_name=run_name,
        description=description or f"V_3d_dem UTK/ORNL U_in={U_in_m_s:.2f} m/s",
        drag_type=drag_type,
        stl_filename=stl_filename,
        nodes=nodes,
    )

    geometry = _load_yaml(PARAMS_DIR / "geometry.yaml")
    dem = geometry["spouted_bed_V_3d_dem"]
    particles = _load_yaml(PARAMS_DIR / "particles.yaml")
    case_V = particles["case_V"]

    # Decide solids model from the coarsening factor.  cgp_factor = 1 is
    # pure DEM; > 1 is coarse-grained DEM (CGP) with cgp_stat_wt(1) set.
    cgp_factor = float(cgp_factor)
    if cgp_factor < 1.0:
        raise ValueError(f"cgp_factor must be >= 1 (got {cgp_factor})")
    if cgp_factor == 1.0:
        solids_model = "DEM"
        cgp_stat_wt = 1.0
    else:
        solids_model = "CGP"
        cgp_stat_wt = cgp_factor

    # Real-particle count for provenance (not an MFiX input; just logged in
    # the manifest and in notes/day11_summary.md).
    rho_p = float(case_V["rho_p_kg_per_m3"])
    d_p = float(case_V["d_p_m"])
    inventory_kg = float(case_V["inventory_g"]) * 1.0e-3
    v_particle_m3 = (math.pi / 6.0) * d_p ** 3
    n_real = inventory_kg / (rho_p * v_particle_m3)
    n_parcels_est = n_real / cgp_stat_wt

    # Strip PIC-only keys so the template never renders them for DEM.
    pic_only_keys = (
        "friction_model",          # PIC uses 'SCHAEFFER' label on its
                                    # continuum-stress branch; DEM template
                                    # does not reference it (see template).
        "fric_exp_pic",
        "fric_non_sing_fac",
        "psfac_fric_pic",
        "mppic_coeff_en1",
        "mppic_coeff_en_wall",
        "mppic_coeff_et_wall",
        "mppic_velfac_coeff",
        "ic_pic_const_statwt",
        "bc_pic_mi_const_statwt",
    )
    for key in pic_only_keys:
        params.pop(key, None)

    # Add DEM closure fields.
    params.update({
        "solids_model": solids_model,          # 'DEM' or 'CGP'
        "cgp_factor": float(cgp_factor),       # for the manifest
        "cgp_stat_wt": float(cgp_stat_wt),     # 1.0 if DEM
        "kn": float(dem["kn"]),
        "kn_w": float(dem["kn_w"]),
        "kt_fac": float(dem["kt_fac"]),
        "kt_w_fac": float(dem["kt_w_fac"]),
        "mew": float(dem["mew"]),
        "mew_w": float(dem["mew_w"]),
        "des_en_input_pp": float(dem["des_en_input_pp"]),
        "des_en_input_pw": float(dem["des_en_input_pw"]),
        "des_etat_fac": float(dem["des_etat_fac"]),
        "des_etat_w_fac": float(dem["des_etat_w_fac"]),
        "des_interp_scheme": str(dem["des_interp_scheme"]),
        "des_epg_clip": float(dem["des_epg_clip"]),
        "dem_dt_s": float(dem["dem_dt_s"]),
        "n_real_particles": float(n_real),
        "n_parcels_est": float(n_parcels_est),
    })
    return params


def _cli() -> None:
    p = argparse.ArgumentParser(description="Render an MFiX .mfx case from a Jinja2 template.")
    p.add_argument(
        "--kind",
        choices=["fb_sweep", "V_2d", "V_3d", "V_3d_dem"],
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
    elif args.kind == "V_3d_dem":
        params = build_V_3d_dem_params(
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

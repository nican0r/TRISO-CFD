"""Render the Day-8 Case V geometry schematic to results/fig_d8_V_2d_geometry.png.

Reads geometry + Case V inventory from params/*.yaml, derives the static
bed height via src.make_case._bed_height_from_inventory, and draws the
2D planar cut-cell slab with IC regions, BCs, and the four monitor
locations annotated.
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.patches as mp  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import yaml  # noqa: E402

from src.make_case import _bed_height_from_inventory

ROOT = Path(__file__).resolve().parent.parent


def _load() -> tuple[dict, dict, dict]:
    geom = yaml.safe_load((ROOT / "params" / "geometry.yaml").read_text())
    particles = yaml.safe_load((ROOT / "params" / "particles.yaml").read_text())
    return geom["bed_V"], geom["spouted_bed_V_2d"], particles["case_V"]


def render_V_2d_schematic(out_path: Path) -> Path:
    bed, v2d, case_V = _load()
    R_c = 0.5 * float(bed["D_c_m"])
    R_i = 0.5 * float(bed["D_i_m"])
    half_angle = float(v2d["cone_half_angle_deg"])
    H_dom = float(v2d["H_dom_m"])
    tan_t = math.tan(math.radians(half_angle))
    cone_top = (R_c - R_i) / tan_t

    H_static = _bed_height_from_inventory(
        mass_kg=float(case_V["inventory_g"]) * 1e-3,
        rho_p_kg_m3=float(case_V["rho_p_kg_per_m3"]),
        eps_s_bed=float(case_V["eps_s_bed"]),
        R_i_m=R_i,
        cone_half_angle_deg=half_angle,
    )
    eps_s_bed = float(case_V["eps_s_bed"])

    fig, ax = plt.subplots(figsize=(7.0, 7.5))
    ax.add_patch(
        mp.Rectangle((-R_c, 0), 2 * R_c, H_dom, facecolor="#f0f0f0",
                     edgecolor="black", linewidth=0.8, zorder=1)
    )
    # IC 2 packed bed
    ax.add_patch(
        mp.Rectangle((-R_c, 0), 2 * R_c, H_static, facecolor="#ffcc80",
                     alpha=0.65, edgecolor="none", zorder=2,
                     label=f"IC2 packed bed, eps_s={eps_s_bed}")
    )
    # IC 1 freeboard
    ax.add_patch(
        mp.Rectangle((-R_c, H_static), 2 * R_c, H_dom - H_static,
                     facecolor="#90caf9", alpha=0.35, edgecolor="none",
                     zorder=2, label="IC1 freeboard (gas)")
    )
    # IC 3 pre-opened spout channel
    ax.add_patch(
        mp.Rectangle((-R_i, 0), 2 * R_i, H_static, facecolor="#c5e1a5",
                     alpha=0.6, edgecolor="none", zorder=3,
                     label=f"IC3 pre-opened spout (|x|<={R_i*1e3:.1f} mm)")
    )
    # MI orifice bar
    ax.add_patch(
        mp.Rectangle((-R_i, -0.004), 2 * R_i, 0.004, facecolor="red", zorder=4,
                     label=f"MI orifice, |x|<={R_i*1e3:.1f} mm")
    )
    # PO top
    ax.add_patch(
        mp.Rectangle((-R_c, H_dom), 2 * R_c, 0.004, facecolor="green",
                     zorder=4, label="PO top (full width)")
    )
    # Cone walls (Y_CONE quadric)
    cone_z = np.linspace(0, cone_top, 50)
    cone_r_right = R_i + cone_z * tan_t
    ax.plot(cone_r_right, cone_z, "k-", linewidth=1.6, zorder=5,
            label=f"Y_CONE wall ({2*half_angle:.0f} deg cone, resolved cut-cell)")
    ax.plot(-cone_r_right, cone_z, "k-", linewidth=1.6, zorder=5)
    ax.plot([R_c, R_c], [cone_top, H_dom], "k-", linewidth=1.6, zorder=5)
    ax.plot([-R_c, -R_c], [cone_top, H_dom], "k-", linewidth=1.6, zorder=5)
    ax.axvline(0, color="black", linestyle=":", linewidth=0.6, alpha=0.5)
    ax.text(0, H_dom * 1.005, "axis (x=0)", va="bottom", ha="center", fontsize=8)
    ax.axhline(H_static, color="black", linewidth=0.5, linestyle=":")
    ax.text(R_c * 0.95, H_static * 1.03,
            f"H_static={H_static*1000:.1f} mm (H0/Dc={H_static/(2*R_c):.3f})",
            fontsize=8, ha="right")

    # Monitor markers
    ax.plot(0.0, 0.0, "r*", ms=10, zorder=6)
    ax.text(R_i * 1.3, -0.008, "M1 P_inlet (1 kHz)", fontsize=8)
    ax.plot(0.0, H_static, "b*", ms=10, zorder=6)
    ax.text(R_c * 0.1, H_static + 0.005, "M2 P_top_bed (1 kHz)", fontsize=8)
    ax.plot(0.0, 0.5 * H_dom, "g*", ms=10, zorder=6)
    ax.text(R_c * 0.3, 0.5 * H_dom,
            "M4 centerline v_g (line)", fontsize=8)

    ax.set_xlim(-R_c * 1.3, R_c * 1.3)
    ax.set_ylim(-0.02, H_dom * 1.03)
    ax.set_aspect("equal")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m] (axial)")
    cells_total = int(v2d["imax"]) * int(v2d["jmax"]) * int(v2d["kmax"])
    dx = (2.0 * R_c) / int(v2d["imax"])
    orifice_cells = (2.0 * R_i) / dx
    ax.set_title(
        f"V_2d (Day 8): ORNL/UTK cold mockup, 2D cartesian cut-cell + Y_CONE quadric\n"
        f"D_c={2*R_c*1e3:.1f} mm, D_i={2*R_i*1e3:.1f} mm, 60 deg cone; "
        f"grid {v2d['imax']} x {v2d['jmax']} x {v2d['kmax']} = {cells_total} cells; "
        f"orifice spans {orifice_cells:.1f} cells",
        fontsize=10,
    )
    ax.legend(loc="upper right", fontsize=7, framealpha=0.92)
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    return out_path


def main() -> None:
    out = ROOT / "results" / "fig_d8_V_2d_geometry.png"
    print("wrote", render_V_2d_schematic(out))


if __name__ == "__main__":
    main()

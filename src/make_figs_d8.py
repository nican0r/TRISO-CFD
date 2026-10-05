"""Render Day-8 geometry screenshots for visual confirmation.

- results/fig_d8_V_3d_geometry.png : 3D vessel inner-wall STL (surface of
  revolution of the cone + cylinder), drawn as an isometric wire+shaded
  projection using matplotlib's mpl_toolkits.mplot3d so the renderer stays
  entirely in Agg (no GL context required).
- results/fig_d8_V_2d_geometry.png : 2D axisymmetric (r, z) domain diagram
  with IC regions, BCs, and the three monitor points annotated.

Both figures encode geometry-level inputs only (no time averaging), so no
averaging-window caption is required per CLAUDE.md rule 6.
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.patches as mp  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import yaml  # noqa: E402
ROOT = Path(__file__).resolve().parent.parent


def _load_geometry() -> tuple[dict, dict, dict, dict]:
    geom = yaml.safe_load((ROOT / "params" / "geometry.yaml").read_text())
    return (
        geom["bed_0p076m"],
        geom["spouted_bed_V_2d"],
        geom["spouted_bed_V_3d"],
        geom["spouted_bed_V_initial"],
    )


def _read_stl_triangles(stl_path: Path) -> np.ndarray:
    """Return an (n, 3, 3) array of triangle vertices from an ASCII STL."""
    text = stl_path.read_text()
    pat = re.compile(
        r"facet normal[^\n]*\n\s*outer loop"
        r"\s*vertex\s+([-+eE.0-9]+)\s+([-+eE.0-9]+)\s+([-+eE.0-9]+)"
        r"\s*vertex\s+([-+eE.0-9]+)\s+([-+eE.0-9]+)\s+([-+eE.0-9]+)"
        r"\s*vertex\s+([-+eE.0-9]+)\s+([-+eE.0-9]+)\s+([-+eE.0-9]+)",
        re.MULTILINE,
    )
    flat = np.array(
        [[float(x) for x in m] for m in pat.findall(text)], dtype=float
    )
    return flat.reshape(-1, 3, 3)


def render_V_3d_screenshot(out_path: Path) -> Path:
    """Side-by-side figure: axial (x-y) slice of the STL and a top-down (x-z) ring.

    An isometric 3D view is misleading at the vessel's 10:1 axial:radial aspect
    ratio, so instead we show the two orthogonal projections that uniquely
    confirm the geometry: the axial profile (cone + cylinder) and the top
    cross-section (circular column).
    """
    bed, _v2d, v3d, _init = _load_geometry()
    R_c = 0.5 * bed["D_c_m"]
    R_i = 0.5 * bed["D_i_m"]
    H_static = v3d["H_static_m"]
    H_dom = v3d["H_dom_m"]

    tris = _read_stl_triangles(ROOT / "cases" / "V_3d" / "geometry.stl")
    cells_total = int(v3d["imax"]) * int(v3d["jmax"]) * int(v3d["kmax"])

    fig, (axL, axR) = plt.subplots(
        1, 2, figsize=(10, 7.5), gridspec_kw={"width_ratios": [2, 1]}
    )

    # --- Left: axial profile (project all facets to x-y) ---
    for tri in tris:
        axL.plot(
            tri[:, 0].tolist() + [tri[0, 0]],
            tri[:, 1].tolist() + [tri[0, 1]],
            color="#888", linewidth=0.25, alpha=0.5,
        )
    axL.plot([-R_i, R_i], [0, 0], color="red", linewidth=3,
             label=f"MI orifice (|x|<={R_i*1e3:.2f} mm)")
    axL.plot([-R_c, R_c], [H_dom, H_dom], color="green", linewidth=2,
             label=f"PO top (y={H_dom:.2f} m)")
    axL.plot([-R_c, R_c], [H_static, H_static], color="gold", linewidth=1.5,
             linestyle="--", label=f"Static bed top (y={H_static:.2f} m)")
    axL.set_xlim(-R_c * 1.3, R_c * 1.3)
    axL.set_ylim(-0.02, H_dom * 1.03)
    axL.set_aspect("equal")
    axL.set_xlabel("x [m]")
    axL.set_ylabel("y [m] (axial)")
    axL.set_title("V_3d: axial profile (x-y projection of STL)")
    axL.legend(loc="upper right", fontsize=8, framealpha=0.9)
    axL.grid(alpha=0.3)

    # --- Right: top view at the cylinder cross-section ---
    theta = np.linspace(0, 2 * np.pi, 120)
    axR.plot(R_c * np.cos(theta), R_c * np.sin(theta),
             color="black", linewidth=1.2, label=f"Column r=R_c={R_c*1e3:.1f} mm")
    axR.fill(R_i * np.cos(theta), R_i * np.sin(theta),
             color="red", alpha=0.5, label=f"Orifice r<={R_i*1e3:.2f} mm")
    # overlay STL top rim as scatter
    top_mask = tris[:, :, 1].max(axis=1) > H_dom - 1e-9
    top_pts = tris[top_mask].reshape(-1, 3)
    axR.scatter(top_pts[:, 0], top_pts[:, 2], s=3, c="#555", alpha=0.6,
                label="STL top-rim vertices")
    axR.set_xlim(-R_c * 1.25, R_c * 1.25)
    axR.set_ylim(-R_c * 1.25, R_c * 1.25)
    axR.set_aspect("equal")
    axR.set_xlabel("x [m]")
    axR.set_ylabel("z [m]")
    axR.set_title("V_3d: top view (x-z)")
    axR.legend(loc="upper right", fontsize=7, framealpha=0.9)
    axR.grid(alpha=0.3)

    fig.suptitle(
        f"V_3d (Day 8): cylinder + 60 deg cone + orifice (surface of revolution)\n"
        f"D_c = {bed['D_c_m']*1e3:.1f} mm, D_i = {bed['D_i_m']*1e3:.2f} mm, "
        f"H_dom = {H_dom:.2f} m; grid "
        f"{v3d['imax']}x{v3d['jmax']}x{v3d['kmax']} = {cells_total} cells, "
        f"{len(tris)} STL facets",
        fontsize=10,
    )
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    return out_path


def render_V_2d_schematic(out_path: Path) -> Path:
    bed, v2d, _v3d, init = _load_geometry()
    R_c = 0.5 * bed["D_c_m"]
    R_i = 0.5 * bed["D_i_m"]
    H_static = v2d["H_static_m"]
    H_dom = v2d["H_dom_m"]
    cone_top = v2d["cone_top_y_m"]

    fig, ax = plt.subplots(figsize=(7.0, 7.5))
    # Domain bounding box, -R_c <= x <= +R_c
    ax.add_patch(
        mp.Rectangle((-R_c, 0), 2 * R_c, H_dom, facecolor="#f0f0f0",
                     edgecolor="black", linewidth=0.8, zorder=1)
    )
    # IC2 packed bed (full width)
    ax.add_patch(
        mp.Rectangle((-R_c, 0), 2 * R_c, H_static, facecolor="#ffcc80", alpha=0.65,
                     edgecolor="none", zorder=2,
                     label=f"IC2 packed bed, eps_s={init['eps_s_bed']}")
    )
    # IC1 freeboard
    ax.add_patch(
        mp.Rectangle((-R_c, H_static), 2 * R_c, H_dom - H_static,
                     facecolor="#90caf9", alpha=0.35, edgecolor="none", zorder=2,
                     label="IC1 freeboard (gas)")
    )
    # MI orifice bar across |x| <= R_i
    ax.add_patch(
        mp.Rectangle((-R_i, -0.004), 2 * R_i, 0.004, facecolor="red", zorder=4,
                     label=f"MI orifice, |x|<={R_i*1e3:.2f} mm")
    )
    # PO top (full width)
    ax.add_patch(
        mp.Rectangle((-R_c, H_dom), 2 * R_c, 0.004, facecolor="green", zorder=4,
                     label="PO top (full width)")
    )
    # Resolved cone walls as SOLID lines on both sides (Y_CONE quadric wall)
    cone_z = np.linspace(0, cone_top, 20)
    cone_r_right = R_i + cone_z * np.tan(np.radians(30.0))
    ax.plot(cone_r_right, cone_z, "k-", linewidth=1.6, zorder=5,
            label="Y_CONE quadric wall (60 deg cone, resolved)")
    ax.plot(-cone_r_right, cone_z, "k-", linewidth=1.6, zorder=5)
    # Cylinder straight walls above the cone
    ax.plot([R_c, R_c], [cone_top, H_dom], "k-", linewidth=1.6, zorder=5)
    ax.plot([-R_c, -R_c], [cone_top, H_dom], "k-", linewidth=1.6, zorder=5)
    # Axis reference
    ax.axvline(0, color="black", linestyle=":", linewidth=0.6, alpha=0.5)
    ax.text(0, H_dom * 1.005, "axis (x=0)", va="bottom", ha="center", fontsize=8)
    ax.axhline(H_static, color="black", linewidth=0.5, linestyle=":")
    ax.text(R_c * 0.95, H_static * 1.03, f"H_static={H_static} m",
            fontsize=8, ha="right")
    # Monitor markers (symmetric about x=0 for the schematic — probes still
    # average over the full strip in the .mfx)
    ax.plot(0.0, 0.0, "r*", ms=10, zorder=6)
    ax.text(R_i * 1.3, -0.01, "M1 P_inlet", fontsize=8)
    ax.plot(0.0, H_static, "b*", ms=10, zorder=6)
    ax.text(R_c * 0.1, H_static + 0.005, "M2 P_top_bed", fontsize=8)
    ax.plot(0.0, H_static * 0.5, "m*", ms=10, zorder=6)
    ax.text(R_i * 1.3, H_static * 0.5, "M3 spout_mid (ep_s, v_s)", fontsize=8)

    ax.set_xlim(-R_c * 1.3, R_c * 1.3)
    ax.set_ylim(-0.015, H_dom * 1.03)
    ax.set_aspect("equal")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m] (axial)")
    cells_total = int(v2d["imax"]) * int(v2d["jmax"]) * int(v2d["kmax"])
    dx = (2.0 * R_c) / int(v2d["imax"])
    orifice_cells = (2.0 * R_i) / dx
    ax.set_title(
        f"V_2d (Day 8 rev 4): 2D Cartesian cut-cell + Y_CONE quadric (no_k)\n"
        f"grid {v2d['imax']} x {v2d['jmax']} x {v2d['kmax']} = {cells_total} cells; "
        f"orifice spans {orifice_cells:.1f} cells across the diameter",
        fontsize=10,
    )
    ax.legend(loc="upper right", fontsize=7, framealpha=0.92)
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    return out_path


def main() -> None:
    results = ROOT / "results"
    v2d_png = render_V_2d_schematic(results / "fig_d8_V_2d_geometry.png")
    print("wrote", v2d_png)
    v3d_png = render_V_3d_screenshot(results / "fig_d8_V_3d_geometry.png")
    print("wrote", v3d_png)


if __name__ == "__main__":
    main()

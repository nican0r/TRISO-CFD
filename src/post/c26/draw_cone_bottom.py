"""Schematic of the three C26 cone-bottom / gas-source layouts (x-y section at z = 0).

Draws, on the 2.1 mm grid, the open apex (flat cap at y = 0, 2 x 2 source on
the 24-cell grid), the closed apex (2 x 2 source in the tip row, 24-cell grid,
the control inlet) and the full-face bore (one source cell walled in, 25-cell
grid, the single-cell inlet).  Dimensions come from params/geometry.yaml
(bed_C26, spouted_bed_C26_3d); see implementation-steps/03-geometry-and-mesh.md.

    MPLBACKEND=Agg .venv/bin/python src/post/c26/draw_cone_bottom.py
"""
import argparse
import math
import sys
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))


def load_geometry():
    g = yaml.safe_load((REPO / "params" / "geometry.yaml").read_text())
    bed, mesh = g["bed_C26"], g["spouted_bed_C26_3d"]
    dx = bed["D_c_m"] / mesh["n_cells_across_D"]
    return dict(
        dx=1e3 * dx,
        R_t=1e3 * 0.5 * bed["D_throat_m"],
        tan=math.tan(math.radians(0.5 * bed["cone_angle_deg"])),
        delta=1e3 * dx * mesh["stl_fullface_wall_offset_frac"],
    )


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=REPO / "implementation-steps" / "figures" / "cone_bottom_variants.svg")
    a = ap.parse_args(argv)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon, Rectangle

    p = load_geometry()
    dx, R_t, tan = p["dx"], p["R_t"], p["tan"]
    x_lim, y_lo, y_hi = 3 * dx, -1.5 * dx, 3.3 * dx
    r = lambda y: R_t + y * tan                          # cone wall radius [mm]
    y_apex = -R_t / tan
    a_bore = 0.5 * dx + p["delta"]                      # bore half-width
    y_meet = (a_bore - R_t) / tan                       # cone meets bore wall (mid-side)

    def solid(wall):                                    # wall: right-hand wall points, bottom to top
        right = wall + [(r(y_hi), y_hi), (x_lim, y_hi), (x_lim, y_lo)]
        left = [(-x, y) for x, y in reversed(wall)]
        return [(-x_lim, y_lo), (-x_lim, y_hi), (-r(y_hi), y_hi)] + left[:-1] + right

    panels = [
        dict(title="open apex (flat cap)\n2 × 2 source, 24-cell grid",
             wall=[(0.0, 0.0), (R_t, 0.0)], src=(-dx, 0.0, 2 * dx), x0=0.0,
             note="runs to 2026-10-23"),
        dict(title="closed apex\n2 × 2 tip source, 24-cell grid",
             wall=[(0.0, y_apex)], src=(-dx, -dx, 2 * dx), x0=0.0,
             note="control inlet (_apex_v232)"),
        dict(title="full-face bore\none source cell, 25-cell grid",
             wall=[(0.0, -dx - p["delta"]), (a_bore, -dx - p["delta"]), (a_bore, y_meet)],
             src=(-0.5 * dx, -dx, dx), x0=0.5 * dx,
             note="single-cell inlet (_fullface_bore_mom_nointerp_v232)"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(11, 4.6), facecolor="white")
    for ax, pn in zip(axes, panels):
        ax.set_facecolor("white")
        x_s, y_s, w = pn["src"]
        ax.add_patch(Rectangle((x_s, y_s), w, dx, color="#f4a259", zorder=1))
        ax.add_patch(Polygon(solid(pn["wall"]), closed=True, facecolor="#c9c9c9", edgecolor="none", zorder=2))
        wall = pn["wall"] + [(r(y_hi), y_hi)]
        ax.plot([x for x, _ in wall], [y for _, y in wall], color="k", lw=1.4, zorder=4)
        ax.plot([-x for x, _ in wall], [y for _, y in wall], color="k", lw=1.4, zorder=4)
        k = math.ceil(x_lim / dx) + 1
        for i in range(-k, k + 1):
            ax.axvline(pn["x0"] + i * dx, color="#7a7a7a", lw=0.5, ls=(0, (3, 3)), zorder=3)
        for j in range(-2, 5):
            ax.axhline(j * dx, color="#7a7a7a", lw=0.5, ls=(0, (3, 3)), zorder=3)
        ax.axvline(0.0, color="#1f5fa8", lw=0.9, ls="-.", zorder=3)
        ax.annotate("", xy=(0.0, y_s + dx + 2.2), xytext=(0.0, y_s + dx + 0.2),
                    arrowprops=dict(arrowstyle="-|>", color="#b5431c", lw=1.6), zorder=5)
        ax.set_xlim(-x_lim, x_lim)
        ax.set_ylim(y_lo, y_hi)
        ax.set_aspect("equal")
        ax.set_title(pn["title"], fontsize=10)
        ax.set_xlabel("x [mm]")
        ax.text(0.5, -0.27, pn["note"], transform=ax.transAxes, ha="center", fontsize=8)
    axes[0].set_ylabel("y [mm]")
    fig.text(0.5, 0.01,
             f"Section through the axis (blue). Grey: solid; orange: source cell(s), fluid part; dashed: "
             f"{dx:g} mm grid. Cone wall 30° from the axis, throat radius {R_t:.2f} mm at y = 0.",
             ha="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    a.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(a.out, facecolor="white")
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()

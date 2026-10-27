"""Source-comparison map metrics, native PIV engine against PIVlab (2026-10-27).

The c02 point-source runs were scored with the native engine
(src/post/piv_c26.py); the PIVlab re-run (notes/c26_pivlab_rerun.md) covered
the control only.  This applies the metrics that separated the sources
(notes/c26_conry2026_case.md "Sweep source choice") to every map available:
fountain reach (row-max |v| >= 0.2 m/s), band ratios to the experiment
(row-max |v|, rows 0.70-0.95 / 0.5-0.7 / < 0.5 of the map height), map mean
and row-mean by fifth over C26's model (Fig. 4b, 12 slpm), CC against the
experiment and the spout-peak x offset in the bottom fifth (|v|-weighted x of
the cells >= half the row max).  The native rows reproduce the notes' values.
Window t = 2.4-3.0 s for every map.

    .venv/bin/python notes/c26_map_metrics_pivlab_vs_native.py
"""
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from src.post.piv_c26 import load_fig4  # noqa: E402

REACH_V = 0.2          # m/s, point_3s fountain_reach threshold
X_HALF_MM = 25.2       # map half-width (R_c), Fig. 4 grid

MAPS = {
    "control native": "results/c26_piv_map_c02_apex_v232.csv",
    "control PIVlab, C26-style frames + defaults":
        "results/pivlab_export/c02_apex_v232/pivlab_c26img/map_fig4grid.csv",
    "control PIVlab, C26-style frames, CLAHE off":
        "results/pivlab_export/c02_apex_v232/pivlab_c26img_noclahe/map_fig4grid.csv",
    "point mass-only native": "results/c26_piv_map_c02_point_v232.csv",
    "point + momentum native": "results/c26_piv_map_c02_point_mom_v232.csv",
    "2x2 + momentum native": "results/c26_piv_map_c02_apex_mom_v232.csv",
}


def load(rel):
    """Map [m/s], row 0 at the bottom (the CSVs store the top row first)."""
    return np.flipud(np.loadtxt(REPO / rel, delimiter=",", comments="#"))


def main():
    exp, mod = load_fig4("exp", 12), load_fig4("model", 12)
    n, ncol = mod.shape
    f = (np.arange(n) + 0.5) / n
    x = (np.arange(ncol) + 0.5) / ncol * 2 * X_HALF_MM - X_HALF_MM

    def reach(m):
        hit = np.nonzero(np.nanmax(m, axis=1) >= REACH_V)[0]
        return (hit.max() + 0.5) / n if hit.size else 0.0

    def band(m, lo, hi):
        r = (f >= lo) & (f < hi)
        return np.nanmean(np.nanmax(m[r], 1)) / np.nanmean(np.nanmax(exp[r], 1))

    def fifths(m):
        return [np.nanmean(m[(f >= a) & (f < a + .2)]) / np.nanmean(mod[(f >= a) & (f < a + .2)])
                for a in (0, .2, .4, .6, .8)]

    def offset(m):
        xs = []
        for row in m[f < 0.2]:
            k = row >= 0.5 * np.nanmax(row)
            xs.append(np.nansum(row[k] * x[k]) / np.nansum(row[k]))
        return float(np.mean(xs))

    print(f"{'map':45s} reach  bands .7-.95/.5-.7/<.5  mean/model  fifths/model (bottom->top)"
          "      CC vs exp  offset [mm]")
    for name, rel in MAPS.items():
        m = load(rel)
        b = (band(m, .7, .95), band(m, .5, .7), band(m, 0, .5))
        print(f"{name:45s} {reach(m):.2f}   {b[0]:.2f}/{b[1]:.2f}/{b[2]:.2f}"
              f"            {np.nanmean(m) / np.nanmean(mod):.2f}       "
              + " ".join(f"{v:.2f}" for v in fifths(m))
              + f"   {np.corrcoef(m.ravel(), exp.ravel())[0, 1]:.3f}    {offset(m):+.1f}")
    print(f"C26 model reach {reach(mod):.2f}, offset {offset(mod):+.1f} mm; "
          f"experiment reach {reach(exp):.2f}, offset {offset(exp):+.1f} mm")


if __name__ == "__main__":
    main()

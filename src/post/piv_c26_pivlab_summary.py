"""Collect the PIVlab re-run scores of the C26 PIV benchmark into one table and
one figure.

Reads <export_root>/<run>/<variant>/score.json written by
src/post/piv_c26_pivlab_score.m (MATLAB; PIVlab fields from
src/post/piv_c26_pivlab_run.m) and, for the side-by-side column, the in-repo
engine references <export_root>/<run>/reference_native.json written by
src.post.piv_c26_export_frames.  Writes

    results/c26_pivlab_rerun.csv         one row per run x variant
    results/fig_c26_pivlab_maps.png      Fig. 4-style maps: C26 experiment,
                                         C26 model (both digitized), PIVlab
                                         default, PIVlab without CLAHE
    results/fig_c26_pivlab_wall.png      wall velocity per run and variant
                                         next to C26 Fig. 6d

    .venv/bin/python -m src.post.piv_c26_pivlab_summary [--runs c01_apex_v232 ...]

All values SI.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
import yaml

REPO = Path(__file__).resolve().parents[2]
EXPORT_ROOT = REPO / "results" / "pivlab_export"
VARIANTS = ("pivlab_default", "pivlab_noclahe", "pivlab_ces25", "pivlab_ces25_noclahe",
            "pivlab_c26img", "pivlab_c26img_noclahe", "pivlab_c26res", "pivlab_c26res_noclahe")
VARIANT_LABEL = {
    "pivlab_default": "PIVlab defaults (CLAHE on)",
    "pivlab_noclahe": "PIVlab defaults, CLAHE off",
    "pivlab_ces25": "CES25 binarization + defaults",
    "pivlab_ces25_noclahe": "CES25 binarization, CLAHE off",
    "pivlab_c26img": "C26-style frames (binary disks, 8-bit)\n+ PIVlab defaults",
    "pivlab_c26img_noclahe": "C26-style frames, CLAHE off",
    "pivlab_c26res": "C26-style frames at 0.05 mm/px,\n64/32/16 px + PIVlab defaults",
    "pivlab_c26res_noclahe": "C26-style frames at 0.05 mm/px,\n64/32/16 px, CLAHE off",
}
METRICS = ("CC", "SSI", "MSE_m2_s2")


def load_table3(path=REPO / "params" / "sweep_C26.yaml"):
    """C26 Table 3 (CC, SSI, MSE) per case id (c01-c03) from the sweep file."""
    y = yaml.safe_load(path.read_text())
    t3 = y.get("piv_targets_table3") or {}
    return {str(k): {m: float(v[m]) for m in METRICS if m in v} for k, v in t3.items() if isinstance(v, dict)}


def load_wall_targets(path=REPO / "params" / "piv_c26.yaml"):
    y = yaml.safe_load(path.read_text())["wall_targets"]
    return ({float(k): float(v) for k, v in y["exp_m_s"].items()},
            {float(k): float(v) for k, v in y["model_m_s"].items()})


def collect(export_root, runs, variants=VARIANTS):
    """Rows (dicts) for every run x variant with a score.json."""
    table3 = load_table3()
    wall_exp, wall_model = load_wall_targets()
    rows, maps = [], {}
    for run in runs:
        d = Path(export_root) / run
        ref = {}
        f = d / "reference_native.json"
        if f.exists():
            ref = json.loads(f.read_text())
        for v in variants:
            sf = d / v / "score.json"
            if not sf.exists():
                continue
            s = json.loads(sf.read_text())
            flow = float(s["flow_slpm"])
            row = {"run": run, "case": s["case"], "flow_slpm": flow, "variant": v,
                   "window_s": f"{s['window_t0_s']:.3f}-{s['window_t1_s']:.3f}",
                   "n_pairs": s["n_pairs_spout"], "pivlab_version": s.get("pivlab_version", ""),
                   "clahe": s["pivlab_settings"]["pre"]["clahe"],
                   "binarize": s["pivlab_settings"]["pre"]["binarize"],
                   "map_mean_m_s": s["map_mean_m_s"], "map_max_m_s": s["map_max_m_s"],
                   "spout_rejected_frac": s["spout_rejected_frac"]}
            for m in METRICS:
                row[m] = s["metrics_vs_exp"][m]
                row[f"{m}_c26model_same_pipeline"] = s["metrics_c26_model_vs_exp"][m]
                row[f"{m}_table3"] = table3.get(str(s["case"]), {}).get(m, np.nan)
                row[f"{m}_vs_c26model"] = s["metrics_vs_c26_model"][m]
            # time-mean conventions (score.json time_mean_conventions): how the
            # empty windows enter the 300-pair mean; 'interp' is the headline
            for conv in ("skip", "pivlab_mean", "skip_fill"):
                c = (s.get("time_mean_conventions") or {}).get(conv)
                for m in METRICS:
                    row[f"{m}_{conv}"] = c["metrics_vs_exp"][m] if c else np.nan
                    row[f"{m}_{conv}_vs_c26model"] = c["metrics_vs_c26_model"][m] if c else np.nan
                row[f"nan_frac_{conv}"] = c["nan_frac_on_grid"] if c else np.nan
            row["wall_piv_m_s"] = s["wall_piv_m_s"]
            row["wall_piv_sd_m_s"] = s["wall_piv_sd_m_s"]
            row["wall_vectors_in_box"] = s["wall_vectors_in_box"]
            row["wall_exp_m_s"] = wall_exp.get(flow, np.nan)
            row["wall_c26model_m_s"] = wall_model.get(flow, np.nan)
            row["wall_over_exp"] = row["wall_piv_m_s"] / row["wall_exp_m_s"] if row["wall_exp_m_s"] else np.nan
            # side by side: the in-repo native engine on the same frames (not a reference)
            # headline wall velocity: the parcels visible in the box (MFiX velocities,
            # no PIV), as C26 reports its model wall velocity "from particles at the wall"
            row["wall_parcels_m_s"] = ref.get("wall_visible_parcels_m_s", np.nan)
            row["wall_parcels_over_exp"] = row["wall_parcels_m_s"] / row["wall_exp_m_s"] if row["wall_exp_m_s"] else np.nan
            row["repo_native_wall_m_s"] = ref.get("wall_piv_m_s", np.nan)
            for m in METRICS:
                row[f"repo_native_{m}"] = (ref.get("metrics_vs_exp") or {}).get(m, np.nan)
            rows.append(row)
            mf = d / v / "map_fig4grid.csv"
            if mf.exists():
                maps[(run, v)] = np.loadtxt(mf, delimiter=",")
            for conv in ("skip", "pivlab_mean", "skip_fill"):
                mf = d / v / f"map_fig4grid_{conv}.csv"
                if mf.exists():
                    maps[(run, v, conv)] = np.loadtxt(mf, delimiter=",")
    return rows, maps


def write_csv(rows, path):
    if not rows:
        raise SystemExit("no score.json found")
    keys = list(rows[0])
    with Path(path).open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})


def digitized(kind, flow):
    return np.loadtxt(REPO / "data" / "digitized" / f"c26_fig4_{kind}_{flow:g}slpm.csv", delimiter=",", comments="#")


def plot_maps(rows, maps, path, variants=("pivlab_c26img", "pivlab_c26img_noclahe", "pivlab_default", "pivlab_noclahe")):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    p = yaml.safe_load((REPO / "params" / "piv_c26.yaml").read_text())["fig4"]
    vmin, vmax = p["v_min_m_s"], p["v_max_m_s"]
    runs = sorted({(r["run"], r["flow_slpm"]) for r in rows}, key=lambda t: t[1])
    cols = [("C26 experiment\n(Fig. 4a, digitized)", "exp"), ("C26 model\n(Fig. 4b, digitized)", "model")] + \
           [(VARIANT_LABEL[v], v) for v in variants]
    fig, axes = plt.subplots(len(runs), len(cols), figsize=(3.2 * len(cols), 3.3 * len(runs)), squeeze=False,
                             constrained_layout=True)
    by = {(r["run"], r["variant"]): r for r in rows}
    for i, (run, flow) in enumerate(runs):
        for j, (title, key) in enumerate(cols):
            ax = axes[i, j]
            if key in ("exp", "model"):
                a = digitized(key, flow)
                sub = ""
                if key == "model":
                    r = next((by[(run, v)] for v in variants if (run, v) in by), None)
                    if r:
                        sub = (f"\nvs exp: CC {r['CC_c26model_same_pipeline']:.2f} SSI {r['SSI_c26model_same_pipeline']:.2f}"
                               f" MSE {r['MSE_m2_s2_c26model_same_pipeline']:.4f}")
            else:
                a = maps.get((run, key))
                r = by.get((run, key))
                sub = f"\nvs exp: CC {r['CC']:.2f} SSI {r['SSI']:.2f} MSE {r['MSE_m2_s2']:.4f}" if r else ""
            if a is None:
                ax.set_axis_off()
                continue
            im = ax.imshow(a, cmap="jet", vmin=vmin, vmax=vmax, aspect="auto", extent=(-25.2, 25.2, 0, 39.8))
            ax.set_title((f"{run}, {flow:g} slpm: " if j == 0 else "") + title + sub, fontsize=7)
            ax.set_xlabel("x [mm]", fontsize=7)
            if j == 0:
                ax.set_ylabel("height above annulus surface [mm]", fontsize=7)
            ax.tick_params(labelsize=6)
    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.6, label="|time-mean velocity| [m/s]")
    r0 = rows[0]
    fig.text(0.5, -0.03,
             f"C26 PIV benchmark re-run with PIVlab {r0['pivlab_version']}. Time averages over t = {r0['window_s']} s "
             f"({r0['n_pairs']} consecutive pairs 2 ms apart, after the 2.4 s startup transient). Maps: |mean vector| on the "
             f"digitized Fig. 4a grid; metrics (corrcoef, ssim, immse) vs the digitized experiment panel of the same flow.",
             ha="center", fontsize=7, wrap=True)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_timemean(rows, maps, path, variant="pivlab_c26img"):
    """C26 model map next to our map under each time-mean convention (how the
    empty windows enter the 300-pair mean), for one PIVlab variant."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    p = yaml.safe_load((REPO / "params" / "piv_c26.yaml").read_text())["fig4"]
    vmin, vmax = p["v_min_m_s"], p["v_max_m_s"]
    runs = sorted({(r["run"], r["flow_slpm"]) for r in rows if r["variant"] == variant}, key=lambda t: t[1])
    if not runs:
        return
    by = {(r["run"], r["variant"]): r for r in rows}
    cols = [("C26 model (Fig. 4b, digitized)", None),
            ("interp: per-pair interpolation, then mean\n(headline)", "interp"),
            ("skip: original vectors only, never-valid = 0", "skip"),
            ("skip + inpaint of never-valid points", "skip_fill"),
            ("PIVlab temporal mean rule:\n> 75 % interpolated -> NaN (= 0)", "pivlab_mean")]
    fig, axes = plt.subplots(len(runs), len(cols), figsize=(3.2 * len(cols), 3.3 * len(runs)), squeeze=False,
                             constrained_layout=True)
    for i, (run, flow) in enumerate(runs):
        r = by[(run, variant)]
        for j, (title, conv) in enumerate(cols):
            ax = axes[i, j]
            if conv is None:
                a = digitized("model", flow)
                sub = (f"\nvs exp: CC {r['CC_c26model_same_pipeline']:.2f} SSI {r['SSI_c26model_same_pipeline']:.2f}"
                       f" MSE {r['MSE_m2_s2_c26model_same_pipeline']:.4f}")
            elif conv == "interp":
                a = maps.get((run, variant))
                sub = f"\nvs exp: CC {r['CC']:.2f} SSI {r['SSI']:.2f} MSE {r['MSE_m2_s2']:.4f}"
            else:
                a = maps.get((run, variant, conv))
                sub = (f"\nvs exp: CC {r[f'CC_{conv}']:.2f} SSI {r[f'SSI_{conv}']:.2f} MSE {r[f'MSE_m2_s2_{conv}']:.4f}"
                       f"; NaN {100 * r[f'nan_frac_{conv}']:.0f} %")
            if a is None:
                ax.set_axis_off()
                continue
            im = ax.imshow(a, cmap="jet", vmin=vmin, vmax=vmax, aspect="auto", extent=(-25.2, 25.2, 0, 39.8))
            ax.set_title((f"{run}, {flow:g} slpm: " if j == 0 else "") + title + sub, fontsize=7)
            ax.tick_params(labelsize=6)
    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.6, label="|time-mean velocity| [m/s]")
    r0 = rows[0]
    fig.text(0.5, -0.03, f"Time-mean conventions for the {VARIANT_LABEL[variant].replace(chr(10), ' ')} run, "
             f"t = {r0['window_s']} s; see notes/c26_pivlab_rerun.md.", ha="center", fontsize=7)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_wall(rows, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    runs = sorted({(r["run"], r["flow_slpm"]) for r in rows}, key=lambda t: t[1])
    variants = [v for v in VARIANTS if any(r["variant"] == v for r in rows)]
    by = {(r["run"], r["variant"]): r for r in rows}
    runs = [rf for rf in runs if all((rf[0], v) in by for v in variants[:1])]
    fig, ax = plt.subplots(figsize=(8, 5.5))
    x = np.arange(len(runs))
    w = 0.8 / (len(variants) + 4)
    ax.bar(x - (len(variants) / 2 + 1) * w, [by[(run, variants[0])]["wall_parcels_m_s"] for run, _ in runs], w,
           color="k", label="parcels in the box (model velocities, headline; as C26 Fig. 6)")
    for k, v in enumerate(variants):
        vals = [by.get((run, v), {}).get("wall_piv_m_s", np.nan) for run, _ in runs]
        ax.bar(x + (k - len(variants) / 2) * w, vals, w, label=VARIANT_LABEL[v])
    ax.bar(x + (len(variants) / 2) * w, [by[(run, variants[0])]["repo_native_wall_m_s"] for run, _ in runs], w,
           color="0.6", label="repo native engine (same frames)")
    ax.bar(x + (len(variants) / 2 + 1) * w, [by[(run, variants[0])]["wall_c26model_m_s"] for run, _ in runs], w,
           color="tab:red", alpha=0.6, label="C26 model (Fig. 6d)")
    ax.bar(x + (len(variants) / 2 + 2) * w, [by[(run, variants[0])]["wall_exp_m_s"] for run, _ in runs], w,
           color="tab:blue", alpha=0.6, label="C26 experiment (Fig. 6d)")
    ax.set_xticks(x, [f"{run}\n{flow:g} slpm" for run, flow in runs], fontsize=8)
    ax.set_ylabel("mean wall velocity, downward + [m/s]")
    ax.set_yscale("symlog", linthresh=1e-3)
    ax.grid(axis="y", alpha=0.3, which="both")
    ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2)
    r0 = rows[0]
    ax.set_title(f"Wall velocity, ~4 mm box (C26 Fig. 1c): model parcels vs PIVlab {r0['pivlab_version']} on images; time mean over t = {r0['window_s']} s",
                 fontsize=9)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--export-root", type=Path, default=EXPORT_ROOT)
    ap.add_argument("--runs", nargs="+", default=["c01_apex_v232", "c02_apex_v232", "c03_apex_v232"])
    ap.add_argument("--out-dir", type=Path, default=REPO / "results")
    ap.add_argument("--tag", default="")
    ap.add_argument("--map-variants", nargs="+", default=["pivlab_c26img", "pivlab_c26img_noclahe", "pivlab_default", "pivlab_noclahe"],
                    help="PIVlab variants drawn as columns of the maps figure")
    ap.add_argument("--timemean-variant", default="pivlab_c26img", help="variant drawn in the time-mean figure")
    a = ap.parse_args(argv)
    rows, maps = collect(a.export_root, a.runs)
    csv_path = a.out_dir / f"c26_pivlab_rerun{a.tag}.csv"
    write_csv(rows, csv_path)
    plot_maps(rows, maps, a.out_dir / f"fig_c26_pivlab_maps{a.tag}.png", variants=tuple(a.map_variants))
    plot_wall(rows, a.out_dir / f"fig_c26_pivlab_wall{a.tag}.png")
    plot_timemean(rows, maps, a.out_dir / f"fig_c26_pivlab_timemean{a.tag}.png", variant=a.timemean_variant)
    for r in rows:
        print(f"{r['run']:14s} {r['variant']:22s} map {r['map_mean_m_s']:.3f} m/s | CC {r['CC']:.3f} SSI {r['SSI']:.3f} "
              f"MSE {r['MSE_m2_s2']:.4f} | Table 3 {r['CC_table3']:.3f} {r['SSI_table3']:.3f} {r['MSE_m2_s2_table3']:.4f} "
              f"| C26 model same pipeline {r['CC_c26model_same_pipeline']:.3f} {r['SSI_c26model_same_pipeline']:.3f} "
              f"{r['MSE_m2_s2_c26model_same_pipeline']:.4f} | wall parcels {r['wall_parcels_m_s']:.2e}, PIVlab {r['wall_piv_m_s']:.2e} (exp {r['wall_exp_m_s']:.2e}, "
              f"C26 model {r['wall_c26model_m_s']:.2e}, repo native {r['repo_native_wall_m_s']:.2e})")
    print(f"wrote {csv_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

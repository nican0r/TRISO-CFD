"""C26 Fig. 4 side by side (2026-10-29; implementation-steps/07-measurement-and-scoring.md).
Our synthetic-PIV spout maps next to C26's experiment (Fig. 4a) and model
(Fig. 4b), on C26's colour scale.  The maps and scores are the ones
src/post/piv_c26.py wrote to results/ (c26_piv_map_<name>.csv, c26_piv_<name>.csv);
the wall trace of its per-case figure is left out.

A case is cNN or cNN<suffix>, e.g. c03_p90 for the 90th-percentile map bottom.

Flow sweep, one row per case (C26 has a Fig. 4 panel at each flow):
    MPLBACKEND=Agg .venv/bin/python src/post/c26/c26_fig4_compose.py fullface_bore_mom_nointerp_v232 c01 c02 c03_p90
One-factor sweep at one flow, one row: C26's panels once, then each case
(C26 ran only its own charge and bare particles, so only the case at C26's
conditions carries scores):
    MPLBACKEND=Agg .venv/bin/python src/post/c26/c26_fig4_compose.py fullface_bore_mom_nointerp_v232 c08_p90 c02 c09 --sweep charge
Control configuration (closed apex, 2x2-cell mass-only tip source):
    MPLBACKEND=Agg .venv/bin/python src/post/c26/c26_fig4_compose.py apex_v232 c01 c02 c03 --label "2x2-cell control inlet"
--no-c26 drops C26's panels and the scores from the one-row layout (cases may then differ in flow):
    MPLBACKEND=Agg .venv/bin/python src/post/c26/c26_fig4_compose.py apex_v232 c01 c02 c03 --sweep flow --no-c26 --label "2x2-cell control inlet"
    MPLBACKEND=Agg .venv/bin/python src/post/c26/c26_fig4_compose.py fullface_bore_mom_nointerp_v232 c04 c05 c06 c07 --sweep diameter --no-c26
"""
import argparse
import csv
import json
import math
import re
import sys
import textwrap
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
from src.post.piv_c26 import load_fig4, load_params  # noqa: E402

RESULTS = REPO / "results"


def load_case(spec, tag, p):
    case, suffix = spec[:3], spec[3:]
    name = f"{case}_{tag}{suffix}"
    with (RESULTS / f"c26_piv_{name}.csv").open() as f:
        row = next(csv.DictReader(f))
    path = RESULTS / f"c26_piv_map_{name}.csv"
    header = path.open().readline()
    pct = re.search(r"p(\d+) of annulus", header)
    flow = float(row["flow_slpm"])
    if float(row["panel_flow_slpm"]) != flow:
        raise ValueError(f"{name}: no C26 Fig. 4 panel at {flow:g} slpm")
    m = json.loads(((REPO / row["run_dir"]).parent / "manifest.json").read_text())
    ours = np.flipud(np.loadtxt(path, delimiter=",", comments="#"))
    H = float(row["map_height_mm"])
    return dict(name=name, case=case, flow=flow, row=row, ours=ours,
                exp=load_fig4("exp", flow), model=load_fig4("model", flow),
                pct=float(pct.group(1)) if pct else float(p["spout_map"]["bed_surface_percentile"]),
                inventory_g=1e3 * m["inventory_kg"], H_static_mm=1e3 * m["H_static_m"],
                stage=m.get("coating_stage") or "bare", d_p_um=1e6 * m["d_p_m"],
                ext=(-0.5 * H * ours.shape[1] / ours.shape[0], 0.5 * H * ours.shape[1] / ours.shape[0], 0, H))


def scored(r, key=""):
    return math.isfinite(float(r["row"]["CC" + key]))


def scores(r, key=""):
    row = r["row"]
    return (f"vs exp: CC {float(row['CC' + key]):.3f}  SSI {float(row['SSI' + key]):.3f}  "
            f"MSE {float(row['MSE_m2_s2' + key]):.4f}")


def conditions(r):
    return f"{r['flow']:g} slpm, {r['inventory_g']:g} g, {r['stage']} {r['d_p_um']:.0f} um"


def engine(r):
    return r["row"].get("engine") or "native"   # rows scored before the engine column were native


def window(r):
    row = r["row"]
    return (f"{r['case']}: {engine(r)} PIV engine, {row['n_pairs']} frame pairs over t = {float(row['window_t0_s']):.3f}"
            f"-{float(row['window_t1_s']):.3f} s, map from bed surface {float(row['bed_surface_mm']):.1f} mm "
            f"(p{r['pct']:g} of annulus parcel y; static bed {r['H_static_mm']:.1f} mm), "
            f"{float(row['map_height_mm']):.1f} mm tall")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tag", help="run tag after cNN_, e.g. fullface_bore_mom_nointerp_v232")
    ap.add_argument("cases", nargs="+", help="cNN or cNN<suffix>, e.g. c01 c03_p90")
    ap.add_argument("--sweep", help="one-row layout for a one-factor sweep at one flow; names the figure")
    ap.add_argument("--label", help="configuration named in the title, e.g. '2x2-cell control inlet'")
    ap.add_argument("--no-c26", action="store_true", help="one-row layout: our maps only, no C26 panels or scores")
    ap.add_argument("--out", type=Path,
                    help="default results/fig_c26_fig4_compose[_<sweep>]_<tag>.png")
    a = ap.parse_args(argv)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    p = load_params()
    vmin, vmax = p["fig4"]["v_min_m_s"], p["fig4"]["v_max_m_s"]
    runs = [load_case(s, a.tag, p) for s in a.cases]

    def draw(ax, m, title, ext, ylabel=None):
        ax.set_title(title, fontsize=8)
        ax.set_xlabel("x [mm]")
        if ylabel:
            ax.set_ylabel(ylabel)
        return ax.imshow(m, origin="lower", extent=ext, cmap="jet", vmin=vmin, vmax=vmax)

    if a.sweep:
        flows = {r["flow"] for r in runs}
        if len(flows) > 1 and not a.no_c26:
            raise ValueError("--sweep with C26 panels needs every case at one flow (one Fig. 4 panel per flow)")
        mixed = len({engine(r) for r in runs}) > 1
        ref = next((r for r in runs if scored(r)), runs[0])
        n_c26 = 0 if a.no_c26 else 2
        width = 3.4 * (len(runs) + n_c26) + 1
        fig, axs = plt.subplots(1, len(runs) + n_c26, figsize=(width, 4.0), squeeze=False)
        axs = axs[0]
        if n_c26:
            ref_cond = conditions(ref) if scored(ref) else f"{ref['flow']:g} slpm"
            model_title = "C26 model (Fig. 4b, digitized)"
            if scored(ref, "_c26model_same_pipeline"):
                model_title += f", {ref_cond}\n{scores(ref, '_c26model_same_pipeline')}"
            draw(axs[0], ref["exp"], f"C26 experiment (Fig. 4a, digitized)\n{ref_cond}", ref["ext"])
            draw(axs[1], ref["model"], model_title, ref["ext"])
        for ax, r in zip(axs[n_c26:], runs):
            title = f"this run, {r['case']}: {conditions(r)}"
            if mixed:
                title += f"\n{engine(r)} PIV engine"
            if scored(r) and n_c26:
                title += f"\n{scores(r)}"
            im = draw(ax, r["ours"], title, r["ext"])
        axs[0].set_ylabel("height above map bottom [mm]")
        cbar = (0.925, 0.21 + 0.04 * (not n_c26), 0.008, 0.62)
        adjust = dict(left=0.6 / width, right=0.91, top=0.83, bottom=0.21 + 0.04 * (not n_c26), wspace=0.28)
        suptitle = (f"C26 spout maps, {a.sweep} sweep" + (f" at {ref['flow']:g} slpm" if len(flows) == 1 else "") + ": "
                    + ("this model (synthetic PIV)" if a.no_c26 else "C26 experiment and model, this model"))
    else:
        fig, axs = plt.subplots(len(runs), 3, figsize=(12.6, 3.5 * len(runs) + 1.2), squeeze=False)
        for i, r in enumerate(runs):
            panels = [("C26 experiment (Fig. 4a, digitized)", r["exp"]),
                      (f"C26 model (Fig. 4b, digitized)\n{scores(r, '_c26model_same_pipeline')}", r["model"]),
                      (f"this run ({r['name']})\n{scores(r)}", r["ours"])]
            for j, (title, m) in enumerate(panels):
                im = draw(axs[i, j], m, title, r["ext"],
                          f"{r['case']}, {r['flow']:g} slpm\nheight above map bottom [mm]" if j == 0 else None)
                if i < len(runs) - 1:
                    axs[i, j].set_xlabel("")
        cbar = (0.91, 0.3, 0.012, 0.5)
        adjust = dict(left=0.07, right=0.89, top=0.93, bottom=0.1, hspace=0.35, wspace=0.22)
        suptitle = "C26 spout maps: experiment, C26 model and this model"
    fig.colorbar(im, cax=fig.add_axes(cbar), label="velocity magnitude [m/s]")
    fig.suptitle(suptitle + (f" ({a.label})" if a.label else ""))
    c26 = not (a.sweep and a.no_c26)
    fig.text(0.01, 0.005, textwrap.fill(
             "Time-averaged particle velocity magnitude above the bed. "
             + ("C26 panels digitized from C26 Fig. 4 (p.6). " if c26 else "")
             + f"This run: synthetic PIV, |time mean of the PIV vectors| ({'; '.join(map(window, runs))})."
             + (" Scores: CC / SSI / MSE [m2/s2] against the digitized experiment, same pipeline for both models."
                if c26 else ""),
             int(18 * fig.get_figwidth())),   # ~18 characters per inch at 7 pt
             fontsize=7)
    fig.subplots_adjust(**adjust)
    out = a.out or RESULTS / f"fig_c26_fig4_compose{'_' + a.sweep if a.sweep else ''}_{a.tag}.png"
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()

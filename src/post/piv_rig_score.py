"""Score a PIV rig run against params/piv_rig_acceptance.yaml.

Fountain (dense top, p50/p90/p99 of parcel heights above the bed surface), spout
closure duty d20, parcel inventory and bore entry, inlet health, then the stage's
outcome.  The metric definitions are the C26 review's; --surface-r-min, --spout-vg
and --w-over-a override the rig's values so the script can be checked on C26 runs.
check_smoke.py still covers the solver's finish and dt.

    .venv/bin/python src/post/piv_rig_score.py cases/PIV_rig_3d/r01_smoke05/aws_out --stage smoke_r01
    .venv/bin/python src/post/piv_rig_score.py cases/PIV_rig_3d/r01/aws_out --stage run_r01_3s
    .venv/bin/python src/post/piv_rig_score.py cases/PIV_rig_3d/r02/aws_out --stage run_r02_3s
    .venv/bin/python src/post/piv_rig_score.py cases/PIV_rig_3d/r03/aws_out --stage run_r03_3s
    .venv/bin/python src/post/piv_rig_score.py cases/PIV_rig_3d/q1/aws_out --stage run_q1_3s
    .venv/bin/python src/post/piv_rig_score.py cases/PIV_rig_3d/r01_step20/aws_out --stage step_down_20
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import yaml

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / ".claude/skills/mfix-smoke-validate"))
from check_smoke import cell_centers, read_pvd, read_vtk_xml  # noqa: E402

ACCEPTANCE = REPO / "params/piv_rig_acceptance.yaml"


def frames_in(run, pvd, t0, t1):
    """(t, path) of a .pvd's frames in [t0, t1] that exist on disk.  A restarted
    run's .pvd also lists the source run's frames, which are not in its dir."""
    if not (run / pvd).is_file():
        return []
    out = []
    for t, f in read_pvd(run / pvd):
        if t0 - 1e-6 <= t <= t1 + 1e-6 and (run / f).is_file() and not (out and t <= out[-1][0]):
            out.append((t, run / f))
    return out


def fountain(run, t0, t1, metrics):
    """Parcel heights above the bed surface, pooled over the window's frames."""
    frames = frames_in(run, "PARTICLES_PIV.pvd", t0, t1)[::metrics["parcel_frame_stride"]]
    if not frames:
        return None
    heights, surfaces, n_above = [], [], []
    for _, path in frames:
        xyz = read_vtk_xml(path)["coordinates"]
        r = np.hypot(xyz[:, 0], xyz[:, 2])
        y_s = np.percentile(xyz[r >= metrics["surface_r_min_m"], 1], metrics["surface_percentile"])
        above = xyz[xyz[:, 1] > y_s, 1] - y_s
        heights.append(above)
        surfaces.append(y_s)
        n_above.append(above.size)
    h = np.concatenate(heights)
    bin_m = metrics["dense_bin_m"]
    counts, edges = np.histogram(h, bins=np.arange(0.0, h.max() + 2 * bin_m, bin_m))
    dense_top = edges[np.nonzero(counts >= metrics["dense_peak_frac"] * counts.max())[0][-1] + 1]
    return {
        "n_frames": len(frames),
        "surface_mm": float(np.median(surfaces)) * 1e3,
        "n_above_per_frame": float(np.mean(n_above)),
        "p50_mm": float(np.percentile(h, 50)) * 1e3,
        "p90_mm": float(np.percentile(h, 90)) * 1e3,
        "p99_mm": float(np.percentile(h, 99)) * 1e3,
        "dense_top_mm": float(dense_top) * 1e3,
    }


def spout_shares(run, t0, t1, dx, q, metrics):
    """Spout gas share at share_y_m for each gas frame: the flow through that
    row's cells with v_y above spout_vg_m_s, over the inlet flow."""
    frames = frames_in(run, "BACKGROUND.pvd", t0, t1)
    if not frames:
        return np.array([])
    d0 = read_vtk_xml(frames[0][1])
    pts, conn, offs = d0["coordinates"], d0["connectivity"], d0["offsets"]
    st = np.r_[0, offs[:-1]]
    y = (np.add.reduceat(pts[conn], st, axis=0) / (offs - st)[:, None])[:, 1]
    band = np.abs(y - metrics["share_y_m"]) < 0.5 * dx
    shares = []
    for _, path in frames:
        d = read_vtk_xml(path)
        vy, ep = d["Gas_Velocity"][band, 1], d["EP_G"][band]
        fast = vy > metrics["spout_vg_m_s"]
        shares.append((vy[fast] * ep[fast] * dx * dx).sum() / q)
    return np.array(shares)


def inventory(run, t0, t_after):
    """First and last parcel counts, and the lowest parcel y after t_after."""
    frames = frames_in(run, "PARTICLES.pvd", t0, 1e9)
    if not frames:
        return None
    first = read_vtk_xml(frames[0][1])["coordinates"].shape[0]
    last, y_min = first, np.inf
    for t, path in frames:
        xyz = read_vtk_xml(path)["coordinates"]
        last = xyz.shape[0]
        if t >= t_after:
            y_min = min(y_min, xyz[:, 1].min())
    return {"first": first, "last": last, "y_min_mm": float(y_min) * 1e3}


def source_cells(xyz, dx):
    """Boolean masks of the axis cell in the source row and the one directly above.
    The cut cells beside them have centers at about 0.63 dx, outside 0.5 dx."""
    axis = np.hypot(xyz[:, 0], xyz[:, 2]) < 0.5 * dx
    source = axis & (xyz[:, 1] < 0.0) & (xyz[:, 1] >= -dx)
    above = axis & (xyz[:, 1] >= 0.0) & (xyz[:, 1] < dx)
    if source.sum() != 1 or above.sum() != 1:
        raise SystemExit(f"expected one source and one exit cell, got {source.sum()} and {above.sum()}")
    return source, above


def inlet_health(run, t0, t1, dx, ps_v_g, h_static):
    """Over the window's gas frames: v_y / ps_v_g and EP_G in the cell above the
    source, and P_g in that cell and in the source cell minus the bed-top slab."""
    frames = frames_in(run, "BACKGROUND.pvd", t0, t1)
    if not frames:
        return None
    xyz = cell_centers(read_vtk_xml(frames[0][1]))
    source, above = source_cells(xyz, dx)
    top = (xyz[:, 1] >= h_static) & (xyz[:, 1] < h_static + dx)
    ratio, ep, dp, dp_source = [], [], [], []
    for _, path in frames:
        d = read_vtk_xml(path)
        ratio.append(d["Gas_Velocity"][above, 1][0] / ps_v_g)
        ep.append(d["EP_G"][above][0])
        p_top = d["P_G"][top].mean()
        dp.append(d["P_G"][above][0] - p_top)
        dp_source.append(d["P_G"][source][0] - p_top)
    return {"n_frames": len(frames), "ratio": float(np.mean(ratio)), "ep_g": float(np.mean(ep)),
            "dp_pa": float(np.mean(dp)), "dp_source_pa": float(np.mean(dp_source))}


def health_failures(limits, inlet, inv, w_over_a_pa=None):
    """Names of the failed health checks measured here.  Missing data is N/A, not a pass."""
    if inlet is None or inv is None:
        return ["N/A (missing frames)"]
    wa = w_over_a_pa or limits["w_over_a_pa"]
    lo, hi = limits["exit_ratio_band"]
    plo, phi = limits["pressure_band_w_over_a"]
    failed = []
    if not lo <= inlet["ratio"] <= hi:
        failed.append("inlet_ratio")
    if not plo <= inlet["dp_pa"] / wa <= phi:
        failed.append("pressure")
    if abs(inv["first"] - limits["seeded_parcels"]) > limits["seeded_tolerance"] * limits["seeded_parcels"]:
        failed.append("seeded")
    if (inv["first"] - inv["last"]) / inv["first"] > limits["parcel_loss_max"]:
        failed.append("parcel_loss")
    if inv["y_min_mm"] < 0.0:
        failed.append("no_bore_entry")
    return failed


def outcome(stage, name, fnt, d20, trend):
    """The stage's result word, from the acceptance file's thresholds."""
    if fnt is None or d20 is None:
        return "N/A (missing frames)"
    top = fnt["dense_top_mm"]
    if name == "smoke_r01":
        ok = d20 < stage["regime_d20_max"] and top >= stage["regime_dense_top_min_mm"]
        return "REGIME PASS" if ok else "REGIME FAIL"
    if name.startswith("run_"):
        ok = top >= stage["spouts_dense_top_min_mm"] and d20 < stage["spouts_d20_max"]
        return "SPOUTS" if ok else "NO SPOUT"
    if trend is not None and trend > stage["trend_max"]:
        return "in_between (still changing)"
    if top >= stage["held_dense_top_min_mm"] and d20 < stage["held_d20_max"]:
        return "held"
    if top < stage["collapsed_dense_top_max_mm"]:
        return "collapsed" if d20 >= stage["collapsed_d20_min"] else "internal_jet"
    return "in_between"


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("run", type=Path, help="run output dir (the one holding the .pvd files)")
    p.add_argument("--stage", required=True, choices=("smoke_r01", "run_r01_3s", "run_r02_3s", "run_r03_3s",
                            "run_q1_3s", "run_q2_3s", "run_q3_3s", "step_down_20"))
    p.add_argument("--window", type=float, nargs=2, default=None, metavar=("T0", "T1"),
                   help="fountain window [s] (default: the stage's window_s)")
    p.add_argument("--d20-window", type=float, nargs=2, default=None, metavar=("T0", "T1"),
                   help="d20 window [s] (default: the stage's d20_window_s, else window_s)")
    p.add_argument("--surface-r-min", type=float, default=None,
                   help="override metrics.surface_r_min_m [m], for checks on C26 runs")
    p.add_argument("--spout-vg", type=float, default=None,
                   help="override metrics.spout_vg_m_s [m/s], for checks on C26 runs")
    p.add_argument("--w-over-a", type=float, default=None,
                   help="override health.w_over_a_pa [Pa], for checks on C26 runs")
    args = p.parse_args()

    acc = yaml.safe_load(ACCEPTANCE.read_text())
    metrics, stage = dict(acc["metrics"]), acc[args.stage]
    if args.surface_r_min is not None:
        metrics["surface_r_min_m"] = args.surface_r_min
    if args.spout_vg is not None:
        metrics["spout_vg_m_s"] = args.spout_vg
    manifest = args.run / "manifest.json"
    if not manifest.is_file():
        manifest = args.run.parent / "manifest.json"
    m = json.loads(manifest.read_text())

    t0, t1 = args.window or stage["window_s"]
    d0, d1 = args.d20_window or stage.get("d20_window_s", (t0, t1))
    fnt = fountain(args.run, t0, t1, metrics)
    mid = 0.5 * (t0 + t1)
    halves = [fountain(args.run, a, b, metrics) for a, b in ((t0, mid), (mid, t1))]
    trend = None
    if all(halves) and max(h["dense_top_mm"] for h in halves) > 0:
        tops = [h["dense_top_mm"] for h in halves]
        trend = abs(tops[0] - tops[1]) / max(tops)
    shares = spout_shares(args.run, d0, d1, m["dx_m"], m["Q_actual_m3_s"], metrics)
    d20 = float(np.mean(shares < metrics["closed_share"])) if shares.size else None
    inv = inventory(args.run, 0.0, 0.1)
    h0, h1 = stage["health_window_s"]
    inlet = inlet_health(args.run, h0, h1, m["dx_m"], m["ps_v_g_m_s"], m["H_static_m"])
    wa = args.w_over_a or acc["health"]["w_over_a_pa"]

    print(f"{args.run}: {m.get('flow_slpm')} slpm, stage {args.stage}")
    print(f"  fountain window {t0:g}-{t1:g} s; d20 window {d0:g}-{d1:g} s")
    if fnt:
        print("  surface {surface_mm:.1f} mm; parcels above per frame {n_above_per_frame:.0f}; "
              "p50 {p50_mm:.1f} / p90 {p90_mm:.1f} / p99 {p99_mm:.1f} mm; "
              "dense top {dense_top_mm:.0f} mm ({n_frames} frames)".format(**fnt))
    if trend is not None:
        print(f"  dense top by half: {halves[0]['dense_top_mm']:.0f} / {halves[1]['dense_top_mm']:.0f} mm "
              f"(change {trend:.2f})")
    if shares.size:
        print(f"  spout gas share at {metrics['share_y_m'] * 1e3:g} mm: mean {shares.mean():.2f}, "
              f"min {shares.min():.2f}, max {shares.max():.2f}; "
              f"d20 {d20:.2f} ({int((shares < metrics['closed_share']).sum())}/{shares.size} frames)")
    if inv:
        print(f"  parcels first / last {inv['first']} / {inv['last']} "
              f"({100 * (inv['first'] - inv['last']) / inv['first']:.3f} % lost); "
              f"lowest y after 0.1 s {inv['y_min_mm']:+.2f} mm")
    if inlet:
        print(f"  inlet, {h0:g}-{h1:g} s ({inlet['n_frames']} gas frames): above the source v_g / ps_v_g "
              f"{inlet['ratio']:.2f}, EP_G {inlet['ep_g']:.3f}, minus bed top {inlet['dp_pa']:+.0f} Pa = "
              f"{inlet['dp_pa'] / wa:+.2f} W/A; source cell minus bed top {inlet['dp_source_pa']:+.0f} Pa "
              f"(reported only)")
    failed = health_failures(acc["health"], inlet, inv, wa)
    print(f"  HEALTH: {'PASS' if not failed else 'FAIL ' + ', '.join(failed)} (finish and dt: check_smoke.py)")
    print(f"  OUTCOME: {outcome(stage, args.stage, fnt, d20, trend)}")


if __name__ == "__main__":
    main()

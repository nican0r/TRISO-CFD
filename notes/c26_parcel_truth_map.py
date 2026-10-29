"""Parcel-truth spout map for a C26 run: the time-mean of the parcel velocities
in the Fig. 4 box taken straight from the PARTICLES_PIV frames (no PIV), on
the pipeline's PIV grid and time-mean convention, then on the Fig. 4 grid.
Used to separate "how we map" from "what the parcels do"
(notes/c26_conry2026_case.md, "Map speed excess: mapping or physics", 2026-10-29).

Variants (maps on the Fig. 4 grid, scored like the pipeline):
  all_interp : every parcel in the box (volume average), per-frame bin mean,
               empty bins inpainted per frame (the pipeline's 'interpolate')
  vis_interp : only parcels the camera sees (own >= 1 pixel, nearest wins)
  pix_interp : camera-visible parcels weighted by the pixels they own
  vis_skip   : visible parcels, bins averaged over occupied frames only
  vis_zero   : visible parcels, averaged over all frames (empty = 0)
  magmean    : visible parcels, time mean of |v| (no cancellation)
  cap_X      : vis_interp with parcels faster than X m/s treated as unseen
               (a stand-in for a PIV velocity ceiling)
Also printed: parcel speed percentiles in the box, the axis-column closure
statistics, the fountain height above the annulus surface and the launch
speeds at the spout root.  All values SI; the printout uses mm and m/s.

    MPLBACKEND=Agg .venv/bin/python notes/c26_parcel_truth_map.py \
        cases/C26_3d/<run>/aws_out results/c26_piv_map_<run>.csv [--pct 90] [--tag _p90]
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from src.post import piv_c26 as P  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run", type=Path, help="run dir with manifest.json and PARTICLES_PIV.pvd")
    ap.add_argument("piv_map", type=Path, help="the pipeline's Fig. 4-grid map CSV for the same run")
    ap.add_argument("--pct", type=float, default=None, help="bed-surface percentile (default: params)")
    ap.add_argument("--tag", default="")
    ap.add_argument("--out", type=Path, default=None, help="output dir (default: results/parcel_truth)")
    ap.add_argument("--stride", type=int, default=1)
    a = ap.parse_args(argv)
    out_dir = a.out or REPO / "results" / "parcel_truth"
    out_dir.mkdir(parents=True, exist_ok=True)

    p = P.load_params()
    if a.pct is not None:
        p["spout_map"]["bed_surface_percentile"] = a.pct
    _, read_vtk_xml = P._vtk_readers()
    t0 = time.time()
    g = P.run_geometry(a.run, p, stride=a.stride)
    frames, R_c, px = g["frames"], g["R_c"], g["px"]
    y_s, y_map, H, spout, nrow, ncol, cols = (g[k] for k in ("y_s", "y_map", "H", "spout", "nrow", "ncol", "cols_xyz"))
    run_name = a.run.parent.name if a.run.name == "aws_out" else a.run.name
    print(f"{run_name}: {len(frames)} frames, bed surface {1e3*y_s:.1f} mm, map {1e3*y_map:.1f}+{1e3*H:.1f} mm")

    gx, gy = spout.x, spout.y                      # PIV grid centres [m]
    ny, nx = len(gy), len(gx)
    xe = np.r_[-R_c, 0.5 * (gx[1:] + gx[:-1]), R_c]
    ye = np.r_[y_map, 0.5 * (gy[1:] + gy[:-1]), y_map + H]

    acc = {}

    def new():
        return dict(su=np.zeros((ny, nx)), sv=np.zeros((ny, nx)), n=np.zeros((ny, nx)), nany=0)

    def add(key, su, sv, occ):
        d = acc.setdefault(key, new())
        d["su"] += su; d["sv"] += sv; d["n"] += occ; d["nany"] += 1

    def frame_mean(x, y, u, v, w=None):
        """per-bin weighted mean vector of one frame; occupancy mask"""
        w = np.ones_like(u) if w is None else w
        i = np.clip(np.searchsorted(ye, y, side="right") - 1, 0, ny - 1)
        j = np.clip(np.searchsorted(xe, x, side="right") - 1, 0, nx - 1)
        W = np.zeros((ny, nx)); U = np.zeros((ny, nx)); V = np.zeros((ny, nx))
        np.add.at(W, (i, j), w); np.add.at(U, (i, j), w * u); np.add.at(V, (i, j), w * v)
        with np.errstate(invalid="ignore", divide="ignore"):
            return U / W, V / W, W > 0

    def accumulate(key, mu, mv, occ, conv):
        if not occ.any():
            return
        if conv == "interp":
            add(key, P.inpaint_nans(np.where(occ, mu, np.nan)), P.inpaint_nans(np.where(occ, mv, np.nan)), occ)
        else:
            add(key, np.where(occ, mu, 0.0), np.where(occ, mv, 0.0), occ)

    speeds_all, speeds_vis, axis_vy, axis_n, h99, hmax, launch = [], [], [], [], [], [], []
    caps = (0.3, 0.4, 0.6)
    for k, (t, f) in enumerate(frames):
        d = read_vtk_xml(f)
        xyz, vel, dia = d["coordinates"][:, cols], d["Velocity"][:, cols], d["Diameter"]
        inbox = (np.abs(xyz[:, 0]) <= R_c) & (xyz[:, 1] >= y_map) & (xyz[:, 1] <= y_map + H)
        x, y, u, v = xyz[inbox, 0], xyz[inbox, 1], vel[inbox, 0], vel[inbox, 1]
        sp = np.hypot(u, v)
        _, own = P.render(xyz, dia, (spout.x0, spout.y0), spout.shape, px)
        own = own[own >= 0]
        ids, npix = np.unique(own, return_counts=True)
        vis = np.zeros(len(xyz), bool); vis[ids] = True
        pix = np.zeros(len(xyz)); pix[ids] = npix
        visb, pixb = vis[inbox], pix[inbox]
        speeds_all.append(sp); speeds_vis.append(sp[visb])
        mu, mv, occ = frame_mean(x, y, u, v); accumulate("all_interp", mu, mv, occ, "interp")
        mu, mv, occ = frame_mean(x[visb], y[visb], u[visb], v[visb])
        accumulate("vis_interp", mu, mv, occ, "interp"); accumulate("vis_skip", mu, mv, occ, "skip")
        z = acc.setdefault("vis_zero", new())
        z["su"] += np.where(occ, mu, 0); z["sv"] += np.where(occ, mv, 0); z["n"] += occ; z["nany"] += 1
        mm, _, occm = frame_mean(x[visb], y[visb], sp[visb], np.zeros(visb.sum()))
        accumulate("magmean", mm, np.zeros_like(mm), occm, "skip")
        mu, mv, occ = frame_mean(x[visb], y[visb], u[visb], v[visb], pixb[visb]); accumulate("pix_interp", mu, mv, occ, "interp")
        for c in caps:
            s = visb & (sp <= c)
            mu, mv, occ = frame_mean(x[s], y[s], u[s], v[s]); accumulate(f"cap_{c}", mu, mv, occ, "interp")
        # fountain height above the annulus surface, launch speeds at the spout root
        above = xyz[:, 1] > y_s
        if above.any():
            hh = xyz[above, 1] - y_s
            h99.append(np.percentile(hh, 99)); hmax.append(hh.max())
        root = above & (xyz[:, 1] < y_s + 3e-3) & (np.hypot(xyz[:, 0], xyz[:, 2]) < 10e-3) & (vel[:, 1] > 0)
        launch.append(vel[root, 1])
        # axis column just above the map bottom: closure statistics
        ax = inbox & (np.abs(xyz[:, 0]) < 2.5e-3) & (xyz[:, 1] <= y_map + 10e-3)
        axis_vy.append(vel[ax, 1].mean() if ax.any() else np.nan); axis_n.append(int(ax.sum()))
        if k % 50 == 0:
            print(f"  frame {k}/{len(frames)} t={t:.3f} parcels in box {inbox.sum()} visible {visb.mean():.2f} ({time.time()-t0:.0f} s)", flush=True)

    exp, model = P.load_fig4("exp", int(g["pflow"])), P.load_fig4("model", int(g["pflow"]))
    # the pipeline's map CSV is stored top row first; load_fig4 maps are bottom row first
    piv_map = np.flipud(np.loadtxt(a.piv_map, delimiter=",", comments="#"))
    assert piv_map.shape == exp.shape, (piv_map.shape, exp.shape)

    def fifths(m):  # bottom -> top, on bottom-first maps
        e = np.linspace(0, m.shape[0], 6).astype(int)
        return np.array([m[e[i]:e[i + 1]].mean() for i in range(5)])

    def report(name, m):
        met = P.map_metrics(exp, m, p["metrics"])
        sel = (piv_map > 0.05) | (m > 0.05)
        r = dict(name=name, mean=float(m.mean()), max=float(m.max()),
                 mean_over_exp=float(m.mean() / exp.mean()), mean_over_model=float(m.mean() / model.mean()),
                 fifths_over_exp_bottom_to_top=(fifths(m) / fifths(exp)).round(2).tolist(),
                 CC=met["CC"], SSI=met["SSI"], MSE=met["MSE_m2_s2"], CC_vs_pivmap=P.pearson_cc(piv_map, m),
                 mean_over_pivmap=float(m.mean() / piv_map.mean()),
                 fifths_over_pivmap_bottom_to_top=(fifths(m) / fifths(piv_map)).round(2).tolist(),
                 med_ratio_to_pivmap=float(np.median(m[sel] / np.maximum(piv_map[sel], 1e-6))))
        print(f"{name:12s} mean {r['mean']:.3f} max {r['max']:.3f} | /exp {r['mean_over_exp']:.2f} /C26model "
              f"{r['mean_over_model']:.2f} | CC {r['CC']:.3f} SSI {r['SSI']:.3f} MSE {r['MSE']:.4f} | vs PIV map: CC "
              f"{r['CC_vs_pivmap']:.3f} mean ratio {r['mean_over_pivmap']:.2f} fifths/piv "
              f"{r['fifths_over_pivmap_bottom_to_top']} | fifths/exp (bottom->top) {r['fifths_over_exp_bottom_to_top']}")
        return r

    out = {"run": run_name, "flow": g["flow"], "bed_surface_mm": 1e3 * y_s, "n_frames": len(frames),
           "exp_mean": float(exp.mean()), "model_mean": float(model.mean()), "maps": {}}
    print(f"\nexp mean {exp.mean():.3f} max {exp.max():.3f}; C26 model mean {model.mean():.3f} max {model.max():.3f}")
    out["maps"]["piv_map"] = report("piv_map", piv_map)
    for key, d in acc.items():
        if key in ("vis_skip", "magmean"):
            mean_map = np.hypot(d["su"] / np.maximum(d["n"], 1), d["sv"] / np.maximum(d["n"], 1))
        else:
            mean_map = np.hypot(d["su"] / max(d["nany"], 1), d["sv"] / max(d["nany"], 1))
        m = P.fig4_grid_map(mean_map, spout, y_map, H, R_c, nrow, ncol)
        out["maps"][key] = report(key, m)
        np.savetxt(out_dir / f"{run_name}{a.tag}_{key}.csv", m, delimiter=",", fmt="%.5f",
                   header=f"{run_name}{a.tag} parcel-truth map '{key}' [m/s], Fig. 4 grid, bottom row first")
    occ_frac = acc["vis_interp"]["n"] / acc["vis_interp"]["nany"]
    print(f"\nbin occupancy (visible parcels): mean {occ_frac.mean():.2f}; bins occupied in <10 % of frames {np.mean(occ_frac < 0.1):.2f}")

    sa, sv_ = np.concatenate(speeds_all), np.concatenate(speeds_vis)
    dtf = np.median(np.diff([t for t, _ in frames]))
    print(f"parcels in the box per frame: all {len(sa)/len(frames):.0f}, visible {len(sv_)/len(frames):.0f}")
    for nm, s in (("all", sa), ("visible", sv_)):
        pr = np.percentile(s, [50, 90, 99])
        print(f"  |v| {nm}: p50 {pr[0]:.2f} p90 {pr[1]:.2f} p99 {pr[2]:.2f} m/s; > 0.4 m/s {np.mean(s>0.4):.3f}; "
              f"> 0.8 m/s {np.mean(s>0.8):.3f}   [dt {1e3*dtf:.1f} ms, {1e3*px:.2f} mm/px: 0.4 m/s = {0.4*dtf/px:.0f} px per pair]")
    out["speed"] = {"p50_vis": float(np.percentile(sv_, 50)), "p90_vis": float(np.percentile(sv_, 90)),
                    "frac_gt_0p4_vis": float(np.mean(sv_ > 0.4)), "frac_gt_0p8_vis": float(np.mean(sv_ > 0.8))}
    av = np.array(axis_vy); ok = np.isfinite(av)
    print(f"axis column (|x| < 2.5 mm, first 10 mm above the map bottom): mean parcels {np.mean(axis_n):.1f}; "
          f"net up {np.mean(av[ok]>0):.2f}, net down {np.mean(av[ok]<0):.2f} of frames; mean v_y {np.nanmean(av):+.3f} m/s, "
          f"mean |v_y| {np.nanmean(np.abs(av)):.3f}, |mean| / mean|.| = {abs(np.nanmean(av))/np.nanmean(np.abs(av)):.2f}")
    out["axis"] = {"frac_up": float(np.mean(av[ok] > 0)), "mean_vy": float(np.nanmean(av)), "mean_abs_vy": float(np.nanmean(np.abs(av)))}
    L = np.concatenate(launch)
    print(f"fountain height above the annulus surface: time-median p99 {1e3*np.median(h99):.1f} mm, max {1e3*np.median(hmax):.1f} mm "
          f"(panel height {1e3*H:.1f} mm)")
    print(f"launch v_y at the spout root (r < 10 mm, 0-3 mm above the surface, rising): p50 {np.percentile(L,50):.2f} "
          f"p90 {np.percentile(L,90):.2f} p99 {np.percentile(L,99):.2f} m/s")
    out["fountain"] = {"h99_mm": 1e3 * float(np.median(h99)), "hmax_mm": 1e3 * float(np.median(hmax)),
                       "launch_p50": float(np.percentile(L, 50)), "launch_p90": float(np.percentile(L, 90))}
    (out_dir / f"{run_name}{a.tag}.json").write_text(json.dumps(out, indent=1))
    print(f"done in {time.time()-t0:.0f} s; outputs in {out_dir}")


if __name__ == "__main__":
    main()

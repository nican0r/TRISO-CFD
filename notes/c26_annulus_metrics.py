"""Annulus / spout-regime metrics for a C26 run over its PIV window.

    .venv/bin/python notes/c26_annulus_metrics.py cases/C26_3d/<run>/aws_out [...]

From the BACKGROUND (gas) frames in the manifest's PIV window:
  - spout gas share at y = 10, 20, 30 mm: upward volumetric gas flux through
    the cells with v_g,y > 2 m/s in a one-cell band at that height, over Q;
  - bed mean EP_G between y = 2 mm and the bed surface - 3 mm, and the
    fraction of those cells with EP_G > 0.7 (bubbles / voids);
  - mid-bed interstitial gas v_y and EP_G more than 6 mm from the cone wall.
The bed surface is the time-median over the PARTICLES_PIV frames of the 95th
percentile of parcel y over ALL parcels: the definition src/post/piv_c26.py
used until 2026-10-23 (it now takes the annulus parcels only, r >= 15 mm).
Kept here unchanged so the pre-registered interp_ab baselines stand.

Written 2026-10-23 for the interpolation A/B (c02_momentum24_nointerp
against c02_momentum24); criteria in params/c26_acceptance.yaml interp_ab.
Pure numpy VTK parsing (no pyvista): readers from the smoke-validate skill.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / ".claude/skills/mfix-smoke-validate"))
from check_smoke import read_pvd, read_vtk_xml  # noqa: E402

SPOUT_VG_M_S = 2.0          # a cell belongs to the spout when v_g,y exceeds this
HEIGHTS_MM = (10.0, 20.0, 30.0)


def cell_centers(d):
    pts, conn, offs = d["coordinates"], d["connectivity"], d["offsets"]
    starts = np.r_[0, offs[:-1]]
    return np.add.reduceat(pts[conn], starts, axis=0) / (offs - starts)[:, None]


def frames(run, name, t0, t1):
    return [(t, run / f) for t, f in read_pvd(run / f"{name}.pvd") if t0 <= t <= t1]


def bed_surface_m(run, t0, t1, percentile=95):
    ys = [np.percentile(read_vtk_xml(f)["coordinates"][:, 1], percentile)
          for _, f in frames(run, "PARTICLES_PIV", t0, t1)]
    return float(np.median(ys))


def analyse(run):
    run = Path(run)
    m = json.loads((run / "manifest.json").read_text())
    t0, t1 = float(m["piv_t_start_s"]), float(m["tstop_s"]) + 1e-4
    dx, Q = float(m["dx_m"]), float(m["Q_actual_m3_s"])
    R_t, R_c = float(m["R_t_m"]), float(m["R_c_m"])
    tanh, h_cone = np.tan(np.radians(float(m["cone_half_angle_deg"]))), float(m["h_cone_m"])
    y_s = bed_surface_m(run, t0, t1)
    acc = {}

    def add(k, v):
        acc.setdefault(k, []).append(float(v))

    bg = frames(run, "BACKGROUND", t0, t1)
    for _, f in bg:
        d = read_vtk_xml(f)
        c = cell_centers(d)
        epg, vy = d["EP_G"], d["Gas_Velocity"][:, 1]
        y, r = c[:, 1], np.hypot(c[:, 0], c[:, 2])
        for h in HEIGHTS_MM:
            band = np.abs(y - h * 1e-3) < 0.5 * dx
            q = vy[band] * epg[band] * dx * dx
            add(f"spout_share_{h:g}mm", q[vy[band] > SPOUT_VG_M_S].sum() / Q)
        bed = (y > 2e-3) & (y < y_s - 3e-3)
        add("bed_epg_mean", epg[bed].mean())
        add("bed_void_frac_epg_gt_0.7", (epg[bed] > 0.7).mean())
        rw = np.where(y < h_cone, R_t + y * tanh, R_c)
        mid = (np.abs(y - 0.5 * y_s) < dx) & (rw - r > 6e-3) & (r > 4e-3)
        add("midbed_inner_annulus_vg_y", vy[mid].mean())
        add("midbed_inner_annulus_epg", epg[mid].mean())
    out = {"run": run.parent.name if run.name == "aws_out" else run.name,
           "n_gas_frames": len(bg), "window_s": [t0, float(m["tstop_s"])],
           "bed_surface_mm": 1e3 * y_s}
    out.update({k: float(np.mean(v)) for k, v in acc.items()})
    return out


if __name__ == "__main__":
    for r in sys.argv[1:]:
        res = analyse(r)
        print(json.dumps(res, indent=1))

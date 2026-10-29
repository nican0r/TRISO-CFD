"""Matched-time parcel metrics (notes/c26_conry2026_case.md, "Bore interpolation A/B at matched times", 2026-10-29).
Like-for-like parcel metrics at matched times for the LINEAR_HAT bore smoke,
the NONE bore 3 s run and the control.  Frames: whichever PVD covers the
requested window (PARTICLES_PIV at 2 ms for the smoke, PARTICLES at 50 ms
for the 3 s runs).  Metrics per frame, then the median over the window.
"""
import sys, json
from pathlib import Path
import numpy as np
REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / ".claude/skills/mfix-smoke-validate"))
from check_smoke import read_pvd, read_vtk_xml  # noqa

R_T, R_C, TAN30 = 1.190625e-3, 25.2e-3, np.tan(np.radians(30))

def wall_r(y):
    return np.minimum(R_T + y * TAN30, R_C)

def frames_in(run, t0, t1, pvd):
    out = []
    for t, f in read_pvd(run / pvd):
        if t0 - 1e-6 <= t <= t1 + 1e-6 and not (out and t <= out[-1][0]):
            out.append((t, run / f))
    return out

def metrics(d):
    xyz, vel = d["coordinates"], d["Velocity"]
    x, y, z = xyz.T; r = np.hypot(x, z); vy = vel[:, 1]; sp = np.linalg.norm(vel, axis=1)
    ann = r >= 15e-3
    y_s = np.percentile(y[ann], 95) if ann.any() else np.percentile(y, 95)
    above = y > y_s
    h = y[above] - y_s
    m = dict(y_s=y_s, n_above=int(above.sum()),
             h99=np.percentile(h, 99) if above.any() else 0.0, hmax=h.max() if above.any() else 0.0,
             h50=np.percentile(h, 50) if above.any() else 0.0,
             sp90_above=np.percentile(sp[above], 90) if above.any() else np.nan,
             frac_above_gt04=np.mean(sp[above] > 0.4) if above.any() else np.nan)
    root = above & (y < y_s + 3e-3) & (r < 10e-3) & (vy > 0)
    m["launch_p90"] = np.percentile(vy[root], 90) if root.sum() > 5 else np.nan
    ax = (r < 2.5e-3) & (y > y_s) & (y < y_s + 10e-3)
    m["axis_vy"] = vy[ax].mean() if ax.any() else np.nan
    sp20 = (r < 1.5e-3) & (y > 18e-3) & (y < 22e-3)
    m["axis_vp_20"] = vy[sp20].mean() if sp20.any() else np.nan
    # spout solids up-flux proxy: parcels rising faster than 0.2 m/s at y 18-22 mm, any r
    m["n_rising_20"] = int(((y > 18e-3) & (y < 22e-3) & (vy > 0.2)).sum())
    # wall layer: within 4 mm of the cone/column wall, below the surface, above 10 mm
    wl = (r > wall_r(y) - 4e-3) & (y < y_s) & (y > 10e-3)
    m["wall_vy_down"] = -vy[wl].mean() if wl.any() else np.nan
    # annulus mid-bed packing proxy: mean nearest-neighbour spacing is costly; use parcel count density in a ring
    ring = (r > 12e-3) & (r < 18e-3) & (y > 0.4 * y_s) & (y < 0.6 * y_s)
    m["n_ring"] = int(ring.sum())
    return m

def run_window(run, pvd, t0, t1, stride=1):
    fr = frames_in(run, t0, t1, pvd)[::stride]
    rows = [metrics(read_vtk_xml(f)) for _, f in fr]
    keys = rows[0].keys()
    med = {k: float(np.nanmedian([r[k] for r in rows])) for k in keys}
    med["n_frames"] = len(rows)
    med["axis_net_down_frac"] = float(np.mean([r["axis_vy"] < 0 for r in rows if np.isfinite(r["axis_vy"])]))
    return med, rows, [t for t, _ in fr]

if __name__ == "__main__":
    C = REPO / "cases/C26_3d"
    runs = {
        "LINEAR_HAT bore smoke": (C / "c02_fullface_bore_mom_v232_smoke05/aws_out", "PARTICLES_PIV.pvd"),
        "NONE bore 3 s": (C / "c02_fullface_bore_mom_nointerp_v232/aws_out", "PARTICLES.pvd"),
        "control (2x2 mass-only, LINEAR_HAT)": (C / "c02_apex_v232/aws_out", "PARTICLES.pvd"),
    }
    windows = [(0.30, 0.50), (0.15, 0.30)]
    fmt = "{:36s} {:>9s} {:>6s} {:>6s} {:>6s} {:>6s} {:>7s} {:>7s} {:>7s} {:>7s} {:>7s} {:>7s} {:>6s} {:>6s}"
    print(fmt.format("run", "window", "nfr", "y_s", "h50", "h99", "hmax", "n_abv", "sp90", ">0.4", "launch", "vp20", "netdn", "wall"))
    summary = {}
    for name, (run, pvd) in runs.items():
        for (t0, t1) in windows:
            stride = 5 if pvd == "PARTICLES_PIV.pvd" else 1   # smoke: every 10 ms
            med, rows, ts = run_window(run, pvd, t0, t1, stride)
            summary[f"{name} {t0}-{t1}"] = med
            print(fmt.format(name, f"{t0:.2f}-{t1:.2f}", str(med["n_frames"]), f"{1e3*med['y_s']:.1f}", f"{1e3*med['h50']:.1f}",
                             f"{1e3*med['h99']:.1f}", f"{1e3*med['hmax']:.1f}", f"{med['n_above']:.0f}", f"{med['sp90_above']:.2f}",
                             f"{med['frac_above_gt04']:.2f}", f"{med['launch_p90']:.2f}", f"{med['axis_vp_20']:.2f}",
                             f"{med['axis_net_down_frac']:.2f}", f"{1e3*med['wall_vy_down']:.1f}"))
    # time histories of the fountain height (h99) at 50 ms for the 3 s runs and the smoke
    print("\nfountain h99 [mm] vs time:")
    for name, (run, pvd) in runs.items():
        pv = "PARTICLES.pvd" if (run / "PARTICLES.pvd").exists() else pvd
        fr = frames_in(run, 0.0, 3.0, pv)
        if pv == "PARTICLES_PIV.pvd":
            fr = fr[::25]
        hs = [(t, 1e3 * metrics(read_vtk_xml(f))["h99"]) for t, f in fr]
        print(f"  {name}: " + " ".join(f"{t:.2f}:{h:.0f}" for t, h in hs))
    Path(sys.argv[1] if len(sys.argv) > 1 else "/dev/null").write_text(json.dumps(summary, indent=1))

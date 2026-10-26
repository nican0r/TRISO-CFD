"""Spout-closure cycle of a C26 run: is the annulus gas-supported, and when?

    .venv/bin/python notes/c26_spout_cycle.py cycle cases/C26_3d/<run>/aws_out [t0 t1]
    .venv/bin/python notes/c26_spout_cycle.py core  cases/C26_3d/<run>/aws_out t0 t1
    .venv/bin/python notes/c26_spout_cycle.py plug  cases/C26_3d/<run>/aws_out t_plug t_back
    .venv/bin/python notes/c26_spout_cycle.py jet   cases/C26_3d/<run>/aws_out t0 t1
    .venv/bin/python notes/c26_spout_cycle.py lattice cases/C26_3d/<run>/aws_out t1 t2 ...

cycle: per BACKGROUND frame,
  - spout gas share at y = 10, 20, 30 mm (as notes/c26_annulus_metrics.py:
    upward gas flux through cells with v_g,y > 2 m/s in a one-cell band, / Q);
  - annulus support ratio S = (-dP/dy) / (eps_s (rho_p - rho_g) g), from the
    row means of P_G and EP_G over annulus cells (r > 4.5 mm and more than one
    cell inside the cone wall, so no cut cells).  S ~ 1: the gas carries the
    annulus weight; S > 1: the bed is being lifted; S ~ 0: no gas support.
  Prints the time-mean S(y), the per-frame series, and the closure duty
  d20 = fraction of frames with spout share at 20 mm below 0.3.
core: EP_G and gas v_y of the four axis cells (r < 1.6 mm), row by row.
plug: the parcels in the spout core (r < 2.1 mm, 8 < y < 16 mm) at t_plug,
  traced by Particle_ID to where they were at t_back (PARTICLES_PIV frames).
jet: time-mean, row by row, of the fastest axis cell's gas v_y and the axis
  cells' EP_G (r < 1.6 mm: the centred cell on the 27-cell layout, the four
  axis cells on the 24/26-cell layouts).
lattice: memory of the seeded HEXA lattice (PARTICLES frames).  Each parcel's
  t = 0 neighbours (distance < 1.25 d_parcel); at t, the fraction of those
  pairs still within 1.25 d_parcel (1 = lattice intact, ~0 = re-mixed) and
  the median displacement, by region at time t.

Written 2026-10-24 for the annulus-fluidization review
(notes/c26_conry2026_case.md, "Annulus fluidization is a spout-closure
cycle").  Pure numpy VTK parsing (no pyvista): readers from the
smoke-validate skill.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / ".claude/skills/mfix-smoke-validate"))
from check_smoke import read_pvd, read_vtk_xml  # noqa: E402

G = 9.81
SPOUT_VG_M_S = 2.0          # a cell belongs to the spout when v_g,y exceeds this
CLOSED_SHARE = 0.3          # spout share at 20 mm below this = spout closed
R_ANNULUS_MIN_M = 4.5e-3    # annulus cells start outside the first ring of cells
R_AXIS_MAX_M = 1.6e-3       # the four axis cells (centres at r = 1.48 mm)


def setup(run):
    run = Path(run)
    m = json.loads((run.parent / "manifest.json").read_text())
    bg = read_pvd(run / "BACKGROUND.pvd")
    d0 = read_vtk_xml(run / bg[0][1])
    pts, conn, offs = d0["coordinates"], d0["connectivity"], d0["offsets"]
    st = np.r_[0, offs[:-1]]
    c = np.add.reduceat(pts[conn], st, axis=0) / (offs - st)[:, None]
    y, r = c[:, 1], np.hypot(c[:, 0], c[:, 2])
    tanh = np.tan(np.radians(m["cone_half_angle_deg"]))
    rw = np.where(y < m["h_cone_m"], m["R_t_m"] + y * tanh, m["R_c_m"])
    row = np.floor(y / m["dx_m"] + 1e-9).astype(int)
    return run, m, bg, y, r, rw, row


def cycle(run, t0=0.0, t1=1e9):
    run, m, bg, y, r, rw, row = setup(run)
    dx, Q = m["dx_m"], m["Q_actual_m3_s"]
    drho = m["rho_p_kg_m3"] - m["rho_g_kg_m3"]
    ann = (r > R_ANNULUS_MIN_M) & (rw - r > dx)
    rows = np.arange(1, 22)
    yc = (rows + 0.5) * dx
    pv = read_pvd(run / "PARTICLES_PIV.pvd")
    pt = np.array([t for t, _ in pv])
    rec = []
    for t, f in [(t, f) for t, f in bg if t0 <= t <= t1]:
        d = read_vtk_xml(run / f)
        ep, p, vy = d["EP_G"], d["P_G"], d["Gas_Velocity"][:, 1]
        k = int(np.argmin(np.abs(pt - t)))
        ys = np.nan
        if abs(pt[k] - t) < 0.01:      # bed surface: p95 of annulus parcels (r >= 15 mm)
            q = read_vtk_xml(run / pv[k][1])["coordinates"]
            ys = np.percentile(q[np.hypot(q[:, 0], q[:, 2]) >= 15e-3, 1], 95)
        sel = [ann & (row == j) for j in rows]
        P = np.array([p[s].mean() if s.any() else np.nan for s in sel])
        E = np.array([ep[s].mean() if s.any() else np.nan for s in sel])
        with np.errstate(divide="ignore", invalid="ignore"):
            S = -np.gradient(P, dx) / ((1.0 - E) * drho * G)
        share = {}
        for h in (10, 20, 30):
            band = np.abs(y - h * 1e-3) < 0.5 * dx
            qf = vy[band] * ep[band] * dx * dx
            share[h] = qf[vy[band] > SPOUT_VG_M_S].sum() / Q
        rec.append(dict(t=t, ys=ys, S=S, E=E, share=share))
    ysm = np.nanmedian([q["ys"] for q in rec])
    mid = (yc > 8e-3) & (yc < ysm - 6e-3)
    with np.errstate(invalid="ignore"):
        Sm = np.nanmean([q["S"] for q in rec], axis=0)
        Em = np.nanmean([q["E"] for q in rec], axis=0)
    print(f"{run.parent.name}: {len(rec)} gas frames, t = {rec[0]['t']:.3f}-{rec[-1]['t']:.3f} s, "
          f"bed surface {ysm * 1e3:.1f} mm")
    print(" y[mm]  S_mean  EP_G_annulus")
    for j in range(len(rows)):
        if np.isfinite(Sm[j]) and yc[j] < ysm:
            print(f" {yc[j] * 1e3:5.1f}  {Sm[j]:6.2f}  {Em[j]:6.3f}")
    print("\n t[s]   share10 share20 share30  S_mid  EP_G_annulus_mid")
    s20, Smid, Emid = [], [], []
    for q in rec:
        s20.append(q["share"][20])
        Smid.append(np.nanmean(q["S"][mid]))
        Emid.append(np.nanmean(q["E"][mid]))
        print(f" {q['t']:.3f}   {q['share'][10]:5.2f}  {q['share'][20]:5.2f}  {q['share'][30]:5.2f}"
              f"   {Smid[-1]:5.2f}  {Emid[-1]:6.3f}")
    s20, Smid, Emid = map(np.array, (s20, Smid, Emid))
    closed = s20 < CLOSED_SHARE
    print(f"\nclosure duty d20 = {closed.mean():.2f} ({closed.sum()}/{len(closed)} frames); "
          f"S_mid closed {Smid[closed].mean():.2f} (max {Smid.max():.2f}), "
          f"open {Smid[~closed].mean():.2f}; EP_G_mid closed {Emid[closed].mean():.3f}, "
          f"open {Emid[~closed].mean():.3f}; corr(share20, S_mid) {np.corrcoef(s20, Smid)[0, 1]:+.2f}")
    # Cycle frequency: dominant FFT peak of the mean-removed S_mid series.
    # Only meaningful for gas frames every <= 10 ms (50 ms frames alias the
    # ~10-20 Hz cycle); resolution 1 / window.
    dt = np.diff([q["t"] for q in rec])
    if len(rec) >= 8 and dt.max() <= 0.0101:
        spec = np.abs(np.fft.rfft(Smid - Smid.mean())) ** 2
        freq = np.fft.rfftfreq(len(Smid), d=float(np.median(dt)))
        k = 1 + int(np.argmax(spec[1:]))
        print(f"S_mid dominant frequency {freq[k]:.1f} Hz (resolution {freq[1]:.1f} Hz)")


def core(run, t0, t1):
    run, m, bg, y, r, rw, row = setup(run)
    axis = r < R_AXIS_MAX_M
    rows = np.arange(0, 17)
    print("columns y[mm] =", " ".join(f"{(j + 0.5) * m['dx_m'] * 1e3:4.0f}" for j in rows))
    for t, f in [(t, f) for t, f in bg if t0 <= t <= t1]:
        d = read_vtk_xml(run / f)
        ep, vy = d["EP_G"], d["Gas_Velocity"][:, 1]
        cells = [axis & (row == j) for j in rows]
        print(f"{t:.3f} EP_G", " ".join(f"{ep[s].mean():4.2f}" if s.any() else "  - " for s in cells))
        print(f"      v_y ", " ".join(f"{vy[s].mean():4.1f}" if s.any() else "  - " for s in cells))


def plug(run, t_plug, t_back):
    run, m, *_ = setup(run)
    tanh = np.tan(np.radians(m["cone_half_angle_deg"]))
    pv = read_pvd(run / "PARTICLES_PIV.pvd")
    T = np.array([t for t, _ in pv])

    def frame(t):
        d = read_vtk_xml(run / pv[int(np.argmin(np.abs(T - t)))][1])
        return d["Particle_ID"].astype(int), d["coordinates"], d["Velocity"]

    idp, xp, vp = frame(t_plug)
    idb, xb, vb = frame(t_back)
    rp = np.hypot(xp[:, 0], xp[:, 2])
    inplug = (rp < 2.1e-3) & (xp[:, 1] > 8e-3) & (xp[:, 1] < 16e-3)
    pos = {i: k for k, i in enumerate(idb)}
    k = np.array([pos[i] for i in idp[inplug] if i in pos])
    rb, yb = np.hypot(xb[k, 0], xb[k, 2]), xb[k, 1]
    gap = m["R_t_m"] + yb * tanh - rb
    vr = (xb[k, 0] * vb[k, 0] + xb[k, 2] * vb[k, 2]) / np.maximum(rb, 1e-9)
    print(f"{run.parent.name}: {inplug.sum()} parcels in the spout core at t = {t_plug} s "
          f"(mean v_y {vp[inplug, 1].mean():+.2f} m/s); at t = {t_back} s they were at")
    print(f"  r median {np.median(rb) * 1e3:.1f} mm (10-90 %: {np.percentile(rb, 10) * 1e3:.1f}-"
          f"{np.percentile(rb, 90) * 1e3:.1f}), y median {np.median(yb) * 1e3:.1f} mm "
          f"({np.percentile(yb, 10) * 1e3:.1f}-{np.percentile(yb, 90) * 1e3:.1f})")
    print(f"  already in the core {(rb < 2.1e-3).mean():.0%}; within 4 mm of the wall {(gap < 4e-3).mean():.0%}; "
          f"above 16 mm {(yb > 16e-3).mean():.0%}; v_r median {np.median(vr):+.3f} m/s, "
          f"v_y median {np.median(vb[k, 1]):+.3f} m/s")


def jet(run, t0, t1):
    run, m, bg, y, r, rw, row = setup(run)
    axis = r < R_AXIS_MAX_M
    rows = np.arange(0, 12)
    sel = [axis & (row == j) for j in rows]
    frames_ = [(t, f) for t, f in bg if t0 <= t <= t1]
    V = np.zeros(len(rows))
    E = np.zeros(len(rows))
    for _, f in frames_:
        d = read_vtk_xml(run / f)
        vy, ep = d["Gas_Velocity"][:, 1], d["EP_G"]
        V += [vy[s].max() if s.any() else np.nan for s in sel]
        E += [ep[s].mean() if s.any() else np.nan for s in sel]
    V /= len(frames_)
    E /= len(frames_)
    print(f"{run.parent.name}: {len(frames_)} gas frames, t = {t0}-{t1} s")
    print("y [mm]   " + " ".join(f"{(j + 0.5) * m['dx_m'] * 1e3:5.1f}" for j in rows))
    print("v_y max  " + " ".join(f"{v:5.1f}" for v in V))
    print("EP_G     " + " ".join(f"{e:5.2f}" for e in E))


def lattice(run, *times):
    from scipy.spatial import cKDTree
    run = Path(run)
    m = json.loads((run.parent / "manifest.json").read_text())
    tanh = np.tan(np.radians(m["cone_half_angle_deg"]))
    fr = read_pvd(run / "PARTICLES.pvd")
    T = np.array([t for t, _ in fr])

    def load(t):
        d = read_vtk_xml(run / fr[int(np.argmin(np.abs(T - t)))][1])
        o = np.argsort(d["Particle_ID"])
        return d["Particle_ID"][o].astype(int), d["coordinates"][o], d["Diameter"][o]

    id0, x0, dia = load(0.0)
    cut = 1.25 * float(np.median(dia))
    pairs = np.array(sorted(cKDTree(x0).query_pairs(cut)))
    print(f"{run.parent.name}: {len(id0)} parcels, {2 * len(pairs) / len(id0):.1f} t = 0 neighbours each")

    def regions(x):
        yy, rr = x[:, 1], np.hypot(x[:, 0], x[:, 2])
        gap = m["R_t_m"] + yy * tanh - rr
        return {"wall layer y<20mm": (gap < 2e-3) & (yy < 20e-3),
                "wall layer y>=20mm": (gap < 2e-3) & (yy >= 20e-3),
                "inner annulus": (gap >= 2e-3) & (rr > R_ANNULUS_MIN_M),
                "all": np.ones(len(yy), bool)}

    print(" t [s]  " + "  ".join(f"{k:>28}" for k in regions(x0)))
    for t in times:
        idt, xa, _ = load(t)
        if not np.isin(idt, id0).all():
            raise ValueError(f"parcels appeared between t = 0 and t = {t}")
        # parcels deleted since t = 0 (outlet, or the 23.2 wall loss) drop
        # out, with their pairs
        alive = np.isin(id0, idt)
        xt = np.full_like(x0, np.nan)
        xt[alive] = xa                    # both sorted by Particle_ID
        keep = np.linalg.norm(xt[pairs[:, 0]] - xt[pairs[:, 1]], axis=1) < cut
        disp = np.linalg.norm(xt - x0, axis=1)
        out = []
        for s in regions(xt).values():
            s = s & alive
            pk = s[pairs[:, 0]] & s[pairs[:, 1]]
            out.append(f"keep {keep[pk].mean():.2f} disp {np.median(disp[s]) * 1e3:4.1f} mm")
        print(f" {t:5.2f}  " + "  ".join(f"{c:>28}" for c in out))


if __name__ == "__main__":
    cmd, run, *a = sys.argv[1:]
    {"cycle": cycle, "core": core, "plug": plug, "jet": jet,
     "lattice": lattice}[cmd](run, *map(float, a))

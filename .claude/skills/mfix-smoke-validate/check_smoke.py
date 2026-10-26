#!/usr/bin/env python
"""Automated Phase-2/3 smoke checks for an MFiX spouted-bed run directory.

Implements the quantitative checks added to SKILL.md after the Day-10 R4
failure (orifice clogged with parcels, orifice v_g ~3x bc_v_g, wrong-sign dP,
dt collapse after the bed's first lift-and-collapse at t ~ 30-40 ms):

  S0  run duration       reached t >= --min-tstop
  S1  dt collapse        last-third median dt vs healthy-baseline median dt
  S2  dt stationarity    last-third vs middle-third median dt
  S3  throat v ratio     mean v_g / bc_v_g in [lo, hi], stub-top and exit layers
  S4  throat voidage     mean EP_G > --orifice-epg-min, exit layer
                         (--exit-voidage-info: report only, for a point source
                         at a closed cone apex with no throat)
  S5  stub parcels       parcels below the orifice plane == 0 in the metric window
                         (S5b FLAG if any entered only during start-up;
                         --skip-stub-parcels for inlets without a stub)
  S6  dP sign/magnitude  mean(P_inlet - P_top) > 0 and within band x W/A
  S7  dP noise           std(dP) <= --dp-std-max-fac x W/A        (FLAG only)
  S8  pressure spectrum  fraction of dP power above 100 Hz        (FLAG only)
  S9  run-dir hygiene    no VTK frames on disk that no .pvd in the dir lists

VTK files are parsed with numpy (raw appended, big-endian, UInt64 headers) --
pyvista/VTK hangs in dlopen on the project laptop, so it is not imported.

Defaults are Case V (ORNL/UTK cold mockup, Day-10 R4 deck); override per case.

Usage:
    .venv/bin/python .claude/skills/mfix-smoke-validate/check_smoke.py <run_dir> [options]

Exit status: 0 if no FAIL, 1 otherwise.  FLAG lines do not change exit status
but must be explained in the smoke report.
"""
import argparse
import math
import re
import sys
from pathlib import Path

import numpy as np

G = 9.80665  # m/s^2, standard gravity (matches gravity_y in the Case V deck)


# ----------------------------------------------------------------------------
# Readers
# ----------------------------------------------------------------------------
_DTYPES = {'Float32': '>f4', 'Float64': '>f8', 'Int32': '>i4', 'Int64': '>i8', 'UInt8': 'u1'}
# Block-size prefix of each appended array.  MFiX 26.1.2 writes
# header_type="UInt64"; 23.2 omits the attribute, so the VTK default UInt32.
_HEADER_DTYPES = {'UInt32': '>u4', 'UInt64': '>u8'}


def read_vtk_xml(path):
    """Return {array_name: ndarray} for an MFiX raw-appended .vtp/.vtu file."""
    raw = Path(path).read_bytes()
    i = raw.index(b'<AppendedData encoding="raw">')
    head = raw[:i].decode('latin1')
    base = raw.index(b'_', i) + 1
    m = re.search(r'<VTKFile[^>]*\bheader_type="(\w+)"', head)
    hdt = np.dtype(_HEADER_DTYPES[m.group(1) if m else 'UInt32'])
    out = {}
    pat = (r'<DataArray type="(\w+)" Name="(\w+)"(?:\s+NumberOfComponents="(\d+)")?'
           r'\s+format="appended" offset="\s*(\d+)\s*"')
    for m in re.finditer(pat, head):
        typ, name, nc, off = m.group(1), m.group(2), int(m.group(3) or 1), int(m.group(4))
        p = base + off
        nbytes = int(np.frombuffer(raw[p:p + hdt.itemsize], hdt)[0])
        p += hdt.itemsize
        a = np.frombuffer(raw[p:p + nbytes], _DTYPES[typ]).astype(
            np.float64 if typ.startswith('Float') else np.int64)
        out[name] = a.reshape(-1, nc) if nc > 1 else a
    return out


def cell_centers(d):
    pts, conn, offs = d['coordinates'], d['connectivity'], d['offsets']
    starts = np.r_[0, offs[:-1]]
    return np.add.reduceat(pts[conn], starts, axis=0) / (offs - starts)[:, None]


def read_pvd(path):
    """[(time, filename), ...] from a .pvd collection."""
    if not path.exists():
        return []
    txt = path.read_text(errors='replace')
    return [(float(t), f) for t, f in re.findall(r'timestep="([^"]+)"[^>]*file="([^"]+)"', txt)]


def read_txt_log(path):
    """(t, dt, nit) arrays from MFiX <RUN>.TXT per-step screen log."""
    t, dt, nit = [], [], []
    rx = re.compile(r'\s+([0-9.]+)\s+([0-9.]+E[-+]\d+)\s+\S+\s+(\d+)\s')
    for line in path.open(errors='replace'):
        m = rx.match(line)
        if m:
            t.append(float(m[1])); dt.append(float(m[2])); nit.append(int(m[3]))
    return np.array(t), np.array(dt), np.array(nit)


def read_monitor(path):
    t, v = [], []
    for line in path.open(errors='replace'):
        parts = line.split(',')
        if len(parts) == 2:
            try:
                t.append(float(parts[0])); v.append(float(parts[1]))
            except ValueError:
                pass
    return np.array(t), np.array(v)


def find_one(run_dir, pattern, explicit):
    if explicit:
        return run_dir / explicit
    hits = sorted(run_dir.glob(pattern))
    return hits[0] if hits else None


# ----------------------------------------------------------------------------
# Reporting
# ----------------------------------------------------------------------------
RESULTS = []


def report(code, name, status, detail):
    RESULTS.append(status)
    print(f'  [{status:4s}] {code} {name:22s} {detail}')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('run_dir', type=Path)
    ap.add_argument('--txt', help='per-step solver log (default: the *.TXT with a "DT (s)" header)')
    ap.add_argument('--p-inlet-csv', help='inlet-pressure monitor CSV (default: *P_INLET*.csv)')
    ap.add_argument('--p-top-csv', help='freeboard/top-of-bed pressure monitor CSV (default: *P_TOP*.csv)')
    # --- geometry / operating point (Case V defaults) ---
    ap.add_argument('--bc-v-g', type=float, default=30.0,
                    help='MI superficial gas velocity [m/s]; Case V bc_v_g = 30 (PLAN.md §1.3, S1/S2 30 m/s jet)')
    ap.add_argument('--orifice-y', type=float, default=0.0,
                    help='orifice plane height [m]; Case V stub ends at y = 0 (params/geometry.yaml spouted_bed_V_3d)')
    ap.add_argument('--orifice-r', type=float, default=2.0e-3,
                    help='orifice radius R_i [m]; Case V D_i = 4 mm (S2; INL/CON-07-12569 cross-check)')
    ap.add_argument('--dx', type=float, default=1.0e-3,
                    help='cell size at the orifice [m]; Case V R4 dx_uniform_m = 1 mm (params/geometry.yaml)')
    ap.add_argument('--bed-mass', type=float, default=54.5e-3,
                    help='solids inventory [kg]; Case V 54.5 g ZrO2 (S1/S2)')
    ap.add_argument('--bed-top-r', type=float,
                    default=2.0e-3 + 32.2e-3 * math.tan(math.radians(30.0)),
                    help='bed free-surface radius [m] used for W/A; Case V = R_i + H_static*tan(30 deg) '
                         '(60 deg cone, H_static 32.2 mm, params/geometry.yaml) = 20.6 mm -> W/A ~ 401 Pa')
    # --- time windows ---
    ap.add_argument('--t-response', type=float, default=0.05,
                    help='end of the bed\'s first lift-and-collapse [s]; metrics use t >= this. Case V R4: '
                         'dP spike +20 kPa and dt collapse at 30-40 ms (SKILL.md R4 timeline)')
    ap.add_argument('--t-baseline', type=float, default=0.025,
                    help='end of the healthy pre-response window used as the dt baseline [s]. Case V R4: dt '
                         'was 0.7-2.9e-4 s for t < 30 ms, then collapsed to ~9e-6 s (SKILL.md R4 timeline)')
    ap.add_argument('--min-tstop', type=float, default=0.06,
                    help='minimum simulated time a smoke must reach [s]; must exceed --t-response (SKILL.md §2)')
    # --- thresholds (SKILL.md Phase 2/3) ---
    ap.add_argument('--dt-collapse-max', type=float, default=5.0,
                    help='FAIL if median dt(pre-response) / median dt(last third) exceeds this (SKILL.md §2)')
    ap.add_argument('--dt-stationary-max', type=float, default=2.0,
                    help='FAIL if median dt(middle third) / median dt(last third) exceeds this (SKILL.md §2)')
    ap.add_argument('--v-ratio-lo', type=float, default=0.5, help='orifice v_g/bc_v_g lower bound (SKILL.md §3.1)')
    ap.add_argument('--v-ratio-hi', type=float, default=2.0, help='orifice v_g/bc_v_g upper bound (SKILL.md §3.1)')
    ap.add_argument('--orifice-epg-min', type=float, default=0.9,
                    help='min mean EP_G in orifice cells (SKILL.md §3.6; R4 clogged throat had 0.56-0.62)')
    ap.add_argument('--dp-band', type=float, nargs=2, default=[0.2, 3.0],
                    help='PASS band for mean dP as multiples of W/A (SKILL.md §3.7)')
    ap.add_argument('--dp-std-max-fac', type=float, default=3.0,
                    help='FLAG if std(dP) exceeds this multiple of W/A (SKILL.md §3.7; R4 was 11x)')
    ap.add_argument('--hf-cut', type=float, default=100.0,
                    help='FLAG if >50%% of dP power lies above this frequency [Hz] (SKILL.md §4)')
    ap.add_argument('--skip-stub-parcels', action='store_true',
                    help='skip S5/S5b for inlets with no stub pipe, where y < --orifice-y is '
                         'the bed itself (C26 point source at a closed cone apex: parcels in '
                         'the apex cut cells are not a clog; implementation-steps/case-c26.md)')
    ap.add_argument('--exit-voidage-info', action='store_true',
                    help='report S4 as INFO instead of grading it. For a point source at a '
                         'closed cone apex (C26) there is no throat: the exit layer is the jet '
                         'base, where annulus solids are entrained into the spout, so EP_G < 0.9 '
                         'there is not a clog. S3 (v ratio) remains the clog guard '
                         '(implementation-steps/case-c26.md)')
    a = ap.parse_args()

    rd = a.run_dir
    if not a.t_baseline < a.t_response < a.min_tstop:
        sys.exit('need --t-baseline < --t-response < --min-tstop')
    w_over_a = a.bed_mass * G / (math.pi * a.bed_top_r ** 2)
    print(f'MFiX smoke checks: {rd}')
    print(f'  bc_v_g={a.bc_v_g} m/s  R_i={a.orifice_r*1e3:.2f} mm  orifice y={a.orifice_y} m  '
          f'W/A={w_over_a:.0f} Pa  metric window t >= {a.t_response} s')

    # ---------------- solver log ----------------
    txt = Path(a.txt) if a.txt else None
    if txt is None:
        for p in sorted(rd.glob('*.TXT')):
            with p.open(errors='replace') as fh:
                if 'DT (s)' in fh.read(4000):
                    txt = p
                    break
    if txt is None or not txt.exists():
        report('S0', 'run duration', 'FAIL', 'no per-step *.TXT log found')
        t = dt = nit = np.array([])
    else:
        t, dt, nit = read_txt_log(txt)
    if len(t):
        t_end = t[-1]
        report('S0', 'run duration', 'PASS' if t_end >= a.min_tstop else 'FAIL',
               f'reached t={t_end:.4f} s (min {a.min_tstop} s; must run past first bed response)')
        print('        window          steps  median dt   median NIT  p90 NIT')
        edges = np.linspace(0, t_end, 7)
        for lo, hi in zip(edges[:-1], edges[1:]):
            s = (t >= lo) & (t < hi + 1e-12)
            if s.any():
                print(f'        {lo:6.3f}-{hi:6.3f} s  {s.sum():6d}  {np.median(dt[s]):.2e}   '
                      f'{np.median(nit[s]):5.0f}      {np.percentile(nit[s], 90):5.0f}')
        pre = t < a.t_baseline
        third = t_end / 3
        mid = (t >= third) & (t < 2 * third)
        late = t >= 2 * third
        if t_end <= a.t_response or not pre.any() or not late.any():
            report('S1', 'dt collapse', 'FAIL', 'run too short to compare pre-response vs late dt')
            report('S2', 'dt stationarity', 'FAIL', 'run too short')
        else:
            r1 = np.median(dt[pre]) / np.median(dt[late])
            report('S1', 'dt collapse', 'PASS' if r1 <= a.dt_collapse_max else 'FAIL',
                   f'median dt t<{a.t_baseline}s = {np.median(dt[pre]):.2e}, last third = '
                   f'{np.median(dt[late]):.2e}  (drop {r1:.1f}x; max {a.dt_collapse_max}x)')
            r2 = np.median(dt[mid]) / np.median(dt[late])
            report('S2', 'dt stationarity', 'PASS' if r2 <= a.dt_stationary_max else 'FAIL',
                   f'median dt middle third / last third = {r2:.2f} (max {a.dt_stationary_max})')
    else:
        t_end = 0.0

    # ---------------- VTK frames ----------------
    bg = [(tt, rd / f) for tt, f in read_pvd(rd / 'BACKGROUND.pvd')]
    pa = [(tt, rd / f) for tt, f in read_pvd(rd / 'PARTICLES.pvd')]
    # Every collection in the dir counts as listing its frames: the glob below
    # also matches extra particle regions (e.g. C26's PARTICLES_PIV_*.vtp,
    # listed in PARTICLES_PIV.pvd), which are not stale.
    listed = {f for pvd in rd.glob('*.pvd') for _, f in read_pvd(pvd)}
    on_disk = {p.name for p in list(rd.glob('BACKGROUND_*.vtu')) + list(rd.glob('PARTICLES_*.vtp'))}
    stale = sorted(on_disk - listed)
    report('S9', 'run-dir hygiene', 'PASS' if not stale else 'FLAG',
           'only .pvd-listed frames on disk' if not stale else
           f'{len(stale)} VTK frames NOT in .pvd (dir reused/overwritten: {stale[0]} ...); ignored')

    # Two throat layers (SKILL.md §3.1/§3.6): the last stub cell row below the
    # orifice plane and the first row above it.  In a straight, solids-free stub
    # of radius R_i, continuity forces v_g = bc_v_g there, so any departure is
    # solids lodged in the pipe (R4: parcels in the stub, v_g >> bc_v_g).
    layers = {'stub-top': (a.orifice_y - a.dx, a.orifice_y), 'exit': (a.orifice_y, a.orifice_y + a.dx)}
    vr = {k: [] for k in layers}
    ep = {k: [] for k in layers}
    for tt, f in bg:
        if tt < a.t_response:
            continue
        g = read_vtk_xml(f)
        c = cell_centers(g)
        rr = np.hypot(c[:, 0], c[:, 2])
        for k, (y0, y1) in layers.items():
            m = (rr < a.orifice_r) & (c[:, 1] >= y0) & (c[:, 1] < y1)
            if m.any():
                vr[k].append(g['Gas_Velocity'][m, 1].mean() / a.bc_v_g)
                ep[k].append(g['EP_G'][m].mean())
    for k in layers:
        if not vr[k]:
            st = 'FAIL' if k == 'exit' else 'FLAG'
            report('S3', f'v ratio ({k})', st, f'no cells in layer at t >= {a.t_response} s '
                   '(stub not inside VTK region?)')
            if k == 'exit':
                report('S4', f'voidage ({k})', st, 'no cells')
            continue
        rmean = float(np.mean(vr[k]))
        ok = a.v_ratio_lo <= rmean <= a.v_ratio_hi
        report('S3', f'v ratio ({k})', 'PASS' if ok else 'FAIL',
               f'mean v_g/bc_v_g = {rmean:.2f} over {len(vr[k])} frames (range {min(vr[k]):.2f}-{max(vr[k]):.2f}; '
               f'band [{a.v_ratio_lo}, {a.v_ratio_hi}]; >hi = clogged/choked throat, <lo = blocked BC)')
        if k == 'stub-top':
            continue   # S4 stub-top removed 2026-10-17 (SKILL.md §3.6)
        emean = float(np.mean(ep[k]))
        st = 'INFO' if a.exit_voidage_info else ('PASS' if emean > a.orifice_epg_min else 'FAIL')
        report('S4', f'voidage ({k})', st,
               f'mean EP_G = {emean:.3f} (range {min(ep[k]):.3f}-{max(ep[k]):.3f}; min {a.orifice_epg_min})')

    if a.skip_stub_parcels:
        print('  [SKIP] S5 stub parcels         --skip-stub-parcels: no stub pipe below the inlet')
    else:
        below = []
        for tt, f in pa:
            p = read_vtk_xml(f)
            below.append((tt, int((p['coordinates'][:, 1] < a.orifice_y).sum())))
        # Pass/fail uses the metric window only (ground rule 1); parcels that enter
        # the stub during start-up and are cleared by the jet are a FLAG to explain.
        win = [n for tt, n in below if tt >= a.t_response]
        pre = [n for tt, n in below if tt < a.t_response]
        per_frame = ", ".join(f"{tt:.2f}s:{n}" for tt, n in below)
        if not below:
            report('S5', 'stub parcels', 'FAIL', 'no PARTICLES frames listed in PARTICLES.pvd')
        elif not win:
            report('S5', 'stub parcels', 'FAIL', f'no PARTICLES frames at t >= {a.t_response} s')
        else:
            report('S5', 'stub parcels', 'PASS' if max(win) == 0 else 'FAIL',
                   f'max parcels below orifice plane at t >= {a.t_response} s = {max(win)} over '
                   f'{len(win)} frames (per frame: {per_frame})')
            if pre and max(pre) > 0:
                report('S5b', 'stub parcels startup', 'FLAG',
                       f'max {max(pre)} parcels below orifice plane at t < {a.t_response} s '
                       f'(start-up only; cleared by the window -- explain in the report)')

    # ---------------- pressure monitors ----------------
    pin = find_one(rd, '*P_INLET*.csv', a.p_inlet_csv)
    ptop = find_one(rd, '*P_TOP*.csv', a.p_top_csv)
    if pin is None or ptop is None:
        report('S6', 'dP sign/magnitude', 'FAIL', 'inlet/top pressure monitor CSVs not found')
        return
    ti, pi_ = read_monitor(pin)
    tp, pt = read_monitor(ptop)
    n = min(len(ti), len(tp))
    ti, dp = ti[:n], pi_[:n] - pt[:n]
    s = ti >= a.t_response
    if s.sum() < 10:
        report('S6', 'dP sign/magnitude', 'FAIL', f'fewer than 10 monitor samples after t={a.t_response} s')
        return
    x = dp[s]
    mean, std = x.mean(), x.std()
    lo, hi = a.dp_band
    if mean <= 0:
        st, why = 'FAIL', 'NEGATIVE: inlet below freeboard cannot support the bed -- do not rationalize'
    elif not (lo * w_over_a <= mean <= hi * w_over_a):
        st, why = 'FAIL', f'outside [{lo}, {hi}] x W/A'
    else:
        st, why = 'PASS', 'positive, order W/A'
    report('S6', 'dP sign/magnitude', st,
           f'mean dP = {mean:+.0f} Pa = {mean / w_over_a:+.2f} x W/A over t=[{ti[s][0]:.3f}, {ti[s][-1]:.3f}] s; {why}')
    report('S7', 'dP noise', 'PASS' if std <= a.dp_std_max_fac * w_over_a else 'FLAG',
           f'std dP = {std:.0f} Pa = {std / w_over_a:.1f} x W/A (flag above {a.dp_std_max_fac}x)')
    dtm = np.median(np.diff(ti[s]))
    xu = np.interp(np.arange(ti[s][0], ti[s][-1], dtm), ti[s], x)
    if len(xu) >= 64:
        tt_ = np.arange(len(xu))
        xd = xu - np.polyval(np.polyfit(tt_, xu, 1), tt_)
        P = np.abs(np.fft.rfft(xd)) ** 2
        f = np.fft.rfftfreq(len(xu), dtm)
        band = f >= 1.0
        hf = P[f >= a.hf_cut].sum() / P[band].sum()
        fpk = f[band][np.argmax(P[band])]
        report('S8', 'pressure spectrum', 'PASS' if hf <= 0.5 else 'FLAG',
               f'window {len(xu) * dtm:.3f} s (df={1 / (len(xu) * dtm):.1f} Hz); peak {fpk:.1f} Hz; '
               f'{hf:.0%} of power > {a.hf_cut:.0f} Hz. Not frequency validation (SKILL.md §4).')
    else:
        report('S8', 'pressure spectrum', 'N/A', 'window too short for any spectrum')


if __name__ == '__main__':
    main()
    n_fail = RESULTS.count('FAIL')
    print(f'\nOVERALL: {"NO-GO" if n_fail else "GO (pending manual phases)"} -- '
          f'{n_fail} FAIL, {RESULTS.count("FLAG")} FLAG')
    sys.exit(1 if n_fail else 0)

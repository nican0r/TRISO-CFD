"""C26 step-3 spout checks on a run dir (criteria pre-registered 2026-10-19; notes/c26_conry2026_case.md).

Geometry, static bed height, parcel size and the averaging window are read
from the run dir's manifest.json, so the same checks apply to every case of
params/sweep_C26.yaml.  The C1/C2/D thresholds compare with C26's PIV and
model (YSZ at 10-14 slpm, 60 g: c01-c03); for the coating-stage and charge
cases (c04-c09) the same quantities are reported without a threshold.

    .venv/bin/python notes/c26_step3_checks.py cases/C26_3d/c01/aws_out
"""
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / '.claude/skills/mfix-smoke-validate'))
from check_smoke import read_vtk_xml, read_pvd, read_monitor  # noqa: E402

PIV_CASES = ('c01', 'c02', 'c03')   # C26 Table 3 / Fig. 6-7 comparisons
YS = (10e-3, 20e-3, 30e-3)          # line-probe heights [m] (C26 Fig. 7)
FOUNTAIN_ABOVE = 10e-3              # fountain: parcels above H_static + this [m]
FOUNTAIN_MIN_PARCELS = 20


def load_case(run):
    """Per-case constants from manifest.json (older manifests lack d_parcel_m)."""
    m = json.loads((Path(run) / 'manifest.json').read_text())
    d_cg = m.get('d_parcel_m', m['d_p_m'] * m['cgp_stat_wt'] ** (1 / 3))
    return {
        'case_id': str(m['case_id']),
        't0': float(m['piv_t_start_s']),
        't1': float(m['tstop_s']) + 1e-4,
        'tstop': float(m['tstop_s']),
        'R_T': float(m['R_t_m']), 'R_C': float(m['R_c_m']),
        'TAN': np.tan(np.radians(m['cone_half_angle_deg'])),
        'H_STATIC': float(m['H_static_m']), 'DX': float(m['dx_m']),
        'D_CG': float(d_cg),
    }


def _dp_spectrum(tp, dp, a, b):
    m = (tp >= a) & (tp < b)
    if m.sum() < 8:
        return None
    tu = np.arange(tp[m][0], tp[m][-1], np.median(np.diff(tp[m])))
    s = np.interp(tu, tp[m], dp[m])
    s -= s.mean()
    P = np.abs(np.fft.rfft(s * np.hanning(len(s)))) ** 2
    fr = np.fft.rfftfreq(len(s), tu[1] - tu[0])
    k = 1 + P[1:].argmax()
    return {'window': (a, min(b, tp[m][-1])), 'peak_hz': fr[k], 'df_hz': fr[1],
            'top3_hz': fr[1 + np.argsort(P[1:])[::-1][:3]],
            'mean_pa': dp[m].mean(), 'std_pa': dp[m].std()}


def step3_metrics(run):
    """Step-3 metrics A-D for one run dir, over the manifest's PIV window."""
    run = Path(run)
    c = load_case(run)
    T0, T1, H, DX, D_CG = c['t0'], c['t1'], c['H_STATIC'], c['DX'], c['D_CG']

    def r_wall(y):
        return np.minimum(c['R_T'] + y * c['TAN'], c['R_C'])

    frames = [(t, f) for t, f in read_pvd(run / 'PARTICLES_PIV.pvd') if T0 <= t <= T1]
    fountain, n_frame, wall1, wall2 = [], [], [], []
    xbins = np.arange(-c['R_C'] - DX / 2, c['R_C'] + DX, DX)   # axis-centred bins
    xc = 0.5 * (xbins[1:] + xbins[:-1])
    prof = {y: [np.zeros(len(xc)), np.zeros(len(xc)), np.zeros(len(xc))] for y in YS}
    for t, f in frames:
        d = read_vtk_xml(run / f)
        x, y, z = d['coordinates'].T
        vy = d['Velocity'][:, 1]
        r = np.hypot(x, z)
        fountain.append((y > H + FOUNTAIN_ABOVE).sum())
        n_frame.append(len(y))
        gap = r_wall(y) - r
        for band, out in ((D_CG, wall1), (2 * D_CG, wall2)):
            m = (gap < band) & (y < H)
            out.append(-vy[m].mean() if m.any() else np.nan)
        for yy in YS:
            m = (np.abs(y - yy) < DX / 2) & (np.abs(z) < D_CG)
            idx = np.digitize(x[m], xbins) - 1
            ok = (idx >= 0) & (idx < len(xc))
            s = np.bincount(idx[ok], weights=vy[m][ok], minlength=len(xc))
            n = np.bincount(idx[ok], minlength=len(xc))
            has = n > 0
            mean = np.where(has, s / np.maximum(n, 1), 0)
            prof[yy][0] += mean
            prof[yy][1] += mean ** 2
            prof[yy][2] += has

    fountain = np.array(fountain)
    tp, pin = read_monitor(next(run.glob('*P_INLET.csv')))
    tt, ptop = read_monitor(next(run.glob('*P_TOP_BED.csv')))
    dp = pin - np.interp(tp, tt, ptop)

    probes = {}
    ax = len(xc) // 2
    for yy in YS:
        s, s2, n = prof[yy]
        mean = np.where(n > 0, s / np.maximum(n, 1), np.nan)
        sd = np.sqrt(np.maximum(np.where(n > 0, s2 / np.maximum(n, 1), np.nan) - mean ** 2, 0))
        occ = np.where(n >= 0.5 * max(len(frames), 1))[0]
        probes[yy] = {'r_wall': r_wall(yy), 'axis_v': mean[ax], 'axis_sd': sd[ax],
                      'outer_x': (xc[occ[0]], xc[occ[-1]]) if len(occ) else (np.nan, np.nan),
                      'outer_v': (mean[occ[0]], mean[occ[-1]]) if len(occ) else (np.nan, np.nan),
                      'profile': mean[occ]}
    w1, w2 = np.array(wall1), np.array(wall2)
    return {
        'case': c, 'n_frames': len(frames),
        'frame_t': (frames[0][0], frames[-1][0]) if frames else (np.nan, np.nan),
        'fountain_frac': (fountain >= FOUNTAIN_MIN_PARCELS).mean() if len(fountain) else np.nan,
        'fountain_median': np.median(fountain) if len(fountain) else np.nan,
        'fountain_min': fountain.min() if len(fountain) else np.nan,
        # share of all parcels in the fountain (comparable across parcel counts)
        'fountain_share_median': np.median(fountain / np.array(n_frame)) if len(fountain) else np.nan,
        'B1': _dp_spectrum(tp, dp, T0, T1),
        'B2': _dp_spectrum(tp, dp, 0.5, T1) if T1 > 0.6 else None,
        'probes': probes,
        'wall_dcg': np.nanmean(w1) if len(w1) else np.nan,
        'wall_dcg_sd': np.nanstd(w1) if len(w1) else np.nan,
        'wall_2dcg': np.nanmean(w2) if len(w2) else np.nan,
    }


def main(run):
    res = step3_metrics(run)
    c = res['case']
    print(f"case {c['case_id']}: H_static {c['H_STATIC']*1e3:.1f} mm, parcel d {c['D_CG']*1e3:.3f} mm")
    print(f"PIV frames in [{c['t0']}, {c['tstop']}] s: {res['n_frames']}  "
          f"(t {res['frame_t'][0]:.3f}-{res['frame_t'][1]:.3f})")
    print(f"\nA fountain: frames with >={FOUNTAIN_MIN_PARCELS} parcels above "
          f"{1e3*(c['H_STATIC']+FOUNTAIN_ABOVE):.1f} mm: {res['fountain_frac']*100:.0f}%  "
          f"(median {res['fountain_median']:.0f}, min {res['fountain_min']})")
    for lab in ('B1', 'B2'):
        b = res[lab]
        if b is None:
            print(f'{lab} dP PSD: window too short, skipped')
            continue
        a, e = b['window']
        print(f"{lab} dP PSD t=[{a},{e:.1f}] s: peak {b['peak_hz']:.1f} Hz (df {b['df_hz']:.2f} Hz); "
              f"top-3 {np.round(b['top3_hz'], 1)}; mean dP {b['mean_pa']:.0f} Pa, std {b['std_pa']:.0f} Pa")
    if c['case_id'] not in PIV_CASES:
        print(f"\nC/D thresholds are C26 PIV comparisons for {', '.join(PIV_CASES)}; "
              f"values for {c['case_id']} are informational (no experiment).")
    print('\nC line probe (time-mean v_y of parcels, |z| < d_cg, 2.1 mm x-bins; frames with parcels):')
    for yy, p in res['probes'].items():
        print(f"  y={yy*1e3:.0f} mm  r_wall={p['r_wall']*1e3:.1f} mm  axis v={p['axis_v']:+.3f} "
              f"(sd {p['axis_sd']:.3f}) m/s  outer bins x={p['outer_x'][0]*1e3:+.1f}/"
              f"{p['outer_x'][1]*1e3:+.1f} mm: {p['outer_v'][0]:+.4f}/{p['outer_v'][1]:+.4f} m/s")
        print('     ' + ' '.join(f'{v:+.2f}' for v in p['profile']))
    print(f"\nD wall (downward +): band d_cg: mean {res['wall_dcg']:.2e} m/s (sd over frames "
          f"{res['wall_dcg_sd']:.2e}); band 2 d_cg: {res['wall_2dcg']:.2e} m/s")


if __name__ == '__main__':
    main(Path(sys.argv[1]))

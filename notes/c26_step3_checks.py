"""C26 step-3 spout checks on a baseline run dir (criteria pre-registered 2026-10-19; notes/c26_conry2026_case.md)."""
import sys
from pathlib import Path
import numpy as np

REPO = Path('/Users/nelsonpereira/Documents/GitHub/Personal_Projects/Nuclear/TRISO-CFD')
sys.path.insert(0, str(REPO / '.claude/skills/mfix-smoke-validate'))
from check_smoke import read_vtk_xml, read_pvd, read_monitor  # noqa

RUN = Path(sys.argv[1])
T0, T1 = 2.4, 3.0001
R_T, R_C, TAN = 1.19e-3, 25.2e-3, np.tan(np.radians(30))
H_STATIC, DX = 34.8e-3, 2.1e-3
D_CG = 411.6e-6 * 10 ** (1 / 3)


def r_wall(y):
    return np.minimum(R_T + y * TAN, R_C)


frames = [(t, f) for t, f in read_pvd(RUN / 'PARTICLES_PIV.pvd') if T0 <= t <= T1]
print(f'PIV frames in [{T0}, 3.0] s: {len(frames)}  (t {frames[0][0]:.3f}-{frames[-1][0]:.3f})')

fountain, wall1, wall2 = [], [], []
xbins = np.arange(-R_C - DX / 2, R_C + DX, DX)   # axis-centred bins (mesh cells are centred on the axis)
xc = 0.5 * (xbins[1:] + xbins[:-1])
YS = (10e-3, 20e-3, 30e-3)
prof = {y: [np.zeros(len(xc)), np.zeros(len(xc)), np.zeros(len(xc))] for y in YS}  # sum v, sum v^2, n-frames
for t, f in frames:
    d = read_vtk_xml(RUN / f)
    x, y, z = d['coordinates'].T
    vy = d['Velocity'][:, 1]
    r = np.hypot(x, z)
    fountain.append((y > H_STATIC + 10e-3).sum())
    gap = r_wall(y) - r
    for band, out in ((D_CG, wall1), (2 * D_CG, wall2)):
        m = (gap < band) & (y < H_STATIC)
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
print(f'\nA fountain: frames with >=20 parcels above {1e3*(H_STATIC+10e-3):.1f} mm: '
      f'{(fountain >= 20).mean()*100:.0f}%  (median {np.median(fountain):.0f}, min {fountain.min()})')

tp, pin = read_monitor(next(RUN.glob('*P_INLET.csv')))
tt, ptop = read_monitor(next(RUN.glob('*P_TOP_BED.csv')))
dp = pin - np.interp(tp, tt, ptop)
for a, b, lab in ((2.4, 3.0001, 'B1'), (0.5, 3.0001, 'B2')):
    m = (tp >= a) & (tp < b)
    tu = np.arange(tp[m][0], tp[m][-1], np.median(np.diff(tp[m])))
    s = np.interp(tu, tp[m], dp[m]); s -= s.mean()
    P = np.abs(np.fft.rfft(s * np.hanning(len(s)))) ** 2
    fr = np.fft.rfftfreq(len(s), tu[1] - tu[0])
    k = 1 + P[1:].argmax()
    top3 = fr[1 + np.argsort(P[1:])[::-1][:3]]
    print(f'{lab} dP PSD t=[{a},{min(b,3.0)}] s: peak {fr[k]:.1f} Hz (df {fr[1]:.2f} Hz); top-3 {np.round(top3,1)}; '
          f'mean dP {dp[m].mean():.0f} Pa, std {dp[m].std():.0f} Pa')

print('\nC line probe (time-mean v_y of parcels, |z| < d_cg, 2.1 mm x-bins; frames with parcels):')
for yy in YS:
    s, s2, n = prof[yy]
    mean = np.where(n > 0, s / np.maximum(n, 1), np.nan)
    sd = np.sqrt(np.maximum(np.where(n > 0, s2 / np.maximum(n, 1), np.nan) - mean ** 2, 0))
    occ = np.where(n >= 0.5 * len(frames))[0]
    ax = len(xc) // 2
    print(f'  y={yy*1e3:.0f} mm  r_wall={r_wall(yy)*1e3:.1f} mm  axis v={mean[ax]:+.3f} (sd {sd[ax]:.3f}) m/s  '
          f'outer bins x={xc[occ[0]]*1e3:+.1f}/{xc[occ[-1]]*1e3:+.1f} mm: {mean[occ[0]]:+.4f}/{mean[occ[-1]]:+.4f} m/s')
    print('     ' + ' '.join(f'{v:+.2f}' for v in mean[occ]))

w1, w2 = np.array(wall1), np.array(wall2)
print(f'\nD wall (downward +): band d_cg: mean {np.nanmean(w1):.2e} m/s (sd over frames {np.nanstd(w1):.2e}); '
      f'band 2 d_cg: {np.nanmean(w2):.2e} m/s')

"""Summarize the C26 one-factor-at-a-time sweep (params/sweep_C26.yaml).

For every case with output, runs the C26 step-3 metrics
(notes/c26_step3_checks.py::step3_metrics) over the manifest's PIV window
(the last 0.6 s, t = [2.4, 3.0] s for tstop = 3 s), adds the inventory and
mid-bed centerline v_g monitors, and writes

  results/c26_sweep_summary.csv              one row per case
  results/fig_c26_sweep_<block>.png          one figure per factor block

Blocks are found from the sweep: each case differs from base.ofat_baseline
(c02) in one factor (flow_slpm, coating_stage or inventory_g); the baseline
appears in every block.  Every caption states the averaging window
(CLAUDE.md rule 6).  A case's run dir is <root>/<id><suffix>/aws_out when
that holds a manifest.json, else <root>/<id><suffix>.

    MPLBACKEND=Agg .venv/bin/python notes/c26_sweep_summary.py
    MPLBACKEND=Agg .venv/bin/python notes/c26_sweep_summary.py --suffix _smoke --tag smoke
"""
import argparse
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / 'notes'))
from c26_step3_checks import step3_metrics  # noqa: E402
from check_smoke import read_monitor  # noqa: E402  (path set by c26_step3_checks)
from src.make_case_c26 import load_sweep  # noqa: E402
from src.correlations import G_STD  # noqa: E402
from smoke_validate_c26 import read_monitor3  # noqa: E402

FACTORS = ('flow_slpm', 'coating_stage', 'inventory_g')
BLOCK_X = {'flow_slpm': ('flow', 'gas flow [slpm]'),
           'coating_stage': ('diameter', 'particle diameter d_p [um] (coating stage)'),
           'inventory_g': ('charge', 'seeded bed charge [g] (inventory monitor, t = 0)')}
G = G_STD         # gravity [m/s^2], as the deck's gravity_y = -9.80665
STARTUP_S = 1.0   # windows starting earlier include the startup transient (C26 needs >= 3 s)


def run_dir(root, cid, suffix):
    d = Path(root) / f'{cid}{suffix}'
    return d / 'aws_out' if (d / 'aws_out' / 'manifest.json').is_file() else d


def case_row(cid, run):
    res = step3_metrics(run)
    c = res['case']
    m = json.loads((run / 'manifest.json').read_text())
    inv = read_monitor3(next(run.glob('*SOLIDS_INVENTORY.csv')))
    tv, vg = read_monitor(next(run.glob('*CENTERLINE_VG.csv')))
    win = (tv >= c['t0']) & (tv <= c['t1'])
    b1 = res['B1'] or {}
    # bed weight per bed-surface area, as check_smoke.py's W/A
    r_top = min(m['R_t_m'] + m['H_static_m'] * math.tan(math.radians(m['cone_half_angle_deg'])),
                m['R_c_m'])
    w_over_a = m['inventory_kg'] * G / (math.pi * r_top ** 2)
    return {
        'case': cid, 'run_dir': str(run.relative_to(REPO)) if run.is_relative_to(REPO) else str(run),
        'coating_stage': m.get('coating_stage') or 'bare',
        'd_p_um': m['d_p_m'] * 1e6, 'rho_p_kg_m3': m['rho_p_kg_m3'],
        'charge_g': m['inventory_kg'] * 1e3,
        # MFiX lattice seeding quantizes the charge (c09: one HEXA row short, -4.9 %)
        'seeded_mass_g': inv[0, 1] * 1e3, 'seeded_parcels': inv[0, 2],
        'H_static_mm': m['H_static_m'] * 1e3,
        'flow_slpm': m['flow_slpm'], 'U_nozzle_m_s': m['U_nozzle_m_s'],
        'cgp_stat_wt': m['cgp_stat_wt'], 'd_parcel_mm': c['D_CG'] * 1e3,
        'window_t0_s': c['t0'], 'window_t1_s': c['tstop'], 'n_piv_frames': res['n_frames'],
        'dP_mean_Pa': b1.get('mean_pa', np.nan), 'dP_std_Pa': b1.get('std_pa', np.nan),
        'dP_over_WA': b1.get('mean_pa', np.nan) / w_over_a,
        'dP_psd_peak_Hz': b1.get('peak_hz', np.nan),
        'fountain_frac_frames': res['fountain_frac'], 'fountain_median_parcels': res['fountain_median'],
        'fountain_share_pct': 100.0 * res['fountain_share_median'],
        'axis_vp_10mm_m_s': res['probes'][10e-3]['axis_v'],
        'axis_vp_20mm_m_s': res['probes'][20e-3]['axis_v'],
        'axis_vp_30mm_m_s': res['probes'][30e-3]['axis_v'],
        'wall_vp_down_m_s': res['wall_dcg'],
        'centerline_vg_mid_m_s': vg[win].mean() if win.any() else np.nan,
        'inventory_end_over_start': inv[-1, 1] / inv[0, 1],
    }


def blocks(sweep):
    cases = {c['id']: c for c in sweep['cases']}
    base_id = sweep['base']['ofat_baseline']
    base = cases[base_id]
    out = {f: [base_id] for f in FACTORS}
    for cid, c in cases.items():
        if cid == base_id:
            continue
        diff = [f for f in FACTORS if c.get(f) != base.get(f)]
        if len(diff) != 1:
            raise ValueError(f'{cid} differs from {base_id} in {diff}, not one factor')
        out[diff[0]].append(cid)
    return out


def plot_block(factor, ids, rows, base_id, tag, out_dir):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    have = [rows[i] for i in ids if i in rows]
    if len(have) < 2:
        return None
    xkey = {'flow_slpm': 'flow_slpm', 'coating_stage': 'd_p_um', 'inventory_g': 'seeded_mass_g'}[factor]
    have.sort(key=lambda r: r[xkey])
    x = np.array([r[xkey] for r in have])
    name, xlabel = BLOCK_X[factor]
    fig, axs = plt.subplots(2, 2, figsize=(9, 6.5), sharex=True)
    panels = (('dP_mean_Pa', 'mean dP (inlet - bed top) [Pa]', 'dP_std_Pa'),
              ('dP_psd_peak_Hz', 'dP PSD peak [Hz]', None),
              ('fountain_share_pct', 'parcels in fountain [% of bed, median]', None),
              ('centerline_vg_mid_m_s', 'mid-bed centerline v_g [m/s]', None))
    for ax, (key, ylabel, err) in zip(axs.flat, panels):
        y = np.array([r[key] for r in have], dtype=float)
        e = np.array([r[err] for r in have], dtype=float) if err else None
        ax.errorbar(x, y, yerr=e, marker='o', capsize=3)
        for r, xi, yi in zip(have, x, y):
            lab = r['case'] + (f" {r['coating_stage']}" if factor == 'coating_stage' else '')
            ax.annotate(lab, (xi, yi), textcoords='offset points', xytext=(4, 4), fontsize=7)
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.3)
    for ax in axs[1]:
        ax.set_xlabel(xlabel)
    t0 = min(r['window_t0_s'] for r in have)
    t1 = max(r['window_t1_s'] for r in have)
    smoke = ' SMOKE: window includes the startup transient, not for reporting.' if t0 < STARTUP_S else ''
    fig.suptitle(f'C26 sweep, {name} block (baseline {base_id})')
    fig.text(0.01, 0.005,
             f'Time averages and PSD over t = [{t0:g}, {t1:g}] s (last {t1 - t0:.1f} s, C26 PIV window); '
             f'dP error bars = std over the window.{smoke}', fontsize=7)
    fig.tight_layout(rect=(0, 0.03, 1, 0.96))
    path = out_dir / f'fig_c26_sweep_{name}{"_" + tag if tag else ""}.png'
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--root', type=Path, default=REPO / 'cases' / 'C26_3d')
    ap.add_argument('--suffix', default='', help="run-dir suffix, e.g. '_smoke'")
    ap.add_argument('--tag', default='', help='appended to output file names')
    ap.add_argument('--out-dir', type=Path, default=REPO / 'results')
    a = ap.parse_args()

    sweep = load_sweep()
    rows = {}
    for c in sweep['cases']:
        run = run_dir(a.root, c['id'], a.suffix)
        if not (run / 'PARTICLES_PIV.pvd').is_file():
            print(f"{c['id']}: no output in {run}, skipped")
            continue
        rows[c['id']] = case_row(c['id'], run)
        r = rows[c['id']]
        print(f"{c['id']}: window [{r['window_t0_s']:g}, {r['window_t1_s']:g}] s, dP {r['dP_mean_Pa']:.0f} Pa, "
              f"peak {r['dP_psd_peak_Hz']:.1f} Hz, fountain {r['fountain_share_pct']:.2f}% of parcels, "
              f"v_g mid {r['centerline_vg_mid_m_s']:.2f} m/s")
    if not rows:
        print('no case output found')
        return 1
    a.out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = a.out_dir / f'c26_sweep_summary{"_" + a.tag if a.tag else ""}.csv'
    keys = list(next(iter(rows.values())).keys())
    with csv_path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for cid in sorted(rows):
            w.writerow(rows[cid])
    print(f'wrote {csv_path}')
    for factor, ids in blocks(sweep).items():
        path = plot_block(factor, ids, rows, sweep['base']['ofat_baseline'], a.tag, a.out_dir)
        print(f'wrote {path}' if path else f'{BLOCK_X[factor][0]} block: fewer than 2 cases with output, no figure')
    return 0


if __name__ == '__main__':
    sys.exit(main())

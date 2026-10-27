"""Summarize the C26 one-factor-at-a-time sweep (params/sweep_C26.yaml).

For every case with output, runs the C26 step-3 metrics
(notes/c26_step3_checks.py::step3_metrics) over the manifest's PIV window
(the last 0.6 s, t = [2.4, 3.0] s for tstop = 3 s), adds the inventory and
mid-bed centerline v_g monitors, three one-factor-at-a-time comparison
quantities (added 2026-10-27, "what to compare c04-c09 with c02 on")

  1. solids circulation: parcels crossing the mid static-bed plane upward
     (spout) and downward (annulus) in the PIV window, tracked by
     Particle_ID with a one-cell hysteresis band around the plane, as a
     mass rate [g/s], the bed turnover time [s] and the fraction of the
     bed that passed up through the plane in the window;
  2. spouting margin: U_sup / U_ms (Mathur & Gishler 1955, column
     cross-section, src/correlations.py), U_sup / U_mf (Wen & Yu) and Ar,
     with the case's own d_p, rho_p and H_static;
  3. bed-only pressure drop: P_G averaged over the annulus cells of the
     two rows y in [2 dx, 4 dx] (4.2-8.4 mm) outside the jet (r > 2 dx)
     minus the slab at the static bed top (y in [H_static, H_static + dx]),
     from the BACKGROUND VTU frames in the window, divided by W/A.  The two
     bottom rows are excluded on purpose: the source cell reads the
     source's acceleration pressure (the inlet monitor's dP, kept as
     dP_over_WA, is that probe value, not a bed value) and the row above it
     is inside the 45 m/s jet, where the static pressure is ~500 Pa below
     the annulus.  From row 2 up the jet and annulus pressures agree within
     a few Pa.  The excluded cone tip below 4.2 mm holds ~1 % of the bed;
     a ratio > 1 is a packed bed held by the cone wall (below minimum
     spouting), < 1 a spouting bed whose annulus carries part of the weight;

and writes

  results/c26_sweep_summary.csv              one row per case
  results/fig_c26_sweep_<block>.png          one figure per factor block

Blocks are found from the sweep: each case differs from base.ofat_baseline
(c02) in one factor (flow_slpm, coating_stage or inventory_g); the baseline
appears in every block.  Every caption states the averaging window
(CLAUDE.md rule 6).  A case's run dir is <root>/<id><suffix>/aws_out when
that holds a manifest.json, else <root>/<id><suffix>.

    MPLBACKEND=Agg .venv/bin/python notes/c26_sweep_summary.py --suffix _massonly   # production runs
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
from src.correlations import (G_STD, archimedes, umf_wen_yu,  # noqa: E402
                              ums_mathur_gishler)
from check_smoke import read_vtk_xml, read_pvd, cell_centers  # noqa: E402
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


def plane_transits(frames, y_plane, band):
    """Count parcels transiting a horizontal plane, by Particle_ID, with hysteresis.

    frames : iterable of (ids, y) arrays, in time order
    A parcel is 'below' once y < y_plane - band and 'above' once
    y > y_plane + band; a below->above change is one up-transit, above->below
    one down-transit.  Parcels inside the band keep their last state, so a
    parcel jiggling about the plane (a bubbling bed) is not counted.  Parcels
    first seen inside the band start with no state.  Deleted parcels simply
    stop updating.
    Returns (n_up, n_down, n_ids_up, n_ids_seen): transit counts, number of
    distinct parcels with >= 1 up-transit, number of distinct parcels seen.
    """
    state = {}
    n_up = n_down = 0
    ids_up = set()
    seen = set()
    for ids, y in frames:
        ids = np.asarray(ids).astype(np.int64)
        y = np.asarray(y, dtype=float)
        seen.update(ids.tolist())
        new = np.where(y < y_plane - band, -1, np.where(y > y_plane + band, 1, 0))
        for i, st in zip(ids.tolist(), new.tolist()):
            if st == 0:
                continue
            prev = state.get(i)
            if prev is not None and prev != st:
                if st == 1:
                    n_up += 1
                    ids_up.add(i)
                else:
                    n_down += 1
            state[i] = st
    return n_up, n_down, len(ids_up), len(seen)


def circulation(run, c, m_parcel_kg, bed_mass_kg):
    """Item 1: solids circulation through the mid static-bed plane in the PIV window."""
    T0, T1, DX = c['t0'], c['t1'], c['DX']
    frames = [(t, f) for t, f in read_pvd(run / 'PARTICLES_PIV.pvd') if T0 <= t <= T1]
    if len(frames) < 2:
        return {}

    def gen():
        for _, f in frames:
            d = read_vtk_xml(run / f)
            yield d['Particle_ID'], d['coordinates'][:, 1]
    y_plane = 0.5 * c['H_STATIC']
    n_up, n_down, n_ids_up, n_seen = plane_transits(gen(), y_plane, DX / 2)
    dt = frames[-1][0] - frames[0][0]
    g_s = m_parcel_kg * n_up / dt                   # spout solids mass rate [kg/s]
    return {
        'circ_plane_mm': y_plane * 1e3, 'circ_window_s': dt,
        'circ_up_transits': n_up, 'circ_down_transits': n_down,
        'circ_rate_g_s': g_s * 1e3,
        'circ_closure_down_over_up': n_down / n_up if n_up else np.nan,
        'turnover_time_s': bed_mass_kg / g_s if g_s > 0 else np.inf,
        'bed_frac_cycled_window': n_ids_up / n_seen if n_seen else np.nan,
    }


def spouting_margin(m):
    """Item 2: U_sup against U_ms (Mathur-Gishler), U_mf (Wen-Yu) and Ar, from the manifest."""
    d_p, rho_p = m['d_p_m'], m['rho_p_kg_m3']
    rho, mu = m['rho_g_kg_m3'], m['mu_g_pa_s']
    D_c, D_i, H = m['D_c_m'], 2.0 * m['R_t_m'], m['H_static_m']
    u_ms = ums_mathur_gishler(d_p, rho_p, rho, D_c, D_i, H)   # on the column cross-section
    u_mf = umf_wen_yu(d_p, rho_p, rho, mu)
    return {
        'U_sup_m_s': m['U_sup_m_s'],
        'U_ms_mg_m_s': u_ms, 'U_over_U_ms': m['U_sup_m_s'] / u_ms,
        'U_mf_wenyu_m_s': u_mf, 'U_over_U_mf': m['U_sup_m_s'] / u_mf,
        'Ar': archimedes(d_p, rho_p, rho, mu),
    }


def bed_dp(run, c, w_over_a):
    """Item 3: bed-only dP from the VTU gas field: annulus at y = 2-4 dx minus the bed-top slab."""
    T0, T1, DX, H = c['t0'], c['t1'], c['DX'], c['H_STATIC']
    frames = [(t, f) for t, f in read_pvd(run / 'BACKGROUND.pvd') if T0 <= t <= T1]
    vals = []
    for _, f in frames:
        d = read_vtk_xml(run / f)
        x, y, z = cell_centers(d).T
        lo = (y > 2 * DX) & (y < 4 * DX) & (np.hypot(x, z) > 2 * DX)
        hi = (y > H) & (y < H + DX)
        if lo.any() and hi.any():
            vals.append(d['P_G'][lo].mean() - d['P_G'][hi].mean())
    if not vals:
        return {}
    v = np.array(vals)
    # 13 frames at 50 ms under-sample the ~12 Hz spout pulsation (aliased): the
    # mean is unbiased but its standard error is std/sqrt(n), ~+-170 Pa on c02.
    return {'dP_bed_n_frames': len(v), 'dP_bed_mean_Pa': v.mean(), 'dP_bed_std_Pa': v.std(),
            'dP_bed_sem_Pa': v.std() / np.sqrt(len(v)),
            'dP_bed_over_WA': v.mean() / w_over_a, 'dP_bed_over_WA_sem': v.std() / np.sqrt(len(v)) / w_over_a}


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
        'W_over_A_Pa': w_over_a,
        **spouting_margin(m),
        **bed_dp(run, c, w_over_a),
        **circulation(run, c, inv[0, 1] / inv[0, 2], inv[0, 1]),
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
    fig, axs = plt.subplots(3, 3, figsize=(13, 10), sharex=True)
    panels = (('dP_mean_Pa', 'dP, inlet monitor - bed top [Pa]', 'dP_std_Pa'),
              ('dP_bed_over_WA', 'bed dP / (W/A) [-]\n(annulus y = 4-8 mm - bed top)', 'dP_bed_over_WA_sem'),
              ('dP_psd_peak_Hz', 'dP PSD peak [Hz]', None),
              ('fountain_share_pct', 'fountain [% of parcels, median]', None),
              ('centerline_vg_mid_m_s', 'mid-bed axis v_g [m/s]', None),
              ('U_over_U_ms', 'U_sup / U_ms [-]\n(Mathur-Gishler, column area)', None),
              ('circ_rate_g_s', 'solids up through y = H/2 [g/s]', None),
              ('turnover_time_s', 'bed turnover time [s]', None),
              ('bed_frac_cycled_window', 'bed fraction up through y = H/2\nin the window [-]', None))
    for ax, (key, ylabel, err) in zip(axs.flat, panels):
        y = np.array([r.get(key, np.nan) for r in have], dtype=float)
        e = np.array([r[err] for r in have], dtype=float) if err else None
        if key == 'turnover_time_s':
            y = np.where(np.isfinite(y), y, np.nan)
        ax.errorbar(x, y, yerr=e, marker='o', capsize=3)
        for r, xi, yi in zip(have, x, y):
            lab = r['case'] + (f" {r['coating_stage']}" if factor == 'coating_stage' else '')
            ax.annotate(lab, (xi, yi), textcoords='offset points', xytext=(4, 4), fontsize=7)
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.3)
    for ax in axs[-1]:
        ax.set_xlabel(xlabel)
    t0 = min(r['window_t0_s'] for r in have)
    t1 = max(r['window_t1_s'] for r in have)
    smoke = ' SMOKE: window includes the startup transient, not for reporting.' if t0 < STARTUP_S else ''
    fig.suptitle(f'C26 sweep, {name} block (baseline {base_id})')
    fig.text(0.01, 0.005,
             f'Time averages and PSD over t = [{t0:g}, {t1:g}] s (last {t1 - t0:.1f} s, C26 PIV window); '
             f'dP error bars = std over the window; bed dP error bars = standard error of the 13 VTU frames (50 ms).\n'
             f'Circulation: Particle_ID transits of the plane y = H_static/2 with a dx/2 hysteresis band, 2 ms frames; '
             f'U_ms on the column cross-section.{smoke}', fontsize=7)
    fig.tight_layout(rect=(0, 0.04, 1, 0.96))
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
              f"v_g mid {r['centerline_vg_mid_m_s']:.2f} m/s, U/U_ms {r['U_over_U_ms']:.2f}, "
              f"bed dP/(W/A) {r.get('dP_bed_over_WA', np.nan):.2f}, "
              f"circulation {r.get('circ_rate_g_s', np.nan):.2f} g/s, "
              f"turnover {r.get('turnover_time_s', np.nan):.1f} s")
    if not rows:
        print('no case output found')
        return 1
    a.out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = a.out_dir / f'c26_sweep_summary{"_" + a.tag if a.tag else ""}.csv'
    keys = []
    for r in rows.values():
        keys += [k for k in r if k not in keys]
    with csv_path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys, restval='')
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

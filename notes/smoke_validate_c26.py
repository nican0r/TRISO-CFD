"""C26 smoke validation for any case of params/sweep_C26.yaml.

Runs the mfix-smoke-validate checker with arguments built from the run dir's
manifest.json (operating point, geometry, bed charge), the way the c01 smokes
were checked (results/run_log.csv 2026-10-17/19: --skip-stub-parcels
--exit-voidage-info), then the C26 step-2 checks of
implementation-steps/case-c26.md step 2:

  P1  point-source cell volume (run.log) is the axis cell, not a sliver;
  P2  seeded parcels within 10 % of the manifest's n_parcels_est; parcels
      leave only through the top pressure outlet (tracked by Particle_ID),
      at most 0.5 % of the seeded mass (amended 2026-10-19, user: the c02
      smoke lost 3 ballistic fountain parcels via the PO at y = 210 mm);
  P3  axial spout channel, PIV frames t = 0.1 s to tstop, slab
      y = 25.2 +/- 1.05 mm: (i) time-mean upward solids mass flux through
      r < 3.15 mm >= 5 m_dot_g; (ii) the centroid of the upward solids
      mass flux lies within 1.05 mm of the axis;
  P4  time-mean centerline v_g at mid-bed (Monitor 4), t >= 0.1 s, >= 5 m/s;
  P5  settled bed height (95th percentile of parcel y, PIV frames
      t = 0.04-0.08 s) reported next to H_static (INFO).

    .venv/bin/python notes/smoke_validate_c26.py cases/C26_3d/c07_smoke

P1-P2 (with check_smoke) gate inlet health and numerics.  P3-P4 say whether
the bed spouts; a spout is required only for the C26 PIV cases (manifest
piv_targets_table3 set: c01-c03).  For the coating-stage and charge cases at
fixed 12 slpm, no spout is reported as a result (user decision 2026-10-19).
Mass-only point source (manifest ps_momentum false; the geometry.yaml default
since 2026-10-19, user decision): check_smoke S3 (exit v_g / U_nozzle band) and P4 (5 m/s) were set
on momentum-source runs, so they are reported, not gated; the spout is judged
on P3.
Prints OVERALL: GO / NO-GO; exits 1 on NO-GO.
"""
import json
import math
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
SKILL = REPO / '.claude/skills/mfix-smoke-validate'
sys.path.insert(0, str(SKILL))
from check_smoke import read_vtk_xml, read_pvd, read_monitor  # noqa: E402

# Thresholds: implementation-steps/case-c26.md step 2 (pre-registered 2026-10-17).
T_RESPONSE = 0.1             # start of the evaluated window [s]
PS_VOLUME_MIN_M3 = 5.0e-9    # axis cell (row 0) 7.0e-9; corner sliver 3.3e-10
SEED_SHORTFALL_MAX = 0.10    # recipe Rule 8
OUTLET_LOSS_MAX = 0.005      # mass fraction that may leave via the top PO (user, 2026-10-19)
SLAB_Y = 25.2e-3             # spout-flux slab centre [m] (12 dx)
CORE_R = 3.15e-3             # spout core radius [m] (1.5 dx)
FLUX_MIN_MDOT = 5.0          # (i) upward solids flux >= 5 m_dot_g
CENTROID_MAX = 1.05e-3       # (ii) half a cell [m]
VG_MID_MIN = 5.0             # P4 [m/s]
SETTLE_WINDOW = (0.04, 0.08)  # [s] after the initial collapse (r2 log)


def read_monitor3(path):
    """Type-101 inventory monitor: rows 'time, pmass, pcount'."""
    rows = []
    for line in Path(path).open(errors='replace'):
        parts = line.split(',')
        if len(parts) == 3:
            try:
                rows.append([float(p) for p in parts])
            except ValueError:
                pass
    return np.array(rows)


RESULTS = []   # (group, status); group 'inlet' (P1-P2) or 'spout' (P3-P4)


def report(code, name, status, detail, group='inlet'):
    RESULTS.append((group, status))
    print(f'  [{status:4s}] {code} {name:22s} {detail}')


def lost_parcel_exits(run, m):
    """(lost, via_outlet): parcels missing between consecutive PIV frames, and
    how many of them reach within dx/2 of the pressure outlet at y = H_dom
    when their last position is carried on at their last v_y for one frame
    interval (ballistic fountain parcels, c02 smoke 2026-10-19)."""
    lost, top, prev = 0, 0, None
    for t, f in read_pvd(run / 'PARTICLES_PIV.pvd'):
        d = read_vtk_xml(run / f)
        ids = d['Particle_ID'].astype(np.int64)
        if prev is not None:
            pid, py, pvy, pt = prev
            gone = ~np.isin(pid, ids)
            reach = py[gone] + np.maximum(pvy[gone], 0.0) * (t - pt)
            lost += int(gone.sum())
            top += int((reach >= m['H_dom_m'] - m['dx_m'] / 2).sum())
        prev = (ids, d['coordinates'][:, 1], d['Velocity'][:, 1], t)
    return lost, top


def run_check_smoke(run, m):
    bed_top_r = min(m['R_t_m'] + m['H_static_m'] * math.tan(math.radians(m['cone_half_angle_deg'])),
                    m['R_c_m'])
    cmd = [sys.executable, str(SKILL / 'check_smoke.py'), str(run),
           '--bc-v-g', f"{m['U_nozzle_m_s']:.6g}",
           '--orifice-r', f"{m['R_nozzle_m']:.6g}",
           '--dx', f"{m['dx_m']:.6g}",
           '--bed-mass', f"{m['inventory_kg']:.6g}",
           '--bed-top-r', f'{bed_top_r:.6g}',
           '--t-response', f'{T_RESPONSE:g}',
           '--min-tstop', f"{0.999 * m['tstop_s']:.6g}",
           '--skip-stub-parcels', '--exit-voidage-info']
    print('$ ' + ' '.join(cmd[1:]))
    proc = subprocess.run(cmd, capture_output=True, text=True)
    print(proc.stdout.rstrip())
    if proc.stderr.strip():
        print(proc.stderr.rstrip())
    if proc.returncode == 0:
        return True
    # Mass-only point source (ps_v_g = 0): S3's exit v_g / bc_v_g band assumes
    # a momentum source (c01_smoke_r4_massonly: 0.39/0.24), so S3 is reported,
    # not gated (mass-only is the C26 default since 2026-10-19, user decision).
    fails = [ln.split(']', 1)[1].split()[0] for ln in proc.stdout.splitlines()
             if ln.strip().startswith('[FAIL]')]
    if m.get('ps_momentum') is False and fails and all(f == 'S3' for f in fails):
        print('  (mass-only source: S3 v ratio reported, not gated)')
        return True
    return False


def main(run):
    run = Path(run)
    m = json.loads((run / 'manifest.json').read_text())
    print(f"== C26 smoke validation: {m['case_id']} ({m['feedstock_label']}), "
          f"{m['flow_slpm']:g} slpm, {m['inventory_kg']*1e3:g} g, run dir {run}")
    cs_ok = run_check_smoke(run, m)

    print('\nC26 step-2 checks (implementation-steps/case-c26.md step 2):')
    # P1 point-source cell volume
    log = (run / 'run.log').read_text(errors='replace') if (run / 'run.log').is_file() else ''
    mv = re.search(r'Point source setup information:\s*PS\s+Volume\s+1\s+([0-9.E+\-]+)', log)
    if mv:
        vol = float(mv.group(1))
        report('P1', 'PS cell volume', 'PASS' if vol >= PS_VOLUME_MIN_M3 else 'FAIL',
               f'{vol:.3e} m^3 (axis cell 7.0e-9; sliver ~3.3e-10)')
    else:
        report('P1', 'PS cell volume', 'FAIL', 'no "Point source setup information" in run.log')

    # P2 seeding and inventory
    inv_csv = next(run.glob('*SOLIDS_INVENTORY.csv'), None)
    inv = read_monitor3(inv_csv) if inv_csv else np.empty((0, 3))
    if len(inv):
        n0, n1, m0, m1 = inv[0, 2], inv[-1, 2], inv[0, 1], inv[-1, 1]
        short = 1.0 - n0 / m['n_parcels_est']
        lost, top = lost_parcel_exits(run, m) if n1 < n0 else (0, 0)
        loss = 1.0 - m1 / m0
        ok = (short <= SEED_SHORTFALL_MAX and lost == n0 - n1 == top
              and loss <= OUTLET_LOSS_MAX)
        report('P2', 'seeding / inventory', 'PASS' if ok else 'FAIL',
               f"seeded {n0:.0f} vs est {m['n_parcels_est']:.0f} ({-short*100:+.1f}%), "
               f"{m0*1e3:.2f} g vs charge {m['inventory_kg']*1e3:g} g; end {n1:.0f} parcels, "
               f"{m1*1e3:.2f} g at t={inv[-1, 0]:.3f} s; {n0 - n1:.0f} lost, {top} via the "
               f"top outlet ({loss*100:.3f}% of seeded mass; <= {OUTLET_LOSS_MAX*100:g}%, "
               f"none elsewhere)")
    else:
        report('P2', 'seeding / inventory', 'FAIL', 'no *SOLIDS_INVENTORY.csv monitor data')

    # P3 spout channel from PIV frames
    t_end = m['tstop_s'] + 1e-4
    frames = [(t, f) for t, f in read_pvd(run / 'PARTICLES_PIV.pvd') if T_RESPONSE <= t <= t_end]
    rho_p, dx = m['rho_p_kg_m3'], m['dx_m']
    flux_core, flux_r = [], []
    for _, f in frames:
        d = read_vtk_xml(run / f)
        x, y, z = d['coordinates'].T
        vy = d['Velocity'][:, 1]
        mp = rho_p * math.pi / 6.0 * d['Diameter'] ** 3   # VTK Diameter is the parcel size
        r = np.hypot(x, z)
        sl = (np.abs(y - SLAB_Y) < dx / 2) & (vy > 0)
        q = mp[sl] * vy[sl] / dx                        # upward mass flux through the slab [kg/s]
        flux_core.append(q[r[sl] < CORE_R].sum())
        flux_r.append((q.sum(), (q * x[sl]).sum(), (q * z[sl]).sum()))
    if frames:
        fc = np.mean(flux_core) / m['ps_massflow_g_kg_s']
        qs = np.array(flux_r)
        qsum = qs[:, 0].sum()
        centroid = (math.hypot(qs[:, 1].sum() / qsum, qs[:, 2].sum() / qsum)
                    if qsum > 0 else np.nan)   # distance of the flux centroid from the axis
        report('P3i', 'spout core flux', 'PASS' if fc >= FLUX_MIN_MDOT else 'FAIL',
               f'{fc:.1f} m_dot_g through r < {CORE_R*1e3:.2f} mm at y = {SLAB_Y*1e3:.1f} mm '
               f'({len(frames)} frames, t >= {T_RESPONSE} s; >= {FLUX_MIN_MDOT:g})', group='spout')
        report('P3ii', 'spout flux centroid', 'PASS' if centroid <= CENTROID_MAX else 'FAIL',
               f'{centroid*1e3:.2f} mm (<= {CENTROID_MAX*1e3:.2f})', group='spout')
    else:
        report('P3', 'spout channel', 'FAIL', f'no PIV frames in [{T_RESPONSE}, {m["tstop_s"]}] s', group='spout')

    # P4 mid-bed centerline v_g
    vg_csv = next(run.glob('*CENTERLINE_VG.csv'), None)
    tv, vg = read_monitor(vg_csv) if vg_csv else (np.array([]), np.array([]))
    sel = tv >= T_RESPONSE
    vmid = vg[sel].mean() if sel.any() else np.nan
    # The 5 m/s threshold was set on momentum-source runs; a mass-only source
    # spreads the jet (c01_smoke_r4_massonly: 3.2 m/s with a P3 spout), so
    # for ps_v_g = 0 P4 is reported and the spout is judged on P3.
    p4 = 'PASS' if vmid >= VG_MID_MIN else 'FAIL'
    if m.get('ps_momentum') is False:
        p4 = 'INFO'
    report('P4', 'mid-bed centerline v_g', p4,
           f"{vmid:.2f} m/s at y = {m['H_mid_m']*1e3:.1f} mm (>= {VG_MID_MIN:g}"
           + ('; momentum-source threshold, reported only)' if p4 == 'INFO' else ')'),
           group='spout')

    # P5 settled bed height (report only)
    early = [(t, f) for t, f in read_pvd(run / 'PARTICLES_PIV.pvd')
             if SETTLE_WINDOW[0] <= t <= SETTLE_WINDOW[1]]
    if early:
        y95 = [np.percentile(read_vtk_xml(run / f)['coordinates'][:, 1], 95) for _, f in early]
        report('P5', 'settled bed height', 'INFO',
               f"y95 {min(y95)*1e3:.1f}-{max(y95)*1e3:.1f} mm at t = {SETTLE_WINDOW[0]}-"
               f"{SETTLE_WINDOW[1]} s vs H_static {m['H_static_m']*1e3:.1f} mm")
    else:
        report('P5', 'settled bed height', 'INFO', 'no PIV frames in the settling window')

    # A spout is required only where C26 observed one (PIV cases c01-c03).
    # Elsewhere (coating stages, charges at fixed 12 slpm) no spout is a
    # result to report, not a failure (user decision 2026-10-19).
    n_fail = {g: sum(1 for gg, s in RESULTS if gg == g and s == 'FAIL') for g in ('inlet', 'spout')}
    inlet_ok = cs_ok and n_fail['inlet'] == 0
    spout = n_fail['spout'] == 0
    need_spout = m.get('piv_targets_table3') is not None
    go = inlet_ok and (spout or not need_spout)
    spout_txt = ('spout yes' if spout else
                 f"NO SPOUT ({n_fail['spout']} spout check(s) failed; "
                 + ('required: C26 PIV case)' if need_spout else 'reported as a result, not required)'))
    print(f"\nOVERALL: {'GO' if go else 'NO-GO'} (inlet/numerics {'GO' if inlet_ok else 'NO-GO'}: "
          f"check_smoke {'GO' if cs_ok else 'NO-GO'}, {n_fail['inlet']} P1-P2 FAIL; {spout_txt})")
    return 0 if go else 1


if __name__ == '__main__':
    sys.exit(main(sys.argv[1]))

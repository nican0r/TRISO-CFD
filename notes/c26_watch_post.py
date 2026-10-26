"""Watch C26 AWS Batch jobs, pull each on success, post-process it and score
it against Conry et al. 2026 (C26) and against the same case on the old grid
layout.

    MPLBACKEND=Agg .venv/bin/python notes/c26_watch_post.py \\
        c01_massonly_axisface c02_massonly_axisface c03_massonly_axisface --tag axisface

Per case (a run-dir name under cases/C26_3d):
  1. job id: the last aws-batch row for the case in results/run_log.csv
     (written by src/aws/submit.py);
  2. poll AWS Batch every --interval s until SUCCEEDED or FAILED;
  3. on SUCCEEDED: s3 sync the job's RESULTS_S3 prefix into
     cases/C26_3d/<case>/aws_out (as src/aws/poll.py pull_results) and check
     run.log reached tstop;
  4. post-process; each step's output goes to <run>/post/<step>.log:
       health  notes/smoke_validate_c26.py   inlet/numerics, t >= 0.1 s
       step3   notes/c26_step3_checks.py     A-D over the PIV window
       piv     src.post.piv_c26              Table 3 map metrics, wall box
       split   notes/c26_source_split.py     source-gas up share, PIV window
     and the same step3/split on the old-layout run <id><baseline_suffix>;
  5. score against params/c26_acceptance.yaml and rewrite
     <results-dir>/c26_compare_<tag>.csv (one row per finished case).

It exits when every case is terminal.  It does not edit results/run_log.csv:
marking a run done stays a manual review step (as src/aws/poll.py).  A case
whose <run>/post/score.json exists is not post-processed again, so a
restarted watcher picks up where it stopped.  --local skips AWS and scores
run dirs already on disk (<case>/aws_out if it holds a manifest, else
<case>), e.g. to check the pipeline on a local smoke.
"""
import argparse
import contextlib
import csv
import json
import math
import os
import re
import subprocess
import sys
import time
import traceback
from pathlib import Path

os.environ.setdefault('MPLBACKEND', 'Agg')

import numpy as np  # noqa: E402
import yaml  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / 'notes'))
import c26_source_split  # noqa: E402
import c26_step3_checks  # noqa: E402
from src.correlations import G_STD  # noqa: E402
from src.post import piv_c26  # noqa: E402

CASES = REPO / 'cases' / 'C26_3d'
RUN_LOG = REPO / 'results' / 'run_log.csv'
CRITERIA = REPO / 'params' / 'c26_acceptance.yaml'
TERMINAL = ('SUCCEEDED', 'FAILED')
NAN = float('nan')


def say(msg):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


def job_id_for(case, run_log=RUN_LOG):
    """aws_job_id of the last aws-batch row for this case in run_log.csv."""
    jid = None
    with open(run_log, newline='') as f:
        for row in csv.DictReader(f):
            if row['step'] == 'aws-batch' and row['case'] == case:
                m = re.search(r'aws_job_id=([0-9a-f-]{36})', row['notes'])
                if m:
                    jid = m.group(1)
    if jid is None:
        sys.exit(f'no aws-batch row with a job id for {case} in {run_log}')
    return jid


def run_dir(case):
    d = CASES / case
    return d / 'aws_out' if (d / 'aws_out' / 'manifest.json').exists() else d


def pull(job, out_dir, log):
    """s3 sync the job's results prefix into out_dir (src/aws/poll.py pull_results)."""
    from src.aws import config as AWS
    env = {e['name']: e['value'] for e in job.get('container', {}).get('environment', [])}
    prefix = env.get('RESULTS_S3')
    if not prefix:
        raise RuntimeError(f"job {job['jobId']} has no RESULTS_S3")
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(log, 'w') as f:
        subprocess.run(['aws', 's3', 'sync', prefix, str(out_dir), '--region', AWS.REGION,
                        '--only-show-errors'], stdout=f, stderr=subprocess.STDOUT, check=True)


def reached_t(run):
    """Last 'Simulation time reached' in run.log [s], or None."""
    log = run / 'run.log'
    if not log.exists():
        return None
    hits = re.findall(r'Simulation time reached\s*=\s*([0-9.Ee+-]+)', log.read_text(errors='replace'))
    return float(hits[-1]) if hits else None


def fwhm(x, v, v_axis):
    """Full width at half the axis value of a probe profile [m]: walk out from
    the bin nearest x = 0 until v drops below v_axis / 2, interpolate linearly
    between bin centres.  NaN if the profile never drops that far or v_axis <= 0."""
    x, v = np.asarray(x, float), np.asarray(v, float)
    if not len(x) or not v_axis > 0:
        return NAN
    half = 0.5 * v_axis
    i0 = int(np.argmin(np.abs(x)))
    edges = []
    for step in (-1, 1):
        i = i0
        while 0 <= i + step < len(x) and v[i + step] >= half:
            i += step
        j = i + step
        if not 0 <= j < len(x):
            return NAN
        edges.append(x[i] + (half - v[i]) * (x[j] - x[i]) / (v[j] - v[i]))
    return abs(edges[1] - edges[0])


def w_over_a(m):
    """Bed weight over the cone cross-section at the static bed top [Pa]
    (as notes/smoke_validate_c26.py run_check_smoke)."""
    r = min(m['R_t_m'] + m['H_static_m'] * math.tan(math.radians(m['cone_half_angle_deg'])), m['R_c_m'])
    return m['inventory_kg'] * G_STD / (math.pi * r * r)


def step3_scalars(res, fwhm_y):
    pr = res['probes']
    axis = {y: float(pr[y]['axis_v']) for y in sorted(pr)}
    y_pk = max(axis, key=axis.get)
    c2 = {y: max(abs(float(u)) for u in pr[y]['outer_v']) / abs(axis[y]) for y in axis}
    y_c2 = max(c2, key=c2.get)
    yf = min(pr, key=lambda y: abs(y - fwhm_y))
    return {
        'n_frames': int(res['n_frames']),
        'A_frac': float(res['fountain_frac']), 'A_median': float(res['fountain_median']),
        'B_peak_hz': float(res['B1']['peak_hz']) if res['B1'] else NAN,
        'B2_peak_hz': float(res['B2']['peak_hz']) if res['B2'] else NAN,
        'dp_mean_pa': float(res['B1']['mean_pa']) if res['B1'] else NAN,
        **{f'axis_vp_{y * 1e3:.0f}mm': v for y, v in axis.items()},
        'C1_peak': axis[y_pk], 'C1_y_mm': y_pk * 1e3,
        'C2_worst': c2[y_c2], 'C2_worst_y_mm': y_c2 * 1e3,
        'D': float(res['wall_dcg']), 'D_sd': float(res['wall_dcg_sd']),
        'fwhm_mm': fwhm(pr[yf]['profile_x'], pr[yf]['profile'], axis[yf]) * 1e3,
        'fwhm_y_mm': yf * 1e3,
    }


def read_piv(results_dir, name):
    p = Path(results_dir) / f'c26_piv_{name}.csv'
    if not p.exists():
        return None
    with p.open(newline='') as f:
        row = next(csv.DictReader(f))
    return {k: (float(v) if re.fullmatch(r'[-+0-9.eE]+', v or '') else v) for k, v in row.items()}


def band(v, lo_hi, rtol=1e-9):
    """PASS if lo <= v <= hi.  rtol absorbs floating-point error only, e.g. a
    PSD bin at 10 Hz computed as 9.999999999999991 (k / (N dt))."""
    lo, hi = lo_hi
    tol = rtol * max(abs(lo), abs(hi))
    return 'PASS' if lo - tol <= v <= hi + tol else 'FAIL'


def ab_verdict(d_new, d_base, ab):
    """Grid-layout A/B rule (params/c26_acceptance.yaml grid_layout_ab)."""
    ch = (d_new - d_base) / d_base
    if ch <= -ab['substantial_drop_frac']:
        return ch, 'substantial'
    if abs(ch) < ab['no_effect_frac']:
        return ch, 'no effect'
    return ch, 'inconclusive'


def piv_match(piv):
    """C26-model-floor rule: CC, SSI >= the lower and MSE <= the higher of
    C26's Table 3 and same-pipeline scores."""
    cc = min(piv['CC_table3'], piv['CC_c26model_same_pipeline'])
    ssi = min(piv['SSI_table3'], piv['SSI_c26model_same_pipeline'])
    mse = max(piv['MSE_m2_s2_table3'], piv['MSE_m2_s2_c26model_same_pipeline'])
    ok = {'CC': piv['CC'] >= cc, 'SSI': piv['SSI'] >= ssi, 'MSE': piv['MSE_m2_s2'] <= mse}
    return ok, {'CC': cc, 'SSI': ssi, 'MSE': mse}


class Steps:
    """Run named post-processing steps, each logged to <post>/<name>.log; a
    failing step is logged and returns None so the others still run."""

    def __init__(self, post, label):
        self.post, self.label = post, label
        post.mkdir(parents=True, exist_ok=True)

    def __call__(self, name, fn):
        log, t = self.post / f'{name}.log', time.time()
        try:
            with log.open('w') as f, contextlib.redirect_stdout(f), contextlib.redirect_stderr(f):
                out = fn()
            say(f'{self.label}: {name} done ({time.time() - t:.0f} s)')
            return out
        except (Exception, SystemExit):   # SystemExit: the scripts' sys.exit
            with log.open('a') as f:
                f.write(traceback.format_exc())
            say(f'{self.label}: {name} FAILED, see {log}')
            return None


def health(run, log_dir):
    log = log_dir / 'health.log'
    with log.open('w') as f:
        subprocess.run([sys.executable, str(REPO / 'notes' / 'smoke_validate_c26.py'), str(run)],
                       stdout=f, stderr=subprocess.STDOUT, cwd=REPO)
    lines = [ln for ln in log.read_text().splitlines() if ln.startswith('OVERALL:')]
    return lines[-1].strip() if lines else 'no OVERALL line'


def process(case, run, crit, results_dir):
    """Post-process one finished case and return its score row."""
    post = run / 'post'
    if (post / 'score.json').exists():
        return json.loads((post / 'score.json').read_text())
    m = json.loads((run / 'manifest.json').read_text())
    cid = str(m['case_id'])
    t0, t1 = float(m['piv_t_start_s']), float(m['tstop_s']) + 1e-4
    ab = crit['grid_layout_ab']
    steps = Steps(post, case)
    row = {'case': case, 'case_id': cid, 'flow_slpm': m['flow_slpm'],
           'window_s': f'{t0:.3f}-{m["tstop_s"]:.3f}', 'reached_t_s': reached_t(run)}

    row['health'] = steps('health', lambda: health(run, post))
    s3 = steps('step3', lambda: c26_step3_checks.main(run))
    # repo-relative when run from the repo, so the committed PIV CSV has no home path
    piv_run = run.relative_to(REPO) if Path.cwd().resolve() == REPO else run
    steps('piv', lambda: piv_c26.main([str(piv_run), '--out-dir', str(results_dir)]))
    sp = steps('split', lambda: c26_source_split.main(run, t0, t1))

    base = run_dir(f'{cid}{ab["baseline_suffix"]}')
    bpost = post / 'baseline'
    bsteps = Steps(bpost, f'{case} baseline {base.relative_to(CASES)}')
    b3 = bsteps('step3', lambda: c26_step3_checks.main(base)) if base.exists() else None
    bsp = bsteps('split', lambda: c26_source_split.main(base, t0, t1)) if base.exists() else None

    new = step3_scalars(s3, ab['fwhm_probe_y_m']) if s3 else {}
    old = step3_scalars(b3, ab['fwhm_probe_y_m']) if b3 else {}
    name = run.parent.name if run.name == 'aws_out' else run.name
    piv = read_piv(results_dir, name)
    bpiv = read_piv(REPO / 'results', base.parent.name if base.name == 'aws_out' else base.name)
    wa = w_over_a(m)

    st = crit['step3']
    gated = cid in st['gated_cases']
    ref = '' if gated else ' (ref)'
    # A, B and C2 come from C26's YSZ 60 g observations: they gate C26's own
    # cases (those with Table 3 targets) and are a reference for c04-c09
    c26_case = bool(m.get('piv_targets_table3'))
    ref_c26 = '' if c26_case else ' (ref)'
    if new:
        row.update({
            'A_fountain_frac': new['A_frac'],
            'A': ('PASS' if new['A_frac'] >= st['A_fountain_min_frame_frac'] else 'FAIL') + ref_c26,
            'B_psd_peak_hz': new['B_peak_hz'], 'B': band(new['B_peak_hz'], st['B_psd_peak_hz']) + ref_c26,
            'C1_axis_vp_m_s': new['C1_peak'], 'C1_y_mm': new['C1_y_mm'],
            'C1': band(new['C1_peak'], st['C1_axis_vp_m_s']) + ref,
            'C2_outer_over_axis': new['C2_worst'], 'C2_y_mm': new['C2_worst_y_mm'],
            'C2': ('PASS' if new['C2_worst'] < st['C2_outer_over_axis_max'] else 'FAIL') + ref_c26,
            'D_wall_m_s': new['D'], 'D_wall_sd_m_s': new['D_sd'],
            'D': band(new['D'], st['D_wall_m_s']) + ref,
            'axis_vp_20mm_m_s': new.get('axis_vp_20mm', NAN), 'fwhm_20mm_mm': new['fwhm_mm'],
            'dp_over_WA': new['dp_mean_pa'] / wa,
        })
        d0 = ab['baseline_D_m_s'].get(cid)
        if d0:
            ch, verdict = ab_verdict(new['D'], d0, ab)
            row.update({'D_baseline_preregistered_m_s': d0, 'D_change_frac': ch, 'grid_ab_verdict': verdict})
    if old:
        row.update({'old_D_wall_m_s': old['D'], 'old_C1_axis_vp_m_s': old['C1_peak'],
                    'old_axis_vp_20mm_m_s': old.get('axis_vp_20mm', NAN), 'old_fwhm_20mm_mm': old['fwhm_mm'],
                    'old_A_fountain_frac': old['A_frac'], 'old_B_psd_peak_hz': old['B_peak_hz'],
                    'old_C2_outer_over_axis': old['C2_worst']})
    if sp:
        row.update({'up_share': sp['up_share'], 'closure': sp['closure']})
    if bsp:
        row['old_up_share'] = bsp['up_share']
    if piv and not c26_case:
        # no C26 map or wall target: piv_c26 leaves the scores and targets NaN
        row.update({'map_vs_c26_model': 'n/a (no C26 data for this case)',
                    'map_mean_m_s': piv['map_mean_m_s'], 'empty_window_frac': piv['empty_window_frac'],
                    'wall_piv_m_s': piv['wall_piv_m_s']})
    elif piv:
        ok, floor = piv_match(piv)
        row.update({'CC': piv['CC'], 'SSI': piv['SSI'], 'MSE_m2_s2': piv['MSE_m2_s2'],
                    'CC_floor': floor['CC'], 'SSI_floor': floor['SSI'], 'MSE_ceiling': floor['MSE'],
                    'CC_table3': piv['CC_table3'], 'SSI_table3': piv['SSI_table3'],
                    'MSE_table3': piv['MSE_m2_s2_table3'],
                    'map_vs_c26_model': 'MATCH' if all(ok.values()) else
                    'NO MATCH (' + ', '.join(k for k, v in ok.items() if not v) + ')',
                    'map_mean_m_s': piv['map_mean_m_s'], 'empty_window_frac': piv['empty_window_frac'],
                    'wall_piv_m_s': piv['wall_piv_m_s'], 'wall_exp_m_s': piv['wall_exp_m_s'],
                    'wall_c26model_m_s': piv['wall_c26model_m_s'],
                    'wall_piv_over_exp': piv['wall_piv_m_s'] / piv['wall_exp_m_s'],
                    'wall_piv_over_c26model': piv['wall_piv_m_s'] / piv['wall_c26model_m_s']})
    if bpiv:
        row.update({'old_CC': bpiv['CC'], 'old_SSI': bpiv['SSI'], 'old_MSE_m2_s2': bpiv['MSE_m2_s2'],
                    'old_wall_piv_m_s': bpiv['wall_piv_m_s']})
    row = {k: (float(v) if isinstance(v, (np.floating, np.integer)) else v) for k, v in row.items()}
    (post / 'score.json').write_text(json.dumps(row, indent=1))
    return row


def report(row):
    """One case's scores as text for the watch log."""
    g = row.get
    f = lambda v, fmt: format(v, fmt) if isinstance(v, (int, float)) and v == v else 'n/a'  # noqa: E731
    lines = [f"== {row['case']} ({g('flow_slpm')} slpm), window {g('window_s')} s, reached t = {g('reached_t_s')}",
             f"  health: {g('health')}",
             f"  A fountain {f(g('A_fountain_frac'), '.0%')} of frames: {g('A')}   "
             f"B PSD peak {f(g('B_psd_peak_hz'), '.1f')} Hz: {g('B')}",
             f"  C1 axis v_p {f(g('C1_axis_vp_m_s'), '.3f')} m/s at y = {f(g('C1_y_mm'), '.0f')} mm: {g('C1')}   "
             f"(old layout {f(g('old_C1_axis_vp_m_s'), '.3f')})",
             f"  C2 outer/axis {f(g('C2_outer_over_axis'), '.3f')} at y = {f(g('C2_y_mm'), '.0f')} mm: {g('C2')}",
             f"  D wall {f(g('D_wall_m_s'), '.2e')} m/s: {g('D')}   old layout {f(g('old_D_wall_m_s'), '.2e')}; "
             f"A/B vs {f(g('D_baseline_preregistered_m_s'), '.2e')}: {f(g('D_change_frac'), '+.0%')} -> "
             f"{g('grid_ab_verdict')}",
             f"  spout y = 20 mm: axis {f(g('axis_vp_20mm_m_s'), '.3f')} m/s, FWHM {f(g('fwhm_20mm_mm'), '.1f')} mm "
             f"(old {f(g('old_axis_vp_20mm_m_s'), '.3f')} m/s, {f(g('old_fwhm_20mm_mm'), '.1f')} mm); "
             f"dP/(W/A) {f(g('dp_over_WA'), '.2f')}; up share {f(g('up_share'), '.2f')} "
             f"(old {f(g('old_up_share'), '.2f')})",
             f"  PIV map CC/SSI/MSE {f(g('CC'), '.3f')} / {f(g('SSI'), '.3f')} / {f(g('MSE_m2_s2'), '.4f')} "
             f"vs C26 floor {f(g('CC_floor'), '.3f')} / {f(g('SSI_floor'), '.3f')} / {f(g('MSE_ceiling'), '.4f')}: "
             f"{g('map_vs_c26_model')}   (old {f(g('old_CC'), '.3f')} / {f(g('old_SSI'), '.3f')} / "
             f"{f(g('old_MSE_m2_s2'), '.4f')})",
             f"  PIV wall {f(g('wall_piv_m_s'), '.2e')} m/s = {f(g('wall_piv_over_exp'), '.1f')}x exp, "
             f"{f(g('wall_piv_over_c26model'), '.1f')}x C26 model (old {f(g('old_wall_piv_m_s'), '.2e')})"]
    return '\n'.join(lines)


def write_table(rows, path):
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('cases', nargs='+', help='run-dir names under cases/C26_3d')
    ap.add_argument('--tag', required=True, help='results/c26_compare_<tag>.csv')
    ap.add_argument('--interval', type=int, default=300, help='poll interval [s]')
    ap.add_argument('--results-dir', type=Path, default=REPO / 'results')
    ap.add_argument('--local', action='store_true', help='no AWS: score run dirs already on disk')
    ap.add_argument('--baseline-suffix', default=None,
                    help='old-layout run dir suffix (default: c26_acceptance.yaml grid_layout_ab)')
    a = ap.parse_args(argv)
    crit = yaml.safe_load(CRITERIA.read_text())
    if a.baseline_suffix is not None:
        crit['grid_layout_ab']['baseline_suffix'] = a.baseline_suffix
    a.results_dir.mkdir(parents=True, exist_ok=True)
    out = a.results_dir / f'c26_compare_{a.tag}.csv'
    rows = {}

    def finish(case, run):
        say(f'{case}: post-processing {run.relative_to(REPO)}')
        rows[case] = process(case, run, crit, a.results_dir)
        write_table([rows[c] for c in a.cases if c in rows], out)
        say(f'{case}: scored -> {out.relative_to(REPO) if out.is_relative_to(REPO) else out}\n{report(rows[case])}')

    if a.local:
        for case in a.cases:
            finish(case, run_dir(case))
        return 0

    from src.aws.cli import batch_describe_jobs
    pending = {c: job_id_for(c) for c in a.cases}
    say('watching ' + ', '.join(f'{c} ({j[:8]})' for c, j in pending.items()) + f' every {a.interval} s')
    while pending:
        try:
            jobs = {j['jobId']: j for j in batch_describe_jobs(list(pending.values()))}
        except RuntimeError as e:   # aws CLI failed (expired `aws login`, network); keep watching
            say(f'describe-jobs failed, retrying in {a.interval} s: {e}')
            time.sleep(a.interval)
            continue
        say('status: ' + ', '.join(f'{c} {jobs[j]["status"]}' for c, j in pending.items()))
        for case, jid in list(pending.items()):
            job = jobs[jid]
            if job['status'] not in TERMINAL:
                continue
            del pending[case]
            if job['status'] == 'FAILED':
                say(f"{case}: job FAILED ({job.get('statusReason', '')}); not pulled")
                rows[case] = {'case': case, 'health': f"AWS job FAILED: {job.get('statusReason', '')}"}
                continue
            run = CASES / case / 'aws_out'
            say(f'{case}: SUCCEEDED, pulling into {run.relative_to(REPO)}')
            try:
                pull(job, run, run / 'pull.log')
            except Exception as e:  # keep watching the other cases
                say(f'{case}: pull FAILED ({e}); see {run / "pull.log"}')
                continue
            t = reached_t(run)
            say(f'{case}: run.log reached t = {t} s')
            finish(case, run)
        if pending:
            time.sleep(a.interval)
    write_table([rows[c] for c in a.cases if c in rows], out)
    say(f'all cases terminal; table {out}')
    return 0


if __name__ == '__main__':
    sys.exit(main())

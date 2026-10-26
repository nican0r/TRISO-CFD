"""Choose the source for the c01-c09 3 s sweep (params/c26_acceptance.yaml
sweep_source_choice) from c02_apex_mom_v232 against c02_point_mom_v232.

    .venv/bin/python notes/c26_sweep_source_choice.py            # score, print the choice
    .venv/bin/python notes/c26_sweep_source_choice.py --submit   # ... and submit its 8 decks

Steps:
  1. health gate on cases/C26_3d/c02_apex_mom_v232/aws_out: reached 3.0 s
     (SOLIDS_INVENTORY), parcel loss <= 0.5 %, PIV map finite;
  2. if results/c26_piv_map_c02_apex_mom_v232.csv is missing, post-process
     with notes/c26_watch_post.py --local (tag apex_mom_v232, baseline the
     control c02_apex_v232);
  3. fountain reach and lower-map ratio as registered; the one-cell
     reference is recomputed from its map and must reproduce the yaml value;
  4. print the choice and the submit commands; with --submit run
     src.aws.submit for each case in the registered order, stopping at the
     first failure.

Job status (SUCCEEDED) is checked by the waiter that pulls the run, not here.
Written 2026-10-25 (notes/c26_conry2026_case.md, "Sweep source choice").

SUPERSEDED 2026-10-26: the rule chose one-cell and its sweep was submitted,
then the user cancelled it and ran the sweep in the control configuration
(c0X_apex_v232, 2x2 mass-only), the closest match to C26's model among the
three sources (params/c26_acceptance.yaml sweep_source_choice
.decision_2026_10_26; notes "Sweep source choice", Decision).  Scoring
still works; --submit refuses so the one-cell decks cannot be resubmitted
by accident.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import yaml

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from src.post.piv_c26 import load_fig4  # noqa: E402

RUN = "c02_apex_mom_v232"
REF = "c02_point_mom_v232"
CONTROL = "c02_apex_v232"
ORDER = ["c01", "c03", "c04", "c05", "c06", "c08", "c09", "c07"]
SUFFIX = {"2x2": "_apex_mom_v232", "one-cell": "_point_mom_v232"}
REACH_V = 0.2          # m/s, point_3s fountain_reach threshold
REACH_MIN = 0.85       # sweep_source_choice choose_2x2
LOWER_FRAC = 0.40      # rows below this fraction of the map height
MAX_LOSS = 0.005       # parcel loss fraction, health gate
FLOW_SLPM = 12


def load_map(run):
    """Synthetic-PIV map [m/s], row 0 at the bottom (the CSV stores the top row first)."""
    return np.flipud(np.loadtxt(REPO / "results" / f"c26_piv_map_{run}.csv", delimiter=",", comments="#"))


def fountain_reach(m):
    rmax = np.nanmax(m, axis=1)
    hit = np.nonzero(rmax >= REACH_V)[0]
    return (hit.max() + 0.5) / m.shape[0] if hit.size else 0.0


def lower_map_ratio(m, model):
    f = (np.arange(model.shape[0]) + 0.5) / model.shape[0]
    low = f < LOWER_FRAC
    return float(np.nanmean(m[low]) / np.nanmean(model[low]))


def health(out):
    inv = sorted(out.glob("C26_*_SOLIDS_INVENTORY.csv"))
    if not inv:
        return False, "no SOLIDS_INVENTORY.csv"
    d = np.loadtxt(inv[0], delimiter=",", comments="#")
    t1, n0, n1 = d[-1, 0], d[0, 2], d[-1, 2]
    loss = (n0 - n1) / n0
    ok = t1 >= 3.0 - 1e-3 and loss <= MAX_LOSS
    return ok, f"t = {t1:.4f} s, parcels {n0:.0f} -> {n1:.0f} ({100 * loss:.2f} % lost)"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--submit", action="store_true", help="submit the chosen branch's decks")
    ap.add_argument("--run", default=RUN, help="run scored as the 2x2 candidate (checks only)")
    a = ap.parse_args(argv)
    if a.submit:
        sys.exit("--submit refused: the sweep_source_choice rule is superseded by "
                 "decision_2026_10_26 (control configuration, c0X_apex_v232, submitted "
                 "2026-10-26); see params/c26_acceptance.yaml")
    if a.submit and a.run != RUN:
        sys.exit(f"--submit decides on {RUN} only")

    rule = yaml.safe_load(open(REPO / "params" / "c26_acceptance.yaml"))["sweep_source_choice"]
    ref = rule["reference_c02_point_mom_v232"]
    model = load_fig4("model", FLOW_SLPM)
    ref_lower = lower_map_ratio(load_map(REF), model)
    if abs(ref_lower - ref["lower_map_ratio"]) > 1e-3:
        sys.exit(f"{REF} lower-map ratio {ref_lower:.4f} does not reproduce the registered "
                 f"{ref['lower_map_ratio']}: stop, pipeline changed")

    out = REPO / "cases" / "C26_3d" / a.run / "aws_out"
    ok, why = health(out)
    print(f"{a.run} health: {why} -> {'PASS' if ok else 'FAIL'}")
    choice, reason = "one-cell", "health gate failed"
    if ok:
        if not (REPO / "results" / f"c26_piv_map_{a.run}.csv").exists():
            subprocess.run([sys.executable, "notes/c26_watch_post.py", a.run, "--tag", "apex_mom_v232",
                            "--local", "--baseline-suffix", "_apex_v232"],
                           cwd=REPO, env={**os.environ, "MPLBACKEND": "Agg"}, check=True)
        m = load_map(a.run)
        if not np.isfinite(m).all():
            print(f"{a.run} health: PIV map has {np.size(m) - np.isfinite(m).sum()} non-finite cells -> FAIL")
        else:
            reach, lower = fountain_reach(m), lower_map_ratio(m, model)
            # compared unrounded against the reference recomputed from its map:
            # the registered 1.297 is 1.2969 rounded, and a tie goes to one-cell
            print(f"{a.run}: fountain reach {reach:.2f} (>= {REACH_MIN}), lower-map ratio {lower:.4f} "
                  f"(< {ref_lower:.4f}, {REF})")
            if reach >= REACH_MIN and lower < ref_lower:
                choice, reason = "2x2", "reach and lower map both qualify"
            else:
                choice, reason = "one-cell", "2x2 does not qualify"
    print(f"choice: {choice} ({reason})")

    cmds = []
    for c in ORDER:
        case = f"{c}{SUFFIX[choice]}"
        d = REPO / "cases" / "C26_3d" / case
        if not (d / f"C26_{case}.mfx").exists():
            sys.exit(f"missing deck {d}/C26_{case}.mfx")
        cmds.append([sys.executable, "-m", "src.aws.submit", "--case-dir", f"cases/C26_3d/{case}",
                     "--mfx", f"C26_{case}.mfx", "--job-name", f"c26-{case.replace('_', '-')}"])
    for cmd in cmds:
        print("  " + " ".join(cmd[1:]))
    if not a.submit:
        return 0
    for cmd in cmds:
        r = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True)
        job = [ln for ln in r.stdout.splitlines() if "job id" in ln]
        print(f"{cmd[4]}: exit {r.returncode} {job[-1] if job else ''}")
        if r.returncode != 0:
            print(r.stdout[-2000:], r.stderr[-2000:])
            return r.returncode
    return 0


if __name__ == "__main__":
    sys.exit(main())

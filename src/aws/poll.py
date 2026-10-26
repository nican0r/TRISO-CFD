"""Poll an AWS Batch MFiX job and pull its results on success.

Usage:
    .venv/bin/python -m src.aws.poll <job-id>                 # one-shot
    .venv/bin/python -m src.aws.poll <job-id> --watch         # loop until done
    .venv/bin/python -m src.aws.poll <job-id> --pull PATH     # on success,
                                                               # s3 sync results
                                                               # into PATH

On terminal state (SUCCEEDED / FAILED), optionally pulls results and prints
log stream handles.  In --watch mode a failing ``aws batch describe-jobs``
call (expired ``aws login`` session, network blip) is logged and retried at
the next interval instead of ending the watch: an overnight watcher must
survive a credential expiry and resume once ``aws login`` is run again
(2026-10-24: a watcher died ~1 h into a 15 h job for exactly this reason). Does not touch results/run_log.csv - the submit step
already appended a 'submitted' row; marking a run done is a manual review
step.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

from src.aws import config as C
from src.aws.cli import batch_describe_jobs


TERMINAL = {"SUCCEEDED", "FAILED"}


def _log(msg: str) -> None:
    print(f"[poll] {msg}", flush=True)


def describe(job_id: str) -> dict:
    jobs = batch_describe_jobs([job_id])
    if not jobs:
        sys.exit(f"job not found: {job_id}")
    return jobs[0]


def summary(j: dict) -> dict:
    out = {
        "jobId": j["jobId"],
        "jobName": j["jobName"],
        "status": j["status"],
        "statusReason": j.get("statusReason", ""),
    }
    if "arrayProperties" in j:
        out["arrayProperties"] = j["arrayProperties"]
    attempts = j.get("attempts", [])
    if attempts:
        a = attempts[-1]
        out["exitCode"] = a.get("container", {}).get("exitCode")
        out["logStreamName"] = a.get("container", {}).get("logStreamName")
    for k in ("createdAt", "startedAt", "stoppedAt"):
        if k in j:
            out[k] = j[k]
    return out


def pull_results(j: dict, out_dir: Path) -> None:
    """Pull the s3 results prefix recorded on the submitting side."""
    env = {e["name"]: e["value"] for e in
            j.get("container", {}).get("environment", [])}
    s3_prefix = env.get("RESULTS_S3") or env.get("RESULTS_S3_ARRAY")
    if not s3_prefix:
        _log("no RESULTS_S3/RESULTS_S3_ARRAY on job; skipping pull")
        return
    out_dir.mkdir(parents=True, exist_ok=True)
    _log(f"s3 sync {s3_prefix} -> {out_dir}")
    subprocess.check_call([
        "aws", "s3", "sync", s3_prefix, str(out_dir),
        "--region", C.REGION,
    ])


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("job_id")
    p.add_argument("--watch", action="store_true")
    p.add_argument("--interval", type=int, default=30, help="seconds")
    p.add_argument("--pull", type=Path, default=None,
                     help="on success, s3 sync results into this directory")
    args = p.parse_args()

    while True:
        try:
            j = describe(args.job_id)
        except RuntimeError as e:          # aws CLI failed (creds, network)
            if not args.watch:
                raise
            _log(f"describe-jobs failed, retrying in {args.interval} s: {e}")
            time.sleep(args.interval)
            continue
        s = summary(j)
        print(json.dumps(s, indent=2, default=str))
        if not args.watch:
            break
        if s["status"] in TERMINAL:
            break
        time.sleep(args.interval)

    if args.pull is not None and j["status"] == "SUCCEEDED":
        pull_results(j, args.pull)


if __name__ == "__main__":
    main()

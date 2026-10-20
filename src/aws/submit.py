"""Submit an MFiX case to AWS Batch.

Workflow:
  1. Caller has already run their day-N driver so a case directory exists
     under cases/<...>/ with the .mfx + STL + manifest.json.
  2. This script `aws s3 sync`s the case up to the project bucket, submits
     a Batch job that points at it, and appends a 'submitted' row to
     results/run_log.csv so provenance is attached before the job even runs.

Single-run usage:
    .venv/bin/python -m src.aws.submit \
        --case-dir cases/V_3d/day10 \
        --mfx V_3d_day10.mfx \
        --job-name d10-pic \
        --vcpu 16 --memory 30000

Array-run usage (Day 13 sweep):
    .venv/bin/python -m src.aws.submit \
        --array-root cases/sweep13 \
        --mfx V_3d_sweep.mfx \
        --job-name d13-sweep \
        --array-size 16 --vcpu 16 --memory 30000

Sizing from the case manifest: if <case>/manifest.json (array mode: case_00's)
has "aws_vcpu", then any of --vcpu, --memory, --nranks and --nodes left unset
come from its aws_vcpu, aws_memory_mb, nranks and nodes_i/j/k.  C26 case
generators write these (4 vCPU / 4 ranks, src/run_C26_3d.py), e.g.:

    .venv/bin/python -m src.aws.submit \
        --case-dir cases/C26_3d/c02 --mfx C26_c02.mfx --job-name c26-c02

Explicit flags win.  Manifests without "aws_vcpu" keep the old defaults
(config.DEFAULT_VCPU / DEFAULT_MEMORY_MB, nranks = vcpu, entrypoint nodes).

In array mode, --array-root must contain subdirectories case_00/, case_01/, ...,
each a complete case dir. Each one is synced to s3://<bucket>/sweep/<name>/case_NN/
and the Batch array job reads AWS_BATCH_JOB_ARRAY_INDEX inside the container.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import subprocess
import sys
from pathlib import Path

import boto3

from src.aws import config as C

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
RUN_LOG = REPO_ROOT / "results" / "run_log.csv"


def _log(msg: str) -> None:
    print(f"[submit] {msg}", flush=True)


def s3_sync(src: Path, dst_url: str) -> None:
    _log(f"s3 sync {src} -> {dst_url}")
    subprocess.check_call([
        "aws", "s3", "sync", str(src), dst_url,
        "--region", C.REGION,
        "--exclude", "*.pyc", "--exclude", "__pycache__/*",
        "--exclude", ".DS_Store",
    ])


def append_run_log_row(row: dict) -> None:
    """Append a row to results/run_log.csv using the existing header order."""
    with RUN_LOG.open("r", newline="") as f:
        header = next(csv.reader(f))
    out = {k: row.get(k, "") for k in header}
    with RUN_LOG.open("a", newline="") as f:
        csv.DictWriter(f, fieldnames=header).writerow(out)
    _log(f"appended run_log.csv row (status={out.get('status')})")


def _nodes_env(nodes: tuple[int, int, int] | None) -> list[dict]:
    """Emit NODESI/NODESJ/NODESK env entries if `nodes` is set.

    docker/entrypoint.sh reads these; absent them it defaults to (1, NP, 1)
    which only divides cleanly if NP divides jmax.  For UTK's 60x220x60
    mesh at NP=32 the default fails -- callers must pass an explicit
    decomposition that divides each axis.
    """
    if nodes is None:
        return []
    ni, nj, nk = nodes
    return [
        {"name": "NODESI", "value": str(int(ni))},
        {"name": "NODESJ", "value": str(int(nj))},
        {"name": "NODESK", "value": str(int(nk))},
    ]


def manifest_sizing(case_dir: Path) -> dict:
    """AWS sizing from case_dir/manifest.json, or {} if it has no aws_vcpu."""
    path = case_dir / "manifest.json"
    if not path.is_file():
        return {}
    m = json.loads(path.read_text())
    if "aws_vcpu" not in m:
        return {}
    out = {"vcpu": int(m["aws_vcpu"])}
    if "aws_memory_mb" in m:
        out["memory"] = int(m["aws_memory_mb"])
    if "nranks" in m:
        out["nranks"] = int(m["nranks"])
    if all(f"nodes_{a}" in m for a in "ijk"):
        out["nodes"] = tuple(int(m[f"nodes_{a}"]) for a in "ijk")
    return out


def submit_single(sess, case_name: str, mfx: str, job_name: str,
                   vcpu: int, memory: int, nranks: int,
                   extra_args: str,
                   nodes: tuple[int, int, int] | None = None) -> str:
    bucket = C.S3_BUCKET
    case_s3 = f"s3://{bucket}/cases/{case_name}/"
    results_s3 = f"s3://{bucket}/results/{case_name}/"
    batch = sess.client("batch")
    _log(f"submit {job_name} (single)")
    resp = batch.submit_job(
        jobName=job_name, jobQueue=C.JOB_QUEUE,
        jobDefinition=C.JOB_DEFINITION,
        containerOverrides={
            "resourceRequirements": [
                {"type": "VCPU", "value": str(vcpu)},
                {"type": "MEMORY", "value": str(memory)},
            ],
            "environment": [
                {"name": "CASE_S3", "value": case_s3},
                {"name": "RESULTS_S3", "value": results_s3},
                {"name": "MFX_FILE", "value": mfx},
                {"name": "NP", "value": str(nranks)},
                {"name": "EXTRA_MFIX_ARGS", "value": extra_args},
                *_nodes_env(nodes),
            ],
        },
    )
    return resp["jobId"]


def submit_array(sess, case_name: str, mfx: str, job_name: str, size: int,
                  vcpu: int, memory: int, nranks: int,
                  extra_args: str,
                  nodes: tuple[int, int, int] | None = None) -> str:
    bucket = C.S3_BUCKET
    case_s3 = f"s3://{bucket}/sweep/{case_name}/"
    results_s3 = f"s3://{bucket}/sweep/{case_name}/results/"
    batch = sess.client("batch")
    _log(f"submit {job_name} (array size={size})")
    resp = batch.submit_job(
        jobName=job_name, jobQueue=C.JOB_QUEUE,
        jobDefinition=C.JOB_DEFINITION,
        arrayProperties={"size": size},
        containerOverrides={
            "resourceRequirements": [
                {"type": "VCPU", "value": str(vcpu)},
                {"type": "MEMORY", "value": str(memory)},
            ],
            "environment": [
                {"name": "CASE_S3_ARRAY", "value": case_s3},
                {"name": "RESULTS_S3_ARRAY", "value": results_s3},
                {"name": "MFX_FILE", "value": mfx},
                {"name": "NP", "value": str(nranks)},
                {"name": "EXTRA_MFIX_ARGS", "value": extra_args},
                *_nodes_env(nodes),
            ],
        },
    )
    return resp["jobId"]


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--case-dir", type=Path,
                   help="single-run case directory (contains the .mfx)")
    g.add_argument("--array-root", type=Path,
                   help="array-run root (contains case_00/, case_01/, ...)")
    p.add_argument("--mfx", required=True, help=".mfx filename inside the case dir")
    p.add_argument("--job-name", required=True)
    p.add_argument("--array-size", type=int, default=None,
                     help="required with --array-root")
    p.add_argument("--vcpu", type=int, default=None,
                     help="default: manifest aws_vcpu, else config.DEFAULT_VCPU")
    p.add_argument("--memory", type=int, default=None,
                     help="MB; default: manifest aws_memory_mb, else "
                          "config.DEFAULT_MEMORY_MB")
    p.add_argument("--nranks", type=int, default=None,
                     help="MPI ranks; default: manifest nranks (if it has "
                          "aws_vcpu), else --vcpu")
    p.add_argument("--nodes", type=str, default=None,
                     metavar="NI,NJ,NK",
                     help=("explicit DMP decomposition passed via "
                           "NODESI/NODESJ/NODESK env vars; required when the "
                           "default (1, NP, 1) axial decomposition does not "
                           "divide the mesh (e.g. UTK 60x220x60 at NP=32)"))
    p.add_argument("--extra-mfix-args", default="")
    p.add_argument("--dry-run", action="store_true",
                     help="skip s3 sync + batch submit; print the plan")
    args = p.parse_args()

    sess = boto3.Session(region_name=C.REGION)
    root = args.case_dir if args.case_dir else args.array_root / "case_00"
    sizing = manifest_sizing(root)
    if sizing:
        _log(f"manifest sizing ({root / 'manifest.json'}): {sizing}")
    args.vcpu = args.vcpu or sizing.get("vcpu", C.DEFAULT_VCPU)
    args.memory = args.memory or sizing.get("memory", C.DEFAULT_MEMORY_MB)
    nranks = args.nranks or sizing.get("nranks", args.vcpu)
    if nranks > args.vcpu:
        sys.exit(f"nranks {nranks} > vcpu {args.vcpu}: ranks would share cores")

    nodes: tuple[int, int, int] | None = None
    if not args.nodes and "nodes" in sizing and args.nranks in (None, sizing.get("nranks")):
        nodes = sizing["nodes"]
    if args.nodes:
        parts = [int(p.strip()) for p in args.nodes.split(",")]
        if len(parts) != 3 or any(p < 1 for p in parts):
            sys.exit(f"--nodes expects NI,NJ,NK positive ints, got {args.nodes!r}")
        if parts[0] * parts[1] * parts[2] != nranks:
            sys.exit(f"--nodes product {parts[0]*parts[1]*parts[2]} != nranks {nranks}")
        nodes = (parts[0], parts[1], parts[2])

    if args.case_dir:
        case_name = args.case_dir.name
        if not args.case_dir.exists():
            sys.exit(f"case dir not found: {args.case_dir}")
        plan = {"mode": "single", "case": case_name, "mfx": args.mfx,
                 "vcpu": args.vcpu, "memory_mb": args.memory, "nranks": nranks,
                 "nodes": list(nodes) if nodes else None}
        _log(json.dumps(plan))
        if args.dry_run:
            return
        s3_sync(args.case_dir, f"s3://{C.S3_BUCKET}/cases/{case_name}/")
        job_id = submit_single(sess, case_name, args.mfx, args.job_name,
                                args.vcpu, args.memory, nranks,
                                args.extra_mfix_args, nodes)
        s3_prefix = f"s3://{C.S3_BUCKET}/results/{case_name}/"
    else:
        if args.array_size is None:
            sys.exit("--array-size is required with --array-root")
        case_name = args.array_root.name
        for i in range(args.array_size):
            sub = args.array_root / f"case_{i:02d}"
            if not sub.exists():
                sys.exit(f"missing array leaf: {sub}")
        plan = {"mode": "array", "case": case_name, "mfx": args.mfx,
                 "size": args.array_size, "vcpu": args.vcpu,
                 "memory_mb": args.memory, "nranks": nranks,
                 "nodes": list(nodes) if nodes else None}
        _log(json.dumps(plan))
        if args.dry_run:
            return
        for i in range(args.array_size):
            sub = args.array_root / f"case_{i:02d}"
            s3_sync(sub, f"s3://{C.S3_BUCKET}/sweep/{case_name}/case_{i:02d}/")
        job_id = submit_array(sess, case_name, args.mfx, args.job_name,
                                args.array_size, args.vcpu, args.memory,
                                nranks, args.extra_mfix_args, nodes)
        s3_prefix = f"s3://{C.S3_BUCKET}/sweep/{case_name}/results/"

    _log(f"job id: {job_id}")
    append_run_log_row({
        "date": dt.date.today().isoformat(),
        "step": "aws-batch",
        "case": case_name,
        "key_params": f"job_name={args.job_name}; vcpu={args.vcpu}; "
                        f"memory_mb={args.memory}; nranks={nranks}",
        "status": "submitted",
        "notes": f"aws_job_id={job_id}; s3={s3_prefix}",
    })


if __name__ == "__main__":
    main()

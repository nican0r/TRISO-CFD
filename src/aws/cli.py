"""Thin wrappers around the AWS CLI for the Batch calls the run tooling makes.

Why the CLI and not boto3: since 2026-10-22 the laptop's credentials come
from the ``aws login`` provider, which boto3 only supports with the optional
``botocore[crt]`` dependency.  The AWS CLI v2 handles that provider natively,
so submit / poll / watch shell out to it (user decision 2026-10-22).
``src/aws/provision.py`` (one-off stack creation) still uses boto3.

Every call is ``aws <service> <op> --region <REGION> --output json`` plus the
op's own arguments; the JSON on stdout is parsed and returned.
"""

from __future__ import annotations

import json
import subprocess

from src.aws import config as C


def aws_json(service: str, op: str, *args: str, region: str = C.REGION) -> dict:
    """Run ``aws <service> <op> ...`` and return its parsed JSON output."""
    cmd = ["aws", service, op, "--region", region, "--output", "json", *args]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"{' '.join(cmd[:3])} failed (exit {res.returncode}): "
                           f"{res.stderr.strip()}")
    return json.loads(res.stdout) if res.stdout.strip() else {}


def batch_submit_job(payload: dict) -> str:
    """``aws batch submit-job --cli-input-json <payload>``; returns the job id.

    ``payload`` uses the API field names (jobName, jobQueue, jobDefinition,
    containerOverrides, arrayProperties, ...), exactly as boto3's
    ``submit_job(**payload)`` would take them.
    """
    resp = aws_json("batch", "submit-job", "--cli-input-json", json.dumps(payload))
    return resp["jobId"]


def batch_describe_jobs(job_ids: list[str]) -> list[dict]:
    """``aws batch describe-jobs --jobs ...``; returns the ``jobs`` list."""
    if not job_ids:
        return []
    return aws_json("batch", "describe-jobs", "--jobs", *job_ids).get("jobs", [])

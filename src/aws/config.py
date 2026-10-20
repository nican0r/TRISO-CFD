"""Project-wide AWS Batch config.

All resource names here are shared across provision / submit / poll so there is
exactly one place to edit if the stack is torn down and re-created.
"""

from __future__ import annotations

import os

REGION = os.environ.get("TRISO_AWS_REGION", "us-east-1")

PROJECT = "triso-cfd"
S3_BUCKET = os.environ.get("TRISO_S3_BUCKET", f"{PROJECT}-mfix")
ECR_REPO = os.environ.get("TRISO_ECR_REPO", f"{PROJECT}/mfix")
IMAGE_TAG = os.environ.get("TRISO_IMAGE_TAG", "26.1.2")

COMPUTE_ENV = f"{PROJECT}-ce"
JOB_QUEUE = f"{PROJECT}-q"
JOB_DEFINITION = f"{PROJECT}-mfix-26-1-2"

INSTANCE_ROLE = f"{PROJECT}-ecs-instance-role"
INSTANCE_PROFILE = f"{PROJECT}-ecs-instance-profile"
JOB_ROLE = f"{PROJECT}-job-role"
EXECUTION_ROLE = f"{PROJECT}-execution-role"
SERVICE_ROLE = f"{PROJECT}-batch-service-role"

INSTANCE_TYPES = ["c7a.4xlarge", "c7a.8xlarge"]
MAX_VCPUS = 64
MIN_VCPUS = 0

DEFAULT_VCPU = 16
DEFAULT_MEMORY_MB = 30000

# C26 (Conry et al. 2026) job sizing: one vCPU per rank, 4 ranks (2x1x2).
# 4 ranks gives the most throughput per vCPU (c01_aws4_timing: 0.234 s/step on
# 4 vCPU vs 0.209 on 16 vCPU / 9 ranks), so a 32-vCPU quota runs 8 jobs.
# Memory: local 4-rank c01 (W=10, 27,782 parcels), measured 2026-10-19:
# steady 0.6-0.7 GB summed over ranks, peak 1.09 GB during mesh/STL
# preprocessing. 6000 MB is ~5x headroom and still fits 4 jobs on a
# c7a.4xlarge (32 GiB). W >= 10 for every case (params/sweep_C26.yaml);
# src/run_C26_3d.py::aws_memory_mb scales this linearly for cases with more
# parcels than c01 (c09, 80 g: ~36k parcels -> 8000 MB).
C26_MEMORY_MB = 6000

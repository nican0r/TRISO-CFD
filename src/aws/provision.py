"""Idempotent provisioning of the AWS Batch stack for MFiX runs.

Creates (or no-ops if already present):
  - S3 bucket               config.S3_BUCKET
  - ECR repository          config.ECR_REPO
  - IAM roles               Batch service, EC2 instance (ECS), job role
  - Instance profile        wraps the EC2 instance role
  - Batch compute env       EC2 on-demand, c7a.4xlarge/8xlarge, BEST_FIT_PROGRESSIVE
  - Batch job queue         wired to the compute env
  - Batch job definition    pointing to <ecr>/<repo>:<tag>, retry=1

Builds and pushes the image (docker/) only with --build
(docker_build_and_push: docker build, tag, ECR login, push).

Usage:
    .venv/bin/python -m src.aws.provision          # provision everything
    .venv/bin/python -m src.aws.provision --build  # also docker build+push
    .venv/bin/python -m src.aws.provision --status # show current state
"""

from __future__ import annotations

import argparse
import json
import subprocess
import time
from pathlib import Path

import boto3
from botocore.exceptions import ClientError

from src.aws import config as C

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DOCKER_DIR = REPO_ROOT / "docker"


def _log(msg: str) -> None:
    print(f"[provision] {msg}", flush=True)


# -- IAM ----------------------------------------------------------------------

BATCH_SERVICE_TRUST = {
    "Version": "2012-10-17",
    "Statement": [{"Effect": "Allow",
                    "Principal": {"Service": "batch.amazonaws.com"},
                    "Action": "sts:AssumeRole"}],
}
EC2_TRUST = {
    "Version": "2012-10-17",
    "Statement": [{"Effect": "Allow",
                    "Principal": {"Service": "ec2.amazonaws.com"},
                    "Action": "sts:AssumeRole"}],
}
ECS_TASKS_TRUST = {
    "Version": "2012-10-17",
    "Statement": [{"Effect": "Allow",
                    "Principal": {"Service": "ecs-tasks.amazonaws.com"},
                    "Action": "sts:AssumeRole"}],
}
S3_RW_POLICY = {
    "Version": "2012-10-17",
    "Statement": [{
        "Effect": "Allow",
        "Action": ["s3:GetObject", "s3:PutObject", "s3:DeleteObject",
                     "s3:ListBucket"],
        "Resource": [f"arn:aws:s3:::{C.S3_BUCKET}",
                       f"arn:aws:s3:::{C.S3_BUCKET}/*"],
    }],
}


def _ensure_role(iam, name: str, trust: dict, managed: list[str],
                 inline: dict | None = None) -> str:
    try:
        role = iam.get_role(RoleName=name)["Role"]
        _log(f"role {name} exists")
    except ClientError as e:
        if e.response["Error"]["Code"] != "NoSuchEntity":
            raise
        _log(f"creating role {name}")
        role = iam.create_role(RoleName=name,
                                 AssumeRolePolicyDocument=json.dumps(trust))["Role"]
    for arn in managed:
        iam.attach_role_policy(RoleName=name, PolicyArn=arn)
    if inline:
        iam.put_role_policy(RoleName=name, PolicyName=f"{name}-inline",
                              PolicyDocument=json.dumps(inline))
    return role["Arn"]


def provision_iam(sess) -> dict:
    iam = sess.client("iam")
    service_arn = _ensure_role(
        iam, C.SERVICE_ROLE, BATCH_SERVICE_TRUST,
        ["arn:aws:iam::aws:policy/service-role/AWSBatchServiceRole"],
    )
    instance_arn = _ensure_role(
        iam, C.INSTANCE_ROLE, EC2_TRUST,
        ["arn:aws:iam::aws:policy/service-role/AmazonEC2ContainerServiceforEC2Role"],
        inline=S3_RW_POLICY,
    )
    # instance profile wraps the EC2 role
    try:
        iam.get_instance_profile(InstanceProfileName=C.INSTANCE_PROFILE)
        _log(f"instance profile {C.INSTANCE_PROFILE} exists")
    except ClientError as e:
        if e.response["Error"]["Code"] != "NoSuchEntity":
            raise
        _log(f"creating instance profile {C.INSTANCE_PROFILE}")
        iam.create_instance_profile(InstanceProfileName=C.INSTANCE_PROFILE)
        iam.add_role_to_instance_profile(InstanceProfileName=C.INSTANCE_PROFILE,
                                           RoleName=C.INSTANCE_ROLE)
        time.sleep(10)  # IAM propagation
    job_arn = _ensure_role(
        iam, C.JOB_ROLE, ECS_TASKS_TRUST, [], inline=S3_RW_POLICY,
    )
    exec_arn = _ensure_role(
        iam, C.EXECUTION_ROLE, ECS_TASKS_TRUST,
        ["arn:aws:iam::aws:policy/service-role/AmazonECSTaskExecutionRolePolicy"],
    )
    return {"service": service_arn, "instance": instance_arn,
            "job": job_arn, "execution": exec_arn}


# -- S3 / ECR -----------------------------------------------------------------

def provision_s3(sess) -> None:
    s3 = sess.client("s3")
    try:
        s3.head_bucket(Bucket=C.S3_BUCKET)
        _log(f"s3 bucket {C.S3_BUCKET} exists")
        return
    except ClientError:
        pass
    _log(f"creating s3 bucket {C.S3_BUCKET}")
    kwargs = {"Bucket": C.S3_BUCKET}
    if C.REGION != "us-east-1":
        kwargs["CreateBucketConfiguration"] = {"LocationConstraint": C.REGION}
    s3.create_bucket(**kwargs)


def provision_ecr(sess) -> str:
    ecr = sess.client("ecr")
    try:
        resp = ecr.describe_repositories(repositoryNames=[C.ECR_REPO])
        uri = resp["repositories"][0]["repositoryUri"]
        _log(f"ecr repo {C.ECR_REPO} exists ({uri})")
    except ClientError as e:
        if e.response["Error"]["Code"] != "RepositoryNotFoundException":
            raise
        _log(f"creating ecr repo {C.ECR_REPO}")
        resp = ecr.create_repository(repositoryName=C.ECR_REPO,
                                      imageScanningConfiguration={"scanOnPush": False})
        uri = resp["repository"]["repositoryUri"]
    return uri


# -- Batch --------------------------------------------------------------------

def provision_compute_env(sess, roles: dict) -> str:
    batch = sess.client("batch")
    ec2 = sess.client("ec2")
    # Default VPC subnets + default security group.
    vpcs = ec2.describe_vpcs(Filters=[{"Name": "isDefault", "Values": ["true"]}])["Vpcs"]
    if not vpcs:
        raise RuntimeError("no default VPC found in region " + C.REGION)
    vpc_id = vpcs[0]["VpcId"]
    subnets = [s["SubnetId"] for s in ec2.describe_subnets(
        Filters=[{"Name": "vpc-id", "Values": [vpc_id]}])["Subnets"]]
    sg = ec2.describe_security_groups(Filters=[
        {"Name": "vpc-id", "Values": [vpc_id]},
        {"Name": "group-name", "Values": ["default"]}])["SecurityGroups"][0]["GroupId"]

    existing = batch.describe_compute_environments(
        computeEnvironments=[C.COMPUTE_ENV])["computeEnvironments"]
    if existing:
        _log(f"compute env {C.COMPUTE_ENV} exists ({existing[0]['status']})")
        return existing[0]["computeEnvironmentArn"]

    _log(f"creating compute env {C.COMPUTE_ENV}")
    resp = batch.create_compute_environment(
        computeEnvironmentName=C.COMPUTE_ENV,
        type="MANAGED",
        state="ENABLED",
        serviceRole=roles["service"],
        computeResources={
            "type": "EC2",
            "allocationStrategy": "BEST_FIT_PROGRESSIVE",
            "minvCpus": C.MIN_VCPUS,
            "maxvCpus": C.MAX_VCPUS,
            "instanceTypes": C.INSTANCE_TYPES,
            "subnets": subnets,
            "securityGroupIds": [sg],
            "instanceRole": f"arn:aws:iam::{sess.client('sts').get_caller_identity()['Account']}:instance-profile/{C.INSTANCE_PROFILE}",
            "tags": {"Project": C.PROJECT},
        },
    )
    return resp["computeEnvironmentArn"]


def provision_job_queue(sess, ce_arn: str) -> str:
    batch = sess.client("batch")
    existing = batch.describe_job_queues(jobQueues=[C.JOB_QUEUE])["jobQueues"]
    if existing:
        _log(f"job queue {C.JOB_QUEUE} exists ({existing[0]['status']})")
        return existing[0]["jobQueueArn"]
    _log("waiting for compute env to become VALID")
    while True:
        ce = batch.describe_compute_environments(
            computeEnvironments=[C.COMPUTE_ENV])["computeEnvironments"][0]
        if ce["status"] == "VALID":
            break
        if ce["status"] == "INVALID":
            raise RuntimeError(f"compute env is INVALID: {ce.get('statusReason')}")
        time.sleep(5)
    _log(f"creating job queue {C.JOB_QUEUE}")
    resp = batch.create_job_queue(
        jobQueueName=C.JOB_QUEUE, state="ENABLED", priority=1,
        computeEnvironmentOrder=[{"order": 1, "computeEnvironment": ce_arn}],
    )
    return resp["jobQueueArn"]


def provision_job_definition(sess, image_uri: str, roles: dict) -> str:
    batch = sess.client("batch")
    _log(f"registering job definition {C.JOB_DEFINITION} -> {image_uri}")
    resp = batch.register_job_definition(
        jobDefinitionName=C.JOB_DEFINITION,
        type="container",
        retryStrategy={"attempts": 1},
        containerProperties={
            "image": image_uri,
            "vcpus": C.DEFAULT_VCPU,
            "memory": C.DEFAULT_MEMORY_MB,
            "jobRoleArn": roles["job"],
            "executionRoleArn": roles["execution"],
            "environment": [],
            "readonlyRootFilesystem": False,
        },
        propagateTags=True,
        tags={"Project": C.PROJECT},
    )
    return resp["jobDefinitionArn"]


# -- Docker build / push ------------------------------------------------------

def docker_build_and_push(image_uri: str) -> None:
    local_tag = f"{C.PROJECT}/mfix:{C.IMAGE_TAG}"
    _log(f"docker build -t {local_tag}")
    subprocess.check_call(
        ["docker", "build", "--platform=linux/amd64", "-t", local_tag,
         str(DOCKER_DIR)])
    _log(f"docker tag {local_tag} -> {image_uri}:{C.IMAGE_TAG}")
    subprocess.check_call(["docker", "tag", local_tag,
                            f"{image_uri}:{C.IMAGE_TAG}"])
    _log("ecr login + push")
    account = image_uri.split(".")[0]
    pw = subprocess.check_output(
        ["aws", "ecr", "get-login-password", "--region", C.REGION])
    subprocess.run(
        ["docker", "login", "--username", "AWS", "--password-stdin",
         f"{account}.dkr.ecr.{C.REGION}.amazonaws.com"], input=pw, check=True)
    subprocess.check_call(["docker", "push", f"{image_uri}:{C.IMAGE_TAG}"])


# -- CLI ----------------------------------------------------------------------

def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--build", action="store_true",
                     help="also build and push the Docker image to ECR")
    p.add_argument("--status", action="store_true",
                     help="print state of the stack and exit")
    args = p.parse_args()

    sess = boto3.Session(region_name=C.REGION)

    if args.status:
        batch = sess.client("batch")
        print(json.dumps({
            "region": C.REGION,
            "bucket": C.S3_BUCKET,
            "ecr_repo": C.ECR_REPO,
            "compute_envs": [c["status"] for c in batch.describe_compute_environments(
                computeEnvironments=[C.COMPUTE_ENV])["computeEnvironments"]],
            "job_queues": [q["status"] for q in batch.describe_job_queues(
                jobQueues=[C.JOB_QUEUE])["jobQueues"]],
        }, indent=2))
        return

    roles = provision_iam(sess)
    provision_s3(sess)
    ecr_uri = provision_ecr(sess)
    if args.build:
        docker_build_and_push(ecr_uri)
    ce_arn = provision_compute_env(sess, roles)
    provision_job_queue(sess, ce_arn)
    provision_job_definition(sess, f"{ecr_uri}:{C.IMAGE_TAG}", roles)
    _log("done")


if __name__ == "__main__":
    main()

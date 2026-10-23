"""Tests for the AWS CLI wrappers (src/aws/cli.py): command shape and JSON
handling, with subprocess mocked so no AWS call is made."""

import json
from types import SimpleNamespace

import pytest

from src.aws import cli
from src.aws import config as C


def _fake_run(stdout="{}", returncode=0, stderr=""):
    calls = []

    def run(cmd, capture_output, text):
        calls.append(cmd)
        return SimpleNamespace(returncode=returncode, stdout=stdout, stderr=stderr)
    return run, calls


def test_submit_job_passes_payload_as_cli_input_json(monkeypatch):
    run, calls = _fake_run(stdout=json.dumps({"jobName": "x", "jobId": "abc-123"}))
    monkeypatch.setattr(cli.subprocess, "run", run)
    payload = {"jobName": "x", "jobQueue": C.JOB_QUEUE, "jobDefinition": C.JOB_DEFINITION,
               "containerOverrides": {"environment": [{"name": "NP", "value": "4"}]}}
    assert cli.batch_submit_job(payload) == "abc-123"
    cmd = calls[0]
    assert cmd[:3] == ["aws", "batch", "submit-job"]
    assert cmd[cmd.index("--region") + 1] == C.REGION
    assert cmd[cmd.index("--output") + 1] == "json"
    assert json.loads(cmd[cmd.index("--cli-input-json") + 1]) == payload


def test_describe_jobs_returns_jobs_list(monkeypatch):
    run, calls = _fake_run(stdout=json.dumps({"jobs": [{"jobId": "a", "status": "RUNNING"}]}))
    monkeypatch.setattr(cli.subprocess, "run", run)
    assert cli.batch_describe_jobs(["a", "b"]) == [{"jobId": "a", "status": "RUNNING"}]
    cmd = calls[0]
    assert cmd[:3] == ["aws", "batch", "describe-jobs"]
    assert cmd[cmd.index("--jobs") + 1:] == ["a", "b"]


def test_describe_jobs_empty_list_makes_no_call(monkeypatch):
    run, calls = _fake_run()
    monkeypatch.setattr(cli.subprocess, "run", run)
    assert cli.batch_describe_jobs([]) == []
    assert calls == []


def test_cli_failure_raises_with_stderr(monkeypatch):
    run, _ = _fake_run(stdout="", returncode=254, stderr="An error occurred (AccessDenied)")
    monkeypatch.setattr(cli.subprocess, "run", run)
    with pytest.raises(RuntimeError, match="AccessDenied"):
        cli.batch_describe_jobs(["a"])


def test_submit_and_poll_modules_do_not_import_boto3():
    import src.aws.submit as submit
    import src.aws.poll as poll
    for mod in (submit, poll):
        assert "boto3" not in mod.__dict__

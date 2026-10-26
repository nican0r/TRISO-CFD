"""Tests for src/aws/poll.py: a failing describe-jobs call must not end a
--watch loop (an expired `aws login` session killed an overnight watcher on
2026-10-24); a one-shot poll still raises."""

import sys

import pytest

from src.aws import poll


def _run_main(monkeypatch, argv, describe_results):
    """Drive poll.main() with describe() returning/raising from a script."""
    script = list(describe_results)

    def fake_describe(job_id):
        item = script.pop(0)
        if isinstance(item, Exception):
            raise item
        return item

    sleeps = []
    monkeypatch.setattr(poll, "describe", fake_describe)
    monkeypatch.setattr(poll.time, "sleep", lambda s: sleeps.append(s))
    monkeypatch.setattr(sys, "argv", ["poll", *argv])
    poll.main()
    return sleeps


def _job(status):
    return {"jobId": "j1", "jobName": "n", "status": status}


def test_watch_retries_after_cli_failure(monkeypatch, capsys):
    sleeps = _run_main(monkeypatch, ["j1", "--watch", "--interval", "7"],
                       [_job("RUNNING"),
                        RuntimeError("aws batch describe-jobs failed (exit 255): "
                                     "Your session has expired"),
                        _job("SUCCEEDED")])
    out = capsys.readouterr().out
    assert "retrying in 7 s" in out and "session has expired" in out
    assert '"status": "SUCCEEDED"' in out
    assert sleeps == [7, 7]      # after RUNNING, after the failure


def test_one_shot_failure_raises(monkeypatch):
    with pytest.raises(RuntimeError, match="expired"):
        _run_main(monkeypatch, ["j1"], [RuntimeError("expired")])

"""Tests for AWS Batch job sizing read from case manifests (src/aws/submit.py)."""

import json

from src.aws.submit import manifest_sizing


def _write(tmp_path, manifest):
    (tmp_path / "manifest.json").write_text(json.dumps(manifest))
    return tmp_path


def test_manifest_with_aws_sizing(tmp_path):
    d = _write(tmp_path, {"aws_vcpu": 4, "aws_memory_mb": 6000, "nranks": 4,
                          "nodes_i": 2, "nodes_j": 1, "nodes_k": 2})
    assert manifest_sizing(d) == {"vcpu": 4, "memory": 6000, "nranks": 4,
                                  "nodes": (2, 1, 2)}


def test_manifest_without_aws_vcpu_is_ignored(tmp_path):
    # Older (V_3d) manifests carry nranks/nodes for local runs only; they must
    # not change the AWS defaults.
    d = _write(tmp_path, {"nranks": 6, "nodes_i": 1, "nodes_j": 6, "nodes_k": 1})
    assert manifest_sizing(d) == {}


def test_missing_manifest(tmp_path):
    assert manifest_sizing(tmp_path) == {}

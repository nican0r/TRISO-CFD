"""Outcome rules of src/post/piv_rig_score.py against params/piv_rig_acceptance.yaml."""
import sys
from pathlib import Path

import numpy as np
import pytest
import yaml

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src/post"))
from piv_rig_score import health_failures, outcome, source_cells  # noqa: E402

ACC = yaml.safe_load((REPO / "params/piv_rig_acceptance.yaml").read_text())


def _top(mm):
    return {"dense_top_mm": mm}


@pytest.mark.parametrize("top, d20, trend, want", [
    (20.0, 0.10, 0.05, "held"),
    (8.0, 0.49, 0.00, "held"),
    (2.0, 0.95, 0.00, "collapsed"),
    (2.0, 0.76, 0.00, "internal_jet"),      # C26 c01 TGS floor 0.50: a jet with no fountain
    (5.0, 0.60, 0.10, "in_between"),
    (20.0, 0.60, 0.10, "in_between"),       # a fountain, but the spout mostly closed
    (20.0, 0.10, 0.40, "in_between (still changing)"),
])
def test_step_down_outcomes(top, d20, trend, want):
    assert outcome(ACC["step_down_20"], "step_down_20", _top(top), d20, trend) == want


@pytest.mark.parametrize("name", ["run_r01_3s", "run_r02_3s", "run_r03_3s",
                                  "run_q1_3s", "run_q2_3s", "run_q3_3s"])
@pytest.mark.parametrize("top, d20, want", [
    (8.0, 0.14, "SPOUTS"),
    (3.0, 0.04, "NO SPOUT"),                # C26 c02 BVK: spout open, no dense fountain
    (2.0, 1.00, "NO SPOUT"),                # C26 c01 BVK
    (43.0, 0.60, "NO SPOUT"),
])
def test_3s_runs_spout(name, top, d20, want):
    assert outcome(ACC[name], name, _top(top), d20, None) == want


@pytest.mark.parametrize("top, d20, want", [
    (3.0, 0.50, "REGIME PASS"),
    (2.9, 0.50, "REGIME FAIL"),
    (20.0, 0.90, "REGIME FAIL"),
])
def test_smoke_regime(top, d20, want):
    assert outcome(ACC["smoke_r01"], "smoke_r01", _top(top), d20, None) == want


def test_missing_frames_are_not_a_result():
    assert outcome(ACC["smoke_r01"], "smoke_r01", None, 0.1, None).startswith("N/A")
    assert outcome(ACC["run_r01_3s"], "run_r01_3s", _top(20.0), None, None).startswith("N/A")


def test_acceptance_flow_targets():
    assert ACC["measured_min_spout"]["flow_lpm"] == 22
    assert ACC["run_r01_3s"]["flow_slpm"] == 24
    assert ACC["run_r02_3s"]["flow_slpm"] == 26.5
    assert ACC["run_r03_3s"]["flow_slpm"] == 30
    assert ACC["step_down_20"]["flow_slpm"] == 20
    assert ACC["step_down_20"]["flow_slpm"] < ACC["measured_min_spout"]["n2_equiv_slpm"]
    # the rig's test flows are the sweep's q cases, all above the measured limit
    sweep = yaml.safe_load((REPO / "params/sweep_piv_rig.yaml").read_text())
    flows = {c["id"]: c["flow_slpm"] for c in sweep["cases"]}
    for cid in ("q1", "q2", "q3"):
        assert ACC[f"run_{cid}_3s"]["flow_slpm"] == flows[cid]
        assert flows[cid] > ACC["measured_min_spout"]["n2_equiv_slpm"]


def test_3s_runs_share_limits_and_windows():
    r01 = ACC["run_r01_3s"]
    for name in ("run_r02_3s", "run_r03_3s", "run_q1_3s", "run_q2_3s", "run_q3_3s"):
        for key in ("window_s", "d20_window_s", "health_window_s",
                    "spouts_dense_top_min_mm", "spouts_d20_max"):
            assert ACC[name][key] == r01[key]


def test_step_down_windows_follow_the_hold_from_the_smoke():
    step, start = ACC["step_down_20"], ACC["smoke_r01"]["window_s"][1]
    end = start + step["hold_s"]
    assert step["window_s"] == [pytest.approx(end - 0.6), pytest.approx(end)]
    assert step["d20_window_s"] == [pytest.approx(end - 1.0), pytest.approx(end)]
    assert step["health_window_s"] == [pytest.approx(start + 0.1), pytest.approx(end)]


def _inlet(ratio=0.97, dp_pa=1077.0):
    return {"ratio": ratio, "dp_pa": dp_pa}


def _inv(first=27023, last=27012, y_min_mm=0.59):
    return {"first": first, "last": last, "y_min_mm": y_min_mm}


@pytest.mark.parametrize("inlet, inv, want", [
    (_inlet(), _inv(), []),                                  # the r01 smoke
    (_inlet(ratio=0.27), _inv(), ["inlet_ratio"]),
    (_inlet(ratio=4.28), _inv(), ["inlet_ratio"]),           # Day-10 R4's clogged throat
    (_inlet(dp_pa=-1483.0), _inv(), ["pressure"]),           # wrong sign is still a FAIL
    (_inlet(dp_pa=3.1 * 513), _inv(), ["pressure"]),
    (_inlet(dp_pa=0.1 * 513), _inv(), ["pressure"]),
    (_inlet(), _inv(first=20000, last=20000), ["seeded"]),
    (_inlet(), _inv(last=26800), ["parcel_loss"]),
    (_inlet(), _inv(y_min_mm=-0.5), ["no_bore_entry"]),
])
def test_health_failures(inlet, inv, want):
    assert health_failures(ACC["health"], inlet, inv) == want


def test_missing_frames_are_not_a_pass():
    assert health_failures(ACC["health"], None, _inv()) == ["N/A (missing frames)"]
    assert health_failures(ACC["health"], _inlet(), None) == ["N/A (missing frames)"]


DX = 2.428571e-3


def test_source_cells_skip_the_cut_cells_beside_the_axis():
    xyz = np.array([[0.0, -0.5 * DX, 0.0],       # source
                    [0.0, 0.5 * DX, 0.0],        # above it
                    [1.53e-3, 0.5 * DX, 0.0],    # cut cell beside it, as on the rig
                    [0.0, 1.5 * DX, 0.0]])       # the row after
    source, above = source_cells(xyz, DX)
    assert source.tolist() == [True, False, False, False]
    assert above.tolist() == [False, True, False, False]


def test_source_cells_refuse_a_missing_cell():
    with pytest.raises(SystemExit):
        source_cells(np.array([[0.0, 0.5 * DX, 0.0]]), DX)

"""Unit tests for the scoring helpers of notes/c26_watch_post.py."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest
import yaml

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "notes"))
import c26_watch_post as W  # noqa: E402

CRIT = yaml.safe_load((REPO / "params" / "c26_acceptance.yaml").read_text())


def test_fwhm_triangle():
    # v = 1 - |x| / 4 mm on 2.1 mm bins: half value at |x| = 2 mm -> 4 mm wide
    x = np.arange(-4, 5) * 2.1e-3
    v = np.clip(1 - np.abs(x) / 4e-3, -0.2, None)
    assert W.fwhm(x, v, 1.0) == pytest.approx(4e-3, rel=1e-9)


def test_fwhm_never_drops_or_no_spout():
    x = np.arange(-2, 3) * 2.1e-3
    assert math.isnan(W.fwhm(x, np.ones(5), 1.0))          # wider than the probe
    assert math.isnan(W.fwhm(x, -np.ones(5), -1.0))        # axis moving down


def test_ab_verdict_thresholds():
    ab = CRIT["grid_layout_ab"]
    assert W.ab_verdict(0.69, 1.0, ab)[1] == "substantial"     # -31 %
    assert W.ab_verdict(0.70, 1.0, ab)[1] == "substantial"     # -30 % exactly
    assert W.ab_verdict(0.80, 1.0, ab)[1] == "inconclusive"    # -20 %
    assert W.ab_verdict(0.90, 1.0, ab)[1] == "no effect"       # -10 %
    assert W.ab_verdict(1.10, 1.0, ab)[1] == "no effect"       # +10 %
    assert W.ab_verdict(1.20, 1.0, ab)[1] == "inconclusive"    # +20 %


def test_piv_match_uses_the_worse_c26_score():
    piv = {"CC_table3": 0.68, "CC_c26model_same_pipeline": 0.39,
           "SSI_table3": 0.61, "SSI_c26model_same_pipeline": 0.74,
           "MSE_m2_s2_table3": 0.0068, "MSE_m2_s2_c26model_same_pipeline": 0.0057,
           "CC": 0.40, "SSI": 0.62, "MSE_m2_s2": 0.0060}
    ok, floor = W.piv_match(piv)
    assert floor == {"CC": 0.39, "SSI": 0.61, "MSE": 0.0068}
    assert all(ok.values())
    ok, _ = W.piv_match({**piv, "SSI": 0.60})
    assert ok == {"CC": True, "SSI": False, "MSE": True}


def test_job_id_for_takes_the_last_aws_row(tmp_path):
    log = tmp_path / "run_log.csv"
    log.write_text(
        "date,step,case,key_params,mesh,dt,files_changed,tests_passed,tests_failed,wall_clock_s,status,notes\n"
        "d,aws-batch,c01_x,,,,,,,,submitted,aws_job_id=11111111-1111-1111-1111-111111111111; s3=a\n"
        "d,case-c26,c01_x,,,,,,,,PASS,aws_job_id=22222222-2222-2222-2222-222222222222\n"
        "d,aws-batch,c01_x,,,,,,,,submitted,aws_job_id=33333333-3333-3333-3333-333333333333; s3=b\n")
    assert W.job_id_for("c01_x", log) == "33333333-3333-3333-3333-333333333333"


def test_w_over_a_matches_check_smoke_c01():
    # check_smoke printed W/A = 413 Pa for c01 (60 g, bed-top r 21.3 mm)
    m = {"R_t_m": 1.19e-3, "H_static_m": 34.8e-3, "cone_half_angle_deg": 30.0,
         "R_c_m": 25.4e-3, "inventory_kg": 0.06}
    assert W.w_over_a(m) == pytest.approx(413, abs=1)


def test_band_absorbs_rounding_only():
    assert W.band(9.999999999999991, (10.0, 20.0)) == "PASS"   # PSD bin at 10 Hz
    assert W.band(9.99, (10.0, 20.0)) == "FAIL"
    assert W.band(20.0, (10.0, 20.0)) == "PASS"

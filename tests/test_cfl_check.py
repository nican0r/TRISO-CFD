"""Unit tests for src/cfl_check.py.

Pinned to hand-worked numbers, per CLAUDE.md rule 4.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from src.cfl_check import (  # noqa: E402
    _load_spout_example,
    cfl_number,
    dt_for_cfl,
)


def test_cfl_number_definition_identity():
    # CFL = u * dt / dx.  Pick u=dx/dt so CFL == 1.
    assert cfl_number(1.0, 1.0, 1.0) == pytest.approx(1.0)
    # u=10 m/s, dx=1e-3 m, dt=1e-4 s -> CFL = 1.0
    assert cfl_number(10.0, 1.0e-3, 1.0e-4) == pytest.approx(1.0)


def test_cfl_number_half_cell_step():
    # dt half of the "just-CFL=1" step should give CFL = 0.5.
    assert cfl_number(20.0, 1.0e-3, 2.5e-5) == pytest.approx(0.5)


def test_dt_for_cfl_inverts_cfl_number():
    # For any (u, dx, target), CFL(u, dx, dt_for_cfl(u, dx, target)) == target.
    u, dx, target = 20.0, 1.0e-3, 0.5
    dt = dt_for_cfl(u, dx, target)
    assert cfl_number(u, dx, dt) == pytest.approx(target, rel=1e-12)


def test_spout_example_dt_hand_calc():
    """Day-3 hand calc — spout example.

    u_max = 20 m/s, dx = 1 mm, CFL <= 0.5
        dt <= CFL * dx / u_max = 0.5 * 1e-3 / 20 = 2.5e-5 s.
    """
    p = _load_spout_example()
    u_max = float(p["u_max_m_per_s"])
    dx = float(p["dx_m"])
    cfl_target = float(p["cfl_target"])
    assert u_max == 20.0
    assert dx == 1.0e-3
    assert cfl_target == 0.5

    dt_max = dt_for_cfl(u_max, dx, cfl_target)
    assert dt_max == pytest.approx(2.5e-5, rel=1e-12)
    # And the corresponding CFL is exactly the target (< 0.5 satisfied at equality).
    assert cfl_number(u_max, dx, dt_max) == pytest.approx(0.5, rel=1e-12)


def test_spout_example_smaller_dt_gives_smaller_cfl():
    # A 10x smaller dt should give a 10x smaller CFL: sanity of the linear relation.
    p = _load_spout_example()
    u_max = float(p["u_max_m_per_s"])
    dx = float(p["dx_m"])
    dt_big = dt_for_cfl(u_max, dx, 0.5)
    dt_small = dt_big / 10.0
    assert cfl_number(u_max, dx, dt_small) == pytest.approx(0.05, rel=1e-12)


@pytest.mark.parametrize(
    "u, dx, dt",
    [
        (0.0, 1.0e-3, 1.0e-5),
        (-1.0, 1.0e-3, 1.0e-5),
        (10.0, 0.0, 1.0e-5),
        (10.0, -1e-3, 1.0e-5),
        (10.0, 1.0e-3, 0.0),
        (10.0, 1.0e-3, -1.0e-5),
    ],
)
def test_cfl_number_rejects_nonpositive(u, dx, dt):
    with pytest.raises(ValueError):
        cfl_number(u, dx, dt)


@pytest.mark.parametrize(
    "u, dx, target",
    [
        (0.0, 1.0e-3, 0.5),
        (10.0, 0.0, 0.5),
        (10.0, 1.0e-3, 0.0),
    ],
)
def test_dt_for_cfl_rejects_nonpositive(u, dx, target):
    with pytest.raises(ValueError):
        dt_for_cfl(u, dx, target)

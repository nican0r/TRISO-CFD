"""CFL sanity-check helpers for MFiX time-step selection.

The 1-D Courant-Friedrichs-Lewy number is

    CFL = u_max * dt / dx                                       (Eq. 1)

with u_max the peak flow (or signal) speed in the region of interest [m/s],
dx the characteristic cell edge length [m], and dt the solver time step [s].
For explicit-in-time transport, CFL <= 1 is a hard stability constraint;
CFL <= 0.5 is a common safety margin used to sanity-check MFiX runs before
launch.

All quantities SI (CLAUDE.md rule 3).  Inputs live in params/numerics.yaml
(CLAUDE.md rule 2); this module contains no hard-coded physics values.
"""

from __future__ import annotations

from pathlib import Path

import yaml


PARAMS_FILE = (
    Path(__file__).resolve().parent.parent / "params" / "numerics.yaml"
)


def cfl_number(u_max: float, dx: float, dt: float) -> float:
    """Return the CFL number for a given velocity, cell size and time step.

    u_max : peak velocity magnitude in the cell [m/s]  (must be > 0)
    dx    : characteristic cell edge length [m]        (must be > 0)
    dt    : solver time step [s]                       (must be > 0)
    Returns CFL [-].
    """
    if u_max <= 0.0:
        raise ValueError(f"u_max must be > 0 m/s, got {u_max}")
    if dx <= 0.0:
        raise ValueError(f"dx must be > 0 m, got {dx}")
    if dt <= 0.0:
        raise ValueError(f"dt must be > 0 s, got {dt}")
    return u_max * dt / dx


def dt_for_cfl(u_max: float, dx: float, cfl_target: float) -> float:
    """Largest dt that satisfies CFL <= cfl_target for the given u_max and dx.

    Inverts Eq. 1:  dt = cfl_target * dx / u_max.

    u_max      : peak velocity magnitude [m/s]         (must be > 0)
    dx         : characteristic cell edge length [m]   (must be > 0)
    cfl_target : upper bound on CFL [-]                (must be > 0)
    Returns dt [s].
    """
    if u_max <= 0.0:
        raise ValueError(f"u_max must be > 0 m/s, got {u_max}")
    if dx <= 0.0:
        raise ValueError(f"dx must be > 0 m, got {dx}")
    if cfl_target <= 0.0:
        raise ValueError(f"cfl_target must be > 0, got {cfl_target}")
    return cfl_target * dx / u_max


def _load_spout_example(path: Path = PARAMS_FILE) -> dict:
    with open(path, "r") as f:
        return yaml.safe_load(f)["spout_cfl_example"]


def _report_spout_example() -> None:
    p = _load_spout_example()
    u_max = float(p["u_max_m_per_s"])
    dx = float(p["dx_m"])
    cfl_target = float(p["cfl_target"])
    dt_max = dt_for_cfl(u_max, dx, cfl_target)
    cfl_at_dt_max = cfl_number(u_max, dx, dt_max)
    print("Spout CFL sanity-check (params/numerics.yaml):")
    print(f"  u_max      = {u_max:.3g} m/s")
    print(f"  dx         = {dx:.3g} m")
    print(f"  CFL target = {cfl_target:.3g}")
    print(f"  => dt_max  = {dt_max:.3e} s   (CFL = {cfl_at_dt_max:.3f})")


if __name__ == "__main__":
    _report_spout_example()

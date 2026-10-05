"""Cantera-based gas property helpers.

Uses the gri30 mechanism (bundled with Cantera). gri30 includes AR, H2, N2, O2
with transport data — sufficient for the coater cold/hot cases in PLAN.md 1.3.
"""

from __future__ import annotations

from typing import Mapping, Tuple

import cantera as ct

_MECH = "gri30.yaml"
_gas: ct.Solution | None = None


def _solution() -> ct.Solution:
    global _gas
    if _gas is None:
        _gas = ct.Solution(_MECH)
    return _gas


def gas_props(
    T_K: float, P_Pa: float, composition: Mapping[str, float]
) -> Tuple[float, float]:
    """Density and dynamic viscosity of a gas mixture at (T, P).

    T_K         : temperature [K]
    P_Pa        : pressure [Pa]
    composition : mapping species -> mole fraction. Species names are
                  case-insensitive here and are upper-cased before being passed
                  to Cantera (gri30 uses e.g. AR, H2, N2, O2).
                  Fractions should sum to 1 but Cantera will normalize.

    Returns (rho [kg/m^3], mu [Pa*s]).
    """
    gas = _solution()
    normalized = {str(k).upper(): float(v) for k, v in composition.items()}
    gas.TPX = float(T_K), float(P_Pa), normalized
    return float(gas.density), float(gas.viscosity)

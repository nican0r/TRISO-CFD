"""Fluid-particle correlations used across the TRISO coater CFD project.

All quantities are SI. Every physical argument carries its unit in the docstring.
CLAUDE.md rule 3 (SI everywhere) and rule 4 (every correlation gets a unit test).
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

G_STD = 9.80665  # m/s^2, standard gravity (CODATA)


def reynolds(u: float, d_p: float, rho: float, mu: float) -> float:
    """Particle Reynolds number.

    u   : slip / relative velocity between particle and fluid [m/s]
    d_p : particle diameter [m]
    rho : fluid density [kg/m^3]
    mu  : fluid dynamic viscosity [Pa*s]
    Returns Re [-].
    """
    return rho * abs(u) * d_p / mu


def cd_schiller_naumann(re):
    """Drag coefficient for an isolated sphere, Schiller-Naumann with Newton cap.

    Cd = 24/Re * (1 + 0.15 Re^0.687)   for  0 < Re <= 1000
    Cd = 0.44                           for  Re >  1000   (Newton regime cap)

    Reference: Schiller & Naumann (1935); Newton-regime constant Cd = 0.44 is the
    standard value cited e.g. in Crowe et al., "Multiphase Flows with Droplets
    and Particles" (2nd ed., 2011), Ch. 4.

    re : Reynolds number [-] (scalar or array-like)
    Returns Cd [-].
    """
    re_arr = np.atleast_1d(np.asarray(re, dtype=float))
    cd = np.empty_like(re_arr)
    zero = re_arr <= 0.0
    high = re_arr > 1000.0
    low = ~(zero | high)
    cd[zero] = np.inf
    cd[high] = 0.44
    cd[low] = 24.0 / re_arr[low] * (1.0 + 0.15 * re_arr[low] ** 0.687)
    return cd.item() if np.ndim(re) == 0 else cd


def terminal_velocity(
    d_p: float,
    rho_p: float,
    rho: float,
    mu: float,
    g: float = G_STD,
) -> float:
    """Terminal settling velocity of a single sphere in a quiescent fluid.

    Force balance at terminal condition (drag = net weight):
        0.5 * rho * u_t^2 * Cd(Re) * (pi/4 d_p^2) = (rho_p - rho) * (pi/6 d_p^3) * g
    which reduces to
        u_t^2 * Cd(Re) = (4/3) * d_p * (rho_p - rho) * g / rho .
    Solved for u_t with scipy.optimize.brentq on the interval [1e-8, 100] m/s
    using Schiller-Naumann Cd (with Newton cap at Re > 1000).

    d_p   : particle diameter [m]
    rho_p : particle density [kg/m^3]
    rho   : fluid density [kg/m^3]
    mu    : fluid dynamic viscosity [Pa*s]
    g     : gravitational acceleration [m/s^2]
    Returns u_t [m/s]. Returns 0.0 if the particle is neutrally or negatively buoyant.
    """
    if rho_p <= rho:
        return 0.0
    rhs = (4.0 / 3.0) * d_p * (rho_p - rho) * g / rho  # = u_t^2 * Cd at terminal

    def residual(u: float) -> float:
        re = reynolds(u, d_p, rho, mu)
        return u * u * cd_schiller_naumann(re) - rhs

    return brentq(residual, 1.0e-8, 100.0, xtol=1.0e-10)


def relaxation_time(d_p: float, rho_p: float, mu: float) -> float:
    """Stokes particle relaxation time.

        tau_p = rho_p * d_p^2 / (18 * mu)   [s]

    Valid strictly in the Stokes limit (Re_p << 1); serves as a reference
    time scale for higher-Re particles as well.

    d_p   : particle diameter [m]
    rho_p : particle density [kg/m^3]
    mu    : fluid dynamic viscosity [Pa*s]
    Returns tau_p [s].
    """
    return rho_p * d_p * d_p / (18.0 * mu)


def archimedes(
    d_p: float, rho_p: float, rho: float, mu: float, g: float = G_STD
) -> float:
    """Archimedes number for a fluid/particle pair.

        Ar = rho * (rho_p - rho) * g * d_p^3 / mu^2   [-]

    Reference: Kunii & Levenspiel, "Fluidization Engineering" (2nd ed., 1991),
    Ch. 3, definition used with the Wen-Yu U_mf correlation.

    d_p   : particle diameter [m]
    rho_p : particle density [kg/m^3]
    rho   : fluid density [kg/m^3]
    mu    : fluid dynamic viscosity [Pa*s]
    g     : gravitational acceleration [m/s^2]
    Returns Ar [-].
    """
    return rho * (rho_p - rho) * g * d_p ** 3 / (mu * mu)


def ergun_dp(
    u: float, d_p: float, rho: float, mu: float, epsilon: float
) -> float:
    """Ergun pressure gradient across a packed bed (per unit bed height).

        dP/L = 150 * mu * (1-eps)^2 * U / (eps^3 * d_p^2)
             + 1.75 * rho * (1-eps) * U^2 / (eps^3 * d_p)          [Pa/m]

    Reference: Ergun, S. (1952), "Fluid flow through packed columns",
    Chem. Eng. Prog. 48(2), 89-94.

    u       : superficial gas velocity [m/s]
    d_p     : particle diameter [m]
    rho     : gas density [kg/m^3]
    mu      : gas dynamic viscosity [Pa*s]
    epsilon : bed voidage [-], 0 < eps < 1
    Returns dP/L [Pa/m] (positive for upward flow through a static bed).
    """
    if not 0.0 < epsilon < 1.0:
        raise ValueError(f"epsilon must be in (0, 1), got {epsilon}")
    one_minus = 1.0 - epsilon
    visc = 150.0 * mu * one_minus * one_minus * u / (epsilon ** 3 * d_p * d_p)
    inert = 1.75 * rho * one_minus * u * u / (epsilon ** 3 * d_p)
    return visc + inert


def umf_wen_yu(
    d_p: float, rho_p: float, rho: float, mu: float, g: float = G_STD
) -> float:
    """Minimum fluidization velocity via the Wen & Yu (1966) correlation.

        Re_mf = sqrt(33.7^2 + 0.0408 * Ar) - 33.7
        U_mf  = Re_mf * mu / (rho * d_p)                             [m/s]

    Reference: Wen, C.Y. & Yu, Y.H. (1966), "A generalized method for
    predicting the minimum fluidization velocity", AIChE J. 12(3), 610-612.

    d_p   : particle diameter [m]
    rho_p : particle density [kg/m^3]
    rho   : gas density [kg/m^3]
    mu    : gas dynamic viscosity [Pa*s]
    g     : gravitational acceleration [m/s^2]
    Returns U_mf [m/s].
    """
    ar = archimedes(d_p, rho_p, rho, mu, g)
    re_mf = np.sqrt(33.7 * 33.7 + 0.0408 * ar) - 33.7
    return re_mf * mu / (rho * d_p)


def tgs_drag_factor(epsilon: float, re: float) -> float:
    """Tenneti-Garg-Subramaniam dimensionless drag force F(phi, Re).

        F = F_isol / eps^3 + F_phi + F_phi_Re
        F_isol   = 1 + 0.15 Re^0.687
        F_phi    = 5.81 phi / eps^3 + 0.48 phi^(1/3) / eps^4
        F_phi_Re = phi^3 Re (0.95 + 0.61 phi^3 / eps^2)
        phi      = 1 - eps

    Reference: Tenneti, Garg & Subramaniam (2011), "Drag law for monodisperse
    gas-solid systems using particle-resolved direct numerical simulation of
    flow past fixed assemblies of spheres", Int. J. Multiphase Flow 37(9),
    1072-1092. Coded as in MFiX 26.1.2 `model/drag_gs.f` DRAG_TGS.

    epsilon : gas volume fraction [-], in (0, 1]
    re      : slip Reynolds number eps * rho * |u_g - v_p| * d_p / mu [-]
    Returns F [-]; F = 1 + 0.15 Re^0.687 for an isolated sphere (eps = 1).
    """
    if not 0.0 < epsilon <= 1.0:
        raise ValueError(f"epsilon must be in (0, 1], got {epsilon}")
    phi = 1.0 - epsilon
    f = (1.0 + 0.15 * re ** 0.687) / epsilon ** 3
    f += 5.81 * phi / epsilon ** 3 + 0.48 * phi ** (1.0 / 3.0) / epsilon ** 4
    f += phi ** 3 * re * (0.95 + 0.61 * phi ** 3 / epsilon ** 2)
    return f


def tgs_drag_coefficient(
    epsilon: float, u_slip: float, d_p: float, rho: float, mu: float
) -> float:
    """TGS drag coefficient per unit particle volume, as MFiX DES applies it.

        DgA = 18 mu eps / d_p^2 * F(phi, Re),   Re = eps rho |u_slip| d_p / mu

    The drag force on one particle is DgA * V_p * u_slip (MFiX model A,
    `des/drag_gp_des.f`: F_gp = DgA * PVOL). In the isolated-sphere limit this is
    3 pi mu d_p u_slip (1 + 0.15 Re^0.687), i.e. Schiller-Naumann.

    epsilon : gas volume fraction [-]
    u_slip  : interstitial gas velocity minus particle velocity [m/s]
    d_p     : particle diameter used by the drag law [m]
    rho     : gas density [kg/m^3]
    mu      : gas dynamic viscosity [Pa*s]
    Returns DgA [kg/(m^3 s)].
    """
    re = epsilon * reynolds(u_slip, d_p, rho, mu)
    return 18.0 * mu * epsilon / (d_p * d_p) * tgs_drag_factor(epsilon, re)


def umf_tgs(
    d_p: float, rho_p: float, rho: float, mu: float, epsilon: float,
    g: float = G_STD,
) -> float:
    """Minimum fluidization velocity of a uniform bed under the TGS drag law.

    Uniform suspended bed in MFiX DES model A (particles at rest, gas
    interstitial velocity U/eps). Gas momentum: eps G = eps rho g + n F_gp u_g;
    particle: F_gp u_g = V_p (rho_p g - G), with G = -dP/dy. Eliminating G
    gives the per-particle-volume balance

        DgA(eps, U/eps) * U / eps = eps (rho_p - rho) g

    solved for the superficial velocity U with brentq.

    d_p     : particle diameter used by the drag law [m]
    rho_p   : particle density [kg/m^3]
    rho     : gas density [kg/m^3]
    mu      : gas dynamic viscosity [Pa*s]
    epsilon : bed voidage at minimum fluidization [-]
    g       : gravitational acceleration [m/s^2]
    Returns U_mf [m/s] (superficial).
    """
    rhs = epsilon * (rho_p - rho) * g

    def residual(u: float) -> float:
        u_g = u / epsilon
        return tgs_drag_coefficient(epsilon, u_g, d_p, rho, mu) * u_g - rhs

    return brentq(residual, 1.0e-9, 100.0, xtol=1.0e-12)


def geldart_class(d_p: float, rho_p: float, rho: float) -> str:
    """Rule-based Geldart (1973) powder classifier.

    Boundaries used (documented in the code below; SI units throughout):

        Group C  : cohesive fines,      d_p < 30 um
        Group D  : spoutable/large,     d_p^2 * (rho_p - rho) > 1e-3 kg/m
        Group A  : aeratable fines,     d_p * (rho_p - rho)   < 0.225 kg/m^2
        Group B  : bubbling / sand-like otherwise

    Check order: C -> D -> A -> B. The A/B and B/D lines come directly from
    Geldart's original chart lines (approx. 225 um at 1 g/cm^3 for A/B, and
    1000 um at 1 g/cm^3 for B/D). See:
        Geldart, D. (1973), "Types of gas fluidization",
        Powder Technol. 7, 285-292.

    d_p   : particle diameter [m]
    rho_p : particle density [kg/m^3]
    rho   : gas density [kg/m^3]
    Returns one of 'A', 'B', 'C', 'D'.
    """
    delta_rho = rho_p - rho
    if delta_rho <= 0.0:
        raise ValueError("rho_p must exceed rho for Geldart classification")
    if d_p < 30.0e-6:
        return "C"
    if d_p * d_p * delta_rho > 1.0e-3:
        return "D"
    if d_p * delta_rho < 0.225:
        return "A"
    return "B"


def ums_mathur_gishler(
    d_p: float,
    rho_p: float,
    rho: float,
    D_c: float,
    D_i: float,
    H: float,
    g: float = G_STD,
) -> float:
    """Minimum spouting velocity, Mathur & Gishler (1955) correlation.

        U_ms = (d_p / D_c) * (D_i / D_c)^(1/3)
             * sqrt( 2 * g * H * (rho_p - rho) / rho )               [m/s]

    Reference: Mathur, K.B. & Gishler, P.E. (1955), "A technique for
    contacting gases with coarse solid particles", AIChE J. 1(2), 157-164;
    also reproduced in Mathur & Epstein, "Spouted Beds" (1974), Eq. 2.1.

    d_p   : particle diameter [m]
    rho_p : particle density [kg/m^3]
    rho   : gas density [kg/m^3]
    D_c   : column (cylindrical section) diameter [m]
    D_i   : inlet-orifice diameter [m]
    H     : static (settled) bed height [m]
    g     : gravitational acceleration [m/s^2]
    Returns U_ms [m/s], defined on the column cross-section (pi/4 * D_c^2).
    """
    return (d_p / D_c) * (D_i / D_c) ** (1.0 / 3.0) * np.sqrt(
        2.0 * g * H * (rho_p - rho) / rho
    )


def coated_particle_stages(
    d_core_m: float,
    rho_core_kg_m3: float,
    layers,
) -> list[tuple[str, float, float]]:
    """Diameter and mean density after each coating layer of a layered sphere.

    Each layer is a concentric shell of uniform thickness and density, so the
    particle after layer k has

        d_k   = d_(k-1) + 2 t_k                                          [m]
        rho_k = (rho_(k-1) d_(k-1)^3 + rho_layer,k (d_k^3 - d_(k-1)^3)) / d_k^3
                                                                    [kg/m^3]

    (mass-weighted by shell volume).  This is how the stage densities of a
    TRISO particle follow from its layer specification, e.g. Marshall et al.
    (2017), Table 1 (425 um kernel -> 855 um, 3.03 g/cm^3).

    d_core_m       : core (kernel) diameter [m]
    rho_core_kg_m3 : core density [kg/m^3]
    layers         : iterable of (name, thickness [m], density [kg/m^3]) in
                     deposition order
    Returns [(name, d [m], rho [kg/m^3]), ...] starting with ("core", d_core,
    rho_core) and then one entry per layer, cumulative.
    """
    if d_core_m <= 0.0 or rho_core_kg_m3 <= 0.0:
        raise ValueError("core diameter and density must be > 0")
    d = float(d_core_m)
    mass_d3 = float(rho_core_kg_m3) * d ** 3   # rho * d^3 is proportional to mass
    out = [("core", d, float(rho_core_kg_m3))]
    for name, t_m, rho_layer in layers:
        if t_m <= 0.0 or rho_layer <= 0.0:
            raise ValueError(f"layer {name!r}: thickness and density must be > 0")
        d_new = d + 2.0 * float(t_m)
        mass_d3 += float(rho_layer) * (d_new ** 3 - d ** 3)
        d = d_new
        out.append((str(name), d, mass_d3 / d ** 3))
    return out

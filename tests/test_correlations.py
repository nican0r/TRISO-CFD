"""Unit tests for src/correlations.py.

Each test pins a correlation to a hand-worked or published number, per CLAUDE.md rule 4.
"""

import math
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from src.correlations import (  # noqa: E402
    G_STD,
    archimedes,
    cd_schiller_naumann,
    coated_particle_stages,
    ergun_dp,
    geldart_class,
    relaxation_time,
    reynolds,
    terminal_velocity,
    tgs_drag_coefficient,
    tgs_drag_factor,
    umf_tgs,
    umf_wen_yu,
    ums_mathur_gishler,
)


def test_reynolds_definition():
    # Re = rho * u * d / mu. Pick numbers that give Re=1.
    assert reynolds(1.0, 1.0, 1.0, 1.0) == pytest.approx(1.0)
    # Air-like: rho=1.2, u=5 m/s, d=500e-6 m, mu=1.8e-5 Pa*s -> Re = 166.667
    assert reynolds(5.0, 500e-6, 1.2, 1.8e-5) == pytest.approx(166.6667, rel=1e-4)


def test_cd_schiller_naumann_at_re_1():
    # Cd(Re=1) = 24 * (1 + 0.15) = 27.6
    assert cd_schiller_naumann(1.0) == pytest.approx(27.6, rel=1e-6)


def test_cd_schiller_naumann_newton_cap():
    # Above Re = 1000 the correlation is capped at Cd = 0.44.
    assert cd_schiller_naumann(1000.001) == 0.44
    assert cd_schiller_naumann(5000.0) == 0.44
    # Just below the cap we still follow Schiller-Naumann.
    assert cd_schiller_naumann(999.0) == pytest.approx(
        24.0 / 999.0 * (1 + 0.15 * 999.0 ** 0.687), rel=1e-9
    )


def test_cd_schiller_naumann_vectorized():
    cd = cd_schiller_naumann(np.array([1.0, 10.0, 100.0, 2000.0]))
    assert cd[0] == pytest.approx(27.6, rel=1e-6)
    assert cd[3] == 0.44
    # monotonic decrease over the SN branch
    assert cd[0] > cd[1] > cd[2]


def test_terminal_velocity_500um_6000kgm3_air_20C():
    """Hand-worked reference — see notes/day1_handcalc.md.

    Air at 20 C (from Cantera gri30, N2:0.79 / O2:0.21):
      rho = 1.199 kg/m^3, mu = 1.830e-5 Pa*s.
    500 um, 6000 kg/m^3 sphere -> u_t = 6.44 m/s (iterated by hand, +/- 0.05).
    """
    u_t = terminal_velocity(500e-6, 6000.0, 1.199, 1.830e-5)
    assert u_t == pytest.approx(6.44, abs=0.1)


def test_terminal_velocity_zero_when_neutrally_buoyant():
    assert terminal_velocity(500e-6, 1.0, 1.0, 1.8e-5) == 0.0
    assert terminal_velocity(500e-6, 0.5, 1.0, 1.8e-5) == 0.0


def test_terminal_velocity_recovers_stokes_for_tiny_particle():
    # For d_p = 20 um in water-like fluid the Stokes limit u_t = g d^2 (rho_p-rho) / (18 mu)
    # should agree with brentq solution to better than 1%.
    d_p = 20e-6
    rho_p, rho, mu = 2500.0, 1.0, 1.8e-5
    stokes = G_STD * d_p ** 2 * (rho_p - rho) / (18.0 * mu)
    u_t = terminal_velocity(d_p, rho_p, rho, mu)
    assert u_t == pytest.approx(stokes, rel=0.05)


def test_relaxation_time_stokes():
    # tau_p = rho_p d^2 / (18 mu). d=100 um, rho_p=2500, mu=1.8e-5:
    # tau = 2500 * 1e-8 / (18 * 1.8e-5) = 0.077160...
    tau = relaxation_time(1e-4, 2500.0, 1.8e-5)
    assert tau == pytest.approx(2500 * 1e-8 / (18 * 1.8e-5), rel=1e-9)
    assert math.isclose(tau, 0.077160, rel_tol=1e-3)


# ---------------------------------------------------------------------------
# Day 2 additions
# ---------------------------------------------------------------------------

def test_archimedes_500um_6000_air_20C():
    # 500 um, 6000 kg/m^3, air-like (rho=1.2, mu=1.8e-5):
    #   Ar = 1.2 * (6000-1.2) * 9.80665 * (5e-4)^3 / (1.8e-5)^2
    #      = 1.2 * 5998.8 * 9.80665 * 1.25e-10 / 3.24e-10
    #      ~= 27245  (hand calc)
    ar = archimedes(500e-6, 6000.0, 1.2, 1.8e-5)
    expected = 1.2 * (6000.0 - 1.2) * G_STD * (500e-6) ** 3 / (1.8e-5) ** 2
    assert ar == pytest.approx(expected, rel=1e-9)
    assert ar == pytest.approx(27245.0, rel=5e-3)


def test_ergun_dp_hand_calc():
    # eps=0.4, d_p=1e-3, U=0.1 m/s, rho=1.2, mu=1.8e-5.
    # viscous : 150 * 1.8e-5 * 0.36 * 0.1 / (0.064 * 1e-6) = 1518.75 Pa/m
    # inertial: 1.75 * 1.2 * 0.6 * 0.01 / (0.064 * 1e-3)   = 196.875 Pa/m
    # total   = 1715.625 Pa/m
    dpl = ergun_dp(0.1, 1.0e-3, 1.2, 1.8e-5, 0.4)
    assert dpl == pytest.approx(1715.625, rel=1e-9)


def test_ergun_dp_rejects_bad_voidage():
    with pytest.raises(ValueError):
        ergun_dp(0.1, 1.0e-3, 1.2, 1.8e-5, 0.0)
    with pytest.raises(ValueError):
        ergun_dp(0.1, 1.0e-3, 1.2, 1.8e-5, 1.0)


def test_umf_wen_yu_500um_6000_air_20C():
    # Same case as above. With Ar = 27245:
    #   Re_mf = sqrt(33.7^2 + 0.0408*27245) - 33.7
    #         = sqrt(1135.69 + 1111.60) - 33.7
    #         = sqrt(2247.29) - 33.7 = 47.406 - 33.7 = 13.706
    #   U_mf = 13.706 * 1.8e-5 / (1.2 * 5e-4) = 0.4112 m/s
    u_mf = umf_wen_yu(500e-6, 6000.0, 1.2, 1.8e-5)
    assert u_mf == pytest.approx(0.411, abs=0.01)


def test_umf_wen_yu_stokes_limit_small_ar():
    # In the small-Ar (Stokes) limit, sqrt(a^2 + x) - a -> x/(2a), so
    #   Re_mf -> 0.0408 * Ar / (2 * 33.7) = Ar / 1651.96
    # -> U_mf -> (rho_p - rho) * g * d_p^2 / (1651.96 * mu).
    # Pick tiny d_p to be sure we sit in the Stokes-like branch.
    d_p, rho_p, rho, mu = 30e-6, 2500.0, 1.2, 1.8e-5
    u_mf = umf_wen_yu(d_p, rho_p, rho, mu)
    stokes_like = (rho_p - rho) * G_STD * d_p ** 2 / (1651.96 * mu)
    assert u_mf == pytest.approx(stokes_like, rel=0.02)


def test_tgs_drag_factor_hand_calc():
    # eps = 0.42, phi = 0.58, Re = 10:
    #   F_isol/eps^3 = (1 + 0.15 * 10^0.687) / 0.42^3 = 1.729610 / 0.074088 = 23.3454
    #   5.81 phi/eps^3 = 3.3698 / 0.074088                               = 45.4838
    #   0.48 phi^(1/3)/eps^4 = 0.48 * 0.833955 / 0.0311170              = 12.8643
    #   phi^3 Re (0.95 + 0.61 phi^3/eps^2)
    #     = 0.195112 * 10 * (0.95 + 0.61 * 0.195112 / 0.1764)           =  3.1700
    #   F = 84.8635
    assert tgs_drag_factor(0.42, 10.0) == pytest.approx(84.8635, rel=1e-4)


def test_tgs_drag_factor_isolated_sphere_is_schiller_naumann():
    # phi = 0: every voidage term vanishes, leaving 1 + 0.15 Re^0.687
    assert tgs_drag_factor(1.0, 50.0) == pytest.approx(1.0 + 0.15 * 50.0 ** 0.687)


def test_tgs_drag_coefficient_isolated_sphere_force():
    # DgA * V_p * u at eps = 1 must equal Cd(Re) * 0.5 rho u^2 * pi d^2 / 4
    d_p, rho, mu, u = 411.6e-6, 1.633, 2.302e-5, 0.8
    force = tgs_drag_coefficient(1.0, u, d_p, rho, mu) * math.pi / 6.0 * d_p ** 3 * u
    re = reynolds(u, d_p, rho, mu)
    expected = cd_schiller_naumann(re) * 0.5 * rho * u * u * math.pi / 4.0 * d_p ** 2
    assert force == pytest.approx(expected, rel=1e-12)


def test_tgs_drag_factor_rejects_bad_voidage():
    with pytest.raises(ValueError):
        tgs_drag_factor(0.0, 1.0)
    with pytest.raises(ValueError):
        tgs_drag_factor(1.2, 1.0)


def test_umf_tgs_stokes_limit():
    # Re -> 0: F -> F0 = (1 + 5.81 phi)/eps^3 + 0.48 phi^(1/3)/eps^4, and
    #   18 mu eps F0 / d^2 * U/eps = eps (rho_p - rho) g
    #   -> U = eps (rho_p - rho) g d^2 / (18 mu F0)
    d_p, rho_p, rho, mu, eps = 10e-6, 2500.0, 1.2, 1.8e-5, 0.42
    phi = 1.0 - eps
    f0 = (1.0 + 5.81 * phi) / eps ** 3 + 0.48 * phi ** (1.0 / 3.0) / eps ** 4
    stokes = eps * (rho_p - rho) * G_STD * d_p ** 2 / (18.0 * mu * f0)
    assert umf_tgs(d_p, rho_p, rho, mu, eps) == pytest.approx(stokes, rel=1e-3)


def test_geldart_class_group_D_triso_surrogate():
    # 500 um ZrO2 surrogate at 6000 kg/m^3 in Ar (~1.66 kg/m^3) -> Group D
    assert geldart_class(500e-6, 6000.0, 1.66) == "D"


def test_geldart_class_group_D_uo2_kernel():
    # 500 um UO2 at 10 800 kg/m^3 -> Group D (even more clearly)
    assert geldart_class(500e-6, 10800.0, 1.66) == "D"


def test_geldart_class_group_A_fcc():
    # ~60 um FCC catalyst at 1500 kg/m^3 in air -> Group A
    #   d_p * dRho = 6e-5 * ~1499 = 0.0899 < 0.225 -> A
    assert geldart_class(60e-6, 1500.0, 1.2) == "A"


def test_geldart_class_group_B_sand():
    # 300 um sand at 2500 kg/m^3 in air -> Group B
    #   d_p^2 * dRho = 9e-8 * 2499 = 2.25e-4 < 1e-3 -> not D
    #   d_p   * dRho = 3e-4 * 2499 = 0.750 > 0.225 -> not A
    assert geldart_class(300e-6, 2500.0, 1.2) == "B"


def test_geldart_class_group_C_fines():
    # 10 um fines -> Group C by the d_p<30 um rule
    assert geldart_class(10e-6, 2000.0, 1.2) == "C"


def test_geldart_class_rejects_neutral_buoyancy():
    with pytest.raises(ValueError):
        geldart_class(500e-6, 1.0, 1.0)


def test_ums_mathur_gishler_hand_calc():
    # ORNL/UTK cold-mockup geometry: D_c=0.050, D_i=0.004.
    # Case D particles: d_p=500e-6 m, rho_p=6050 kg/m^3.
    # Ar at 20 C: rho ~= 1.66 kg/m^3.  H = 0.10 m.
    #
    #   d_p/D_c        = 500e-6 / 0.050     = 1.0000e-2
    #   (D_i/D_c)^(1/3)= (0.004/0.050)^(1/3) = (0.08)^(1/3) = 0.43089
    #   inside sqrt    = 2 * 9.80665 * 0.10 * (6050-1.66)/1.66
    #                  = 1.96133 * 3643.579 = 7147.01
    #   sqrt(...)      = 84.5401
    #   U_ms           = 1.0000e-2 * 0.43089 * 84.5401 = 0.3643 m/s
    u_ms = ums_mathur_gishler(500e-6, 6050.0, 1.66, 0.050, 0.004, 0.10)
    assert u_ms == pytest.approx(0.3643, rel=5e-3)


def test_ums_mathur_gishler_scales_as_sqrt_H():
    # Analytical scaling: U_ms proportional to sqrt(H).
    u1 = ums_mathur_gishler(500e-6, 6050.0, 1.66, 0.050, 0.004, 0.05)
    u2 = ums_mathur_gishler(500e-6, 6050.0, 1.66, 0.050, 0.004, 0.20)
    assert u2 / u1 == pytest.approx(2.0, rel=1e-9)


# Marshall et al. (2017), Table 1: nominal TRISO layers on a 425 um kernel.
MARSHALL_2017_LAYERS = [
    ("buffer", 100e-6, 1000.0),
    ("IPyC", 40e-6, 1900.0),
    ("SiC", 35e-6, 3200.0),
    ("OPyC", 40e-6, 1900.0),
]


def test_coated_particle_stages_marshall_2017_table1():
    # Published (Marshall 2017 Table 1): kernel 425 um, 11.1 g/cm^3 ->
    #   +buffer 625 um 4.18, +IPyC 705 um 3.49, +SiC 775 um 3.42, +OPyC 855 um 3.03
    stages = coated_particle_stages(425e-6, 11100.0, MARSHALL_2017_LAYERS)
    assert [s[0] for s in stages] == ["core", "buffer", "IPyC", "SiC", "OPyC"]
    d_um = [s[1] * 1e6 for s in stages[1:]]
    rho_gcc = [s[2] / 1000.0 for s in stages[1:]]
    assert d_um == pytest.approx([625.0, 705.0, 775.0, 855.0], abs=1e-6)
    assert rho_gcc == pytest.approx([4.18, 3.49, 3.42, 3.03], abs=0.01)


def test_coated_particle_stages_ysz_core_hand_calc():
    # C26 YSZ core (411.6 um, 6000 kg/m^3) with the same layers.
    #   +buffer: d = 611.6 um; f = (411.6/611.6)^3 = 0.30481
    #            rho = 0.30481*6000 + 0.69519*1000 = 2524.0 kg/m^3
    #   +IPyC  : d = 691.6 um; f = (611.6/691.6)^3 = 0.69157
    #            rho = 0.69157*2524.0 + 0.30843*1900 = 2331.6
    #   +SiC   : d = 761.6 um; f = (691.6/761.6)^3 = 0.74883
    #            rho = 0.74883*2331.6 + 0.25117*3200 = 2549.7
    #   +OPyC  : d = 841.6 um; f = (761.6/841.6)^3 = 0.74108
    #            rho = 0.74108*2549.7 + 0.25892*1900 = 2381.5
    stages = coated_particle_stages(411.6e-6, 6000.0, MARSHALL_2017_LAYERS)
    assert [s[1] * 1e6 for s in stages] == pytest.approx(
        [411.6, 611.6, 691.6, 761.6, 841.6], abs=1e-6)
    assert [s[2] for s in stages] == pytest.approx(
        [6000.0, 2524.0, 2331.6, 2549.7, 2381.5], abs=0.5)


def test_coated_particle_stages_conserves_mass():
    # Total mass = core mass + sum of shell masses.
    stages = coated_particle_stages(411.6e-6, 6000.0, MARSHALL_2017_LAYERS)
    vol = lambda d: math.pi / 6.0 * d ** 3  # noqa: E731
    m = 6000.0 * vol(411.6e-6)
    for (_, d0, _), (_, _, layer_rho), (_, d1, rho1) in zip(
            stages[:-1], MARSHALL_2017_LAYERS, stages[1:]):
        m += layer_rho * (vol(d1) - vol(d0))
        assert rho1 * vol(d1) == pytest.approx(m, rel=1e-12)


def test_coated_particle_stages_rejects_bad_input():
    with pytest.raises(ValueError):
        coated_particle_stages(0.0, 6000.0, [])
    with pytest.raises(ValueError):
        coated_particle_stages(411.6e-6, 6000.0, [("buffer", -1e-6, 1000.0)])

"""MFiX 23.2 gas-solids drag closures evaluated at the C26 bed states, relative
to the deck's TGS.  Each function is a line-for-line port of the routine in
~/mamba/envs/mfix-23.2/share/mfix/src/model/drag_gs.f (DRAG_TGS, DRAG_WEN_YU,
DRAG_GIDASPOW, DRAG_BVK, DRAG_KOCH_HILL, DRAG_SYAM_OBRIEN); the returned
quantity is lDgA, the drag coefficient per unit solids volume, so ratios
between closures at one state are the ratios of the interphase momentum
exchange.  Diameters: MFiX 23.2 CGP evaluates drag with the real particle
diameter (DPM = DPA = d_p; notes "CGP drag check").  Gas: the deck's argon.
Written 2026-10-29 for the drag-closure test design
(notes/c26_conry2026_case.md, "Bore interpolation A/B at matched times; drag
closure test").  Not a correlations.py entry (no unit test): a port for a
sensitivity table, not a model input.

    .venv/bin/python notes/c26_drag_closure_table.py
"""
import numpy as np
from scipy.optimize import brentq

MU_G, RHO_G = 2.302324e-05, 1.632918   # deck mu_g0, ro_g0 [Pa s, kg/m3] (argon)
D_P, RHO_S, G = 411.6e-6, 6000.0, 9.81  # YSZ, params/feedstock


def wen_yu(epg, v, d=D_P):
    ropg = RHO_G * epg; re = d * v * ropg / MU_G
    cd = 24 / re * (1 + 0.15 * re ** 0.687) if re <= 1000 else 0.44
    return 0.75 * cd * v * ropg * epg ** -2.65 / d


def gidaspow(epg, v, d=D_P):
    if epg <= 0.8:
        return 150 * (1 - epg) * MU_G / (epg * d * d) + 1.75 * RHO_G * v / d
    return wen_yu(epg, v, d)


def bvk(epg, v, d=D_P):
    phis = 1 - epg; ropg = RHO_G * epg; re = d * v * ropg / MU_G
    a1, a2, a3, a4, a5 = 1 / epg ** 2, 1 / epg, epg, epg ** 2, 1 / re
    fst = 18 * MU_G * a3 / d ** 2
    F = 10 * phis * a1 + a4 * (1 + 1.5 * np.sqrt(phis))
    F += 0.017208 * re * a1 * (a2 + 3 * a3 * phis + 8.4 * a5 ** 0.343) / (1 + 10 ** (3 * phis) * a5 ** (0.5 + 2 * phis))
    return F * fst


def tgs(epg, v, d=D_P):
    phis = 1 - epg; ropg = RHO_G * epg; re = d * v * ropg / MU_G
    fst = 18 * MU_G * epg / d ** 2; fiso = 1 + 0.15 * re ** 0.687
    F = fiso / epg ** 3 + 5.81 * phis / epg ** 3 + 0.48 * phis ** (1 / 3) / epg ** 4
    F += phis ** 3 * re * (0.95 + 0.61 * phis ** 3 / (1 - phis) ** 2)
    return F * fst


def koch_hill(epg, v, d=D_P):
    phis = 1 - epg; ropg = RHO_G * epg; re = 0.5 * d * v * ropg / MU_G
    fst = 18 * MU_G * epg * epg / d ** 2; ww = np.exp(-10 * (0.4 - phis) / phis)
    if 0.01 < phis < 0.4:
        F0 = ((1 - ww) * (1 + 3 * np.sqrt(phis / 2) + 135 / 64 * phis * np.log(phis) + 17.14 * phis)
              / (1 + 0.681 * phis - 8.48 * phis ** 2 + 8.16 * phis ** 3) + ww * 10 * phis / (1 - phis) ** 3)
    else:
        F0 = 10 * phis / (1 - phis) ** 3
    F1 = np.sqrt(2 / phis) / 40 if phis <= 0.1 else 0.11 + 5.1e-4 * np.exp(11.6 * phis)
    if phis < 0.4:
        F2 = ((1 - ww) * (1 + 3 * np.sqrt(phis / 2) + 135 / 64 * phis * np.log(phis) + 17.89 * phis)
              / (1 + 0.681 * phis - 11.03 * phis ** 2 + 15.41 * phis ** 3) + ww * 10 * phis / (1 - phis) ** 3)
    else:
        F2 = 10 * phis / (1 - phis) ** 3
    F3 = 0.9351 * phis + 0.03667 if phis < 0.0953 else 0.0673 + 0.212 * phis + 0.0232 / (1 - phis) ** 5
    re_t2 = (F3 + np.sqrt(F3 * F3 - 4 * F1 * (F0 - F2))) / (2 * F1)
    F = F0 + F1 * re * re if re <= re_t2 else F2 + F3 * re
    return F * fst


def syam_obrien(epg, v, d=D_P, c1=0.8, d1=2.65):
    re = d * v * RHO_G / MU_G; A = epg ** 4.14; B = c1 * epg ** 1.28 if epg <= 0.85 else epg ** d1
    vrm = 0.5 * (A - 0.06 * re + np.sqrt(3.6e-3 * re * re + 0.12 * re * (2 * B - A) + A * A))
    x = re / vrm; cd = (0.63 + 4.8 / np.sqrt(x)) ** 2      # Dalla Valle, C_DSXRE_DV
    return 0.75 * MU_G * epg * cd * x / (vrm * d * d)


CLOSURES = {"TGS (deck)": tgs, "WEN_YU": wen_yu, "GIDASPOW": gidaspow, "BVK": bvk,
            "KOCH_HILL": koch_hill, "SYAM_OBRIEN": syam_obrien}
# (label, EP_g, |slip| [m/s]): annulus EP_g 0.44 from the bore runs' EP_G_mid;
# spout core EP_g 0.93-0.96 and 13-30 m/s from the jet profiles; fountain dilute
STATES = [("annulus packed", 0.44, 0.25), ("annulus, slow", 0.44, 0.10), ("annulus loosened", 0.55, 0.30),
          ("spout core", 0.95, 15.0), ("fountain", 0.99, 2.0)]

if __name__ == "__main__":
    print("lDgA relative to TGS (EP_g, |slip| m/s):")
    print(f"{'closure':14s}" + "".join(f"{s[0] + f' ({s[1]}, {s[2]})':>28s}" for s in STATES))
    for nm, f in CLOSURES.items():
        print(f"{nm:14s}" + "".join(f"{f(e, v) / tgs(e, v):28.2f}" for _, e, v in STATES))
    print("TGS lDgA [kg/m3/s]:", ", ".join(f"{tgs(e, v):.3g}" for _, e, v in STATES))
    print("\nsuperficial gas velocity a packed annulus (EP_g 0.44) passes at the fluidization pressure gradient"
          " (lDgA * slip = (rho_s - rho_g) g) [m/s]:")
    for nm, f in CLOSURES.items():
        vrel = brentq(lambda v: f(0.44, v) * v - (RHO_S - RHO_G) * G, 1e-3, 5.0)
        print(f"  {nm:14s} {0.44 * vrel:.3f}")
    print("  (U_sup at 12 slpm is 0.100 m/s; Wen-Yu U_mf at EP_mf in the sweep summary is 0.238 m/s)")

"""Unit tests for the comparison quantities of notes/c26_sweep_summary.py."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "notes"))
import c26_sweep_summary as S  # noqa: E402


def _frames(tracks):
    """tracks: {id: [y per frame]} -> [(ids, y), ...] with NaN = parcel absent."""
    n = max(len(v) for v in tracks.values())
    out = []
    for k in range(n):
        ids, ys = [], []
        for i, tr in tracks.items():
            if k < len(tr) and not math.isnan(tr[k]):
                ids.append(i)
                ys.append(tr[k])
        out.append((np.array(ids), np.array(ys)))
    return out


def test_plane_transits_counts_up_and_down_once():
    # parcel 1 rises through the plane and comes back down; parcel 2 only rises
    tracks = {1: [-3, -1, 0.0, 1, 3, 1, -1, -3], 2: [-3, -3, -2, 2, 3, 3, 3, 3]}
    n_up, n_down, n_ids_up, n_seen = S.plane_transits(_frames(tracks), 0.0, 1.5)
    assert (n_up, n_down, n_ids_up, n_seen) == (2, 1, 2, 2)


def test_plane_transits_ignores_jitter_in_band():
    # a parcel oscillating inside the +-band never changes state
    tracks = {7: [-3, -1, 1, -1, 1, -1, 1]}
    assert S.plane_transits(_frames(tracks), 0.0, 1.5)[:2] == (0, 0)
    # ... but a parcel that fully crosses the band does
    tracks = {7: [-3, -1, 1, 2, 1, -1, -2]}
    assert S.plane_transits(_frames(tracks), 0.0, 1.5)[:2] == (1, 1)


def test_plane_transits_starting_in_band_and_deleted_parcels():
    # parcel 3 is first seen inside the band and then goes up: no transit
    # parcel 4 goes up and is then deleted (23.2 cut-cell loss): one transit
    tracks = {3: [0.5, 1, 3, 3], 4: [-3, 3, math.nan, math.nan]}
    n_up, n_down, n_ids_up, n_seen = S.plane_transits(_frames(tracks), 0.0, 1.5)
    assert (n_up, n_down, n_ids_up, n_seen) == (1, 0, 1, 2)


def test_spouting_margin_c02_hand_worked():
    # C26 c02 manifest values: YSZ 411.6 um, 6000 kg/m3, Ar 1.633 kg/m3,
    # 2.302e-5 Pa s, D_c 50.4 mm, D_i 2.381 mm, H 34.84 mm, U_sup 0.1002 m/s
    m = {'d_p_m': 411.6e-6, 'rho_p_kg_m3': 6000.0, 'rho_g_kg_m3': 1.6329,
         'mu_g_pa_s': 2.3023e-5, 'D_c_m': 0.0504, 'R_t_m': 0.001190625,
         'H_static_m': 0.034844, 'U_sup_m_s': 0.100249}
    r = S.spouting_margin(m)
    # U_ms = (d/Dc) (Di/Dc)^(1/3) sqrt(2 g H (rho_p - rho)/rho)
    u_ms = (411.6e-6 / 0.0504) * (0.00238125 / 0.0504) ** (1 / 3) * math.sqrt(
        2 * 9.80665 * 0.034844 * (6000 - 1.6329) / 1.6329)
    assert r['U_ms_mg_m_s'] == pytest.approx(u_ms, rel=1e-6)
    assert r['U_ms_mg_m_s'] == pytest.approx(0.1479, abs=5e-4)     # hand value: 0.008167 * 0.3614 * 50.11
    assert r['U_over_U_ms'] == pytest.approx(0.100249 / u_ms, rel=1e-6)
    assert r['Ar'] == pytest.approx(1.6329 * (6000 - 1.6329) * 9.80665 * 411.6e-6 ** 3 / 2.3023e-5 ** 2, rel=1e-6)
    # a spouted bed runs below U_mf on the column cross-section (annulus not fluidized);
    # Wen-Yu for 412 um YSZ in Ar: Ar = 1.26e4 -> Re_mf = 6.94 -> U_mf = 0.238 m/s
    assert r['U_mf_wenyu_m_s'] == pytest.approx(0.238, abs=2e-3)
    assert r['U_over_U_mf'] < 1.0


def test_spouting_margin_larger_lighter_particle_lowers_the_margin():
    base = {'d_p_m': 411.6e-6, 'rho_p_kg_m3': 6000.0, 'rho_g_kg_m3': 1.6329,
            'mu_g_pa_s': 2.3023e-5, 'D_c_m': 0.0504, 'R_t_m': 0.001190625,
            'H_static_m': 0.034844, 'U_sup_m_s': 0.100249}
    opyc = dict(base, d_p_m=841.6e-6, rho_p_kg_m3=2381.5, H_static_m=0.049193)
    # Mathur-Gishler: U_ms grows with d_p and H, falls with rho_p; net rises for OPyC
    assert S.spouting_margin(opyc)['U_over_U_ms'] < S.spouting_margin(base)['U_over_U_ms']

"""Unit tests for src.post.piv_c26 (synthetic PIV and C26 map metrics) and
src.post.digitize_c26_fig4 (colormap inversion, digitized Fig. 4 files)."""

from __future__ import annotations

import numpy as np
import pytest

from src.post import piv_c26 as P
from src.post.digitize_c26_fig4 import fill_nan, invert_colormap

PX = 1.0e-4                       # pixel size [m]
WINDOWS = [64, 32, 32]            # params/piv_c26.yaml piv.windows_px
KW = dict(windows_px=WINDOWS, min_std=1e-3, uod_threshold=2.0, uod_eps=0.1)


def _shifted_pair(shift_px, n=600, seed=0, size=200):
    """Two renders of random spheres (0.6 mm, 6 px) moved by shift_px = (dx, dy)."""
    rng = np.random.default_rng(seed)
    xyz = np.column_stack([rng.uniform(0, size * PX, n), rng.uniform(0, size * PX, n),
                           rng.uniform(0, 0.01, n)])
    d = np.full(n, 0.6e-3)
    a, _ = P.render(xyz, d, (0.0, 0.0), (size, size), PX)
    b, _ = P.render(xyz + np.array([shift_px[0] * PX, shift_px[1] * PX, 0.0]), d, (0.0, 0.0),
                    (size, size), PX)
    return a, b


def _grid(size=200, margin=64, step=16):
    c = np.arange(margin, size - margin + 1, step)
    return np.meshgrid(c, c, indexing="ij")


def test_render_single_sphere_profile():
    """R = 5 px sphere centred on pixel (10, 10): intensity sqrt(1 - (r/R)^2)."""
    xyz = np.array([[10.5 * PX, 10.5 * PX, 0.0]])
    img, own = P.render(xyz, np.array([10 * PX]), (0.0, 0.0), (21, 21), PX)
    assert img[10, 10] == pytest.approx(1.0)
    assert img[10, 14] == pytest.approx(np.sqrt(1 - 16 / 25))    # 0.6
    assert img[14, 10] == pytest.approx(0.6)                      # rows follow y
    assert img[10, 15] == pytest.approx(0.0, abs=1e-6)            # r = R: rim
    assert img[10, 16] == 0.0 and own[10, 16] == -1               # outside
    assert own[10, 10] == 0


def test_render_nearest_sphere_wins():
    """Camera at z -> -inf: the sphere with the smaller z hides the other."""
    xyz = np.array([[1e-3, 1e-3, 0.0], [1e-3, 1e-3, -1e-3]])
    _, own = P.render(xyz, np.array([0.8e-3, 0.8e-3]), (0.0, 0.0), (20, 20), PX)
    assert own[9, 9] == 1


def test_correlate_integer_shift_sign():
    """B = A moved by (+3 rows, -2 cols) gives displacement (+3, -2)."""
    rng = np.random.default_rng(1)
    A = rng.random((1, 32, 32))
    B = np.roll(A, (3, -2), axis=(1, 2))
    dr, dc = P.correlate(A, B)
    assert dr[0] == pytest.approx(3.0, abs=0.05)
    assert dc[0] == pytest.approx(-2.0, abs=0.05)


@pytest.mark.parametrize("shift", [(2.3, -1.7), (-0.4, 0.6), (0.05, -0.02)])
def test_piv_recovers_subpixel_shift(shift):
    """Mean bias < 0.03 px, including wall-like shifts of a few hundredths of
    a pixel (C26 wall velocity 1e-3 m/s = 0.02 px per frame at 0.1 mm/px)."""
    a, b = _shifted_pair(shift)
    rows, cols = _grid()
    dr, dc, empty = P.piv_pair(a, b, rows, cols, **KW)
    assert not empty.any()
    assert dc.mean() == pytest.approx(shift[0], abs=0.03)
    assert dr.mean() == pytest.approx(shift[1], abs=0.03)


def test_piv_multipass_recovers_large_shift():
    """12 px is too large for a lone 32 px window; the 64 px pass offsets it."""
    a, b = _shifted_pair((0.0, 12.0))
    rows, cols = _grid()
    dr, dc, _ = P.piv_pair(a, b, rows, cols, **KW)
    assert np.median(dr) == pytest.approx(12.0, abs=0.1)
    assert np.median(dc) == pytest.approx(0.0, abs=0.1)


def test_piv_static_and_empty_images_give_zero():
    a, _ = _shifted_pair((0.0, 0.0))
    rows, cols = _grid()
    dr, dc, _ = P.piv_pair(a, a, rows, cols, **KW)
    assert np.abs(dr).max() < 0.01 and np.abs(dc).max() < 0.01
    z = np.zeros_like(a)
    dr, dc, empty = P.piv_pair(z, z, rows, cols, **KW)
    assert empty.all() and not dr.any() and not dc.any()


def test_median_test_flags_isolated_spike():
    dr = np.zeros((5, 5))
    dr[2, 2] = 5.0
    out, mr, _ = P.median_test(dr, np.zeros((5, 5)), np.ones((5, 5), bool), 2.0, 0.1)
    assert out[2, 2] and out.sum() == 1
    assert mr[2, 2] == 0.0


def test_map_metrics_hand_values():
    a = np.array([[0.0, 1.0], [2.0, 3.0]])
    assert P.pearson_cc(a, 2 * a + 1) == pytest.approx(1.0)
    assert P.pearson_cc(a, -a) == pytest.approx(-1.0)
    assert P.mse(np.zeros((2, 2)), np.full((2, 2), 0.1)) == pytest.approx(0.01)
    rng = np.random.default_rng(2)
    m = rng.random((20, 20))
    assert P.ssim(m, m) == pytest.approx(1.0)
    # constant maps: sigma terms vanish, SSIM = (2ab + C1) / (a^2 + b^2 + C1)
    # with C1 = (0.01 * 1)^2: (0.04 + 1e-4) / (0.05 + 1e-4) = 0.800399...
    s = P.ssim(np.full((10, 10), 0.1), np.full((10, 10), 0.2))
    assert s == pytest.approx(0.0401 / 0.0501, rel=1e-9)


def test_make_view_grid_inside_box():
    v = P.make_view(-0.01, 0.01, 0.0, 0.008, PX, 64, 16)
    assert v.shape == (80 + 128, 200 + 128)
    assert v.x.min() > -0.01 and v.x.max() < 0.01
    assert v.y.min() > 0.0 and v.y.max() < 0.008
    assert np.allclose(np.diff(v.x), 16 * PX)


def _matlab_jet(t):
    """MATLAB jet, RGB 0-255: piecewise linear, breaks at 1/8, 3/8, 5/8, 7/8.
    (Matplotlib's jet has a flat blue stretch at t = 0.11-0.125, so it cannot
    be inverted there; C26 Fig. 4 uses MATLAB's.)"""
    t = np.asarray(t)[..., None]
    lo = np.array([-1.5, -0.5, 0.5])          # r, g, b rising edges: 4t + lo
    hi = np.array([4.5, 3.5, 2.5])            # falling edges: -4t + hi
    return np.clip(np.minimum(4 * t + lo, -4 * t + hi), 0, 1) * 255


def test_invert_colormap_jet_round_trip():
    lut = _matlab_jet(np.linspace(0, 1, 1000))
    values = np.linspace(4e-4, 0.4, 1000)
    v = np.linspace(4e-4, 0.4, 50).reshape(5, 10)
    rgb = _matlab_jet((v - 4e-4) / (0.4 - 4e-4))
    rgb[0, 0] = (128, 128, 128)                    # grey: off the colormap
    out = invert_colormap(rgb, lut, values, max_dist=30)
    assert np.isnan(out[0, 0])
    ok = ~np.isnan(out)
    assert np.abs(out[ok] - v[ok]).max() < 1.5 * (values[1] - values[0])
    assert not np.isnan(fill_nan(out)).any()


@pytest.mark.parametrize("flow", [10, 12, 14])
def test_digitized_fig4_maps(flow):
    """Committed digitized maps: experiment and model on one grid, in the
    colorbar range."""
    exp, model = P.load_fig4("exp", flow), P.load_fig4("model", flow)
    assert exp.shape == model.shape and exp.shape[1] == 100
    for a in (exp, model):
        assert a.min() >= 4e-4 - 1e-4 and a.max() <= 0.4 + 1e-4

"""Unit tests for src.post.psd (Welch dominant-frequency estimator)."""

from __future__ import annotations

import numpy as np
import pytest

from src.post.psd import dominant_frequency, welch_psd


def test_psd_recovers_pure_sine():
    """Welch PSD of a sine at f0 recovers f0 within the frequency resolution.

    8 s of a 20 Hz sine sampled at 1 kHz -> df = fs/nperseg; recovered peak
    should be within one bin of 20.0 Hz.
    """
    fs = 1000.0
    t = np.arange(0.0, 8.0, 1.0 / fs)
    f0 = 20.0
    sig = np.sin(2 * np.pi * f0 * t)
    res = welch_psd(sig, fs_hz=fs)
    df = res.freqs_hz[1] - res.freqs_hz[0]
    assert abs(res.dominant_hz - f0) <= df


def test_psd_rejects_dc():
    """A DC-dominated signal must still return a non-DC peak at f0."""
    fs = 1000.0
    t = np.arange(0.0, 4.0, 1.0 / fs)
    f0 = 50.0
    sig = 100.0 + 1.0 * np.sin(2 * np.pi * f0 * t)
    f = dominant_frequency(sig, fs_hz=fs)
    assert 40.0 < f < 60.0


def test_psd_picks_louder_tone_when_two_present():
    """Two tones: f1 at amplitude 1, f2 at amplitude 0.1 — pick f1."""
    fs = 1000.0
    t = np.arange(0.0, 8.0, 1.0 / fs)
    sig = np.sin(2 * np.pi * 25.0 * t) + 0.1 * np.sin(2 * np.pi * 100.0 * t)
    f = dominant_frequency(sig, fs_hz=fs)
    assert 20.0 < f < 30.0


def test_psd_rejects_bad_input():
    with pytest.raises(ValueError):
        welch_psd(np.zeros(10), fs_hz=1000.0)   # too short
    with pytest.raises(ValueError):
        welch_psd(np.zeros(100), fs_hz=-1.0)    # bad fs
    with pytest.raises(ValueError):
        welch_psd(np.zeros((10, 10)), fs_hz=1000.0)  # not 1-D

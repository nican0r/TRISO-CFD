"""Welch power-spectral-density helper for inlet-pressure time series.

Returns the dominant frequency peak of the PSD (strictly > 0 Hz, so the
DC bin is ignored) using scipy.signal.welch. All quantities SI.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple

import numpy as np
from scipy.signal import welch


@dataclass
class PSDResult:
    freqs_hz: np.ndarray
    psd: np.ndarray
    dominant_hz: float
    fs_hz: float


def welch_psd(
    signal: np.ndarray,
    fs_hz: float,
    nperseg: int | None = None,
    window: str = "hann",
) -> PSDResult:
    """Welch PSD estimator; returns (freqs, psd, dominant_frequency, fs).

    signal     : 1-D array of samples (any physical unit) [-]
    fs_hz      : sampling rate [Hz]
    nperseg    : Welch segment length (default: min(len/4, 1024))
    window     : scipy.signal window name (default 'hann')

    The dominant frequency is argmax(PSD) restricted to f > 0 (ignore DC).
    """
    sig = np.asarray(signal, dtype=float)
    if sig.ndim != 1:
        raise ValueError(f"signal must be 1-D, got shape {sig.shape}")
    if fs_hz <= 0.0:
        raise ValueError(f"fs_hz must be > 0, got {fs_hz}")
    if sig.size < 16:
        raise ValueError(f"signal must have >=16 samples, got {sig.size}")
    if nperseg is None:
        nperseg = min(max(sig.size // 4, 64), 1024)
    nperseg = min(nperseg, sig.size)

    freqs, psd = welch(sig, fs=fs_hz, nperseg=nperseg, window=window)
    mask = freqs > 0.0
    if not np.any(mask):
        raise RuntimeError("no positive-frequency bins in PSD")
    k = int(np.argmax(psd[mask]))
    dom = float(freqs[mask][k])
    return PSDResult(freqs_hz=freqs, psd=psd, dominant_hz=dom, fs_hz=fs_hz)


def dominant_frequency(signal: np.ndarray, fs_hz: float) -> float:
    """Convenience: just the dominant frequency [Hz]."""
    return welch_psd(signal, fs_hz).dominant_hz

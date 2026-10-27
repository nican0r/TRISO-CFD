"""Synthetic PIV for Case C26 (Conry et al. 2026, Table 3 and Fig. 6d).

C26 scores its model by running PIVlab on side-view images of the simulated
bed and comparing them with PIV of the experiment (p.5 3.1, p.9 3.2.1).
This module does the same on MFiX PARTICLES_PIV frames:

1. render each frame as a side-view camera image: orthographic along z,
   opaque parcels at their VTK diameter, nearest surface wins;
2. run windowed FFT cross-correlation on consecutive frames: multi-pass with
   sub-pixel window offset, 3-point Gaussian sub-pixel peak, normalized median
   test (Westerweel & Scarano 2005) after each pass;
3. spout map: |time mean of the PIV vectors| over the region above the bed
   (C26 Fig. 4), resampled onto the digitized Fig. 4 grid and scored with
   CC, SSI and MSE;
4. wall velocity: mean vertical PIV velocity in a 4 mm box at mid-bed on the
   bed wall (C26 Fig. 1c), downward positive (C26 Fig. 6).

The same metrics are computed between the digitized experiment and C26's own
model map (Fig. 4a vs 4b), so a run is compared with C26's model through one
pipeline as well as with Table 3.  Settings and sources: params/piv_c26.yaml;
design notes: notes/c26_conry2026_case.md ("Synthetic PIV").

    MPLBACKEND=Agg .venv/bin/python -m src.post.piv_c26 cases/C26_3d/c01/aws_out

All values SI (CLAUDE.md rule 3); image quantities are in pixels [px].
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import yaml
from scipy import ndimage, sparse
from scipy.sparse.linalg import spsolve
from scipy.interpolate import RegularGridInterpolator

REPO = Path(__file__).resolve().parents[2]
PARAMS = REPO / "params" / "piv_c26.yaml"
DIGITIZED = REPO / "data" / "digitized"


def load_params(path: Path = PARAMS) -> dict:
    with Path(path).open() as f:
        return yaml.safe_load(f)


# ----------------------------------------------------------------------------
# Rendering
# ----------------------------------------------------------------------------
def render(xyz, diam, origin, shape, pixel_m):
    """Side-view image of spheres seen from z -> -inf.

    xyz      : (N, 3) centres [m]; image columns follow x, rows follow y (row 0
               at the bottom, so rows increase upward)
    diam     : (N,) diameters [m]
    origin   : (x0, y0) of the image's lower-left corner [m]
    shape    : (ny, nx) [px]
    pixel_m  : pixel size [m]

    Returns (img, owner): intensity sqrt(1 - (r/R)^2) of the nearest sphere
    surface per pixel (0 = background) and that sphere's index (-1 = none).
    """
    ny, nx = shape
    img = np.zeros(shape)
    owner = np.full(shape, -1)
    u = (xyz[:, 0] - origin[0]) / pixel_m - 0.5    # pixel j centre is u = j
    v = (xyz[:, 1] - origin[1]) / pixel_m - 0.5
    R = 0.5 * np.asarray(diam) / pixel_m
    idx = np.flatnonzero((u > -R - 1) & (u < nx + R) & (v > -R - 1) & (v < ny + R))
    if idx.size == 0:
        return img, owner
    u, v, R, z = u[idx], v[idx], R[idx], xyz[idx, 2]
    m = int(np.ceil(R.max()))
    di, dj = (a.ravel() for a in np.meshgrid(np.arange(-m, m + 1), np.arange(-m, m + 1), indexing="ij"))
    I = np.rint(v).astype(int)[:, None] + di
    J = np.rint(u).astype(int)[:, None] + dj
    r2 = ((I - v[:, None]) ** 2 + (J - u[:, None]) ** 2) / R[:, None] ** 2
    ok = (r2 < 1.0) & (I >= 0) & (I < ny) & (J >= 0) & (J < nx)
    k = np.broadcast_to(np.arange(idx.size)[:, None], I.shape)[ok]
    h = np.sqrt(1.0 - r2[ok])
    depth = z[k] - R[k] * pixel_m * h          # sphere surface facing the camera
    flat = I[ok] * nx + J[ok]
    order = np.lexsort((depth, flat))          # per pixel, nearest first
    fs = flat[order]
    sel = order[np.r_[True, fs[1:] != fs[:-1]]]
    img.flat[flat[sel]] = h[sel]
    owner.flat[flat[sel]] = idx[k[sel]]
    return img, owner


# ----------------------------------------------------------------------------
# PIV
# ----------------------------------------------------------------------------
def _windows(img, r0, c0, w):
    """(K, w, w) windows with top-left corners at fractional (r0, c0) [px],
    bilinear interpolation, zero outside the image."""
    k = np.arange(w, dtype=float)
    r = np.broadcast_to(r0[:, None, None] + k[None, :, None], (r0.size, w, w))
    c = np.broadcast_to(c0[:, None, None] + k[None, None, :], (r0.size, w, w))
    return ndimage.map_coordinates(img, [r.ravel(), c.ravel()], order=1,
                                   mode="constant").reshape(r0.size, w, w)


def _subpixel(cm, c0, cp):
    """Peak offset from three samples at -1, 0, +1: Gaussian fit, parabolic
    where a sample is not positive."""
    with np.errstate(divide="ignore", invalid="ignore"):
        pos = (cm > 0) & (c0 > 0) & (cp > 0)
        lm, l0, lp = (np.log(np.where(pos, a, 1.0)) for a in (cm, c0, cp))
        g = (lm - lp) / (2.0 * (lm - 2.0 * l0 + lp))
        p = (cm - cp) / (2.0 * (cm - 2.0 * c0 + cp))
    s = np.where(pos & np.isfinite(g), g, p)
    return np.clip(np.nan_to_num(s), -0.5, 0.5)


def correlate(A, B):
    """Displacement (dr, dc) [px] of each window in B relative to A.

    A, B: (K, w, w) windows.  Circular FFT cross-correlation of the
    mean-subtracted windows; the peak C(s) = sum_x A(x) B(x + s) sits at the
    displacement.
    """
    K, w, _ = A.shape
    A = A - A.mean(axis=(1, 2), keepdims=True)
    B = B - B.mean(axis=(1, 2), keepdims=True)
    C = np.fft.irfft2(np.conj(np.fft.rfft2(A)) * np.fft.rfft2(B), s=(w, w))
    C = np.fft.fftshift(C, axes=(1, 2))
    pr, pc = np.divmod(C.reshape(K, -1).argmax(axis=1), w)
    n = np.arange(K)
    out = []
    for p, q, axis in ((pr, pc, 0), (pc, pr, 1)):
        inner = (p > 0) & (p < w - 1)
        pm, pp = np.clip(p - 1, 0, w - 1), np.clip(p + 1, 0, w - 1)
        if axis == 0:
            s = _subpixel(C[n, pm, q], C[n, p, q], C[n, pp, q])
        else:
            s = _subpixel(C[n, q, pm], C[n, q, p], C[n, q, pp])
        out.append(p - w // 2 + np.where(inner, s, 0.0))
    return out[0], out[1]


def _neighbours(a):
    p = np.pad(a, 1, constant_values=np.nan)
    ny, nx = a.shape
    return np.stack([p[1 + i:1 + i + ny, 1 + j:1 + j + nx]
                     for i in (-1, 0, 1) for j in (-1, 0, 1) if i or j])


def median_test(dr, dc, valid, threshold, eps):
    """Normalized median test on the 3x3 neighbourhood (Westerweel & Scarano
    2005).  Returns (outlier mask, neighbour median of dr, of dc); invalid
    points are left out of every neighbourhood."""
    res2 = np.zeros(dr.shape)
    meds = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)   # all-NaN neighbourhoods
        for d in (dr, dc):
            nb = _neighbours(np.where(valid, d, np.nan))
            med = np.nanmedian(nb, axis=0)
            rm = np.nanmedian(np.abs(nb - med), axis=0)
            res2 += ((d - med) / (rm + eps)) ** 2
            meds.append(med)
    out = valid & (np.sqrt(res2) > threshold)            # NaN (no neighbours) -> kept
    return out, meds[0], meds[1]


def inpaint_nans(a):
    """Fill the NaNs of a 2-D field by solving the discrete Laplace equation
    with the valid values as boundary conditions (PIVlab's "interpolate
    missing data", inpaint_nans).  Harmonic fields (e.g. a linear ramp) are
    recovered exactly; an all-NaN field is returned unchanged."""
    a = np.asarray(a, float)
    bad = np.isnan(a)
    if not bad.any() or bad.all():
        return a.copy()
    ny, nx = a.shape
    pos = np.full(a.shape, -1, int)
    pos[bad] = np.arange(bad.sum())
    rows, cols, vals = [], [], []
    rhs = np.zeros(bad.sum())
    for k, (i, j) in enumerate(zip(*np.nonzero(bad))):
        nb = [(i + di, j + dj) for di, dj in ((-1, 0), (1, 0), (0, -1), (0, 1))
              if 0 <= i + di < ny and 0 <= j + dj < nx]
        rows.append(k); cols.append(k); vals.append(float(len(nb)))
        for ii, jj in nb:
            if bad[ii, jj]:
                rows.append(k); cols.append(pos[ii, jj]); vals.append(-1.0)
            else:
                rhs[k] += a[ii, jj]
    lap = sparse.csc_matrix((vals, (rows, cols)), shape=(rhs.size, rhs.size))
    out = a.copy()
    out[bad] = spsolve(lap, rhs)
    return out


EMPTY_WINDOW_MODES = ("zero", "skip", "interpolate")


def piv_pair(a, b, rows, cols, windows_px, min_std, uod_threshold, uod_eps):
    """Displacement field [px] from image a to image b.

    rows, cols : integer arrays (same 2-D shape) of window centres [px]
    windows_px : window size per pass, largest first.  Each pass moves the A
                 window back and the B window forward by half the previous
                 estimate (sub-pixel, bilinear) and correlates the residual,
                 so repeated last passes converge where the peak fit is
                 unbiased (PIVlab's iterative window offset, without warping)
    Returns (dr, dc, empty): dr along rows (up), dc along columns; empty marks
    windows with intensity std < min_std in either frame (displacement 0).
    """
    shape = rows.shape
    r, c = rows.ravel().astype(float), cols.ravel().astype(float)
    wmax = max(windows_px)
    dr = np.zeros(r.size)
    dc = np.zeros(r.size)
    empty = np.zeros(r.size, bool)
    for w in windows_px:
        pr, pc = np.clip(dr, -wmax, wmax), np.clip(dc, -wmax, wmax)
        WA = _windows(a, r - pr / 2 - w // 2, c - pc / 2 - w // 2, w)
        WB = _windows(b, r + pr / 2 - w // 2, c + pc / 2 - w // 2, w)
        empty = (WA.std(axis=(1, 2)) < min_std) | (WB.std(axis=(1, 2)) < min_std)
        sr, sc = correlate(WA, WB)
        dr = np.where(empty, 0.0, pr + sr)
        dc = np.where(empty, 0.0, pc + sc)
        out, mr, mc = median_test(dr.reshape(shape), dc.reshape(shape),
                                  ~empty.reshape(shape), uod_threshold, uod_eps)
        out = out.ravel()
        dr = np.where(out, np.nan_to_num(mr.ravel()), dr)
        dc = np.where(out, np.nan_to_num(mc.ravel()), dc)
    return dr.reshape(shape), dc.reshape(shape), empty.reshape(shape)


# ----------------------------------------------------------------------------
# Views
# ----------------------------------------------------------------------------
@dataclass
class View:
    """A camera image region and its PIV grid.

    The image covers the target box plus a margin of one pass-1 window on each
    side, so every window near the box edge sees particles.  Grid points
    (rows, cols) are window centres inside the box at the PIV step.
    """
    x0: float
    y0: float
    shape: tuple
    pixel_m: float
    rows: np.ndarray
    cols: np.ndarray

    @property
    def x(self):
        return self.x0 + (self.cols[0] + 0.5) * self.pixel_m

    @property
    def y(self):
        return self.y0 + (self.rows[:, 0] + 0.5) * self.pixel_m


def bed_surface(xyz, percentile, r_min):
    """Bed-surface height [m] from one parcel frame: the given percentile of
    parcel y over the annulus parcels, those at least r_min from the axis
    (x = z = 0).  The fountain parcels standing above the spout are left
    out, as C26 masks the annular bed surface (Fig. 4 caption).  Falls back
    to all parcels when none is outside r_min."""
    xyz = np.asarray(xyz, float)
    sel = np.hypot(xyz[:, 0], xyz[:, 2]) >= r_min
    if not sel.any():
        sel = np.ones(len(xyz), bool)
    return float(np.percentile(xyz[sel, 1], percentile))


def make_view(x_lo, x_hi, y_lo, y_hi, pixel_m, margin_px, step_px):
    nx_box = int(round((x_hi - x_lo) / pixel_m))
    ny_box = int(round((y_hi - y_lo) / pixel_m))

    def centres(n_box):
        n = max(1, n_box // step_px)
        start = margin_px + (n_box - (n - 1) * step_px) // 2
        return start + step_px * np.arange(n)

    rows, cols = np.meshgrid(centres(ny_box), centres(nx_box), indexing="ij")
    shape = (ny_box + 2 * margin_px, nx_box + 2 * margin_px)
    return View(x_lo - margin_px * pixel_m, y_lo - margin_px * pixel_m, shape, pixel_m, rows, cols)


# ----------------------------------------------------------------------------
# Map metrics (C26 p.5: CC, SSI, MSE)
# ----------------------------------------------------------------------------
def pearson_cc(a, b):
    return float(np.corrcoef(np.ravel(a), np.ravel(b))[0, 1])


def mse(a, b):
    return float(np.mean((np.asarray(a) - np.asarray(b)) ** 2))


def ssim(a, b, sigma=1.5, radius=5, k1=0.01, k2=0.03, dynamic_range=1.0):
    """Mean structural similarity (Wang et al. 2004) with a Gaussian window,
    replicate padding, as MATLAB ssim."""
    c1, c2 = (k1 * dynamic_range) ** 2, (k2 * dynamic_range) ** 2

    def f(m):
        return ndimage.gaussian_filter(m, sigma, mode="nearest", truncate=radius / sigma)

    a, b = np.asarray(a, float), np.asarray(b, float)
    ma, mb = f(a), f(b)
    saa, sbb, sab = f(a * a) - ma ** 2, f(b * b) - mb ** 2, f(a * b) - ma * mb
    s = (2 * ma * mb + c1) * (2 * sab + c2) / ((ma ** 2 + mb ** 2 + c1) * (saa + sbb + c2))
    return float(s.mean())


def map_metrics(ref, test, mp):
    return {"CC": pearson_cc(ref, test),
            "SSI": ssim(ref, test, mp["ssim_sigma_px"], mp["ssim_radius_px"], mp["ssim_k1"],
                        mp["ssim_k2"], mp["ssim_dynamic_range"]),
            "MSE_m2_s2": mse(ref, test)}


def fig4_path(kind, flow):
    return DIGITIZED / f"c26_fig4_{kind}_{flow:g}slpm.csv"


def load_fig4(kind, flow):
    """Digitized C26 Fig. 4 map [m/s] ('exp' or 'model'), row 0 at the bottom
    (the CSV stores the figure's top row first)."""
    return np.flipud(np.loadtxt(fig4_path(kind, flow), delimiter=",", comments="#"))


# ----------------------------------------------------------------------------
# Run analysis
# ----------------------------------------------------------------------------
def _vtk_readers():
    sys.path.insert(0, str(REPO / ".claude/skills/mfix-smoke-validate"))
    from check_smoke import read_pvd, read_vtk_xml  # noqa: E402
    return read_pvd, read_vtk_xml


def piv_frames(run, manifest):
    read_pvd, _ = _vtk_readers()
    t0, t1 = float(manifest["piv_t_start_s"]), float(manifest["tstop_s"]) + 1e-4
    frames = []
    for t, f in read_pvd(run / "PARTICLES_PIV.pvd"):
        # MFiX 23.2 can write the tstop frame twice (same time, identical file,
        # c02_point_v232); a pair with dt = 0 would put inf into the map.
        if t0 <= t <= t1 and not (frames and t <= frames[-1][0]):
            frames.append((t, run / f))
    return frames


def run_geometry(run, params=None, stride=1, max_frames=None):
    """Camera views and sampling regions of one run dir over the manifest's
    PIV window (shared by analyse_run and the frame exporter for PIVlab).

    Returns a dict: manifest m, frames [(t, path)], flow, R_c, px, step, wins,
    cols_xyz (parcel columns seen by the camera), y_s (bed surface), y_map
    (map bottom), H (map height), pflow / nrow / ncol (digitized experiment
    panel and its grid), the spout and wall Views, xc / yc / half (wall box),
    in_roi (wall grid points inside the box) and roi_pix (box pixels)."""
    run = Path(run)
    p = params or load_params()
    m = json.loads((run / "manifest.json").read_text())
    _, read_vtk_xml = _vtk_readers()
    frames = piv_frames(run, m)[::stride][:max_frames]
    if len(frames) < 2:
        raise ValueError(f"{run}: fewer than 2 PIV frames in the window")

    cam, piv, sm, wr = p["camera"], p["piv"], p["spout_map"], p["wall_roi"]
    px, step, wins = float(cam["pixel_m"]), int(piv["step_px"]), [int(w) for w in piv["windows_px"]]
    margin = max(wins)
    R_c = float(m["R_c_m"])
    flow = float(m["flow_slpm"])

    # line of sight: camera along z sees (x, y); along x it sees (z, y)
    cols_xyz = {"z": [0, 1, 2], "x": [2, 1, 0]}[cam["view_axis"]]
    # bed surface: time-median over the frames of the percentile of parcel y
    # over the annulus parcels (r >= bed_surface_annulus_r_min_m), outside the
    # spout / fountain footprint; C26 masks "the bed annular surface".
    r_min = float(sm["bed_surface_annulus_r_min_m"])
    y_s = float(np.median([bed_surface(read_vtk_xml(f)["coordinates"], sm["bed_surface_percentile"], r_min)
                           for _, f in frames]))
    y_map = y_s + float(sm["bottom_offset_m"])
    # comparison grid: the digitized experiment panel of the nearest flow
    panels = sorted(p["fig4"]["panels_exp"])
    pflow = min(panels, key=lambda q: abs(q - flow))
    nrow, ncol = np.loadtxt(fig4_path("exp", pflow), delimiter=",", comments="#").shape
    H = 2 * R_c * nrow / ncol
    spout = make_view(-R_c, R_c, y_map, y_map + H, px, margin, step)
    half = 0.5 * float(wr["size_m"])
    yc = float(wr["y_center_frac"]) * y_s
    xc = float(wr["x_center_m"])
    wall = make_view(xc - half - step * px, xc + half + step * px,
                     yc - half - step * px, yc + half + step * px, px, margin, step)
    in_roi = (np.abs(wall.x[None, :] - xc) <= half) & (np.abs(wall.y[:, None] - yc) <= half)
    r_lo, r_hi = (np.array([yc - half, yc + half]) - wall.y0) / px
    c_lo, c_hi = (np.array([xc - half, xc + half]) - wall.x0) / px
    roi_pix = (slice(int(r_lo), int(np.ceil(r_hi))), slice(int(c_lo), int(np.ceil(c_hi))))
    return dict(m=m, frames=frames, flow=flow, R_c=R_c, px=px, step=step, wins=wins, cols_xyz=cols_xyz,
                y_s=y_s, y_map=y_map, H=H, pflow=pflow, nrow=nrow, ncol=ncol, spout=spout, wall=wall,
                xc=xc, yc=yc, half=half, in_roi=in_roi, roi_pix=roi_pix)


def fig4_grid_map(mean_map, spout, y_map, H, R_c, nrow, ncol):
    """Resample a map on the spout View's grid onto the digitized Fig. 4 grid
    (nrow x ncol over x in [-R_c, R_c], height H above y_map)."""
    gy, gx = spout.y - y_map, spout.x
    ci = (np.arange(nrow) + 0.5) * H / nrow
    cj = -R_c + (np.arange(ncol) + 0.5) * 2 * R_c / ncol
    Y, X = np.meshgrid(np.clip(ci, gy[0], gy[-1]), np.clip(cj, gx[0], gx[-1]), indexing="ij")
    return RegularGridInterpolator((gy, gx), mean_map)(np.stack([Y, X], -1))


def analyse_run(run, params=None, stride=1, max_frames=None, on_pair=None):
    """Synthetic PIV of one run dir over the manifest's PIV window.
    on_pair: optional callback(info) after each frame pair, with info =
    {"i", "t0", "t1", "spout": (img_a, img_b, dr, dc, empty), "wall": (...)};
    dr / dc are the fields as accumulated (spout: after the empty-window
    policy) [px per pair].  Used by the PIVlab frame exporter."""
    run = Path(run)
    p = params or load_params()
    _, read_vtk_xml = _vtk_readers()
    g = run_geometry(run, p, stride, max_frames)
    m, frames, flow, R_c, px, step, wins = (g[k] for k in ("m", "frames", "flow", "R_c", "px", "step", "wins"))
    cols_xyz, y_s, y_map, H, pflow, nrow, ncol = (g[k] for k in ("cols_xyz", "y_s", "y_map", "H", "pflow",
                                                                  "nrow", "ncol"))
    spout, wall, xc, yc, in_roi, roi_pix = (g[k] for k in ("spout", "wall", "xc", "yc", "in_roi", "roi_pix"))
    cam, piv, sm = p["camera"], p["piv"], p["spout_map"]

    kw = dict(windows_px=wins, min_std=float(piv["min_window_std"]),
              uod_threshold=float(piv["uod_threshold"]), uod_eps=float(piv["uod_eps_px"]))
    mode = str(sm["empty_windows"])
    if mode not in EMPTY_WINDOW_MODES:
        raise ValueError(f"spout_map.empty_windows must be one of {EMPTY_WINDOW_MODES}, got {mode!r}")
    su = np.zeros(spout.rows.shape)
    sv = np.zeros(spout.rows.shape)
    n_full = np.zeros(spout.rows.shape)        # pairs in which a window had particles
    n_any = 0                                  # pairs with particles in any window
    t_pair, wall_piv, wall_vis, empty_frac = [], [], [], []
    prev = None
    for n_pair, (t, f) in enumerate(frames, start=-1):
        d = read_vtk_xml(f)
        xyz, dia = d["coordinates"][:, cols_xyz], d["Diameter"]
        img_s, _ = render(xyz, dia, (spout.x0, spout.y0), spout.shape, px)
        img_w, own_w = render(xyz, dia, (wall.x0, wall.y0), wall.shape, px)
        if prev is not None:
            tp, ps, pw, own_p, vel_p = prev
            dt = t - tp
            dr, dc, emp = piv_pair(ps, img_s, spout.rows, spout.cols, **kw)
            if mode == "interpolate":
                # PIVlab-like: a window without particles gets no vector of
                # its own and is filled from its neighbours in this pair, so
                # every window carries a vector in every pair that has any
                # particles.  A pair with none anywhere contributes nothing.
                if emp.all():
                    dr = dc = np.zeros(spout.rows.shape)
                else:
                    dr = inpaint_nans(np.where(emp, np.nan, dr))
                    dc = inpaint_nans(np.where(emp, np.nan, dc))
                    n_any += 1
            su += dc * px / dt
            sv += dr * px / dt
            empty_frac.append(emp.mean())
            n_full += ~emp
            dr_w, dc_w, emp_w = piv_pair(pw, img_w, wall.rows, wall.cols, **kw)
            wall_piv.append(-(dr_w[in_roi] * px / dt).mean())      # downward +
            if on_pair is not None:
                on_pair({"i": n_pair, "t0": tp, "t1": t, "spout": (ps, img_s, dr, dc, emp),
                         "wall": (pw, img_w, dr_w, dc_w, emp_w)})
            o = own_p[roi_pix].ravel()
            o = o[o >= 0]
            wall_vis.append(-vel_p[o, 1].mean() if o.size else np.nan)
            t_pair.append(0.5 * (tp + t))
        prev = (t, img_s, img_w, own_w, d["Velocity"])

    n = len(t_pair)
    # empty windows add 0 to su, sv: 'zero' averages over all pairs, 'skip'
    # over the pairs in which the window had particles, 'interpolate' (filled
    # per pair, see above) over the pairs with particles anywhere
    n_avg = {"zero": n, "skip": np.maximum(n_full, 1), "interpolate": max(n_any, 1)}[mode]
    mean_map = np.hypot(su / n_avg, sv / n_avg)        # |time mean vector| [m/s]
    grid_map = fig4_grid_map(mean_map, spout, y_map, H, R_c, nrow, ncol)

    res = {
        "case": str(m["case_id"]), "run_dir": str(run), "flow_slpm": flow,
        "window_t0_s": frames[0][0], "window_t1_s": frames[-1][0], "n_pairs": n,
        "bed_surface_m": y_s, "map_bottom_m": y_map, "map_height_m": H, "panel_flow_slpm": pflow,
        "view_axis": cam["view_axis"], "empty_windows": mode,
        "pixel_m": px, "windows_px": wins, "step_px": step,
        "empty_window_frac": float(np.mean(empty_frac)),
        "map_piv_grid": mean_map, "map": grid_map,
        "wall_t_s": np.array(t_pair), "wall_piv_m_s": np.array(wall_piv),
        "wall_visible_m_s": np.array(wall_vis), "wall_roi_center_m": (xc, yc),
        "targets_table3": m.get("piv_targets_table3"),
    }
    if res["targets_table3"] and pflow == flow:
        exp, model = load_fig4("exp", flow), load_fig4("model", flow)
        res["fig4_exp"], res["fig4_model"] = exp, model
        res["metrics"] = map_metrics(exp, grid_map, p["metrics"])
        res["metrics_c26_model"] = map_metrics(exp, model, p["metrics"])
        res["wall_exp_m_s"] = p["wall_targets"]["exp_m_s"].get(int(flow))
        res["wall_c26_model_m_s"] = p["wall_targets"]["model_m_s"].get(int(flow))
    return res


# ----------------------------------------------------------------------------
# Output
# ----------------------------------------------------------------------------
def summary_row(res):
    row = {k: res[k] for k in ("case", "run_dir", "flow_slpm", "window_t0_s", "window_t1_s",
                               "n_pairs", "panel_flow_slpm", "empty_window_frac")}
    row.update({"bed_surface_mm": 1e3 * res["bed_surface_m"], "map_bottom_mm": 1e3 * res["map_bottom_m"],
                "map_height_mm": 1e3 * res["map_height_m"], "view_axis": res["view_axis"],
                "empty_windows": res["empty_windows"],
                "pixel_mm": 1e3 * res["pixel_m"], "windows_px": "/".join(map(str, res["windows_px"])),
                "step_px": res["step_px"],
                "map_mean_m_s": float(res["map"].mean()), "map_max_m_s": float(res["map"].max()),
                "wall_piv_m_s": float(np.mean(res["wall_piv_m_s"])),
                "wall_piv_sd_m_s": float(np.std(res["wall_piv_m_s"])),
                "wall_visible_parcels_m_s": float(np.nanmean(res["wall_visible_m_s"]))})
    t3 = res.get("targets_table3") or {}
    for key in ("CC", "SSI", "MSE_m2_s2"):
        row[key] = res.get("metrics", {}).get(key, np.nan)
        row[f"{key}_c26model_same_pipeline"] = res.get("metrics_c26_model", {}).get(key, np.nan)
        row[f"{key}_table3"] = t3.get(key, np.nan)
    row["wall_exp_m_s"] = res.get("wall_exp_m_s", np.nan)
    row["wall_c26model_m_s"] = res.get("wall_c26_model_m_s", np.nan)
    return row


def caption(res):
    return (f"Synthetic PIV, {res['n_pairs']} frame pairs over t = [{res['window_t0_s']:.3f}, "
            f"{res['window_t1_s']:.3f}] s (time averages over this window). Camera along {res['view_axis']}, "
            f"{1e3 * res['pixel_m']:.2f} mm/px; windows {'/'.join(map(str, res['windows_px']))} px, "
            f"step {res['step_px']} px. Map: |mean vector| (empty windows: {res['empty_windows']}) "
            f"from y = {1e3 * res['map_bottom_m']:.1f} mm "
            f"(bed surface {1e3 * res['bed_surface_m']:.1f} mm), {1e3 * res['map_height_m']:.1f} mm tall. "
            f"Wall: 4 mm box at y = {1e3 * res['wall_roi_center_m'][1]:.1f} mm, downward +.")


def plot(res, path, p):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    f4 = p["fig4"]
    vmin, vmax = f4["v_min_m_s"], f4["v_max_m_s"]
    R_mm = 0.5e3 * res["map_height_m"] * res["map"].shape[1] / res["map"].shape[0]
    ext = (-R_mm, R_mm, 0, 1e3 * res["map_height_m"])
    panels = [("this run (synthetic PIV)", res["map"], res.get("metrics"))]
    if "fig4_exp" in res:
        panels = [(f"C26 experiment, {res['flow_slpm']:g} slpm (Fig. 4a, digitized)", res["fig4_exp"], None),
                  ("C26 model (Fig. 4b, digitized)", res["fig4_model"], res["metrics_c26_model"])] + panels
    fig = plt.figure(figsize=(4.2 * len(panels), 7.2))
    gs = fig.add_gridspec(2, len(panels), height_ratios=(1.25, 1))
    for i, (title, a, met) in enumerate(panels):
        ax = fig.add_subplot(gs[0, i])
        im = ax.imshow(a, origin="lower", extent=ext, cmap="jet", vmin=vmin, vmax=vmax)
        if met:
            title += f"\nvs exp: CC {met['CC']:.3f}  SSI {met['SSI']:.3f}  MSE {met['MSE_m2_s2']:.4f}"
        ax.set_title(title, fontsize=8)
        ax.set_xlabel("x [mm]")
        if i == 0:
            ax.set_ylabel("height above map bottom [mm]")
    fig.colorbar(im, cax=fig.add_axes((0.9, 0.56, 0.012, 0.32)), label="velocity magnitude [m/s]")
    ax = fig.add_subplot(gs[1, :])
    t, w = res["wall_t_s"], res["wall_piv_m_s"]
    ax.plot(t, w, lw=0.6, color="0.6", label="synthetic PIV, 4 mm box")
    k = min(len(w), max(1, int(round(0.05 / np.median(np.diff(t))))))
    ax.plot(t[k - 1:], np.convolve(w, np.ones(k) / k, "valid"), lw=2, color="k",
            label=f"moving average 0.05 s; mean {np.mean(w):.2e} m/s")
    ax.axhline(np.nanmean(res["wall_visible_m_s"]), color="tab:green", ls=":",
               label=f"visible parcels in box, mean {np.nanmean(res['wall_visible_m_s']):.2e} m/s")
    if res.get("wall_exp_m_s"):
        ax.axhline(res["wall_exp_m_s"], color="tab:blue", ls="--",
                   label=f"C26 experiment mean {res['wall_exp_m_s']:.2e} m/s (Fig. 6d)")
        ax.axhline(res["wall_c26_model_m_s"], color="tab:red", ls="--",
                   label=f"C26 model mean {res['wall_c26_model_m_s']:.2e} m/s (Fig. 6d)")
    ax.set_xlabel("t [s]")
    ax.set_ylabel("wall v, downward + [m/s]")
    ax.set_yscale("symlog", linthresh=1e-3)
    ax.legend(fontsize=7, loc="upper right")
    ax.grid(alpha=0.3)
    fig.suptitle(f"C26 {res['case']} synthetic PIV vs C26 Fig. 4 / Fig. 6d")
    fig.text(0.01, 0.005, caption(res), fontsize=7, wrap=True)
    fig.subplots_adjust(left=0.06, right=0.88, top=0.9, bottom=0.11, hspace=0.35, wspace=0.25)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run", type=Path, help="run dir with manifest.json and PARTICLES_PIV.pvd")
    ap.add_argument("--out-dir", type=Path, default=REPO / "results")
    ap.add_argument("--tag", default="", help="appended to output file names")
    ap.add_argument("--stride", type=int, default=1, help="use every n-th frame (checks only)")
    ap.add_argument("--max-frames", type=int, default=None, help="cap on frames (checks only)")
    a = ap.parse_args(argv)
    p = load_params()
    t_start = time.time()
    res = analyse_run(a.run, p, a.stride, a.max_frames)
    run_name = a.run.parent.name if a.run.name == "aws_out" else a.run.name   # c01, c01_massonly
    name = f"{run_name}{'_' + a.tag if a.tag else ''}"
    a.out_dir.mkdir(parents=True, exist_ok=True)
    row = summary_row(res)
    with (a.out_dir / f"c26_piv_{name}.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(row))
        w.writeheader()
        w.writerow(row)
    np.savetxt(a.out_dir / f"c26_piv_map_{name}.csv", np.flipud(res["map"]), delimiter=",", fmt="%.5f",
               header=f"C26 {res['case']} synthetic-PIV spout map [m/s], top row first, x from -R_c to R_c. "
                      + caption(res))
    plot(res, a.out_dir / f"fig_c26_piv_{name}.png", p)
    print(caption(res))
    for k, v in row.items():
        print(f"  {k}: {v}")
    print(f"wrote {a.out_dir}/c26_piv_{name}.csv, c26_piv_map_{name}.csv, fig_c26_piv_{name}.png "
          f"({time.time() - t_start:.0f} s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

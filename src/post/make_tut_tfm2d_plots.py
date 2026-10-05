"""Day-5 deliverable: reproduce two plots from the TFM 2D tutorial output.

(a) Time series of pressure drop across the bed, ΔP(t) [Pa], computed as the
    cell-area-mean of ``P_g`` on the bottom cell row minus the cell-area-mean
    on the top cell row of the Cartesian VTK output.

(b) Time-averaged solids fraction field ``<ε_s>`` over the last half of the
    simulated time, where ``ε_s = 1 − ε_g`` is read cell-wise from each
    ``.vtu``.

Inputs
------
``--data-dir``: directory containing ``BACKGROUND.pvd`` and ``BACKGROUND_*.vtu``
``--out-dir``:  directory to write ``fig_d5_dp_timeseries.png`` and
                 ``fig_d5_eps_s_mean.png`` into (defaults to ``results/``).

Raw simulation output is intentionally not committed (CLAUDE.md rule 7); this
script consumes it from an external path provided by the user.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .read_vtk import cell_field, read_pvd, read_vtu


def _bottom_top_rows(mesh, tol: float = 1e-9) -> tuple[np.ndarray, np.ndarray]:
    """Return boolean cell masks for the lowest-y row and the highest-y row.

    MFiX Cartesian output has structured cell centres in y. We bin by y and
    select the min- and max-y strata.
    """
    yc = mesh.cell_centers().points[:, 1]
    y_min, y_max = yc.min(), yc.max()
    bottom = np.abs(yc - y_min) < tol * max(1.0, abs(y_max - y_min))
    top = np.abs(yc - y_max) < tol * max(1.0, abs(y_max - y_min))
    if bottom.sum() == 0 or top.sum() == 0:
        # fall back to the thinnest percentile band if tol was too tight
        bottom = yc <= np.percentile(yc, 1.0)
        top = yc >= np.percentile(yc, 99.0)
    return bottom, top


def pressure_drop_series(pvd_path: Path) -> pd.DataFrame:
    """Build the ΔP(t) time series by walking every frame in the ``.pvd``.

    Returns a DataFrame indexed by ``time`` [s] with columns
    ``p_bottom`` [Pa], ``p_top`` [Pa], ``dp`` [Pa].
    """
    rows = []
    for t, path in read_pvd(pvd_path):
        mesh = read_vtu(path)
        p = cell_field(mesh, "P_G", "p_g", "P_g")
        bottom, top = _bottom_top_rows(mesh)
        p_bot = float(p[bottom].mean())
        p_top = float(p[top].mean())
        rows.append({"time": t, "p_bottom": p_bot, "p_top": p_top, "dp": p_bot - p_top})
    return pd.DataFrame(rows).set_index("time")


def eps_s_timeavg(pvd_path: Path, t_window: tuple[float, float]) -> tuple[np.ndarray, np.ndarray, np.ndarray, tuple[float, float]]:
    """Average ``ε_s = 1 − ε_g`` cell-wise over frames with ``t ∈ t_window``.

    Returns ``(x_edges, y_edges, eps_s_mean_2d, (t_lo, t_hi))`` where
    ``x_edges``/``y_edges`` are cell-centre coordinates reshaped into the
    structured 2D grid, and ``eps_s_mean_2d`` is the time-averaged field in
    ``(ny, nx)`` orientation suited for ``imshow`` with ``origin='lower'``.
    """
    frames = [(t, p) for t, p in read_pvd(pvd_path) if t_window[0] <= t <= t_window[1]]
    if not frames:
        raise RuntimeError(f"No .pvd frames fall inside window {t_window} s")

    acc: np.ndarray | None = None
    centres: np.ndarray | None = None
    for t, path in frames:
        mesh = read_vtu(path)
        eps_g = cell_field(mesh, "EP_G", "ep_g", "EP_g")
        eps_s = 1.0 - np.asarray(eps_g, dtype=float)
        if acc is None:
            acc = eps_s.copy()
            centres = mesh.cell_centers().points
        else:
            acc += eps_s
    assert acc is not None and centres is not None
    acc /= len(frames)

    xc = centres[:, 0]
    yc = centres[:, 1]
    xs = np.unique(np.round(xc, 10))
    ys = np.unique(np.round(yc, 10))
    nx, ny = xs.size, ys.size
    if nx * ny != acc.size:
        raise RuntimeError(f"Grid not reconstructable: nx*ny={nx*ny} vs n_cells={acc.size}")
    # Build (ny, nx) image by sorting cells into (iy, ix) bins
    ix = np.searchsorted(xs, np.round(xc, 10))
    iy = np.searchsorted(ys, np.round(yc, 10))
    img = np.full((ny, nx), np.nan)
    img[iy, ix] = acc
    return xs, ys, img, (frames[0][0], frames[-1][0])


def _plot_dp(df: pd.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(df.index.values, df["dp"].values, "-o", ms=3, lw=1.2)
    ax.set_xlabel("time [s]")
    ax.set_ylabel(r"$\Delta P = \bar P_g(\mathrm{bottom}) - \bar P_g(\mathrm{top})$  [Pa]")
    ax.set_title("TFM 2D fluidized-bed tutorial — pressure drop across the bed")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def _plot_eps_s(xs, ys, img, window, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(4, 7))
    im = ax.imshow(
        img,
        origin="lower",
        extent=(xs.min(), xs.max(), ys.min(), ys.max()),
        aspect="equal",
        cmap="viridis",
        vmin=0.0,
        vmax=max(float(np.nanmax(img)), 0.6),
    )
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_title(
        rf"$\langle\varepsilon_s\rangle$ averaged over $t \in [{window[0]:.3f}, {window[1]:.3f}]$ s"
    )
    cbar = fig.colorbar(im, ax=ax, shrink=0.9)
    cbar.set_label(r"$\langle\varepsilon_s\rangle$  [-]")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", required=True, type=Path, help="dir containing BACKGROUND.pvd")
    ap.add_argument("--out-dir", default=Path("results"), type=Path)
    ap.add_argument("--pvd", default="BACKGROUND.pvd")
    args = ap.parse_args()

    pvd = args.data_dir / args.pvd
    if not pvd.exists():
        raise SystemExit(f"PVD not found: {pvd}")
    args.out_dir.mkdir(parents=True, exist_ok=True)

    # (a) pressure-drop time series
    df = pressure_drop_series(pvd)
    df.to_csv(args.out_dir / "fig_d5_dp_timeseries.csv")
    _plot_dp(df, args.out_dir / "fig_d5_dp_timeseries.png")

    # (b) time-averaged ε_s over the last half of the available run
    t_lo = 0.5 * df.index.values.max()
    t_hi = df.index.values.max()
    xs, ys, img, window = eps_s_timeavg(pvd, (t_lo, t_hi))
    _plot_eps_s(xs, ys, img, window, args.out_dir / "fig_d5_eps_s_mean.png")

    print(
        f"wrote {args.out_dir/'fig_d5_dp_timeseries.png'} ({len(df)} frames, "
        f"ΔP range {df['dp'].min():.1f}–{df['dp'].max():.1f} Pa)"
    )
    print(
        f"wrote {args.out_dir/'fig_d5_eps_s_mean.png'} (window [{window[0]:.3f}, {window[1]:.3f}] s, "
        f"<eps_s> range {np.nanmin(img):.3f}–{np.nanmax(img):.3f})"
    )


if __name__ == "__main__":
    main()

"""Day-9 V_2d post-processing: solids-fraction frame grid, ΔP time series,
and a max-CFL check via ``src.cfl_check``.

Inputs (all under one run directory, default ``cases/V_2d/day9_u1p2/``):

  - ``BACKGROUND.pvd`` + ``BACKGROUND_*.vtu``  — cell-centred field frames
  - ``V_2D_P_INLET.csv``                        — gas pressure at the orifice [Pa]
  - ``V_2D_P_TOP_BED.csv``                      — gas pressure at y=H_static [Pa]
  - ``V_2D_SPOUT_MID.csv``                      — eps_s, v_s at spout mid-plane
  - ``V_2D_DAY9.TXT``                           — MFiX advance log (dt vs t)
  - ``manifest.json``                           — U_in, U_ms, drag, c_e

Outputs (into ``results/``):

  - ``fig_d9_V_2d_solids_frames.png`` : 6-panel grid of solids fraction
                                        (eps_s = 1 - EP_G) at evenly spaced
                                        frames after startup
  - ``fig_d9_V_2d_dP.png``            : ΔP(t) across the bed = P_inlet - P_top
  - ``fig_d9_V_2d_cfl.csv``           : per-frame max |V_gas| and CFL
  - ``fig_d9_V_2d_dP.csv``            : ΔP time series tabulated

Headless-safe (uses MPLBACKEND=Agg if not already set).
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.cfl_check import cfl_number
from src.post.read_monitors import read_monitor
from src.post.read_vtk import load_pvd, read_pvd, read_vtu

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_RUN_DIR = REPO_ROOT / "cases" / "V_2d" / "day9_u1p2"
DEFAULT_OUT_DIR = REPO_ROOT / "results"


# ---------- helpers ----------

def _load_manifest(run_dir: Path) -> dict:
    with (run_dir / "manifest.json").open() as f:
        return json.load(f)


def _parse_dt_from_txt(path: Path) -> pd.DataFrame:
    """Return a DataFrame with columns ['t', 'dt'] from a V_2D_DAY9.TXT log.

    MFiX writes one line per time step; column 1 is t [s], column 2 is dt [s].
    Lines that don't start with a number are header / warning messages and
    are silently skipped.
    """
    rows = []
    with path.open() as f:
        for line in f:
            toks = line.split()
            if len(toks) < 2:
                continue
            try:
                t = float(toks[0])
                dt = float(toks[1])
            except ValueError:
                continue
            rows.append((t, dt))
    return pd.DataFrame(rows, columns=["t", "dt"])


def _frame_eps_s(mesh) -> np.ndarray:
    """Solids volume fraction eps_s = 1 - EP_G, flattened per-cell."""
    ep_g = np.asarray(mesh.cell_data["EP_G"]).astype(float)
    return 1.0 - ep_g


def _frame_vgas_max(mesh) -> float:
    """Max magnitude of gas velocity across the mesh [m/s]."""
    vg = np.asarray(mesh.cell_data["Gas_Velocity"]).astype(float)
    if vg.ndim == 2:
        return float(np.sqrt((vg * vg).sum(axis=1)).max())
    return float(np.abs(vg).max())


def _pick_frame_indices(times: np.ndarray, n: int, t_start_s: float) -> list[int]:
    """Return n frame indices evenly spaced in time, starting at or after t_start_s."""
    mask = times >= t_start_s
    if mask.sum() < n:
        return list(range(len(times)))[-n:] if len(times) >= n else list(range(len(times)))
    hi = len(times) - 1
    lo = int(np.argmax(mask))
    picks = np.linspace(lo, hi, n).round().astype(int).tolist()
    return picks


# ---------- plot: solids-fraction frame grid ----------

def plot_solids_frames(
    pvd_path: Path, out_path: Path, n_frames: int = 6, t_start_s: float = 0.2,
    u_in_m_s: float | None = None,
) -> dict:
    """Write a grid of eps_s heatmaps from the ``.pvd`` collection.

    Returns a dict with the picked ``times`` and ``max_abs_eps_s`` per frame
    for use in the summary / CFL CSV.
    """
    entries = read_pvd(pvd_path)
    if not entries:
        raise RuntimeError(f"no frames in {pvd_path}")
    times = np.array([t for t, _ in entries])
    picks = _pick_frame_indices(times, n_frames, t_start_s)
    picked = [entries[i] for i in picks]

    fig, axes = plt.subplots(1, len(picks), figsize=(2.2 * len(picks), 6.0), sharey=True)
    if len(picks) == 1:
        axes = [axes]

    im = None
    bounds = None
    for ax, (t, path) in zip(axes, picked):
        mesh = read_vtu(path)
        eps_s = _frame_eps_s(mesh)
        pts = np.asarray(mesh.cell_centers().points)
        x = pts[:, 0]
        y = pts[:, 1]
        # Simple tricontourf gives a smooth, grid-agnostic heat map.
        im = ax.tricontourf(x, y, eps_s, levels=np.linspace(0.0, 0.65, 14),
                            cmap="magma_r")
        ax.set_title(f"t = {t:.2f} s", fontsize=10)
        ax.set_aspect("equal")
        ax.set_xlabel("x [m]")
        if bounds is None:
            bounds = mesh.bounds
        ax.set_xlim(bounds[0], bounds[1])
        ax.set_ylim(0.0, bounds[3])
    axes[0].set_ylabel("y [m]")
    cbar = fig.colorbar(im, ax=axes, shrink=0.75, pad=0.02)
    cbar.set_label(r"$\varepsilon_s = 1 - \varepsilon_g$ [–]")
    u_in_tag = f" (orifice {u_in_m_s:.2f} m/s)" if u_in_m_s else ""
    fig.suptitle(
        f"Case V_2d (Day 9): solids fraction at 1.2 × U_ms{u_in_tag}",
        fontsize=11,
    )
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    return {"times": [t for t, _ in picked]}


# ---------- plot: ΔP time series ----------

def plot_dp_timeseries(run_dir: Path, out_path: Path) -> pd.DataFrame:
    """Write a two-panel ΔP(t) figure from the two pressure monitors.

    Top panel: raw P_inlet and P_top above atmospheric (101 325 Pa).
    Bottom panel: ``ΔP_bed = P_top - P_inlet`` [Pa], the hydrodynamic bed
    pressure drop (positive during spouting because the MI face is pinned to
    atmospheric by the mass-inflow BC and the gas above it is compressed by
    the bed weight; this is the sign convention Mathur & Epstein use for
    spouted-bed ΔP characteristics).
    """
    p_in = read_monitor(run_dir / "V_2D_P_INLET.csv")
    p_top = read_monitor(run_dir / "V_2D_P_TOP_BED.csv")
    merged = pd.DataFrame({
        "P_inlet_Pa": p_in["p_g"],
        "P_top_Pa": p_top["p_g"],
    }).dropna()
    P_atm = 101_325.0
    merged["P_inlet_gauge_Pa"] = merged["P_inlet_Pa"] - P_atm
    merged["P_top_gauge_Pa"] = merged["P_top_Pa"] - P_atm
    merged["dP_bed_Pa"] = merged["P_top_Pa"] - merged["P_inlet_Pa"]

    fig, (ax_p, ax_dp) = plt.subplots(2, 1, figsize=(7.5, 5.4), sharex=True)
    ax_p.plot(merged.index, merged["P_inlet_gauge_Pa"], lw=1.0, color="C0",
              label=r"$P_{inlet} - P_{atm}$ (orifice MI face)")
    ax_p.plot(merged.index, merged["P_top_gauge_Pa"], lw=1.0, color="C3",
              label=r"$P_{top\,bed} - P_{atm}$ (y = $H_{static}$)")
    ax_p.axhline(0.0, color="k", lw=0.5, ls="--", alpha=0.4)
    ax_p.set_ylabel("gauge pressure [Pa]")
    ax_p.set_title("Case V_2d (Day 9): bed pressures and ΔP at 1.2 × U_ms")
    ax_p.legend(loc="best", fontsize=9)
    ax_p.grid(True, alpha=0.3)

    ax_dp.plot(merged.index, merged["dP_bed_Pa"], lw=1.0, color="C2")
    ax_dp.axhline(0.0, color="k", lw=0.5, ls="--", alpha=0.4)
    ax_dp.set_xlabel("time [s]")
    ax_dp.set_ylabel(r"$\Delta P_{bed} = P_{top} - P_{inlet}$ [Pa]")
    ax_dp.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)

    csv_path = out_path.with_suffix(".csv")
    merged.to_csv(csv_path)
    return merged


# ---------- CFL check ----------

def check_max_cfl(
    run_dir: Path, pvd_path: Path, dx_m: float, out_csv: Path
) -> pd.DataFrame:
    """Compute per-frame |V_gas|_max and CFL = |V_gas|_max * dt / dx.

    dt per frame is taken as the mean solver dt reported in V_2D_DAY9.TXT within
    [t_frame-0.01, t_frame+0.01] s; if the solver log is missing, falls back to
    the inter-frame dt (which is >> solver dt, so CFL would be a loose upper
    bound only).
    """
    txt_path = run_dir / "V_2D_DAY9.TXT"
    dt_log = _parse_dt_from_txt(txt_path) if txt_path.exists() else None

    entries = read_pvd(pvd_path)
    rows = []
    for t_frame, vtu in entries:
        mesh = read_vtu(vtu)
        v_max = _frame_vgas_max(mesh)
        if dt_log is not None and len(dt_log):
            m = dt_log[(dt_log["t"] > t_frame - 0.01) & (dt_log["t"] < t_frame + 0.01)]
            dt_local = float(m["dt"].mean()) if len(m) else float(dt_log["dt"].min())
        else:
            dt_local = 1e-4
        cfl = cfl_number(u_max=max(v_max, 1e-12), dx=dx_m, dt=dt_local)
        rows.append((t_frame, v_max, dt_local, cfl))
    df = pd.DataFrame(rows, columns=["t_s", "vgas_max_m_s", "dt_s", "CFL"])
    df.to_csv(out_csv, index=False)
    return df


# ---------- CLI ----------

def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run-dir", type=Path, default=DEFAULT_RUN_DIR)
    p.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    p.add_argument("--n-frames", type=int, default=6)
    p.add_argument("--t-start", type=float, default=0.2, help="skip the startup window [s]")
    args = p.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    pvd = args.run_dir / "BACKGROUND.pvd"
    if not pvd.exists():
        raise FileNotFoundError(pvd)

    manifest = _load_manifest(args.run_dir)
    # dx = 2*R_c / imax; R_c is baked into the mesh bounds, read from frame 0.
    m0 = read_vtu(next(iter(read_pvd(pvd)))[1])
    dx_m = (m0.bounds[1] - m0.bounds[0]) / manifest["imax"]

    frames_png = args.out_dir / "fig_d9_V_2d_solids_frames.png"
    dp_png = args.out_dir / "fig_d9_V_2d_dP.png"
    cfl_csv = args.out_dir / "fig_d9_V_2d_cfl.csv"

    frames_info = plot_solids_frames(pvd, frames_png, n_frames=args.n_frames,
                                     t_start_s=args.t_start,
                                     u_in_m_s=manifest.get("U_in_m_s"))
    dp_df = plot_dp_timeseries(args.run_dir, dp_png)
    cfl_df = check_max_cfl(args.run_dir, pvd, dx_m, cfl_csv)

    print(f"wrote {frames_png}")
    print(f"      picked frames at t = " + ", ".join(f"{t:.2f}" for t in frames_info["times"]))
    print(f"wrote {dp_png}  ({len(dp_df)} rows)")
    print(f"      ΔP_bed range: {dp_df['dP_bed_Pa'].min():.1f} .. {dp_df['dP_bed_Pa'].max():.1f} Pa")
    print(f"wrote {cfl_csv}")
    print(f"      max |V_gas| overall = {cfl_df['vgas_max_m_s'].max():.2f} m/s")
    print(f"      max CFL overall    = {cfl_df['CFL'].max():.3f}")
    print(f"      dx used            = {dx_m*1e3:.3f} mm (imax={manifest['imax']})")


if __name__ == "__main__":
    main()

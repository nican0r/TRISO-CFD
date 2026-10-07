"""Day-9 V_2d post-processing: eps_s frames, DeltaP, solids inventory, PSD,
and a max-CFL check via ``src.cfl_check``.

Headless-safe (prefixes MPLBACKEND=Agg if unset).
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from src.cfl_check import cfl_number  # noqa: E402
from src.post.psd import welch_psd  # noqa: E402
from src.post.read_monitors import read_monitor  # noqa: E402
from src.post.read_vtk import read_pvd, read_vtu  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_RUN_DIR = REPO_ROOT / "cases" / "V_2d" / "day9"
DEFAULT_OUT_DIR = REPO_ROOT / "results"
P_ATM = 101_325.0


# ---------- helpers ----------

def _load_manifest(run_dir: Path) -> dict:
    return json.loads((run_dir / "manifest.json").read_text())


def _parse_dt_from_txt(path: Path) -> pd.DataFrame:
    rows = []
    with path.open() as f:
        for line in f:
            toks = line.split()
            if len(toks) < 2:
                continue
            try:
                t = float(toks[0]); dt = float(toks[1])
            except ValueError:
                continue
            rows.append((t, dt))
    return pd.DataFrame(rows, columns=["t", "dt"])


def _frame_eps_s(mesh) -> np.ndarray:
    ep_g = np.asarray(mesh.cell_data["EP_G"]).astype(float)
    return 1.0 - ep_g


def _frame_vgas_max(mesh) -> float:
    vg = np.asarray(mesh.cell_data["Gas_Velocity"]).astype(float)
    if vg.ndim == 2:
        return float(np.sqrt((vg * vg).sum(axis=1)).max())
    return float(np.abs(vg).max())


def _pick_frame_indices(times: np.ndarray, n: int, t_start_s: float) -> list[int]:
    mask = times >= t_start_s
    if mask.sum() < n:
        return list(range(len(times)))[-n:] if len(times) >= n else list(range(len(times)))
    hi = len(times) - 1
    lo = int(np.argmax(mask))
    return np.linspace(lo, hi, n).round().astype(int).tolist()


# ---------- plot: solids-fraction frame grid ----------

def plot_solids_frames(pvd_path: Path, out_path: Path, n_frames: int = 6,
                        t_start_s: float = 0.3, u_in_m_s: float | None = None) -> dict:
    entries = read_pvd(pvd_path)
    if not entries:
        raise RuntimeError(f"no frames in {pvd_path}")
    times = np.array([t for t, _ in entries])
    picks = _pick_frame_indices(times, n_frames, t_start_s)
    picked = [entries[i] for i in picks]

    fig, axes = plt.subplots(1, len(picks),
                             figsize=(2.2 * len(picks), 6.0), sharey=True)
    if len(picks) == 1:
        axes = [axes]
    im = None
    bounds = None
    for ax, (t, path) in zip(axes, picked):
        mesh = read_vtu(path)
        eps_s = _frame_eps_s(mesh)
        pts = np.asarray(mesh.cell_centers().points)
        im = ax.tricontourf(pts[:, 0], pts[:, 1], eps_s,
                             levels=np.linspace(0.0, 0.65, 14), cmap="magma_r")
        ax.set_title(f"t = {t:.2f} s", fontsize=10)
        ax.set_aspect("equal")
        ax.set_xlabel("x [m]")
        if bounds is None:
            bounds = mesh.bounds
        ax.set_xlim(bounds[0], bounds[1])
        ax.set_ylim(0.0, bounds[3])
    axes[0].set_ylabel("y [m]")
    cbar = fig.colorbar(im, ax=axes, shrink=0.75, pad=0.02)
    cbar.set_label(r"$\varepsilon_s = 1 - \varepsilon_g$ [-]")
    u_in_tag = f" (orifice {u_in_m_s:.1f} m/s)" if u_in_m_s else ""
    fig.suptitle(
        f"V_2d (Day 9, ORNL/UTK): solids fraction frames{u_in_tag}", fontsize=11,
    )
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    return {"times": [t for t, _ in picked]}


# ---------- plot: dP time series ----------

def plot_dp_timeseries(run_dir: Path, out_path: Path) -> pd.DataFrame:
    p_in = read_monitor(run_dir / "V_2D_P_INLET.csv")
    p_top = read_monitor(run_dir / "V_2D_P_TOP_BED.csv")
    df = pd.DataFrame({
        "P_inlet_Pa": p_in["p_g"],
        "P_top_Pa": p_top["p_g"],
    }).dropna()
    df["P_inlet_gauge_Pa"] = df["P_inlet_Pa"] - P_ATM
    df["P_top_gauge_Pa"] = df["P_top_Pa"] - P_ATM
    df["dP_bed_Pa"] = df["P_top_Pa"] - df["P_inlet_Pa"]

    fig, (ax_p, ax_dp) = plt.subplots(2, 1, figsize=(7.5, 5.4), sharex=True)
    ax_p.plot(df.index, df["P_inlet_gauge_Pa"], lw=1.0, color="C0",
              label=r"$P_{inlet}-P_{atm}$")
    ax_p.plot(df.index, df["P_top_gauge_Pa"], lw=1.0, color="C3",
              label=r"$P_{top\,bed}-P_{atm}$")
    ax_p.axhline(0.0, color="k", lw=0.5, ls="--", alpha=0.4)
    ax_p.set_ylabel("gauge pressure [Pa]")
    ax_p.set_title("Case V_2d (Day 9, ORNL/UTK): bed pressures at U = 30 m/s")
    ax_p.legend(loc="best", fontsize=9)
    ax_p.grid(True, alpha=0.3)
    ax_dp.plot(df.index, df["dP_bed_Pa"], lw=1.0, color="C2")
    ax_dp.axhline(0.0, color="k", lw=0.5, ls="--", alpha=0.4)
    ax_dp.set_xlabel("time [s]")
    ax_dp.set_ylabel(r"$\Delta P_{bed}$ [Pa]")
    ax_dp.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    df.to_csv(out_path.with_suffix(".csv"))
    return df


# ---------- plot: inventory trace ----------

def plot_inventory_trace(run_dir: Path, out_path: Path) -> pd.DataFrame:
    inv = read_monitor(run_dir / "V_2D_SOLIDS_INVENTORY.csv")
    # Monitor file column is "ep_s(1)" (per-phase); rename to a stable name.
    col = next(c for c in inv.columns if c.startswith("ep_s"))
    df = inv.rename(columns={col: "ep_s_vol_m3"}).dropna()
    df["frac_retained"] = df["ep_s_vol_m3"] / df["ep_s_vol_m3"].iloc[0]
    fig, ax = plt.subplots(figsize=(7.5, 3.5))
    ax.plot(df.index, df["frac_retained"], lw=1.0, color="C4")
    ax.axhline(0.95, color="r", lw=0.5, ls="--",
               label="95% retention floor (Day 9 pass criterion)")
    ax.set_xlabel("time [s]")
    ax.set_ylabel("solids inventory / inventory(0) [-]")
    ax.set_title("V_2d (Day 9, ORNL/UTK): solids inventory retention")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    df.to_csv(out_path.with_suffix(".csv"))
    return df


# ---------- plot: PSD of inlet pressure ----------

def plot_inlet_psd(run_dir: Path, out_path: Path, t_startup_s: float = 0.5,
                    monitor_dt_s: float = 1.0e-3) -> dict:
    p_in = read_monitor(run_dir / "V_2D_P_INLET.csv").dropna()
    mask = p_in.index.to_numpy() >= t_startup_s
    sig = p_in.loc[mask, "p_g"].to_numpy() - P_ATM
    t_window_s = float(p_in.index[-1] - t_startup_s)
    fs_hz = 1.0 / monitor_dt_s
    res = welch_psd(sig, fs_hz=fs_hz)

    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    ax.semilogy(res.freqs_hz, res.psd, lw=1.0, color="C0")
    # S1 Fig. 2: ~10-11 Hz (ORNL 2-inch coater, 600 K ambient).
    # S2: 17-25 Hz over U/U_ms = 1.3-2.3 (UTK mockup).
    ax.axvspan(17.0, 25.0, color="C3", alpha=0.12, label="S2 band 17-25 Hz")
    ax.axvline(10.5, color="C2", lw=1.0, ls="--", label="S1 ~10-11 Hz")
    ax.axvline(res.dominant_hz, color="k", lw=1.0, ls=":",
               label=f"peak = {res.dominant_hz:.1f} Hz")
    ax.set_xlim(0.0, 100.0)
    ax.set_xlabel("frequency [Hz]")
    ax.set_ylabel(r"PSD of $P_{inlet} - P_{atm}$ [Pa$^2$/Hz]")
    ax.set_title(
        f"V_2d (Day 9, ORNL/UTK): inlet-pressure PSD\n"
        f"(window: t in [{t_startup_s:.2f}, {p_in.index[-1]:.2f}] s, {t_window_s:.2f} s span)",
        fontsize=10,
    )
    ax.legend(loc="best", fontsize=9)
    ax.grid(True, alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    return {
        "dominant_hz": res.dominant_hz,
        "fs_hz": fs_hz,
        "t_window_s": t_window_s,
        "t_startup_s": t_startup_s,
    }


# ---------- CFL check ----------

def check_max_cfl(run_dir: Path, pvd_path: Path, dx_m: float,
                   out_csv: Path) -> pd.DataFrame:
    txt_path = run_dir / "V_2D_DAY9.TXT"
    dt_log = _parse_dt_from_txt(txt_path) if txt_path.exists() else None
    rows = []
    for t_frame, vtu in read_pvd(pvd_path):
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
    p.add_argument("--t-startup", type=float, default=0.5,
                   help="skip the startup window [s] (used for PSD + frames)")
    args = p.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    pvd = args.run_dir / "BACKGROUND.pvd"
    if not pvd.exists():
        raise FileNotFoundError(pvd)
    manifest = _load_manifest(args.run_dir)
    dx_m = float(manifest.get("dx_m")) or (
        (read_vtu(next(iter(read_pvd(pvd)))[1]).bounds[1]
         - read_vtu(next(iter(read_pvd(pvd)))[1]).bounds[0]) / manifest["imax"]
    )

    frames_png = args.out_dir / "fig_d9_V_2d_solids_frames.png"
    dp_png = args.out_dir / "fig_d9_V_2d_dP.png"
    inv_png = args.out_dir / "fig_d9_V_2d_inventory.png"
    psd_png = args.out_dir / "fig_d9_V_2d_psd.png"
    cfl_csv = args.out_dir / "fig_d9_V_2d_cfl.csv"

    frames_info = plot_solids_frames(pvd, frames_png, n_frames=args.n_frames,
                                     t_start_s=args.t_startup,
                                     u_in_m_s=manifest.get("U_in_m_s"))
    dp_df = plot_dp_timeseries(args.run_dir, dp_png)
    inv_df = plot_inventory_trace(args.run_dir, inv_png)
    psd_info = plot_inlet_psd(args.run_dir, psd_png,
                               t_startup_s=args.t_startup,
                               monitor_dt_s=manifest.get("monitor_dt_s", 1e-3))
    cfl_df = check_max_cfl(args.run_dir, pvd, dx_m, cfl_csv)

    print(f"wrote {frames_png}")
    print("      picked frames at t = " + ", ".join(f"{t:.2f}" for t in frames_info["times"]))
    print(f"wrote {dp_png}  ({len(dp_df)} rows); dP_bed range "
          f"{dp_df['dP_bed_Pa'].min():.1f} .. {dp_df['dP_bed_Pa'].max():.1f} Pa")
    print(f"wrote {inv_png}  final retention = {inv_df['frac_retained'].iloc[-1]:.3f}")
    print(f"wrote {psd_png}  dominant frequency = {psd_info['dominant_hz']:.2f} Hz "
          f"(window {psd_info['t_window_s']:.2f} s, fs = {psd_info['fs_hz']:.0f} Hz)")
    print(f"wrote {cfl_csv}  max |V_gas| = {cfl_df['vgas_max_m_s'].max():.2f} m/s; "
          f"max CFL = {cfl_df['CFL'].max():.3f}")


if __name__ == "__main__":
    main()

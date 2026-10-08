"""Day-10 V_3d post-processing: time-averaged centerline v_g vs height
overlay (S1 Fig. 1), inlet-pressure Welch PSD with S2's 17-25 Hz band
and S1's ~10-11 Hz line, plus a max-CFL check via ``src.cfl_check``.

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
DEFAULT_RUN_DIR = REPO_ROOT / "cases" / "V_3d" / "day10"
DEFAULT_OUT_DIR = REPO_ROOT / "results"
DEFAULT_DIGITIZED = REPO_ROOT / "data" / "digitized" / "s1_fig1_centerline_vg.csv"
P_ATM = 101_325.0


def _load_manifest(run_dir: Path) -> dict:
    return json.loads((run_dir / "manifest.json").read_text())


def _manual_cell_centers(mesh) -> np.ndarray:
    """Return cell centroids as a (n_cells, 3) ndarray WITHOUT using
    PyVista's ``cell_centers()`` (which hangs on MFiX cut-cell VTUs on
    this box).  Centroids are the mean of each cell's corner points,
    read directly from ``mesh.points`` and ``mesh.cells``.  The cells
    array is flat: ``[n_pts_0, pt0_0, pt0_1, ..., n_pts_1, pt1_0, ...]``.
    """
    pts = np.asarray(mesh.points)
    cells = np.asarray(mesh.cells)
    centers = np.zeros((mesh.n_cells, 3))
    idx = 0
    cid = 0
    while idx < len(cells):
        n = int(cells[idx])
        pt_ids = cells[idx + 1:idx + 1 + n]
        centers[cid] = pts[pt_ids].mean(axis=0)
        cid += 1
        idx += 1 + n
    return centers


def plot_centerline_vg(pvd_path: Path, out_path: Path,
                       digitized_csv: Path | None,
                       r_select_m: float,
                       t_startup_s: float = 1.0,
                       u_in_m_s: float | None = None) -> dict:
    """Time-averaged centerline v_g vs y, overlaid against S1 Fig. 1.

    The mesh is static across frames, so cell centres and the near-axis
    mask are computed once on the first frame and reused for every
    subsequent frame.  Cells with sqrt(x^2 + z^2) <= r_select_m form the
    central column; v_g is averaged within each y-row at each frame, then
    time-averaged over frames after the startup window.

    MFiX cut-cell VTUs report centroids shifted into the fluid portion of
    the cell, so the near-axis cells (nominal centre at r <= 1 mm) have
    centroid r > R_i.  Caller should set r_select_m >= R_i_m + dx_min_m
    (orifice radius + near-axis cell size) to capture the inner column.
    """
    entries = read_pvd(pvd_path)
    if not entries:
        raise RuntimeError(f"no frames in {pvd_path}")
    frames = [(t, p) for t, p in entries if t >= t_startup_s]
    if not frames:
        raise RuntimeError(
            f"no frames past t_startup = {t_startup_s} s in {pvd_path}"
        )
    # Precompute static-mesh selection ONCE from the first post-startup frame.
    mesh0 = read_vtu(frames[0][1])
    pts0 = _manual_cell_centers(mesh0)
    r0 = np.sqrt(pts0[:, 0] ** 2 + pts0[:, 2] ** 2)
    sel = r0 <= r_select_m
    if not sel.any():
        raise RuntimeError("no cells within r_select_m of the axis")
    y_sel = pts0[sel, 1]
    # Group cells into y-rows (same mesh across frames).
    order = np.argsort(y_sel)
    y_sorted = y_sel[order]
    row_starts = [0]
    for i in range(1, len(y_sorted)):
        if (y_sorted[i] - y_sorted[row_starts[-1]]) > 1e-6:
            row_starts.append(i)
    row_starts.append(len(y_sorted))
    y_bins = np.array([
        float(np.mean(y_sorted[row_starts[k]:row_starts[k + 1]]))
        for k in range(len(row_starts) - 1)
    ])
    sel_idx = np.where(sel)[0][order]

    acc = np.zeros(len(y_bins))
    n_frames = 0
    for _, path in frames:
        mesh = read_vtu(path)
        vg = np.asarray(mesh.cell_data["Gas_Velocity"]).astype(float)
        vy = vg[:, 1] if vg.ndim == 2 else vg
        vy_sorted = vy[sel_idx]
        for k in range(len(y_bins)):
            acc[k] += float(np.mean(
                vy_sorted[row_starts[k]:row_starts[k + 1]]
            ))
        n_frames += 1
    v_mean = acc / n_frames
    t_window = (float(frames[0][0]), float(frames[-1][0]))

    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    ax.plot(v_mean, y_bins, lw=1.8, color="C0", label="V_3d simulation (time-avg)")
    digitized_label = None
    if digitized_csv is not None and digitized_csv.exists():
        df = pd.read_csv(digitized_csv, comment="#")
        digitized_label = f"S1 Fig. 1 (digitized, placeholder)"
        ax.plot(df["v_g_m_s"], df["y_m"], "o--", color="C3", lw=1.2, ms=4,
                label=digitized_label)
    if u_in_m_s:
        ax.axvline(u_in_m_s, color="k", lw=0.5, ls=":", alpha=0.5,
                   label=f"orifice jet U_in = {u_in_m_s:.0f} m/s")
    ax.axhline(0.0, color="k", lw=0.3, alpha=0.4)
    ax.set_xlabel("centerline axial gas velocity v_g [m/s]")
    ax.set_ylabel("height above orifice y [m]")
    ax.set_title(
        f"V_3d (Day 10, ORNL/UTK): centerline v_g vs S1 Fig. 1\n"
        f"(time-avg window: t in [{t_window[0]:.2f}, {t_window[1]:.2f}] s, "
        f"{n_frames} frames; r_select = {r_select_m*1000:.1f} mm)",
        fontsize=10,
    )
    ax.legend(loc="best", fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    # Also save the numeric trace for the record.
    pd.DataFrame({"y_m": y_bins, "v_g_m_s": v_mean}).to_csv(
        out_path.with_suffix(".csv"), index=False,
    )
    return {
        "n_frames": n_frames,
        "t_window_s": t_window,
        "r_select_m": r_select_m,
    }


def plot_inlet_psd(run_dir: Path, out_path: Path, t_startup_s: float = 1.0,
                   monitor_dt_s: float = 1.0e-3) -> dict:
    p_in = read_monitor(run_dir / "V_3D_P_INLET.csv").dropna()
    mask = p_in.index.to_numpy() >= t_startup_s
    sig = p_in.loc[mask, "p_g"].to_numpy() - P_ATM
    t_window_s = float(p_in.index[-1] - t_startup_s)
    fs_hz = 1.0 / monitor_dt_s
    res = welch_psd(sig, fs_hz=fs_hz)

    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    ax.semilogy(res.freqs_hz, res.psd, lw=1.0, color="C0")
    ax.axvspan(17.0, 25.0, color="C3", alpha=0.12, label="S2 band 17-25 Hz")
    ax.axvline(10.5, color="C2", lw=1.0, ls="--", label="S1 ~10-11 Hz")
    ax.axvline(res.dominant_hz, color="k", lw=1.0, ls=":",
               label=f"peak = {res.dominant_hz:.1f} Hz")
    ax.set_xlim(0.0, 100.0)
    ax.set_xlabel("frequency [Hz]")
    ax.set_ylabel(r"PSD of $P_{inlet} - P_{atm}$ [Pa$^2$/Hz]")
    ax.set_title(
        f"V_3d (Day 10, ORNL/UTK): inlet-pressure PSD\n"
        f"(window: t in [{t_startup_s:.2f}, {p_in.index[-1]:.2f}] s, "
        f"{t_window_s:.2f} s span)",
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


def plot_dp_and_inventory(run_dir: Path, dp_out: Path, inv_out: Path) -> dict:
    p_in = read_monitor(run_dir / "V_3D_P_INLET.csv")
    p_top = read_monitor(run_dir / "V_3D_P_TOP_BED.csv")
    df = pd.DataFrame({
        "P_inlet_gauge_Pa": p_in["p_g"] - P_ATM,
        "P_top_gauge_Pa": p_top["p_g"] - P_ATM,
        "dP_bed_Pa": p_top["p_g"] - p_in["p_g"],
    }).dropna()
    fig, (ax_p, ax_dp) = plt.subplots(2, 1, figsize=(7.5, 5.4), sharex=True)
    ax_p.plot(df.index, df["P_inlet_gauge_Pa"], lw=1.0, color="C0",
              label=r"$P_{inlet}-P_{atm}$")
    ax_p.plot(df.index, df["P_top_gauge_Pa"], lw=1.0, color="C3",
              label=r"$P_{top\,bed}-P_{atm}$")
    ax_p.axhline(0.0, color="k", lw=0.5, ls="--", alpha=0.4)
    ax_p.set_ylabel("gauge pressure [Pa]")
    ax_p.set_title("Case V_3d (Day 10, ORNL/UTK): bed pressures at U = 30 m/s")
    ax_p.legend(loc="best", fontsize=9)
    ax_p.grid(True, alpha=0.3)
    ax_dp.plot(df.index, df["dP_bed_Pa"], lw=1.0, color="C2")
    ax_dp.axhline(0.0, color="k", lw=0.5, ls="--", alpha=0.4)
    ax_dp.set_xlabel("time [s]")
    ax_dp.set_ylabel(r"$\Delta P_{bed}$ [Pa]")
    ax_dp.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(dp_out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    df.to_csv(dp_out.with_suffix(".csv"))

    inv = read_monitor(run_dir / "V_3D_SOLIDS_INVENTORY.csv")
    col = next(c for c in inv.columns if c.startswith("ep_s"))
    inv_df = inv.rename(columns={col: "ep_s_vol_m3"}).dropna()
    inv_df["frac_retained"] = inv_df["ep_s_vol_m3"] / inv_df["ep_s_vol_m3"].iloc[0]
    fig, ax = plt.subplots(figsize=(7.5, 3.5))
    ax.plot(inv_df.index, inv_df["frac_retained"], lw=1.0, color="C4")
    ax.axhline(0.95, color="r", lw=0.5, ls="--", label="95% retention floor")
    ax.set_xlabel("time [s]")
    ax.set_ylabel("solids inventory / inventory(0) [-]")
    ax.set_title("V_3d (Day 10, ORNL/UTK): solids inventory retention")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()
    fig.savefig(inv_out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    inv_df.to_csv(inv_out.with_suffix(".csv"))
    return {
        "dP_min": float(df["dP_bed_Pa"].min()),
        "dP_max": float(df["dP_bed_Pa"].max()),
        "retention_final": float(inv_df["frac_retained"].iloc[-1]),
    }


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


def _frame_vgas_max_from_cells(mesh) -> float:
    """Max |V_gas| per frame without using PyVista's cell_centers()."""
    vg = np.asarray(mesh.cell_data["Gas_Velocity"]).astype(float)
    if vg.ndim == 2:
        return float(np.sqrt((vg * vg).sum(axis=1)).max())
    return float(np.abs(vg).max())


def check_max_cfl(run_dir: Path, pvd_path: Path, dx_min_m: float,
                   out_csv: Path, run_name: str) -> pd.DataFrame:
    txt_path = run_dir / f"{run_name.upper()}.TXT"
    dt_log = _parse_dt_from_txt(txt_path) if txt_path.exists() else None
    rows = []
    for t_frame, vtu in read_pvd(pvd_path):
        mesh = read_vtu(vtu)
        v_max = _frame_vgas_max_from_cells(mesh)
        if dt_log is not None and len(dt_log):
            m = dt_log[(dt_log["t"] > t_frame - 0.01) & (dt_log["t"] < t_frame + 0.01)]
            dt_local = float(m["dt"].mean()) if len(m) else float(dt_log["dt"].min())
        else:
            dt_local = 1e-4
        cfl = cfl_number(u_max=max(v_max, 1e-12), dx=dx_min_m, dt=dt_local)
        rows.append((t_frame, v_max, dt_local, cfl))
    df = pd.DataFrame(rows, columns=["t_s", "vgas_max_m_s", "dt_s", "CFL"])
    df.to_csv(out_csv, index=False)
    return df


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run-dir", type=Path, default=DEFAULT_RUN_DIR)
    p.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    p.add_argument("--digitized", type=Path, default=DEFAULT_DIGITIZED)
    p.add_argument("--t-startup", type=float, default=1.0,
                   help="skip the startup window [s] for time-avg and PSD")
    args = p.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    pvd = args.run_dir / "BACKGROUND.pvd"
    if not pvd.exists():
        raise FileNotFoundError(pvd)
    manifest = _load_manifest(args.run_dir)
    dx_min_m = float(manifest["dx_min_m"])
    # r_select for the centerline column: R_i (orifice radius) + 1.5 x
    # near-axis cell size.  Day-10 R1 uses a uniform 2 mm mesh, so the
    # near-axis cells on either side of the axis have centroids at r up
    # to ~1.4 mm; r_select = R_i + 1.5*dx_min_m = 2 + 3 = 5 mm captures
    # the central column.  MFiX cut-cell VTUs report centroids shifted
    # into the fluid portion of each cell, so the nominal-centerline
    # cells have centroid r > R_i.  R_i = 0.002 m (params/geometry.yaml
    # bed_V).
    r_select_m = 0.002 + 1.5 * dx_min_m

    centerline_png = args.out_dir / "fig_d10_V_3d_centerline_vg.png"
    psd_png = args.out_dir / "fig_d10_V_3d_psd.png"
    dp_png = args.out_dir / "fig_d10_V_3d_dP.png"
    inv_png = args.out_dir / "fig_d10_V_3d_inventory.png"
    cfl_csv = args.out_dir / "fig_d10_V_3d_cfl.csv"

    centerline_info = plot_centerline_vg(
        pvd, centerline_png, args.digitized,
        r_select_m=r_select_m, t_startup_s=args.t_startup,
        u_in_m_s=manifest.get("U_in_m_s"),
    )
    psd_info = plot_inlet_psd(args.run_dir, psd_png,
                              t_startup_s=args.t_startup,
                              monitor_dt_s=manifest.get("monitor_dt_s", 1e-3))
    dp_inv_info = plot_dp_and_inventory(args.run_dir, dp_png, inv_png)
    cfl_df = check_max_cfl(args.run_dir, pvd, dx_min_m, cfl_csv,
                            run_name=manifest["run_name"])
    print(f"wrote {centerline_png}  "
          f"({centerline_info['n_frames']} frames, "
          f"window {centerline_info['t_window_s']})")
    print(f"wrote {psd_png}  dominant frequency = {psd_info['dominant_hz']:.2f} Hz "
          f"(window {psd_info['t_window_s']:.2f} s)")
    print(f"wrote {dp_png}  dP range "
          f"{dp_inv_info['dP_min']:.1f} .. {dp_inv_info['dP_max']:.1f} Pa")
    print(f"wrote {inv_png}  final retention = {dp_inv_info['retention_final']:.3f}")
    print(f"wrote {cfl_csv}  max |V_gas| = {cfl_df['vgas_max_m_s'].max():.2f} m/s; "
          f"max CFL = {cfl_df['CFL'].max():.3f}")


if __name__ == "__main__":
    main()

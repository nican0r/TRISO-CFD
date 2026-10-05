"""Day 6-7 fluidization-curve post-processing.

For each fb_sweep run directory, read the VTK output, compute ΔP(t) across the
bed (same recipe as src/post/make_tut_tfm2d_plots.py), time-average over the
last ``t_avg`` seconds, and assemble a CSV + the day's figure.

Entry points:
    time_averaged_dp(data_dir, t_avg_s) -> (dp_mean_pa, dp_series)
    collect_sweep(sweep_root, out_csv, t_avg_s) -> pandas.DataFrame
    plot_fluidization_curve(csv, out_png, ...) -> None
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from src.post.read_vtk import cell_field, load_pvd

# ---------- ΔP(t) time series ------------------------------------------------


def _dp_from_frame(mesh) -> float:
    """ΔP = mean(P_g bottom row) - mean(P_g top row) [Pa] for a 2D structured grid."""
    pg = cell_field(mesh, "P_G", "p_g")
    centers = mesh.cell_centers().points
    yc = centers[:, 1]
    ys = np.unique(np.round(yc, 10))
    y_bot = ys[0]
    y_top = ys[-1]
    tol = 1e-9
    bot = np.isclose(yc, y_bot, atol=tol)
    top = np.isclose(yc, y_top, atol=tol)
    return float(np.mean(pg[bot]) - np.mean(pg[top]))


def _pvd_in(data_dir: Path) -> Path:
    cands = sorted(Path(data_dir).glob("*.pvd"))
    if not cands:
        raise FileNotFoundError(f"No .pvd found in {data_dir}")
    # Prefer the Background collection written by VTK output 1.
    for p in cands:
        if "background" in p.stem.lower():
            return p
    return cands[0]


def time_averaged_dp(data_dir: Path, t_avg_s: float) -> tuple[float, pd.Series]:
    """Return (mean ΔP [Pa] over the last ``t_avg_s`` seconds, ΔP(t) series)."""
    pvd = _pvd_in(Path(data_dir))
    frames = load_pvd(pvd)
    if not frames:
        raise RuntimeError(f"No frames in {pvd}")
    times = np.array([f.time for f in frames])
    dps = np.array([_dp_from_frame(f.mesh) for f in frames])
    series = pd.Series(dps, index=pd.Index(times, name="t [s]"), name="dp [Pa]")
    t_end = float(times[-1])
    mask = times >= (t_end - float(t_avg_s))
    if mask.sum() < 2:
        # Not enough frames in the window; fall back to the last-half window
        mask = times >= 0.5 * t_end
    return float(dps[mask].mean()), series


# ---------- sweep collector --------------------------------------------------


def collect_sweep(
    sweep_root: Path, out_csv: Path, t_avg_s: float
) -> pd.DataFrame:
    """Walk every ``u*/`` directory under ``sweep_root`` and tabulate ΔP̄.

    Each run directory is expected to carry a ``manifest.json`` written by the
    runner (U [m/s], ratio [-], mesh tag).  Runs with no .pvd yet are skipped
    (status 'pending'); runs that error out are recorded with ΔP̄=NaN.
    """
    rows = []
    for run_dir in sorted(Path(sweep_root).iterdir()):
        if not run_dir.is_dir():
            continue
        manifest = run_dir / "manifest.json"
        if not manifest.exists():
            continue
        meta = json.loads(manifest.read_text())
        try:
            dp_mean, series = time_averaged_dp(run_dir, t_avg_s)
            status = "ok"
            n_frames = len(series)
            t_end = float(series.index[-1])
        except FileNotFoundError:
            dp_mean, status, n_frames, t_end = (float("nan"), "pending", 0, 0.0)
        except Exception as e:  # pragma: no cover - runtime failure path
            dp_mean, status, n_frames, t_end = (float("nan"), f"error:{e}", 0, 0.0)
        rows.append(
            dict(
                tag=meta["tag"],
                ratio=meta["ratio"],
                U_m_s=meta["U"],
                mesh=meta["mesh"],
                imax=meta["imax"],
                jmax=meta["jmax"],
                dp_mean_Pa=dp_mean,
                n_frames=n_frames,
                t_end_s=t_end,
                status=status,
            )
        )
    df = pd.DataFrame(rows).sort_values(["mesh", "ratio"]).reset_index(drop=True)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    return df


# ---------- figure -----------------------------------------------------------


def _ergun_curve(u_arr, d_p, rho_g, mu_g, eps, H):
    from src.correlations import ergun_dp
    return np.array([ergun_dp(u, d_p, rho_g, mu_g, eps) * H for u in u_arr])


def _intersection_u(u_arr, dp_curve, dp_plateau):
    """Smallest U at which ``dp_curve`` crosses ``dp_plateau`` from below."""
    diff = dp_curve - dp_plateau
    for i in range(1, len(u_arr)):
        if diff[i - 1] <= 0.0 <= diff[i]:
            u0, u1 = u_arr[i - 1], u_arr[i]
            d0, d1 = diff[i - 1], diff[i]
            if d1 == d0:
                return 0.5 * (u0 + u1)
            return float(u0 - d0 * (u1 - u0) / (d1 - d0))
    return float("nan")


def plot_fluidization_curve(
    csv: Path,
    out_png: Path,
    *,
    umf_wen_yu_m_s: float,
    bed_weight_pa: float,
    d_p_m: float,
    rho_g_kg_m3: float,
    mu_g_pa_s: float,
    eps_packed: float,
    H_static_m: float,
    t_avg_s: float,
) -> dict:
    """Produce fig_d7_fluidization_curve.png and return a summary dict.

    The simulated ΔP̄ vs U is plotted alongside:
      - Ergun ΔP = (dP/L)_Ergun * H_static  (packed-bed branch, U <= U_mf)
      - bed weight per area = (ρ_p − ρ_g) * ε_s * g * H_static  (plateau branch)
    U_mf,sim is the Ergun / plateau crossover drawn from the simulated points
    (nearest linear interpolation).  The returned dict carries this plus the
    plateau ΔP̄ (mean of points above U_mf,wen-yu), used as the Done-when
    acceptance read-out.
    """
    import matplotlib.pyplot as plt

    df = pd.read_csv(csv)
    base = df[df["mesh"] == "base"].sort_values("U_m_s")

    u_dense = np.linspace(
        max(1e-4, float(base["U_m_s"].min()) * 0.5),
        float(base["U_m_s"].max()) * 1.2,
        300,
    )
    ergun_pa = _ergun_curve(
        u_dense, d_p_m, rho_g_kg_m3, mu_g_pa_s, eps_packed, H_static_m
    )
    plateau = np.full_like(u_dense, bed_weight_pa)

    # Intersection of simulated points with the plateau line.
    u_sim = base["U_m_s"].to_numpy()
    dp_sim = base["dp_mean_Pa"].to_numpy()
    good = np.isfinite(dp_sim)
    umf_sim = _intersection_u(u_sim[good], dp_sim[good], bed_weight_pa)

    fig, ax = plt.subplots(figsize=(6.0, 4.5))
    ax.plot(u_dense, ergun_pa, "--", color="tab:gray", label="Ergun (packed)")
    ax.plot(u_dense, plateau, ":", color="k", label=f"bed weight/area = {bed_weight_pa:.0f} Pa")
    ax.axvline(umf_wen_yu_m_s, color="tab:red", lw=1.0, label=f"U_mf,Wen-Yu = {umf_wen_yu_m_s:.3f} m/s")
    if np.isfinite(umf_sim):
        ax.axvline(umf_sim, color="tab:blue", lw=1.0, ls="-.",
                   label=f"U_mf,sim ≈ {umf_sim:.3f} m/s")
    ax.plot(base["U_m_s"], base["dp_mean_Pa"], "o", color="tab:blue", label="sim (base mesh)")
    ax.set_xlabel("Inlet superficial velocity U [m/s]")
    ax.set_ylabel("Time-averaged ΔP across bed [Pa]")
    ax.set_title(
        "Day 7 fluidization curve — Case D (500 µm, 6000 kg/m³) in air @ 20 °C\n"
        f"ΔP averaged over last {t_avg_s:.1f} s of each run"
    )
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="best")
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)

    # Plateau check: mean ΔP̄ over U >= U_mf,Wen-Yu
    sup_mask = good & (u_sim >= umf_wen_yu_m_s)
    plateau_mean = float(np.mean(dp_sim[sup_mask])) if sup_mask.any() else float("nan")

    return dict(
        umf_sim_m_s=umf_sim,
        umf_wen_yu_m_s=umf_wen_yu_m_s,
        bed_weight_pa=bed_weight_pa,
        plateau_sim_pa=plateau_mean,
        plateau_rel_err=(
            (plateau_mean - bed_weight_pa) / bed_weight_pa
            if np.isfinite(plateau_mean)
            else float("nan")
        ),
    )


def _cli() -> None:
    from src.correlations import umf_wen_yu
    from src.gasprops import gas_props
    from src.make_case import REPO_ROOT, _load_yaml

    p = argparse.ArgumentParser(description="Collect fb_sweep ΔP and plot the fluidization curve.")
    p.add_argument("--sweep-root", type=Path, default=REPO_ROOT / "cases" / "fb_sweep")
    p.add_argument("--csv-out", type=Path, default=REPO_ROOT / "results" / "fb_sweep_dp.csv")
    p.add_argument("--fig-out", type=Path, default=REPO_ROOT / "results" / "fig_d7_fluidization_curve.png")
    args = p.parse_args()

    particles = _load_yaml(REPO_ROOT / "params" / "particles.yaml")
    geometry = _load_yaml(REPO_ROOT / "params" / "geometry.yaml")
    gases = _load_yaml(REPO_ROOT / "params" / "gases.yaml")
    case_D = particles["case_D"]
    fb = particles["fb_sweep"]
    gas_case = gases["cases"][fb["gas_case"]]
    rho_g, mu_g = gas_props(gas_case["T_K"], gases["pressure_Pa"], gas_case["composition"])

    umf = umf_wen_yu(case_D["d_p_m"], case_D["rho_p_kg_per_m3"], rho_g, mu_g)
    g = 9.80665
    bed_weight = (case_D["rho_p_kg_per_m3"] - rho_g) * float(fb["eps_s_bed"]) * g * float(fb["H_static_m"])
    eps_packed = 1.0 - float(fb["eps_s_bed"])
    t_avg = float(fb["t_avg_s"])

    df = collect_sweep(args.sweep_root, args.csv_out, t_avg)
    print(df.to_string(index=False))

    if df["dp_mean_Pa"].notna().any():
        summary = plot_fluidization_curve(
            args.csv_out,
            args.fig_out,
            umf_wen_yu_m_s=umf,
            bed_weight_pa=bed_weight,
            d_p_m=case_D["d_p_m"],
            rho_g_kg_m3=rho_g,
            mu_g_pa_s=mu_g,
            eps_packed=eps_packed,
            H_static_m=float(fb["H_static_m"]),
            t_avg_s=t_avg,
        )
        print(json.dumps(summary, indent=2))
    else:
        print("No completed runs yet; figure not produced.")


if __name__ == "__main__":
    _cli()

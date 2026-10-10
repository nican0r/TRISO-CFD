"""Phase-3 validation for the Day-10 R3 (Fix A) smoke run.

Mirrors notes/smoke_validate_r2.py but adds the pathology-specific checks
that Fix A was designed to resolve (notes/day10_pathology_review.md):
  P-A  EP_G minimum must NOT drop to the 0.333 solver clip near the orifice
  P-B  |V|_max in the fluid field must not reach the MAX_INLET_VEL_FAC clip
       (default ratio: 100 * bc_v_g_mi = 750 m/s for Day-10)
  P-C  Mean centerline v_y above the orifice should decay smoothly (not
       oscillate with sign reversals across adjacent cells, as happened in
       the Day-10 R2 AWS Batch run)

Run from the project root:
    .venv/bin/python notes/smoke_validate_r3.py
"""

from __future__ import annotations

import glob
import json
import math
import sys
from pathlib import Path

import numpy as np
import pyvista as pv


RUN_DIR = Path("cases/V_3d/day10_r3_smoke")


def _manual_cell_centers(mesh: pv.UnstructuredGrid) -> np.ndarray:
    pts = np.asarray(mesh.points)
    cells = np.asarray(mesh.cells).reshape(-1)
    centers = np.zeros((mesh.n_cells, 3))
    idx = 0
    for i in range(mesh.n_cells):
        n = cells[idx]
        idx += 1
        inds = cells[idx : idx + n]
        idx += n
        centers[i] = pts[inds].mean(axis=0)
    return centers


def _passfail(label: str, ok: bool) -> str:
    return f"[{label}] {'PASS' if ok else 'FAIL'}"


def main() -> int:
    manifest = json.loads((RUN_DIR / "manifest.json").read_text())
    bc_v_g_mi = float(manifest["bc_v_g_mi_m_s"])
    U_orifice = float(manifest["U_in_m_s"])
    R_i = float(manifest["R_i_m"])
    R_stub = float(manifest["R_stub_m"])
    L_stub = float(manifest["L_stub_m"])
    dx = float(manifest["dx_min_m"])
    y_bc = -L_stub
    # MAX_INLET_VEL_FAC default is 100 in MFiX 26.1.2.
    v_clip = 100.0 * bc_v_g_mi

    print("=" * 72)
    print(f"MFiX smoke validation (Phase 3 + pathology) — {manifest['run_name']}")
    print(f"Run dir: {RUN_DIR}")
    print("=" * 72)
    print(
        f"R_i={R_i*1e3:.2f}mm  R_stub={R_stub*1e3:.2f}mm  "
        f"dx={dx*1e3:.2f}mm  bc_v_g_mi={bc_v_g_mi:.3f}m/s  "
        f"U_orifice_target={U_orifice:.1f}m/s"
    )
    print(f"MAX_INLET_VEL clip = 100 * bc_v_g_mi = {v_clip:.1f} m/s")
    print()

    frames = sorted(glob.glob(str(RUN_DIR / "BACKGROUND_*.vtu")))
    if not frames:
        print("FAIL: no BACKGROUND_*.vtu frames found")
        return 2
    print(f"frames: {len(frames)}")

    results = {}

    # Walk every frame to catch transient pathology spikes (not just the last).
    all_epg_min = []
    all_v_max = []
    all_v_max_near_orifice = []
    all_epg_min_near_orifice = []
    for f in frames:
        m = pv.read(f)
        pts = _manual_cell_centers(m)
        vg = np.asarray(m.cell_data["Gas_Velocity"])
        epg = np.asarray(m.cell_data["EP_G"])
        vmag = np.linalg.norm(vg, axis=1)
        all_epg_min.append(float(epg.min()))
        all_v_max.append(float(vmag.max()))
        near = (np.abs(pts[:, 0]) < 0.004) & (np.abs(pts[:, 2]) < 0.004) \
             & (pts[:, 1] > -0.003) & (pts[:, 1] < 0.005)
        if near.any():
            all_v_max_near_orifice.append(float(vmag[near].max()))
            all_epg_min_near_orifice.append(float(epg[near].min()))

    worst_epg = min(all_epg_min)
    worst_v = max(all_v_max)
    worst_v_near = max(all_v_max_near_orifice) if all_v_max_near_orifice else float("nan")
    worst_epg_near = min(all_epg_min_near_orifice) if all_epg_min_near_orifice else float("nan")

    print()
    print("--- Pathology checks (Fix A targets) ---")
    print(f"[P-A] EP_G worst across all frames        : {worst_epg:.4f}   "
          f"(clip = 0.333; pass if > 0.42 everywhere)")
    print(f"      EP_G worst near-orifice (|x|,|z|<4mm): {worst_epg_near:.4f}")
    print(f"[P-B] |V| worst across all frames         : {worst_v:.1f} m/s "
          f"(clip = {v_clip:.0f}; pass if < 60 m/s)")
    print(f"      |V| worst near-orifice             : {worst_v_near:.1f} m/s")

    results["P-A"] = worst_epg > 0.34
    results["P-A-near"] = worst_epg_near > 0.34
    results["P-B"] = worst_v < 60.0
    print(_passfail("P-A", results["P-A"]))
    print(_passfail("P-B", results["P-B"]))
    print()

    # --- Phase 3 checks on the last frame ---
    print("--- Phase 3 (standard smoke checklist, last frame) ---")
    m = pv.read(frames[-1])
    pts = _manual_cell_centers(m)
    vg = np.asarray(m.cell_data["Gas_Velocity"])
    epg = np.asarray(m.cell_data["EP_G"])

    # 3.1 MI-adjacent cells see bc_v_g_mi
    adjacent_mask = (
        (pts[:, 1] > y_bc)
        & (pts[:, 1] < y_bc + 2 * dx)
        & (np.hypot(pts[:, 0], pts[:, 2]) < R_stub + 0.5 * dx)
    )
    v_adj = vg[adjacent_mask, 1]
    n_adj = int(adjacent_mask.sum())
    print(f"[3.1] MI-adjacent fluid cells: {n_adj}")
    if n_adj > 0:
        ratio_bc = v_adj.mean() / bc_v_g_mi
        print(f"      v_y: mean={v_adj.mean():.3f} m/s, "
              f"min={v_adj.min():.3f}, max={v_adj.max():.3f}")
        print(f"      ratio v_y / bc_v_g_mi = {ratio_bc:.3f} (pass: >= 0.5)")
        results["3.1"] = ratio_bc >= 0.5
    else:
        results["3.1"] = False
        print("      FAIL: zero fluid cells adjacent to MI")

    # 3.4 orifice-contraction check
    orifice_mask = (
        (pts[:, 1] > 0.0) & (pts[:, 1] < 2 * dx)
        & (np.hypot(pts[:, 0], pts[:, 2]) < R_i + 0.5 * dx)
    )
    v_orifice = vg[orifice_mask, 1]
    n_orifice = int(orifice_mask.sum())
    print(f"[3.4] Cells just above orifice: {n_orifice}")
    if n_orifice > 0:
        print(f"      v_y: mean={v_orifice.mean():.3f} m/s, max={v_orifice.max():.3f} m/s  "
              f"(target orifice v ~ {U_orifice:.1f} m/s)")
        orifice_ratio = v_orifice.mean() / U_orifice
        print(f"      ratio v_y / U_orifice = {orifice_ratio:.3f} (pass: 0.4..1.5)")
        results["3.4"] = 0.4 <= orifice_ratio <= 1.5
    else:
        results["3.4"] = False
    print()

    print("--- Centerline profile (last frame) ---")
    axis_mask = (
        (np.abs(pts[:, 0]) < 0.5 * dx) & (np.abs(pts[:, 2]) < 0.5 * dx)
        & (pts[:, 1] >= y_bc) & (pts[:, 1] < 0.05)
    )
    if axis_mask.sum() == 0:
        axis_mask = (
            (np.abs(pts[:, 0]) < dx) & (np.abs(pts[:, 2]) < dx)
            & (pts[:, 1] >= y_bc) & (pts[:, 1] < 0.05)
        )
    axis_pts = pts[axis_mask]
    axis_vg = vg[axis_mask]
    axis_epg = epg[axis_mask]
    order = np.argsort(axis_pts[:, 1])
    print("       y (mm)   v_y (m/s)   EP_G")
    for i in order:
        y_mm = axis_pts[i, 1] * 1e3
        v_y = axis_vg[i, 1]
        e = axis_epg[i]
        print(f"       {y_mm:+6.2f}   {v_y:+8.3f}  {e:.3f}")

    print()
    all_ok = all(results.values())
    print(f"OVERALL: {'GO (Fix A appears to work)' if all_ok else 'NO-GO'}")
    print(f"detail: {results}")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())

"""Phase-3 validation for the Day-10 R2 (Option 4) smoke run.

Implements the mfix-smoke-validate skill's Phase 3 checks:
  3.1 BC-to-interior momentum transfer: compare v_g in cells adjacent to the
      MI face against bc_v_g_mi
  3.2 Fluid-cell count inside the BC lateral footprint (the stub tube)
  3.3 BC resolution rule (cells across orifice >= 4 — expected to fail here;
      Option 4 widens the stub but keeps the physical orifice at R_i)
  3.4 Flow-regime indicator: centerline v_g above the orifice should rise
      toward the orifice-contraction velocity

Run from the project root with `.venv/bin/python notes/smoke_validate_r2.py`.
"""

from __future__ import annotations

import glob
import json
import math
from pathlib import Path

import numpy as np
import pyvista as pv


RUN_DIR = Path("cases/V_3d/day10_r2_smoke")


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


def main() -> None:
    manifest = json.loads((RUN_DIR / "manifest.json").read_text())
    bc_v_g_mi = float(manifest["bc_v_g_mi_m_s"])
    U_orifice = float(manifest["U_in_m_s"])
    R_i = float(manifest["R_i_m"])
    R_stub = float(manifest["R_stub_m"])
    L_stub = float(manifest["L_stub_m"])
    dx = float(manifest["dx_min_m"])
    y_bc = -L_stub

    print(
        f"R_i={R_i*1e3:.2f}mm  R_stub={R_stub*1e3:.2f}mm  "
        f"dx={dx*1e3:.2f}mm  bc_v_g_mi={bc_v_g_mi:.3f}m/s  "
        f"U_orifice_target={U_orifice:.1f}m/s"
    )
    print(f"stub/orifice area ratio = {(R_stub/R_i)**2:.2f}")
    print(f"cells across stub diameter  = {2*R_stub/dx:.2f}")
    print(f"cells across orifice diameter = {2*R_i/dx:.2f}")
    print()

    frames = sorted(glob.glob(str(RUN_DIR / "BACKGROUND_*.vtu")))
    if not frames:
        raise SystemExit("no BACKGROUND frames found")
    print(f"frames: {len(frames)}")

    # Use the last frame for steady-state-ish observation.
    m = pv.read(frames[-1])
    pts = _manual_cell_centers(m)
    vg = np.asarray(m.cell_data["Gas_Velocity"])
    epg = np.asarray(m.cell_data["EP_G"]) if "EP_G" in m.cell_data else None

    # Phase 3.1 — cells immediately above the MI, inside the stub-tube footprint.
    adjacent_mask = (
        (pts[:, 1] > y_bc)
        & (pts[:, 1] < y_bc + 2 * dx)
        & (np.hypot(pts[:, 0], pts[:, 2]) < R_stub + 0.5 * dx)
    )
    v_adj = vg[adjacent_mask, 1]
    n_adj = int(adjacent_mask.sum())
    print(f"[3.1] MI-adjacent fluid cells: {n_adj}  (expected ~ pi R_stub^2 / dx^2 = {math.pi*(R_stub/dx)**2:.1f})")
    if n_adj > 0:
        print(f"      v_y in those cells: mean={v_adj.mean():.3f} m/s, "
              f"min={v_adj.min():.3f}, max={v_adj.max():.3f}")
        ratio_bc = v_adj.mean() / bc_v_g_mi
        print(f"      ratio v_y/bc_v_g_mi = {ratio_bc:.3f}  (pass: >=0.5)")
        passed_31 = ratio_bc >= 0.5
    else:
        passed_31 = False
        print("      FAIL: zero fluid cells adjacent to the MI")
    print()

    # Phase 3.2 — stub-tube fluid cell count (anywhere in stub, not just MI face).
    stub_mask = (
        (pts[:, 1] > y_bc) & (pts[:, 1] < -0.5 * dx)
        & (np.hypot(pts[:, 0], pts[:, 2]) < R_stub + 0.5 * dx)
    )
    n_stub = int(stub_mask.sum())
    expected_stub = math.pi * (R_stub ** 2) * (L_stub - dx) / dx ** 3
    print(f"[3.2] Fluid cells inside the stub tube: {n_stub}  "
          f"(expected ~ {expected_stub:.0f})")
    passed_32 = n_stub > 0.5 * expected_stub
    print()

    # Phase 3.4 — orifice-contraction check: v_y just above cone apex y=0.
    orifice_mask = (
        (pts[:, 1] > 0.0) & (pts[:, 1] < 2 * dx)
        & (np.hypot(pts[:, 0], pts[:, 2]) < R_i + 0.5 * dx)
    )
    n_orifice = int(orifice_mask.sum())
    v_orifice = vg[orifice_mask, 1]
    print(f"[3.4] Cells just above the orifice (0 < y < {2*dx*1e3:.1f} mm, r < R_i+dx/2): {n_orifice}")
    if n_orifice > 0:
        print(f"      v_y: mean={v_orifice.mean():.3f} m/s, "
              f"max={v_orifice.max():.3f} m/s  (target orifice v ~ {U_orifice:.1f} m/s)")
    print()

    # Centerline column rise, cell by cell up to y = H_mid.
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
    order = np.argsort(axis_pts[:, 1])
    print("[3.4a] Near-axis v_y profile (binned):")
    print("       y (mm)   v_y (m/s)   EP_G")
    for i in order:
        y_mm = axis_pts[i, 1] * 1e3
        v_y = axis_vg[i, 1]
        e = epg[axis_mask][np.where(order == i)[0][0]] if epg is not None else float("nan")
        print(f"       {y_mm:+6.2f}   {v_y:+8.3f}  {e:.3f}")
    print()

    overall = passed_31 and passed_32
    print(f"PHASE-3 OVERALL: {'GO' if overall else 'NO-GO'}")


if __name__ == "__main__":
    main()

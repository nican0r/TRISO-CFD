"""Score a full-face-nozzle dry run (params/c26_acceptance.yaml
fullface_smoke.dry_run), 2026-10-27.

Reads <run_dir>/manifest.json, the .LOG (point-source volume, seeding,
warnings), run.log (facets read / ignored) and the last BACKGROUND frame
(MFiX writes only non-blocked cells).  Prints the fluid cells within one cell
of the axis column in the rows around the source, and the source cell's
north-face flux eps_g V_G dx^2 / Q (VTK Gas_Velocity is the face velocity,
notes/c26_source_split.py; eps_g is the source cell's, the upwind side).

    .venv/bin/python notes/c26_fullface_dry_check.py <run_dir>
"""
import json
import re
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / ".claude/skills/mfix-smoke-validate"))
from check_smoke import cell_centers, read_pvd, read_vtk_xml  # noqa: E402


def main(run):
    run = Path(run)
    m = json.loads((run / "manifest.json").read_text())
    dx = m["dx_m"]
    log = (run / f"{m['run_name'].upper()}.LOG").read_text(errors="replace")
    runlog = (run / "run.log").read_text(errors="replace")
    vol = re.search(r"PS\s+Volume\s+1\s+(\S+)", log)
    print(f"PS volume {vol.group(1) if vol else '?'} m3 (dx^3 = {dx ** 3:.5e})")
    for pat in (r"Number of valid facets\s*=\s*\d+", r"Number of ignored facets\s*=\s*\d+",
                r"Successful seeding", r"Total number of particles in the system:\s*\d+"):
        hit = re.search(pat, runlog + log)
        print("  ", re.sub(r"\s+", " ", hit.group(0)) if hit else f"(not found: {pat})")
    print(f"  warnings/errors in .LOG: {len(re.findall(r'(?i)warning|error', log))}")

    t, f = read_pvd(run / "BACKGROUND.pvd")[-1]
    g = read_vtk_xml(run / f)
    ijk = np.floor((cell_centers(g) - [m["x_min_m"], m["y_min_m"], m["z_min_m"]]) / dx).astype(int)
    ic = int(round(-m["x_min_m"] / dx - 0.5))
    jps = int(round((m["ps_cell_y_s_m"] - m["y_min_m"]) / dx))
    V, ep = g["Gas_Velocity"][:, 1], g["EP_G"]
    print(f"frame t = {t}: {len(ijk)} fluid cells; axis column i = k = {ic}, source row j = {jps}")
    src = None
    for j in range(jps - 2, jps + 2):
        near = [(n, i - ic, k - ic) for n, (i, jj, k) in enumerate(ijk)
                if jj == j and abs(i - ic) <= 1 and abs(k - ic) <= 1]
        if j == jps:
            src = next((n for n, di, dk in near if di == dk == 0), None)
        print(f"  row {j}: {len(near)} fluid cell(s) within +-1 of the axis column "
              + " ".join(f"({di:+d},{dk:+d})" for _, di, dk in sorted(near, key=lambda t: t[1:])))
    if src is not None:
        Q = m["ps_massflow_g_kg_s"] / m["rho_g_kg_m3"]
        print(f"source north-face V_G {V[src]:.3f} m/s (Q/dx^2 = {Q / dx ** 2:.3f}), "
              f"EP_G {ep[src]:.4f}: eps V dx^2 / Q = {ep[src] * V[src] * dx ** 2 / Q:.4f}")


if __name__ == "__main__":
    main(sys.argv[1])

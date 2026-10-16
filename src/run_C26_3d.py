"""Generate Case C26 run directories (Conry et al. 2026, Table 2 matrix).

For each selected case in params/sweep_C26.yaml, writes:

  cases/C26_3d/<id>/C26_<id>.mfx      rendered deck (cases/C26_3d/template_dem.mfx.j2)
  cases/C26_3d/<id>/geometry_0001.stl closed cone + column STL (gitignored)
  cases/C26_3d/<id>/manifest.json     run provenance
  cases/C26_3d/<id>/run.sh            MFiX runner (serial or DMP)

and cases/C26_3d/sweep_index.json listing every generated case.  It does
NOT launch MFiX.  Per CLAUDE.md rules 11 and 13, smoke-run the base case
(c01) locally and validate it with the mfix-smoke-validate skill before
running anything else.

Examples:

  .venv/bin/python -m src.run_C26_3d --case c01 --tstop 0.3 --run-dir-suffix _smoke
  .venv/bin/python -m src.run_C26_3d --case all
  .venv/bin/python -m src.run_C26_3d --case c01 --dry-run
"""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

from src.make_case import REPO_ROOT, render_case
from src.make_case_c26 import build_C26_3d_params, load_sweep
from src.make_stl import build_cone_column_stl_implicit
from src.run_V_3d_day10 import _parse_nodes, write_runner

TEMPLATE = REPO_ROOT / "cases" / "C26_3d" / "template_dem.mfx.j2"
OUT_ROOT = REPO_ROOT / "cases" / "C26_3d"

# Keys too bulky or redundant for the manifest.
_MANIFEST_SKIP = {"dx_m_list", "dy_m_list", "dz_m_list"}


def _select(cases: list[dict], spec: str) -> list[dict]:
    if spec == "all":
        return list(cases)
    wanted = [s.strip() for s in spec.split(",") if s.strip()]
    by_id = {c["id"]: c for c in cases}
    missing = [w for w in wanted if w not in by_id]
    if missing:
        raise SystemExit(f"unknown case id(s) {missing}; known: {sorted(by_id)}")
    return [by_id[w] for w in wanted]


def _build_stl(params: dict, path: Path) -> int:
    build_cone_column_stl_implicit(
        R_c_m=params["R_c_m"],
        R_t_m=params["R_t_m"],
        cone_half_angle_deg=params["cone_half_angle_deg"],
        H_dom_m=params["H_dom_m"],
        dx_sample_m=params["dx_sample_m"],
        out_path=path,
        top_extend_m=params["top_extend_m"],
    )
    with path.open("r") as f:
        return sum(1 for ln in f if "facet normal" in ln)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--case", default="all",
                   help="case id(s) from params/sweep_C26.yaml, comma-separated, or 'all'")
    p.add_argument("--tstop", type=float, default=None,
                   help="override simulated end time [s] (default: sweep base.tstop_s)")
    p.add_argument("--nodes", type=_parse_nodes, default=(1, 1, 1), metavar="NI,NJ,NK",
                   help="DMP decomposition; must divide imax/jmax/kmax (27/104/27)")
    p.add_argument("--out-root", type=Path, default=OUT_ROOT)
    p.add_argument("--run-dir-suffix", default="",
                   help="appended to each run dir name, e.g. '_smoke'")
    p.add_argument("--no-stl", action="store_true",
                   help="render decks and manifests only (skip the STL build)")
    p.add_argument("--dry-run", action="store_true",
                   help="print manifests; write nothing")
    args = p.parse_args()

    sweep = load_sweep()
    selected = _select(sweep["cases"], args.case)
    index = []
    stl_cache: Path | None = None
    facet_count = None

    for case in selected:
        name_suffix = args.run_dir_suffix
        params = build_C26_3d_params(
            case, tstop_s=args.tstop, run_name=f"C26_{case['id']}{name_suffix}",
            nodes=args.nodes,
        )
        run_dir = args.out_root / f"{case['id']}{name_suffix}"
        mfx_path = run_dir / f"{params['run_name']}.mfx"
        manifest = {k: v for k, v in params.items() if k not in _MANIFEST_SKIP}
        manifest.update({
            "case": "C26 (Conry et al. 2026, Nucl. Eng. Des. 456:115009)",
            "template": str(TEMPLATE.relative_to(REPO_ROOT)),
            "mfx": mfx_path.name,
            "stl": params["stl_filename"],
            "is_base_case": case["id"] == sweep["base"]["base_case"],
            "piv_targets_table3": sweep.get("piv_targets_table3", {}).get(case["id"]),
            "notes": "notes/c26_conry2026_case.md",
        })
        if args.dry_run:
            print(json.dumps(manifest, indent=2))
            continue

        run_dir.mkdir(parents=True, exist_ok=True)
        stl_path = run_dir / params["stl_filename"]
        if not args.no_stl:
            # Geometry is identical for every case: build once, then copy.
            if stl_cache is None:
                facet_count = _build_stl(params, stl_path)
                stl_cache = stl_path
            else:
                shutil.copy2(stl_cache, stl_path)
            manifest["facet_count"] = facet_count
        render_case(TEMPLATE, params, mfx_path)
        (run_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        write_runner(run_dir, mfx_path.name, params["nranks"])
        index.append({
            "id": case["id"], "run_dir": str(run_dir.relative_to(REPO_ROOT)),
            "feedstock": params["feedstock"], "flow_slpm": params["flow_slpm"],
            "solids_model": params["solids_model"], "cgp_stat_wt": params["cgp_stat_wt"],
            "mew_w": params["mew_w"], "U_nozzle_m_s": params["U_nozzle_m_s"],
            "dt_max_s": params["dt_max_s"], "n_parcels_est": params["n_parcels_est"],
            "tstop_s": params["tstop_s"],
        })
        print(f"wrote {run_dir}")

    if not args.dry_run and index and not args.run_dir_suffix:
        # Merge into the existing index so a single-case regeneration does
        # not drop the other cases.
        index_path = args.out_root / "sweep_index.json"
        merged = {}
        if index_path.is_file():
            merged = {row["id"]: row for row in json.loads(index_path.read_text())}
        merged.update({row["id"]: row for row in index})
        rows = [merged[k] for k in sorted(merged)]
        index_path.write_text(json.dumps(rows, indent=2) + "\n")


if __name__ == "__main__":
    main()

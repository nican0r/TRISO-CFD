"""Generate Case C26 run directories (flow / coating-stage / charge sweep).

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
  # non-uniform partition: narrow slabs over the cone, where the parcels are
  .venv/bin/python -m src.run_C26_3d --case c01 --gridmap-i 9,3,3,3,9 --gridmap-k 9,3,3,3,9

With neither --nodes nor --gridmap-*, cases get DEFAULT_GRIDMAP: 4 ranks
(2x1x2, i = k = [13,14]), the laptop split and the best throughput per vCPU on
AWS (notes/c26_conry2026_case.md, "AWS run cost"). The manifest's aws_vcpu
(= nranks) and aws_memory_mb (src/aws/config.py C26_MEMORY_MB) size the Batch
job; src/aws/submit.py reads them when --vcpu/--memory are omitted.
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
from pathlib import Path

from src.make_case import REPO_ROOT, render_case
from src.make_case_c26 import build_C26_3d_params, load_sweep, write_gridmap
from src.make_stl import build_cone_column_stl_implicit
from src.run_V_3d_day10 import _parse_nodes, write_runner
from src.aws.config import C26_MEMORY_MB

TEMPLATE = REPO_ROOT / "cases" / "C26_3d" / "template_dem.mfx.j2"
OUT_ROOT = REPO_ROOT / "cases" / "C26_3d"

# Default DMP partition (cells per rank along i, j, k): 2x1x2 = 4 ranks.
# 27 cells do not split evenly in two, hence a gridmap rather than --nodes.
DEFAULT_GRIDMAP = ([13, 14], [104], [13, 14])

# Parcel count of the c01 memory probe behind C26_MEMORY_MB (src/aws/config.py,
# results/run_log.csv 2026-10-19 "c01 memory probe").
C26_MEMORY_PROBE_PARCELS = 27782


def aws_memory_mb(n_parcels_est: float) -> int:
    """Batch memory [MB]: C26_MEMORY_MB scaled linearly with the parcel count
    above the probed c01 deck, rounded up to 1000 MB.  Conservative (memory
    was 0.6-1.1 GB at 27,782 parcels); revisit with each case's smoke RSS."""
    scale = max(1.0, float(n_parcels_est) / C26_MEMORY_PROBE_PARCELS)
    return int(math.ceil(C26_MEMORY_MB * scale / 1000.0) * 1000)


# Keys too bulky or redundant for the manifest.
_MANIFEST_SKIP = {"dx_m_list", "dy_m_list", "dz_m_list"}


def _parse_sizes(spec: str) -> list[int]:
    try:
        sizes = [int(s) for s in spec.split(",") if s.strip()]
    except ValueError as e:
        raise argparse.ArgumentTypeError(f"gridmap sizes must be ints: {e}")
    if not sizes or min(sizes) < 1:
        raise argparse.ArgumentTypeError(f"gridmap sizes must be >= 1; got {spec!r}")
    return sizes


def _repo_relative(path: Path) -> str:
    """Path relative to the repo root, or absolute for an --out-root outside it."""
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


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
    p.add_argument("--nodes", type=_parse_nodes, default=None, metavar="NI,NJ,NK",
                   help="DMP decomposition; must divide imax/jmax/kmax (27/104/27). "
                        "Default (no --nodes, no --gridmap-*): DEFAULT_GRIDMAP, 4 ranks; "
                        "1,1,1 for serial")
    for axis in "ijk":
        p.add_argument(f"--gridmap-{axis}", type=_parse_sizes, default=None, metavar="N,N,...",
                       help=f"non-uniform DMP partition: cells per rank along {axis} "
                            "(writes gridmap.dat; overrides --nodes; an omitted axis is "
                            "one rank)")
    p.add_argument("--stl-sample-fraction", type=float, default=None, metavar="F",
                   help="STL SDF sampling step as a fraction of dx (default: "
                        "geometry.yaml dx_sample_fraction); facet count ~ 1/F^2")
    ps = p.add_mutually_exclusive_group()
    ps.add_argument("--ps-mass-only", action="store_true",
                    help="point source adds gas mass only (ps_v_g = 0); the "
                         "geometry.yaml ps_momentum default since 2026-10-19")
    ps.add_argument("--ps-momentum", action="store_true",
                    help="point source carries the nozzle momentum "
                         "(ps_v_g = U_nozzle), for A/Bs against the default")
    p.add_argument("--epg-clip", type=float, default=None, metavar="EPG",
                   help="DES_EPG_CLIP override (gas-fraction floor); default "
                        "geometry.yaml des_epg_clip")
    p.add_argument("--interp-scheme", choices=("LINEAR_HAT", "SQUARE_DPVM"), default=None,
                   help="DES_INTERP_SCHEME override; SQUARE_DPVM also writes "
                        "DES_INTERP_WIDTH (geometry.yaml des_interp_width_cells * dx)")
    p.add_argument("--flow-slpm", type=float, default=None, metavar="Q",
                   help="override the case's gas flow [slpm] for a one-off check "
                        "outside the sweep matrix (e.g. c07 spout check); "
                        "requires --run-dir-suffix")
    p.add_argument("--out-root", type=Path, default=OUT_ROOT)
    p.add_argument("--run-dir-suffix", default="",
                   help="appended to each run dir name, e.g. '_smoke'")
    p.add_argument("--no-stl", action="store_true",
                   help="render decks and manifests only (skip the STL build)")
    p.add_argument("--dry-run", action="store_true",
                   help="print manifests; write nothing")
    args = p.parse_args()
    if args.flow_slpm is not None and not args.run_dir_suffix:
        p.error("--flow-slpm needs --run-dir-suffix (keeps the sweep's run dir and index)")

    sweep = load_sweep()
    selected = _select(sweep["cases"], args.case)
    gridmap = None
    if any(getattr(args, f"gridmap_{a}") for a in "ijk"):
        gridmap = tuple(getattr(args, f"gridmap_{a}") or [] for a in "ijk")
    elif args.nodes is None:
        gridmap = DEFAULT_GRIDMAP
    nodes = args.nodes or (1, 1, 1)
    index = []
    stl_cache: Path | None = None
    facet_count = None

    for case in selected:
        name_suffix = args.run_dir_suffix
        if args.flow_slpm is not None:
            case = {**case, "flow_slpm": args.flow_slpm}
        params = build_C26_3d_params(
            case, tstop_s=args.tstop, run_name=f"C26_{case['id']}{name_suffix}",
            nodes=nodes, gridmap=gridmap,
            dx_sample_fraction=args.stl_sample_fraction,
            ps_momentum=True if args.ps_momentum else (False if args.ps_mass_only else None),
            epg_clip=args.epg_clip,
            interp_scheme=args.interp_scheme,
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
            "is_ofat_baseline": case["id"] == sweep["base"]["ofat_baseline"],
            "flow_slpm_override": args.flow_slpm is not None,
            "piv_targets_table3": sweep.get("piv_targets_table3", {}).get(case["id"]),
            "notes": "notes/c26_conry2026_case.md",
            # AWS Batch sizing, read by src/aws/submit.py: one vCPU per rank.
            "aws_vcpu": params["nranks"],
            "aws_memory_mb": aws_memory_mb(params["n_parcels_est"]),
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
        gridmap_path = run_dir / "gridmap.dat"
        if params["gridmap"] is not None:
            write_gridmap(run_dir, params["gridmap"])
        elif gridmap_path.exists():
            # MFiX reads any gridmap.dat in the run dir; never leave a stale one.
            gridmap_path.unlink()
        (run_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        write_runner(run_dir, mfx_path.name, params["nranks"])
        index.append({
            "id": case["id"], "run_dir": _repo_relative(run_dir),
            "feedstock": params["feedstock"], "coating_stage": params["coating_stage"],
            "d_p_m": params["d_p_m"], "rho_p_kg_m3": params["rho_p_kg_m3"],
            "inventory_kg": params["inventory_kg"], "H_static_m": params["H_static_m"],
            "flow_slpm": params["flow_slpm"],
            "solids_model": params["solids_model"], "cgp_stat_wt": params["cgp_stat_wt"],
            "d_parcel_m": params["d_parcel_m"], "U_nozzle_m_s": params["U_nozzle_m_s"],
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

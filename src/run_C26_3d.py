"""Generate Case C26 run directories (flow / coating-stage / charge sweep).

For each selected case in params/sweep_C26.yaml, writes:

  cases/C26_3d/<id>/C26_<id>.mfx      rendered deck (cases/C26_3d/template_dem.mfx.j2)
  cases/C26_3d/<id>/geometry_0001.stl closed cone + column STL (gitignored)
  cases/C26_3d/<id>/manifest.json     run provenance
  cases/C26_3d/<id>/run.sh            MFiX runner (serial or DMP; for an MFiX
                                      release with no local DMP install, i.e.
                                      the 23.2 default, it runs the AWS image
                                      locally with docker)

It does NOT launch MFiX.  Per CLAUDE.md rules 11 and 13, smoke-run a case
locally and validate it with notes/smoke_validate_c26.py before submitting
its long run.

Decks are for MFiX 23.2 by default (user decision 2026-10-24: all new runs;
numerics.yaml mfix_version); --mfix-version 26.1.2 re-renders the earlier
runs' release.  Use a fresh --run-dir-suffix for new runs: the unsuffixed
c01-c03 dirs and the existing *_smoke dirs are run provenance
(results/run_log.csv), and a run dir holding a case for another release is
never overwritten.

Examples:

  .venv/bin/python -m src.run_C26_3d --case c04 --run-dir-suffix _massonly
  .venv/bin/python -m src.run_C26_3d --case c02 --tstop 0.3 --run-dir-suffix _smoke_r2
  .venv/bin/python -m src.run_C26_3d --case all --dry-run
  # C26's grid layout (axis on a cell face, 2x2-cell point source)
  .venv/bin/python -m src.run_C26_3d --case c01 --axis-on-face --run-dir-suffix _massonly_axisface
  # non-uniform partition: narrow slabs over the cone, where the parcels are
  .venv/bin/python -m src.run_C26_3d --case c01 --gridmap-i 9,3,3,3,9 --gridmap-k 9,3,3,3,9 \
      --run-dir-suffix _gridmap_test
  # an MFiX 26.1.2 deck, as for the runs before 2026-10-24 (local native solver)
  .venv/bin/python -m src.run_C26_3d --case c01 --mfix-version 26.1.2 --run-dir-suffix _v26

With neither --nodes nor --gridmap-*, cases get default_gridmap(): 4 ranks
(2x1x2, i = k = [13,14] on the 27-cell box, [13,13] with --axis-on-face), the
laptop split and the best throughput per vCPU on AWS
(notes/c26_conry2026_case.md, "AWS run cost"). The manifest's aws_vcpu
(= nranks) and aws_memory_mb (src/aws/config.py C26_MEMORY_MB) size the Batch
job; src/aws/submit.py reads them when --vcpu/--memory are omitted.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
from pathlib import Path

from src.make_case import REPO_ROOT, render_case
from src.make_case_c26 import (
    MFIX_VERSIONS, STL_BUILDERS, build_C26_3d_params, default_gridmap, load_sweep,
    write_gridmap,
)
from src.make_stl import build_cone_column_stl_implicit
from src.make_stl_c26_apex import build_cone_column_stl_apex
from src.make_stl_c26_fullface import build_cone_column_stl_fullface
from src.run_V_3d_day10 import _parse_nodes, write_runner
from src.aws.config import C26_MEMORY_MB, PROJECT
from src.aws.submit import manifest_mfix_version

TEMPLATE = REPO_ROOT / "cases" / "C26_3d" / "template_dem.mfx.j2"
OUT_ROOT = REPO_ROOT / "cases" / "C26_3d"

# The release behind the laptop's mamba env (src/run_V_3d_day10.py
# MFIX_SOLVER_*); any other release runs from its Docker image.
LOCAL_MFIX_VERSION = "26.1.2"

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


def _select(cases: list[dict], spec: str) -> list[dict]:
    if spec == "all":
        return list(cases)
    wanted = [s.strip() for s in spec.split(",") if s.strip()]
    by_id = {c["id"]: c for c in cases}
    missing = [w for w in wanted if w not in by_id]
    if missing:
        raise SystemExit(f"unknown case id(s) {missing}; known: {sorted(by_id)}")
    return [by_id[w] for w in wanted]


def _md5(path: Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()


def check_run_dir_release(run_dir: Path, mfix_version: str) -> None:
    """Refuse to render into a run dir that holds a case for another MFiX
    release (manifest ``mfix_version``; none = 26.1.2): those dirs are run
    provenance, and the 23.2 default must not overwrite the 26.1.2 ones."""
    if not (run_dir / "manifest.json").is_file():
        return
    old = manifest_mfix_version(run_dir)
    if old != mfix_version:
        raise SystemExit(f"{run_dir} holds an MFiX {old} case (run provenance); render "
                         f"the MFiX {mfix_version} deck under a new --run-dir-suffix")


def write_docker_runner(run_dir: Path, mfx_name: str, params: dict) -> Path:
    """run.sh for a release with no local DMP install: runs the AWS image
    ``triso-cfd/mfix:<mfix_version>`` (docker/) on the run dir in the
    entrypoint's LOCAL_SMOKE mode, with the deck's own DMP decomposition.
    The entrypoint writes run.log itself (tee), so stdout is discarded."""
    image = f"{PROJECT}/mfix:{params['mfix_version']}"
    script = run_dir / "run.sh"
    lines = [
        "#!/usr/bin/env bash",
        "# Auto-generated by src/run_C26_3d.py",
        f"# MFiX {params['mfix_version']} (no local DMP install): image {image}, "
        "built by src.aws.provision --build",
        f"# nranks = {params['nranks']}",
        "set -euo pipefail",
        f'cd "{run_dir.resolve()}"',
        f'docker run --rm --platform=linux/amd64 -e LOCAL_SMOKE=1 -e MFX_FILE="{mfx_name}" '
        f'-e NP={params["nranks"]} -e NODESI={params["nodes_i"]} '
        f'-e NODESJ={params["nodes_j"]} -e NODESK={params["nodes_k"]} '
        f'-v "$PWD":/work {image} > /dev/null',
        "",
    ]
    script.write_text("\n".join(lines))
    script.chmod(0o755)
    return script


# Facts of the apex builder carried into the manifest (src/make_stl_c26_apex.py).
_APEX_FACTS = ("n_theta", "sagitta_Rc_m", "sagitta_Rt_m", "n_bands_cone", "n_bands_column",
               "open_area_y0_m2", "tip_volume_m3")
# ... and of the full-face nozzle builder (src/make_stl_c26_fullface.py).
_FULLFACE_FACTS = ("n_theta", "sagitta_Rc_m", "n_bands_cone", "n_bands_column",
                   "nozzle_face_area_m2", "bore_area_m2", "sliver_fraction",
                   "protrusion_fraction", "neighbour_fluid_fraction", "y_curve_min_m",
                   "y_curve_max_m", "y_collar_top_m")


def _build_stl(params: dict, path: Path) -> dict:
    """Write the vessel STL for ``params['stl_builder']``; return manifest facts
    (facet_count, stl_md5, and the apex builder's geometry numbers)."""
    builder = params["stl_builder"]
    if builder == "implicit":
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
            facts = {"facet_count": sum(1 for ln in f if "facet normal" in ln)}
    elif builder == "parametric_apex":
        _, apex = build_cone_column_stl_apex(
            R_c_m=params["R_c_m"],
            R_t_m=params["R_t_m"],
            cone_half_angle_deg=params["cone_half_angle_deg"],
            H_dom_m=params["H_dom_m"],
            out_path=path,
            top_extend_m=params["top_extend_m"],
            sagitta_max_m=params["stl_apex_sagitta_max_m"],
            band_height_m=params["stl_apex_band_height_m"],
        )
        if abs(apex["y_apex_m"] - params["y_apex_m"]) > 1e-12:
            raise AssertionError("apex height differs between builder and params")
        facts = {"facet_count": apex["facet_count"]}
        facts.update({f"stl_{k}": apex[k] for k in _APEX_FACTS})
    elif builder == "parametric_fullface":
        _, ff = build_cone_column_stl_fullface(
            R_c_m=params["R_c_m"],
            R_t_m=params["R_t_m"],
            cone_half_angle_deg=params["cone_half_angle_deg"],
            H_dom_m=params["H_dom_m"],
            out_path=path,
            top_extend_m=params["top_extend_m"],
            sagitta_max_m=params["stl_apex_sagitta_max_m"],
            band_height_m=params["stl_apex_band_height_m"],
            cell_half_m=0.5 * params["dx_m"],
            nozzle_y_s_m=params["ps_cell_y_s_m"],
            wall_offset_m=params["stl_fullface_wall_offset_m"],
        )
        if (abs(ff["y_floor_m"] - params["y_nozzle_floor_m"]) > 1e-12
                or abs(ff["half_width_m"] - params["nozzle_half_width_m"]) > 1e-12):
            raise AssertionError("nozzle geometry differs between builder and params")
        facts = {"facet_count": ff["facet_count"]}
        facts.update({f"stl_{k}": ff[k] for k in _FULLFACE_FACTS})
    else:
        raise ValueError(f"unknown stl_builder {builder!r}")
    facts["stl_md5"] = _md5(path)
    return facts


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--case", default="all",
                   help="case id(s) from params/sweep_C26.yaml, comma-separated, or 'all'")
    p.add_argument("--tstop", type=float, default=None,
                   help="override simulated end time [s] (default: sweep base.tstop_s)")
    p.add_argument("--nodes", type=_parse_nodes, default=None, metavar="NI,NJ,NK",
                   help="DMP decomposition; must divide imax/jmax/kmax (27/104/27, "
                        "26/104/26 with --axis-on-face). Default (no --nodes, no "
                        "--gridmap-*): default_gridmap(), 4 ranks; 1,1,1 for serial")
    for axis in "ijk":
        p.add_argument(f"--gridmap-{axis}", type=_parse_sizes, default=None, metavar="N,N,...",
                       help=f"non-uniform DMP partition: cells per rank along {axis} "
                            "(writes gridmap.dat; overrides --nodes; an omitted axis is "
                            "one rank)")
    p.add_argument("--stl-builder", choices=STL_BUILDERS, default=None,
                   help="vessel STL builder (default geometry.yaml stl_builder = implicit). "
                        "parametric_apex: analytic cone closed at its apex y = -R_t/tan30 "
                        "(src/make_stl_c26_apex.py) with the point source in the tip row "
                        "y in [-dx, 0]; the deck is otherwise unchanged. "
                        "parametric_fullface: the apex cone with a one-cell square bore "
                        "cut into its tip around the axis cell of the odd 25-cell box "
                        "(src/make_stl_c26_fullface.py), the "
                        "source open only through its whole top face, ps_v_g = Q/dx^2 "
                        "with --ps-momentum; not with --axis-on-face")
    p.add_argument("--stl-sample-fraction", type=float, default=None, metavar="F",
                   help="implicit builder: STL SDF sampling step as a fraction of dx "
                        "(default: geometry.yaml dx_sample_fraction); facet count ~ 1/F^2")
    ps = p.add_mutually_exclusive_group()
    ps.add_argument("--ps-mass-only", action="store_true",
                    help="point source adds gas mass only (ps_v_g = 0); the "
                         "geometry.yaml ps_momentum default since 2026-10-19")
    ps.add_argument("--ps-momentum", action="store_true",
                    help="point source carries the nozzle momentum "
                         "(ps_v_g = U_nozzle), for A/Bs against the default")
    p.add_argument("--ps-point", action="store_true",
                   help="write the point source as a point on the axis (C26 p.4), "
                        "which MFiX maps to ONE cell (the -x/-z cell of the source "
                        "row) instead of the 2x2 axis cells; needs --axis-on-face "
                        "(geometry.yaml ps_point)")
    p.add_argument("--epg-clip", type=float, default=None, metavar="EPG",
                   help="DES_EPG_CLIP override (gas-fraction floor); default "
                        "geometry.yaml des_epg_clip")
    p.add_argument("--interp-scheme", choices=("LINEAR_HAT", "SQUARE_DPVM", "NONE"), default=None,
                   help="DES_INTERP_SCHEME override; SQUARE_DPVM also writes "
                        "DES_INTERP_WIDTH (geometry.yaml des_interp_width_cells * dx); "
                        "NONE is MFiX's default (no interpolation: des_interp_on and "
                        "des_interp_mean_fields are then written .False.)")
    p.add_argument("--vtk-dt-background", type=float, default=None, metavar="DT",
                   help="cell-field VTK interval [s] (default numerics.yaml vtk_dt_background_s)")
    p.add_argument("--axis-on-face", action="store_true",
                   help="C26's grid layout: an even cell count (26) with the axis on "
                        "the edge of the four central cells, and the point source over "
                        "those 2x2 cells (geometry.yaml axis_on_cell_face); default "
                        "is the odd 27-cell box with a one-cell source on the axis")
    p.add_argument("--flow-slpm", type=float, default=None, metavar="Q",
                   help="override the case's gas flow [slpm] for a one-off check "
                        "outside the sweep matrix (e.g. c07 spout check); "
                        "requires --run-dir-suffix")
    p.add_argument("--mfix-version", choices=MFIX_VERSIONS, default=None,
                   help="MFiX release the deck is for (default numerics.yaml "
                        "mfix_version = 23.2, the AWS image; its run.sh runs "
                        "the image with docker). 26.1.2 re-renders the earlier "
                        "runs' release (the four 26.1.2-only keywords, local "
                        "native run.sh)")
    p.add_argument("--out-root", type=Path, default=OUT_ROOT)
    p.add_argument("--run-dir-suffix", default="",
                   help="appended to each run dir name, e.g. '_smoke'")
    p.add_argument("--no-stl", action="store_true",
                   help="render decks and manifests only (skip the STL build)")
    p.add_argument("--dry-run", action="store_true",
                   help="print manifests; write nothing")
    args = p.parse_args()
    if args.flow_slpm is not None and not args.run_dir_suffix:
        p.error("--flow-slpm needs --run-dir-suffix (keeps the sweep's run dir)")
    fullface = args.stl_builder == "parametric_fullface"
    if fullface and args.axis_on_face:
        p.error("--stl-builder parametric_fullface needs the odd layout; drop --axis-on-face")
    # None = geometry.yaml axis_on_cell_face; the full-face nozzle forces the odd box.
    axis_on_face = True if args.axis_on_face else (False if fullface else None)

    sweep = load_sweep()
    selected = _select(sweep["cases"], args.case)
    gridmap = None
    if any(getattr(args, f"gridmap_{a}") for a in "ijk"):
        gridmap = tuple(getattr(args, f"gridmap_{a}") or [] for a in "ijk")
    elif args.nodes is None:
        gridmap = default_gridmap(axis_on_face)
    nodes = args.nodes or (1, 1, 1)
    stl_cache: Path | None = None
    stl_facts: dict | None = None

    for case in selected:
        name_suffix = args.run_dir_suffix
        if args.flow_slpm is not None:
            case = {**case, "flow_slpm": args.flow_slpm}
        params = build_C26_3d_params(
            case, tstop_s=args.tstop, run_name=f"C26_{case['id']}{name_suffix}",
            nodes=nodes, gridmap=gridmap,
            dx_sample_fraction=args.stl_sample_fraction,
            ps_momentum=True if args.ps_momentum else (False if args.ps_mass_only else None),
            ps_point=True if args.ps_point else None,
            epg_clip=args.epg_clip,
            interp_scheme=args.interp_scheme,
            axis_on_face=axis_on_face,
            vtk_dt_background_s=args.vtk_dt_background,
            stl_builder=args.stl_builder,
            mfix_version=args.mfix_version,
        )
        run_dir = args.out_root / f"{case['id']}{name_suffix}"
        check_run_dir_release(run_dir, params["mfix_version"])
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
                stl_facts = _build_stl(params, stl_path)
                stl_cache = stl_path
            else:
                shutil.copy2(stl_cache, stl_path)
            manifest.update(stl_facts)
        render_case(TEMPLATE, params, mfx_path)
        gridmap_path = run_dir / "gridmap.dat"
        if params["gridmap"] is not None:
            write_gridmap(run_dir, params["gridmap"])
        elif gridmap_path.exists():
            # MFiX reads any gridmap.dat in the run dir; never leave a stale one.
            gridmap_path.unlink()
        (run_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        if params["mfix_version"] == LOCAL_MFIX_VERSION:
            write_runner(run_dir, mfx_path.name, params["nranks"])
        else:
            write_docker_runner(run_dir, mfx_path.name, params)
        print(f"wrote {run_dir}")


if __name__ == "__main__":
    main()

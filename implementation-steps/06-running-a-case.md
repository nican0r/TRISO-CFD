# 6. Running a case

The workflow is: parameters → case directory → local dry run → smoke →
long run on AWS Batch → pull → score → log. The rules behind it are in
[.claude/CLAUDE.md](../.claude/CLAUDE.md): validate the base case first,
dry-run locally before any cloud submission, generous timeouts, and every
run recorded in `results/run_log.csv`.

## 1. Build the case

All physical inputs are in `params/*.yaml`, each with units and a source.
`src/run_C26_3d.py` reads a case row from `params/sweep_C26.yaml`, builds
the parameters (`src/make_case_c26.py`), writes the STL (file 03) and
renders `cases/C26_3d/template_dem.mfx.j2`. Output:

```
cases/C26_3d/<case><suffix>/
  C26_<case><suffix>.mfx   rendered deck (gas properties from Cantera)
  geometry_0001.stl        vessel STL (gitignored; regenerated)
  manifest.json            every parameter, mesh, sizing and provenance
  run.sh                   runner (docker for 23.2)
```

The two primary sweeps were rendered with:

```bash
.venv/bin/python -m pytest tests/
# control inlet: closed apex, 2x2 mass-only tip source, 24-cell grid
.venv/bin/python -m src.run_C26_3d --case c02 --axis-on-face \
    --stl-builder parametric_apex --interp-scheme LINEAR_HAT --run-dir-suffix _apex_v232
# single-cell inlet: full-face bore, one source cell with momentum, 25-cell grid
.venv/bin/python -m src.run_C26_3d --case c02 --stl-builder parametric_fullface \
    --ps-momentum --run-dir-suffix _fullface_bore_mom_nointerp_v232
```

(`--interp-scheme LINEAR_HAT` is needed for the control now that NONE is
the default. `--case all` renders c01–c09.)

- **Use a new `--run-dir-suffix` for every new run.** Committed run dirs
  are provenance, and regenerating one rewrites its deck.
- **`run.sh` changes into the absolute repo case directory.** Running a
  copied `run.sh` runs the production deck in the repo, not the copy.

## 2. Dry run and smoke

1. **Dry run** (1 ms, local docker, 4 ranks, ~30 s): checks that the deck
   parses under 23.2, that the point-source volume is what the layout
   expects (e.g. dx³ for the bore), that the source is walled in, that
   ~27,760 parcels seed for 60 g, and that the STL md5 and gridmap are the
   intended ones. Diff the deck against the reference deck; only the
   intended keywords may differ.
2. **Smoke** (0.3–0.5 s): `src/post/c26/smoke_validate_c26.py <run_dir>`
   wraps the `mfix-smoke-validate` skill's `check_smoke.py` with
   C26-specific checks:
   - point-source volume;
   - inventory: parcels may leave only through the top outlet, ≤ 0.5 % of
     the seeded mass;
   - the jet on the axis;
   - an axial spout channel: upward solids flux ≥ 5 ṁ_g within r < 3.15 mm,
     centroid ≤ 1.05 mm from the axis;
   - no dt collapse.

   The spout checks gate only c01–c03. c04–c09 may legitimately not spout.
   A 0.3 s smoke tests inlet health and numerics only. The bed needs ~1 s to
   settle a fountain, and C26 averages only after 2.4 s.
3. **Pre-register.** Before a run is scored, its pass/fail criteria and
   decision rule are written into `params/c26_acceptance.yaml` (one block
   per test, e.g. `fullface_3s_nointerp`, `drag_closure_smoke`). Results
   are judged against what was written first. Outcomes outside the
   registered branches go to the user.

Known false alarms: the inlet-pressure monitor sits in the source cell on
the closed-apex and bore decks. It reads the jet's acceleration pressure,
so `smoke_validate` reports NO-GO on S6 for every such run. That check is
not meaningful there.

## 3. Submit, poll, pull

```bash
.venv/bin/python -m src.aws.submit --case-dir cases/C26_3d/<run> \
    --mfx C26_<run>.mfx --job-name c26-<run> --dry-run   # drop --dry-run to submit
.venv/bin/python -m src.aws.poll <job-id> --watch --pull cases/C26_3d/<run>/aws_out
```

- `submit` syncs the case to `s3://triso-cfd-mfix/cases/<run>/`, submits
  it to Batch, and appends a `submitted` row to `results/run_log.csv`
  before the job starts. vCPU and memory come from `manifest.json`.
- The Batch stack (queue `triso-cfd-q`, c7a on-demand instances, image
  `triso-cfd/mfix:23.2` in ECR) is described in `src/aws/provision.py`.
- Output lands in `s3://triso-cfd-mfix/results/<run>/`. Raw output is
  never committed (CLAUDE.md rule 7). Superseded runs' local `aws_out/`
  may have been deleted to save disk; re-pull them from S3.
- Credentials come from `aws login`, which the venv's boto3 cannot read,
  so submit and poll shell out to the AWS CLI (`src/aws/cli.py`).
- Jobs have no timeout and never retry. Budget 12–30 h per 3 s case
  (file 05).

## 4. Score

`src/post/c26/c26_watch_post.py` polls a list of jobs, pulls each when it
succeeds, and runs the scoring steps:

- health (`smoke_validate_c26.py`);
- step-3 spout checks (`c26_step3_checks.py`);
- synthetic PIV (`src.post.piv_c26`);
- source-gas split (`c26_source_split.py`).

It then writes `results/c26_compare_<tag>.csv`. With `--local` it scores
run dirs already on disk. Example (single-cell sweep, compared with the
control runs):

```bash
MPLBACKEND=Agg .venv/bin/python src/post/c26/c26_watch_post.py \
    c0{1..9}_fullface_bore_mom_nointerp_v232 --tag fullface_nointerp_v232 --baseline-suffix _apex_v232
MPLBACKEND=Agg .venv/bin/python src/post/c26/c26_sweep_summary.py \
    --suffix _fullface_bore_mom_nointerp_v232 --tag fullface_nointerp_v232
# PIVlab re-score (MATLAB R2026b + PIVlab 3.14), file 07
.venv/bin/python -m src.post.piv_c26_export_frames cases/C26_3d/c02_fullface_bore_mom_nointerp_v232/aws_out
src/post/piv_c26_pivlab_batch.sh c02_fullface_bore_mom_nointerp_v232 pivlab_c26img
.venv/bin/python -m src.post.piv_c26_pivlab_summary --tag _fullface_nointerp_v232 \
    --runs c01_fullface_bore_mom_nointerp_v232 c02_fullface_bore_mom_nointerp_v232 \
           c03_fullface_bore_mom_nointerp_v232_p90 --map-variants pivlab_c26img pivlab_c26img_noclahe
# Fig. 4-style composite and the cone-bottom schematic
MPLBACKEND=Agg .venv/bin/python src/post/c26/c26_fig4_compose.py \
    fullface_bore_mom_nointerp_v232 c01 c02 c03_p90
MPLBACKEND=Agg .venv/bin/python src/post/c26/draw_cone_bottom.py
```

All headless plotting needs `MPLBACKEND=Agg`. The measurement methods are
in file 07.

## 5. Log

Each run gets a row in `results/run_log.csv`: case, key parameters, mesh,
dt, wall-clock time and status. `submit` writes the `submitted` row;
marking a run done is a manual review step after scoring. Edit the log
line-wise and diff it. Rewriting it with a CSV writer can truncate it on a
malformed row.

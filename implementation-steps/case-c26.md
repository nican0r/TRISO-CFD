# Case C26 — Conry et al. 2026 benchtop spout-fluidized bed

Task: reproduce the MFiX CFD-CGDEM model of Conry, Dal Forno Chuahy, Oyedeji,
Finney, Heldt & López-Honorato (2026), *Nucl. Eng. Des.* 456: 115009, on the
project's MFiX 26.1.2, then benchmark it against the paper's PIV data.
Design, parameter sources and recipe deviations: `notes/c26_conry2026_case.md`.

## Already in place (branch `case-c26-conry2026`)

- Params: `bed_C26`, `spouted_bed_C26_3d`, `spouted_bed_C26_3d_dem`
  (`params/geometry.yaml`); `case_C26_ysz`, `case_C26_ss`
  (`params/particles.yaml`); `ar_25C`, `slpm_reference_C26`
  (`params/gases.yaml`); `spouted_bed_C26_3d_numerics`
  (`params/numerics.yaml`); Table 2 matrix `params/sweep_C26.yaml`.
- Builder `src/make_case_c26.py`, STL `src/make_stl.py::build_cone_column_stl_implicit`,
  template `cases/C26_3d/template_dem.mfx.j2`, driver `src/run_C26_3d.py`.
- Tests: `tests/test_make_case_c26.py`, `tests/test_make_stl_c26.py`.

## Steps

1. **Regenerate locally** (the decks embed Cantera gas properties and the
   STL is gitignored): `.venv/bin/python -m pytest tests/` then
   `.venv/bin/python -m src.run_C26_3d --case all`.
2. **Base-case smoke (CLAUDE.md rule 11, 13).**
   `.venv/bin/python -m src.run_C26_3d --case c01 --tstop 0.3 --run-dir-suffix _smoke`,
   run `cases/C26_3d/c01_smoke/run.sh`, validate with the
   `mfix-smoke-validate` skill. Case-specific checks:
   - run.log "Point source setup information" shows a non-zero PS volume
     (no `PS_VOLUME == ZERO` abort); the apex cell is fluid.
   - Seeded parcels ≈ manifest `n_parcels_est` (shortfall ≤ 10 %, recipe Rule 8).
   - Settled bed height after the initial collapse: report it next to
     H_static = 34.8 mm (60° cone) and the paper's ~25 mm annotation.
   - A spout forms at 10 slpm (fountain above the bed by the end of the smoke);
     no dt collapse; CFL < 2; NIT < 30.
   If the apex cell is blocked, set `ps_row_offset_cells: 1` first.
   Do not widen the inlet (recipe Rule 3 history).
3. **Baseline runs c01–c03** (tstop 3 s). Log each in `results/run_log.csv`.
4. **Synthetic PIV post-processor** (`src/post/piv_c26.py`, not written yet):
   project `Particles_PIV` frames onto the camera plane, render them, run the
   same windowed cross-correlation as PIVlab (~4 mm windows, C26 Fig. 1c),
   time-average 300 frames, and compute CC, SSI and MSE against the
   experimental maps. Digitize C26 Fig. 4 maps and wall v(t) into
   `data/digitized/`.
5. **Sensitivity cases c04–c12** after c01–c03 pass.

## Done when

- c01 smoke passes the step-2 checks, and the result is recorded in `results/run_log.csv`.
- c01–c03 reach t = 3.0 s without `DT < DT_MIN`.
- CC, SSI and MSE for c01–c03 are reported next to C26 Table 3
  (0.68/0.82/0.83, 0.61/0.77/0.77, 0.0068/0.0042/0.0049 m²/s²) with the
  averaging window stated (CLAUDE.md rule 6). No tuning to match.
- The settled-bed-height discrepancy (35 vs ~25 mm) is stated plainly, with the
  cone-angle implication.
- All tests pass.

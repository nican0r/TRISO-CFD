# Case C26 — Conry et al. 2026 benchtop spout-fluidized bed

Task: reproduce the MFiX CFD-CGDEM model of Conry, Dal Forno Chuahy, Oyedeji,
Finney, Heldt & López-Honorato (2026), *Nucl. Eng. Des.* 456: 115009, on the
project's MFiX 26.1.2, then benchmark it against the paper's PIV data.
Design, parameter sources and recipe deviations: `notes/c26_conry2026_case.md`.

## Already in place (branch `case-c26-conry2026`)

- Params: `bed_C26`, `spouted_bed_C26_3d`, `spouted_bed_C26_3d_dem`
  (`params/geometry.yaml`); `case_C26_ysz`
  (`params/particles.yaml`); `ar_25C`, `slpm_reference_C26`
  (`params/gases.yaml`); `spouted_bed_C26_3d_numerics`
  (`params/numerics.yaml`); `triso_coating_layers` (`params/particles.yaml`);
  case matrix `params/sweep_C26.yaml` (C26 Table 2 cases 1–3 plus the
  coating-stage / charge sweep c04–c09; notes "Coating-stage / charge sweep").
- Builder `src/make_case_c26.py`, STL `src/make_stl.py::build_cone_column_stl_implicit`,
  template `cases/C26_3d/template_dem.mfx.j2`, driver `src/run_C26_3d.py`.
- Tests: `tests/test_make_case_c26.py`, `tests/test_make_stl_c26.py`.
- Checks: `notes/smoke_validate_c26.py` (step 2), `notes/c26_step3_checks.py`
  (step 3), `notes/c26_sweep_summary.py` (step 5).

## Steps

1. **Generate decks locally** (the decks embed Cantera gas properties and the
   STL is gitignored): `.venv/bin/python -m pytest tests/` then
   `.venv/bin/python -m src.run_C26_3d --case <ids> --run-dir-suffix <new suffix>`.
   Use a new suffix: the committed run dirs (unsuffixed c01–c03, `*_massonly`,
   `*_smoke*`) are run provenance, and regenerating one rewrites its deck.
   The production decks are `cases/C26_3d/<id>_massonly` (2026-10-19).
2. **Base-case smoke (CLAUDE.md rule 11, 13).** Done: c01_smoke_r3 GO
   (2026-10-17); mass-only c01_smoke_r4_massonly (2026-10-19). To repeat:
   `.venv/bin/python -m src.run_C26_3d --case c01 --tstop 0.3 --run-dir-suffix _smoke_r7`,
   run its `run.sh`, validate with `notes/smoke_validate_c26.py` (which runs
   the `mfix-smoke-validate` skill's `check_smoke.py`). The 0.3 s smoke tests **inlet health and
   numerics only**. C26 needs ≥ 3.0 s to pass the startup transient (p.3,
   §2.2.2), so a spout is judged in step 3, not here. Case-specific checks:
   - run.log "Point source setup information" shows the **axis cell's** fluid
     volume (≈ 7.0e-9 m³ for row 0; 9.26e-9 m³ for a full cell). A volume of
     ~3.3e-10 m³ means the source landed in a corner sliver cut cell (first
     c01 smoke, 2026-10-16).
   - Inlet jet on the axis (BACKGROUND frames, t ≥ 0.1 s): the fastest
     first-row gas cell is the axis cell, with v_g within 0.5–2× U_nozzle
     (37.4 m/s); no off-axis cell above 2× U_nozzle.
   - Axial spout channel (all PIV frames, t = 0.1 s to tstop; slab
     y = 25.2 ± 1.05 mm):
     - (i) time-mean upward solids mass flux through r < 3.15 mm ≥ 5 ṁ_g;
     - (ii) the centroid radius of the upward solids mass flux is ≤ 1.05 mm
       (half a cell).
     Both are time means, so the result does not depend on the pulse phase.
     Thresholds were written before evaluation (2026-10-17). First smoke
     (misplaced source): 19 ṁ_g, 5.8 mm → fails (ii). r2: 39 ṁ_g, 0.08 mm.
     Supporting check (threshold chosen after seeing both runs): the time-mean
     of the centerline v_g monitor at mid-bed is ≥ 5 m/s. Uniform cone flow
     would give ~0.9 m/s interstitial; the runs gave 0.5 vs 13.9 m/s.
     This replaces an earlier "axis gas momentum flux at 10 mm ≥ 0.5 ṁ·U"
     test, which is unsound for a jet into a bed: the bed weight is ~58 ṁU,
     and in r2 the gas hands ~0.65 ṁU to spout particles by 11.5 mm.
   - Seeded parcels ≈ manifest `n_parcels_est` (shortfall ≤ 10 %, recipe Rule 8);
     Monitor 3 (particle mass/count) loses parcels only through the top
     pressure outlet (Particle_ID tracking), ≤ 0.5 % of the seeded mass.
     (Amended 2026-10-19 by user decision: was "stays at the seeded value";
     the c02 smoke lost 3 ballistic fountain parcels via the PO.)
   - Point source: mass-only (`ps_v_g = 0`) is the default since 2026-10-19
     (user decision; `--ps-momentum` for A/Bs). For mass-only decks the
     `check_smoke` S3 exit-velocity band and the mid-bed v_g ≥ 5 m/s check
     were set on momentum runs, so they are reported, not gated, and the
     spout is judged on the axial spout-channel checks.
   - Spout checks (axial spout channel, mid-bed v_g) gate only C26's PIV
     cases c01–c03. For the coating-stage / charge cases at fixed 12 slpm,
     no spout is reported as a result (user decision 2026-10-19).
   - Settled bed height after the initial collapse: report it next to
     H_static = 34.8 mm (60° cone) and the C26 Fig. 1b "~25 mm" schematic
     annotation.
   - No dt collapse; CFL < 2; NIT < 30.
   The `check_smoke.py` stub checks (S3–S5) assume Case V's stub pipe. For
   C26, run it with `--skip-stub-parcels --exit-voidage-info` and use the
   inlet-jet checks above instead. (`--exit-voidage-info` added 2026-10-19,
   after the c01 baseline's exit-layer EP_G of 0.891; see SKILL.md §3.6.) Do not widen the inlet
   (recipe Rule 3 history).
3. **Baseline runs c01–c03** (tstop 3 s). Log each in `results/run_log.csv`.
   Status: c01 with the momentum source reached 3 s (health GO; step-3
   A/B PASS, C1/D FAIL); c01–c03_massonly were submitted 2026-10-19.
   Spout check, over the last 0.6 s (C26's PIV window): at 10 slpm C26 calls
   the regime "near minimum spouting" with "approximately 10 spouting events
   every 0.5 s" (§3.1, §2.2.2), so expect a pulsed spout, not a steady one.
   - Fountain parcels above the bed surface in most PIV frames.
   - Pulsation of order 10–20 Hz in ΔP.
   - Vertical particle velocity across the bed at 10 slpm peaking near
     0.5 m/s on the axis and ~0 at the walls (C26 Fig. 7, model).
   - Mean downward wall velocity of order 1e-3 m/s (C26 Fig. 6d: experiment
     0.93e-3, their model 3.3e-3).
   Batch jobs carry no timeout (no `attemptDurationSeconds`); budget
   ~17–21 h per case at 4 ranks from the smoke's s/step (C26 needed 167 h
   for 3 s at W = 10, Table 5).
4. **Synthetic PIV post-processor** (`src/post/piv_c26.py`, not written yet):
   project `Particles_PIV` frames onto the camera plane, render them, run the
   same windowed cross-correlation as PIVlab (~4 mm windows, C26 Fig. 1c),
   time-average 300 frames, and compute CC, SSI and MSE against the
   experimental maps. Digitize C26 Fig. 4 maps and wall v(t) into
   `data/digitized/`.
5. **Coating-stage / charge sweep c04–c09** (one factor at a time around c02:
   YSZ, 60 g, 12 slpm, W = 10). It replaces C26's wall-friction,
   statistical-weight and steel cases.
   - c04–c07 add the TRISO layers buffer, IPyC, SiC and OPyC to the YSZ core.
     The parcels grow to 1.81 mm in 2.1 mm cells.
   - c08 and c09 use 40 and 80 g charges.
   - **Smokes first** (rule 13): c02, then c07 (largest parcel), then c04,
     then c09 (most parcels, memory check). Each is checked with
     `notes/smoke_validate_c26.py`.
   - **Full runs:** 3 s runs only after the smokes pass. Then run
     `notes/c26_sweep_summary.py --suffix _massonly` over the last 0.6 s.
   - Status: smokes done (notes "Local smokes and criteria changes");
     c04–c09_massonly submitted 2026-10-19/20.
   - **Reporting:** there is no experiment for c04–c09, so results are
     relative sensitivities. A bed that doesn't spout at 12 slpm is reported
     as is.

## Done when

- c01 smoke passes the step-2 checks, and the result is recorded in
  `results/run_log.csv`. (Done: c01_smoke_r3, 2026-10-17.)
- c01–c03 reach t = 3.0 s without `DT < DT_MIN`.
- CC, SSI and MSE for c01–c03 are reported next to C26 Table 3
  (0.68/0.82/0.83, 0.61/0.77/0.77, 0.0068/0.0042/0.0049 m²/s²) with the
  averaging window stated (CLAUDE.md rule 6). No tuning to match.
- The settled-bed-height discrepancy (35 vs ~25 mm) is stated plainly, with the
  cone-angle implication. (Done: notes, "Geometry decisions" item 3.)
- All tests pass.

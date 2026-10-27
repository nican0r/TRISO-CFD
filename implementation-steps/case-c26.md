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
- Closed-apex parametric STL (2026-10-24): `src/make_stl_c26_apex.py`
  (analytic surface of revolution, cone continued to its apex at
  y = −R_t/tan 30° = −2.062 mm, 8,480 facets, inward normals checked per
  facet), selected with `run_C26_3d --stl-builder parametric_apex`
  (`geometry.yaml stl_builder`, default `implicit` unchanged); the point
  source then fills the tip row y ∈ [−dx, 0]. Tests
  `tests/test_make_stl_c26_apex.py`; 1 ms dry run c02_apex_dry scored
  against five pre-registered criteria in the notes ("Closed apex,
  parametric STL"): 3 pass, 2 fail on the literal numbers with the cause
  identified (MFiX's polyhedral cut-cell volume of the tip, 1.95e-9 m³ vs
  the 3.06e-9 m³ cone; 8 sliver cells of 1.01 % cell volume under the 1 %
  small-cell cutoff). Not a production default; no AWS submission.
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
   c01–c03_massonly reached 3 s on 2026-10-20 (22.7–23.5 h).
   - Health: GO for c01 and c02. c03 is NO-GO on P3ii because its spout is
     off-centre.
   - Step 3: C1 PASS and D FAIL in all three; A fails for c01 only; C2 fails
     at y = 10 mm in all three.
   - All mass-only runs have the spout leaning towards −x. Gas-only tests
     show the bias exists without a bed and is not the domain split or the
     STL. The likely cause is MFiX 26.1.2's unlimited U interpolation factor
     (`get_alpha.f`); the patched-solver test has not been run.
   - Details: notes, "c01–c03 mass-only 3 s results".
   - C26 grid layout (`--axis-on-face`), c01–c03_massonly_axisface, reached
     3 s on 2026-10-21. A, B, C1 and C2 pass at all three flows. D still
     fails, and changes by −20 / −4 / 0 %. The PIV map gives NO MATCH at all
     three flows. Health is GO for c01 and c02; c03 is NO-GO on S6 (ΔP
     3.09 W/A). Details: notes, "Axis-face 3 s results".
   - First-principles review of the velocity gap (2026-10-22; notes,
     "First-principles review of the velocity gap"). The "missing retarding
     force in the fountain" hypothesis is not supported: drag is ≤ 5 % of
     weight on fountain parcels and descending parcels decelerate relative
     to free fall. The mass-only runs are in a bubbling regime (jet at
     1/3–1/4 of the nozzle momentum flux); the 27-cell momentum c01
     over-drives the spout. The default layout is now exactly 24 cells
     across (`geometry.yaml` margin 0, axis on the cell corner, 2 × 2
     source; no physics change against 26 cells with a blocked margin).
     c01–c09_momentum24 (momentum source on the 24-cell layout) were
     submitted 2026-10-22; c05 and c07 reached 3 s on 2026-10-23 (14–16 h),
     the rest were still running. Not pulled or scored yet.
   - 2026-10-23: c01–c03_momentum24 scored. The momentum source changed
     nothing but the pressure on the 24-cell grid (the four source cells
     are enclosed, so continuity sets the 11–13 m/s jet; the source type
     is closed as a candidate). The PIV convention was changed to C26's
     (`piv_c26.yaml empty_windows: interpolate`, PIVlab-like per-frame
     fill) and every `results/c26_piv_*` output was regenerated: c02 now
     MATCHES C26's model (CC 0.85, SSI 0.84, MSE 0.0037); c03 is close on
     SSI; c01 does not match. D is unchanged (0.016–0.03 m/s, 5–7× C26's
     model): a fluidized, internally circulating annulus. Interpolation
     A/B `c02_momentum24_nointerp` (des_interp_scheme NONE) submitted
     2026-10-23 with pre-registered criteria (`c26_acceptance.yaml`
     interp_ab). Details: notes, "momentum24 c01–c03 results" and
     "Interpolation A/B".
   - 2026-10-23 (later): map bottom moved to the annulus surface
     (`piv_c26.yaml bed_surface_annulus_r_min_m: 15 mm`; C26 masks the
     annular bed surface and the all-parcel percentile was fountain-inflated
     by 1.6–6.3 mm). All `results/c26_piv_*` regenerated again: c01 now at
     the CC floor (match), c02 misses the MSE floor by 0.0001 (CC 0.83, SSI
     0.85), c03 unchanged (no match). A 0.3 s AWS smoke of the no-interp
     deck (`c02_momentum24_nointerp_smoke03`, gas VTK every 5 ms) was
     submitted to read the annulus-fluidization metrics early; criteria in
     `c26_acceptance.yaml interp_ab_smoke`. Details: notes, "Map bottom
     moved to the annulus surface" and "Interpolation A/B, 0.3 s smoke".
   - 2026-10-24: smoke pulled and scored: NO EFFECT (bed EP_G 0.571 vs
     0.571 baseline, void fraction 0.228 vs 0.25; heave persists; jet
     above the bed no stronger). Interpolation is not what fluidizes the
     annulus; the 3 s job (e7395c56…) still to be scored under
     `interp_ab` for the wall D. c08_momentum24 SUCCEEDED, not pulled.
   Spout check, over the last 0.6 s (C26's PIV window): at 10 slpm C26 calls
   the regime "near minimum spouting" with "approximately 10 spouting events
   every 0.5 s" (§3.1, §2.2.2), so expect a pulsed spout, not a steady one.
   - Fountain parcels above the bed surface in most PIV frames.
   - Pulsation of order 10–20 Hz in ΔP.
   - Vertical particle velocity across the bed at 10 slpm peaking near
     0.5 m/s on the axis and ~0 at the walls (C26 Fig. 7, model).
   - Mean downward wall velocity of order 1e-3 m/s (C26 Fig. 6d: experiment
     0.93e-3, their model 3.3e-3).
   Batch jobs carry no timeout (no `attemptDurationSeconds`). Budget
   ~23–24 h per case at 4 ranks, as measured on c01–c03_massonly (the
   momentum24 c05 / c07 took 16.3 / 14.3 h). The
   earlier 17–21 h, taken from a 10 ms window, was too low. (C26 needed
   167 h for 3 s at W = 10, Table 5.)
4. **Synthetic PIV** (`src/post/piv_c26.py`, written 2026-10-20). It
   renders `PARTICLES_PIV` frames as a side-view camera image and runs
   multi-pass FFT cross-correlation on consecutive frames. It writes:
   - the time-averaged map above the bed, scored with CC, SSI and MSE against
     the digitized C26 Fig. 4a;
   - the wall velocity in C26 Fig. 1c's 4 mm box, compared with Fig. 6d.

   Settings and sources are in `params/piv_c26.yaml`. C26 gives no PIV
   settings, so they are assumed. Conry et al. 2025 (CES 306:121197), which
   C26 cites for the process, gives the processing steps but no PIVlab
   settings either; the comparison is in `notes/c26_conry2026_case.md`
   ("Conry et al. 2025 PIV method"). The wall box sits at 0.64 of the bed
   height (C26 Fig. 1c, measured 2026-10-25). The "~4 mm" in Fig. 1c is the wall
   sampling box, not the PIV window (this step's earlier text said
   otherwise). The Fig. 4 maps come from `src/post/digitize_c26_fig4.py`
   (`data/digitized/c26_fig4_*.csv`). Fig. 6 time series are not digitized;
   the Fig. 6d means are given in the text (p.7).

   Run: `MPLBACKEND=Agg .venv/bin/python -m src.post.piv_c26 <run_dir>`
   (about 2 min per run). Report each run next to Table 3 **and** next to
   C26's own model scored through the same pipeline (Fig. 4b vs 4a): the
   digitized figure does not reproduce Table 3 (notes, "Synthetic PIV").
   Status (2026-10-20): c01 (momentum) and c01–c03_massonly are scored
   (notes, "Synthetic PIV"). The map placement and the handling of empty
   windows are assumptions that move CC a lot. For the mass-only runs, SSI
   and MSE stay far from C26's model under every variant tested.
   - 2026-10-27: **re-run with PIVlab 3.14** (MATLAB R2026b) on c01–c03
     `_apex_v232`, C26's settings + PIVlab defaults, scoring in MATLAB
     (`corrcoef`/`ssim`/`immse`): `src/post/piv_c26_pivlab_batch.sh`,
     `results/c26_pivlab_rerun.csv`, `fig_c26_pivlab_{maps,wall}.png`,
     note `notes/c26_pivlab_rerun.md`.  Pure defaults on the shaded 16-bit
     renders are ruined by CLAHE on the black background (fixed spurious
     0.43 m/s vector in every empty window; c02 CC 0.10); with CLAHE off,
     c02 scores CC 0.846 / SSI 0.850 / MSE 0.0028 (C26 model through the
     same code 0.872 / 0.740 / 0.0022), c03 0.633 / 0.771 / 0.0074, c01
     0.464 / 0.706 / 0.0049.  Frames made as C26's Fig. 4d shows (binary
     white parcel disks on black, 8-bit) + untouched PIVlab defaults
     (`pivlab_c26img`): empty windows go NaN and are interpolated as in C26;
     c02 CC 0.778 / SSI 0.825 / MSE 0.0048, c03 0.559 / 0.738 / 0.0085, c01
     0.273 / 0.647 / 0.0070; CLAHE off on the same frames 0.856 / 0.871 /
     0.0025, 0.611 / 0.772 / 0.0078, 0.354 / 0.691 / 0.0058.  Wall PIV on
     binary frames is not usable (merged disks; C26 Fig. 8 describes the
     same failure); the headline wall velocity is the parcels in the box,
     as C26 reports its model (2.09 / 1.53 / 3.04e-2 m/s), which PIVlab on
     the shaded frames reproduces within 2 %.
     Checks of why C26's model map is speckled and ours smooth (notes
     "Why is our map smooth"): their vector spacing measured on Fig. 4d is
     2.0 parcel diameters, same as ours; a finer-grid run (0.05 mm/px,
     64/32/16) stays smooth (CLAHE off: CC 0.870 / SSI 0.874 / MSE 0.0026,
     equal to C26's own model); only a mean over original vectors ("skip")
     reproduces the speckle, but it leaves blanks their maps do not have.
     Headline convention unchanged (per-pair interpolation, then mean).
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
   - Status: smokes done (notes "Local smokes and criteria changes").
     c04–c09_massonly and c04/c05/c06/c08_massonly_axisface reached 3 s
     and were scored on 2026-10-22; c07/c09_massonly_axisface were still
     running. Results are in the notes, "Coating-stage / charge 3 s
     results":
     - The charge block (c08, c09) spouts. Its wall velocity has the same
       gap to C26 as c02.
     - Coated beds stop spouting at 12 slpm from c05 on. c04 spouts weakly.
     - c06 shows an unexplained bed expansion.
     - c05 and c06 are health NO-GO on P2, because losses before the PIV
       window cannot be traced.
   - **Reporting:** there is no experiment for c04–c09, so results are
     relative sensitivities. A bed that doesn't spout at 12 slpm is reported
     as is.

## Done when

- c01 smoke passes the step-2 checks, and the result is recorded in
  `results/run_log.csv`. (Done: c01_smoke_r3, 2026-10-17.)
- c01–c03 reach t = 3.0 s without `DT < DT_MIN`. (Done: c01–c03_massonly,
  2026-10-20; min dt 6.6e-6 s.)
- CC, SSI and MSE for c01–c03 are reported next to C26 Table 3
  (0.68/0.82/0.83, 0.61/0.77/0.77, 0.0068/0.0042/0.0049 m²/s²) and next to
  C26's model scored through the same pipeline, with the averaging window
  stated (CLAUDE.md rule 6). No tuning to match.
- The settled-bed-height discrepancy (35 vs ~25 mm) is stated plainly, with the
  cone-angle implication. (Done: notes, "Geometry decisions" item 3.)
- All tests pass.
- 2026-10-24: closed-apex local 0.05 s smoke passed gate J (jet 2.2–2.5× the 2 × 2 momentum source at 5–10 mm; column asymmetric, one quadrant cell); AWS 0.3 s smoke c26-c02-apex-smoke03 submitted (job 51bb89b1); 3 s decision per params/c26_acceptance.yaml apex_smoke.
- 2026-10-24: AWS 0.3 s apex smoke pulled and scored (notes "AWS smoke result and decision"): gate NOT met (J 1.51× at 5 mm, 1.27× at 10 mm; exit v_p +73 %, PIV map mean +56 %, bed EP_G 0.545); the spout runs up one quadrant cell (a bistable selection, different quadrant locally and on AWS); 3 s run c02_apex generated but NOT launched. Odd-grid (axis-centred) closed apex fails at setup: the single tip cell has no fluid corner node and is blocked.
- 2026-10-24: user decision: 3 s closed-apex run launched on MFiX 23.2 (c02_apex_v232, job 38c9dc92, c26-c02-apex-v232) after a 1 ms 23.2 docker dry run (mesh identical to 26.1.2); predictions pre-registered in params/c26_acceptance.yaml apex_3s; score at 2.4–3.0 s when pulled.
- 2026-10-25: c02_apex_v232 scored at 2.4–3.0 s: all three `apex_3s` predictions met (map mean 0.98× C26 model, CC 0.914, bed EP_G 0.561); first c02 map to MATCH all three C26 floors; wall PIV 0.0198 m/s (3.9× C26 model). Caveat: MFiX 23.2 deletes slow wall parcels at one cut-cell corner (−67 of 27 765; absent under 26.1.2 on the same STL). Next lever: the annulus.
- 2026-10-25: interp_ab (c02_momentum24_nointerp, 3 s) scored: D wall −50 %, annulus metrics unchanged, CC 0.783 → INCONCLUSIVE by the pre-registered rule.
- 2026-10-25/26: source A/B on the closed-apex 23.2 layout (notes "Point-source 3 s results", "Sweep source choice"): control (2×2 mass-only) map vs C26 model CC 0.907 / MSE 0.0015, reach 0.69; one-cell + momentum 0.879 / 0.0058, reach 0.99 but 30 % fast in the lower map and 5 mm lean; 2×2 + momentum 0.849 / 0.0025, reach 0.51. The pre-registered rule (momentum branches only) picked one-cell and its sweep was submitted 10:45 EDT.
- 2026-10-26: user decision: the c01–c09 sweep runs in the control configuration. The eight one-cell jobs were terminated ~11:04 EDT (not pulled); c0X_apex_v232 (c01, c03, c04, c05, c06, c08, c09, c07; MFiX 23.2, 3 s) submitted ~11:08 EDT as c26-c0X-apex-v232. Score at 2.4–3.0 s through the c02_apex_v232 pipeline. Rule record: params/c26_acceptance.yaml sweep_source_choice.decision_2026_10_26.
- 2026-10-24: `src/aws/poll.py --watch` and `notes/c26_watch_post.py` retry a failing describe-jobs call (expired `aws login` session killed the first watcher); `tests/test_aws_poll.py`.
- 2026-10-27: control-configuration sweep pulled and scored (notes "Control-configuration sweep results"; results/c26_compare_apex_v232.csv, c26_sweep_summary_apex_v232.csv, fig_c26_sweep_*_apex_v232.png). All nine at 3.0 s, parcel loss ≤ 0.29 % except c08 0.57 %. Only c02 matches C26 (CC 0.914); c01 CC 0.356 and c03 CC 0.634 / MSE 0.0075 do not (short fountains; both worse than on the open apex). c04–c07 do not spout at 12 slpm. Charge trend physical (fountain share 4.1 / 2.3 / 1.2 %, D 2.8 / 2.1 / 1.8e-2 m/s). ΔP/(W/A) 8–14.6 is the inlet monitor in the source cell on the closed apex; move it one row up in the next decks.
- 2026-10-27: sweep comparison quantities added to `notes/c26_sweep_summary.py` (notes "Sweep comparison quantities"): solids circulation through y = H_static/2 by Particle_ID (c02 26.5 g/s, turnover 2.3 s; c04–c07 14.2 / 5.7 / 3.7 / 0 g/s; charge block flat 25–28 g/s so turnover 1.6 / 2.3 / 2.7 s), U_sup/U_ms Mathur–Gishler (spouting cases ≥ 0.56, non-spouting ≤ 0.54; predicted flows for c02's margin: c04 13.6, c05 15, c06 17, c07 18.4 slpm, to pre-register), bed ΔP from the VTU annulus at 4–8 mm (tracks charge 582 / 804 / 885 Pa where the inlet monitor was flat; ratio to W/A 1.4–2.0 everywhere, not a regime discriminator in the cone). `tests/test_c26_sweep_summary.py` (5 tests).

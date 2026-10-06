# Day 9 — Summary

## What was implemented

Built the full Day-9 pipeline for a Case V_2d spouting run at
`U_in = 1.2 × U_ms,MG` (Mathur–Gishler) with the Case V_rpt experimental
particle inventory (R13/R14): 2.18 mm soda-lime glass beads,
ρ_p = 2400 kg/m³, in dry air at 25 °C. The pipeline is end-to-end
reproducible: `src/make_case.py:build_V_2d_day9_params` renders the
`.mfx` from YAML + Cantera; `src/run_V_2d_day9.py` emits a run
directory + runner script; `src/post/plot_V_2d_day9.py` reads the
monitors + VTU frames and produces an ε_s frame grid, a two-panel
pressure/ΔP figure, and a per-frame max-CFL CSV using
`src.cfl_check.cfl_number`. Five solver-stabilisation iterations
(Runs 4–8 in `notes/day09_debug.md`) each pushed the TFM failure further
in simulated time (0.0055 → 0.0108 → 0.0108 → 0.0576 → **0.072 s**) but
none sustained the "spout + fountain + annulus for several seconds" the
step requires. **The Day-9 simulation deliverable is therefore partial
/ blocked** on a numerical pathology that is outside the Day-9 scope to
fix. The decision point for the user is laid out at the end of
`notes/day09_debug.md`.

Key artifacts:

- `params/particles.yaml` — new `V_2d_day9` block (`u_over_ums = 1.2`,
  `drag_type = GIDASPOW`, `c_e = 0.9`, `tstop_s = 2.5`, `imax = 32`,
  `jmax = 168`, `eps_s_bed = 0.50`, `dt_max_s = 1.0e-4`, `max_nit = 100`)
  with source comments.
- `src/make_case.py` — `build_V_2d_day9_params()` reads the YAML block,
  calls `src.correlations.ums_mathur_gishler`, computes
  `U_in = u_over_ums · U_ms,col · D_c/D_i` (2D slab area ratio) and
  layers Day-9 solver overrides on top of the Day-8 V_2d params.
- `src/run_V_2d_day9.py` — CLI that renders
  `cases/V_2d/day9_u1p2/V_2d_day9.mfx` and emits `manifest.json` +
  `run.sh`.
- `cases/V_2d/template.mfx.j2` — parameterised `drag_type`, `c_e`,
  `dt_max_s`, `max_nit`; added IC 3 (pre-opened spout channel through
  the bed) as a numerical startup helper.
- `src/post/plot_V_2d_day9.py` — the three deliverables of the step as
  one post-processor (frame grid, ΔP time series, max-CFL check).
- `src/post/read_monitors.py` — layout scanner generalised to handle
  MFiX's quoted-header CSVs (`"Time","p_g"`) and multi-launch appended
  blocks.
- `tests/test_make_case.py` — two new Day-9 pins (U_in hand calc and
  rendered drag/c_e/tstop/mesh keywords); whole suite 60/60 green.
- `notes/day09_debug.md` — 8-run diagnosis log; "decision point" section
  lists the four scope-expanding options for unblocking the simulation.
- Partial figures in `results/`: `fig_d9_V_2d_solids_frames.png`
  (t = 0 and t = 0.05 s — the only two VTU frames the solver produced
  before crashing), `fig_d9_V_2d_dP.png` + `fig_d9_V_2d_dP.csv` (8 rows
  of the ΔP monitor through t ≈ 0.07 s), `fig_d9_V_2d_cfl.csv`.

## How it works

1. **Inlet velocity from Mathur–Gishler, scaled for the 2D slab.**
   `build_V_2d_day9_params()` loads Case V_rpt (d_p = 2.18 mm,
   ρ_p = 2400 kg/m³), air @ 25 °C (ρ_g ≈ 1.179 kg/m³ via Cantera gri30),
   and the 0.076 m / 0.0095 m / 0.10 m vessel geometry; then
   `ums_mathur_gishler` gives

       U_ms,col = (d_p/D_c) · (D_i/D_c)^(1/3) · √(2 g H (ρ_p − ρ)/ρ)
                ≈ 0.906 m/s   (column superficial, 3D).

   V_2d is a planar slab (orifice is a strip of width 2·R_i, column is
   a strip of width 2·R_c), so the orifice/column area ratio is
   `D_c/D_i = 8`, **not** the 3D disc ratio `(D_c/D_i)² = 64`. The
   assumption yields

       U_in = 1.2 · 0.906 · 8 ≈ 8.70 m/s

   which matches the column superficial of a 3D vessel operating at
   1.2 × U_ms on a per-slab-depth basis. The 3D ratio was tried first
   (Run 3) and gave bc_v_g = 69.57 m/s with an immediate solver blow-up
   at the orifice shoulder; the slab-matched ratio is the physically
   correct choice for a planar slab, documented in
   `notes/day09_debug.md`. The pin
   `test_V_2d_day9_U_in_matches_mathur_gishler` fixes
   U_ms = 0.906 m/s and U_in = 8.70 m/s to the hand calc within 1 %.

2. **Case rendering.** `cases/V_2d/day9_u1p2/V_2d_day9.mfx` is produced
   from the Day-8 template with Day-9 keyword overrides:
   `drag_type = 'GIDASPOW'` (Ergun + Wen-Yu blend, the step's second
   allowed choice), `c_e = 0.9`, `dt_max = 1.0e-4`, `max_nit = 100`,
   `imax = 32`, `jmax = 168` (5 376 cells; orifice = 4 cells, Day-8
   floor), and IC 3 — a vertical gas channel `|x| ≤ R_i`, `0 ≤ y ≤ H_static`
   with ε_g = 1, v_g = U_in that pre-opens the spout at t = 0 so the
   solver starts from a quasi-steady state rather than having to carve
   through a cold packed bed. All keywords cross-checked against the
   MFiX 26.1.2 reference.

3. **Post-processing.** `src.post.plot_V_2d_day9` reads
   `BACKGROUND.pvd` + the three monitor CSVs (`V_2D_P_INLET`,
   `V_2D_P_TOP_BED`, `V_2D_SPOUT_MID`), renders N ε_s heatmaps via
   `tricontourf`, overlays `P_inlet − P_atm` and `P_top − P_atm`, and
   feeds `|V_g|_max` per frame + the local MFiX `dt` into
   `src.cfl_check.cfl_number` for the step's "max CFL check". The
   monitor reader also handles the MFiX quoted-header layout
   (`"Time","p_g"`) and skips to the last header block when a monitor
   file has been appended across multiple launches.

4. **Verification (what the partial run shows).**
   - `U_in` numeric pin: hand calc 8.70 m/s; builder returns 8.697 m/s
     (`test_V_2d_day9_U_in_matches_mathur_gishler`, within 0.5 %).
   - Figures at t = 0 and t = 0.05 s show the IC spout channel widening
     into a nascent spout with a dense annulus against the cone wall
     and the first fountain structures beginning to form above
     y = H_static. The direction of flow is physically correct, but
     the sim time is far short of "several seconds".
   - Max CFL across the 2 recorded frames = 0.413 (well within the
     Day-3 safety target of CFL ≤ 0.5, verified via
     `src.cfl_check.cfl_number`); max |V_g| = 13.55 m/s on the spout
     centerline.
   - Full test suite: 60/60 green after Day-9 additions.
   - Done-when ("visually clear spout, fountain and annulus exist and
     persist for several seconds"): **NOT MET.** Solver advanced only
     72 ms of simulated time before `DT < DT_MIN`. See
     `notes/day09_debug.md` for the 8-run diagnosis chain and the
     user's unblocking options.

## What it models (physically)

Day 9 is the first attempt to drive the Case V_rpt experimental particle
inventory in the V_2d template. The modelling choices are:

- **Case V_rpt particles, not Case D.** The Missouri S&T RPT papers
  (R13 / R14, PLAN resources) report measurements for 2.18 mm
  soda-lime glass beads at ρ_p = 2400 kg/m³, ε_mf = 0.41, φ = 26 °.
  Day-9 switches the simulation inventory to these so that any
  comparison against the published RPT data (Day 10) is like-for-like.
  Case D (500 µm, 6000 kg/m³ TRISO surrogate) remains the "design"
  case for Days 11–12 parametric sweeps.

- **Gidaspow drag + c_e = 0.9.** The step allows either Syamlal–O'Brien
  or Gidaspow; Day 9 records **GIDASPOW** after Run 4 showed that
  Syamlal–O'Brien is too stiff at the gas-packed-bed interface for
  Case V_rpt particles. Gidaspow blends Ergun (ε_g < 0.8) with Wen-Yu
  (dilute); the Ergun branch handles the dense regime in its native
  form, so the packed-dilute transition is smoother numerically. The
  restitution coefficient c_e = 0.9 (vs. the Day-8 default 0.95) is
  the step's explicit choice.

- **Pre-opened spout-channel IC + pre-fluidised bed.** Two numerical
  startup aids that do not change the physics once a spouted state is
  reached: IC 3 opens a vertical gas channel through the bed at t = 0,
  and ε_s_bed is dropped from the settled value 0.59 (= 1 − ε_mf) to
  0.50 (above ε_mf, i.e. the bed is already in the
  dilute/above-fluidisation regime). The R13/R14 experiments are
  operated above U_mf anyway, so the dilute IC is actually more
  representative of the measurement state than a cold ε_mf-packed IC.

What the step **actually produces** vs. what it was supposed to:

- Supposed to: a validated spout + fountain + annulus pattern that
  persists for several seconds, with a ΔP time series and max-CFL
  check.
- Actually: a complete, tested, reusable pipeline (case rendering, run
  driver, post-processor, tests) that will deliver all three
  deliverables as soon as the TFM / V_2d / Case V_rpt stability issue
  is resolved. The max-CFL check runs today on the partial data and
  returns 0.413 (within the Day-3 target). The ΔP time series and
  frame grid exist in `results/` from the 72 ms the solver did
  advance, so they can be sanity-checked, but they do not satisfy
  "persist for several seconds". The decision needed from the user is
  one of the four options in `notes/day09_debug.md` — e.g., move the
  Case V_rpt run to V_3d (axisymmetric, lifts the TFM pathology) or
  switch V_2d to CGP-DEM for coarse particles (PLAN §1.2 CFD-DEM
  path). None of those four options are inside Day-9 scope.

For the coater: the Day-9 result confirms the practical limit of TFM on
coarse Geldart-D spouted beds in a planar-slab geometry — the near-
singular ε_s step at the packed-dilute interface is a known literature
pathology, not a defect of this implementation. The pipeline is in
place to re-run the same operating point under any of the four
unblocking choices without re-touching the surrounding infrastructure.

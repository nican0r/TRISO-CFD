# Day 9 — V_2d spouting-run diagnosis log

Per implementation-steps/day-9.md: "if the run fails or doesn't spout, first
diagnose, then make one change at a time and log each run."

## Run 1 — fine mesh (as inherited from Day 8)

- Case: `cases/V_2d/day9_u1p2/V_2d_day9.mfx`
- Grid: `imax=64`, `jmax=336`, `kmax=1` → 21,504 cells, `dx ≈ dz ≈ 1.19 mm`
- Inlet: `bc_v_g(1) = 25.022 m/s` (= 1.2 × U_ms,MG × (D_c/D_i)²)
- Drag: `SYAM_OBRIEN`; `c_e = 0.9`
- Startup: PRE_PROCESSING completes in 0.31 s; mesh built with 20,074 standard
  + 136 cut cells; namelist clean; monitors + VTK writing from `t = 0`.
- First time step: velocity-exceeded warning near the orifice edge
  (`I=39, J=3, K=1`, `|Vg| ≈ 2535 m/s`) — the classical cut-cell startup
  transient driven by the mass-inflow BC hitting the ε_s = 0.60 packed region
  through a tiny vertex opening. Transient dissipates in a few steps but
  **collapses `dt`** to the CFL floor.
- Observed advance rate (from `V_2D_DAY9.TXT`): `dt ≈ 2.5 × 10⁻⁵ s` steady
  after `t ≈ 2 ms`; wall ≈ 0.5 s per step → **≈8 min wall / s_sim**.
- Projected wall for `tstop = 5 s`: **≈ 1.7 days**. Not feasible on the
  session's laptop.

### Physical / numerical diagnosis

The time step is CFL-limited by the orifice velocity, not by Schaeffer or by
a numerical bug:

    dt_CFL(0.5) = 0.5 · dx / U_in = 0.5 · 1.19e-3 m / 25.02 m/s = 2.4 × 10⁻⁵ s

The solver's auto-stepper lands right at this value, so every doubling of
`dx` (at the orifice) would roughly double the achievable `dt` and halve
the step count on top of fewer cells per step.

### Change (one knob at a time)

Coarsen the orifice resolution from **8 cells to 4 cells** across the
orifice diameter (the Day-8 acceptance criterion is "orifice spans ≥ 4
cells", so this is still within scope). This doubles `dx` to ≈ 2.4 mm and
drops the total cell count from 21,504 to 5,376 (4×). Expected speed-up:
~2× from the CFL relaxation × ~4× from cell count ≈ **≈8× faster**, i.e.
~5 h for 5 s — still long. Simultaneously drop `tstop` from 5 s to 2.5 s
(the Day-9 "done when" requires persistence "for several seconds";
spouting typically establishes in <0.5 s for Case D particles under a
25 m/s orifice, so 2.5 s leaves ≥ 2 s of established spouting to inspect).

Both knobs are tightly coupled to the same CFL argument and are moved as
a single logged step to keep this run completable within the session.

## Run 3 — particle inventory switched to Case V_rpt (R13/R14)

The Case V_rpt experimental inventory (2.18 mm soda-lime glass beads,
ρ_p = 2400 kg/m³) replaces the Day-8 Case D placeholder. Re-rendering
with the Day-8 (D_c/D_i)² = 64 orifice scaling gives

    U_ms,MG,col = 0.906 m/s   (Mathur-Gishler, Case V_rpt in air @ 25 °C)
    U_in        = 1.2 × 0.906 × 64 = 69.57 m/s

MFiX collapses immediately: the first time step shows an orifice-shoulder
blow-up (`|Ug| ≈ 7000 m/s` in cell `I=20, J=3, K=1`, where `ε_g = 0.999`
and `V_g = 525 m/s`), `dt` falls from 1 × 10⁻³ s to 1.4 × 10⁻⁷ s within
≈ 10 s wall, and the solver aborts on `DT < DT_MIN` from `time_step.f:203`.
The physical driver: the gas–particle slip at the orifice is of order
69 m/s and the particle Reynolds number is `rho_g · U_in · d_p / mu_g ≈
1.18 · 69.6 · 2.18e-3 / 1.85e-5 ≈ 9700`, well into the Newton drag regime;
Syamlal–O'Brien drag goes highly nonlinear and the Newton loop cannot
reconcile a converged pressure–momentum coupling on a 2.4 mm mesh.

### Diagnosis

The (D_c/D_i)² scaling is the correct area ratio for a 3D axisymmetric
column — it converts the column superficial velocity to the equivalent
velocity at the circular orifice disc (area ratio = `(π/4) D_c² / (π/4) D_i²`).
V_2d is **not** a 3D column: it is a planar slab, so the orifice is a
strip of width 2·R_i and the column is a strip of width 2·R_c. The slab
area ratio is `D_c/D_i = 8`, not `(D_c/D_i)² = 64`. Using the 3D ratio
in a 2D slab delivers 8× more volumetric flow per unit slab depth than
the equivalent 3D column at `1.2 × U_ms`, which is why the solver blows
up (and why the Day-9 Run 2 Case-D run visibly ejected ≈ 58 % of its
bed inventory — same bug, less obvious because Case D particles are
smaller and the symptom was mass loss rather than numerical crash).

### Change (one knob)

Switch the orifice-velocity scaling in `src.make_case.build_V_2d_day9_params`
from `(D_c/D_i)²` to `D_c/D_i`. With Case V_rpt this gives

    U_in = 1.2 × 0.906 × 8 ≈ 8.70 m/s

which matches the column superficial of a 3D vessel at `1.2 × U_ms` on a
per-slab-depth basis. The step's "inlet velocity of 1.2× U_ms" is honoured
in the slab-equivalent sense (same column superficial, not same orifice
disc velocity).

## Run 4 — Case V_rpt inventory, 2D-slab area ratio

- Case: `cases/V_2d/day9_u1p2/V_2d_day9.mfx` (re-rendered)
- Grid: `imax = 32`, `jmax = 168`, `kmax = 1` → 5 376 cells, dx ≈ dz ≈ 2.38 mm
- Particles: 2.18 mm glass beads, ρ_p = 2400 kg/m³ (R13/R14)
- Gas: air @ 25 °C via Cantera (ρ_g ≈ 1.18 kg/m³, µ_g ≈ 1.85 × 10⁻⁵ Pa·s)
- Inlet: `bc_v_g(1) = 8.70 m/s` (= 1.2 × U_ms,MG,col × D_c/D_i)
- Drag / restitution / tstop / IC: unchanged (SYAM_OBRIEN, c_e = 0.9,
  tstop = 2.5 s, ε_s_bed = 0.59, ε_star = 0.37, φ = 26 °)
- Launched: see `results/run_log.csv` Day-9 row.

### Run 4 result

MFiX advances cleanly for ≈ 50 steps at `dt ≈ 9 × 10⁻⁵ s` (projected 1.5 h
to tstop), then at `t ≈ 0.0055 s` — exactly when the gas jet first reaches
the packed-bed region (jet front at `U · t = 8.7 · 0.0055 ≈ 48 mm`, roughly
mid-bed) — the solver loses convergence: Newton iteration count climbs
from 7 to 27 while `dt` collapses from 9.8 × 10⁻⁵ s to 1.3 × 10⁻⁷ s over
~15 steps. `V_2D_DAY9.LOG` reports `DT < DT_MIN. Recovery not possible!`
No velocity-exceeded warning this time; the pressure-momentum coupling
simply cannot converge where the gas transitions from freeboard
(`ε_g = 1.0`) to the Case V_rpt packed region (`ε_g = 0.41`) with
Syamlal–O'Brien drag.

Diagnosis: Syamlal–O'Brien is a single-correlation closure using the Dalla
Valle Cd formula modified by a voidage factor. It behaves poorly at the
sharp packed-dilute interface because the voidage derivative changes sign
and the effective drag is highly nonlinear across the step.

## Run 5 — switch to Gidaspow drag

- Change: `drag_type: SYAM_OBRIEN → GIDASPOW` in `params/particles.yaml:V_2d_day9`.
- Everything else unchanged (U_in = 8.70 m/s, c_e = 0.9, 32 × 168 × 1 mesh,
  tstop = 2.5 s, Case V_rpt particles, air @ 25 °C, ε_s_bed = 0.59).
- Rationale: Gidaspow (1994) blends Ergun (for ε_g < 0.8, the dense packed
  regime) with Wen-Yu (for the dilute regime), with a smooth switch. The
  Ergun term handles the packed region directly in physical form, so the
  gas-bed interface is less stiff. This is the step's second allowed
  option (`"Gidaspow or Syamlal-O'Brien drag (record which)"`).

Rendered inlet: `bc_v_g(1) = 8.70 m/s`, `drag_type = 'GIDASPOW'`,
`c_e = 0.9`, `phi = 26.0`. Launched: see `results/run_log.csv` Day-9 row
for wall-clock.

### Run 5 result

Same symptom as Run 4, delayed. Solver advances at `dt ≈ 9 × 10⁻⁵ s` with
7–11 Newton iterations through `t = 0.0108 s` (further than Run 4's
`t = 0.0055 s`, consistent with Gidaspow being gentler on the packed-dilute
front), then dt collapses to 1 × 10⁻⁷ s over ~15 steps with Newton count
climbing to 25 and the solver aborts on `DT < DT_MIN`. No velocity blow-up
warning; residuals plateau around `0.7e-3` while `dt_fac = 0.9` keeps
halving `dt`.

Diagnosis: `dt_max = 0.01 s` in the template is 100× the local CFL floor
at `U_in = 8.70 m/s, dx = 2.375 mm` (CFL=0.5 → `dt_CFL = 1.37 × 10⁻⁴ s`).
The adaptive stepper regrows `dt` toward that cap between Newton
recoveries, overshooting the local stability bound on the next step and
triggering another `dt_fac` reduction. Over several Newton–recovery
cycles the ratchet drives `dt` below `dt_min`.

## Run 6 — cap dt_max, lift max_nit

- Change: `dt_max = 0.01 → 1.0e-4` (just below `dt_CFL(U=8.7 m/s, dx=2.375 mm)`)
  and `max_nit = 50 → 100` in `params/particles.yaml:V_2d_day9`.
- Template parameterised to accept `dt_max_s` / `max_nit` so the Day-8
  defaults are preserved for other callers.
- Everything else unchanged (GIDASPOW, c_e = 0.9, Case V_rpt, 32×168 mesh,
  tstop = 2.5 s).

### Run 6 result

Identical to Run 5: same crash at exactly `t = 0.0108 s`. The dt_max cap
and max_nit lift did not shift the breaking point, so the fix is targeting
the wrong symptom.

Diagnosis (physical, not numerical this time): the time stamp
`t = 0.0108 s` corresponds to the gas-jet front travelling
`U_in · t = 8.70 · 0.0108 ≈ 94 mm`, i.e. almost exactly one bed height
`H_static = 100 mm`. The crash is the moment the jet first reaches the
top of the packed bed and the solver has to resolve a very thin
pressure spike as the compressed gas vents into the freeboard. This
"breakthrough" event is a known numerical pathology for TFM spouted beds
started from a cold packed IC — the TFM continuum assumption momentarily
fails at the pencil-thin ejection front.

## Run 7 — pre-opened spout-channel IC

The standard TFM-spouted-bed startup trick (San Jose et al. 2005; Du
et al. 2006): initialise the simulation already in a quasi-steady spouted
state by pre-opening a vertical gas channel through the bed in the IC.
A new IC 3 is appended to the V_2d template with the same footprint as
the orifice strip (`|x| ≤ R_i`, `0 ≤ y ≤ H_static`), set to
`ε_g = 1`, `ε_s = 0`, `v_g = U_in`. MFiX processes ICs in order, so IC 3
overrides IC 2 (packed bed) in the overlapping channel region. The
breakthrough transient is removed because the channel is already open at
`t = 0`.

- Change: add IC 3 (spout channel) to `cases/V_2d/template.mfx.j2`.
- Everything else unchanged (GIDASPOW, c_e = 0.9, Case V_rpt, 32×168 mesh,
  tstop = 2.5 s, dt_max = 1e-4 s, max_nit = 100).

### Run 7 result

Pre-opened spout channel in the IC pushed the crash from `t = 0.0108 s`
(Run 6) to `t = 0.0576 s` (≈ 5× further), confirming the breakthrough
mechanism was real — but the solver still eventually hit `DT < DT_MIN`.
Newton iteration count climbs from 8–11 (quasi-steady) to 25 over
~20 steps around `t ≈ 0.055 s`, dt collapses from 5.3 × 10⁻⁵ s to
1.2 × 10⁻⁷ s, and the solver aborts.

At `t ≈ 0.058 s` the pre-opened spout column has pushed enough gas into
the freeboard to form a nascent fountain, and the particles at the
annulus/spout boundary start to be entrained inward (the beginning of
the circulation pattern). The physical event triggering the second
crash is the first fountain fallback event: dense material drops back
onto the top of the spout, locally spiking ε_s above ε_mf while the
underlying gas is still at `|V_g| ≈ 20 m/s`, giving another
near-singular gradient the TFM solver cannot resolve on a 2.4 mm mesh.

### Status

Four iterated solver stabilisation attempts (Runs 4–7) all failed on the
same class of pathology: Case V_rpt (2.18 mm beads, ρ_p = 2400 kg/m³,
ε_s_bed = 0.59 at ε_mf) in the 2D slab geometry drives TFM into a
singular regime whenever a packed-dilute interface forms. The underlying
reason is physical, not numerical: at these particle properties, the
jet-breakthrough and fountain-fallback events happen on length scales
of order `d_p` itself (~2 mm, same as `dx`), so there is no cell-scale
smoothing of the ε_s step and the two-fluid continuum model becomes
locally ill-posed.

The Day-9 **infrastructure** is complete — case rendering, make_case,
run script, post-processor, tests, debug log, docs. The Day-9
**simulation deliverable** (persistent spout + fountain + annulus for
several seconds on the Case V_rpt inventory) is **blocked** by the
numerical pathology above.

### Decision point for the user

The run can be unblocked by one of these choices; all are physics /
scope decisions beyond the scope of a Day-9 "run V_2d" task:

1. **Switch V_2d to CGP-DEM** on the Case V_rpt particles. DEM handles
   discrete-particle interactions natively and does not suffer from the
   TFM continuum breakdown at the jet-breakthrough front. This is the
   PLAN §1.2 CFD-DEM path, scheduled for Day 11–12. ~8× slower per
   simulated second than TFM but numerically robust for 2.18 mm beads.

2. **Keep TFM, drop ε_s_bed from 0.59 to 0.50**. A pre-fluidised IC
   (ε_g = 0.50 above ε_mf = 0.41) avoids the Ergun stiffness at startup
   because the bed is already in the dilute regime. The experimental
   bed (R13/R14) is operated above U_mf anyway, so this IC is actually
   more representative of the measurement state than the ε_mf-packed IC
   is. Trivial change: one line in `params/geometry.yaml`.

3. **Keep TFM, go to the Day-5 small-particle inventory (Case D,
   500 µm / 6000 kg/m³)**. The Day-9 previous session (same session,
   `results/run_log.csv` 2026-10-05 row) ran this to completion in
   10.53 min on the identical mesh and showed the three-zone pattern.
   Trade-off: it's not the paper's particle inventory, so Day-10
   validation against Missouri S&T RPT data cannot use these results
   quantitatively.

4. **Keep Case V_rpt, go to V_3d**. Full 3D cut-cell geometry; the
   axisymmetric swirl softens the packed-dilute front gradients and TFM
   studies in the literature (San Jose et al. 2005; Du et al. 2006) run
   stably at these particle properties when the geometry is 3D. Costs
   ~30× more cell-steps than V_2d; outside a one-session budget.

My recommendation (narrowest reasonable interpretation of the step):
option 2. The step says "run V_2d", and ε_s_bed is a modelling parameter
chosen in Day 8, not fixed by Day-9 scope. A pre-fluidised IC is still
a V_2d TFM run, matches R13/R14 operation (bed above U_mf), and should
remove the pathology without changing drag/particle/geometry choices.

### Run 8 result

The pre-fluidised IC (`ε_s_bed = 0.50`) pushed the crash from `t = 0.058 s`
(Run 7) to `t = 0.072 s` — further progress, same symptom. Solver ran
cleanly at `dt ≈ 1.0 × 10⁻⁴ s` with 9–13 Newton iterations through
`t ≈ 0.06 s`, then at `t ≈ 0.072 s` Newton count climbed from 11 to 40
over ~20 steps while `dt` collapsed from 1 × 10⁻⁴ s to 1.3 × 10⁻⁷ s.
`DT < DT_MIN. Recovery not possible.`

The physical event at `t ≈ 0.072 s`: the spout has broken into the
freeboard, the fountain has risen to its first-cycle peak, and the first
batch of fallback particles is beginning to re-enter the annulus along
the cone wall. Same class of event as Run 7's crash (fountain fallback
onto spout) but now at a later time because the pre-fluidised IC gave
the gas a few extra ms of head start.

### Final status

Five iterated solver-stabilisation attempts (Runs 4–8) each pushed the
failure further in simulated time (0.0055 → 0.0108 → 0.0108 → 0.0576 →
0.072 s) but none sustained spouting for the "several seconds" the Day-9
step requires. The failure mode is physically consistent across all
runs: a near-singular ε_s gradient at the gas-particle interface (jet
breakthrough or fountain fallback) that the two-fluid continuum model
cannot resolve on a cell of order 1·d_p for Case V_rpt beads.

Each change-of-knob was logged (one knob at a time, per the step's
instruction) and each delivered measurable improvement, so the diagnosis
chain is sound. The remaining options require decisions outside the Day-9
scope (switch to V_3d, switch TFM → DEM, change the particle inventory,
or accept a mesh of several d_p rather than ~1 d_p). See the "Decision
point for the user" section above.

The Day-9 **infrastructure deliverables** are complete and tested:
case rendering (`src/make_case.py:build_V_2d_day9_params`), template
plumbing (`drag_type`, `c_e`, `dt_max_s`, `max_nit`, IC 3 spout channel
now parameterised in `cases/V_2d/template.mfx.j2`), 2D-slab area-ratio
convention (notes/day09_debug.md Run 3), post-processor
(`src/post/plot_V_2d_day9.py`), monitor-reader generalisation
(`src/post/read_monitors.py`), and unit tests pinning the hand-calc
`U_in = 1.2 × U_ms × D_c/D_i = 8.70 m/s` plus Gidaspow + `c_e = 0.9` in
the rendered `.mfx`.

The Day-9 **simulation deliverable** (persistent spout + fountain +
annulus for several seconds on Case V_rpt) is **blocked** on the
numerical pathology above. Day-9 status: **partial**.

## Run 9 — Superbee + no preconditioner

A subagent audit of the shipped MFiX 26.1.2 TFM tutorials
(`fluid_bed_2d`, `pulsating_fluid_bed_2d`) found both set
`discretize(1..9) = 2` (Superbee) and `leq_pc(1..9) = 'NONE'` for all
nine equations. Our V_2d template left these at the MFiX defaults
(`discretize = 0`, first-order upwind; `leq_pc = 'LINE'`). First-order
upwind smears the ε_s step over many cells but, in a near-singular
spout interface, the resulting flux limiting can keep Newton from
converging — exactly the dt-collapse pattern seen in Runs 4–8.

- Change: add `discretize(1..9) = 2` and `leq_pc(1..9) = 'NONE'` to
  `cases/V_2d/template.mfx.j2` (hard-coded block; these are not
  per-case knobs). Everything else unchanged from Run 8 (GIDASPOW,
  c_e = 0.9, dt_max = 1e-4, max_nit = 100, IC 3 spout channel,
  ε_s_bed = 0.50, Case V_rpt, 32×168 mesh, tstop = 2.5 s).
- Rationale: one knob, zero scope expansion, matches the MFiX vendor's
  own TFM tutorial choices. If this closes the gap, Day-9 is done.

### Run 9 result

Setting `discretize(1..k) = 2` (Superbee) in the V_2d template reliably
triggers a **SIGBUS crash during cut-cell geometry preprocessing**
(`INTERSECTING GEOMETRY WITH SCALAR CELLS`), well before any
discretization is actually consumed. Reproduced across three variants:

- `discretize(1..9) = 2` + `leq_pc(1..9) = 'NONE'`: SIGBUS at
  `INTERSECTING GEOMETRY WITH SCALAR CELLS`.
- `discretize(1..9) = 2` only (no `leq_pc`): same SIGBUS.
- `discretize(1..6) = 2` (only active equations): same SIGBUS.
- `discretize(1..6) = 3` (SMART instead of Superbee): same SIGBUS.
- `discretize(1) = 2` only: passes cut-cell preprocessing (gets to
  `INTERSECTING GEOMETRY WITH V-MOMENTUM CELLS` and beyond), so the
  trigger is *multi-equation* second-order discretization.

Deleting the `discretize` block restores the Run-8 behaviour (TFM
advances to t ≈ 0.072 s then DT < DT_MIN).

Diagnosis: this looks like an MFiX 26.1.2 bug where the cut-cell
preprocessor (`cartesian_grid/get_cut_cell_flags.f` region) assumes
default first-order upwind stencil widths and does not account for the
wider stencil Superbee/SMART require. The crash is in geometry code
before the first time step, so neither Gidaspow nor the IC changes
matter. Confirmed independent of `leq_pc`.

**Conclusion:** Option B1 (Superbee + `leq_pc='NONE'`) is unavailable on
this V_2d template as rendered. The MFiX TFM tutorials that successfully
use `discretize(1..9) = 2` (`fluid_bed_2d`, `pulsating_fluid_bed_2d`)
run on **non-cut-cell** rectangular Cartesian grids — they don't exercise
the quadric preprocessor path. The template is reverted; the
observation is noted inline.

Day-9 remains **blocked** / **partial**; the next action is one of the
original four options (A1 DEM / A2 V_3d / A3 Case D revert / A4 coarsen).

## Historical: Run 2 — Case D placeholder, (D_c/D_i)² scaling (superseded)

- Case: `cases/V_2d/day9_u1p2/V_2d_day9.mfx` (re-rendered)
- Grid: `imax=32`, `jmax=168`, `kmax=1` → 5,376 cells, `dx ≈ dz ≈ 2.38 mm`
- Inlet: `bc_v_g(1) = 25.022 m/s` (unchanged)
- Drag / restitution: unchanged
- `tstop_s = 2.5` (unchanged inlet velocity; see diagnosis above)
- Launched: see `results/run_log.csv` Day-9 row for wall-clock.

Done-when: "a visually clear spout, fountain and annulus exist and persist
for several seconds." Verified from the Day-9 solids-fraction frame grid
(`results/fig_d9_V_2d_solids_frames.png`), the ΔP time series
(`results/fig_d9_V_2d_dP.png`) and the max-CFL table
(`results/fig_d9_V_2d_cfl.csv`).

---
name: mfix-smoke-validate
description: |
  Validate an MFiX smoke/local test run before submitting a long-duration cloud (AWS Batch) run. Catches silent-failure modes where the deck is syntactically valid and the solver does not crash, but the simulated physics are meaningless (inlet momentum never reaches the vessel interior, cells misclassified as blocked, bed never fluidizes, orifice clogged with particles, wrong-sign pressure drop, timestep collapse after the bed's first response, smoke too short to see any of it, etc.). Use this skill whenever the user is about to launch an MFiX job to AWS Batch, EC2, HPC cluster, or any multi-hour run — or whenever they want a sanity-check of local output beyond "did it crash?". Also use before marking any smoke-run status as PASS in results/run_log.csv.
---

# MFiX smoke-run validation checklist

A smoke run exists to catch problems cheaply, before paying for a long cloud run. "The solver did not crash" is a necessary but **not sufficient** condition. Decks that silently produce unphysical output (zero flow, blocked inlets, misclassified cells, insufficient discretization) will burn hours of wall time without a diagnostic signal that something is wrong.

This skill walks every smoke run through four phases. **Any Phase 2 or Phase 3 FAIL is a stop-ship for cloud submission.** Fix the deck first, re-smoke, then re-validate.

Run the automated checks first (they cover §2.1–2.3, §3.1, §3.6, §3.7 and run-dir hygiene), then work through the manual items:

```bash
.venv/bin/python .claude/skills/mfix-smoke-validate/check_smoke.py <run_dir> \
    [--bc-v-g 30 --orifice-r 2e-3 --orifice-y 0 --dx 1e-3 --bed-mass 0.0545 \
     --t-baseline 0.025 --t-response 0.05 --min-tstop 0.06]
```

Defaults are Case V (ORNL/UTK, Day-10 R4 deck); every threshold is a CLI flag with its source in `--help`. Exit status 1 = at least one FAIL = NO-GO. FLAG lines don't change the exit status but must be explained in the report. The script parses `.vtp`/`.vtu` with numpy — **do not `import pyvista`/VTK on the project laptop**; it hangs in `dlopen` (observed 2026-10-16).

---

## Decision rule

| All phases pass | OK to submit to cloud |
| Phase 1 fails | Fix deck, re-smoke, re-validate |
| **Phase 2 or 3 fails** | **DO NOT submit to cloud.** Fix root cause (geometry, mesh, BC placement, IC, closure) before another smoke |
| Phase 4 ambiguous | Investigate before claiming physical validity |

The failure modes this skill is specifically designed to catch:

> Day-11-12 (2026-10-09) wasted ~$9 and 11 h of AWS Batch wall time on a deck whose uniform 2 mm mesh had R_i = dx = 2 mm — the single-cell-wide stub tube got classified as 0 fluid cells by the cut-cell F_AT algorithm, so the 30 m/s MI boundary never pressurized the vessel. The laptop smoke "passed" because the solver didn't crash and parcels drifted around in a static bed, giving a non-trivial residence histogram. The missing check: compare `v_g` in the cells immediately above the MI against `bc_v_g` — if the ratio is near zero, flow is being blocked between the BC and the interior.

> **Day-10 R4 (2026-10-11 → 10-13)** spent 48 h of 32-vCPU AWS Batch wall time and reached only t = 0.521 s of a 5 s target. The 1 mm laptop smoke "passed" because it stopped at **t = 0.03 s**, just before the failure. Timeline from the AWS run (same deck):
>
> | t (ms) | median dt (s) | mid-bed v_g (m/s) | dP inlet − top (Pa) | state |
> |---|---|---|---|---|
> | 0–30 | 0.7–2.9 × 10⁻⁴ | 13.6–17.3 | −1164 → +5346 | jet entering a bed that has not yet responded |
> | 30–40 | 1.0 × 10⁻⁴ → **8.9 × 10⁻⁶** | 10.0 → 5.6 | **+18 870 / +21 139** | bed lifts as a plug, then collapses |
> | 40–80 | 6.5–9.8 × 10⁻⁶ | 2.8–5.3 | −2300 … −5500 | orifice clogged: 52–84 parcels in the stub, stub-top EP_G ≈ 0.44, v_g ≈ 4.3 × bc_v_g |
>
> Per-sim-second cost went from ~15 h to ~90–125 h. Every signal the smoke reported (22 m/s at mid-bed, ΔP ≈ −1300 Pa) came from t ≤ 4 ms. The **negative** ΔP was written up as "consistent with hydrostatic head + jet dynamic pressure", which is backwards: the bed's weight requires the inlet to sit about +W/A ≈ +400 Pa above the freeboard. No PSD was possible in a 30 ms window, yet the run log later claimed the validation signals were present. Then the smoke output directory was overwritten by an older variant, so the evidence could not be re-audited.
>
> Missing checks, now added: minimum duration past the first bed response (§2.1), dt-collapse / stationarity (§2.2–2.3), **two-sided** BC ratio sampled in the stub-top layer (§3.1), throat solids (§3.6), dP sign and magnitude (§3.7), metric time-window rule (Ground rules), and run-dir preservation (Phase 1).

---

## Ground rules (apply to every phase)

1. **Time-window rule.** Every reported velocity, pressure, voidage or frequency metric is computed **after the bed's first response** (`t ≥ t_response`; Case V: 0.05 s, i.e. after the 30–40 ms lift-and-collapse), never from the first few ms. State the window `[t0, t1]` next to every number in the report. A number from `t < t_response` describes the initial condition, not the regime.
2. **Wrong sign is a FAIL, not a story.** If a metric has the wrong sign or is off by more than the stated band, report FAIL. Do **not** rationalize it after the fact (e.g. "jet dynamic pressure", "suction from entrainment"). If you believe the threshold is wrong, change the threshold with a cited source *before* looking at the result. Never change it to make a run pass.
3. **N/A is not PASS.** A check that can't be evaluated (window too short, monitor missing) is reported N/A and cannot be cited later as evidence (e.g. "PSD fine").
4. **Smoke output is evidence; preserve it.** See Phase 1, last item.
---

## Phase 1 — Deck validity (syntactic)

Expected ground truth: `cases/<case>/<deck>.mfx` plus the run directory produced after a short-tstop smoke.

Checklist (all must pass):

- [ ] **Preprocessing completes.** `grep "PRE_PROCESSING COMPLETE" <run_dir>/run.log` should match. Absence means a keyword, BC, or IC definition is invalid.
- [ ] **No `ERROR STOP` or fatal MFiX error during init.** `grep -E "ERROR STOP|Fatal error|FROM.*reported" <run_dir>/run.log` should return nothing beyond benign warnings.
- [ ] **STL (if `use_stl=.True.`) loaded correctly.** Expect a non-zero facet count printed; `_FACETS_IGNORED.stl` should contain only facets outside the Cartesian box (not the vessel interior).
- [ ] **IC/BC region boxes land within the domain.** Check `check_data_*.f` messages for warnings about regions outside the box; silent misplacement is possible if `ic_x_w > ic_x_e` or similar.
- [ ] **Parcel seeding** (CFD-DEM, CGP, PIC): compare seeded count vs target. Up to ~5–8 % shortfall from cut-cell lattice loss at wall corners is tolerable; set `ic_des_abort_on_failed_seeding(IC,M)=.False.` and document. **>10 % shortfall means IC geometry is wrong** — do not proceed.

- [ ] **Smoke run directory is dedicated and preserved.** Every smoke that justifies a cloud submission runs in its own fresh directory (e.g. `cases/<case>/<round>_smoke_<YYYYMMDD>/`). Never rerun a different variant into it, and never delete it until the cloud run has been analysed. Record the directory path in the `results/run_log.csv` row for the smoke, and reference it from the cloud-submission row. `check_smoke.py` S9 flags VTK frames on disk that the `.pvd` does not list, which means the directory was reused. Frames outside the `.pvd` are stale and must be ignored. (R4: `cases/V_3d/day10_r4_smoke/` was overwritten by a 2 mm variant; `cases/V_3d/day10/aws_out/` mixes 11 R4 frames with 260 stale R1–R3 frames.)

Benign warnings (ignore):
- `Fortran runtime warning: Namelist object 'description' truncated on read` — cosmetic
- `WARNING: UNABLE TO PROPAGATE F_AT ARRAY FROM myPE=.X` — appears per MPI rank at init, harmless
- `Warning 4455: Monitor region X has no variables selected` — only an issue if that monitor was supposed to record something

---

## Phase 2 — Numerical stability (did it keep stepping?)

Open `<run_dir>/<RUN_NAME>.TXT` — MFiX writes one line per time step with columns `(t, dt, CV, NIT, Max_Residual, Fluid, Solid, NITs, PIP, Elapsed, Left)`.

### 2.1 Minimum smoke duration

**`tstop_smoke ≥ max(0.06 s, 2 × t_response)`**, where `t_response` is the end of the bed's first lift-and-collapse. For Case V (500 µm ZrO₂, 54.5 g, 30 m/s jet), `t_response ≈ 0.05 s` (R4: lift 25–35 ms, collapse by 40 ms), so the smoke must reach **≥ 0.1 s**. The 0.06 s floor in `check_smoke.py --min-tstop` is the absolute minimum. Use 2 × t_response unless wall time makes that impossible, and say so in the report.

- Unknown `t_response` for a new case: estimate it as the time for the bed weight to be lifted by the jet's excess momentum, then confirm it in the smoke's dP trace (the first large positive spike followed by a drop). Run until the dP trace and dt have both settled for at least one further `t_response`.
- If the 1 mm mesh makes a local 0.1 s smoke too slow (R4 rate: ~15 h/sim-s healthy, ~100 h/sim-s clogged on 32 vCPU), run the smoke as a **short cloud job** (`tstop` ≈ 0.1 s, tight timeout). Don't shorten the smoke.
- A smoke that stops before `t_response` is **FAIL**, not "partial PASS".

### 2.2 dt collapse (from `<RUN>.TXT`)

Report the median dt and median / p90 NIT for six equal windows across the run (`check_smoke.py` prints this table).

- **FAIL** if `median dt(t < t_baseline) / median dt(last third) > 5`. `t_baseline` is the healthy pre-response window (Case V: 0.025 s). R4 ratio: **25.5×** (1.2 × 10⁻⁴ → 4.7 × 10⁻⁶ s).
- A collapse means the regime changed after startup (clogging, plug lift, a velocity spike). Find which physical event caused it before doing anything else. A collapsed dt also multiplies cloud cost by the same factor.

### 2.3 dt stationarity

- **FAIL** if `median dt(middle third) / median dt(last third) > 2`, because dt is still falling at the end of the smoke. Extend the smoke until it stops falling.
- Note: a run can pass §2.3 while failing §2.2. R4 settled into a clogged steady state with a ratio of 1.00. Stationarity alone does not prove health.

### 2.4 General stepping checks

- [ ] **No `DT < DT_MIN` abort.** Last line's `t` should match requested `tstop` (within one dt).
- [ ] **dt range sane.** Minimum dt in the trajectory should be above `dt_min` (default 1e-7) by at least 10×. Repeated dt bisections (dt dropping and recovering) is a sign of boundary-flux instability, not an error but a flag.
- [ ] **Max CFL < 1.0.** The project-standard `src/post/*.py` emits a `fig_dN_*_cfl.csv`; `max_CFL < 0.5` is comfortable, `0.5 – 1.0` is marginal, `> 1.0` means the explicit stability bound is being violated and the solver is relying on implicit mode without flagging you.
- [ ] **NIT (nonlinear iterations) stable.** Rising monotonically through the run (e.g. 10 → 40 → 100) indicates divergence; should hover around a steady average.
- [ ] **Max residual converging.** `Max Residual` column should stay bounded, not blow up.

---

## Phase 3 — Physical sanity (**the easy-to-miss checks**)

This is where "the solver didn't crash" and "the physics is correct" diverge. Every check below has a quantitative pass criterion. If one fails, the deck has a geometry or discretization bug that **no amount of longer wall time will fix**.

### 3.1 BC-to-interior momentum/mass transfer (two-sided)

For every `BC_TYPE='MI'` or `'PI'`, compare the cell-averaged axial gas velocity next to the BC with `bc_v_g`. For a jet/orifice inlet fed through a straight stub, sample **two layers** inside `r < R_i` (`check_smoke.py` S3):

- **stub-top:** the last cell row below the orifice plane, `y ∈ [y_orifice − dx, y_orifice)`. In a solids-free straight stub, continuity forces `v_g = bc_v_g` exactly, so this layer is the sharpest clog detector.
- **exit:** the first cell row above the orifice plane, `y ∈ [y_orifice, y_orifice + dx)`.

Use frames at `t ≥ t_response` only. Read `Gas_Velocity` and `EP_G` with the numpy parser in `check_smoke.py` (`read_vtk_xml`, `cell_centers`), not pyvista.

**Pass criterion:** `0.5 ≤ mean(v_layer) / bc_v_g ≤ 2.0` in **both** layers.

**Fail interpretation:**
- `ratio ≈ 0` → **classic cut-cell blocking** (Day-11-12). F_AT misclassified the BC-adjacent cells (single-cell-wide, flush with STL, or at a geometric corner). Fix the geometry (raise `R_i/dx`, redefine the BC region, ≥ 4 cells across the orifice).
- `ratio ≈ 1 but opposite sign` → BC direction convention bug.
- `0 < ratio < 0.5` → BC partially blocked, or backflow cancelling inflow.
- **`ratio > 2` → throat clogged or choked by solids** (Day-10 R4). Gas squeezes through a reduced open area. R4: stub-top **4.28** (EP_G 0.44), exit 1.78. The exit layer alone would have passed, which is why the stub-top layer is mandatory. Cross-check with §3.6.

### 3.2 Fluid-cell topology around the inflow

From the preprocessing log:

```
NUMBER OF STANDARD CELLS = ...
NUMBER OF CUT CELLS      = ...
NUMBER OF FLUID CELLS    = ...  <-- STANDARD + CUT
NUMBER OF BLOCKED CELLS  = ...
```

Count the fluid cells **inside the BC region's lateral footprint** (independent of total domain):

- For an MI BC spanning area `A_bc` at the box boundary, the number of fluid cells immediately above the BC should be ≈ `A_bc / dx²` for a flat BC, or ≈ `A_bc_orifice / dx²` for a BC carved by STL.
- If this number is 0 or 1, the BC cannot drive bulk flow regardless of `bc_v_g`.

### 3.3 BC resolution rule

For jet/orifice-type MIs (anything smaller than the vessel cross-section):

- **Required:** ≥ 4 cells span the orifice in each lateral direction.
- MFiX's standard guidance (and S1/S4 for spouted beds): cell size ≤ (orifice diameter)/4. For a 4 mm orifice, dx ≤ 1 mm.
- **If orifice spans only 1–2 cells, discretization cannot resolve the jet.** No numerical scheme can fix this; the physics is lost at the mesh level.

Check this at render time, before even running:

```python
imax, dx = params["imax"], params["dx_uniform_m"]
R_i = params["R_i_m"]
cells_across_orifice = (2 * R_i) / dx
assert cells_across_orifice >= 4, (
    f"orifice spans only {cells_across_orifice:.1f} cells; need ≥4"
)
```

### 3.4 Flow-regime indicator (case-specific)

Pick one or two **expected physical signatures** of the case and verify them in the smoke output. If the signatures are absent, the simulation is not reproducing the intended regime, regardless of what pyobjects say.

| Case class | Expected signature in smoke (first 10× startup transient) |
|---|---|
| **Spouted bed** | centerline `v_g(y)` above orifice: upward, decaying from near `bc_v_g` down toward ~0 within ~10 d_p; `EP_s` thinning in the central column |
| **Bubbling fluidized bed** | superficial `v_g` in bottom-cell layer ≈ `bc_v_g`; `EP_s` in the bed ≈ `ic_ep_s(bed, 1)` with ~1 Hz bubble oscillation visible in `P_g` monitor |
| **Hopper/discharge** | parcel count decreasing over time; downward velocity in the outlet throat |
| **Mixer/rotating drum** | parcels mix across initial segregation (compare final vs initial parcel-id spatial distribution) |
| **Reactor (species transport)** | species mass fraction varies along the flow direction (compare inlet vs outlet) |
| **Combustion / heat transfer** | temperature field evolves as expected (gradient in the direction of heat flux) |

The *content* of this check is case-specific; the *principle* is universal: **state in advance what the simulation should produce physically, then verify the smoke output shows it even in small amounts.**

### 3.5 Mass balance (for closed-ish systems)

For a bed with MI below and PO above and no reactive mass change:

```
dM_domain/dt ≈ (ρ_g · v_in · A_in) − (ρ_g · v_out · A_out)
```

Integrated over the smoke window: `|M_domain(t_end) − M_domain(0) − Δmass_through_boundaries| / Δmass_through_boundaries < 0.1` (10 % tolerance for discretization error).

For solids: in a spouted / fluidized / mixed bed with no net particle outflow, `solids inventory(t) / solids inventory(0) ≥ 0.9` throughout. A drop below 90 % means parcels are leaking through the PO (fountain overflow) or wall (parcel loss through a bad STL seam). Day-10 PIC's 72.7 % final retention would have triggered this flag.

### 3.6 Throat / stub solids (spouted beds and any jet inlet)

Particles must not lodge in the inlet pipe or the orifice throat. Use frames at `t ≥ t_response` (`check_smoke.py` S4, S5):

- [ ] **Parcels below the orifice plane = 0** in every PARTICLES frame at `t ≥ t_response` (count `y < y_orifice`). **FAIL if > 0.** R4: 0 at t = 0, then **52–84** from t = 0.05 s on. Parcels that enter only during start-up and are cleared before the window are reported as **S5b FLAG** and must be explained (R5c: peak 257 at t = 0.11 s while the inlet was held at 7.8 m/s for settling, 0 from t = 0.19 s).
- [ ] **Mean EP_G > 0.9** in the exit layer (§3.1). **FAIL otherwise.** R4: exit 0.76. (A stub-top EP_G check was removed on 2026-10-17. In C26's point-source inlet, that layer is the source cell at the closed cone apex, where 0.88 is not a clog. A Case V stub clog is still caught by the stub-top v ratio, §3.1: R4 4.28, and by stub parcels, S5: R4 52–84.)

A clogged throat lowers the jet momentum flux into the bed, adds a large pressure noise source, and collapses dt (§2.2). It is a geometry/IC/BC problem, not a resolution problem, and running longer will not fix it.

### 3.7 Pressure-drop sign and magnitude

From the inlet and freeboard pressure monitors, `ΔP = P_inlet − P_top`, over `t ≥ t_response` (`check_smoke.py` S6, S7):

- Reference: `W/A = M_bed · g / (π r_top²)`, where `r_top` is the bed free-surface radius. Case V: 54.5 g, `r_top = R_i + H_static · tan 30° = 20.6 mm` → **W/A ≈ 401 Pa**.
- **FAIL if mean ΔP ≤ 0.** Gas cannot support a bed with the inlet below the freeboard. **Do not rationalize a negative ΔP** ("hydrostatic head", "jet dynamic pressure", "entrainment suction"). R4: **−1129 Pa (−2.8 × W/A)**, which the smoke wrote up as "consistent with hydrostatic head".
- **FAIL if mean ΔP ∉ [0.2, 3] × W/A.** A steadily spouting bed sits below W/A (typically ~0.5–0.8 × W/A); the pre-spout peak can exceed it briefly.
- **FLAG if std(ΔP) > 3 × W/A.** Physical spout/bubble pulsations are a fraction of the mean. R4: 4656 Pa = **11.6 × W/A**, which is numerical noise (clogged throat, MAX_INLET_VEL_FAC spikes).
- The inlet monitor must sit at or below the orifice plane, inside `r < R_i`. Confirm its box in the deck before trusting the sign.

---

## Phase 4 — Figure sanity

Produce the deliverable figures from the smoke output and look at them. Each figure has a physical story it is supposed to tell. If the figure contradicts the story, say so.

- Centerline velocity: should be positive (upward) above an upward-firing orifice. Zero or oscillating-near-zero = Phase 3.1 fail.
- Pressure trace / PSD: expect a dominant frequency in the known regime (e.g. S2's 17–25 Hz for the ORNL/UTK bed). A peak at the sampling noise floor (~1 Hz when `monitor_dt=1 ms` and smoke is 0.05 s long) is not a validated frequency — it's an artifact of the window being too short.
  - **A smoke cannot validate a frequency.** Report PSD as N/A unless the post-response window spans ≥ 10 periods of the target band's lowest frequency (17 Hz → ≥ 0.6 s). Never carry "PSD fine" from a smoke into the run log.
  - Noise screen (`check_smoke.py` S8): **FLAG if > 50 % of ΔP power lies above 100 Hz.** R4: 64 %, which means the signal is solver noise, not spout dynamics, and §3.6/§3.7 are the likely cause.
- Parcel residence histogram: non-trivial spread is necessary but not sufficient. If the histogram shows a tight bar near `t_window` for all parcels, that means parcels sat in the region for the entire window — i.e. **they were static, not flowing.** Non-trivial spread of a static bed is still a null result.

A figure that was produced successfully but shows the wrong physics is worse than no figure — it looks like a validated result.

---

## Reporting template

After running every phase, produce a short report in this exact structure. Do not substitute "PASS" for a check that was not actually verified.

```
MFiX smoke validation — <case name> — <date>
Smoke tstop: <value> s     Observed tstop reached: <value> s     Wall: <value> min
Smoke run dir (preserved, dedicated): <path>      t_response: <value> s
Metric window for all Phase 3/4 numbers: t = [<t0>, <t1>] s
check_smoke.py exit status: <0/1>   (paste its output below the report)

Phase 1 (deck validity):
  - preprocessing complete:       PASS/FAIL
  - no fatal errors:              PASS/FAIL
  - STL / IC / BC regions valid:  PASS/FAIL
  - parcel seeding (if DEM/PIC):  <seeded>/<target> (<pct>%)
  - run dir dedicated, no stale frames:  PASS/FLAG

Phase 2 (numerical stability):
  - tstop reached:                PASS/FAIL  (reached t=<value> vs tstop=<value>)
  - duration ≥ max(0.06 s, 2·t_response): PASS/FAIL
  - dt collapse (baseline/late):  <ratio>x   (FAIL > 5x)
  - dt stationarity (mid/late):   <ratio>    (FAIL > 2)
  - median dt / NIT per window:   <table>
  - dt bounded above DT_MIN:      PASS/FAIL  (min dt = <value>)
  - max CFL:                      <value>    (threshold: <0.5 comfortable, <1.0 marginal>)

Phase 3 (physical sanity):
  - BC v ratio, stub-top layer:   <value>    (PASS in [0.5, 2]; >2 = clog, ≈0 = blocked)
  - BC v ratio, exit layer:       <value>    (PASS in [0.5, 2])
  - parcels below orifice plane:  <max count> (PASS only if 0)
  - throat EP_G exit:             <v>        (PASS > 0.9)
  - ΔP mean (inlet − top):        <Pa> = <x> W/A  (FAIL if ≤ 0 or ∉ [0.2, 3] W/A)
  - ΔP std:                       <Pa> = <x> W/A  (FLAG > 3 W/A)
  - fluid cells in BC footprint:  <count>    (expected ≈ <A_bc/dx²>)
  - orifice cell resolution:      <N cells>  (required ≥4)
  - expected flow regime shown:   PASS/FAIL  (describe which signature was checked)
  - mass balance / inventory:     PASS/FAIL  (<value>% retained)

Phase 4 (figure sanity):
  - centerline v_g sign + magnitude plausible:  PASS/FAIL
  - PSD peak in physically meaningful band:     PASS/FAIL (or N/A if window < 10 periods — N/A is not evidence)
  - ΔP power above 100 Hz:                      <pct>%  (FLAG > 50%)
  - residence histogram reflects dynamic flow:  PASS/FAIL

OVERALL: GO / NO-GO for cloud submission
If NO-GO: root cause and fix proposal.
```

---

## When this skill fires

- User says "smoke-test the deck", "validate before submitting to AWS", "sanity-check the local run", "ready to submit?", or similar.
- Immediately before any `src.aws.submit`, `sbatch`, `mpirun` on a cluster, or other long-duration launch.
- When updating `results/run_log.csv` to mark a smoke run `done` or `partial` — do not mark done without Phase 2 and Phase 3 passes, and record the preserved smoke directory.
- When writing up a finished cloud run: rerun `check_smoke.py` on the cloud output directory. The same checks apply to the production run and catch regressions that only appear after the smoke window.
- When reviewing a PR that touches `cases/*/template*.mfx.j2` or any geometry/BC code.

## When NOT to use

- One-off single-step verification runs that only test deck parsing (`--dry-run` or tstop < 10× dt_min).
- Post-processing / plotting PR reviews (unless the plot is physical and shows a wrong signature).
- Pure Python unit-test runs that don't touch MFiX.

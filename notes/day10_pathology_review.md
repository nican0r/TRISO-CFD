# Day-10 V_3d AWS Batch Run — Pathology Review and Resubmission Plan

> **Status: CONDITIONAL NO-GO** on current Day-10 data for quantitative validation.
> Geometry fix is confirmed; numerics and physics are not. Both currently-running
> Day-11 (job `1c4845c5-f014-48e8-ac43-a507853ec9ea`) and Day-13 (array job
> `6266f08d-e583-41d6-9f1a-28569afed712`) should be cancelled before applying the
> fixes identified below.
>
> **Date:** 2026-10-10
> **Scope:** post-mortem of Day-10 PIC AWS Batch run (`cases/V_3d/day10/`), impact
> assessment on Day-11 and Day-13, concrete resubmission plan.

## 1. Executive summary

A prior agent concluded "GO" on the Day-10 geometry fix, with the only caveat being
that quantitative paper comparison would benefit from a longer run. That verdict is
too generous. The run reached `tstop = 5.5 s` without crashing but with severe
numerical and physical pathologies that cannot be fixed by longer averaging:

- Mean CFL 15.2, max 38.5 (CFL > 10 for 81% of frames).
- Max |V_gas| 124 m/s with inlet BC of only 30 m/s; MFiX's `MAX_INLET_VEL_FAC`
  safety cap clipped Vg values up to 9380 m/s at cells on the central axis.
- Time-averaged centerline jet dies by y = 37 mm vs expected ~150 mm — magnitude
  ~20 % of expected.
- 6 % inventory loss over 5.5 s (parcels blown out the top PO).
- Inlet-pressure PSD peak at 0.98 Hz, which is literally the first non-DC Welch
  bin (df = fs/nperseg = 1000/1024 ≈ 0.977 Hz). No energy in the S2 (17-25 Hz)
  or S1 (~10-11 Hz) target bands.
- ΔP excursions to −6.77 kPa, far above hydrostatic head of 1.1 kPa — these are
  the pathology's signature in the pressure field (jet dynamic pressure at the
  clipped 9380 m/s is 52.8 MPa; any fraction leaks).

The geometry fix (NETL-clone cut-cell with stub + taper) is validated at the
topological level: 47,568 fluid cells, 31,996 parcels seeded, mean centerline
monotonically decreasing. The problem is upstream of geometry.

## 2. Root cause

**Primary (dominant, ~95 % of pathology):** parcels seeded directly over the
open orifice at `t = 0`, combined with `des_interp_scheme = 'LINEAR_HAT'`.

- `IC_2` seeds the packed bed at `ic_ep_s = 0.57` over the entire box footprint
  from `y = 0` to `y = H_static = 32 mm`. The 4-mm orifice at `y = 0` sits
  directly below seeded parcels.
- LINEAR_HAT drag projection smears parcel mass across neighbouring cells with
  ~2-cell support. Within the first few ms, parcels in the center of IC_2 fall
  through the orifice, and LINEAR_HAT projects their mass into the taper
  cut-cells at `y ≈ -1 mm`, `r ≈ 1.4 mm`.
- EP_G in those cells drops to **0.33** — below close-packing, which is
  physically impossible from real parcel density and can only be a projection
  artifact.
- Mass balance then forces gas velocity in those cells from the expected
  13 m/s (plug-flow ratio `U_in · (R_stub/R(y))²`) to ~100 m/s. The semi-implicit
  SIMPLE loop, running at CFL ~15 due to uncapped `dt_max`, under-converges the
  momentum–drag coupling and produces instantaneous overshoots to the 750 m/s
  MAX_INLET_VEL_FAC clip ceiling (and above — Vg = 9380 m/s was observed).

**Secondary (amplifier):** `dt_max` is set to 1e-3 s in the current deck, which
is 30× the CFL-safe value for `U_in = 30 m/s` and `dx = 2 mm`. At CFL ~15 the
Picard iterations saturate before converging a sharp mass-flux imbalance, so
each spike gets amplified ~7× instead of being damped per step.

### What was ruled out

- **Stub/taper geometry alone:** without parcels (smoke run `cases/V_3d/day10_r2_taper_smoke/`),
  centerline velocity ratio is 0.998 of expected. The taper geometry is correct.
- **MI boundary condition:** stub interior shows uniform `v_y ≈ 7.5 m/s`; the MI
  is 19 mm below the trouble zone.
- **Step inlet:** `bc_v_g = 7.5 m/s` is well sub-sonic; no shock. SIMPLE handles
  it in 6-11 NIT within 3 steps.
- **F_AT / master-cell artifact:** 9,512 cut cells / 51,432 blocked / 47,568
  fluid is self-consistent with the SDF; no cell-flag warnings.
- **Welch resolution:** `nperseg = 1024` on 4,500-sample signal resolves bins at
  16.6, 17.6, …, 24.4, 25.4 Hz. The 0.98 Hz peak is lack of band energy, not
  parameter choice.

## 3. The two fixes

### Fix A — Quadrant IC_2 + gas-only orifice column (primary, correctness)

Replace the single rectangular IC_2 with four strips plus a gas-only IC_5 column
over the orifice, matching NETL's shipped `tutorials/pic/spouted_bed_3d/spouted_bed_pic_3d.mfx`.

- Exclude `|x| ≤ R_i + dx = 4 mm` and `|z| ≤ R_i + dx = 4 mm` from the parcel-seeded
  region; keep IC_2a-d as four packed-bed quadrant strips (ep_s = 0.57).
- Add IC_5: `|x|,|z| ≤ R_i + dx`, `0 ≤ y ≤ H_static`, gas-only (`ic_ep_g = 1`,
  `ic_u_g = ic_v_g = ic_w_g = 0`).
- Expected parcel count drops from ~32,000 to ~30,600 (−4 %), within Day-10
  validation band.
- Files to edit: `src/make_case.py` `build_V_3d_params()` (IC_2 construction);
  renders into `cases/V_3d/day10/V_3d_day10.mfx:302-319`.

### Fix B — Tighten `dt_max` from 1e-3 to 5e-5 s (secondary, robustness)

The current `dt_max = 1e-3 s` is 30× the CFL-safe value for the inlet. Tightening
to `dt_max = 5e-5 s` (CFL ~1 at the paper-expected 13 m/s plug-flow spout) prevents
the semi-implicit solver from amplifying any future transient 7× before the Picard
iterations can correct it.

- Files to edit: `cases/V_3d/template.mfx.j2:43-45` and `cases/V_3d/template_dem.mfx.j2:31`.
- Alternative if post-A CFL still > 1: drop further to 1e-5 s.

**Fix A is strictly required (correctness).** Fix B alone would reduce the spike
magnitude but not remove the source — the solver would converge on wrong physics.
Both fixes should land together.

## 4. Impact on currently-running Day-11 and Day-13

Both AWS Batch jobs inherit the same deck structure as Day-10 (identical IC_2,
same `LINEAR_HAT`, same uncapped `dt_max`), **and** the pathology directly
contaminates their primary deliverables.

### Day-11 (CFD-CGP, job `1c4845c5-f014-48e8-ac43-a507853ec9ea`)

Headline deliverables (`notes/day11_summary.md`, `src/post/residence.py`):

| Deliverable | Status | Why |
|---|---|---|
| Spout residence histogram | **fatally contaminated** | Region defined as `r ≤ 2 mm, 0 ≤ y ≤ 32 mm` — exactly where Day-10's jet dies at y=37 mm instead of reaching the full spout |
| Fountain residence histogram | **fatally contaminated** | Fountain = `r ≤ 25 mm, 32 mm < y ≤ 200 mm`. Jet can't eject parcels above the bed → fountain empty by construction. Smoking gun: `results/fig_d11_V_3d_dem_residence.csv:3` already shows `fountain n_particles = 0` from the smoke run |
| Circulation cycle time | **fatally contaminated** | Derived from the two residence percentiles |
| Centerline v_g vs height | contaminated | Same artifact as Day-10 |
| Inlet PSD peak | likely contaminated | Same drag/EP_G coupling as Day-10 |
| Bulk ΔP, inventory | partially unaffected | Weight-dominated, but supporting trace only |

**Verdict:** Letting the run finish produces an empty fountain histogram. The
headline deliverable cannot be recovered by longer averaging.

### Day-13 (12-case OFAT sweep, array job `6266f08d-e583-41d6-9f1a-28569afed712`)

The sweep varies U_in, gas composition, ρ_p, drag law, en coefficient at fixed
geometry (`R_i = 2 mm` everywhere per `cases/D_sweep13/manifest.json`). Despite
the "D_sweep" directory name, no geometry axis is swept.

Headline deliverables (`notes/day13_summary.md:106-113`, `src/post/plot_V_3d_day13.py`):

| Axis | Status | Why |
|---|---|---|
| Gas composition (hot Ar, Ar/H₂) | **fatally confounded** | ρ_g drops 4× → LINEAR_HAT contraction amplifies more at same mass flux. Spurious ρ_g dependence is dimensionally identical to the real physics being measured. This is the flagship axis per `notes/day13_summary.md:150-155` |
| Inlet-pressure PSD peak | contaminated across all cases | All lock onto 0.98 Hz artifact bin |
| Mid-bed centerline ⟨v_g⟩ | contaminated | Mid-bed `y ≈ 16 mm` is inside the dead-spout region |
| ⟨ΔP_bed⟩ absolute values | partially contaminated | Weight-dominated, but shifted by artifact |
| Density axis trend | **least contaminated** | ρ_p scales bed weight, orthogonal to artifact |
| Drag axis trend | qualitatively usable | Both drag laws share LINEAR_HAT; artifact common-mode cancels |
| Restitution axis | contaminated | Real effect drowned by artifact |

**Verdict:** The gas-axis fan plot would look like data but encode a ρ_g artifact
that mimics the physics. Worse than no data because it is hard to spot as wrong
without this investigation in hand.

## 5. Recommendation

**Cancel both AWS Batch jobs now.** Apply Fix A and Fix B together to both
templates, then resubmit Day-10, Day-11, and Day-13 from the fixed baseline.

Files to change (one commit, four files):

1. `src/make_case.py` — `build_V_3d_params()`, split IC_2 into quadrant strips
   + add IC_5 gas-only column.
2. `cases/V_3d/template.mfx.j2` — tighten `dt_max` to 5e-5 s.
3. `cases/V_3d/template_dem.mfx.j2` — tighten `dt_max` to 5e-5 s.
4. `tests/test_make_case.py` — assert four-quadrant IC_2 layout + IC_5 gas-only.

Rerun order:

- **Day-10 PIC** first (~2 h on AWS, `tstop = 7 s` so the Welch window covers
  df ≤ 0.5 Hz for S1/S2 band discrimination). Confirm max CFL < 2, max |V_gas| < 50 m/s,
  centerline matches S1 shape within tolerance, PSD peak lands in 10-25 Hz band.
- **Day-11 CGP** only after Day-10 is validated.
- **Day-13 sweep** only after Day-11 confirms the residence histogram is non-empty.

## 6. Additional post-processing follow-ups (not blockers)

Lower-priority items that improve the next run but are not required for the fix:

- **Shrink inlet pressure monitor patch** from ±R_i = ±2 mm (4×4 mm² square,
  21 % outside the orifice circle) to ±R_i/√2 ≈ ±1.41 mm, OR relocate to mid-stub
  plenum. Template line `cases/V_3d/template.mfx.j2:249-258`.
- **Welch tuning in `src/post/plot_V_3d_day10.py:164`**: pass explicit
  `nperseg = 2048, noverlap = 1024, detrend = 'linear'` once `tstop ≥ 7 s`.
- **Audit the digitized S1 Fig. 1 trace** (`data/digitized/s1_fig1_centerline_vg.csv`)
  — currently flagged as placeholder with "UNRESOLVED provenance" (is it an
  experimental PIV/LDV trace or a TFM simulation output?). Must be settled
  before any quantitative claim.

## 7. References

- **Day-10 run artifacts:** `cases/V_3d/day10/aws_out/` (120 BACKGROUND VTUs,
  3 monitor CSVs, run.log, V_3D_DAY10.TXT).
- **Post-processing output:** `results/fig_d10_V_3d_*.{png,csv}`.
- **AWS job records:** `results/run_log.csv` rows tagged `aws-batch,day10`,
  `aws-batch,day11`, `aws-batch,day13`.
- **Related commits:** `28cf021` (stub→orifice taper), `6e3e2c4` (stub radius
  decouple), `3346697` (Day-10 PIC baseline promotion).
- **NETL reference:** `/Users/nelsonpereira/mamba/envs/mfix-26.1.2/share/mfix/templates/tutorials/pic/spouted_bed_3d/spouted_bed_pic_3d.mfx`
  — canonical quadrant IC pattern.

## 8. Investigation audit trail

Subagents launched during this review (available via SendMessage for follow-up):

- `a0a7b218e6117b097` — root-cause investigation (identified LINEAR_HAT +
  parcel-over-orifice coupling).
- `aec5e3a60c1ba1296` — fix evaluation (confirmed A + B resolves PSD, CFL, ΔP
  issues; identified Welch nperseg + monitor patch follow-ups).
- `a98c94e787028c8fb` — Day-11/13 structural inheritance (confirmed identical
  IC_2, LINEAR_HAT, dt_max across all three days).
- `ac1385435bf3fac2c` — deliverable-relevance (confirmed Day-11 residence and
  Day-13 gas-axis are fatally contaminated by the artifact).

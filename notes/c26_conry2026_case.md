# Case C26 — Conry et al. 2026 ORNL benchtop rig: design note (CLAUDE.md rule 9)

Written 2026-10-16, BEFORE any code for this case. It explains what the new
case is, why it needs its own blocks instead of reusing Case V, and every place
it deviates on purpose from `.claude/skills/mfix-spouted-bed-recipe`.

## Source

Conry B., Dal Forno Chuahy F., Oyedeji O., Finney C.E.A., Heldt R.,
López-Honorato E. (2026). "Development of a coupled experimental–computational
approach for engineering optimization of spout-fluidized bed particle coating
systems." *Nucl. Eng. Des.* 456: 115009. DOI 10.1016/j.nucengdes.2026.115009.
Cited below as **C26**. The rig is also described in Conry et al. (2025),
*Chem. Eng. Sci.* 306: 121197 (**C25**).

C26 runs MFiX 23.2 CFD–CGDEM and benchmarks it against PIV of a room-temperature
3D glass spouted bed. The project's MFiX is 26.1.2. Every keyword below was
checked against the installed 26.1.2 source (see "Keyword verification").

## Why a separate case and not a Case V variant

Case V reproduces the UTK/ORNL cold mockup (S1 = ORNL/TM-2006/520, S2 = Zhou et
al. 2005). C26 is a different rig with different everything:

| | Case V | C26 | C26 source |
|---|---|---|---|
| Column ID | 50 mm | 50.4 mm in the model (24 × 2.1 mm); text 50 mm; C25 50.8 mm | §2.2.1, §2.1; C25 Fig. 1 |
| Cone | 60° | 60° | §2.1 |
| Inlet | 4 mm orifice, stub pipe, MI | single 2.38 mm (3/32 in) nozzle, **gas point source** | §2.1, §2.2.1 |
| Mesh | 1 mm uniform | 2.1 mm uniform, 24 cells across D, 100 vertical | §2.2.1 |
| Gas | air 25 °C | 100 % Ar, 298 K | §2.2.2 |
| Flow | 30 m/s jet | 10/12/14 slpm (YSZ), 22/26/30 slpm (SS) | §2.1, Table 2 |
| Solids | 500 µm ZrO₂, 6050, mono | YSZ 6000 kg/m³, measured PSD mean 411.6 µm; SS 8460, 784.9 µm | Table 1, §2.2.2 |
| Inventory | 54.5 g | 60 g | §2.2.2 |
| Drag | SYAM_OBRIEN | TGS | §2.2.1 |
| e_n p-p / p-w | 0.9 / 0.9 | 0.5 / 0.5 | Table 1 |
| μ p-p / p-w | 0.3 / 0.3 | 0.1 / 0.1 (p-w swept 0.1–0.7) | Table 1, Table 2 |
| k_n, k_t/k_n, η_t/η_n | 1000, 2/7, 0.5 | 1000, 0.285, 0.5 | Table 1 |
| CGP weight | 10 (5 in later runs) | 10 baseline; 5/3/1; SS at 1 | Table 2 |
| Validation | centerline v_g (S1 Fig. 1), ΔP PSD (S2) | side-view PIV velocity maps (CC/SSI/MSE), wall v(t) | §3, Table 3 |

Mixing these into `bed_V` / `case_V` would silently change Case V. New blocks
are added instead (`bed_C26`, `spouted_bed_C26_3d`, `spouted_bed_C26_3d_dem`,
`case_C26_ysz`, `triso_coating_layers`, `ar_25C`, `spouted_bed_C26_3d_numerics`,
`params/sweep_C26.yaml`). Case V files are untouched.

## Deliberate deviations from the spouted-bed recipe

The recipe was written for Case V, and its own decision table says to deviate
only with a pathology-review-grade justification. Here are the justifications.

### Recipe Rules 1, 3, 4, 5 (4 cells across the orifice; straight stub; MI over the full box face) → replaced by a gas point source

- C26's mesh has D_nozzle / dx = 2.38 / 2.1 = **1.13 cells**. On a stub-pipe MI
  that is the Day-10 R1 failure mode (R_i/dx = 1 → 0 fluid cells in the stub,
  run_log 2026-10-09). Rule 1's floor of 4 cells would need dx ≤ 0.595 mm
  (≈ 85 × 360 × 85 = 2.6 M cells, ~35× the C26 mesh). That is no longer the
  model the paper benchmarked.
- C26 avoids the problem the same way NETL's
  `tutorials/dem/loop_seal_pseudo_2d` does: a **gas point source**
  (`ps_*` keywords) "centered with the reactor and placed at the bottom of the
  cone" (§2.2.1). The authors chose dx ≈ nozzle diameter because "the cell
  size near the nozzle need to be similar to the nozzle size, otherwise when
  the flow mass is added to the domain, calculated velocities will not be
  correct" (§2.2.1, p.3).
- In 26.1.2 the point source adds `ps_massflow_g` spread over the fluid volume
  of its cells (`set_ps.f`). It adds momentum `ps_massflow_g × |v| × v̂`, where
  |v| is the magnitude of (`ps_u_g`, `ps_v_g`, `ps_w_g`) (`source_v_g.f:984-987`,
  `set_ps.f:CALC_PS_VEL_MAG`). Since 2026-10-19 the source is mass-only by
  default (`ps_momentum: false`, `ps_v_g = 0`; user decision, see "AWS
  production runs: mass-only point source"). With `--ps-momentum`, `ps_v_g`
  = nozzle exit velocity Q / (π D_n²/4) and the jet carries ṁ·U_n; that was
  the design until 2026-10-19 and is kept for A/Bs.
- **Axis placement.** 24 cells across D puts the vessel axis on a cell face.
  To put it at a cell centre, the box uses an odd cell count (27 = 24 across D
  plus a 1.5-cell margin each side). The point source is then one cell at
  (0, dy/2, 0). One cell matches C26's "cell ≈ nozzle" argument. A 2 × 2
  footprint would cut the cell velocity ~4×.
- Gas point sources work with DEM/CGP solids in 26.1.2. Only *solids* point
  sources are rejected for DEM/CGP (`check_point_sources.f:443-451`).

### Recipe Rule 6 (4-quadrant IC with a central gas column) → one loose bed region

- R5 (`notes/day10_r5_ic_review.md`) found that the quadrant pattern copies
  NETL's draft-tube bore, which this geometry does not have. With a closed cone
  apex and no stub, there is also no orifice for seeded particles to fall
  through. That was the original Fix-A motivation.
- C26 §2.2.2: "Particles were initialized at the bottom of the reactor and
  followed the particle size distribution". We seed one region over the full
  cross-section from y = 0 to H_ic at a loose ε_s,ic = 0.30 (assumed), as in
  R5's loose-IC approach. The bed settles under gravity while the point source
  flows.

### Recipe Rule 12 (SYAM_OBRIEN) → TGS

C26 §2.2.1 uses TGS (Tenneti, Garg & Subramaniam 2011). It is a valid
26.1.2 `drag_type` value for DES (`init_namelist.f` valid list;
`des/drag_gp_des.f:221`). `drag_c1/drag_d1` are Syamlal–O'Brien parameters and
are not emitted.

### Recipe Rule 11 (DEM closure) → C26 Table 1 values

μ = 0.1 (not 0.3), e_n = 0.5 (not 0.9), k_t/k_n = 0.285 (C26 rounding of 2/7).

### Recipe Rule 9 (dt_max ≤ 0.5·dx/U_in) → paper step, capped by Rule 9

dt_max = min(1.0e-5 s, 0.5 · dx / U_nozzle), from `params/numerics.yaml`.
C26 p.9 (§3.2.1): "a single time step was 0.01 ms". That is the binding value
for every case in `params/sweep_C26.yaml` (U_nozzle ≤ 52.4 m/s, so the Rule-9
bound is ≥ 2.0e-5 s). Rule 9 would bind only above ~28 slpm. This deliberately does not copy the V_3d decks' 1e-3 s ceiling,
which R5 flagged as a Rule-9 violation. (Changed 2026-10-17. The first c01
smoke used the Rule-9 value alone, 2.8e-5 s.)

### Unchanged from the recipe

Uniform Cartesian mesh (Rule 2), IC 1 whole-domain gas (Rule 7), abort-off
seeding (Rule 8), LINEAR_HAT + mean fields (Rule 13), params in yaml (Rule 14),
PO at y = H_dom with the STL top cap 5 mm above it (Rule 19), inward STL
normals (Rule 20), and the R5 overlapping-seam SDF union (no membrane).

## Geometry decisions the paper does not settle (all marked "assumed" in yaml)

1. **Throat diameter.** Unknown. We assume the cone narrows to the nozzle
   diameter (D_throat = 2.38 mm) and closes with a flat cap at y = 0. Gas
   enters only through the point source.
2. **Cone angle.** 60° included, per the text (half-angle 30°). Cone height
   = (R_c − R_t)/tan 30° = 41.6 mm.
3. **Bed-height inconsistency.** The C25/C26 Fig. 1 annotation "~25 mm" is the
   *static bed height* (C25 caption: "dimensions of the cone width, cone
   angle, and bed height"). 60 g of YSZ at ε_s = 0.57 in a 60° cone gives
   **~35 mm**. ~25 mm is only consistent with a ~90° cone. Every manifest
   logs `H_static_m` and the 90° cross-check `H_static_if_90deg_m`. The first
   smoke settled at 34.5 mm (2026-10-17), matching the 60° estimate.
4. **Column height / outlet.** Domain top = 100 × 2.1 mm = 0.210 m (§2.2.1),
   with a pressure outlet there.
5. **Below-cap margin.** 4 blocked cell rows under the y = 0 cap. That keeps
   the cap off the box face (cf. Rule 19) and makes jmax = 104 divisible by
   2/4/8 for DMP.
6. **slpm reference state.** 25 °C, 1 atm (Alicat default "standard"; the MFC
   is an Alicat per C25 §2.2). That equals the 298 K gas state, so the
   standard and actual volumetric flows coincide.

## Validation target (what "done" means for this case)

C26 reports no gas velocity or pressure. Its quantitative targets are:

- **Table 3**, time-averaged side-view velocity-magnitude maps (spout region),
  YSZ cases 1–3: CC = 0.6817 / 0.8170 / 0.8286; SSI = 0.6086 / 0.7687 / 0.7657;
  MSE = 0.0068 / 0.0042 / 0.0049 m²/s².
- Wall v-component time series (Fig. 1d and the wall-velocity figures).

Both come from PIV on side-view frames. A fair comparison renders simulated
particles into the camera plane and runs the same PIV. That post-processor is
**not** part of this change (see open items). The deck writes the data it
needs: particle VTK with ids every 2 ms for 300 frames over the last 0.6 s
(`vtk_t_start`/`vtk_t_end`, verified keywords), matching §2.2.2.

## Keyword verification (CLAUDE.md rule 1)

Checked 2026-10-16 against the installed source at
`~/mamba/envs/mfix-26.1.2/share/mfix/src/model/`:

- `ps_x_w/x_e/y_s/y_n/z_b/z_t`, `ps_massflow_g`, `ps_u_g/v_g/w_g`:
  `init_namelist.f` 3931-4022; checks in `check_data/check_point_sources.f`;
  usage pattern from `templates/tutorials/dem/loop_seal_pseudo_2d`.
- `drag_type = 'TGS'`: `init_namelist.f` valid-value list;
  `check_data/check_solids_common_all.f:270`.
- `ic_psd_type ('NORMAL')`, `ic_psd_mean_dp`, `ic_psd_stdev`, `ic_psd_min_dp`,
  `ic_psd_max_dp`: `init_namelist.f:2215-2275`. In CGP mode MFiX scales them
  by `cgp_stat_wt^(1/3)` internally (`check_solids_common_discrete.f:190-194`),
  so the deck carries **real** particle sizes.
- `vtk_t_start`, `vtk_t_end`: `cartesian_grid/cartesian_grid_init_namelist.f:976-992`.
- `monitor_pmass`, `monitor_pcount` (type-101 Monitor 3, added 2026-10-17):
  `init_namelist.f` 4715/4722; `check_output_control.f` 888-892; usage from
  `tutorials/dem/hopper_3d`.
- `des_interp_width` (emitted only with `--interp-scheme SQUARE_DPVM`,
  2026-10-19): `des/des_init_namelist.f` 1021-1040.
- All other keywords are reused unchanged from `cases/V_3d/template_dem.mfx.j2`.
- Cross-check: all 135 keyword names in the rendered 2026-10-16 decks were
  grepped against every `*init_namelist*.f` in the installed source; none
  are missing. `tests/test_make_case_c26.py::test_deck_keywords_are_verified`
  renders c01/c04/c07/c09 and keeps each deck within the V_3d-template set
  plus the list above (default decks, so not `des_interp_width`).

## First c01 smoke and full-text review (2026-10-17)

- **Point-source placement bug.** The first smoke (AWS, 0.3 s, run_log
  2026-10-16) gave the PS as a zero-width point at the axis-cell centre.
  `get_ps.f`/`calc_cell.f` round each coordinate to a face, and a cell centre
  is the tie. The source landed in the (−,−) corner sliver cut cell
  (PS volume 3.27e-10 m³). It drove a 150–170 m/s wall-hugging jet that lost
  most of its momentum within ~10 mm. The axis cell (7.0e-9 m³, 76 % fluid)
  got ~3 m/s. Fix: the PS region is now the axis cell's faces.
- **Full text (experimental-computational-modeling-approach.pdf)**:
  - dt = 0.01 ms (p.9).
  - ≥ 3.0 s startup before averaging (p.3).
  - 10 slpm is "near minimum spouting", with ~10 spouting events per 0.5 s.
  - A single point source "centered with the reactor and placed at the
    bottom of the cone", fed by one MFC via a T-junction (Fig. 1b), so no
    auxiliary fluidizing gas.
  - The simulated cone closes to a point (Fig. 8e–f).
  - PSD variance 149.1 µm² (σ = 12.2 µm) confirmed.
  - Fig. 1b's "~25 mm (1 in.)" is a schematic annotation (C25 Fig. 1 draws
    the same bracket for 15 g SiC). The first smoke settled at 34.5 mm,
    matching the 60° estimate.
  - Fig. 8e's simulated cone looks wider than 60° by eye. Measure it from a
    high-resolution extraction before changing geometry; the text says 60°.
- Monitor 3 is now a type-101 particle mass/count sum (the type-11
  `monitor_ep_s` monitor was dropped for DEM, Warning 4455).

## Run cost: DMP partition and STL facet density (2026-10-17)

The second smoke (`c01_smoke_r2`, 3×1×3 = 9 ranks on 16 vCPU) ran at
0.52 s/step. CGDEM took 0.40 s of that and the fluid solve 0.09 s, which puts
the 3 s baseline at ~44 h. Two ways to cut that were tested. All runs are in
`results/run_log.csv`.

1. **Non-uniform DMP (`gridmap.dat`) on more ranks: slower, not adopted.**
   - Support was added: `--gridmap-i/-j/-k` in `src/run_C26_3d.py`, and
     `write_gridmap` plus the `gridmap=` argument in `src/make_case_c26.py`.
     The run directory gets MFiX's `gridmap.dat` (`dmp_modules/gridmap_mod.f`),
     and nodes are taken from the gridmap lengths.
   - Check run `c01_gridmap_check`: 0.03 s, i = k = [9,3,3,3,9], 25 ranks,
     32 vCPU.
   - Physics matched r2: inlet pressure within 10 Pa, centerline v_g within
     2 %, identical bed heights. So MFiX handles uneven i/k slabs despite its
     "j-direction only" source comment.
   - CGDEM got *slower*: 0.52 s/step against 0.40 s in r2. The busiest rank's
     parcel count fell 8× without any gain, so particle load per rank is not
     the bottleneck.
   - The support is kept. The c01 baseline used uniform 3×1×3; since
     2026-10-19 the default is the 4-rank gridmap i = k = [13,14] (see "AWS
     run cost").
2. **Coarser STL: 3.5× faster DEM, adopted.**
   - `dx_sample_fraction` (`params/geometry.yaml`) changed from 0.25 to 0.5.
     The marching-cubes sampling step goes from 0.525 to 1.05 mm, and the
     STL from 311,640 to 77,088 facets.
   - The STL does **not** set the fluid mesh: 27×104×27 at dx = 2.1 mm, dt
     and CFL are unchanged. It only sets how the wall cuts cells, and the
     facets parcels collide with. DEM wall contact is checked against the
     facets binned in each cell, so cut cells went from ~32 facets each to ~8.
   - Timing pair `c01_stl_q_timing` vs `c01_stl_h_timing`: local, 4 ranks,
     0.01 s, decks identical apart from the STL.
     - CGDEM 0.140 → 0.040 s/step; fluid 0.100 → 0.090 s/step.
     - Standard/cut/fluid/blocked cell counts identical
       (34615/8949/43564/32252).
     - PS volume identical (7.013e-9 m³).
     - Seeded 27,782 vs 27,787 parcels. The gap is fixed from step 1, so it
       comes from seeding, not loss.
   - Validation smoke `c01_smoke_r3`: local, 0.3 s, same deck as r2 except
     the STL. `check_smoke.py` gave GO with 0 FAIL / 0 FLAG. Every step-2
     criterion matched r2 over t = 0.1–0.3 s:

     | metric (t = 0.1–0.3 s) | r2 (dx/4) | r3 (dx/2) |
     |---|---|---|
     | upward solids flux r < 3.15 mm, y = 25.2 mm (≥ 5 ṁ_g) | 38.8 ṁ_g | 40.9 ṁ_g |
     | flux centroid (≤ 1.05 mm) | 0.08 mm | 0.05 mm |
     | fastest row-0 cell | axis, 1.12–1.24 U | axis, 1.15–1.25 U |
     | centerline v_g, y = 17–19 mm | 13.0–14.7 m/s | 11.9–15.3 m/s |
     | ΔP / W/A, PSD peak | 1.79, 15 Hz | 1.69, 14.9 Hz |
     | solids inventory | constant | constant |

   - Wall roughness: at dx/2 a facet edge is ~1 mm against 0.41 mm
     particles. Checked 2026-10-19 after c01's wall velocity came out off:
     r2 (dx/4) and r3 (dx/2) give 1.61e-2 and 1.79e-2 m/s, within each run's
     own spread (see "Not the STL" below).
   - `--stl-sample-fraction` in `src/run_C26_3d.py` overrides the yaml value
     for such tests.
   - c01–c03 were regenerated with the new STL (3×1×3, 3 s).
   - The dx/4 c01 baseline (job f476de9a…) was cancelled after 1.4 h and
     resubmitted as `c26-c01-baseline-stl05`
     (57ca0039-97ea-4722-8932-2bf48430a8c0). It should take ~19 h if the
     local speedup carries over to AWS.

## c01 baseline result and step-3 spout checks (2026-10-19)

**The run.** c01 reached t = 3.000 s in 22.3 h on 16 vCPU (job 57ca0039),
which is 0.264 s/step:
- CGDEM took 0.15–0.19 s/step, a 2.3× gain over the dx/4 STL on AWS (3.5×
  locally). Fluid took 0.09 s/step.
- dt stayed at 1e-5 s; 12 % of steps were briefly reduced, to as low as
  6.6e-6 s.
- No parcels were lost: 27,782 parcels and 60.95 g at t = 3 s.

**Health check** (`check_smoke.py`, t = 0.1–3.0 s): GO.
- ΔP = 1.66 W/A, std 0.7 W/A.
- S4 exit-layer EP_G was 0.891, below 0.9. S4 is now report-only for C26
  (`--exit-voidage-info`).
  - The decision was made after seeing this value, and that is stated in
    SKILL.md §3.6.
  - In C26 the exit layer is the jet base at a closed apex, not a throat.
  - The value showed no trend over 3 s: 0.879–0.898 per 0.5 s window.
  - The exit v_g ratio was 0.97, so the inlet is not choked.

**Step-3 criteria.** Thresholds were written before evaluation
(2026-10-19). The window is the last 0.6 s, t = [2.4, 3.0] s, 300 PIV
frames; the script is noted in run_log.

| criterion | threshold | c01 | |
|---|---|---|---|
| A fountain: ≥ 20 parcels above H_static + 10 mm | ≥ 50 % of frames | 100 % (median 499) | PASS |
| B ΔP PSD peak, t = 2.4–3.0 s / 0.5–3.0 s | 10–20 Hz | 18.3 / 18.4 Hz | PASS |
| C1 Fig. 7 axis v_p (line probe, axis-centred 2.1 mm bin) at peak height | 0.26–1.06 m/s | 1.75 m/s at y = 10 mm (1.22 at 20, 0.82 at 30) | **FAIL** |
| C2 outermost bin \|v_p\| | < 0.1 × axis | 0.007–0.020 m/s | PASS |
| D wall v_p (downward +, parcels within d_cg of the wall, y < H_static) | 0.3e-3–10e-3 m/s | 1.11e-2 m/s (2·d_cg band: 1.10e-2) | **FAIL** |

- B agrees with C26's "~10 spouting events every 0.5 s" (20 Hz).
- **C1 binning.** The first evaluation put a bin edge on the axis by mistake.
  It gave 0.97 m/s peak, which would have PASSED. Centring the bins on the
  axis, as the mesh is, changed the result to FAIL. The corrected value
  stands.
- **D in context.** The value is 3.4× C26's model (3.31e-3) and 12× their
  experiment (0.933e-3). It is close to their W = 1 model (13.9e-3,
  Fig. 8c).
- **Not the STL.** Over t = 0.1–0.3 s, the r2 (dx/4) and r3 (dx/2) smokes
  give 1.61e-2 and 1.79e-2 m/s. Each run's own half-to-half spread is
  larger than that gap.

**Reading.** Both failures point the same way. Our bed circulates ~2–3×
faster than C26's model: a faster spout core means faster annulus return
flow at the wall. Candidate causes at the time (1 was then A/B-tested below
and is being run as the 3 s c01–c03 mass-only runs; 3 is reviewed under
"Gas-fraction floor"):
1. **Point-source momentum.** Our source injects gas at `ps_v_g =
   U_nozzle`, adding a momentum flux ṁU at the apex.
   - C26 §2.2.1 justifies the 2.1 mm cell because, with a point source,
     "when the flow mass is added to the domain, calculated velocities will
     not be correct" unless the cell is nozzle-sized.
   - That suggests a mass-only source (`ps_v_g = 0`), where the jet speed
     comes from continuity through the nozzle-sized cell.
   - Test: c01 with `ps_v_g = 0`, compared on C1/D.
2. **Probe definitions.** C26 states neither the Fig. 7 probe height nor
   how it averages, nor how its "line-averaged" wall velocity is formed. A
   wider or Eulerian average would lower C1's peak. That would not explain
   D.
3. **MFiX 23.2 vs 26.1.2.** Changes in cut-cell wall contact are possible.
   Not testable without the old version.

### Point-source momentum A/B (c01_smoke_r4_massonly, 2026-10-19)

- **Paper.** C26 never states a point-source velocity (subagent review of
  the full PDF).
- **MFiX defaults.** The GUI defaults ps_u_g/ps_v_g/ps_w_g to 0 for a
  mass-flow PS (`mfixgui/pss.py:964-980`). The solver then zeroes undefined
  components, which gives no momentum source
  (`check_point_sources.f:277-281`, `source_v_g.f:964-987`).
- **Reading.** So a GUI user following §2.2.1 gets a mass-only source. The
  p.3 cell-size argument also reads most naturally that way: a 2.1 mm cell
  face gives 37.8 m/s by continuity, against the nozzle's 37.4 m/s. Not
  proven.
- **PIV method ruled out for D.** C26 Fig. 8c gives PIV/parcel ratios of
  0.83–1.49 (W = 10: 2.9 vs 3.2e-3 m/s). PIV-like variants of our wall
  average (front sector, a Fig. 1c-like 4 mm ROI, a gridded front face) give
  0.96–1.7e-2 m/s.
- **A/B setup.** Identical to r3 except ps_v_g = 0. Rule (pre-registered):
  - substantial if wall and peak axis v_p both fall ≥ 30 %;
  - ruled out if both change < 15 %;
  - otherwise inconclusive.

| t = 0.1–0.3 s | r3 (ps_v_g = U_nozzle) | r4 (mass only) |
|---|---|---|
| peak axis v_p (y = 10 mm) | 1.56 m/s | 0.42 m/s (−73 %) |
| wall v_p (downward) | 1.79e-2 m/s | 1.93e-2 m/s (+8 %) |
| axis-cell v_g / U_nozzle | 1.15–1.25 | 0.45–0.53 |
| axis gas momentum flux at y = 3 mm | ~0.8–1.1 ṁU | 0.03–0.08 ṁU |
| centerline v_g mid-bed | 11.9–15.3 m/s | 2.5–3.3 m/s |
| ΔP / W/A, PSD peak | 1.69, 14.9 Hz | 2.70, 10 Hz |

- **Verdict: inconclusive** under the rule.
- PS momentum sets the spout core. With mass only, the axis v_p is close to
  C26 Fig. 7's ~0.53 m/s, and source gas leaves the apex cell partly
  sideways.
- ~~PS momentum does not explain the 3× wall velocity.~~ Retracted
  2026-10-19. The wall number cannot decide this. A circulation check
  (parcel flux through y = 20 mm, same frames) shows:
  - **r3:** spout net up +4 500 /s against annulus net down −10 200 /s. The
    bed is still compacting from the loose IC, so the early wall value
    includes settling, not just circulation. c01 at 2.4–3.0 s is balanced
    (+5 180 / −5 285 /s).
  - **r4:** spout net up only 14 % below r3 (+3 900 /s): the jet is slower
    but carries nearly as many parcels. Annulus net down is only −2 560 /s,
    and there are 24 fountain parcels against r3's 360.
  - **r4 wall band:** mean |v_y| is 6.0e-2 m/s, 2× r3. The frame-to-frame sd
    of the band mean (5.9e-2) is 3× the mean (1.93e-2). The "+8 %" is a
    small net of large up/down motions. It is not evidence of faster
    circulation.
  - **C26 also reports pulsation.** Its model's wall velocity pulses between
    positive and negative at 10 slpm; the experiment does not (p.7, Fig. 6).
    C26 links this to bubble formation, "spout incoherence" (p.10). It
    still calls the regime spouting (~10 events / 0.5 s) and says the wall
    pulsations are small relative to the core.
  - **c01 at 2.4–3.0 s looks like that:** wall mean upward in 8 % of frames,
    0.05 s moving average −5.5e-3 to +2.4e-2 m/s, PSD 18 Hz.
  - **r4 pulses far harder:** 46 % of frames upward, moving average
    −4.5e-2 to +8.2e-2 m/s, 24 fountain parcels. It is also still in the
    startup window.
  - A steady 3 s mass-only run scored at 2.4–3.0 s is the only valid test.
- r4 fails `check_smoke` S3 (v ratio 0.39/0.24) and the mid-bed v_g ≥ 5 m/s
  supporting check. Those were written for a momentum source; the mass-only
  criteria were set afterwards (see "Mass-only smoke criteria").
- `ps_momentum` stayed `true` at the time. It became `false` (mass-only) on
  2026-10-19 by user decision; see "AWS production runs: mass-only point
  source".
- The 0.1–0.3 s window is early. c01's wall v_p fell from ~1.8e-2 (0.1–0.3 s)
  to 1.1e-2 (2.4–3.0 s), so only a 3 s mass-only run can be scored against
  Figs. 6–7.

### Gas-fraction floor (`des_epg_clip`) and interpolation (2026-10-19)

- **MFiX version.** A subagent compared wall contact in MFiX 23.2 and 26.1.2.
  It used the 23.2 keyword docs, the 23.2–26.1 release notes, and a ~21.x copy
  of the wall-contact routine. Every wall and contact keyword and default is
  unchanged. The only relevant code change, per-particle wall damping in 25.1,
  is about ±4 % at our 3 % size spread. Wall contact is unlikely to explain D.
  The 23.2 source itself was not available, so this is not a direct diff.
- **The floor is active in most of the bed.** `des_epg_clip = 0.42` is our
  assumption (copied from Case V); MFiX's default is no floor. In c01 (13
  BACKGROUND frames, t = 2.4–3.0 s) EP_G sits at 0.42 in 61 % of cells within
  one cell of the wall and 63 % of annulus cells, but 2 % in the spout.
- **No floor cannot run on this mesh.** Error 1100 (EP_G ≤ 0,
  `des/calc_epg_des.f`) hits sliver cut cells (VOL 1.0–1.8 % of a full cell,
  no parcel centred in them) that receive neighbouring parcel volume. One
  parcel is 3.7e-10 m³, bigger than the cell.
  - LINEAR_HAT, no floor: fails at t = 0.
  - SQUARE_DPVM, `des_interp_width` = 2.1 mm (the C26 cell size; width
    assumed, not reported), no floor: fails at t = 0.0048 s as the bed
    settles. A 3 ms dry run had passed on the loose IC.
  - These slivers sit just above MFiX's `TOL_SMALL_CELL` (default 0.01).
  - So C26 either set a floor or had no such slivers. Asked in the author
    email.
- **Floor 0.30 A/B (`c01_smoke_r6_clip030`).** This is the r3 deck with the
  floor lowered to 0.30 and LINEAR_HAT kept, so one variable changes against
  r3. The criteria were pre-registered: a wall v_p change ≥ 30 % is
  substantial, < 15 % means no effect. Metrics are over PIV frames
  t = 0.1–0.3 s.
  - **Health:** check_smoke GO, 0 FAIL / 0 FLAG. dP +728 Pa = 1.76 W/A;
    S3 ratios 1.15 / 0.94.
  - **Wall v_p:** 1.79e-2 → 2.06e-2 m/s (+15.1 %), which falls in the
    **inconclusive** band. It moved up, away from the C26 model value of
    3.31e-3. The r3 halves alone differ by 1.5× (2.17e-2 / 1.41e-2), so a
    15 % change is within the noise of a 0.2 s window.
  - **Peak axis v_p:** +3 %.
  - **The floor does move:** at-floor cells drop from 62 % to 28 % in the
    wall band and from 56 % to 0 % in the annulus. Median EP_G goes from
    0.420 to 0.402 (wall) and 0.417 (annulus).
  - **Run time:** 4882 s against 4655 s for r3.
  - **Reading:** lowering the floor clearly changes the bed packing but does
    not slow the wall layer. The floor is not a credible cause of the gap to
    C26: 3.4× their model in the paper's window, t = 2.4–3.0 s, from the
    3 s c01 run.
- `--epg-clip` and `--interp-scheme` in `src/run_C26_3d.py` set these per run.
  `des_interp_width_cells` (geometry.yaml) gives the DPVM width.

## AWS run cost: open performance issues (2026-10-19)

Per step, t = 3–10 ms, same deck (dx/2 STL):

| run | ranks | fluid | CGDEM | total |
|---|---|---|---|---|
| laptop M4 (`c01_stl_h_timing`) | 4, i = k = [13,14] | 0.090 s | 0.040 s | 0.130 s |
| AWS c7a, c01 baseline | 9, uniform 3×1×3 | 0.090 s | 0.120 s | 0.210 s |
| AWS c7a, `c01_aws9_timing` | 9, i = k = [11,5,11] | 0.090 s | 0.110 s | 0.200 s |

So AWS runs this case ~1.6–1.7× slower per step than the laptop, even with
more ranks. Over a 3 s run that is 22.3 h on AWS (0.264 s/step) against ~13 h
locally at r3's whole-run rate (0.155 s/step).

- **Not the particle load per rank.** Uniform 3×1×3 puts 47 % of the parcels
  (13,077 of 27,782) on the centre rank. The 11/5/11 split cut the busiest
  rank's parcels 2.6× and gained only 1.05× (pre-registered threshold 1.15×),
  so it was not adopted and the presets were removed. The earlier
  `c01_gridmap_check` (8× fewer parcels on the busiest rank, 25 ranks) was
  slower too.
- **Suspected causes** (1 tested, see "Test 1 result"; 2 and 3 untested):
  1. **Per-rank overhead grows with rank count.** CGDEM went from 0.40 s
     (9 ranks) to 0.52 s (25 ranks) with the dx/4 STL. Candidates are ghost
     particle exchange each DEM substep (~5 per fluid step) and STL facet
     work at subdomain boundaries. Test: AWS with the laptop's 4-rank split
     (i = k = [13,14]) for 0.01 s. ≈ 0.13 s/step means fewer ranks is the
     fix; ≈ 0.20 s/step means hardware or container.
  2. **Per-core speed.** The fluid solve takes 0.09 s/step on both machines,
     even though the busiest AWS rank has fewer cells (7,792 vs 12,056). That
     implies ~1.5× slower cores on c7a than the M4 (inferred, not
     benchmarked).
  3. **Container MPI settings.** `docker/entrypoint.sh` runs mpirun with
     `--bind-to none`, which overrides the image's core binding, and sets the
     shared-memory/TCP transport (`btl self,vader,tcp`). Worth checking only
     if test 1 points to hardware or container.
- **Test 1 result (`c01_aws4_timing`, job 743edb60).** AWS with the laptop's
  split (2×1×2, i = k = [13,14]): fluid 0.14, CGDEM 0.09 s/step (medians,
  t = 3–10 ms). The whole window ran at 0.234 s/step, slower than 9 ranks
  (0.209). So fewer ranks is not the fix; the gap is mostly hardware or
  container.
  - Same decomposition, AWS against laptop: fluid 1.56× slower, CGDEM 2.25×.
  - Rank count on AWS, 4 → 9: CGDEM 0.09 → 0.12 (exchange overhead; MFiX
    runs `DES_PAR_EXCHANGE` every DEM substep when numPEs > 1,
    `des/des_time_march.f:430`), fluid 0.14 → 0.09. The two roughly cancel.
  - DEM loses more than fluid (2.25× against 1.56×), which points at
    container MPI (test 3) on top of per-core speed. Not yet separated.
- **Cost (before the change below).** 16 vCPU were requested for 9 ranks, so
  7 sat idle.
- **Adopted sizing (2026-10-19): 4 ranks on 4 vCPU, 6000 MB.**
  - Per run, 4 ranks is ~12 % slower than 9 (0.234 against 0.209 s/step),
    but it uses a quarter of the vCPUs.
  - On the 32-vCPU quota that means 8 jobs at once instead of 2. The nine
    runs c01–c09 then take 2 waves (8 at once, the ninth queued) instead of
    5.
  - Memory (local 4-rank c01, W = 10): steady 0.6–0.7 GB summed over ranks,
    peak 1.09 GB during mesh/STL preprocessing.
  - `src/run_C26_3d.py` now defaults to this split and writes
    `aws_vcpu`/`aws_memory_mb` to the manifest; `src/aws/submit.py` reads
    them. c02 and c03 were regenerated with it. c01 keeps its 3×1×3 deck,
    which produced the baseline.
  - Every case is W = 10. `src/run_C26_3d.py::aws_memory_mb` scales the
    probe linearly by parcel count (c09, 36.4k parcels: 8000 MB). The local
    0.3 s smokes peaked at 285–315 MB RSS (c04/c07 mass-only, c07_q18).
- **Where runs go:** smokes run locally (4 ranks). The 3 s runs
  c01–c09_massonly run on AWS at 4 vCPU each (submitted 2026-10-19/20),
  accepting ~1.7× the laptop wall time.

## Coating-stage / charge sweep (replaces Table 2 cases 4–12, 2026-10-19)

This project's sweep no longer reproduces C26's sensitivity cases 4–12, which
varied wall friction, statistical weight and steel feedstock. The matrix in
`params/sweep_C26.yaml` instead varies gas flow, particle diameter (TRISO
coating stage) and bed height (bed charge), one factor at a time around
**c02**: YSZ, 60 g, 12 slpm. It builds on the C26 model only, with no
correlations or operating points from Case V/D. The decisions are the user's,
2026-10-19.

| id | block | particle | d_p µm | ρ_p kg/m³ | charge g | H_static mm | Q slpm | parcel d mm | dx/d_parcel | parcels |
|---|---|---|---|---|---|---|---|---|---|---|
| c01 | flow | YSZ | 411.6 | 6000 | 60 | 34.8 | 10 | 0.89 | 2.37 | 27.3k |
| c02 | **baseline** | YSZ | 411.6 | 6000 | 60 | 34.8 | 12 | 0.89 | 2.37 | 27.3k |
| c03 | flow | YSZ | 411.6 | 6000 | 60 | 34.8 | 14 | 0.89 | 2.37 | 27.3k |
| c04 | diameter | YSZ + buffer | 611.6 | 2524 | 60 | 47.9 | 12 | 1.32 | 1.59 | 19.8k |
| c05 | diameter | + IPyC | 691.6 | 2332 | 60 | 49.7 | 12 | 1.49 | 1.41 | 14.8k |
| c06 | diameter | + SiC | 761.6 | 2550 | 60 | 47.7 | 12 | 1.64 | 1.28 | 10.2k |
| c07 | diameter | + OPyC (full TRISO) | 841.6 | 2381 | 60 | 49.2 | 12 | 1.81 | 1.16 | 8.1k |
| c08 | charge | YSZ | 411.6 | 6000 | 40 | 30.2 | 12 | 0.89 | 2.37 | 18.2k |
| c09 | charge | YSZ | 411.6 | 6000 | 80 | 38.6 | 12 | 0.89 | 2.37 | 36.4k |

**c01–c03** are C26 Table 2 cases 1–3. They keep their PIV targets, and
their deck physics is unchanged: a regenerated deck differs from the
existing one only in comment and `description` lines, so the c01 3 s run
stays valid. **c04–c09 are not C26 cases 4–9.**

- **Particles.**
  - The base particle is C26's YSZ (C26 calls it a TRISO-kernel surrogate).
  - Each diameter level adds the next TRISO layer, cumulatively and in
    deposition order:

    | Layer | Thickness | Density |
    |---|---|---|
    | buffer | 100 µm | 1000 kg/m³ |
    | IPyC | 40 µm | 1900 kg/m³ |
    | SiC | 35 µm | 3200 kg/m³ |
    | OPyC | 40 µm | 1900 kg/m³ |

    Values are from Marshall et al. 2017, Table 1, in `particles.yaml`
    `triso_coating_layers`. The AGR-5/6/7 Table 8 spec is given in comments.
  - Stage d and ρ come from `src/correlations.py::coated_particle_stages`.
    Its unit test reproduces Marshall's published stage values for a 425 µm
    UCO kernel.
  - **Assumed:** layers are uniform, so the YSZ PSD σ (12.21 µm) carries
    through coating unchanged.
- **Flow.**
  - The flow block is C26's own 10 / 12 / 14 slpm.
  - The diameter and charge blocks run at a fixed 12 slpm, the user's
    "optimal" value. It is also C26 Table 3's best case (c02: SSI 0.7687,
    MSE 0.0042).
  - Flow is deliberately **not** scaled with stage or charge.
- **Statistical weight.**
  - W = 10 for every case, and **W ≥ 10 is a lower bound** (the builder
    raises for W < 10).
  - Larger particles therefore get larger parcels, d·10^(1/3). For c04–c07
    that means dx/d_parcel falls from 1.59 to 1.16, against c01's 2.37.
  - **Risk:** expect noisier voidage and possible gas-fraction-floor or
    sliver-cell trouble (Error 1100 history above).
  - **Fallback:** if a smoke fails, the C26 builder's own
    `--interp-scheme SQUARE_DPVM` is the first candidate. It is proposed
    before it is applied (rule 9).
- **Fixed parameters.**
  - `mew_w` is fixed at the Table 1 baseline of 0.1 (`geometry.yaml` key
    renamed from `mew_w_base`).
  - All other closures are unchanged. Every case keeps the paper's
    dt_max = 1e-5 s, since U_nozzle is at most 52.4 m/s.
- **Diameter confounds bed height.**
  - At a fixed 60 g charge, coating lowers the particle density, so the
    static bed rises from 34.8 mm to 47.7–49.7 mm. That is above the
    41.6 mm cone top.
  - So the diameter block also varies bed height and bed weight per unit
    area. This is inherent to real layers at a fixed charge.
  - The SiC stage is denser than the IPyC stage, so c06's bed is slightly
    shallower than c05's.
- **Charge block.** In the 60° cone, H grows only as m^(1/3), so 40 / 60 /
  80 g gives 30.2 / 34.8 / 38.6 mm (about ±12 %).
- **Seeded charge** (1 ms dry runs, 2026-10-19, inventory monitor at t = 0):

  | case | c04 | c05 | c06 | c07 | c08 | c09 |
  |---|---|---|---|---|---|---|
  | seeded mass | 60.5 g | 58.8 g | 61.1 g | 60.1 g | 39.6 g | **76.1 g** |
  | error vs charge | +0.9 % | −2.0 % | +1.9 % | +0.1 % | −0.9 % | **−4.9 %** |

  For comparison, c01–c03 seed 61.0 g. Every case is within the 10 %
  seeding rule.
  - **Cause of c09's shortfall:** it is about one HEXA lattice row short
    at the top of the seeded region. That row sits in the column, above
    the cone top, where one row is about 1,750 parcels.
  - **Status:** accepted as is (user decision 2026-10-19). The charge block
    is 39.6 / 61.0 / 76.1 g as seeded.
  - **Reporting:** the sweep summary uses the seeded mass, not the
    nominal charge.
- **Spouting at 12 slpm is not guaranteed.** The deeper coated beds and the
  80 g bed may spout weakly or not at all. That is a result to report, not a
  reason to change the flow.

### Local smokes and criteria changes (2026-10-19)

0.3 s, 4 ranks, `notes/smoke_validate_c26.py`; t = 0.1–0.3 s. Rows in
`results/run_log.csv`.

| case | s/step | core flux (≥ 5 ṁ_g) | centroid (≤ 1.05 mm) | mid-bed v_g (≥ 5 m/s) | ΔP / (W/A) | lost parcels (all via top PO) | verdict |
|---|---|---|---|---|---|---|---|
| c02 (YSZ, 60 g) | 0.169 | 34.3 | 0.04 mm | 17.5 m/s | 1.36 | 3 (0.011 %) | GO, spout |
| c07 (full TRISO) | 0.121 | 1.6 | 1.43 mm | 1.64 m/s | 0.22 | 0 | GO (inlet), **no spout** |
| c09 (YSZ, 80 g → 76.0 g seeded) | 0.161 | 35.8 | 0.48 mm | 15.2 m/s | 1.38 | 2 (0.006 %) | GO, spout |
| c07_q18 (full TRISO, 18 slpm) | 0.159 | 16.1 | 0.39 mm | 18.0 m/s | −2.00 | 0 | **NO-GO** (S6), spout |
| c04, mass-only (YSZ + buffer) | 0.132 | 19.2 | 0.41 mm | 2.25 m/s (reported) | 2.05 | 0 | GO, spout |
| c07, mass-only (full TRISO) | 0.113 | 0.2 | 3.53 mm | 0.73 m/s (reported) | 2.19 | 0 | GO (inlet), **no spout** |

The first four rows use the momentum source, the last two the mass-only
default (criteria below).

- **Outlet loss (criterion amended, user decision).** The c02 smoke lost 3
  parcels. Particle_ID tracking shows each was last seen 1–9 mm below the
  pressure outlet at y = H_dom = 210 mm (C26's domain top), moving up at up
  to 4.7 m/s: ballistic fountain parcels, not wall leaks. c01 (10 slpm) lost
  none in 3 s. P2 was "the inventory stays at the seeded value". It is now:
  parcels may leave only through the top outlet (Particle_ID tracking), at
  most 0.5 % of the seeded mass; any other loss fails. The domain stays at
  210 mm, as in the paper.
- **c07 does not spout at 12 slpm.** The inlet is healthy: PS in the axis
  cell, exit v_g 1.19 U, inventory constant. The jet only opens a cavity:
  axis EP_G is > 0.78 up to y ≈ 22 mm and 0.42–0.47 above. Axis v_g falls
  from 13.6 m/s at 19.9 mm to 0.2–0.9 m/s at 28–47 mm.
  - The bed comes to rest: parcels moving > 5 cm/s fall from 69 % to 2.5 %,
    and there is no fountain.
  - Upward core flux decays from 4.3 to 0.9 ṁ_g over the run.
  - The bed is below minimum spouting: the particles are 2× larger and
    2.5× lighter, and the bed is 41 % deeper.
  - The `des_epg_clip` floor (0.42) is active in 71–72 % of dense-bed cells
    in both c02 and c07, so it is not specific to c07.
  - A 0.3 s smoke cannot exclude a spout later (C26: ≥ 3 s).
- **Spout checks are now results outside the PIV cases (user decision).**
  The validator splits its verdict:
  - P1–P2 plus `check_smoke` gate inlet health and numerics.
  - P3–P4 report whether the bed spouts. A spout is required only for c01–c03
    (manifest `piv_targets_table3`), where C26 observed one.
  - The misplaced-source c01 smoke is still NO-GO through P1, and
    c01_smoke_r3 is still GO.
- **c07 spout check at 18 slpm (pre-registered before running).** The user
  kept c04–c07 at 12 slpm and added one check run:
  `run_C26_3d --case c07 --flow-slpm 18 --ps-momentum --run-dir-suffix _q18_smoke`
  (momentum source, the default at the time).
  - **Why 18 slpm:** a rough Mathur–Gishler scaling, used only for this
    choice and not in the model, gives c07 about 1.5× bare YSZ's minimum
    spouting flow. So 18 = 1.5 × 12 slpm puts c07 at about c02's excess over
    minimum spouting.
  - **Time step:** U_nozzle 67.4 m/s still gives dt_max = 1e-5 s.
  - **Spout criterion:** P3i ≥ 5 ṁ_g and P4 ≥ 5 m/s, the same thresholds.
  - **If it spouts:** c07's no-spout at 12 slpm is physical, and 1.81 mm
    parcels can carry a spout in this model.
  - **If it doesn't spout:** the run is inconclusive. Either c07's minimum
    spouting flow is above 18 slpm, or 1.81 mm parcels in 2.1 mm cells hit a
    numerical limit. The next step would be proposed then.
  - **Result (`c07_q18_smoke`, 2026-10-19): it spouts.** P3i 16.1 ṁ_g,
    P3ii 0.39 mm, P4 18.0 m/s, so by the rule above c07's no-spout at
    12 slpm is physical. `check_smoke` S6 failed (mean ΔP −590 Pa =
    −2.00 W/A): a static-pressure depression in the 67 m/s momentum source
    cell, not an unsupported bed (the bed above row 0 is pressure-supported).
    It does not apply to the mass-only decks.

### AWS production runs: mass-only point source (user decision 2026-10-19)

- **Decision.** The sweep's 3 s runs use the mass-only point source
  (`--ps-mass-only`, `ps_v_g = 0`), matching the c01–c03 runs already on AWS
  (`c26-c0X-massonly`, submitted 2026-10-19).
  - **Project default changed (user decision, same day).**
    `params/geometry.yaml` now has `ps_momentum: false`, so every new C26
    deck is mass-only.
    - `--ps-momentum` restores the momentum source for A/Bs.
    - `--ps-mass-only` is kept and now equals the default.
    - The existing momentum artifacts are left as they are and not
      regenerated: `cases/C26_3d/c01` (the completed 3 s momentum run and its
      `aws_out`), `c02` and `c03` (momentum decks, never run), and the
      `c02`/`c07`/`c09`/`c07_q18` smokes.
  - The decks differ from the momentum decks only in `ps_v_g`, the run name
    and the monitor names.
  - Run dirs are `cases/C26_3d/<id>_massonly`, and job names
    `c26-<id>-massonly`.
  - The sweep summary reads them with `--suffix _massonly`.
- **Caveat.** The c01–c03 mass-only runs are the A/B test of the point-source
  momentum (see "Point-source momentum A/B"). Mass-only was adopted before
  their 3 s results, by user decision. The c04–c09 runs follow the c02
  mass-only baseline, whichever source scores better against the C26 PIV.
- **Mass-only smoke criteria** (set before any mass-only c04–c09 result).
  - `check_smoke` S3 (exit v_g / U_nozzle band 0.5–2) and P4 (mid-bed
    v_g ≥ 5 m/s) were set on momentum-source runs. For `ps_momentum: false`
    they are reported, not gated. On c01_smoke_r4_massonly, S3 gave 0.39
    and P4 3.2 m/s, while P3 shows a spout (37.8 ṁ_g, 0.18 mm).
  - All other inlet and numerics checks gate as before, and the spout is
    judged on P3.
  - The momentum smokes (c02, c07, c09) keep their verdicts.
- **Automation** (`smoke_chain5.sh`, scratchpad; user: queue on AWS if the
  smokes pass).
  1. Generate `c04`–`c09_massonly` (tstop 3 s, 4 ranks).
  2. Run a 1 ms serial dry run of each in scratch (rule 13).
  3. Run a mass-only c04 smoke (`c04_massonly_smoke`, 0.3 s).
  4. Submit via `src/aws/submit.py`, sizing from each manifest, with these
     gates:
     - the c04 mass-only smoke must be GO, or nothing is submitted;
     - a case whose dry run fails is held;
     - c05–c07 are held unless c07 at 18 slpm (momentum) is GO and spouts.
  - **Outcome.** The c04 mass-only smoke was GO and spouts, so c04, c08 and
    c09 were submitted (2026-10-19/20). `c07_q18_smoke` was NO-GO on S6 but
    spouts, so c05–c07 were held; they were submitted on 2026-10-20 at the
    user's direction, c07 last (lowest priority, no spout at 12 slpm). All
    nine `c26-<id>-massonly` jobs are submitted.
- **Momentum smokes still apply.** They cover parcel size up to 1.81 mm
  (c07), parcel count up to 34,640 (c09), and outlet loss. These do not
  depend on the source model.
- **Timeout (rule 14).** The Batch job definition sets no
  `attemptDurationSeconds`, so jobs are not time-limited.
  - Expected time: 300k steps at about 0.2–0.25 s/step (AWS 4-rank c01
    0.234), roughly 17–21 h per case.
  - Capacity: 9 jobs × 4 vCPU = 36 vCPU. With a 32-vCPU quota, one job
    waits in RUNNABLE until a slot frees.
- **Tools.**
  - Smoke validation is `notes/smoke_validate_c26.py <run_dir>`. It runs
    `check_smoke.py` with arguments from the manifest, plus the step-2
    checks: PS volume, seeding/inventory, spout-core flux and centroid,
    mid-bed v_g, and settled height.
  - Step-3 metrics come from `notes/c26_step3_checks.py`, which now reads
    the manifest; the C/D thresholds apply to c01–c03 only.
  - The sweep summary is `notes/c26_sweep_summary.py`, which writes
    `results/c26_sweep_summary.csv` and one figure per block. The averaging
    window is in each caption.
- **AWS memory.** `aws_memory_mb` scales linearly with parcels above the
  c01 probe, so c09 gets 8000 MB.

## Open items (not done here)

- Point-source momentum and wall velocity: asked the C26 authors for the
  .mfx/STL, the PS velocity, the apex geometry and probe definitions
  (2026-10-19); no reply yet. The 3 s mass-only c01–c03 runs (submitted
  2026-10-19) are the A/B: score them with `notes/c26_step3_checks.py`
  (C1/D, t = [2.4, 3.0] s) when they return. The MFiX 23.2 vs 26.1.2
  wall-contact review is done (see "Gas-fraction floor"); the apex geometry
  is still open.
- Fig. 8e cone angle: measure it from a high-resolution extraction before
  changing geometry (the text says 60°).
- AWS per-step cost: the 4-rank AWS test is done (the gap is hardware or
  container, not rank count; 4 ranks adopted). Per-core speed and the
  container MPI settings (`--bind-to none`, btl) are untested.
- Synthetic-PIV post-processor and digitized C26 Fig. 4 / wall-velocity data.
- If the apex cut cell misbehaves, the fallback is to move the point source up
  one row (`ps_row_offset_cells: 1`). It should not be a widened inlet
  (recipe Rule 3 history).

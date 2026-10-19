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
`case_C26_ysz`, `case_C26_ss`, `ar_25C`, `spouted_bed_C26_3d_numerics`,
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
  `set_ps.f:CALC_PS_VEL_MAG`). We set `ps_v_g` = nozzle exit velocity
  Q / (π D_n²/4). The jet then carries the physical momentum flux ṁ·U_n.
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
for the YSZ cases. Rule 9 binds only for stainless steel at 30 slpm
(9.4e-6 s). This deliberately does not copy the V_3d decks' 1e-3 s ceiling,
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
   **~35 mm**. ~25 mm is only consistent with a ~90° cone. The builder logs
   both numbers in every manifest. The first smoke must report the settled
   bed height so this can be judged.
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
- All other keywords are reused unchanged from `cases/V_3d/template_dem.mfx.j2`.
- Cross-check: all 135 keyword names in the rendered c01/c09/c12 decks were
  grepped against every `*init_namelist*.f` in the installed source; none
  are missing. `tests/test_make_case_c26.py::test_deck_keywords_are_verified`
  keeps the deck within the V_3d-template set plus the list above.

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
   - The support is kept for later use. The baseline uses uniform 3×1×3.
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

   - Still to check: wall roughness. At dx/2 a facet edge is ~1 mm against
     0.41 mm particles. Compare the wall-velocity target (Fig. 6d, ~1e-3 m/s)
     against a dx/4 run if c01's wall velocity comes out off.
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
flow at the wall. Candidate causes (not yet tested; nothing changed):
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
- PS momentum does not explain the 3× wall velocity.
- r4 fails `check_smoke` S3 (v ratio 0.39/0.24) and the mid-bed v_g ≥ 5 m/s
  supporting check. Those were written for a momentum source, so a
  mass-only baseline would need its own inlet criteria before adoption.
- `ps_momentum` stays `true`.
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
- **Next A/B.** `c01_smoke_r6_clip030` lowers the floor to 0.30 with
  LINEAR_HAT, one variable against r3. Criteria were pre-registered (wall v_p
  change ≥ 30 % substantial, < 15 % no effect).
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
- **Suspected, not yet tested:**
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
- **Cost.** 16 vCPU are requested for 9 ranks, so 7 sit idle.
- **Until resolved:** C26 smokes run locally (4 ranks). Long runs go to AWS
  only when the laptop can't be spared, accepting ~1.7× the wall time.

## Open items (not done here)

- Point-source momentum and wall velocity: asked the C26 authors for the
  .mfx/STL, the PS velocity, the apex geometry and probe definitions
  (2026-10-19). Options pending their reply:
  - a 3 s mass-only c01 (~22 h);
  - looking for the wall-velocity cause elsewhere (MFiX 23.2 vs 26.1.2
    cut-cell wall contact, apex geometry).

- AWS per-step cost: the dx/2 baseline ran at 0.264 s/step overall (0.210
  over t = 3–10 ms). See "AWS run cost: open performance issues" for the
  untested causes and the 4-rank AWS test.
- Synthetic-PIV post-processor and digitized C26 Fig. 4 / wall-velocity data.
- If the apex cut cell misbehaves, the fallback is to move the point source up
  one row (`ps_row_offset_cells: 1`). It should not be a widened inlet
  (recipe Rule 3 history).

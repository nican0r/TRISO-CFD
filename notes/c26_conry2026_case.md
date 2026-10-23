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
  the design until 2026-10-19 and is kept for A/Bs. The nine
  `c0X_momentum24` runs (2026-10-22) use it on the 24-cell layout; see
  "First-principles review of the velocity gap" below.
- **Axis placement.** Since 2026-10-22 the box is exactly 24 cells across
  (= D_c, as C26 p.3 states; `n_margin_cells_xz: 0`, `axis_on_cell_face:
  true`). The vessel axis sits on the corner shared by the four central
  cells and the point source covers those 2 × 2 cells (±2.1 mm in x and z).
  Layout history: every run to 2026-10-20 used 27 cells (margin 1, odd
  count, axis at a cell centre, one source cell at (0, dy/2, 0), which
  matched C26's "cell ≈ nozzle" argument); the `_massonly_axisface` runs
  (2026-10-21) used 26 (margin 1, axis on the corner). 26 with an all-blocked
  margin and 24 with none have the same fluid cells (42,520; 8,608 cut),
  the same STL and the same source; only the box extents and the DMP gridmap
  ([12,12] instead of [13,13]) change, so there is no physics difference.
  See "First-principles review of the velocity gap" below.
- Gas point sources work with DEM/CGP solids in 26.1.2. Only *solids* point
  sources are rejected for DEM/CGP (`check_point_sources.f:443-451`).

### Recipe Rule 6 (4-quadrant IC with a central gas column) → one bed IC region

- R5 (`notes/day10_r5_ic_review.md`) found that the quadrant pattern copies
  NETL's draft-tube bore, which this geometry does not have. With a closed cone
  apex and no stub, there is also no orifice for seeded particles to fall
  through. That was the original Fix-A motivation.
- C26 §2.2.2: "Particles were initialized at the bottom of the reactor and
  followed the particle size distribution". We seed one IC region over the
  full cross-section from y = 0 to H_ic, with ε_s,ic = 0.30 (assumed). H_ic
  is the height the charge would fill at 0.30. The bed settles under gravity
  while the point source flows.
- **The seeded bed is not loose (checked 2026-10-20).** In MFiX 26.1.2,
  `ic_ep_s` sets only the parcel count, `floor(ic_ep_s · IC_VOL / P_VOL)`
  (`des/generate_particles_mod.f:731`). The parcels go on a HEXA lattice
  spaced at (1 + `ic_des_spacing`) × the largest parcel diameter (L849–856),
  filled from y = 0 up until the count is reached (L941–1033).
  - So the bed starts at the lattice packing, not at 0.30: about 0.50 for
    YSZ (c01–c03, c08, c09) and 0.56 for c07. In the c02 smoke at t = 0 the
    median nearest-neighbour distance is 1.014 mm (= 1.05 × 0.966 mm) and
    the interior ε_s is 0.51.
  - The top of the IC region starts empty. Seeded bed tops at t = 0 (99th
    percentile of parcel y): c02 38.1 mm (H_ic 43.7), c04 51.5 (66.8),
    c07 51.9 (69.1), c09 40.6 (49.3).
  - The lattice then compacts under gravity to an interior ε_s of 0.59–0.62
    within ~40–80 ms. The settled height is set by the seeded charge and
    that packing. H_ic enters no check; the monitors and thresholds use
    H_static.

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
particles into the camera plane and runs the same PIV. `src/post/piv_c26.py`
does this (see "Synthetic PIV" below). The deck writes the data it needs:
particle VTK with ids every 2 ms for 300 frames over the last 0.6 s
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
    bed is still compacting from the seeded lattice, so the early wall value
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
    settles. A 3 ms dry run had passed on the seeded IC.
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
  - **Cause of c09's shortfall:** the parcel count, not the lattice.
    MFiX counts parcels from `IC_VOL`, the summed volume of whole fluid
    cells in the IC index range (`GET_IC_VOLUME`,
    `des/generate_particles_mod.f:2623`), so `ic_y_n = H_ic` acts as if
    snapped to the nearest cell face. c09's H_ic of 49.3 mm (23.48 cells)
    snaps down to 48.3 mm. That drops ~1 mm of column, ~3.6 g at
    ε_s 0.30.
  - The other cases' errors have the same cause. Their signs follow the
    snap direction: up for c01–c04, c06 and c07; down for c05 and c08. The
    magnitudes also carry the cut-cell volume error against the nominal
    cone.
  - Corrected 2026-10-20: this note first blamed one missing HEXA lattice
    row. The lattice fills until the count is reached, so it is the count
    that is short.
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
    0.234), roughly 17–21 h per case. **Too low:** c01–c03 took 22.7–23.5 h
    (see "c01–c03 mass-only 3 s results"). The 0.234 s/step came from a
    10 ms window, and the per-step cost grows as the bed develops.
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

## c01–c03 mass-only 3 s results (2026-10-20)

**Runs.** The AWS jobs `c26-c01/c02/c03-massonly` ran on 4 ranks (2×1×2),
4 vCPU each on a shared c7a.4xlarge. All three reached t = 3.000 s with
mfix exit 0.

| case | flow | wall time | steps | steps with dt < 1e-5 s (min dt) |
|---|---|---|---|---|
| c01 | 10 slpm | 23.52 h | 300,653 | 1.8 % (6.6e-6 s) |
| c02 | 12 slpm | 22.71 h | 300,711 | 2.0 % (7.3e-6 s) |
| c03 | 14 slpm | 23.21 h | 300,917 | 2.5 % (6.6e-6 s) |

- No parcels were lost: each run holds 27,782 parcels and 60.97 g at
  t = 3 s.
- Outputs are in `cases/C26_3d/<id>_massonly/aws_out` (gitignored).
- Each run is logged in `results/run_log.csv` (2026-10-20).

**Health** (`notes/smoke_validate_c26.py`, t ≥ 0.1 s).
- **c01 and c02 are GO.**
- **c03 is NO-GO on P3ii only.** Its spout-flux centroid is 1.18 mm off the
  axis, against a 1.05 mm limit. The validator calls this "no spout", but a
  spout is present (P3i 33.1 ṁ_g, a fountain in every frame); it is
  off-centre (next section).
- **All three:**
  - `check_smoke` GO; ΔP 2.58 / 2.65 / 2.77 W/A; S8 peak 10.7–11.0 Hz;
  - S3 exit ratio 0.35–0.36 (mass-only, reported, not gated);
  - P1 7.013e-9 m³; P2 PASS with 0 parcels lost.

**Step-3 checks.** The window is t = [2.4, 3.0] s with 300 PIV frames; the ΔP
spectrum is also taken over 0.5–3.0 s. Thresholds were pre-registered on
2026-10-19 (table in "c01 baseline result").

| criterion | threshold | c01 | c02 | c03 | c01 momentum |
|---|---|---|---|---|---|
| A fountain | ≥ 50 % of frames | **45 % (median 12) FAIL** | 96 % (336) PASS | 100 % (1088) PASS | 100 % PASS |
| B ΔP PSD peak, 2.4–3.0 / 0.5–3.0 s | 10–20 Hz | 11.7 / 11.2 PASS | 10.0 / 10.8 PASS | 10.0 / 10.4 PASS | 18.3 / 18.4 PASS |
| C1 axis v_p, peak of y = 10/20/30 mm | 0.26–1.06 m/s | 0.42 (y = 20) PASS | 0.47 (20) PASS | 0.53 (30) PASS | **1.75 FAIL** |
| C2 outermost bin \|v_p\| | < 0.1 × axis, each height | **FAIL at y = 10** (0.057 vs 0.038) | **FAIL at 10** (0.074 vs 0.041) | **FAIL at 10** (0.076 vs 0.041) | PASS |
| D wall v_p (downward +) | 0.3e-3–10e-3 m/s | **2.23e-2 FAIL** | **2.92e-2 FAIL** | **3.39e-2 FAIL** | **1.11e-2 FAIL** |

- **B:** c02 and c03 sit on the 10 Hz edge. The 0.6 s window resolves
  1.67 Hz.
- **C2:** the pre-registration names no probe height. The momentum run
  passed at every height, so the strict reading (every height) is used.
  - All three pass at y = 20 and 30 mm.
  - The y = 10 mm failure is the outermost bin, x = ±6.3 mm, in the cone
    bottom, where the wall is at r = 7.0 mm. It is annulus return flow next
    to the wall.
- **C1 and D bands** come from C26's 10 slpm figures (Fig. 7, Fig. 6d). For
  c02 and c03 they are only a reference.
- **D pulses.** The sd over frames exceeds the mean in all three runs
  (3.3e-2 / 4.6e-2 / 6.3e-2 m/s).
- **Not assessed:** the Table 3 PIV maps (step 4, post-processor not
  written).

**Point-source A/B scored on the 3 s runs.** The rule was pre-registered
under "Point-source momentum A/B": substantial only if both wall and peak axis
v_p fall by ≥ 30 %.
- Peak axis v_p went from 1.75 to 0.42 m/s (−76 %).
- Wall v_p went from 1.11e-2 to 2.23e-2 m/s (+101 %).
- **Verdict: inconclusive** under the rule.
- **Reading:** source momentum explains the excess axis speed (C1) but not
  the excess wall speed (D). Without it, D doubles.
- **Caveat:** the momentum c01 ran on 9 ranks (3×1×3); the mass-only runs
  used 4 (2×1×2).

### Spout off-centre in the mass-only runs (diagnosis 2026-10-20)

**Observation.** In all three mass-only runs the spout leans towards −x.
- **Particles:** the upward particle-flux centroid at y = 25.2 mm,
  t = 2.4–3.0 s, sits at:

  | case | x (mm) | z (mm) |
  |---|---|---|
  | c01 | −0.84 | −0.05 |
  | c02 | −0.99 | +0.02 |
  | c03 | −1.17 | +0.21 |

  It is the same in every 0.1 s sub-window. The momentum c01 run gives
  0.10 mm with no persistent offset. Momentum smokes on the same 4-rank split
  are centred (|x| ≤ 0.04 mm).
- **Gas:** the upward gas-flux centroid (BACKGROUND frames, t = 0.5–3.0 s,
  y = 3–18 mm) sits at x = −0.6 to −1.4 mm, with z within 0.1 mm. The
  momentum c01 run is within 0.03 mm.

**Where it starts: the source cell.**
- **Reading the data.** BACKGROUND `Gas_Velocity` holds the staggered face
  values. The +x and +z neighbours read exactly 0, because their far faces
  are wall. The face fluxes out of the source cell close its mass balance:
  for c02, 2.007e-4 m³/s × 1.633 kg/m³ = 3.277e-4 kg/s, against
  ṁ_g = 3.266e-4 kg/s (+0.3 %).
- **c02 face speeds** (time-mean outflow from the source cell, t = 2.0–3.0 s):

  | face | +x (E) | −x (W) | +z (T) | −z (B) | top (N) |
  |---|---|---|---|---|---|
  | speed, m/s | 6.91 | 11.19 | 10.53 | 10.66 | 15.75 |
  | open area, mm² | 3.34 | 3.34 | 3.34 | 3.34 | 4.41 |

- **Most of the gas leaves sideways:** about two-thirds of the source gas
  exits through the side faces. The +x face carries 38 % less than the −x
  face; the two z faces match.
- **The E/W ratio is fixed from the start.**
  - u_E/|u_W| is 0.58–0.65 in every time window from the first output frame
    (0.05 s). Per frame it ranges 0.37–0.75 and is never ≥ 1.
  - The z ratio w_T/|w_B| is 0.99–1.07 (per frame 0.88–1.47).
- **Every mass-only run shows it:**

  | run | u_E/\|u_W\| |
  |---|---|
  | r4 smoke (c01) | 0.65 |
  | c04 smoke | 0.69 |
  | c07 smoke (does not spout) | 0.65 |
  | c01 / c02 / c03, 3 s | 0.61 / 0.62 / 0.62 |

**What it is not.**
- **Not the scalar cut-cell geometry.** In the BACKGROUND mesh, the source
  cell's four side faces have equal open area (3.34 mm²). Its four side
  neighbours have equal volume (1.853 mm³).
- **Not a bed instability.** The ratio is present at 0.05 s and never
  reverses. It has the same size in a bed that does not spout (c07), and z
  stays symmetric.
- **The momentum source hides it.** The imposed 44–53 m/s jet carries the
  gas up, and the side faces carry < 2 m/s inward.

**Candidates considered before testing.**
1. **Staggered U-momentum cut cells.** MFiX cuts the U, V and W cells
   separately, and they are not written to VTK. The STL comes from
   marching cubes: its vertices are mirror-symmetric, but its facet
   diagonals are not.
   - Below y = 6.3 mm, 230 of 334 facets have no mirror image under
     x → −x, and a similar number under z → −z.
   - Wall area below 6.3 mm is 56.56 mm² at x < 0 and 59.18 mm² at x > 0;
     in z it is 57.57 and 58.16 mm².
2. **Linear-solver sweep direction.** The deck uses MFiX's default `leq_*`
   settings.
3. **Domain split.** Rank boundaries lie on the source cell's −x and −z
   faces (gridmap 13|14 in both i and k).

**Gas-only tests (2026-10-20, local).**
- **Setup.** The c02 deck (12 slpm, mass-only) with no solids
  (`ic_ep_g(2)` 1.0, `ic_ep_s(2,1)` 0.0), BACKGROUND every 5 ms, 0.05 s.
  Nothing else in the deck changed (diffed against `c02_massonly`).
- **Decision rule**, set before the runs:
  - E/W = u_E/|u_W| and T/B = w_T/|w_B|, the source-cell outflow ratios,
    averaged over t = 0.01–0.05 s;
  - bias if a ratio is < 0.85, symmetric if it is 0.9–1.1.

| run | E/W | T/B | reading |
|---|---|---|---|
| as is, 4 ranks | 0.663 | 0.997 | the bias exists without a bed (production c02 with bed: 0.62) |
| serial | 0.664 | 0.997 | identical, so not the domain split |
| STL rotated 90° about y, (x, y, z) → (z, y, −x) | u_W = 0.00, u_E = 10.51 m/s | 0.980 | the bias did not move to z |

- All three runs are steady from t = 0.005 s, with 0 particles and 0
  errors.
- In the rotated run, the x-face flows change with the STL triangulation,
  but z stays symmetric in both orientations. So the STL alone is not the
  cause. Something in MFiX treats x differently from z, and that something
  is sensitive to the STL details.

**Code finding.** In MFiX 26.1.2, `src/model/cartesian_grid/get_alpha.f`
computes the cut-cell interpolation factors that weight face values between
neighbouring momentum cells.
- **The V and W factors are limited; the U factor is not.**
  - `Theta_Vn` (l. 538–544) and `Theta_Wt` (l. 954–960) are reset to 0.5
    when they fall outside (0, 1).
  - `Theta_Ue` (l. 107–108) has no such reset.
- **In this mesh** (`MESH_STATS.LOG`, 4-rank run):

  | factor (cut cells) | range | share outside (0, 1) |
  |---|---|---|
  | `Theta_Ue` (U) | 0.5–7.39 | about 21 % are above 1.19 |
  | `Theta_Wt` (W) | 0.5–0.975 | none (99.7 % at 0.5) |
  | `Theta_Vn` (V) | 0.5 | none |

  A factor above 1 means the east-face value is extrapolated, with a
  negative weight on the neighbour. That error is one-sided in x, which
  breaks the +x/−x mirror symmetry.
- **The rotated mesh gives the same `Theta_Ue` range** (max 7.3858), and W
  is unchanged. So the U-only extrapolation is tied to the x direction in
  the code, not to the STL.
- **The output reading holds.** VTK `Gas_Velocity` is the raw
  U_G / V_G / W_G (`vtk_out.f` l. 197), so the face velocities above are
  not an output artefact.
- **Status: consistent with every result, but not proven causal.** The
  test is a custom solver built from a copy of `get_alpha.f` with the same
  reset added for `Theta_Ue`, MFiX's standard user-source mechanism. Rerun
  the gas-only test with it: if this is the cause, E/W should go to ≈ 1.
  - Whether the reset is missing deliberately is not known. It has not
    been raised with the MFiX developers.
  - Possible wider effect (untested): U cut cells line the whole cone and
    column wall, so near-wall x-momentum is affected everywhere, not only
    at the source. The momentum-source runs use the same code.

**Effect on results so far.**
- Every mass-only run, c01–c09, carries this lean.
- It is at most 0.6 cells, along the PIV camera's horizontal axis. So it
  shifts the time-mean spout in the Table 3 maps.
- c03 fails P3ii because of it.

## Synthetic PIV (2026-10-20)

**Method.** `src/post/piv_c26.py`; settings and sources in
`params/piv_c26.yaml`; tests in `tests/test_piv_c26.py`.
- **Rendering.** Each `PARTICLES_PIV` frame becomes a side view along z at
  0.1 mm/px. Parcels are opaque spheres at their VTK (parcel) diameter, and
  the nearest surface wins.
- **PIV.** Consecutive frames are correlated: 299 pairs, 2 ms apart.
  - FFT cross-correlation, windows 64 / 32 / 32 px (6.4 / 3.2 / 3.2 mm),
    vectors every 16 px (1.6 mm).
  - Between passes, the windows are offset by the previous estimate to
    sub-pixel accuracy.
  - 3-point Gaussian peak fit; normalized median test after each pass.
  - An empty window counts as zero displacement.
- **Spout map.** The map value is the magnitude of the time-mean vector.
  - It spans x in [−R_c, R_c]. It runs up from the bed surface (time median
    of the 95th percentile of parcel y) by 2 R_c × the panel aspect ratio:
    39.8 mm at 10 and 12 slpm.
  - It is resampled onto the digitized Fig. 4 grid and scored with CC, SSI
    (MATLAB ssim defaults, L = 1) and MSE.
- **Wall velocity.** Mean vertical PIV velocity in a 4 mm box at x = 0,
  mid-bed (C26 Fig. 1c), downward positive.
- **Settings are assumed.** C26 gives no camera or PIV settings; it cites
  Conry et al. 2025, which we don't have. The "~4 mm" in Fig. 1c is the wall
  sampling box, not the interrogation window.

**Checks.**
- **Sub-pixel bias.** Integer window offsets left a bias of ~0.1 px toward
  the nearest whole pixel for 9 px parcel images. The sub-pixel offset brings
  it below 0.01 px. Tests cover shifts of 0.05–12 px.
- **Against the parcels on c01.** The PIV wall velocity is 1.24e-2 m/s; the
  parcels visible in the same box move at 1.18e-2 m/s (+5 %).
- **Table 3 not reproduced from the figure.** The digitized Fig. 4 maps do
  not reproduce Table 3 (`data/digitized/README.md`). C26 model vs experiment
  through this pipeline gives CC 0.39 / 0.87 / 0.84, against Table 3's
  0.68 / 0.82 / 0.83. So each run is reported next to Table 3 **and** next to
  C26's model scored through the same pipeline.

**c01, momentum source, 3 s run.** Window t = [2.400, 2.998] s, 299 pairs.
Output: `results/c26_piv_c01.csv`, `results/fig_c26_piv_c01.png`.

| 10 slpm | CC | SSI | MSE [m²/s²] | wall v [m/s] |
|---|---|---|---|---|
| c01, this pipeline | 0.214 | 0.589 | 0.0162 | 1.24e-2 |
| C26 model, same pipeline (Fig. 4b vs 4a) | 0.390 | 0.742 | 0.0057 | — |
| C26 model, Table 3 / Fig. 6d | 0.682 | 0.609 | 0.0068 | 3.31e-3 |
| C26 experiment, Fig. 6d | | | | 0.933e-3 |

- **Map shape matches.** Our map has C26's structure: two fast lobes of
  falling particles at the sides and a slow centre, where up and down motion
  in the line of sight cancels.
- **Our lobes differ.** They are faster (peak 0.45 against ~0.32 m/s), lower
  (~5 against 10–12 mm above the map bottom) and further out (|x| ≈ 15
  against 11–12 mm).
- **Wall velocity.** 13× the experiment and 3.7× C26's model. This agrees
  with step-3 D (1.11e-2 m/s from parcels), so the PIV method does not
  explain D.

**Sensitivity.** The same c01 output, with one assumption changed at a time:

| variant | CC | SSI | MSE [m²/s²] |
|---|---|---|---|
| base: map bottom at the bed surface (35.4 mm), camera along z | 0.214 | 0.589 | 0.0162 |
| camera along x | 0.230 | 0.572 | 0.0164 |
| map bottom +3.4 mm (the mask strip height) | 0.095 | 0.540 | 0.0168 |
| map bottom −3.4 mm | 0.615 | 0.672 | 0.0080 |

- **Camera axis.** Metrics move by ≤ 0.02 and the wall velocity by 1 %, so
  the 0.6 s average is converged enough for this comparison.
- **Map bottom dominates CC.** Lowering the map by 3.4 mm triples CC, to
  0.62, above C26's model through the same pipeline (0.39).
  - It brings the slow bed top into the bottom rows, like the dark band along
    the bottom of the experiment map, and moves our lobes up.
  - SSI and MSE stay worse than C26's model in every variant.
- **Default kept.** The map bottom stays at the bed surface, chosen before
  scoring. It is not moved to the best-scoring variant (no tuning to match).
- **Still needed:** a placement rule that does not depend on the score.
  Sources: the authors (the 2026-10-19 email is still open) or a
  measurement from the Fig. 4 snapshots.
- **Wall velocity** does not depend on the map position.

**c01–c03 mass-only (production source).** Window t = [2.400, 2.998] s, 299
pairs, default settings. Output: `results/c26_piv_c0X_massonly.csv` and
`results/fig_c26_piv_c0X_massonly.png`.

| run | flow | CC | SSI | MSE [m²/s²] | C26 model, same pipeline (CC / SSI / MSE) | Table 3 | wall PIV [m/s] | C26 exp / model wall |
|---|---|---|---|---|---|---|---|---|
| c01_massonly | 10 | −0.134 | 0.125 | 0.0171 | 0.390 / 0.742 / 0.0057 | 0.682 / 0.609 / 0.0068 | 3.07e-2 | 0.93e-3 / 3.31e-3 |
| c02_massonly | 12 | 0.130 | 0.132 | 0.0348 | 0.872 / 0.740 / 0.0022 | 0.817 / 0.769 / 0.0042 | 3.56e-2 | 1.83e-3 / 5.05e-3 |
| c03_massonly | 14 | −0.067 | 0.169 | 0.0283 | 0.838 / 0.741 / 0.0033 | 0.829 / 0.766 / 0.0049 | 3.85e-2 | 2.55e-3 / 6.95e-3 |

- **Almost nothing moves above the bed surface on average.**
  - The map means are 0.015 / 0.021 / 0.024 m/s, against 0.11 / 0.18 /
    0.16 m/s in the experiment.
  - 85–88 % of the map's windows are empty in a given pair, against 36 % for
    the momentum c01.
  - Only a low dome reaches ~8 mm above the bed surface (36.6 / 39.5 /
    43.3 mm). The step-3 fountain counts (A) still pass for c02 and c03:
    the parcels above the bed move up and down in place and largely cancel
    in the time mean. They do not form C26's two side lobes of falling
    particles.
- **Wall velocity.** 33× / 19× / 15× the experiment and 9× / 7× / 5.5× C26's
  model. The box pulses at ~11 Hz between about +0.1 and −0.03 m/s. The PIV
  agrees with the visible parcels within 2 %, and with step-3 D
  (2.2–3.4e-2 m/s from all wall parcels) in order of magnitude.
- **Sensitivity.**
  - Map bottom ±3.4 mm: CC −0.14 to 0.18, SSI ≤ 0.20, MSE ≥ 0.015.
  - Empty windows left out of the mean instead of counted as zero
    (`spout_map.empty_windows: skip`):

    | run | CC | SSI | MSE [m²/s²] |
    |---|---|---|---|
    | c01 momentum | 0.195 | 0.589 | 0.0314 |
    | c01_massonly | −0.020 | 0.155 | 0.0173 |
    | c02_massonly | 0.546 | 0.279 | 0.0222 |
    | c03_massonly | 0.281 | 0.274 | 0.0230 |

  - CC depends on these assumptions. SSI (≤ 0.28) and MSE (≥ 0.017) are far
    from C26's model through the same pipeline under every variant.
- **Point-source A/B at 10 slpm, on the C26 PIV.** Momentum c01 scores
  higher than c01_massonly on CC and SSI in every variant (SSI 0.54–0.67
  against 0.12–0.18).
  - MSE is about equal by default (0.016 / 0.017) and favours mass-only with
    empty windows skipped (0.031 / 0.017).
  - The wall velocity is closer to the experiment for momentum (1.24e-2
    against 3.07e-2 m/s), but both are far off.
  - Step 3 pointed the other way on C1: mass-only passes, momentum fails.
  - So neither source reproduces C26's maps. The momentum source gives the
    C26-like fountain structure, though with lobes that are too fast and
    too low.
  - This bears on the 2026-10-19 caveat above ("whichever source scores
    better against the C26 PIV"). It is not a decision; that stays with the
    user.

## Wall-velocity gap: diagnosis and the grid-layout A/B (2026-10-20)

**Is it the measurement? No.** Post-processing of the c01 / c01–c03
mass-only PIV frames, t = [2.4, 3.0] s (scratch scripts, no new runs):
- Synthetic PIV in the wall box equals the visible parcels within 2–5 %
  (c01_massonly 3.07e-2 against 3.12e-2 m/s). Pixel size 0.2 mm, 128 px
  windows and frame pairs 4–10 ms apart all stay within ±5 %.
- The wall layer is not special: at y = 16–20 mm it descends at 0.027 m/s,
  the next three parcel layers at 0.033 m/s. All wall parcels, the visible
  parcels and the 4 mm box give 2.0–3.2e-2 m/s (mass-only c01). C26 Fig. 8c
  also has PIV ≈ cumulative parcel mean at W = 10.
- Box height is a 2–3× effect. At C26's Fig. 1c box height (~0.65 of the
  bed) we are 6.8× (mass-only) and ~3× (momentum) C26's model. Less than
  1.2 % of our visible wall is as slow as C26's Fig. 8e baseline map.
- No start-up drift: 0.5 s block means over 1–3 s stay at 1.7–4.4e-2 m/s.
- Static image content can make PIV read low (glare blended at 60 %: −92 %;
  90 % of the wall layer frozen: −83 %). That could only affect C26's
  experiment, not either model.

**What sets it: circulation.** The wall flow is the annulus return flow.
Upward solids flux is 24 g/s at y = 18 mm and 40–46 g/s at 24 mm
(mass-only); C26's 3.3e-3 m/s implies ~3 g/s by the same mass balance. The
whole bed heaves at ~11 Hz (centre-of-mass ±5–8 cm/s), and the box signal
follows it (r = 0.78–0.84), 3–5× the amplitude of C26's Fig. 6 model trace.

**Leading suspect: the gas inlet cells.** Neither source reproduces C26's
spout (Fig. 7: ~0.53 m/s peak, ~4–5 mm wide, read by eye, probe height not
stated):

| c01, 10 slpm | spout core v_p | spout FWHM, y = 20 mm | wall D / C26 model |
|---|---|---|---|
| momentum source | 1.75 m/s | ~4 mm | 3.4 |
| mass-only source | 0.42 m/s | ~8 mm | 6.7 |

- Momentum: core 3.3× too fast, wall 3.4× too fast, the same ratio.
- Mass-only: ~2/3 of the gas leaves the source cell sideways into the bottom
  of the annulus ("Spout off-centre"); ΔP 2.6–2.8 W/A against 1.66 W/A.
- Second suspect, not separable from this one without a W test: C26's own
  wall velocity rises 3–5× from W = 10 to W = 1 (Fig. 8c), and our W = 10
  momentum run (1.1e-2 m/s) sits on their W = 1 value. MFiX 26.1.2 computes
  parcel drag from the real particle diameter (`des/drag_gp_des.f:114`,
  `DES_CGP_RPR`); 23.2 checked indirectly on 2026-10-22 (very likely the
  same; see "CGP drag check").

**Why the grid layout.** C26's 24 cells across D put the axis on the edge
shared by the four central cells; we use 27 cells so the source is one
centred cell. In our layout the source cell's side neighbours start at
|x| = 1.05 mm, inside the cone already at the cap (R_t = 1.19 mm), so its
side faces are open. In C26's layout the next cells out start at 2.1 mm,
outside the cone in the first row. Open face area of the first-row source
region, straight cone (scratch `source_faces.py`):

| layout | sideways share of the open face area |
|---|---|
| ours (1 cell), flat cap | 79 % (MFiX's cut faces: 75 %, 13.4 of 17.8 mm²) |
| C26 (2 × 2 cells), flat cap | 17 % (top opening ~16 mm²) |
| C26 (2 × 2 cells), cone closed at y = 0 | 0 % (top opening 4.6 mm² ≈ nozzle 4.45 mm²) |

- With our flat cap, C26's layout sends the gas up at ~10–12 m/s: between
  our two sources, a case not yet run. If C26's cone closes to a point, the
  same layout gives a nozzle-speed jet, i.e. our momentum run (3.4×).
- Gas can still turn sideways one row up; only a run shows how much.

**The A/B case.** `run_C26_3d --case c01 --axis-on-face --run-dir-suffix
_massonly_axisface` (geometry.yaml `axis_on_cell_face`). Against
`c01_massonly`, the only non-comment deck changes are the box (26 cells,
x, z = ±27.3 mm), the PS region (±2.1 mm in x and z) and monitors 1 and 4
over it. dx, the STL (identical md5), the cone, the gas and all DEM
settings are unchanged. DMP default `[13, 13]` in i and k: the rank
boundaries pass through the axis, and MFiX sums the PS volume over ranks
(`set_ps.f` l.116).
- **Dry run** (local, 4 ranks, to t = 0.042 s): no errors; seeded 27,778
  parcels, 60.97 g (c01_massonly 27,782, 60.97 g); PS volume 1.196e-8 m³
  over 4 cells (axis cell 7.0e-9).

**Criteria, written before any A/B result.**
- *G, gas-only, 0.05 s, c01 (10 slpm), both layouts* (deck changes as the
  2026-10-20 gas-only tests: no solids, BACKGROUND every 5 ms). Metrics over
  t = 0.01–0.05 s from `notes/c26_source_split.py`:
  - **G1 up share** = net gas flux out of the source cells' north faces / Q.
    Mechanism confirmed if the new layout gives ≥ 0.80 and the old ≤ 0.50.
    If the new layout gives < 0.60, the premise fails: stop before the bed
    smoke.
  - **G2 closure**: signed flux through the whole first-row top plane / Q
    within 1 ± 0.05; otherwise the approximate cut-face areas are not
    trusted and G1 is reported, not scored.
  - **G3 DMP**: new layout, 4 ranks against serial, up share within 0.02.
- *Bed smoke, 0.3 s local, 4 ranks*: `notes/smoke_validate_c26.py` GO (P3
  spout required for c01) before the 3 s AWS run.
- *3 s A/B against c01_massonly*, t = [2.4, 3.0] s, same scripts
  (`notes/c26_step3_checks.py`, synthetic PIV):
  - wall D: **substantial** if it falls ≥ 30 % from 2.23e-2 m/s, **no
    effect** if it changes < 15 %, otherwise inconclusive;
  - reported with it: spout axis v_p and FWHM at y = 20 mm, the source-gas
    up share (BACKGROUND, same script), ΔP/(W/A) and its PSD peak, the PIV
    wall box, and C26's bands (D 0.3e-3–10e-3 m/s, C1 0.26–1.06 m/s).
  - Reading fixed in advance: substantial → the inlet cells are a cause of
    the gap; no effect → the source-cell leak does not set the wall
    velocity, and coarse-graining becomes the leading suspect.
- **G not run** (user decision, 2026-10-20): the three gas-only runs were
  stopped during start-up, before any output was evaluated. The up share is
  measured on the bed runs instead (`notes/c26_source_split.py` on their
  BACKGROUND frames).
- **Plan** (user decision, 2026-10-20): if the bed smoke is GO, queue all
  three flows on AWS: `c01/c02/c03_massonly_axisface` (10 / 12 / 14 slpm),
  each generated with `--axis-on-face`. Against `c0X_massonly`, each deck
  changes only the box, the PS region and monitors 1 and 4; ps_massflow_g
  (2.722e-4 / 3.266e-4 / 3.810e-4 kg/s) and dt_max (1e-5 s) are unchanged.
  The c02/c03 A/Bs use the same rule as c01, against their own mass-only D
  (2.92e-2 / 3.39e-2 m/s).
- **Bed smoke `c01_massonly_axisface_smoke`: GO** (0.3 s, 4 ranks, 1.31 h;
  `notes/smoke_validate_c26.py`, mass-only rules of 2026-10-19). Old-layout
  smoke `c01_smoke_r4_massonly` in brackets; metrics over t ≥ 0.1 s.
  - check_smoke: S0–S2, S9 PASS (dt 1e-5 fixed, NIT median 2). S3 exit
    v_g/U_nozzle 0.28 [0.24–0.39], reported, not gated. S6 ΔP 2.87 W/A
    [2.70]; S7 sd 1.6 W/A; S8 PSD peak 10 Hz [10 Hz].
  - P1 PS 1.196e-8 m³ over 4 cells; P2 27,778 parcels, 0 lost; P3i 42.3 ṁ_g
    [37.8]; P3ii 0.17 mm [0.18]; P4 3.68 m/s [3.2], reported.
  - Source split (`c26_source_split.py`, BACKGROUND t = 0.1–0.3 s, 5
    frames). Up share **1.02 [0.38]**; closure 1.02 [1.08]; cross-axis flux
    −0.008 / −0.005 Q. Upward-flux centroid in rows 2–9 within 0.1 mm of the
    axis [x −0.2 to −0.6 mm]. The source gas now leaves upward, as the
    open-area estimate predicted (G1's ≥ 0.80 threshold, here on the bed
    run, not gas-only).
  - The smoke is a gate, not the A/B: wall D is judged on the 3 s runs.
- **Submitted 2026-10-20** (AWS Batch, 4 vCPU, 6000 MB, 4 ranks, tstop 3 s):
  `c26-c01-massonly-axisface` 4021f349-f487-4cca-b23c-fd30902da70c,
  `c26-c02-massonly-axisface` 413be729-5f0d-46bd-936e-c0f515ac2f09,
  `c26-c03-massonly-axisface` 1540440e-b526-4985-9768-96cae6f121d2.
- **c04–c09 requeued on the new layout** (user decision 2026-10-21; the old
  c04–c09_massonly 3 s runs all reached 3.0 s on AWS and are kept for
  comparison). Decks from `--axis-on-face --run-dir-suffix
  _massonly_axisface`; only the box, PS and monitors 1/4 differ. 1 ms 4-rank
  local dry runs PASS; no 0.3 s smoke (the layout smoked GO on c01).
  Submitted in this order, c07 last because it does not spout at 12 slpm:
  c04 11c096db-634b-48a4-94a7-1ac416adac52, c05
  5b960762-7415-4c08-be54-ca630b0ae8a9, c06
  ada58098-bae1-4d18-aa0c-2a45f7ef95c0, c08
  89e2191c-2860-4a6f-976c-d327a00ab652, c09
  e8969c1f-085f-4f25-a660-c888cf2ed311, c07
  083b9d06-2742-4d46-819f-cc3302cef224.
- **Dry-run pitfall.** Generated `run.sh` files `cd` to the absolute case
  directory. A copied case run through its `run.sh` therefore runs the
  production deck in the repo, not the copy. That happened on 2026-10-21:
  c04 ran to 0.5 s locally before it was caught (logged as operator error).
  Before running a copy through `run.sh`, change the copy's `cd` line to
  point at the copy.

**Scoring the 3 s axis-face runs** (written 2026-10-21, before any axis-face
3 s result). `notes/c26_watch_post.py` polls the Batch jobs. When a job
succeeds it pulls `aws_out`, runs health, step 3, synthetic PIV and the
source split, and does the same for the old-layout run. It writes
`results/c26_compare_<tag>.csv`. Thresholds are in
`params/c26_acceptance.yaml`:
- **Step 3:** A–D as pre-registered on 2026-10-19. C1 and D gate c01 only;
  for c02 and c03 they are a reference.
- **Grid-layout A/B:** the 2026-10-20 rule (D ≥ 30 % lower → substantial,
  < 15 % change → no effect).
- **PIV map, new on 2026-10-21:** the map "matches C26's model" if CC and
  SSI are at least the lower, and MSE at most the higher, of C26's two own
  model scores (Table 3, and Fig. 4b against 4a through this pipeline).
  C26 states no acceptance level, so this floor is our assumption: our map
  must score at least as well as C26's own model.
- **PIV wall:** reported as multiples of C26's experiment and model
  (Fig. 6d). For c01 the step-3 D band gates.

## Axis-face 3 s results (2026-10-21)

`c01/c02/c03_massonly_axisface` (10 / 12 / 14 slpm) all reached t = 3.000 s
on AWS (4 ranks; 22.4 / 23.5 / 23.2 h). They were scored by
`notes/c26_watch_post.py`, giving `results/c26_compare_axisface.csv` and
`results/c26_piv_c0X_massonly_axisface.csv` / `fig_*.png`. Step 3 and the
PIV use the window t = [2.4, 3.0] s (299–300 frames); health uses
t ≥ 0.1 s. Old-layout `c0X_massonly` values are in brackets.

**Not validated.** A, B, C1 and C2 pass at all three flows. The wall
velocity (D) and the PIV map do not match C26 at any flow.

**Health.**
- c01 and c02 are GO. c03 is NO-GO on check_smoke S6 only: mean ΔP is
  3.09 W/A, against the band's 3.0 upper bound.
- ΔP is 2.76 / 2.92 / 3.09 W/A [2.58 / 2.65 / 2.77]. It rose by 0.18–0.32 W/A
  at each flow.
- Cause not established. S6 is meant to catch a clogged or blocked inlet,
  and none of that run's signs are present:
  - dt fixed at 1e-5 s, NIT median 2;
  - S7 sd 2.4 W/A, under the 3.0 flag;
  - 0 parcels lost;
  - spout centred (P3ii 0.34 mm).
- The layout's main change at the inlet is that all the source gas now goes
  up into the spout. The gate is not relaxed; c03's numbers below carry the
  NO-GO.
- Spout-flux centroid (P3ii) 0.20 / 0.36 / 0.34 mm [0.84 / 0.99 / 1.18].
  The −x lean of the mass-only runs is gone. The gas upward-flux centroid
  (c01, y = 5–20 mm) is within 0.2 mm of the axis.

**Step 3** (thresholds of 2026-10-19; C1 and D gate c01 only):

| criterion | threshold | c01 | c02 | c03 |
|---|---|---|---|---|
| A fountain | ≥ 50 % of frames | 88 % PASS [45 % FAIL] | 100 % PASS [96 %] | 100 % PASS [100 %] |
| B ΔP PSD peak | 10–20 Hz | 11.7 PASS [11.7] | 11.7 PASS [10.0] | 10.0 PASS [10.0] |
| C1 peak axis v_p | 0.26–1.06 m/s | 0.61 (y 10 mm) PASS [0.42] | 0.67 (y 20 mm) ref [0.47] | 0.77 (y 20 mm) ref [0.53] |
| C2 outer/axis, worst height | < 0.1 | 0.060 PASS [0.148 FAIL] | 0.076 PASS [0.182 FAIL] | 0.092 PASS [0.188 FAIL] |
| D wall, downward | 0.3e-3–10e-3 m/s | 1.78e-2 FAIL [2.23e-2] | 2.81e-2 ref, above [2.92e-2] | 3.40e-2 ref, above [3.39e-2] |

- **Spout at y = 20 mm.**
  - Axis v_p 0.59 / 0.67 / 0.77 m/s [0.42 / 0.47 / 0.49].
  - FWHM 6.3 / 6.7 / 6.8 mm [7.3 / 7.9 / 5.7].
  - C26 Fig. 7 (by eye): ~0.53 m/s peak, ~4–5 mm wide.
- **Source-gas up share** 1.01 at all three flows [0.36 / 0.35 / 0.34].

**Grid-layout A/B on D** (pre-registered 2026-10-20):
- c01: −20 % → inconclusive.
- c02: −4 % → no effect.
- c03: +0.3 % → no effect.
- The PIV wall box fell 32 / 10 / 0 % [3.07 / 3.56 / 3.85e-2 m/s]. The rule
  is on D, so this is reported only.
- Reading fixed in advance: "no effect" means the source-cell leak does not
  set the wall velocity, and coarse-graining becomes the leading suspect.
  The layout fixed the spout (direction, centring, C2, c01's fountain) but
  not the wall velocity.

**Synthetic PIV**, t = [2.4, 3.0] s, default settings. Map floor = the worse
of C26's two own-model scores (rule of 2026-10-21).

| run | CC | SSI | MSE [m²/s²] | floor CC / SSI / MSE | map | wall PIV [m/s] | × exp | × C26 model |
|---|---|---|---|---|---|---|---|---|
| c01 | −0.107 [−0.134] | 0.129 [0.125] | 0.0172 [0.0171] | 0.390 / 0.609 / 0.0068 | NO MATCH | 2.08e-2 [3.07e-2] | 22 | 6.3 |
| c02 | 0.186 [0.130] | 0.162 [0.132] | 0.0327 [0.0348] | 0.817 / 0.740 / 0.0042 | NO MATCH | 3.20e-2 [3.56e-2] | 17.5 | 6.3 |
| c03 | 0.035 [−0.067] | 0.264 [0.169] | 0.0247 [0.0283] | 0.829 / 0.741 / 0.0049 | NO MATCH | 3.84e-2 [3.85e-2] | 15.1 | 5.5 |

- Map means 0.024 / 0.028 / 0.042 m/s [0.015 / 0.021 / 0.024], against
  0.11 / 0.18 / 0.16 in the experiment. Empty windows: 85 / 83 / 77 %
  [88 / 86 / 85 %].
- The fountain is still a low dome, ~5–10 mm above the bed surface. It has
  no side lobes of falling particles.
- The wall box still pulses at ~11 Hz, between about +0.08 and −0.015 m/s
  (c02).

## CGP drag check: MFiX 23.2 vs 26.1.2 (2026-10-22)

**Question.** After the grid-layout A/B, coarse-graining became the leading
suspect for the wall velocity. The planned first step, before any W run, was
to check whether MFiX 23.2 (C26) computes CGP parcel drag the same way as
26.1.2 (ours). If 23.2 evaluated the drag law at the parcel diameter, its
W = 10 parcels would couple to the gas very differently from ours.

**Method.** `notes/c26_cgp_drag_check.py` (log in scratch):
- Part 1 searches the 26.1.2 source behind `run.sh` for the lines that set
  CGP drag, parcel size and volume, and contact.
- Part 2 evaluates TGS (`src/correlations.py::tgs_drag_coefficient`,
  `umf_tgs`, unit-tested) with the c01 deck's gas properties and W, at both
  diameters.

**26.1.2, from the source** (all checks PASS):
- **Drag uses the real diameter.** The drag law is evaluated at the real
  diameter, `IF(CGDEM) DPM = 2.0d0*DES_CGP_RPR(NP)` (`des/drag_gp_des.f:114`).
  DES_CGP_RPR is the parcel radius / W^(1/3) (`generate_particles_mod.f:1391`).
- **Then multiplied by the parcel volume.** `F_gp = DgA * PVOL`
  (`drag_gp_des.f:254`), with PVOL and PMASS taken from the parcel radius
  (`make_arrays_des.f:311`). Drag per unit parcel volume therefore equals
  that on one real particle, so a parcel is accelerated by the gas like a
  real particle.
- **Our TGS mirrors MFiX's.** For one solids phase DPA = DPM, so TGS's
  argument order (`drag_gp_des.f:222`) cannot matter. MFiX's TGS uses
  `RE = DPA*VREL*ROPg/Mug` and `F_STOKES = 18D0*Mug*EPg/DPM**2`
  (`drag_gs.f:2193`, `2203`); `src/correlations.py` mirrors both.
- **Contacts have no CGP term.** `calc_force_dem.f`,
  `calc_collision_wall_mod.f`, `cfwallcontact.f` and `cffctowall.f` contain
  none. Contacts use the parcel's size and mass, with kn, e and μ as input
  (kn not scaled with W).

**23.2.** No direct diff was possible:
- The source is not reachable here. The conda channel recorded in the
  26.1.2 install returns HTTP 401, even for the installed 26.1.2 package.
- NETL's GitLab exposes no MFiX project publicly.
- The GitHub copies are 2016–2017 MFiX or MFiX-Exa.

The indirect evidence all points the same way:
- **Forum, before 23.2.** MFiX staff, June 2022 (MFiX 22.x, forum topic
  3900, "Source code of mfix-cgp"), quote the same line: "the drag is using
  the real (not scaled) particle diameter: drag_gp_des.f: IF(CGDEM) DPM =
  2.0d0*DES_CGP_RPR(NP)".
- **Forum, after 23.2.** A user found the same line in 25.x (topic 6408,
  2025-10-21).
- **Release notes 23.3–26.1.** They list no change to CGP drag. CGP items:
  - 23.4 GUI statistical-weight fix;
  - 24.2 tangential wall velocity;
  - 24.3 CSV particle input;
  - 24.4 solids time step computed at runtime;
  - 25.2 CGDEM mass-inflow size bug;
  - 25.3 wide-PSD handling.

  The only drag fix is SQP drag (26.1). No 25.4 announcement was found.
- **23.2's own notes.** They include "Fixed call to TGS drag law". That fix
  is in the version C26 ran, and for a single phase the DPM/DPA order is
  immaterial.
- **Docs.** The 23.2 and 26.1 descriptions of MFiX-CGP are identical except
  for capitalization (`about.rst`).
- **C26's own model behaviour.** The parcel diameter would make a gross
  difference (table below), and C26's model spouts near 10 slpm, as ours does.

What the parcel diameter would change (TGS, Ar at 298 K, YSZ 411.6 µm,
W = 10 → parcel 886.8 µm; voidages assumed, `eps_mf_range` in
`params/particles.yaml`):

| ε | U_mf, real d | U_mf, parcel d | ratio | DgA real / parcel at U_mf |
|---|---|---|---|---|
| 0.40 | 0.104 m/s | 0.385 m/s | 3.7 | 4.4 |
| 0.42 | 0.129 m/s | 0.462 m/s | 3.6 | 4.4 |
| 0.45 | 0.175 m/s | 0.595 m/s | 3.4 | 4.3 |

- Single-particle terminal velocity: 4.53 m/s real, 8.88 m/s parcel.
- The c01 bed-top superficial velocity at 10 slpm is 0.117 m/s, i.e.
  U_mf-scale with real-diameter drag and ~1/4 of U_mf with parcel-diameter
  drag.

**Verdict.**
- CGP drag is evaluated at the real diameter in 26.1.2. It very likely is
  in 23.2 too. This rests on the evidence above, not on a diff.
- Contact handling has no CGP-specific code in 26.1.2. Nothing in the
  release notes suggests 23.2 differed.

**Reading.**
- With the same CGP scaling, W, drag law and contact inputs (kn 1000 N/m,
  e 0.5, μ 0.1, C26 Table 1), coarse-graining can explain a gap between a
  W = 10 model and the experiment (C26 Fig. 8c). It cannot by itself explain
  why our W = 10 wall velocity is 5.5–6.3× C26's W = 10 model.
- That gap is in something C26 and we do differently:
  - unstated inputs: point-source momentum, apex geometry, des_epg_clip and
    interpolation, tangential/damping settings, the solids time-step factor;
  - non-CGP version changes: 24.4/25.3 solids time step, 25.1 wall damping
    (reviewed, ±4 %);
  - or the definition of C26's wall velocity.
- C26 also reports its model's wall particles pulsing between positive and
  negative peaks, absent in the experiment (p.7). Ours do too.
- A W < 10 run would not separate the two models and is outside the C26
  scope (W ≥ 10).

**To make it a direct diff:** refresh the MFiX download link (registration
e-mail) and run `diff` on `des/drag_gp_des.f`, `drag_gs.f`
(DRAG_TGS), `make_arrays_des.f`, `generate_particles_mod.f` and
`calc_force_dem.f` between 23.2 and 26.1.2.

## Coating-stage / charge 3 s results (2026-10-22)

**Runs.** All reached t = 3.000 s on AWS (exit 0, 4 ranks, 14.5–23.4 h):
- c04–c09_massonly (old layout);
- c04, c05, c06 and c08_massonly_axisface.

c07 and c09_massonly_axisface were still running when this was written.

**Outputs.**
- Pulled with `aws s3 sync` into `cases/C26_3d/<case>/aws_out`.
- Scored with `notes/c26_watch_post.py --local`, giving
  `results/c26_compare_axisface.csv` (c01–c06, c08) and
  `results/c26_compare_massonly_sweep.csv`.
- Summarized with `notes/c26_sweep_summary.py`, giving
  `results/c26_sweep_summary{,_axisface}.csv` and
  `fig_c26_sweep_{flow,diameter,charge}{,_axisface}.png`.

Step 3, the PIV and the sweep use t = [2.4, 3.0] s; health uses t ≥ 0.1 s.

**What can be compared with C26.**
- C26 measured only YSZ at 60 g, at 10, 12 and 14 slpm (Table 3, Fig. 4a,
  Fig. 6d). c04–c09 have no C26 map or wall target.
  - `src/post/piv_c26.py` leaves their map scores NaN.
  - `notes/c26_watch_post.py` now labels the map "n/a (no C26 data for
    this case)". Before this change it read "NO MATCH" against NaN
    floors, and `read_piv` left the NaN as a string, so `piv_match` could
    not compare it.
  - A, B and C2 (C26's YSZ observations) are now labelled "(ref)" for
    these cases, as C1 and D already were.
- What is left to compare:
  - the C26 12 slpm YSZ 60 g values (wall 1.83e-3 m/s experiment, 5.05e-3
    model; Fig. 4a map mean ~0.18 m/s) as a reference for the 12 slpm
    sweep;
  - C26's spout-regime observations (fountain; 10–20 Hz pulsation).
- The sweep result itself is relative to c02.

**Results, axis-face [old layout].** C26 12 slpm YSZ 60 g reference:
experiment 1.83e-3 m/s, model 5.05e-3 m/s.

| case | particle, charge | health | spout | axis v_p y = 10/20/30 mm [m/s] | C2 | fountain [% parcels] | ΔP [Pa] (W/A) | PSD peak [Hz] | wall PIV [m/s] | wall ÷ C26 exp / model |
|---|---|---|---|---|---|---|---|---|---|---|
| c02 | YSZ 60 g | GO [GO] | yes | 0.65/0.67/0.61 [0.41/0.47/0.45] | 0.076 [0.182] | 1.74 [1.21] | 1202 ± 690 (2.91) [1049 ± 684] | 11.7 [10.0] | 3.20e-2 [3.56e-2] | 17.5 / 6.3 [19 / 7.0] |
| c04 | + buffer 612 µm, 60 g | GO [GO] | yes; weak fountain [off-centre, P3ii 2.58 mm] | 0.59/0.51/0.40 [0.16/0.19/0.20] | 0.061 [0.363] | 0.08 [0.09] | 675 ± 48 (2.29) [592 ± 50] | 13.3 [13.3] | 1.57e-2 [1.73e-2] | 8.6 / 3.1 [9.5 / 3.4] |
| c05 | + IPyC 692 µm, 60 g | **NO-GO** P2 [**NO-GO** P2] | internal jet only | 0.29/0.08/0.04 [0.06/0.04/0.03] | 0.185 [0.900] | 0.84 [0.80] | 617 ± 24 (2.09) [542 ± 23] | 23.3 [11.7] | 8.48e-3 [5.11e-3] | 4.6 / 1.7 [2.8 / 1.0] |
| c06 | + SiC 762 µm, 60 g | **NO-GO** P2 [**NO-GO** P2] | no; bed expanded | 0.04/0.02/0.02 [0.01/0.02/0.04] | 0.406 [0.729] | 8.04 [8.40] | 684 ± 45 (2.32) [618 ± 48] | 11.7 [15.0] | 8.09e-3 [9.08e-3] | 4.4 / 1.6 [5.0 / 1.8] |
| c07 | + OPyC 842 µm, 60 g | running [GO] | [no; packed bed] | [0.11/0.01/0.00] | [n/a] | [0.00] | [596 ± 22] | [15.0] | [8.2e-6] | [~0] |
| c08 | YSZ 40 g (39.6 seeded) | GO [GO] | yes | 0.66/0.66/0.63 [0.38/0.39/0.52] | 0.092 [0.199] | 4.13 [3.87] | 1044 ± 891 (2.90) [969 ± 929] | 11.6 [10.0] | 3.40e-2 [3.84e-2] | 18.6 / 6.7 [21 / 7.6] |
| c09 | YSZ 80 g (76.0 seeded) | running [GO] | [yes, off-centre: P3ii 1.83 mm] | [0.36/0.41/0.42] | [0.168] | [0.45] | [1145 ± 444] | [11.7] | [2.92e-2] | [16 / 5.8] |

**Charge block (YSZ, 12 slpm).**
- c08 (40 g) spouts strongly in both layouts. On the axis-face layout, A
  is 100 %, C1 0.66 m/s, C2 0.092 and mid-bed v_g 6.2 m/s.
- c09 (80 g, old layout) spouts: P3i 37.7 ṁ_g, A 89 %. Its spout is
  off-centre (P3ii 1.83 mm), the old layout's known −x lean, so health
  labels it "NO SPOUT" by the P3ii rule.
- Trends:
  - the fountain share falls with charge (3.9 → 1.2 → 0.45 %, old layout);
  - ΔP rises with charge (969 → 1049 → 1145 Pa);
  - the PSD peak stays at 10–11.7 Hz.
- **The wall velocity is insensitive to charge.** It is 16–21× C26's
  12 slpm experiment and 5.8–7.6× its model at 40, 60 and 80 g, the same
  gap as c02. This adds to the evidence that the gap is not a bed-depth
  effect.

**Diameter block (60 g, 12 slpm).** At a fixed 12 slpm the coated beds
stop spouting.
- **c04 (buffer)** spouts on the axis-face layout. Its fountain is weak:
  A 16 %, 0.08 % of parcels.
- **c05 (IPyC)** has only an internal jet: the axis v_p decays from
  0.29 m/s at 10 mm to 0.04 m/s at 30 mm. Mid-bed v_g is 0.86 m/s, about
  the uniform-percolation value. A passes on ~120 parcels above the
  threshold, which reads as surface bubbling, not a spout fountain. Its
  axis-face PSD peak is 23.3 Hz, outside the 10–20 Hz band.
- **c07 (full TRISO, old layout)** is a packed bed with no motion. Its D
  sits in band only because the bed is static.
- ΔP falls from about 1050–1200 Pa ± 690 (bare, pulsating) to 540–680 Pa
  ± 22–50 (steady).
- The PSD peak of a steady ΔP (c06, c07) is noise, not a spouting
  frequency.
- The lower wall velocity of c05 and c06 (2.8–5.0× C26's 12 slpm
  experiment, 1.0–1.8× its model) comes from a bed that barely circulates. It is **not** better
  agreement with C26.
- **Confound** (design note above): at 60 g the coated beds are 13–15 mm
  deeper (47.7–49.7 mm, above the 41.6 mm cone top). So the loss of the
  spout is not attributable to diameter alone. TGS U_mf (ε 0.42) is 0.129
  m/s for bare YSZ against 0.118 / 0.136 / 0.174 m/s for c04 / c05 / c06,
  so the minimum-fluidization velocity alone does not separate c04 from
  c02.

**c06 anomaly (unexplained).**
- In both layouts the PIV bed surface sits at ~61 mm against H_static
  47.7 mm, and the parcels thin out up to 80 mm. The bed is expanded.
- Its A, fountain share (8 %) and PIV map count that expanded bed, not a
  fountain.
- This is physically unexpected: c06's U_mf (0.174 m/s) is above c05's and
  c04's, and neither of those expands.
- c06 also loses the most parcels (35–46; every traceable one via the top
  PO).
- Same in both layouts, so not the source layout. Candidates: voidage
  interpolation and des_epg_clip with 1.64 mm parcels in 2.1 mm cells
  (dx/d_parcel 1.28). Not tested.

**Health P2 (c05, c06, both layouts): NO-GO, not relaxed.**
- Parcels lost: 6–46 (0.04–0.44 % of the mass, under the 0.5 % outlet
  limit).
- P2 also requires every lost parcel to be traced to the top PO. Tracing
  uses consecutive PARTICLES_PIV frames, which exist only for t ≥ 2.4 s.
- Inside the window every loss is traced to the top PO: c05 axis-face 2/2,
  c06 axis-face 9/9, c06 old 5/5, c05 old 0.
- Before the window, 6–37 losses per run cannot be traced. The loss rate
  is about the same before and during the window (c06 axis-face ~15/s in
  both).
- The PIV-frame parcel counts equal the inventory monitor, and no parcel
  centre is beyond the cone wall in any sampled frame. So these are real
  removals, most likely ballistic exits through the PO, but unproven
  before t = 2.4 s.
- To make P2 decisive for 3 s runs, write a sparse particle VTK
  (positions and ids) over the whole run. Proposed, not applied.

**Grid layout at 12 slpm** (axis-face vs old; no pre-registered A/B for
c04–c08, so reported only):
- The layout centres the spout: P3ii 0.91 / 0.13 mm for c04 / c08, against
  2.58 / 0.85 mm.
- It raises the spout's axis v_p and brings C2 under 0.1 where the old
  layout fails (c04 0.363 → 0.061, c08 0.199 → 0.092).
- D changes by −5 / +24 / −14 / −12 % (c04 / c05 / c06 / c08). As at
  10–14 slpm, the layout fixes the spout, not the wall velocity.

## First-principles review of the velocity gap; momentum24 runs (2026-10-22)

Asked for on 2026-10-22: explain, from the results in `results/` and the
paper alone (without relying on the conclusions above), why our particle
velocities disagree with C26's model and experiment. The starting hypothesis
was that fountain parcels fall faster than in the experiment because a
retarding force is missing on the way down. Scratch scripts (parcel
kinematics by ID across the PIV frames, BACKGROUND cell-centre binning,
spout flux share against height) were run on the 3 s outputs of c01
(momentum, 27 cells), c01–c03_massonly and c01–c03_massonly_axisface. No new
runs were needed for the review.

**The hypothesis is not supported.**
- The terminal velocity of a 412 µm YSZ sphere in air is ≈ 3.3 m/s, so gas
  drag on a fountain parcel moving at 0.3–1 m/s is ≤ 5 % of its weight. No
  gas-phase force the model could be missing is large enough to matter in
  the fountain. MFiX 26.1.2 computes CGP drag with the real particle diameter
  and scales the force by the parcel volume, so a parcel's acceleration
  equals a single particle's ("CGP drag check").
- Tracked by ID, descending fountain parcels decelerate relative to free
  fall by 16–32 % of g in every run (collisions with rising parcels and with
  the bed surface). They are not accelerating too much. In the mass-only
  runs the fountain is too *weak*, not too fast.
- Where the model is too fast is the annulus wall layer (D, 6–30× the
  experiment, in every run) and, in the momentum run only, the spout core
  (C1).

**Two regimes, neither of them C26's.**

| | c0X_massonly (27 and 26 cells) | c01 momentum (27 cells) | C26 model / experiment |
|---|---|---|---|
| regime | bubbling, heaving bed: surface heave ±0.1 m/s at ~11.6 Hz; 24 % of bed cells at EP_G > 0.7 | spouting, 18.3 Hz pulsation | spouting, ~20 Hz |
| mean ΔP | ≈ 2.9 W/A | 1.66 W/A | – |
| gas speed in the source cells | ~11–13 m/s | 30 m/s on the axis at y = 5 mm | nozzle 45 m/s (12 slpm) |
| axis parcel v at y = 10 mm | 0.42 m/s | 1.75 m/s | ~0.5 m/s (Fig. 7, model) |
| fountain top (median) | a few mm; 83–88 % of PIV windows empty | 55 mm | ~20 mm |
| PIV map mean, 10/12/14 slpm | 0.024 / 0.028 / 0.042 m/s | 0.134 (10 slpm) | expt 0.111 / 0.178 / 0.158; model 0.132 / 0.182 / 0.190 |
| wall velocity | 3–4e-2 m/s | 1.24e-2 m/s | model 3.3e-3; expt 0.93e-3; model at W = 1 1.39e-2 |

- The mass-only source releases the 12 slpm into the fluid volume of its
  cells with no momentum, and the cut-cell source region leaks part of it
  sideways. The jet that forms is set by the open face area, ~11–13 m/s
  instead of the 45 m/s nozzle speed, i.e. 1/3–1/4 of the momentum flux.
  That drives a bubbling bed, not a spouted one, so the fountain is weak and
  the time-averaged maps are dominated by empty windows.
- c01's momentum map mean matched C26's by coincidence: concentrated lobes at
  0.35–0.40 m/s over a cold background, against C26's diffuse 0.2–0.3 m/s
  fountain (`results/fig_c26_piv_c01.png`). The earlier statement that the
  momentum source "restores the fountain" was too generous.
- The momentum run over-drives the spout because the imposed ṁU stays in a
  1–2 cell column. A 2.1 mm grid cannot resolve the spout base, so the source
  footprint (one cell on the 27-cell layout, 2 × 2 on 24/26) is the free
  parameter that sets how fast the jet spreads ("smearing"). Nozzle momentum
  spread over the 2 × 2 footprint that C26's own grid implies is the untested
  middle case.
- Wall velocity: the net annulus descent is confined to a ≤ 4 mm layer on
  the 60° smooth wall (μ_w = 0.1, no rolling resistance, gas ≈ 0.85 U_mf
  next to the wall). Nothing in the model holds that layer back. Part of the
  gap to the experiment is probably measurement physics (a tribocharged skin
  of particles adhered to the glass, which PIV sees as slow); the paper's own
  W = 1 model value (1.39e-2 m/s) is close to our momentum run.

**Ranked causes (best guesses).**
1. Mass-only source + 2 × 2 footprint + cut-cell side leakage gives a
   1/3–1/4 momentum-flux jet and a bubbling regime.
2. Wall-layer return flow with nothing to restrain it: wall friction,
   rolling resistance, wall roughness and gas bypass at the wall are all
   candidates; none is set by the paper.
3. Measurement physics in the experiment's 4 mm wall box.
4. Grid resolution of the spout base, which turns the source footprint into
   a tuning parameter.

**Layout change (user decision 2026-10-22).** `params/geometry.yaml`
`n_margin_cells_xz: 1 → 0` with `axis_on_cell_face: true`, so the box is
exactly 24 cells = D_c as C26 p.3 states. Against the 26-cell axis-face
layout the fluid cells (42,520), cut cells (8,608), STL (md5 93f1ea18,
77,088 facets) and the 2 × 2 source are identical; only the blocked margin,
the box extents (x, z = ±25.2 mm) and the DMP gridmap ([12,12]) differ.
There is no physics difference; the change is for consistency with the
paper. `src/make_case_c26.py::default_gridmap` and
`tests/test_make_case_c26.py` were updated (77 tests pass). The odd layout is
still available with `axis_on_cell_face: false` (25 cells at margin 0).

**momentum24 runs.** Nine 3 s runs with `--ps-momentum` on the 24-cell
layout: `cases/C26_3d/c0X_momentum24`, jobs `c26-c0X-momentum24`, 4 ranks
(2×1×2), 6000 MB (c09 8000 MB), submitted 2026-10-22 after local 1 ms dry
runs (`results/run_log.csv` rows carry the job ids). `ps_v_g` = 37.42 /
44.91 / 52.39 m/s for c01 / c02 / c03 and 44.91 m/s (12 slpm) for c04–c09;
everything else is the c0X_massonly_axisface deck on the 24-cell box. c07
(full TRISO, may not spout at 12 slpm) was submitted last. c09 seeds 34,636
parcels = 76.0 g, the known lattice shortfall against 80 g.
- Expectation, written before any result: step-3 and PIV metrics between
  the 27-cell momentum c01 and the mass-only runs; a spout at 12 slpm; a
  fountain top below c01's 55 mm; and D below 1.24e-2 m/s if the footprint
  matters. If D stays at 1.2e-2 m/s the wall layer, not the source, is the
  next thing to work on.
- Status 2026-10-23 11:30 UTC: c05 and c07 reached t = 3.000 s (16.3 h and
  14.3 h wall, 300 PIV frames, results on S3 under `results/c0X_momentum24/`);
  c01–c04, c06 and c09 were running (18 h in); c08 started at ~08:00 UTC when
  c07 freed a slot. Nothing pulled or scored yet. Score with
  `notes/c26_watch_post.py` and `src/post/piv_c26.py`, compare c01–c03 with
  C26 Table 3 and Fig. 6d next to the `_massonly_axisface` rows, and add a
  24-cell column to the regime table above.

**Tooling change (2026-10-22).** The AWS submit, poll and watch scripts now
shell out to the AWS CLI through `src/aws/cli.py` instead of boto3, because
the laptop's credentials come from `aws login` since 2026-10-22; see
`notes/aws_batch.md`, "Credentials".

## Open items (not done here)

- momentum24 runs (2026-10-22): c05 and c07 done, the other seven running
  on 2026-10-23. Pull, score and compare as set out in "First-principles
  review of the velocity gap"; then decide between the mass-only and the
  momentum source on the 24-cell layout for production. The production
  default in `geometry.yaml` is still `ps_momentum: false`.
- Wall-layer restraint (ranked cause 2): wall friction, rolling resistance
  and wall roughness are untested. Pre-register a criterion on D before any
  such run.
- Point-source momentum and wall velocity: asked the C26 authors for the
  .mfx/STL, the PS velocity, the apex geometry and probe definitions
  (2026-10-19); no reply yet. The 3 s mass-only A/B was scored on
  2026-10-20: inconclusive by the rule. The mass-only source fixes C1 but D
  doubles (see "c01–c03 mass-only 3 s results"). The MFiX 23.2 vs 26.1.2
  wall-contact review is done (see "Gas-fraction floor"); the apex geometry
  is still open.
- Wall velocity after the grid-layout A/B (2026-10-21): no effect at 12 and
  14 slpm and inconclusive at 10 slpm.
- CGP drag check (2026-10-22): 26.1.2 uses the real diameter. 23.2 very
  likely does too (indirect evidence; see "CGP drag check"). So the CGP
  implementation does not explain our 5.5–6.3× gap to C26's own model.
  Remaining candidates:
  - inputs C26 does not state: PS momentum, apex, floor/interpolation,
    contact settings, solids dt factor;
  - non-CGP version changes: 24.4/25.3 solids time step;
  - C26's wall-velocity definition.

  A direct 23.2 diff needs a refreshed MFiX download link.
- c03_massonly_axisface S6 (ΔP 3.09 W/A, over the 3.0 bound): cause not
  established. ΔP rose 0.18–0.32 W/A with the new layout at every flow.
- Mass-only spout lean towards −x: the gas-only tests are done (bias
  without a bed, not the domain split, not carried by the STL). The
  axis-face layout removes it (P3ii 0.20–0.36 mm), so the patched-solver test
  of the unlimited `Theta_Ue` in `get_alpha.f` is no longer needed for the
  production runs; see "Spout off-centre in the mass-only runs".
- Fig. 8e cone angle: measure it from a high-resolution extraction before
  changing geometry (the text says 60°).
- AWS per-step cost: the 4-rank AWS test is done (the gap is hardware or
  container, not rank count; 4 ranks adopted). Per-core speed and the
  container MPI settings (`--bind-to none`, btl) are untested.
- Synthetic-PIV assumptions that move CC ("Synthetic PIV"):
  - where the map bottom sits relative to the bed surface (a 3.4 mm shift
    moves c01's CC from 0.21 to 0.62);
  - whether empty windows count as zero in the time mean (c02_massonly CC
    0.13 against 0.55).

  Settle both independently of the score: ask the authors, or use Conry et
  al. 2025's PIVlab settings. C26 Fig. 6 wall v(t) traces are not
  digitized; only the Fig. 6d means are used.
- If the apex cut cell misbehaves, the fallback is to move the point source up
  one row (`ps_row_offset_cells: 1`). It should not be a widened inlet
  (recipe Rule 3 history).
- c07 and c09_massonly_axisface: still running on 2026-10-22. Score them
  with `notes/c26_watch_post.py <case> --tag ... --local` once pulled, and
  rerun `notes/c26_sweep_summary.py --suffix _massonly_axisface --tag
  axisface`.
- c06 bed expansion (both layouts, against the U_mf ordering): cause not
  established. See "Coating-stage / charge 3 s results".
- P2 cannot trace losses before the PIV window in 3 s runs (c05, c06
  NO-GO). Proposed: a sparse whole-run particle VTK. Not applied.

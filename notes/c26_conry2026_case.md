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
  cone" (§2.2.1). The authors chose dx ≈ nozzle diameter so that "when the
  flow mass is added to the domain, calculated velocities will be correct".
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

## Open items (not done here)

- Confirm the AWS per-step time of the dx/2 baseline against the ~0.2
  s/step estimate, from the `.TXT` once results sync at job end.
- Synthetic-PIV post-processor and digitized C26 Fig. 4 / wall-velocity data.
- If the apex cut cell misbehaves, the fallback is to move the point source up
  one row (`ps_row_offset_cells: 1`). It should not be a widened inlet
  (recipe Rule 3 history).

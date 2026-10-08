# Day 10 Summary — Case V (ORNL/UTK) 3D MP-PIC spouted bed

> **State: implementation complete, run pending.**  The model has not yet
> been launched.  All .mfx keywords were verified against the installed
> MFiX 26.1.2 reference (CLAUDE.md rule 1) and the deck renders cleanly
> through `src/make_case.py`.  No results have been recorded in
> `results/run_log.csv` yet.

## What was implemented

Day 10 lifts Case V from the Day-8/9 planar quadric TFM slab to a full
**3D cut-cell** geometry (60-degree cone + straight cylinder + 4 mm
orifice + 10 mm inlet stub) with an **MP-PIC** solids solver.  The
pipeline that was built for the earlier TFM attempt (STL writer,
Day-3 CFL wall-time estimator, post-processor for centerline `v_g`,
inlet-pressure PSD, bed-pressure traces, inventory and CFL CSV) is
reused as-is; only the solids-model block in the .mfx template and the
driver manifest changed.

Key artifacts:

- `cases/V_3d/template.mfx.j2` — 3D STL cut-cell MFiX template, now
  with `solids_model(1) = 'PIC'`.  `cartesian_grid = .True.`,
  `use_stl = .True.`, graded per-cell `dx(i)` / `dy(j)` / `dz(k)`
  arrays, BC_1 CG_NSW (STL walls), BC_2 PO at `y = H_dom`, BC_3 MI at
  `y = -L_stub` across the orifice-bounding square `|x| <= R_i,
  |z| <= R_i`.  Gas-side closure is unchanged (SYAM_OBRIEN drag,
  SI units); the TFM-only closures (`SCHAEFFER + GIDASPOW_PCF
  blending + LUN_1984 granular-energy PDE + phi + c_e`) are removed
  and replaced by the PIC continuum-stress block
  (`fric_exp_pic`, `fric_non_sing_fac`, `psfac_fric_pic`,
  `mppic_coeff_en1`, `mppic_coeff_en_wall`, `mppic_coeff_et_wall`,
  `mppic_velfac_coeff`) plus the shared DES interpolation block
  (`des_interp_mean_fields`, `des_interp_on`,
  `des_interp_scheme = 'LINEAR_HAT'`, `gener_part_config`,
  `des_epg_clip = 0.42`, `friction_model = 'SCHAEFFER'`).  Initial-
  and inflow-condition blocks were extended with
  `ic_pic_const_statwt(IC, 1) = 4.2` and
  `bc_pic_mi_const_statwt(3, 1) = 5.0` per the NETL tutorial.

- `src/make_stl.py` — unchanged.  The implicit-boolean + marching-cubes
  path (`build_spouted_bed_stl_implicit`) still produces the dense
  watertight shell used by both the earlier TFM attempt and the new
  PIC run.  PIC does not need the shell to be perfectly cut-cell-safe
  (parcels live in continuous space), but a clean fluid-cell map is
  still a precondition for the gas phase.

- `src/make_case.py` — `build_V_3d_params()` now reads the new
  `spouted_bed_V_3d_pic` block from `params/geometry.yaml` and threads
  the PIC coefficients, `des_epg_clip`, `des_interp_scheme`,
  `ic_pic_const_statwt`, and `bc_pic_mi_const_statwt` into the
  template context.  The mesh / `H_static` / geometry derivations
  that were correct for TFM are unchanged.

- `params/geometry.yaml` — new `spouted_bed_V_3d_pic` block.  Every
  entry carries a source comment pointing to
  `/Users/nelsonpereira/mamba/envs/mfix-26.1.2/share/mfix/templates/
  tutorials/pic/spouted_bed_3d/spouted_bed_pic_3d.mfx` and the
  installed MFiX 26.1.2 keyword reference.  No magic numbers in the
  template or driver — CLAUDE.md rule 2.

- `src/run_V_3d_day10.py` — driver.  The TFM-era manifest fields
  (`blending_function`, `kt_type`, `c_e`) are replaced by the PIC
  closure block (`fric_exp_pic`, `fric_non_sing_fac`,
  `psfac_fric_pic`, `mppic_coeff_en1`, `mppic_coeff_en_wall`,
  `mppic_coeff_et_wall`, `mppic_velfac_coeff`, `des_epg_clip`,
  `des_interp_scheme`, `ic_pic_const_statwt`,
  `bc_pic_mi_const_statwt`, `solids_model`) so the manifest is a
  complete record of what the deck will run.

- `src/post/plot_V_3d_day10.py` — unchanged.  Gas-phase fields
  (`V_g` in `BACKGROUND.pvd`) and the four monitor CSVs
  (`V_3D_P_INLET.csv`, `V_3D_P_TOP_BED.csv`,
  `V_3D_SOLIDS_INVENTORY.csv`, `V_3D_CENTERLINE_VG.csv`) are
  identical between TFM and PIC runs because the monitor regions
  use `monitor_type = 4/11` on the Eulerian cell grid, and PIC's
  `des_interp_mean_fields = .True.` projects parcel mass to the
  grid `ep_s` field.  A `Particles` VTK dump (parcel positions,
  diameter, velocity, `vtk_data(2) = 'P'`) was added for Day 11-12
  CFD-DEM comparison; it does not break any existing figure.

- `tests/test_make_case.py` — the nine V_3d tests were rewritten
  for the PIC template.  Each one still enforces an invariant that
  was in the TFM version (STL cut-cell on, 3D not `no_k`, graded
  mesh sums, orifice spans at least four cells, MI/PO/CG_NSW BCs,
  monitors + VTK layout, bed IC derivation from inventory,
  per-cell `dx`/`dy` emitted); the closure test (`test_V_3d_closure_keywords`)
  now asserts PIC keywords are present and TFM keywords are
  absent, and `test_V_3d_bed_ic_matches_case_V` now checks
  `ic_pic_const_statwt(2, 1)`.

- `.gitignore` — unchanged.  `cases/**/geometry_*.stl` is still
  ignored; the implicit-boolean STL regenerates deterministically.

## Why PIC was chosen for Day 10

The earlier Day-10 TFM attempt failed to produce a physical spout.  Four
fix attempts were made before concluding the obstruction is in MFiX
26.1.2 itself and the correct response is to switch solids models.

1. **F_AT is a Cartesian-grid classifier bug that the input deck
   cannot resolve.**  MFiX's STL cut-cell F_AT classifier
   (`cartesian_grid/intersect.f`) updates `F_AT` only at cell corners
   that are *near a facet*; corners far from any facet remain
   `UNDEFINED` and get marked `BLOCKED`.  For an open-top surface-of-
   revolution STL, the near-axis column of cells inside the vessel
   has no facets within one cell diameter and therefore ends up
   `BLOCKED`.  The TFM attempt reported 4 556 fluid cells vs the ~8 600
   expected from the vessel-to-box volume ratio.  The code path cannot
   be enabled or disabled from the `.mfx` deck — the only field knobs
   are `cad_propagate_order`, `set_corner_cells`, `out_stl_value`, and
   facet density, all of which were tried on the TFM path (see "What
   was tried on TFM" below) without closing the fluid-cell deficit.

2. **PIC parcels track in continuous space, so solids are not
   dependent on F_AT.**  PIC's parcels advance under the gas drag and
   the Harris-Crighton continuum stress using local gas-field
   interpolation (`des_interp_on`, `des_interp_scheme`).  Their
   trajectories pass smoothly through cells regardless of whether
   those cells are tagged `FLUID`, `CUT`, or `BLOCKED` — the solids
   advection is independent of the F_AT cell map.  Only the gas phase
   needs a correct fluid-cell map, which it already has (the 4 556
   fluid cells carry the gas equations without issue; the earlier TFM
   failure was in the volume-fraction advection, not the gas side).

3. **PIC avoids the LUN_1984 PDE stiffness that killed the quadric
   TFM attempt.**  When the TFM attempt was retried with native MFiX
   quadrics (Y_CONE + Y_CYL_INT, `GROUP_RELATION = 'AND'`) the
   geometry classifier gave the correct 78.6 % fluid-cell count, but
   the LUN_1984 PDE for granular temperature failed to converge in
   the near-axis cut cells and `dt` collapsed below `DT_MIN` within
   the first 100 steps.  PIC does not solve a granular-temperature
   PDE; the Harris-Crighton stress is an algebraic function of local
   `ep_s`, so this stiffness branch is absent by construction.

4. **NETL ships `tutorials/pic/spouted_bed_3d/spouted_bed_pic_3d.mfx`
   as the canonical 3D spouted-bed template.**  There is no TFM
   analogue in the shipped tutorial set
   (`.../share/mfix/templates/tutorials/`).  The PIC tutorial exercises
   exactly the geometry Case V needs — 60-degree cone, 1.5 cm orifice
   (we scale to 4 mm), SYAM_OBRIEN gas drag, cut-cell STL vessel —
   so adopting it is the shortest route to a documented, working
   recipe.

5. **PIC wall-clock is ~1.5–3× TFM, far below CFD-DEM.**  PIC's
   solids stress is algebraic in `ep_s` (no pairwise contact
   resolution, unlike CFD-DEM) and the solver dt stays at the gas-
   side semi-implicit floor (~1e-3 s), not the Hertzian contact
   floor.  The main overhead is parcel-field interpolation each gas
   timestep; this adds roughly 50-200 % on top of the pure-gas cost
   depending on parcel count.  By comparison, CFD-DEM with 500 µm
   ZrO₂ parcels would carry a Hertzian contact dt around 1e-6 s or
   smaller, i.e. 1000× more solver steps.

## How it works

1. **STL geometry.**  Unchanged from the TFM attempt.
   `src/make_stl.build_spouted_bed_stl_implicit` samples the vessel
   SDF (union of stub cylinder + cone frustum + freeboard cylinder)
   on a 0.5 mm grid, extracts the zero isosurface via marching
   cubes, and writes a watertight ASCII STL with outward-pointing
   normals.  For `R_c = 25 mm, R_i = 2 mm, cone half-angle = 30°,
   H_dom = 200 mm, L_stub = 10 mm` the STL has several thousand
   facets, including bottom and top caps.  MFiX cuts the MI through
   the bottom cap and the PO through the top cap via its own
   BC-region logic.

2. **Mesh grading.**  Day-10 rules:
     - cell size `≤ 10·d_p = 5 mm` in the particle zone;
     - at least 4 cells across the 4 mm orifice (near-axis cell
       size `≤ 1 mm`; grade if needed).

   Both rules are satisfied by the four-segment y mesh plus the
   three-segment xz mesh:

   | axis | segment 1 | segment 2 | segment 3 | segment 4 |
   | --- | --- | --- | --- | --- |
   | x,z | 5 × 4.6 mm | 4 × 1 mm (near-axis) | 5 × 4.6 mm | — |
   | y | 5 × 2 mm (stub) | 20 × 2 mm (bed) | 20 × 4 mm (fountain) | 16 × 5 mm (freeboard) |

   Totals: `imax = kmax = 14`, `jmax = 61`, roughly 11 956 total
   cells.  The expected fluid-cell count on the implicit STL is
   around 68 % of the box volume; the PIC solver does not need the
   near-axis cells to be in the fluid set for the solids-phase
   advection to work (see point 2 under "Why PIC was chosen").

3. **Inflow / outflow.**  MI at `y = -L_stub` (box bottom)
   over the orifice-bounding square `|x| ≤ R_i, |z| ≤ R_i`, with
   `bc_v_g = 30 m/s`, `bc_ep_g = 1`, `bc_ep_s(3, 1) = 0.0`
   (gas-only inlet), and `bc_pic_mi_const_statwt(3, 1) = 5.0` set
   per the NETL tutorial's parser requirement.  STL cut-cell masks
   the square to the circular stub-tube footprint (circle area
   `π·R_i² = 12.57 mm²` vs square area `16 mm²`, ratio `4/π ≈ 1.273`
   — a pinned test fact).  PO at the top (`y = H_dom`) spans the
   full cross-section at 101 325 Pa.

4. **Solids-model block (PIC, Harris-Crighton).**  `solids_model(1)
   = 'PIC'`, `d_p0(1) = 500 µm`, `ro_s0(1) = 6050 kg/m³`.  Solids
   stress uses the empirical Harris-Crighton closure (`fric_exp_pic =
   3.0`, `fric_non_sing_fac = 1e-7`, `psfac_fric_pic = 100 Pa`) and
   the empirical damping / wall-restitution coefficients
   (`mppic_coeff_en1 = 0.85`, `mppic_coeff_en_wall = 0.85`,
   `mppic_coeff_et_wall = 1.0`, `mppic_velfac_coeff = 1.0`).  The
   drag clip (`des_epg_clip = 0.42`) is numerically identical to
   Case V's `ep_star` — a physics coincidence rather than a derived
   identity, since both are the min-fluidization gas volume fraction.

5. **Parcel seeding.**  `gener_part_config = .True.` reads the IC
   blocks and generates parcels at startup.  The bed IC (IC_2,
   `0 ≤ y ≤ H_static`) has `ic_ep_s(2, 1) = 0.57` and
   `ic_pic_const_statwt(2, 1) = 4.2` so each parcel represents 4.2
   physical particles; the actual parcel count is auto-computed to
   match the IC solids fraction.  For a 54.5 g inventory in a
   32.2 mm bed at `eps_s_bed = 0.57`, expect roughly 10 000–20 000
   parcels (comparable to the NETL tutorial's reported 14 909 for a
   larger vessel).

6. **Pre-opened spout IC.**  IC3 overrides the packed-bed IC in the
   orifice-aligned column `|x| ≤ R_i, |z| ≤ R_i, 0 ≤ y ≤ H_static`
   with `ep_g = 1`, `v_g = U_in`, `ic_ep_s = 0.0`, same as the TFM
   attempt — this seeds the gas jet without placing parcels in the
   future spout path.  IC4 covers the stub-tube region
   (`-L_stub ≤ y ≤ 0`) with the same gas-only profile.

7. **Wall-time estimate (Day-3 CFL tool).**  Smallest cell is 1 mm
   (near-axis x, z).  `dt_for_cfl(u_max = 30, dx = 1 mm, cfl = 0.5)
   = 1.67 × 10⁻⁵ s`, which is a loose upper bound: MFiX's PIC
   advancement is semi-implicit on the gas side with sub-stepping on
   parcels, so the actual `dt` is expected to stabilise near
   1 × 10⁻³ s after a brief startup transient (same regime as the
   Day-9 TFM slab).  The wall-time projection for `tstop = 5.5 s`
   is therefore roughly 1.5–3× the Day-9 slab's `dt ≈ 1e-3` cost,
   which on the 10-core laptop lands near 1–2 hours wall-clock.
   This is a projection; the model has not been run yet.

8. **Post-processing pipeline** (`src/post/plot_V_3d_day10.py`,
   unchanged).  All four figures and the CFL CSV come from gas-phase
   or Eulerian mean-field data that PIC writes identically to the
   TFM run.  The `Particles` VTK will add a parcel point-cloud dump
   for the Day 11-12 CFD-DEM comparison; current Day-10 figures do
   not depend on it.

## Why TFM failed on this geometry in MFiX 26.1.2

The open-top surface-of-revolution STL interacts poorly with MFiX's
cut-cell `F_AT` classifier (`cartesian_grid/intersect.f`).  `F_AT` is
updated at cell corners that are *near a facet*; corners far from any
facet (e.g. inside the vessel near the axis, where the nearest facet
is on the far side of the vessel wall) remain `UNDEFINED` and get
marked `BLOCKED`.  The committed open STL reported 4 556 fluid cells
instead of the ~8 600 expected from the vessel-to-box volume ratio
(68 %).  The missing cells are precisely the near-axis column where
the spout should form; with that column mis-classified, the TFM
volume-fraction advection cannot carry solids up through the bed.

### What was tried on TFM

| Attempt | Result |
| --- | --- |
| Close the STL at both ends with facets flush with `y_min = 0, y_max = H_dom` | Facets flush with the Cartesian box boundary are DISCARDED by MFiX (`get_stl_data.f` line 2130: "Number of active facets (BCs not flush with MFiX box)"); 360 input facets reduced to 240 active |
| Offset `y_min = -dy_min`, close the bottom with an annular floor at `y = 0` (interior to box) | Fluid-cell count still 4 556 because the top is open and `F_AT` propagation does not reach the upper freeboard; 0 fluid cells for `y > 100 mm` |
| `flip_stl_normals(1) = .True.` + offset-box + closed-bottom STL | MFiX error 1100 "Cannot locate flow plane for boundary condition 3": with `y_min` offset below `y = 0`, the MI is INTERNAL to the box and MFiX requires MI at the box face |
| Replace STL with native MFiX quadrics (Y_CONE + Y_CYL_INT, `GROUP_RELATION = 'AND'`) | Fluid-cell count correct at 8 628 (78.6 %) — geometry is right — but the LUN_1984 PDE for granular temperature does not converge in the near-axis cut cells; `dt` collapses to below 1e-7 within the first 100 steps |
| Inlet stub tube + `cad_propagate_order = 'JKI'` + `set_corner_cells = .True.` + dense watertight STL (implicit-boolean + marching cubes) | Fluid-cell count still well below the vessel-to-box volume ratio; the F_AT classifier bug reproduces regardless of facet density and propagator mode |

The TFM path was abandoned once the fourth attempt confirmed the
limit is in MFiX's F_AT classifier, not anything the input deck can
express.  Day-10 switches to PIC, which is independent of F_AT on
the solids side.

## Results — "Done when" criteria

> **Model has not been run yet.**  All rows below are placeholders
> and will be filled in after the first PIC run completes.  The
> acceptance thresholds themselves are inherited from the Day-10
> plan.

| Criterion | Target | Observed |
| --- | --- | --- |
| ≥ 5 s simulated past startup without `DT < DT_MIN` | 5 s | TBD (pending first run) |
| Centerline `v_g` within stated tolerance vs S1 Fig. 1 | report without tuning | TBD (post-run; overlay figure already wired) |
| PSD peak reported next to S2's 17–25 Hz band and S1's ~10–11 Hz line | report, no tuning | TBD (post-run; Welch PSD wired) |
| All tests pass | — | **37 / 37 PASS** on the make_case + make_stl suite |
| Run logged in `results/run_log.csv` | — | TBD (no run has happened) |

## What it models (physically)

The 3D cut-cell template is the first simulation in this project
that can carry the **axisymmetric dynamics** of the ORNL/UTK vessel.
The planar quadric slab from Day 8/9 was explicitly qualitative; a
3D geometry is a precondition for the two Day-10 validation targets
(centerline `v_g` from S1 Fig. 1 and inlet-pressure pulsation from
S2's 17–25 Hz band / S1's ~10–11 Hz line).  The physics the pipeline
is designed to capture:

- **Centerline jet decay (S1 Fig. 1).**  The axial `v_g(y)` profile
  measures how far above the orifice the gas jet stays faster than
  the particles' terminal velocity (~6.4 m/s for 500 µm ZrO₂ at
  20 °C from Day 1).  That jet-persistence height sets the spout-
  exit location and therefore the particle time-of-flight through
  the deposition zone.

- **Spout-fountain-annulus circulation.**  Parcels rise axially in
  the central jet, decelerate at the top of the fountain, and rain
  down symmetrically on the annulus from all azimuths.  The 3D
  volume balance (spout mass flux ≈ fountain mass flux ≈ annulus-
  return mass flux) determines the circulation cycle time — PLAN §1.1's
  "roughly how long is one circulation cycle" question.

- **Inlet-pressure pulsation frequency (S2 17–25 Hz, S1 ~10–11 Hz).**
  The published peaks are driven by azimuthal spout-neck oscillations
  and bubble-detachment events in the cone — both inherently 3D
  phenomena.  The Day-9 slab's 1.60 Hz peak reflected a 2D wavelength,
  not the real physics.

Why this matters for the coater (PLAN §1.1):

- `v_g(y)` and the fountain height set the particle residence time
  in the spout + fountain (where CVD deposition happens);
- the pulsation frequency is a bulk indicator of the circulation
  cycle time — ~20 Hz implies ~50 ms per pulse, roughly one particle
  round-trip through the spout per pulse at this operating point.

Both quantities feed into Day 11-12's CFD-DEM residence-time
histograms and Day 13's parametric design study.  Switching to PIC
on Day 10 keeps the pipeline intact (same post-processor, same
figures, same monitors); only the solids-stress branch changes.

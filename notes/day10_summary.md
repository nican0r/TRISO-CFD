# Day 10 Summary — Case V (ORNL/UTK) 3D MP-PIC spouted bed

> **State: R1 implementation complete, re-run pending (2026-10-08).**
> First PIC attempt (commit 760d4b0) produced no spout — the gas jet
> collapsed to 4.3 m/s one cell above the orifice. Root cause traced
> to MFiX 26.1.2's F_AT cut-cell classifier mis-classifying near-axis
> cells on the graded fine-cell mesh. R1 re-implements the geometry
> as a clone of NETL's shipped `tutorials/pic/spouted_bed_3d` template
> (uniform 2 mm Cartesian mesh, full-face MI, stub-cylinder
> orifice-equivalent). Not yet re-run.
>
> **Follow-up (not Day-10 scope):** R1 was submitted via AWS Batch
> (`results/run_log.csv` rows `aws-batch,day10`); Day-11 CFD-DEM work
> continues from the R1 geometry — see `src/run_V_3d_day11.py` and
> `cases/V_3d/day11/`. Day-10 PIC was promoted to the baseline for
> Days 11-12 and 13 in commit `3346697`.

## What was implemented

Day 10 lifts Case V from the Day-8/9 planar quadric TFM slab to a full
**3D cut-cell** geometry (60° cone + cylinder + 4 mm orifice-equivalent
+ stub inlet) with an **MP-PIC** solids solver. The STL writer, CFL
wall-time estimator, and post-processor (centerline `v_g`, inlet-pressure
PSD, bed-pressure traces, inventory, CFL CSV) are reused from Day 8/9;
only the solids-model block in the template and the driver manifest
changed.

- `cases/V_3d/template.mfx.j2` — 3D STL cut-cell MFiX template with
  `solids_model(1) = 'PIC'`. `cartesian_grid = .True.`, `use_stl = .True.`,
  uniform 2 mm mesh, BC_1 CG_NSW (STL walls), BC_2 PO at `y = H_dom`,
  BC_3 MI across the full box-bottom face. `cad_propagate_order` and
  `set_corner_cells` are intentionally not emitted (NETL does not set
  them; see template header lines 59-61). Gas-side closure unchanged
  (SYAM_OBRIEN drag, SI units); TFM closures removed and replaced by
  the PIC continuum-stress block (`fric_exp_pic`, `fric_non_sing_fac`,
  `psfac_fric_pic`, `mppic_coeff_en1`, `mppic_coeff_en_wall`,
  `mppic_coeff_et_wall`, `mppic_velfac_coeff`) plus the DES interpolation
  block (`des_interp_mean_fields`, `des_interp_on`,
  `des_interp_scheme = 'LINEAR_HAT'`, `gener_part_config`,
  `des_epg_clip = 0.42`, `friction_model = 'SCHAEFFER'`).

- `src/make_stl.py` — `build_spouted_bed_stl_implicit` now takes a
  `top_extend_m` parameter (added for Day-10 PIC fix set; see run_log row
  `V_3d_day10_pic`) that pushes the top cap above `y = H_dom` so the PO
  plane cuts side-wall facets rather than the top cap, and the facet
  normals were corrected to point inward (into the fluid) in commit
  `88b67d9` (`fix(STL): correct normal direction for MFiX internal flow`).
  R1 driver calls it with `dx_sample = 1 mm` and `top_extend_m = 5 mm`,
  yielding 79 612 facets (see `cases/V_3d/day10/manifest.json`).

- `src/make_case.py` — `build_V_3d_params()` (line 316) reads the new
  `spouted_bed_V_3d_pic` block from `params/geometry.yaml` and threads
  PIC coefficients, `des_epg_clip`, `des_interp_scheme`,
  `ic_pic_const_statwt`, and `bc_pic_mi_const_statwt` into the template.
  The box in `x, z` extends to ±(R_c + 5 mm) = ±30 mm with uniform
  `dx = dy = dz = 2 mm` (`imax = kmax = 30`, `jmax = 110`, 99 000 cells);
  the `y` box spans `[-L_stub, H_dom] = [-20, 200] mm`.

- `params/geometry.yaml` — new `spouted_bed_V_3d` (R1 geometry) and
  `spouted_bed_V_3d_pic` (closure) blocks, every entry sourced to NETL's
  `tutorials/pic/spouted_bed_3d/spouted_bed_pic_3d.mfx` and the MFiX
  26.1.2 keyword reference (CLAUDE.md rule 2). A `spouted_bed_V_3d_dem`
  block has also landed for Day-11-12 (out of scope here).

- `src/run_V_3d_day10.py` — driver, manifest fields updated for PIC
  closure. The driver passes `top_extend_m = 5e-3` to the STL builder
  (`src/run_V_3d_day10.py:146`) and records `dx_sample_m` and
  `facet_count` in `manifest.json`.

- `src/post/plot_V_3d_day10.py` — monitors live on the Eulerian cell
  grid and PIC's `des_interp_mean_fields` projects parcel mass to
  `ep_s`, so all four figures (centerline `v_g`, ΔP, inventory, PSD)
  and the CFL CSV come through identically. The centerline post-processor
  uses a `manual_cell_centers` helper to avoid a PyVista `cell_centers()`
  hang on MFiX cut-cell VTUs (see `src/post/plot_V_3d_day10.py:37-55`).
  A `Particles` VTK dump (`vtk_data = 'P'`, `vtk_part_id = .True.`) was
  added to the template for the Day 11-12 CFD-DEM comparison and is
  consumed by `src/post/residence.py` (new file: per-parcel residence
  time inside caller-supplied spout/fountain region predicates).

- `tests/test_make_case.py` — V_3d test coverage rewritten for PIC and
  substantially expanded (now ~24 V_3d-only tests plus the DEM block for
  Day-11-12). The closure test (`test_V_3d_closure_keywords`,
  `tests/test_make_case.py:401`) asserts PIC keywords are present and
  TFM keywords absent; `test_V_3d_bed_ic_matches_case_V` checks
  `ic_pic_const_statwt(2, 1)`. New STL and residence suites also exist
  (`tests/test_make_stl.py`, `tests/test_residence.py`). Full
  `test_make_case + test_make_stl + test_residence` suite collects
  55 tests at time of writing.

## Why PIC + NETL-clone geometry

Two independent problems forced the change from TFM + custom geometry.

1. **TFM path.** MFiX 26.1.2's cut-cell F_AT classifier mis-classifies
   near-axis cells on open-top surface-of-revolution STLs (the near-axis
   column has no facets within one cell diameter and ends up `BLOCKED`).
   The TFM attempt reported 4 556 fluid cells instead of the ~8 600
   expected and the volume-fraction advection could not carry solids up
   through the spout. Four fix attempts (closing the STL, offset box,
   flipped normals, native MFiX quadrics) all failed — the quadric
   attempt gave the right fluid-cell count but LUN_1984 granular-
   temperature PDE stiffness collapsed `dt` within 100 steps. PIC
   parcels advance in continuous space (independent of F_AT on the
   solids side) and PIC's Harris–Crighton stress is algebraic (no PDE
   stiffness), so PIC removes both failure modes at once. NETL also
   ships a canonical 3D-PIC spouted-bed tutorial but no TFM analogue.

2. **First PIC attempt** (commit 760d4b0, state `V_3d_day10_pic`).
   Produced no spout: F_AT mis-classified the near-axis column for the
   **gas phase** too, so the inlet momentum never saw the fluid cells
   above the MI. The failure mode compounded three ways on the custom
   geometry: a graded fine-cell mesh (1 mm near-axis / 2 mm bed /
   4-5 mm freeboard), a dense 322 k-facet STL, and a tiny
   orifice-aligned square MI (16 mm² island at the box bottom with no
   F_AT seed for the near-axis cells above it). R1 adopts NETL's
   shipped pattern verbatim to avoid all three: uniform 2 mm Cartesian
   mesh, ~80 k-facet STL, straight stub cylinder as orifice-equivalent
   running from cone apex to box bottom, MI defined across the full
   box-bottom face (the STL blocks all cells outside the stub
   cylinder, so only the circular orifice-equivalent footprint sees
   inflow), and `cad_propagate_order` / `set_corner_cells` deleted
   (NETL does not set them).

Reference template:
`/Users/nelsonpereira/mamba/envs/mfix-26.1.2/share/mfix/templates/tutorials/pic/spouted_bed_3d/spouted_bed_pic_3d.mfx`.

## How it works

1. **STL geometry.** `build_spouted_bed_stl_implicit` samples the vessel
   SDF (union of stub cylinder + cone frustum + freeboard cylinder) on
   a 1 mm grid, extracts the zero isosurface via marching cubes, and
   writes a watertight ASCII STL with outward-pointing normals.

2. **Mesh.** Uniform `dx = dy = dz = 2 mm`. Box: `x, z = ±30 mm`
   (R_c + 5 mm margin); `y = [-L_stub, H_dom] = [-20 mm, 200 mm]`.
   `imax = kmax = 30`, `jmax = 110`, 99 000 total cells. Stub tube
   length `L_stub = 20 mm` (10× stub radius for inflow development).

3. **Inflow / outflow.** MI across the full `y = -L_stub` face with
   `bc_v_g = 30 m/s`, `bc_ep_g = 1`, `bc_ep_s(3, 1) = 0.0`, and
   `bc_pic_mi_const_statwt(3, 1) = 5.0`. The STL stub cylinder masks
   the MI down to the circular orifice-equivalent footprint
   (`π R_i² = 12.57 mm²`). PO at `y = H_dom` at 101 325 Pa.

4. **Solids-model block (PIC, Harris–Crighton).** `d_p0(1) = 500 µm`,
   `ro_s0(1) = 6050 kg/m³`. Solids stress uses `fric_exp_pic = 3.0`,
   `fric_non_sing_fac = 1e-7`, `psfac_fric_pic = 100 Pa`; damping /
   wall-restitution `mppic_coeff_en1 = 0.85`,
   `mppic_coeff_en_wall = 0.85`, `mppic_coeff_et_wall = 1.0`,
   `mppic_velfac_coeff = 1.0`. Drag clip `des_epg_clip = 0.42`.

5. **Parcel seeding.** `gener_part_config = .True.`. Bed IC (IC_2) has
   `ic_ep_s(2, 1) = 0.57`, `ic_pic_const_statwt(2, 1) = 4.2`; theoretical
   parcel count for the 54.5 g / 500 µm / 6050 kg/m³ ZrO₂ inventory is
   `N_phys / statwt = 1.376e5 / 4.2 ≈ 32 800` (see geometry.yaml comment
   on `ic_pic_const_statwt`, lines 206-215). The first PIC attempt that
   seeded IC on the whole domain reported 206 838 parcels (run_log row
   `V_3d_day10_pic`); R1 restricts seeding to the bed IC_2 region so the
   count is expected near the 32 k figure.

6. **IC layout.** The R1 template uses three ICs (`cases/V_3d/template.mfx.j2`
   lines 112-188): IC_1 covers the whole box with gas-only / `ic_ep_g = 1`
   (required by MFiX for multi-IC PIC runs; error 1003 fires without it);
   IC_2 is the packed bed (`0 ≤ y ≤ H_static`, `ic_ep_s = 0.57`); IC_3 is
   the inlet stub-tube region (`-L_stub ≤ y ≤ 0`, gas-only, `ic_v_g = 0`).
   The TFM-era pre-opened spout IC was **removed** (parcel-fluid drag
   spike and inconsistent with NETL's reference pattern). Monitors still
   sample the orifice-aligned column via `orifice_x_w_m` / `_e_m` /
   `orifice_z_b_m` / `_t_m` params so near-axis probes work despite the
   full-face MI.

7. **Wall-time estimate.** `dt_for_cfl(u = 30, dx = 2 mm, cfl = 0.5)
   ≈ 3.3×10⁻⁵ s` is a loose upper bound. MFiX's semi-implicit `dt` is
   expected to stabilise near 1×10⁻³ s (same regime as Day 9); for
   `tstop = 5.5 s` wall-clock projects to 1–2 hours on the 10-core
   laptop. Not yet run.

## Results — "Done when" criteria

> **R1 not yet run.** All simulation-dependent rows are placeholders.

| Criterion | Target | Observed |
| --- | --- | --- |
| ≥ 5 s simulated past startup without `DT < DT_MIN` | 5 s | TBD |
| Centerline `v_g` within stated tolerance vs S1 Fig. 1 | report without tuning | TBD (overlay figure wired) |
| PSD peak reported next to S2's 17–25 Hz band and S1's ~10–11 Hz line | report, no tuning | TBD (Welch PSD wired) |
| All tests pass | — | **55 collected** (make_case + make_stl + residence); last reported 37/37 at R1 commit, suite has since grown to include the Day-11 DEM block and residence tests |
| Run logged in `results/run_log.csv` | — | Logged (rows `V_3d_day10` TFM-era partial, `V_3d_day10_pic` first-PIC-attempt, plus two `aws-batch,day10` submissions) |

## What it models (physically)

The 3D cut-cell template is the first case in this project that can
carry the **axisymmetric dynamics** of the ORNL/UTK vessel. Day 8/9's
planar slab was qualitative; 3D is a precondition for the Day-10
validation targets:

- **Centerline jet decay (S1 Fig. 1).** The axial `v_g(y)` profile sets
  how far above the orifice the jet stays faster than the particles'
  terminal velocity, and therefore sets the spout-exit location and
  particle time-of-flight through the deposition zone.
- **Spout–fountain–annulus circulation.** 3D volume balance (spout flux
  ≈ fountain flux ≈ annulus return) determines the circulation cycle
  time — PLAN §1.1's "how long is one circulation cycle".
- **Inlet-pressure pulsation (S2 17–25 Hz, S1 ~10–11 Hz).** Driven by
  azimuthal spout-neck oscillations and bubble detachment in the cone,
  both inherently 3D. Day 9's 1.60 Hz slab peak was a 2D artifact.

Both quantities feed Day 11-12's CFD-DEM residence-time histograms and
Day 13's parametric sweep.

## Known open items

- R1 geometry has not been simulated locally; it was submitted to AWS
  Batch (see `results/run_log.csv` rows tagged `aws-batch,day10`); the
  "done when" validation rows remain TBD until the artifact is pulled
  back from `s3://triso-cfd-mfix/results/day10/`.
- Expected parcel count ≈ 32 800 (theoretical from inventory /
  `ic_pic_const_statwt`). Will be confirmed after the R1 run.
- The committed TFM-era and PIC-first-attempt rows in
  `results/run_log.csv` are historical records and are left in place.
- Figures and CSVs produced by the Day-10 pipeline sit under `results/`
  with the `fig_d10_V_3d_` prefix (`centerline_vg.{png,csv}`, `dP.{png,csv}`,
  `inventory.{png,csv}`, `psd.png`, `cfl.csv`). The current PNGs reflect
  the first-PIC-attempt (not R1).

## References

- **Snider (2001)**, *J. Comp. Phys.* 170: 523–549,
  DOI [10.1006/jcph.2001.6747](https://doi.org/10.1006/jcph.2001.6747) —
  MP-PIC method implemented as MFiX `solids_model = 'PIC'`.
- **Harris & Crighton (1994)**, *J. Fluid Mech.* 266: 243–276,
  DOI [10.1017/S0022112094000996](https://doi.org/10.1017/S0022112094000996) —
  frictional-stress closure (`psfac_fric_pic`, `fric_exp_pic`,
  `fric_non_sing_fac`).
- **NETL (2020)** `tutorials/pic/spouted_bed_3d/spouted_bed_pic_3d.mfx`
  shipped with MFiX 26.1.2 — reference template R1 clones.
- **S1** Collins et al., ORNL/TM-2006/520 — vessel geometry (50 mm
  column, 60° cone, 4 mm orifice), 54.5 g ZrO₂ inventory, centerline
  `v_g(y)` validation trace.
- **S2** Zhou et al., ANS Winter Meeting (2005) — inlet-pressure
  pulsation band (17–25 Hz) overlay reference.

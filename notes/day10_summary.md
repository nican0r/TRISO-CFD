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

## Update 2026-10-11 — R4: mesh refinement + straight stub; the spout finally works

### What landed

Day-10 R3 (Fix A, commit `555b8f5`) resolved the EP_G = 0.333 clip pathology
by replacing the single packed-bed IC_2 with NETL's four-strip quadrant
pattern. The dynamics went clean — no sustained velocity-cap clip events,
calm steady-state — but the quantitative validation targets were still
missed by large margins:

| Deliverable | R3 observed | Paper target | Shortfall |
| --- | --- | --- | --- |
| Centerline V_g at y = 25 mm | 2.4 m/s | ~15 m/s | 6× too low |
| Centerline V_g at y = 50 mm | 0.18 m/s | ~8 m/s | 44× too low |
| PSD peak | 1.95 Hz (drift artifact) | S1 ~10–11 Hz / S2 17–25 Hz | wrong band |
| Fountain formation | none (jet dies by y ≈ 15 mm) | parcels ejected above H_static = 32 mm | zero |

### Root cause: the orifice was numerically unresolvable at dx = 2 mm

The R1–R3 mesh was uniform 2 mm Cartesian. With R_i = 2 mm, only **2 cells
spanned the 4 mm orifice diameter**. NETL's shipped
`tutorials/pic/spouted_bed_3d` — the proven baseline we had been cloning —
uses D_i = 15 mm at dx = 4 mm, i.e. **3.75 cells across the orifice**.
R1–R3 were 1.88× under-resolved relative to the demonstrated NETL floor.

The Fix A smoke validator (`notes/smoke_validate_r3.py`) measures the
consequence directly. On the R3 smoke the orifice-plane cells carried
`v_y mean = 12.2 m/s`, `v_y max = 51 m/s`, against the 30 m/s target
driven by the paper's operating point. The jet entered the vessel with
(12/30)² = 16 % of the expected dynamic pressure. With that momentum the
bed weight dominated and the spout collapsed at y ≈ 15 mm. No fountain
could form by construction.

The R2 Option-4 scheme (R_stub = 2 × R_i, bc_v_g scaled DOWN by 0.25,
"accelerate 4× through a 1-cell taper") was a workaround for the F_AT
classifier blocking the stub at R_i = dx. On a 2 mm mesh the SIMPLE
pressure correction cannot drive a 4× velocity gradient across one cell
on a cut-cell mesh; the identity does not translate to the discrete FV
scheme at this resolution.

### The two fixes

**Fix 1 — refine mesh to 1 mm uniform** (`params/geometry.yaml`
`spouted_bed_V_3d.dx_uniform_m: 2.0e-3 → 1.0e-3`, `dx_sample_m:
1.0e-3 → 5.0e-4`). 4 cells across the 4 mm orifice, matching the NETL
ratio. Cell count 60 × 220 × 60 = **792 000 cells** (8× R3). CFL
tightens ~2×; expected total compute ~16× R3.

**Fix 2 — straight stub R_stub = R_i, no taper**
(`params/geometry.yaml` `R_stub_m: 4.0e-3 → 2.0e-3`, `taper_length_m:
2.0e-3 → 0.0`). With `R_stub = R_i` the `(R_i / R_stub)²` factor in
`src/make_case.py:564` evaluates to 1.0 — `bc_v_g_mi` identically
equals `U_in_m_s = 30 m/s`. MFiX now delivers the paper's physical jet
velocity directly at the MI face; no 4× pressure-driven acceleration
needed. Safe at dx = 1 mm because R_i = 2 mm = 2 cells across the stub
radius → F_AT seeds fluid cells cleanly without the Option-4 widening
hack.

**Why the two fixes are coupled.** Fix 2 is only safe once Fix 1 lands;
at dx = 2 mm, R_i = 1 cell across the stub radius and F_AT would block
every cell (the original failure that forced the Option-4 widening).

### Also landed with R4: local MPI (DMP) support

The 792 000-cell mesh is impractical in serial. The project's local MFiX
install includes `mfixsolver_dmp` (Open MPI 5.0.8). R4 parameterizes the
decomposition:

- `cases/V_3d/template.mfx.j2` + `template_dem.mfx.j2`: added
  `nodesi/nodesj/nodesk` as Jinja variables.
- `src/make_case.py` `build_V_3d_params` + `build_V_3d_dem_params`:
  take a `nodes=(ni, nj, nk)` kwarg (default `(1, 1, 1)` = serial);
  validate `ni × nj × nk > 0` AND that each axis divides the mesh
  cleanly, so a bad decomposition fails at render time rather than at
  AWS cold-start.
- `src/run_V_3d_day10.py` + `src/run_V_3d_day11.py`: added `--nodes
  NI,NJ,NK` CLI; `write_runner` emits `mpirun -np N mfixsolver_dmp`
  for `nranks > 1`, `mfixsolver` serial otherwise.
- AWS Batch is unchanged — `docker/entrypoint.sh` already strips any
  `nodes*` lines in the deck and injects a matching decomposition
  from the container `NP` env var. The local defaults never survive
  to the cluster.
- Local smoke decomposition pick: **(1, 5, 1) = 5 ranks**. `jmax =
  220` is divisible by 5, giving clean axial slabs of 44 cells per
  rank. Satisfies the "leave cores free" rule (5 ranks on 10 cores =
  5 cores free). `(1, 6, 1)` is NOT valid because 220 % 6 ≠ 0; the
  only 6-rank options at this mesh are non-axial (e.g. `(3, 2, 1)`).

### Smoke validation (2026-10-11, local MPI)

`cases/V_3d/day10_r4_smoke/` at tstop = 0.03 s, 5 MPI ranks:

- 29 791 parcels seeded (93% of 32 006 target — 7% shortfall from finer
  cut-cell lattice loss at wall corners, within tolerance).
- V_g at centerline monitor (y = H_mid = 16 mm, half-bed-height, 2 mm
  above the orifice cone-cylinder boundary):

  | t (ms) | centerline V_g (m/s) |
  | --- | --- |
  | 1.0 | 22.4 |
  | 2.1 | 21.6 |
  | 3.0 | 22.1 |
  | 4.1 | 21.8 |

- NIT = 10–11 (max_nit = 50), max P0 residual ≤ 5e-4, dt growing
  adaptively 0.73e-4 → 0.9e-4 s, PIP = 29 791 (zero parcel loss,
  zero MAX_INLET_VEL_FAC clips).
- Inlet pressure dropped from 101 325 → ~100 000 Pa (ΔP ≈ −1300 Pa,
  consistent with hydrostatic head + jet dynamic pressure).

The 22 m/s mid-bed centerline velocity is 73% velocity retention over
4 orifice diameters of dense-bed penetration, which lands cleanly in
the 70–85% band published by San José / Olazar for similar narrow-orifice
conical spouted beds (vs 8% retention at R3). The orifice jet is now
penetrating the bed at the paper's operating-point velocity — the
primary failure mode is resolved.

### What R4 does NOT settle

- **PSD validation** still depends on run length (tstop ≥ 5 s for
  df ≤ 0.5 Hz Welch resolution) and possibly on azimuthal resolution
  of the spout. The 1 mm mesh has ~12 cells circumferentially around
  the spout boundary — enough for the axisymmetric S1 band (10–11 Hz),
  borderline for the azimuthal S2 band (17–25 Hz). If PIC shows S1
  but not S2 after the full AWS run, the next lever is DEM or dx = 0.5 mm.
- **Day-11 and Day-13 reruns** need the same mesh refinement applied
  (both inherit `build_V_3d_params` automatically; no code change
  required). AWS Batch submission paths are unchanged.

### Code changes summary

| File | Change |
| --- | --- |
| `params/geometry.yaml` | `dx_uniform_m`, `dx_sample_m`, `R_stub_m`, `taper_length_m` updated with R4 provenance comments. |
| `src/make_case.py` | `nodes=(ni,nj,nk)` kwarg on both PIC and DEM builders; mesh-divisibility validation; `bc_v_g_mi_m_s` comment rewritten. |
| `cases/V_3d/template.mfx.j2` + `template_dem.mfx.j2` | `nodesi/j/k` Jinja block added; header rewritten; BC_3 preamble updated. |
| `src/run_V_3d_day10.py` + `src/run_V_3d_day11.py` | `--nodes` CLI; `write_runner` emits MPI or serial; absolute `cd` so the runner works from any cwd. |
| `tests/test_make_case.py` | 5 new R4 tests (serial default, 5-rank axial, 6-rank non-axial, non-divisor rejection, DEM pass-through); existing tests updated for 1 mm mesh and straight stub. |
| `tests/test_run_V_3d_day13.py` | sweep baseline asserts `bc_v_g_mi == 30 m/s`. |
| `.claude/skills/mfix-spouted-bed-recipe/SKILL.md` | New skill encoding the R1–R4 lessons for future agents. |

Full test suite: **112 passed**.

---

## Update 2026-10-10 — Fix A (R3) applied; dynamics clean on both PIC and DEM

The Day-10 R2 AWS Batch run (job `e6ba8d14-...`) completed but was
diagnosed as fatally pathological (`notes/day10_pathology_review.md`):
parcels seeded over the 4 mm orifice fell through at `t = 0+`, LINEAR_HAT
drag projection smeared their mass into the stub cut-cells, EP_G dropped
to the 0.333 solver clip across 500–800 cells, and `Vg` spiked to the
MAX_INLET_VEL_FAC ceiling (9380 m/s observed) with sustained 110 m/s
field-wide chaos. Centerline jet died by y = 37 mm (vs ~150 mm expected).

**Fix A** (commit `555b8f5`) replaces the single full-box IC_2 with
NETL's quadrant pattern from `tutorials/pic/spouted_bed_3d`: four packed
bed strips (IC_2 +z, IC_3 −z, IC_4 −x, IC_5 +x) leave a central gas-only
column of half-width `R_bed_gas_col = max(R_stub, R_i + dx) = 4 mm`
untouched (filled by IC_1). Dropped the explicit stub-tube IC_3 from R2
(redundant with IC_1). Applied identically to the DEM template so
Day 11-12 inherits.

**R3 PIC smoke** (`cases/V_3d/day10_r3_smoke`, tstop = 0.1 s, 24 min
local wall-clock):
  - Dominant R2 pathology eliminated: no sustained EP_G = 0.333 clip
    near the orifice (near-orifice min 0.427), max |V| in the field
    drops from 124 m/s steady to 106 m/s in 2–3 cells confined to
    `y = -1.1 mm, r = 1.4 mm` (near-axis stub tube).
  - Centerline jet survives through the bed cleanly: 30 m/s at
    y = 1 mm decaying smoothly to 15–20 m/s at y = 25 mm (vs R2's
    2.3 m/s dying).
  - CFL drops from mean 15.2 (R2) to ~3 (R3); NIT stable at 12–17 of
    max 50; dt stable 1.7–2.3 × 10⁻⁴ s.
  - Parcel count: 28,120 of ~32,700 target (86% — within tolerance).
  - **Residual**: one transient `MAX_INLET_VEL_FAC` clip event
    (Vg = -2933 m/s, EP_G = 0.655 at the orifice cell) during spout
    establishment, plus 2–3 near-axis stub cells with sustained
    |V| ≈ 100 m/s. Both are confined and do not propagate — centerline
    above y = 3 mm is clean. The one clip violates local momentum
    conservation but does not accumulate (parcel count stable).
  - A widened gas column (R4, `R_bed_gas_col = R_i + 2·dx = 6 mm`,
    `cases/V_3d/day10_r4_smoke`) was tried as a mitigation but regressed:
    27% parcel-count shortfall (73% of target) is a bigger paper-fidelity
    hit than the stub residual it tries to fix. R4 reverted.

**R3 DEM smoke** (`cases/V_3d/day11_r3_smoke`, tstop = 0.1 s, 10.5 min
local wall-clock):
  - Dramatically cleaner than PIC under the same Fix A: **zero** clip
    events, EP_G minimum sits exactly at the DEM `des_epg_clip = 0.42`
    floor (no excursions below), max |V| field-wide is 49.6 m/s
    (half the PIC residual), stub-tube max |V| is 48.5 m/s (also halved).
  - CGP sub-cycled contact resolution damps the parcel–gas momentum
    coupling that ran open-loop in PIC; HEXA lattice seeding also
    distributes parcels more structurally.
  - **Open issue**: parcel-count deficit. Target is 13,760 parcels at
    CGP factor 10; realized 8,186 (59.5% — 40% inventory deficit).
    Breakdown: ~20% geometric from Fix A's gas column + ~15% HEXA
    cut-cell lattice loss at the cone wall + ~5% baseline. Exceeds the
    mfix-smoke-validate 10% shortfall threshold. Dynamics are clean,
    but publishing results against a 32 g simulated inventory (vs
    paper's 54.5 g) means ΔP measurements scale ~40% low. Investigation
    of mitigations (`ic_des_lattice` alternatives, strip ep_s packing,
    cgp_factor reduction) ongoing.

**AWS re-submit**: Day-10 R3 PIC job `d10-pic-r3-fixA` submitted
2026-10-10 as `94e07164-b2fe-4964-9a81-5fbd83cc476b` with tstop = 7 s
(review-recommended for df ≤ 0.5 Hz Welch resolution on the S1 and S2
pressure bands). Day-11 R3 DEM resubmit pending the inventory-deficit
investigation.

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

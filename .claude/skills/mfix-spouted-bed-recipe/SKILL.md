---
name: mfix-spouted-bed-recipe
description: |
  Build or modify an MFiX 26.1.2 spouted-bed deck for a conical bed with a narrow orifice (D_c/D_i ~> 10, ORNL/UTK coater class) on the first try. Use whenever the task involves building/modifying an MFiX spouted-bed deck, a cone + narrow-orifice vessel, a conical spouted-bed CFD-DEM or MP-PIC model, or diagnosing "spout fails to form / dies / no fountain", "jet decays too fast", "centerline velocity collapses above orifice", "EP_G = 0.333 clip", "MAX_INLET_VEL_FAC clip", "F_AT blocked cells near axis", "parcel count seeded over orifice". Encodes 4 rounds (R1-R4) of hard-won Day-10 iteration lessons so the next agent does NOT retry proven-failed approaches (graded meshes, 2 mm uniform meshes, R_stub > R_i widening with mass-flux scaling, single-rectangle packed IC, TFM on cut-cell geometries, bare orifice MI without the stub pipe, etc.).
---

# MFiX spouted-bed recipe

Build a working MFiX spouted-bed deck on the **first** try.

This skill encodes the recipe the project converged on after Day-10 R1/R2/R3 PIC
iterations, Day-11 DEM/CGP iterations, Day-13 D-sweep contamination, and the
Day-10 R4 mesh+stub fix. Each rule below has a concrete MFiX keyword, a
"why" (which failed attempt proved it), and a source citation. Follow the rules
as written; do not re-derive them from first principles.

The target geometry class is **ORNL/UTK coater**: cylindrical column with a
60-degree cone transitioning to a narrow 4 mm orifice (D_c/D_i = 12.5), 500 um
ZrO2 particles at 6050 kg/m^3, 54.5 g inventory, U_in ~ 30 m/s at the orifice
plane (U/U_ms ~ 1.9). The recipe generalises to any conical narrow-orifice
spouted bed in MFiX 26.1.2.

---

## Decision rule

| Situation | Action |
|---|---|
| Building new conical narrow-orifice spouted-bed deck in MFiX 26.1.2 | **Apply this skill verbatim.** Deviate from a rule only with a Day-10-pathology-review-grade justification. |
| Modifying an existing project deck (template.mfx.j2 / template_dem.mfx.j2) | Diff against this recipe; any deviation needs a note in `notes/`. |
| User says "spout dies", "no fountain", "jet decays too fast", "ratio v_orifice / U_in << 1" | Jump to the diagnostic matrix (Section 7) first; most failures are Rule 1 (mesh) or Rule 6 (IC). |
| Non-conical geometry (fluidised bed, hopper, rotating drum) | **Do not use this skill.** The quadrant IC, stub pipe, and 4-cells-across-orifice rule are spouted-bed specific. |
| MFiX version != 26.1.2 | Re-verify every keyword against the installed version's `init_namelist.f` / keyword HTML. Keywords named here (BLENDING_FUNCTION, GIDASPOW_PCF, etc.) have changed between versions. |

---

## The 20 hard-won rules

Each rule is written as a concrete directive. "Why" cites the specific failed
iteration that proved it. "Source" cites the exact MFiX file / NETL tutorial /
project note the rule comes from.

### Rule 1 — MESH: `D_i / dx >= 4` cells across the orifice DIAMETER

- **Rule:** Set `dx_uniform_m` so that `(2 * R_i) / dx >= 4`. For `R_i = 2 mm`
  (D_i = 4 mm), use `dx = 1 mm`. If the orifice is 15 mm (NETL tutorial),
  dx = 4 mm is enough.
- **Why:** NETL's shipped `tutorials/pic/spouted_bed_3d` resolves its 15 mm
  orifice with 3.75 cells at dx = 4 mm; that is the demonstrated floor.
  Day 10-13 burned ~16 h of AWS Batch time running `dx = 2 mm` with
  `R_i = 2 mm` (2 cells across D_i) and the smoke validator measured orifice
  `v_y mean = 12.2 m/s` against a 30 m/s target. The jet carried
  (12/30)^2 = 16% of the expected dynamic pressure, could not penetrate the
  bed, and died at y ~ 15 mm with no fountain forming.
- **Finer tradeoff:** `dx = 0.5 mm` (8 cells across D_i, 6.3 M cells) is only
  warranted if the 1 mm smoke still misses the quantitative targets after
  Rules 2-20 are applied. 1 mm is the committed Day-10 R4 baseline (792 000
  cells, ~16x the compute of the 2 mm run but still tractable on AWS Batch).
- **Source:** NETL `tutorials/pic/spouted_bed_3d/spouted_bed_pic_3d.mfx`;
  `notes/day10_summary.md R4`.

### Rule 2 — MESH: UNIFORM Cartesian, not graded

- **Rule:** `dx(i) = dy(j) = dz(k) = dx_uniform_m` for all i, j, k. Do NOT grade
  near-axis vs outer cells.
- **Why:** Day-10 first attempt used a graded mesh (1 mm near-axis / 2 mm bed /
  4.6 mm outer) with a dense STL. MFiX 26.1.2's F_AT cut-cell classifier
  misclassified the near-axis fine cells as BLOCKED. The spout never formed
  because the inlet could not reach any fluid cell path up the axis.
- **Keyword stance:** `cad_propagate_order` and `set_corner_cells` are
  intentionally NOT emitted. NETL's tutorial does not set them either; on
  graded meshes with dense STLs they interact pathologically.
- **Source:** `notes/day10_summary.md` ("R1 - NETL-clone redesign") and
  `cases/V_3d/template.mfx.j2` lines 63-65.

### Rule 3 — STUB: `R_stub = R_i`, NOT widened + mass-flux-scaled

- **Rule:** Set `R_stub_m = R_i_m` (straight stub pipe, same radius as the
  orifice). Set `taper_length_m = 0`. Set `bc_v_g = U_in_m_s` directly at the
  MI, no `(R_i/R_stub)^2` scaling.
- **Why:** Day-10 R2 ("Option 4") set `R_stub = 2 * R_i` and scaled
  `bc_v_g` DOWN by 0.25 so the orifice velocity would recover U_in via 4x
  acceleration across a 1-cell taper. On a 2 mm cut-cell mesh the SIMPLE
  pressure correction cannot drive that gradient; the orifice-plane velocity
  came in at 12.2 m/s not 30 m/s. The widening was a workaround for the
  F_AT-blocks-the-stub problem that only existed on the 2 mm mesh. On the
  1 mm mesh (Rule 1), F_AT seeds the stub as fluid directly and the
  workaround becomes counter-productive.
- **Keyword:** `R_stub_m` (project yaml param) feeds the STL builder; the
  MFiX template keeps `(R_i/R_stub)^2` in the `bc_v_g_mi_m_s` expression as
  a no-op (= 1.0) for backward compatibility.
- **Source:** `notes/day10_summary.md` R4 sections "Why R4 exists" and
  "Fix 2 - straight stub".

### Rule 4 — STUB: provide a straight inlet pipe of length `>= 10 * R_i`

- **Rule:** The stub cylinder extends from `y = -L_stub` to `y = 0` (the cone
  apex). Set `L_stub_m >= 10 * R_i_m`. The project uses `L_stub = 20 mm = 10 * R_i`.
  Do NOT place the MI directly at the cone apex (y = 0).
- **Why:** The MI delivers a flat velocity profile. Without a straight pipe
  below the cone, the profile hits the diverging cone wall immediately and
  the near-wall boundary layer never develops, producing an off-axis velocity
  spike. 10 radii of straight pipe lets the profile relax before the
  divergence.
- **Source:** NETL `tutorials/pic/spouted_bed_3d` (uses L_stub ~ 10 x R_i);
  `params/geometry.yaml` `spouted_bed_V_3d.L_stub_m`.

### Rule 5 — BC: MI spans the FULL box bottom face

- **Rule:**
  ```
  bc_type(3) = 'MI'
  bc_x_w(3)  = x_min_m    (= full x extent)
  bc_x_e(3)  = x_max_m
  bc_y_s(3)  = y_min_m    (= -L_stub)
  bc_y_n(3)  = y_min_m
  bc_z_b(3)  = z_min_m    (= full z extent)
  bc_z_t(3)  = z_max_m
  ```
  The STL stub cylinder masks every cell outside the circular stub footprint
  as BLOCKED, so only `pi R_stub^2` actually sees inflow. Do NOT define a
  tiny orifice-aligned square MI patch.
- **Why:** The project's pre-Day-10 attempt used a 4 mm x 4 mm MI patch
  (16 mm^2 island at the box bottom with no F_AT seed for the near-axis cells
  above it). F_AT could not classify the near-axis column and blocked it.
  NETL's two reference tutorials (PIC spouted_bed_3d, DEM spouted_bed_pseudo_2d)
  both use a full-face MI with STL-defined footprint.
- **Source:** `cases/V_3d/template.mfx.j2` lines 266-285; NETL PIC and DEM
  tutorials.

### Rule 6 — IC: 4-quadrant packed-bed pattern with central gas-only column

- **Rule:** Split the packed-bed IC into four strips that leave a central
  gas-only column of half-width `R_bed_gas_col = max(R_stub, R_i + dx)`
  untouched (filled by the whole-domain gas IC).
  ```
  IC 1: whole domain, ep_g = 1.0 (gas-only background; see Rule 7)
  IC 2: +z strip   (full x,  z in [+R_bed_gas_col, z_max])
  IC 3: -z strip   (full x,  z in [z_min, -R_bed_gas_col])
  IC 4: -x center  (x in [x_min, -R_bed_gas_col], z in [-R_bed_gas_col, +R_bed_gas_col])
  IC 5: +x center  (x in [+R_bed_gas_col, x_max], z in [-R_bed_gas_col, +R_bed_gas_col])
  ```
  All four strips carry `ic_ep_s(k, 1) = eps_s_bed` (= 0.57 for Case V) and
  the same y-range `[0, H_static_m]`.
- **Why:** The pre-Day-10-R3 single-rectangle packed IC seeded parcels
  directly over the 4 mm orifice. At t = 0+, those parcels fell through the
  open orifice; LINEAR_HAT drag projection (Rule 13) smeared their mass with
  ~2-cell support into the stub cells. EP_G dropped to the 0.333 solver clip
  across 500-800 cells, mass balance forced gas velocity from 13 m/s (plug
  flow) up to the MAX_INLET_VEL_FAC ceiling, and 9380 m/s Vg values were
  observed with 110 m/s field-wide chaos. Fix A (Day-10 R3, commit `555b8f5`)
  implements the quadrant pattern verbatim from NETL.
- **Source:** NETL `tutorials/pic/spouted_bed_3d/spouted_bed_pic_3d.mfx`
  IC_2-5; `notes/day10_pathology_review.md` section 3 "Fix A".

### Rule 7 — IC: `ic_ep_g(1) = 1.0` on the whole-domain background IC is REQUIRED

- **Rule:** IC 1 MUST span the entire box (x_min..x_max, y_min..y_max,
  z_min..z_max) with `ic_ep_g(1) = 1.0`, `ic_ep_s(1, 1) = 0.0`,
  `ic_u_g(1) = ic_v_g(1) = ic_w_g(1) = 0`, `ic_p_g(1) = 101325`,
  `ic_t_g(1) = 293.15`. The sub-region ICs (bed quadrants) override within
  their boxes.
- **Why:** MFiX fires error 1003 in
  `check_data/check_ic_common_discrete.f:156` on multi-IC PIC runs if IC_1
  does not span the entire domain with `ep_g = 1`. This is unobvious: the
  parser passes without it, but the solver aborts at the first time step.
- **Source:** `cases/V_3d/template.mfx.j2` lines 115-138; MFiX 26.1.2
  `check_data/check_ic_common_discrete.f`.

### Rule 8 — IC: `ic_des_abort_on_failed_seeding(k, 1) = .False.` per bed IC

- **Rule:** Set `ic_des_abort_on_failed_seeding(k, 1) = .False.` on every bed
  quadrant IC (k = 2..5). For DEM/CGP also set `ic_des_lattice(k, 1) = 'HEXA'`.
- **Why:** Cut-cell wall corners drop some lattice sites. Up to ~5-10%
  parcel-count shortfall is tolerable (geometric loss, not a bug).
  **> 10% shortfall is a signal the IC geometry is wrong** — investigate
  Rule 6's `R_bed_gas_col`, the strip bounds, or the cone-wall cut-cell
  pattern before raising the ceiling.
- **Observed project values:** Day-10 R3 PIC seeded 28 120 of 32 700
  parcels (86%, acceptable); Day-11 R3c DEM/CGP seeded 18 527 of ~27 000
  parcels (67%, documented as inventory deficit — pay-the-price call).
- **Source:** `cases/V_3d/template_dem.mfx.j2` lines 148-220; NETL DEM tutorial.

### Rule 9 — TIMESTEP: `dt_max <= 0.5 * dx / U_in`

- **Rule:** For `dx = 1 mm`, `U_in = 30 m/s`, use `dt_max <= 1.7e-5 s`
  (CFL <= 0.5). `dt_min = 1.0e-7 s`, `dt_fac = 0.9`, `detect_stall = .True.`
  are standard.
- **Why:** The pre-Fix-B deck had `dt_max = 1e-3 s` at `dx = 2 mm`,
  `U_in = 30 m/s` running sustained CFL = 15 (max 38). Any mass-flux imbalance
  was amplified ~7x per step by the semi-implicit SIMPLE loop before the
  Picard iterations could correct it. On a correct deck CFL < 1 keeps NIT
  within its budget; CFL > 2 means the solver is relying on implicit mode
  without a safety signal.
- **Source:** `notes/day10_pathology_review.md` section 3 "Fix B";
  `cases/V_3d/template.mfx.j2` line 48 `dt_max`.

### Rule 10 — PIC CLOSURE: Harris-Crighton + SCHAEFFER with NETL defaults

- **Rule (verbatim):**
  ```
  friction_model      = 'SCHAEFFER'
  fric_exp_pic        = 3.0
  fric_non_sing_fac   = 1.0e-7
  psfac_fric_pic      = 100.0
  mppic_coeff_en1     = 0.85
  mppic_coeff_en_wall = 0.85
  mppic_coeff_et_wall = 1.0
  mppic_velfac_coeff  = 1.0
  des_epg_clip        = 0.42
  ```
- **Why:** All values from NETL `tutorials/pic/spouted_bed_3d/spouted_bed_pic_3d.mfx`.
  These are the proven-stable defaults for a 60-deg cone + narrow orifice +
  SYAM_OBRIEN drag combination. Do not tune individually — if PIC is
  unstable, the problem is upstream (Rule 1, 6, or 9), not the stress
  closure.
- **Source:**
  `/Users/nelsonpereira/mamba/envs/mfix-26.1.2/share/mfix/templates/tutorials/pic/spouted_bed_3d/spouted_bed_pic_3d.mfx`.

### Rule 11 — DEM/CGP CLOSURE: Cundall-Strack soft-sphere, NETL defaults

- **Rule (verbatim):**
  ```
  kn                   = 1000.0        # N/m, particle-particle
  kn_w                 = 1000.0        # N/m, particle-wall
  kt_fac               = 0.2857142857  # = 2/7 (Cundall-Strack)
  kt_w_fac             = 0.2857142857
  mew                  = 0.3           # p-p friction (ceramic/ceramic)
  mew_w                = 0.3           # p-wall friction
  des_en_input(1)      = 0.9           # p-p normal restitution
  des_en_wall_input(1) = 0.9           # p-wall normal restitution
  des_etat_fac         = 0.5
  des_etat_w_fac       = 0.5
  cgp_stat_wt(1)       = 10            # 1 for pure DEM
  ```
- **Why:** Base values from NETL `tutorials/dem/spouted_bed_pseudo_2d/
  spouted_bed_dem_pseudo_2d.mfx`. Project uses `mew = 0.3` (ceramic-on-
  ceramic) vs NETL's 0.1 (glass beads); project uses `des_en_input = 0.9`
  (matches `case_V.c_e` from `params/particles.yaml`) vs NETL's lower value.
  CGP factor 10 cuts the real-particle count 1.376e5 to ~14 k parcels and
  scales Rayleigh dt_c to ~2e-5 s so a 5+ s run fits in the AWS Batch window.
- **Rayleigh dt:** `dt_c = pi * sqrt(m_parcel / kn)`. Target DEM substep
  `dem_dt_s ~ 0.1 * dt_c` (ends up near Rayleigh at CGP = 10).
- **Source:**
  `/Users/nelsonpereira/mamba/envs/mfix-26.1.2/share/mfix/templates/tutorials/dem/spouted_bed_pseudo_2d/spouted_bed_dem_pseudo_2d.mfx`;
  NETL `tutorials/cgp/fluidbed_cgp_3d` for `cgp_stat_wt`.

### Rule 12 — DRAG: `drag_type = 'SYAM_OBRIEN'` as the isothermal baseline

- **Rule:** `drag_type = 'SYAM_OBRIEN'` with `drag_c1 = 0.8` and
  `drag_d1 = 2.65`. Use for Day-10 PIC baseline, Day-11-12 DEM/CGP, and the
  reference arm of any Day-13 drag-axis sweep.
- **When to switch:** `GIDASPOW_BLEND` is appropriate for non-isothermal
  runs (S1 Table 1 recommendation). The Day-13 D-sweep uses both drag laws
  as a sensitivity axis; do NOT switch for the baseline run.
- **Why:** SYAM_OBRIEN is the NETL default across both spouted-bed
  tutorials. It is well-conditioned at the ep_g ~ 0.42-0.58 range where
  the bed sits between IC and spouting.
- **Source:** NETL PIC and DEM tutorials; `cases/V_3d/template.mfx.j2` line 52.

### Rule 13 — INTERP: `LINEAR_HAT` with mean-fields on

- **Rule:**
  ```
  des_interp_mean_fields = .True.
  des_interp_on          = .True.
  des_interp_scheme      = 'LINEAR_HAT'
  gener_part_config      = .True.
  ```
- **Why:** NETL default (both PIC and DEM tutorials). **Known risk:**
  LINEAR_HAT smears parcel mass with ~2-cell support. This is what caused
  Rule 6's EP_G = 0.333 clip when parcels seeded over the orifice. The
  quadrant IC (Rule 6) is specifically what makes LINEAR_HAT safe for this
  geometry class — do NOT revert the quadrant IC in favour of LINEAR_HAT's
  cell-centred alternative (`des_interp_scheme = 'NONE'`); NETL's PIC
  tutorial uses LINEAR_HAT and the quadrant IC together as a matched pair.
- **Source:** NETL PIC tutorial;
  `notes/day10_pathology_review.md` section 2 "Root cause".

### Rule 14 — PARAMS live in `params/*.yaml`, with units + source comment

- **Rule:** Every physical parameter goes into `params/geometry.yaml` or
  `params/particles.yaml` with a unit comment and a `source:` line (paper +
  page, or `assumed`). Python scripts (`src/make_case.py`, `src/run_*.py`)
  must read these values, never hard-code them.
- **Why:** CLAUDE.md rule 2. The project has already added
  `spouted_bed_V_3d`, `spouted_bed_V_3d_pic`, `spouted_bed_V_3d_dem`, and
  `operating_point_V` blocks. New cases/sweeps extend these blocks (e.g.
  `params/sweep_V_3d.yaml` for the D-sweep), they do NOT redefine the
  values inline in Python.
- **Source:** `/Users/nelsonpereira/Documents/GitHub/Personal_Projects/Nuclear/TRISO-CFD/.claude/CLAUDE.md`
  rule 2; `params/geometry.yaml`.

### Rule 15 — KEYWORDS: verify every MFiX keyword against installed reference

- **Rule:** Before using any keyword, check:
  ```
  /Users/nelsonpereira/mamba/envs/mfix-26.1.2/share/mfix/doc/html/reference/
  /Users/nelsonpereira/mamba/envs/mfix-26.1.2/share/mfix/src/model/init_namelist.f
  /Users/nelsonpereira/mamba/envs/mfix-26.1.2/share/mfix/src/cartesian_grid/cartesian_grid_init_namelist.f
  /Users/nelsonpereira/mamba/envs/mfix-26.1.2/share/mfix/src/check_data/check_*.f
  ```
  If a keyword is not in `init_namelist.f`, the parser rejects it even if
  the HTML header mentions it.
- **Known gotchas:**
  - `BLENDING_FUNCTION = 'SIGM_BLEND'` is documented in the HTML namelist
    header but **rejected by the parser** in 26.1.2. Use the drag model
    option `GIDASPOW_PCF` instead (it embeds the blend).
  - `cad_propagate_order` and `set_corner_cells` are valid keywords but
    SHOULD NOT be set (Rule 2); NETL omits them.
  - `BC_PIC_MI_CONST_STATWT(BC, PHASE)` must be set on a PIC MI even when
    `bc_ep_s(BC, 1) = 0` (no parcels enter via that BC) — the parser
    requires the field. See `cases/V_3d/template.mfx.j2` line 280.
  - `IC_PIC_CONST_STATWT(IC, PHASE)` must be set on EVERY IC that could
    hold parcels, including the whole-domain gas IC_1 (set to a sentinel
    value; no parcels seed there because `ic_ep_s(1, 1) = 0`).
- **Source:** CLAUDE.md rule 1; `cases/V_3d/template.mfx.j2` header lines 35-38.

### Rule 16 — ALWAYS smoke-run locally before cloud submission

- **Rule:** Run `tstop = 0.1 - 0.3 s` locally. Validate per the
  `mfix-smoke-validate` skill. Only then submit to AWS Batch / EC2 / HPC.
- **Why:** The project lost ~$9 and 11 h of AWS Batch wall time to decks
  that passed "deck validates" but had zero fluid cells in the stub (R_i /
  dx = 1 cell case). A 10-minute local smoke run with the orifice-velocity
  check and fluid-cell-count check catches this before cloud submission.
  CLAUDE.md rule 13 "Dry-run locally before any cloud submission."
- **Source:** `.claude/skills/mfix-smoke-validate/SKILL.md`;
  `notes/day10_pathology_review.md`.

### Rule 17 — MP-PIC vs CFD-DEM vs TFM

- **Rule:**
  - **MP-PIC** (`solids_model(1) = 'PIC'`): use as the baseline on
    cut-cell spouted-bed geometries in MFiX 26.1.2.
  - **CFD-DEM / CGP** (`solids_model(1) = 'DEM'` or `'CGP'`): use when
    particle-level residence histograms are the deliverable (e.g. coating
    uniformity). Dramatically cleaner dynamics than PIC under the same Fix A
    (0 clip events vs 1-3 transient events) at the cost of ~2x wall time.
  - **TFM**: **do NOT use** on cut-cell spouted-bed geometries. The
    LUN_1984 granular-temperature PDE collapses dt within 100 steps on
    quadric or STL spouted beds (observed Day-8/9 iteration). PIC's
    Harris-Crighton stress is algebraic and removes the PDE stiffness.
- **Why:** Four Day-8/9 TFM fix attempts (closed STL, offset box, flipped
  normals, native quadrics) all failed on the TFM path. PIC succeeded on
  the first NETL-cloned attempt (after Rules 1-6 landed).
- **Source:** `notes/day10_summary.md` section "Why PIC + NETL-clone
  geometry"; `notes/day11_summary.md` section "Clean dynamics across all
  three smokes" for the DEM-vs-PIC quantitative comparison.

### Rule 18 — NETL's PIC tutorial geometry is NOT the UTK/ORNL geometry

- **Rule:** CLONE NETL's **IC pattern, BC pattern, closure parameters, and
  mesh convention** (Rules 1-15 above). DO NOT clone NETL's geometry.
- **Why:** NETL's `tutorials/pic/spouted_bed_3d` has a 15 mm orifice-equivalent
  cylinder that runs the full height of the vessel as a central "spine" —
  gas never sees a narrowing orifice. The UTK/ORNL coater has a 4 mm orifice
  at the cone apex with no spine. The project burned ~1 week of iterations
  assuming NETL's deck transferred verbatim; it does not.
- **What transfers:** mesh resolution ratio (3.75-4+ cells across orifice),
  full-face MI pattern, 4-quadrant IC, Harris-Crighton closure, SCHAEFFER
  friction, SYAM_OBRIEN drag, LINEAR_HAT interpolation.
- **What does NOT transfer:** absolute orifice size, cone angle, vessel
  diameter, inventory, parcel count, `U_in` value, bed height.
- **Source:** `notes/day10_summary.md` section "Why PIC + NETL-clone
  geometry" paragraph on NETL vs UTK geometry.

### Rule 19 — PO at `y = H_dom`, full x,z extent, atmospheric

- **Rule:**
  ```
  bc_type(2) = 'PO'
  bc_x_w(2)  = x_min_m, bc_x_e(2) = x_max_m
  bc_y_s(2)  = H_dom_m, bc_y_n(2) = H_dom_m
  bc_z_b(2)  = z_min_m, bc_z_t(2) = z_max_m
  bc_p_g(2)  = 101325
  ```
  **Separate requirement:** the STL top cap must sit ABOVE `y = H_dom` by at
  least `top_extend_m ~ 5 mm`. The STL builder (`src/make_stl.py
  build_spouted_bed_stl_implicit`) takes this as an argument.
- **Why:** If the STL top cap is at exactly `y = H_dom`, the PO plane cuts
  a horizontal facet instead of the side-wall facets. MFiX's F_AT ray-cast
  then classifies the top cells as BLOCKED, blocking the outflow entirely.
  Push the STL cap 5 mm above the PO plane so the PO cuts only side-wall
  facets.
- **Source:** `notes/day10_summary.md` section "What was implemented" bullet
  on `top_extend_m`; `src/run_V_3d_day10.py:146`.

### Rule 20 — STL NORMALS point INWARD (into the fluid)

- **Rule:** STL facet normals point into the fluid interior, not outward.
  The project's `build_spouted_bed_stl_implicit` enforces inward normals
  (fixed in commit `88b67d9`); agents writing a custom STL path must do
  the same.
- **Why:** MFiX 26.1.2's F_AT classifier uses normal direction to decide
  "fluid" vs "solid". Outward normals invert the classification: the
  cone wall becomes fluid and the vessel interior becomes blocked. The
  project's first STL attempt had outward normals and the whole bed region
  came up as BLOCKED cells.
- **Source:** `src/make_stl.py build_spouted_bed_stl_implicit`; commit
  `88b67d9` "fix(STL): correct normal direction for MFiX internal flow".

---

## Quick parameter reference (copy-pasteable)

Ground truth from `params/geometry.yaml` + `params/particles.yaml` as of
Day-10 R4 commit and later. All SI units.

```yaml
# Geometry (ORNL/UTK Case V, 60-deg cone + 4 mm orifice)
D_c_m: 0.050            # column inside diameter    [m]
D_i_m: 0.004            # orifice diameter          [m]
R_c_m: 0.025            # column radius             [m]
R_i_m: 0.002            # orifice radius            [m]
cone_angle_deg: 60      # 30 deg half-angle

# Mesh (Rule 1 + Rule 2)
dx_uniform_m: 1.0e-3    # uniform Cartesian cell size, dx=dy=dz
                        # 4 cells across D_i = 4 mm orifice
dx_sample_m:  5.0e-4    # STL SDF sampling grid (0.5 * dx)

# Box
x_min_m: -0.030, x_max_m: +0.030   # +/- (R_c + 5 mm margin)
z_min_m: -0.030, z_max_m: +0.030
y_min_m: -0.020                     # = -L_stub
H_dom_m:  0.200                     # ~3 * H_static + fountain headroom

# Stub + taper (Rule 3 + Rule 4)
R_stub_m:        2.0e-3   # = R_i (no widening)
L_stub_m:        2.0e-2   # = 10 * R_i
taper_length_m:  0.0      # no taper (straight pipe)

# Bed layout (Rule 6)
H_static_m:        0.0322         # derived: m_bed / (rho_p * eps_s * A_c)
eps_s_bed:         0.57           # IC packed-bed solids fraction
R_bed_gas_col_m:   0.004          # = max(R_stub, R_i + dx) = R_i + dx

# Operating point
U_in_m_s:          30.0           # orifice-plane velocity
bc_v_g_mi_m_s:     30.0           # = U_in (Rule 3: no scaling)
rho_g_kg_m3:       1.179          # air at 25 C
mu_g_pa_s:         1.854e-5

# Particles (Case V ZrO2)
d_p_m:             500.0e-6
rho_p_kg_m3:       6050
inventory_g:       54.5
c_e:               0.9            # = des_en_input(1) in DEM
ep_star:           0.42           # = des_epg_clip

# Timestep (Rule 9)
dt_s:      1.0e-5
dt_min_s:  1.0e-7
dt_max_s:  1.7e-5                  # <= 0.5 * dx / U_in
dt_fac:    0.9

# PIC closure (Rule 10)
friction_model:       'SCHAEFFER'
fric_exp_pic:         3.0
fric_non_sing_fac:    1.0e-7
psfac_fric_pic:       100.0
mppic_coeff_en1:      0.85
mppic_coeff_en_wall:  0.85
mppic_coeff_et_wall:  1.0
mppic_velfac_coeff:   1.0
des_epg_clip:         0.42
ic_pic_const_statwt:  4.2         # bed parcel stat weight
bc_pic_mi_const_statwt: 5.0       # MI sentinel (no parcels enter)

# DEM/CGP closure (Rule 11)
kn:               1000.0
kn_w:             1000.0
kt_fac:           0.2857142857    # = 2/7
kt_w_fac:         0.2857142857
mew:              0.3
mew_w:            0.3
des_en_input_pp:  0.9
des_en_input_pw:  0.9
des_etat_fac:     0.5
des_etat_w_fac:   0.5
cgp_stat_wt:      10              # 1 for pure DEM
dem_dt_s:         2.0e-5          # 0.1 * Rayleigh dt_c at CGP=10

# Drag + interpolation (Rules 12 + 13)
drag_type:              'SYAM_OBRIEN'
drag_c1:                0.8
drag_d1:                2.65
des_interp_scheme:      'LINEAR_HAT'
des_interp_mean_fields: .True.
des_interp_on:          .True.
gener_part_config:      .True.

# STL (Rule 19 + Rule 20)
top_extend_m:  5.0e-3             # STL cap pushed above H_dom
# Normals: INWARD (see src/make_stl.py build_spouted_bed_stl_implicit)
```

---

## Known failure modes and what to check

Direct symptom -> rule to re-check matrix. If the smoke validator reports any
of these, do NOT raise `tstop` or submit to cloud; fix the deck.

| Observed symptom | Likely rule violated | What to check |
|---|---|---|
| Orifice-plane `v_y mean < 0.5 * U_in` | **Rule 1** (mesh) | `(2 * R_i) / dx >= 4`? If not, refine `dx_uniform_m`. |
| Jet dies at `y << H_static`; no fountain | **Rule 1 + Rule 6** | Mesh first; then quadrant IC. |
| Field-wide `max |V| > 2 * U_in` sustained | **Rule 6** or **Rule 9** | Did parcels seed over the orifice? Is `dt_max` CFL-safe? |
| EP_G = 0.333 clip near the stub-tube cells | **Rule 6** | Quadrant IC missing or `R_bed_gas_col` too small. |
| `MAX_INLET_VEL_FAC` clip events (Vg > 750 m/s) | **Rule 6 + Rule 9** | Both. Fix A (IC) + Fix B (dt_max). |
| CFL > 2 sustained | **Rule 9** | `dt_max <= 0.5 * dx / U_in`? |
| Parcel seeding < 50% of target | **Rule 6 + Rule 8** | Is `R_bed_gas_col` too wide? Is `ic_des_lattice` set? Is strip geometry pinched by cone walls? |
| MFiX error 1003 at init | **Rule 7** | IC 1 must span whole box with `ic_ep_g = 1`. |
| Near-axis cells show as BLOCKED | **Rule 2 + Rule 20** | Graded mesh? Outward STL normals? |
| Zero fluid cells in stub region | **Rule 1 + Rule 3** | `R_stub / dx` and `R_i / dx` must each be >= 2. |
| PO plane blocked (top cells = BLOCKED) | **Rule 19** | STL top cap at or below `H_dom`; push it up by `top_extend_m = 5 mm`. |
| Spout dies at y ~ 15 mm with the mesh "right" | **Rule 3** | Did R_stub get set to 2 * R_i somewhere? Check `R_stub_m` and verify `bc_v_g_mi = U_in` (not scaled). |
| TFM run collapses dt within 100 steps | **Rule 17** | TFM is not viable on cut-cell spouted beds; switch to PIC. |
| Deck parses but solver aborts at step 0 | **Rule 15** | Keyword not in `init_namelist.f`? E.g. `BLENDING_FUNCTION` is doc'd but rejected. |
| PSD peak at the first Welch bin (~0.98 Hz) | time-window too short | Not a deck bug — need longer `tstop` for `df <= 0.5 Hz` resolution. |

---

## Pre-cloud-submission gate

Before any AWS Batch / HPC submission:

1. **Deck must satisfy Rules 1-20 above.** Verify by grep on the rendered
   `.mfx`:
   ```
   grep -E "cartesian_grid|use_stl|drag_type|des_interp_scheme|friction_model" <deck>
   grep -E "bc_type\(3\)|bc_v_g\(3\)|bc_x_w\(3\)|bc_y_s\(3\)" <deck>
   grep -E "ic_ep_g\(1\)|ic_ep_s\(2,1\)|ic_ep_s\(5,1\)" <deck>
   ```

2. **Local smoke run** (`tstop = 0.3 s` is the Day-10 R4 convention). Produce
   `cases/<name>_smoke/run.log` plus BACKGROUND VTUs.

3. **Smoke validation** per `.claude/skills/mfix-smoke-validate/SKILL.md`.
   For a spouted-bed deck specifically, Phase 3 pass criteria are:
   - **Orifice `v_y mean >= 25 m/s`** at a target `U_in = 30 m/s`
     (= 83% of target; `v_y max <= 60 m/s`).
   - **Centerline `v_g` at y = 25 mm >= 10 m/s** (jet survives through the
     bed; vs <= 2.4 m/s observed on the Day-10 R3 pathological run).
   - **Fountain nucleation**: at least one parcel crosses
     `y = H_static_m` by the end of the smoke window.
   - **Field-wide `max |V| < 60 m/s`**, no sustained `EP_G < 0.42`
     anywhere (DEM) or `< 0.333` (PIC).
   - **CFL max < 2**, **NIT < 30** per step.

4. **If any criterion fails, STOP.** Diagnose per Section 7 and re-render
   the deck. Do NOT submit a failing smoke to cloud and hope averaging
   fixes it.

5. **Record** the submission in `results/run_log.csv` with case name, mesh,
   dt, key params, status, and wall-clock (CLAUDE.md rule 8).

---

## References

**Project files (ground truth for this recipe):**
- `/Users/nelsonpereira/Documents/GitHub/Personal_Projects/Nuclear/TRISO-CFD/params/geometry.yaml` (blocks `bed_V`, `spouted_bed_V_3d`, `spouted_bed_V_3d_pic`, `spouted_bed_V_3d_dem`, `operating_point_V`)
- `/Users/nelsonpereira/Documents/GitHub/Personal_Projects/Nuclear/TRISO-CFD/params/particles.yaml` (block `case_V`)
- `/Users/nelsonpereira/Documents/GitHub/Personal_Projects/Nuclear/TRISO-CFD/cases/V_3d/template.mfx.j2` (PIC deck)
- `/Users/nelsonpereira/Documents/GitHub/Personal_Projects/Nuclear/TRISO-CFD/cases/V_3d/template_dem.mfx.j2` (DEM/CGP deck)
- `/Users/nelsonpereira/Documents/GitHub/Personal_Projects/Nuclear/TRISO-CFD/src/make_case.py` (`build_V_3d_params`, `build_V_3d_dem_params`)
- `/Users/nelsonpereira/Documents/GitHub/Personal_Projects/Nuclear/TRISO-CFD/src/make_stl.py` (`build_spouted_bed_stl_implicit`, inward normals)

**Project notes (iteration history; read these before deviating):**
- `/Users/nelsonpereira/Documents/GitHub/Personal_Projects/Nuclear/TRISO-CFD/notes/day10_summary.md` (R1 / R2 / R3 PIC iteration; Fix A)
- `/Users/nelsonpereira/Documents/GitHub/Personal_Projects/Nuclear/TRISO-CFD/notes/day10_pathology_review.md` (EP_G clip + LINEAR_HAT diagnosis)
- `/Users/nelsonpereira/Documents/GitHub/Personal_Projects/Nuclear/TRISO-CFD/notes/day11_summary.md` (DEM/CGP template rationale; R3c inventory deficit)
- `/Users/nelsonpereira/Documents/GitHub/Personal_Projects/Nuclear/TRISO-CFD/notes/day10_summary.md R4` (mesh refinement + straight stub; retires R2 Option 4)
- `/Users/nelsonpereira/Documents/GitHub/Personal_Projects/Nuclear/TRISO-CFD/results/run_log.csv` (every submission + cost)

**NETL reference tutorials (installed MFiX 26.1.2):**
- `/Users/nelsonpereira/mamba/envs/mfix-26.1.2/share/mfix/templates/tutorials/pic/spouted_bed_3d/spouted_bed_pic_3d.mfx` (PIC reference — IC pattern + closure)
- `/Users/nelsonpereira/mamba/envs/mfix-26.1.2/share/mfix/templates/tutorials/dem/spouted_bed_pseudo_2d/spouted_bed_dem_pseudo_2d.mfx` (DEM reference — contact closure)
- `/Users/nelsonpereira/mamba/envs/mfix-26.1.2/share/mfix/templates/tutorials/cgp/fluidbed_cgp_3d/fluidbed_cgp_3d.mfx` (CGP reference — `cgp_stat_wt` usage)

**MFiX 26.1.2 keyword / parser sources:**
- `/Users/nelsonpereira/mamba/envs/mfix-26.1.2/share/mfix/doc/html/reference/` (HTML keyword reference, per section)
- `/Users/nelsonpereira/mamba/envs/mfix-26.1.2/share/mfix/src/model/init_namelist.f` (authoritative keyword list — if it is not here, the parser rejects it)
- `/Users/nelsonpereira/mamba/envs/mfix-26.1.2/share/mfix/src/cartesian_grid/cartesian_grid_init_namelist.f` (cut-cell / STL keywords)
- `/Users/nelsonpereira/mamba/envs/mfix-26.1.2/share/mfix/src/check_data/` (`check_ic_common_discrete.f` for error 1003; `check_bc_*.f` for BC validation)

**Project rules (always apply):**
- `/Users/nelsonpereira/Documents/GitHub/Personal_Projects/Nuclear/TRISO-CFD/.claude/CLAUDE.md`
  rules 1 (no invented keywords), 2 (params in yaml), 3 (SI), 8 (log every
  run), 9 (explain before changing), 10 (use `.venv`), 13 (dry-run locally),
  14 (set generous timeouts).

**Sibling skill:**
- `/Users/nelsonpereira/Documents/GitHub/Personal_Projects/Nuclear/TRISO-CFD/.claude/skills/mfix-smoke-validate/SKILL.md`
  (always run before any cloud submission).

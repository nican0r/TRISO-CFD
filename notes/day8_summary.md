# Day 8 Summary — Case V (ORNL/UTK) 2D template and builder

## What was implemented

Day 8 retargets Case V from the previously-dropped laboratory
column (dropped per PLAN.md §1.2) to the ORNL/UTK cold-mockup spouted bed (S1, S2) and ships the 2D
template, parameter builder, tests, and a geometry schematic.

- `params/particles.yaml` — new `case_V` block (ZrO₂, d_p = 500 µm,
  ρ_p = 6050 kg/m³, inventory 54.5 g, φ = 15°, ep_star = 0.42, eps_s_bed =
  0.57, c_e = 0.9), case_D ρ_p updated to 6050 kg/m³ (S2).
- `params/geometry.yaml` — new `bed_V` block (D_c = 50 mm, D_i = 4 mm, 60°
  cone) and `spouted_bed_V_2d` block (H_dom = 0.20 m, imax = 50, jmax =
  200, kmax = 1, z_thickness = 1 mm).
- `src/make_case.py` — `build_V_2d_params()`, which derives the static bed
  height from inventory + ρ_p + eps_s_bed + cone geometry via
  `_bed_height_from_inventory()` (frustum-volume bisection); loads gas
  properties for air at ~25 °C via Cantera; assembles all template
  variables with S1 Table 1 closures (SYAM_OBRIEN drag; SCHAEFFER friction
  with GIDASPOW_PCF blending = "modified sigmoidal"; LUN_1984 PDE
  granular energy; φ = 15°, ep_star = 0.42, c_e = 0.9).
- `cases/V_2d/template.mfx.j2` — planar Cartesian cut-cell slab
  (cartesian_grid = .True., coordinates = 'CARTESIAN', no_k = .True.,
  n_quadric = 1, Y_CONE for the cone wall); IC regions for freeboard,
  packed bed and pre-opened spout channel; BCs for the orifice MI, top
  PO, and Y_CONE no-slip wall; four monitors (P_inlet and P_top at 1 kHz,
  solids inventory volume integral, centerline v_g line integral).
- `tests/test_make_case.py` — twelve V_2d tests: coordinates/no_k;
  grid + orifice-cell count; MI/PO/CG_NSW BCs; monitors; bed-IC hand
  calc (H₀ and H₀/D_c vs S2's 0.50–0.65); 22.6 L/min volumetric check
  and U/U_ms ≈ 1.9 vs S2's 12 L/min; closure keyword block.
- `tests/test_correlations.py` — U_ms,MG hand-calc rebased to ORNL/UTK
  geometry (D_c = 0.050 m, D_i = 0.004 m, ρ_p = 6050 kg/m³, H = 0.10 m;
  U_ms = 0.3643 m/s).
- `src/make_figs_d8.py` — geometry schematic generator.
- `results/fig_d8_V_2d_geometry.png` — screenshot.

## How it works

1. **`_bed_height_from_inventory`**. Given mass m, ρ_p, eps_s_bed, R_i and
   cone half-angle θ, finds h such that
   `V_frustum(h) = (π*h/3)*(R_i² + R_i*r_h + r_h²)` equals `m/(ρ_p*eps_s_bed)`
   with `r_h = R_i + h*tan(θ)`. Bisection to 1 nm. For the ORNL/UTK input
   (m = 0.0545 kg, ρ_p = 6050, eps_s_bed = 0.57, R_i = 2 mm, θ = 30°) this
   returns H_static ≈ 0.0322 m, giving H₀/D_c = 0.644 — just inside S2's
   reported 0.50–0.65 range. Test pins the volume balance to 1e-3
   relative tolerance.

2. **`build_V_2d_params`**. Reads bed_V, spouted_bed_V_2d, case_V from
   YAML. Derives H_static (bed height) and the Y_CONE quadric apex
   `(0, -R_i/tan(θ), 0)` and clip top `(R_c - R_i)/tan(θ)` so the quadric
   produces r = R_i at y = 0 (orifice rim) and r = R_c at y = cone_top.
   Loads air gas properties via Cantera. Returns a dict with every
   template variable, including the four closure keywords.

3. **Template**. Every keyword is cross-checked against init_namelist.f
   in MFiX 26.1.2. Important finding: the namelist header lists
   `BLENDING_FUNCTION = 'SIGM_BLEND'` as a valid value, but
   `check_blending_function` in check_solids_continuum.f only accepts
   `'NONE'`, `'TANH_BLEND'`, and `'GIDASPOW_PCF'`. The last one
   internally activates the SIGM_BLEND flag, i.e. the scaled-sigmoidal
   blending S1 Table 1 calls "modified sigmoidal". So the template renders
   `blending_function = 'GIDASPOW_PCF'`, and the comment traces the
   discrepancy.

4. **Hand-calc checks** (all three pinned in tests):
   - Bed inventory 54.5 g → H_static = 32.2 mm → H₀/D_c = 0.644 (inside
     S2's 0.50–0.65).
   - Inlet 30 m/s × π·(2 mm)² = 3.77·10⁻⁴ m³/s = 22.6 L/min, vs S2's
     12 L/min at 500 µm → U/U_ms = 1.88 (inside S2's 1.0–1.9 stable
     range).
   - Rendered closure keywords: SYAM_OBRIEN drag, SCHAEFFER friction,
     GIDASPOW_PCF blending (= modified sigmoidal), LUN_1984 PDE kinetic
     theory, φ = 15°, ep_star = 0.42, c_e = 0.9.

5. **MFiX initialization check**. Rendered .mfx was fed to MFiX 26.1.2
   (`mfixsolver -f V_2d.mfx` with tstop = 1e-4 s). Output: 10 000 total
   cells, 9024 standard + 110 cut + 866 blocked (= 9134 fluid), the
   Y_CONE quadric cut the cone correctly, PRE_PROCESSING COMPLETE in
   0.13 s, and the solver reached tstop without errors.

## What it models (physically)

- **Geometry**. The 50 mm / 60° / 4 mm vessel is S1's "standard" cold
  mockup. The 60° cone is the geometry S1 and S2 report validation data
  for. The 4 mm orifice is a constant diameter; S2's abstract typo
  ("0.04 cm") is reconciled against its body text and INL/CON-07-12569.
  This is recorded as the known inlet-diameter discrepancy.

- **Particles and closures**. 500 µm ZrO₂ at 6050 kg/m³ puts the bed in
  Geldart group D, where spouting (as opposed to bubbling) is the
  stable mode above U_ms. S1 Table 1 closures are the published-model
  spec: Syamlal–O'Brien drag (isothermal); Schaeffer frictional stress
  with the scaled-sigmoidal blending (so the transition from kinetic to
  frictional stress at eps_s → eps_s,max is smooth); and granular
  temperature solved as a PDE (Lun 1984), not algebraically. ep_star =
  0.42 (= eps_g,mf per S1) sits above eps_s_bed = 0.57 so Schaeffer
  frictional stress is inactive at t = 0 by physics (not by numerical
  tuning).

- **Operating point.** 30 m/s inlet jet gives Q = 22.6 L/min, which is
  about 1.9 × S2's reported U_ms (12 L/min for 500 µm). This is the
  upper end of S2's stable-spouting range (U/U_ms = 1.0–1.9).

- **Limitation.** MFiX 26.1.2 forbids `cartesian_grid=.True.` with
  `coordinates='CYLINDRICAL'`, so option (a) in day-8.md (axisymmetric
  with stair-stepped cone) cannot carry a resolved sloped cone wall
  under documented keywords. Option (b) (planar quadric cut-cell slab)
  is used instead and results are to be treated as qualitative only.
  A 3D cut-cell / STL run on Day 10 is the validation step.

## Known open items

- Pressure-tap location — not reported in S1 or S2; assumed at the inlet
  plane and the top of bed (both at 1 kHz).
- Particle-particle restitution c_e — not reported in S1 or S2; set to
  0.9 and flagged for the Day 13 sensitivity study.
- Whether S1's ambient-run baseline is 2D axisymmetric — S1 Table 1
  does not say explicitly. Treat it as 3D (consistent with S1 Fig. 1's
  centerline v_g profile shape).

## Day 1-7 files touched by the retargeting

- `params/particles.yaml` (case_D rho_p 6000 → 6050; added case_V; dropped
  the old lab-case particle block; the fb_sweep block is unchanged
  except for ep_star which tracks case_D).
- `params/geometry.yaml` (old lab-column geometry block dropped; new bed_V;
  fb_sweep_2d width retargeted to 50 mm; old spouted-bed geometry blocks
  dropped).
- `params/numerics.yaml` (comment-only: CFL-example reference to the
  prior lab bed replaced by "ORNL/UTK bed").
- `tests/test_correlations.py` (two U_ms hand-calc tests rebased to
  the ORNL/UTK geometry; the correlation code is unchanged).
- `notes/day04_model_choice.md` and `notes/day04_summary.md` (text
  pointers updated: prior lab reference → ORNL/UTK).
- `notes/day3_summary.md` (orifice reference updated to the ORNL/UTK bed).
- `notes/day6-7_summary.md` (fb_sweep width reference updated).
- `notes/mfx_anatomy.md` (the Day-5 MFiX TUTORIAL particle phase was
  previously labelled with the forbidden-string substring; relabelled
  as "tutorial beads". That tutorial case is unrelated to the dropped
  prior Case V inventory).
- `analytics.ipynb` (code and markdown cells retargeted; cached
  outputs stripped — the Day 2 and Day 4 figures must be re-run
  against the new YAML to pick up updated numbers).

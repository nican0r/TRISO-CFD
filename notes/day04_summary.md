# Day 4 — Summary

## What was implemented

A **modelling-method decision** for the two-week plan: which gas–solid approach
(TFM, DEM, CGP-DEM, PIC) is used where, pinned to a concrete particle-count
estimate for the ORNL/UTK 0.050 m bed loaded with Case D (500 µm, 6050
kg/m³ — see `params/particles.yaml` `case_D`) at a static bed height
H = 0.10 m. The particle counts are computed in the notebook so they stay in
sync with the geometry and particle YAMLs — and they are the numbers that make
the TFM-vs-DEM cost argument quantitative rather than hand-wavy.

Key artifacts:

- `notes/day04_model_choice.md` — the one-page decision table (TFM / DEM /
  CGP-DEM / PIC × outputs / cost / limitations / when-used-here).
- `analytics.ipynb` — Day-4 cells that read `params/geometry.yaml` and
  `params/particles.yaml` and print the full-3D and pseudo-2D-slab particle
  counts.
- `params/particles.yaml` — new `day4_count_estimate` block (H_static_m = 0.10,
  eps_s = 0.6, slab_thickness_particles = 5) with source comments pointing at
  `implementation-steps/day-4.md`.

No code was added to `src/`, no new correlations, no simulation runs. This step
is a quantitative scoping decision.

## How it works

Note on geometry: this section was originally written against a pre-retargeting
column (D_c = 76 mm, D_i = 9.5 mm). The project geometry has since been
retargeted to the ORNL/UTK column — **D_c = 50 mm, D_i = 4 mm, 60° cone**
(see `params/geometry.yaml` `bed_V`). The particle counts below have been
recomputed for the current geometry; the old numbers (N_3D ≈ 2.67 × 10⁶,
N_slab ≈ 1.30 × 10⁵) are preserved in the decision table of
`notes/day04_model_choice.md` as a historical note.

1. **Cone frustum geometry.** The bed's cone section runs from the inlet
   orifice (D_i = 4 mm) up to the column (D_c = 50 mm) at an included cone
   angle of 60°, so the half-angle is 30°. The cone's axial height is

   `h_cone = (r_c − r_i) / tan(30°) ≈ 39.84 mm`  [m].

   For H_static = 100 mm, the bed fills the whole cone and extends 60.16 mm
   into the cylindrical section.

2. **Bed volume (SI).**
   - Cone frustum volume: `V_cone = π·h_cone/3 · (r_i² + r_i·r_c + r_c²)
     ≈ 2.84 × 10⁻⁵ m³`.
   - Cylinder volume above the cone: `V_cyl = π·r_c² · (H − h_cone)
     ≈ 1.18 × 10⁻⁴ m³`.
   - Total: **V_bed ≈ 1.47 × 10⁻⁴ m³**.

3. **Solids volume and particle count.** With ε_s = 0.6:
   `V_s = ε_s · V_bed ≈ 8.80 × 10⁻⁵ m³`.
   Single particle: `V_p = π/6 · d_p³ ≈ 6.54 × 10⁻¹¹ m³` for d_p = 500 µm.
   Count: **N_3D = V_s / V_p ≈ 1.34 × 10⁶ particles** (≈ 1.3 million).

4. **Pseudo-2D slab.** Thickness `t = 5·d_p = 2.5 mm`. The in-plane projected
   area through the bed axis is a trapezoid (cone) plus a rectangle (cylinder):

   `A_proj = h_cone·(D_i + D_c)/2 + (H − h_cone)·D_c ≈ 4.08 × 10⁻³ m²`.

   Slab volume `V_slab = A_proj · t ≈ 1.02 × 10⁻⁵ m³`, giving
   **N_slab ≈ 6.7 × 10⁴ particles** — a **20× collapse** from full 3D.

5. **Verification.** The notebook's Day-4 cell reads geometry and particle
   parameters from YAML (no magic numbers; inputs under
   `params/particles.yaml` `day4_count_estimate` plus `params/geometry.yaml`
   `bed_V`) and prints both counts. Any future revision of H, ε_s, or
   geometry is picked up automatically.

## What it models (physically)

Day 4 does not add physics — it decides which physics engine runs for each
downstream question. Three operating decisions in PLAN §1.1 fall out of this:

- **"At what gas flow does the bed spout stably?"** → **TFM** is sufficient:
  pressure drop, bed expansion, and the spout–annulus structure are mean-flow
  observables. TFM's cost is set by mesh and CFL (Day 3: dt ≤ 25 µs at a 1 mm
  jet mesh), not by particle count — so it is the right tool for the Days-6–7
  flow sweep, the Case V validation runs, and the Case D sensitivity sweeps
  (gas T, drag model, particle density). PLAN §1.4 "parametric study" leans
  entirely on TFM.
- **"What is the solids circulation pattern and roughly how long is one
  cycle?"** → **TFM** gives the mean velocity field; a cycle time is `H_bed /
  u_s_annulus_mean`. Individual-particle histories would be nicer but are not
  required for the mean cycle.
- **"How uniform is particle exposure in the spout + fountain?"** → **DEM**,
  but only in a **pseudo-2D slab**. Full 3D DEM would need ≈ 1.34 M particles
  at a τ_c-bounded sub-step — out of reach on a workstation for a Geldart-D
  system. The 20× geometric collapse of the slab lands at ≈ 67 k particles,
  which is in the band where MFiX-DEM is routinely run. This was intended as
  the single CFD-DEM entry in PLAN §1.2.

Update after Days 10–11: the plan evolved. Full-3D **MP-PIC** was ultimately
adopted as the Day-10 baseline (`cases/V_3d/day10/` with the
`spouted_bed_V_3d_pic` closure block in `params/geometry.yaml`), and Days
11–12 added a full-3D **CGP-DEM** variant at statistical weight ~10
(`spouted_bed_V_3d_dem`). The pseudo-2D-DEM slab argument above is still the
right scoping lens, but the actual runs that produced the Case V particle-
level statistics are 3D PIC/CGP-DEM, not slab DEM.

The particle count also tells us what the slab run can and cannot resolve:
with ~67 k particles distributed between spout (narrow, fast), fountain
(sparse, intermittent), and annulus (dense, slow), the fountain is the region
with the fewest statistically-independent trajectories and will dominate the
uncertainty in the per-particle residence-time distribution. That drives the
averaging-window choice in later days (CLAUDE.md rule 6).

# Day 4 — Summary

## What was implemented

A **modelling-method decision** for the two-week plan: which gas–solid approach
(TFM, DEM, CGP-DEM, PIC) is used where, pinned to a concrete particle-count
estimate for the Missouri S&T 0.076 m bed loaded with Case D (500 µm, 6000
kg/m³) at a static bed height H = 0.10 m. The particle counts are computed in
the notebook so they stay in sync with the geometry and particle YAMLs — and
they are the numbers that make the TFM-vs-DEM cost argument quantitative rather
than hand-wavy.

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

1. **Cone frustum geometry.** The bed's cone section runs from the inlet
   orifice (D_i = 9.5 mm) up to the column (D_c = 76 mm) at an included cone
   angle of 60°, so the half-angle is 30°. The cone's axial height is

   `h_cone = (r_c − r_i) / tan(30°) ≈ 57.59 mm`  [m].

   For H_static = 100 mm, the bed fills the whole cone and extends 42.41 mm
   into the cylindrical section.

2. **Bed volume (SI).**
   - Cone frustum volume: `V_cone = π·h_cone/3 · (r_i² + r_i·r_c + r_c²)
     ≈ 9.93 × 10⁻⁵ m³`.
   - Cylinder volume above the cone: `V_cyl = π·r_c² · (H − h_cone)
     ≈ 1.92 × 10⁻⁴ m³`.
   - Total: **V_bed ≈ 2.92 × 10⁻⁴ m³**.

3. **Solids volume and particle count.** With ε_s = 0.6:
   `V_s = ε_s · V_bed ≈ 1.75 × 10⁻⁴ m³`.
   Single particle: `V_p = π/6 · d_p³ ≈ 6.54 × 10⁻¹¹ m³` for d_p = 500 µm.
   Count: **N_3D = V_s / V_p ≈ 2.67 × 10⁶ particles** (≈ 2.67 million).

4. **Pseudo-2D slab.** Thickness `t = 5·d_p = 2.5 mm`. The in-plane projected
   area through the bed axis is a trapezoid (cone) plus a rectangle (cylinder):

   `A_proj = h_cone·(D_i + D_c)/2 + (H − h_cone)·D_c ≈ 5.69 × 10⁻³ m²`.

   Slab volume `V_slab = A_proj · t ≈ 1.42 × 10⁻⁵ m³`, giving
   **N_slab ≈ 1.30 × 10⁵ particles** — a **20× collapse** from full 3D.

5. **Verification.** The notebook's Day-4 cell reads geometry and particle
   parameters from YAML (no magic numbers) and prints both counts; the printed
   values (`2.674e+06` and `1.303e+05`) agree with the hand calc above to the
   third significant figure. Any future revision of H, ε_s, or geometry is
   picked up automatically.

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
  but only in a **pseudo-2D slab**. Full 3D DEM would need ≈ 2.67 M particles
  at a τ_c-bounded sub-step — out of reach on a workstation for a Geldart-D
  system. The 20× geometric collapse of the slab lands at ≈ 130 k particles,
  which is in the band where MFiX-DEM is routinely run. This is the single
  CFD-DEM entry in PLAN §1.2.

CGP-DEM and MP-PIC are **not used** in the two-week plan; they are listed in
the decision table as the natural upgrade paths if the pseudo-2D-DEM run turns
out to be insufficient and HPC becomes available (future work for REPORT.md,
PLAN §1.4).

The particle count also tells us what the DEM slab run can and cannot resolve:
with ~130 k particles distributed between spout (narrow, fast), fountain
(sparse, intermittent), and annulus (dense, slow), the fountain is the region
with the fewest statistically-independent trajectories and will dominate the
uncertainty in the per-particle residence-time distribution. That drives the
averaging-window choice in later days (CLAUDE.md rule 6).

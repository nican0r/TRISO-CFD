# Day 4 — Which gas–solid model, and when

Decision table for picking a modelling approach for the 0.050 m ORNL/UTK
cold-flow spouted bed loaded with Case D (500 µm, 6000 kg/m³, Geldart D) at a
static bed height H = 0.10 m.

Particle-count anchors (computed in `analytics.ipynb`, Day-4 cells):

- Full 3D bed, H = 0.10 m, ε_s = 0.6: **N_3D ≈ 2.67 × 10⁶ particles**
  (V_bed ≈ 2.92 × 10⁻⁴ m³, V_s ≈ 1.75 × 10⁻⁴ m³, V_p ≈ 6.54 × 10⁻¹¹ m³).
- Pseudo-2D slab, thickness t = 5·d_p = 2.5 mm (axial cross-section):
  **N_slab ≈ 1.30 × 10⁵ particles** (collapse factor ≈ 20× vs 3D).

| Axis | TFM (two-fluid, Eulerian–Eulerian) | DEM (resolved particles, CFD-DEM / MFiX-DEM) | CGP-DEM (coarse-grained parcels) | PIC (MP-PIC, e.g. MFiX-PIC) |
|---|---|---|---|---|
| **What it outputs** | Phase-averaged fields: ε_g(x,t), u_g(x,t), ε_s(x,t), u_s(x,t), p(x,t), granular T. One velocity per phase per cell — no individual particle tracks. | Position, velocity, contact force for every real particle; gas fields on the CFD mesh. Trajectories, collisions, residence-time distributions directly. | Same as DEM but each "parcel" lumps k real particles sharing a position/velocity, with contact forces scaled consistently. Trajectories per parcel, not per particle. | Parcels advected by gas with drag + an interpolated solids-stress gradient (Harris–Crighton closure). No pairwise contacts; parcel positions and velocities, cell-averaged solids volume fraction. |
| **Cost for this bed (H = 0.10 m, N_3D ≈ 2.67 M)** | Cheapest — grid-only. For a 1 mm fluid mesh over the bed+freeboard it is O(10⁵–10⁶) cells; dt bounded by CFL (Day 3: ≤ 25 µs at u_jet = 20 m/s). Hours on a workstation per 10 s of physical time. | Most expensive. 2.67 M particles with contact detection at a particle-scale sub-step (τ_c ≪ dt_fluid for Geldart-D). Current MFiX-DEM practice tops out around 10⁵–10⁶ particles on a workstation → full 3D run is marginal / needs HPC. The **pseudo-2D slab (≈ 1.3 × 10⁵ particles) is tractable** and is why we use it. | Middle: pick a parcel size k so N_parcel ≈ 10⁵–10⁶. For k = 20 the full 3D run becomes ≈ 1.3 × 10⁵ parcels — same order as the pseudo-2D DEM slab but keeping the full geometry. Contact mechanics is still resolved, just at parcel scale. | Cheapest particle-based option — pairwise contacts replaced by a solids-stress gradient, so cost scales ~linearly in number of parcels and is largely insensitive to packing. Full 3D with ~10⁵–10⁶ parcels runs in hours on a workstation. |
| **What it can't tell us** | Individual-particle residence times, individual-particle histories through spout/fountain/annulus (the uniformity-of-coating question in PLAN §1.1). Collision statistics. Any observable that needs tagging specific particles. Multiple particle sizes/densities only through extra solid phases. | Nothing fundamental is excluded at the particle level, but at Geldart-D scale a resolved-particle run on a workstation **cannot afford full 3D for this geometry** — a hard practical limit, not a modelling one. | Smooths out collision-scale statistics (two parcels colliding ≠ two particles colliding). Fine for mean flow and bulk circulation, less reliable for narrow RTD tails and segregation at the particle scale. | No real contact mechanics → cannot resolve dense-packing stresses accurately, spout–annulus interface sharpness is weaker, friction angle has to be prescribed. Not a good tool for "does it spout stably at this U?" near the Umsf boundary. |
| **When we use it in this plan** | **Day 5 tutorial (`tut_tfm2d`)** to learn MFiX; **Days 6–7 flow-rate sweep (`fb_sweep`)**; the **Case V validation runs** and the **Case D parametric sweeps** (gas T, particle ρ, drag model) in PLAN §1.4. TFM is the workhorse — one velocity field per phase is all we need for pressure drop, spouting status, mean circulation, and the correlation comparisons. | **One small CFD-DEM run on the pseudo-2D slab** (PLAN §1.2 "a small CFD-DEM run for particle-level statistics") to get per-particle residence-time distributions in spout + fountain. Not used elsewhere — too expensive. | Not used in the two-week plan. Listed here because it is the natural upgrade path if the DEM pseudo-2D slab turns out insufficient and HPC becomes available (noted as future work in REPORT.md). | Not used in the two-week plan. A reasonable fallback if DEM-pseudo-2D is still too slow, but it would weaken the particle-level statistics we actually need from Day-14. |

## Why these cost numbers are what they are

- **DEM cost scales with N_particles and 1/τ_c.** With N_3D ≈ 2.67 M particles and a Hertz-contact time of order 10 µs for 500 µm, 6000 kg/m³ spheres, DEM for the full bed is in the ≥10⁷ particle-step / fluid-step range — out of reach on a single workstation but tractable on HPC. The **pseudo-2D slab brings this down by the 20× geometric collapse factor**, so N_slab ≈ 1.3 × 10⁵ is exactly the band where MFiX-DEM is routinely run.
- **CGP-DEM parcel size k is the one free knob.** Picking k so N_parcel is in [10⁵, 10⁶] recovers DEM-like behaviour at a tractable cost; this is why it is listed as the natural upgrade path. We do not use it in the two-week plan to keep validation apples-to-apples with TFM and the slab-DEM reference.
- **TFM has no particle count — only cells.** Its cost is set by mesh refinement and CFL, not by N_particles, which is why it is the workhorse for sweeps.
- **PIC's dt is not bounded by τ_c** (no pairwise contacts), only by CFL. That is where its speedup comes from, and also why it is weakest at the spouting-onset boundary where contact stresses matter.

Numeric anchors are reproduced in the notebook so any change in geometry, H, or
ε_s is picked up without editing this note.

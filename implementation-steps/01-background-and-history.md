# 1. Background and project history

## The physical problem

TRISO fuel particles are uranium kernels (~0.4–0.5 mm) coated with four
layers (porous carbon buffer, inner pyrolytic carbon, silicon carbide, outer
pyrolytic carbon) by chemical vapor deposition in a **spouted bed**: a
conical vessel with a narrow gas nozzle at its apex. The gas jet carries
particles up a central **spout**, throws them into a **fountain** above the
bed, and they rain back onto the **annulus**, which slides slowly down the
cone wall to the nozzle again. Coating uniformity depends on how evenly
every particle cycles through the hot spout and fountain, so the bed
hydrodynamics (spout stability, circulation rate, residence times) are the
quantities of interest.

This project builds that hydrodynamic model in MFiX (NETL's open-source
multiphase CFD code). It is cold flow only: no chemistry, no heat transfer.
The goal and scope are in [PLAN.md](../PLAN.md).

## Stages

The project ran as numbered days (the original task specs are in git
history). Days 1–9 are foundations; Days 10–13 tried to model the
ORNL/UTK cold-mockup coater ("Case V") and did not succeed; the current
work is the C26 benchmark (files 02–08).

| Stage | What it established | Status |
|---|---|---|
| Day 1 | `src/gasprops.py` (Cantera gas properties), terminal velocity with Schiller–Naumann drag. 500 µm, 6000 kg/m³ in 20 °C air: u_t = 6.44 m/s by hand, 6.442 m/s by code | done, unit-tested |
| Day 2 | `src/correlations.py`: Archimedes, Ergun, Wen–Yu U_mf, Geldart class, Mathur–Gishler U_ms. Surrogate and UO₂ kernels are Geldart D | done, unit-tested |
| Day 3 | `src/cfl_check.py`. 20 m/s on 1 mm cells needs dt ≤ 2.5e-5 s for CFL 0.5. MFiX is semi-implicit, so this is a loose bound, but Day 10 showed an uncapped dt does break spout runs | done |
| Day 4 | Model choice (TFM vs DEM vs coarse-grained DEM vs PIC). 3D bed ≈ 1.3 × 10⁶ particles, too many for plain DEM | superseded at Day 10 |
| Day 5 | Readers for MFiX monitors and VTK; ΔP(t) and time-averaged ε_s on the MFiX 2D TFM tutorial; every tutorial keyword checked against the reference | done |
| Days 6–7 | 2D flat-bottom TFM bubbling bed, U/U_mf = 0.3–2.0 (`cases/fb_sweep/`). Initial ε_s lowered to 0.55 because Schaeffer friction blew up at the IC step at 0.60; the fine-mesh check was dropped for the same reason | run log "partial"; acceptance checks not recorded |
| Day 8 | Case V 2D template: ORNL/UTK mockup (50 mm column, 60° cone, 4 mm orifice, 54.5 g of 500 µm ZrO₂, 30 m/s air jet). MFiX forbids cut cells in cylindrical coordinates, so a planar slab was used | infrastructure only |
| Day 9 | 2D TFM smoke to 3 s: spout, fountain and annulus visible; inventory 94.1 % (floor 95 %); slab PSD not physical | qualitative only |
| Day 10 | 3D Case V: TFM, then MP-PIC, rounds R1–R5 (below) | never validated |
| Days 11–13 | CGP-DEM residence times and the Case D parameter sweep; on branch `case-v-d-day10-13`, not merged | not completed |

## Day 10: what the 3D spouted bed taught

Day 10 is where most of the lasting modeling rules came from. These
became the project's "spouted-bed recipe" (the
`.claude/skills/mfix-spouted-bed-recipe` skill), and C26 departed from
some of them on purpose (file 02).

- **3D cut-cell TFM failed.** On an STL vessel MFiX's cut-cell classifier
  blocked near-axis cells (4,556 fluid cells where ~8,600 were expected),
  and the quadric version collapsed dt through granular-energy PDE
  stiffness. Day 10 switched to MP-PIC. CGP-DEM later proved cleaner still
  under the same setup.
- **Uniform Cartesian meshes only.** Graded meshes and dense STLs trigger
  the same blocking.
- **Inlet.** Mass inlet over the full box bottom with an STL stub pipe
  (length 10 R_i, radius exactly R_i, the physical 30 m/s), at least 4 cells
  across the orifice. Widening the stub and scaling the velocity down, or
  tapering it, kills the jet.
- **STL seams.** SDF primitives that meet exactly on a sampled plane leave a
  disc of facets MFiX treats as a wall (the R5 "membrane" across the
  orifice); primitives must overlap.
- **Initial condition.** Never seed parcels over the orifice. A packed bed
  next to an empty gas column is near-singular under the PIC stress model.
  NETL's empty central column is the bore of a draft tube this vessel does
  not have. Settle the bed first, then ramp the gas.
- **Time step.** dt_max ≤ 0.5·dx/U_in. A 1e-3 s cap on 1 mm cells with a
  30 m/s jet let the solver run at CFL ~30 and killed the jet; capping at
  5e-5 s let it penetrate.

Best result (R5d): the jet penetrated and a fountain erupted, but the
centerline gas velocity was 0.15–0.56× the S1 Pitot data. An R4 claim of a
"matching" ΔP was the wrong sign and is void.

## Why the project moved to the C26 benchmark

Two things ended Case V. First, its validation data was thin. An audit
(2026-10-16, on the side branch) left one measured target, 7 Pitot points
of S1 Fig. 1: the "17–25 Hz" pulsation band is not in S2, and S1 Fig. 2 is
a different rig and gas. Second, the model never matched. Conry et al.
2026 ("C26") instead reports PIV maps of the spout and annulus on a
benchtop TRISO-coater rig. It also describes the MFiX CFD-CGDEM model the
authors built of that rig and reports the model's scores against the
experiment. The model can be rebuilt from that description "as the paper
did" and scored the same way. Case C26 has been the project's focus since
2026-10-16.

Dead ends from this period (do not retry):

- TFM where cells are not ≫ d_p. The Missouri S&T glass-bead case had
  2.18 mm beads and a 9.5 mm orifice, and crashed at 0.072 s.
- 3D cut-cell TFM; graded meshes; a 2 mm mesh on a 4 mm orifice.
- A widened or tapered stub with scaled velocity; a widened gas column; a
  pre-opened spout IC; a packed IC beside empty cells.
- A slab PSD compared with experiment.
- Retuning PIC stress parameters away from NETL defaults.

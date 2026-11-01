# 2. The C26 benchmark and the model

## The experiment

Conry, Dal Forno Chuahy, Oyedeji, Finney, Heldt & López-Honorato (2026),
*Nucl. Eng. Des.* 456:115009 ("C26"), studied the ORNL benchtop
spout-fluidized bed, a cold-flow analogue of a TRISO coater:

- glass vessel: 60° cone (30° half-angle) opening into a 50 mm column,
  one 2.38 mm (3/32 in.) nozzle at the bottom of the cone;
- ultra-high-purity argon at room temperature through one mass-flow
  controller (slpm referenced to 25 °C, 1 atm);
- 60 g of yttria-stabilized zirconia (YSZ) particles, a TRISO-kernel
  surrogate: 411.6 µm mean diameter (normal PSD, σ = 12.2 µm), 6000 kg/m³;
- a side-view high-speed camera (430 fps) and PIVlab processing, giving
  time-averaged particle-velocity maps above the bed (Fig. 4) and the
  particle velocity at the wall (Fig. 6). The processing steps are from the
  same group's earlier paper (Conry et al. 2025, *Chem. Eng. Sci.*
  306:121197).

C26 also built an MFiX 23.2 CFD-CGDEM model and scored its maps against the
experiment with three image metrics: correlation coefficient (CC),
structural similarity (SSI) and mean squared error (MSE) (Table 3). This
project rebuilds that model and scores it the same way.

## The model

| Item | Value | Source |
|---|---|---|
| Solver | MFiX 23.2, CFD-CGDEM (coarse-grained DEM) | C26 (23.2 since 2026-10-24; 26.1.2 before, file 05) |
| Gas | 100 % Ar at 298 K: ρ 1.633 kg/m³, μ 2.302e-5 Pa·s (Cantera) | C26 p.4 |
| Flow | 10 / 12 / 14 slpm; nozzle velocity Q/A_nozzle = 37.4 / 44.9 / 52.4 m/s | C26 Table 2 |
| Particles | YSZ, d_p 411.6 µm, ρ_p 6000 kg/m³, measured PSD | C26 |
| Coarse graining | statistical weight W = 10 for every case: each parcel stands for 10 particles, d_parcel = d_p·W^⅓ = 0.887 mm, ~27,770 parcels for 60 g | C26 baseline W; W < 10 out of scope (user decision) |
| Drag | TGS (Tenneti, Garg & Subramaniam 2011), evaluated at the real particle diameter under CGP | C26 §2.2.1 |
| DEM contacts | spring-dashpot, k_n 1000 N/m, k_t/k_n 0.285, η_t/η_n 0.5, friction 0.1 (particle–particle and particle–wall), restitution 0.5 | C26 Table 1 |
| Mesh | uniform 2.1 mm, 24 cells across the column, 100 cells tall (H = 0.21 m), pressure outlet at the top | C26 §2.2.1 (file 03) |
| Gas inlet | a point source of gas at the bottom of the cone | C26 §2.2.1 (file 04) |
| Time | 3.0 s simulated; dt ≤ 1e-5 s (fluid), DEM sub-step 2.13e-6 s; averages over t = 2.4–3.0 s | C26 p.3 ("minimum flow time of 3.0 s"), p.9 (0.01 ms) |

All values live in `params/*.yaml` with their sources: `bed_C26`,
`spouted_bed_C26_3d` and `spouted_bed_C26_3d_dem` in `geometry.yaml`;
`case_C26_ysz` in `particles.yaml`; `ar_25C` in `gases.yaml`;
`spouted_bed_C26_3d_numerics` in `numerics.yaml`.

### Initial bed

One initial-condition region spans the full cross-section from y = 0 up.
In MFiX DEM, `ic_ep_s` only sets the parcel *count*. Parcels are placed on
a hexagonal lattice at 1.05 × the largest parcel diameter, filled from the
bottom up. The seeded bed therefore starts denser (ε_s ≈ 0.50) and shorter
than the IC box, and settles to ε_s 0.59–0.62 within 40–80 ms. All checks
use the settled static height H_static (34.8 mm for 60 g of YSZ), not the
IC box height.

**Bed-height discrepancy.** C26 Fig. 1 marks the bed as "~25 mm".
60 g of YSZ at ε_s = 0.57 in a 60° cone is ~35 mm, and the first smoke
settled at 34.5 mm. 25 mm would need a ~90° cone. The bracket in Fig. 1 is
schematic: the 2025 paper draws the same bracket for 15 g of SiC. The
model keeps the stated 60° cone and the stated 60 g.

### What is assumed (C26 does not say)

- The cone bottom: the cone narrows to the nozzle diameter (throat
  2.38 mm) and the gas enters only through the point source. How the
  source cell sits in the cone tip is the main open question (files 03–04).
- The point source's size and velocity, until the authors replied
  (2026-10-27; file 04).
- A gas-fraction floor `des_epg_clip = 0.42`. It is required on this mesh
  (file 05).
- The particle–mesh interpolation scheme (file 05).
- Every PIV setting (file 07).

## Deviations from the spouted-bed recipe

The Day-10 recipe (file 01) was written for a 4 mm orifice on 1 mm cells.
C26 sizes its cells *to* the nozzle (dx = 2.1 mm ≈ D_n = 2.38 mm), so
several recipe rules were replaced on purpose:

| Recipe rule | C26 model | Why |
|---|---|---|
| ≥ 4 cells across the orifice; stub pipe with a mass inlet | gas point source | at R_i/dx ≈ 0.6 a stub has no fluid cells; 4 cells across would need dx ≤ 0.6 mm, ~2.6 M cells (35× C26's mesh) |
| 4-quadrant IC with a central gas column | one IC region over the full cross-section | the quadrant IC copied a draft-tube bore this vessel lacks |
| Syamlal–O'Brien drag | TGS | C26 |
| recipe DEM values (e = 0.9, μ = 0.3) | C26 Table 1 (e = 0.5, μ = 0.1) | C26 |
| LINEAR_HAT interpolation | NONE (MFiX default) since 2026-10-27 | file 05 |

Kept from the recipe: the uniform mesh, a gas IC over the whole domain,
inward STL normals, an STL cap 5 mm above the pressure outlet, and the
dt ≤ 0.5·dx/U rule (C26's 1e-5 s is the tighter bound at every flow).

## The case matrix: c01–c09

`params/sweep_C26.yaml`. One factor at a time around **c02** (YSZ, 60 g,
12 slpm, the best-scored C26 case).

| Case | Factor | Bed | d_p [µm] | ρ_p [kg/m³] | Flow [slpm] | H_static [mm] | parcel [mm] |
|---|---|---|---|---|---|---|---|
| c01 | flow (C26 case 1) | YSZ, 60 g | 411.6 | 6000 | 10 | 34.8 | 0.89 |
| c02 | baseline (C26 case 2) | YSZ, 60 g | 411.6 | 6000 | 12 | 34.8 | 0.89 |
| c03 | flow (C26 case 3) | YSZ, 60 g | 411.6 | 6000 | 14 | 34.8 | 0.89 |
| c04 | coating stage | + buffer | 611.6 | 2524 | 12 | 47.9 | 1.32 |
| c05 | coating stage | + IPyC | 691.6 | 2332 | 12 | 49.7 | 1.49 |
| c06 | coating stage | + SiC | 761.6 | 2550 | 12 | 47.7 | 1.64 |
| c07 | coating stage | + OPyC (full TRISO) | 841.6 | 2381 | 12 | 49.2 | 1.81 |
| c08 | charge | YSZ, 40 g | 411.6 | 6000 | 12 | 30.2 | 0.89 |
| c09 | charge | YSZ, 80 g | 411.6 | 6000 | 12 | 38.6 | 0.89 |

- Only c01–c03 have experimental data (C26 Table 2 cases 1–3). c04–c09
  are this project's own cases. They replace C26's wall-friction,
  statistical-weight and steel-feedstock cases and are compared only with
  each other and with c02.
- Coating layers are applied to the YSZ core in deposition order, with
  thicknesses and densities from Marshall et al. 2017 Table 1: buffer
  100 µm at 1000 kg/m³, IPyC 40 µm at 1900, SiC 35 µm at 3200, OPyC 40 µm
  at 1900. The PSD width is assumed to carry through coating.
- Flow stays at 12 slpm and W at 10, so coated parcels are larger (up to
  1.81 mm in 2.1 mm cells) and coated beds are deeper than the 41.6 mm
  cone. Bed depth and particle size change together in the coating block.
  A bed that does not spout at 12 slpm is reported as a result.
- c09 seeds 76.1 g, not 80 g: MFiX snaps the IC top to a cell face.
  Summaries use the seeded mass.

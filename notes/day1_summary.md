# Day 1 — Summary

## What was implemented

Set up the project skeleton and delivered the first quantitative building block: **single-sphere terminal settling velocity** for the TRISO surrogate particles across the cold-flow and hot-gas conditions the coater will be run at.

Key artifacts:

- `src/correlations.py` — drag / settling helpers
- `src/gasprops.py` — temperature-dependent gas properties from Cantera
- `params/particles.yaml`, `params/gases.yaml` — inputs (no magic numbers in code)
- `tests/test_correlations.py` — unit tests pinned to a hand calc
- `analytics.ipynb` — table + plot of u_t vs particle diameter
- `results/fig_d1_terminal_velocity.png` — the Day 1 deliverable figure

## How it works

1. **Gas properties.** `gas_props(T, P, composition)` calls Cantera's gri30 mechanism to return the mixture density ρ [kg/m³] and dynamic viscosity µ [Pa·s] at a given temperature, pressure, and mole-fraction composition. This lets us treat Ar, H₂, N₂, O₂ mixtures at both room temperature and coating temperature (1400 °C) consistently, instead of using tabulated values.

2. **Drag law.** `cd_schiller_naumann(Re)` gives the drag coefficient of an isolated sphere: the standard Schiller–Naumann correlation for Re ≤ 1000, capped at the Newton-regime constant Cd = 0.44 for Re > 1000.

3. **Terminal velocity.** `terminal_velocity(d_p, ρ_p, ρ, µ)` solves the force balance drag = buoyancy-corrected weight, i.e. `u² · Cd(Re) = (4/3) · d_p · (ρ_p − ρ) · g / ρ`, using `scipy.optimize.brentq` on `u ∈ [1e-8, 100] m/s`. Because Cd depends on Re which depends on u, this is an implicit equation — brentq iterates to a fixed root.

4. **Relaxation time.** `relaxation_time(d_p, ρ_p, µ) = ρ_p d_p² / (18 µ)`, the Stokes-limit particle response time. Included as a reference time scale for later CFD-DEM work.

5. **Verification.** A hand calculation (`notes/day1_handcalc.md`) iterates the same force balance manually for a 500 µm, 6000 kg/m³ ZrO₂ surrogate in air at 20 °C, giving u_t ≈ 6.44 m/s. The solver returns 6.442 m/s — pinned as a unit test so any future change to the drag law is caught.

## What it models (physically)

A single, isolated sphere falling through a **quiescent** gas at steady state. u_t is a **material property** of the particle–gas pair — the speed at which drag balances net weight in still fluid. It is not a prediction of anything happening inside the coater.

The only physical variables that appear:

- particle: d_p, ρ_p
- gas: ρ, µ (from Cantera at the chosen T, P, composition)
- gravity: g

The only physical effects included: viscous + form drag (Schiller–Naumann Cd) and buoyancy-corrected gravity. No walls, no other particles, no gas flow field, no temperature gradient, no chemistry.

## What Day 1 actually delivers

One deliverable: `results/fig_d1_terminal_velocity.png` — u_t vs d_p (200–1200 µm) at three densities {2500, 6000, 10 800} kg/m³ across four gas conditions {air 20 °C, Ar 20 °C, Ar 1400 °C, Ar/H₂ 50/50 mol 1400 °C}, plus the underlying numeric table in the notebook.

That is the entire Day-1 output. It shows which particles are heavy/light relative to which gas and where the hot-gas curves sit relative to the cold-gas curves. Nothing else has been computed, and nothing has been compared against a simulation.

## What Day 1 does NOT model

- **No fountain height, no residence time.** Both require the gas velocity field above the bed (specifically the jet decay `u_gas(z)`), which isn't in this model. A uniform upward gas at U > u_t would give an infinite fountain — the real fountain is finite only because the jet decelerates with height as it entrains surrounding gas.
- **No bed / particle-particle interactions.** Hindered settling and packed-bed pressure drop start Day 2 (U_mf, U_ms).
- **No geometry.** Bed diameter, cone angle, inlet orifice — Day 2's `geometry.yaml`.
- **No turbulence, no heat transfer, no chemistry.** Out of scope for the two-week plan (PLAN §1.2).

# Day 3 — Summary

## What was implemented

A small numerical-sanity utility: a **CFL-number calculator** and its inverse, `dt_for_cfl`, so that any MFiX run planned in later days can be pre-screened against the Courant–Friedrichs–Lewy stability bound before launch. Applied to the Day-3 spout example (inlet jet ~20 m/s, cells of 1 mm) it answers the question the step file asks: *what dt does CFL < 0.5 require?* Answer: dt ≤ **2.5 × 10⁻⁵ s** (25 µs).

Key artifacts:

- `src/cfl_check.py` — `cfl_number(u_max, dx, dt)`, `dt_for_cfl(u_max, dx, cfl_target)`, and a `__main__` reporter that prints the spout example from YAML.
- `params/numerics.yaml` — the spout CFL sanity-check inputs (`u_max_m_per_s = 20`, `dx_m = 1e-3`, `cfl_target = 0.5`) with source comments — no magic numbers in scripts (CLAUDE.md rule 2).
- `tests/test_cfl_check.py` — 14 tests pinning the CFL definition, the round-trip inverse, the spout-example hand calc, and the input-validation guards.

## How it works

1. **Definition.** The 1-D CFL number is
   `CFL = u_max · dt / dx`  [-],
   with `u_max` [m/s] the peak signal / flow speed in the cell, `dx` [m] the cell edge length, and `dt` [s] the solver time step. For explicit-in-time advection this must satisfy CFL ≤ 1 for stability; CFL ≤ 0.5 is the safety margin the step file asks about.

2. **`cfl_number(u_max, dx, dt)`** returns the CFL as defined above. It rejects non-positive inputs so accidental sign errors surface as a clear `ValueError` rather than a silently negative CFL.

3. **`dt_for_cfl(u_max, dx, cfl_target)`** inverts the definition: `dt = cfl_target · dx / u_max`. This is the largest time step that meets the target CFL for a given (u_max, dx), and is the practical form used to pick `dt` for an MFiX case.

4. **Inputs come from YAML.** `params/numerics.yaml` holds the spout example values (u_max = 20 m/s, dx = 1 mm, CFL target = 0.5) with source comments (`day-3.md`, plus `assumed` for the order-of-magnitude choices). The `__main__` reporter loads that file and prints the resulting dt so the CLI has no hard-coded physics.

5. **Verification.** Hand calc:  
   `dt_max = 0.5 · 10⁻³ / 20 = 2.5 × 10⁻⁵ s`.  
   The unit test `test_spout_example_dt_hand_calc` pins the tool to that value at relative tolerance 1e-12; a companion test (`test_dt_for_cfl_inverts_cfl_number`) confirms the two functions are exact inverses. Running `python -m src.cfl_check` prints `dt_max = 2.500e-5 s (CFL = 0.500)`, in exact agreement with the hand calc. All 14 tests pass.

## What it models (physically)

The CFL number is not a physical model — it is a numerical stability criterion that says a wave (here the fluid or jet velocity) must not cross more than one cell per time step. For the spouted-bed coater it matters because:

- **The inlet jet sets the ceiling.** Anywhere in the domain, `u_max` is dominated by the spout jet — gas enters through the 9.5 mm orifice of the 0.076 m bed (PLAN §1.3) and accelerates well above the superficial column velocity. A jet of order 20 m/s is a defensible upper-bound estimate to pre-screen dt before any run; once a real case is run, `u_max` should be updated from the solver output.
- **The near-inlet mesh sets the floor.** The smallest cells will sit near the orifice/cone transition where gradients are steep. A 1 mm target there is coarse but representative; finer meshes will demand proportionally smaller dt (linear in `dx`).
- **The 0.5 margin protects the coupled solve.** MFiX runs a coupled fluid + particle system; empirical practice is to sit at CFL ≤ 0.5 to keep the fluid advection stable and to leave headroom for the particle contact / drag sub-steps. With u_max = 20 m/s and dx = 1 mm, that pins dt ≤ 25 µs — an order-of-magnitude constraint the case-builder scripts (starting Day 5) can consult.

Downstream this feeds into two operating decisions: (i) whether a proposed mesh + dt combination is even viable before we spend wall-clock on it, and (ii) how much a mesh refinement (e.g. dropping to 0.5 mm cells near the inlet) or an operating change (e.g. hotter gas → higher jet velocity for the same mass flow) tightens the time-step budget. Nothing in this step touches physics; it only prevents wasting the physics runs that come next.

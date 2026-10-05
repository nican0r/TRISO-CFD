# Day 1 hand calculation — terminal velocity, 500 µm ZrO₂-surrogate in air at 20 °C

Reference point pinned in `tests/test_correlations.py`.

## Inputs (all SI)

- d_p     = 5.00e-4 m
- rho_p   = 6000    kg/m³
- Fluid: dry air at 293.15 K, 101 325 Pa
  - rho   = 1.199   kg/m³   (Cantera gri30, N2:0.79 / O2:0.21)
  - mu    = 1.830e-5 Pa·s
- g       = 9.80665 m/s²

## Force balance

At terminal velocity: drag = net weight (buoyancy-corrected).

    0.5 · rho · u_t² · Cd · (π/4 d_p²) = (rho_p − rho) · (π/6 d_p³) · g
    ⇒  u_t² · Cd  =  (4/3) · d_p · (rho_p − rho) · g / rho

Numerically:

    RHS = (4/3) · 5e-4 · (6000 − 1.199) · 9.80665 / 1.199
        = (4/3) · 5e-4 · 5998.80 · 9.80665 / 1.199
        ≈ 32.70   m²/s²

Drag law: Schiller–Naumann, Cd = (24/Re)(1 + 0.15 Re^0.687) for Re ≤ 1000.

## Iterate

| u [m/s] | Re = ρ u d / µ | Re^0.687 | Cd     | u² · Cd | vs RHS 32.70 |
|--------:|---------------:|---------:|-------:|--------:|-------------:|
|   5.00  |         163.8  |     35.0 |  0.918 |   22.95 |  −9.75       |
|   6.00  |         196.6  |     39.6 |  0.847 |   30.49 |  −2.21       |
|   6.45  |         211.4  |     41.6 |  0.788 |   32.77 |  +0.07       |
|   6.50  |         213.0  |     41.8 |  0.785 |   33.16 |  +0.46       |

**Hand answer: u_t ≈ 6.44 m/s** (bracket ±0.05).

The test in `tests/test_correlations.py` asserts `terminal_velocity(...) ≈ 6.44 m/s`
with absolute tolerance 0.1 m/s to cover small differences in rho/mu when air
is queried at run time from Cantera (v3.2 gri30 vs. table lookups).

## Sanity: Stokes limit relaxation time

τ_p = rho_p · d_p² / (18 · µ) = 6000 · 2.5e-7 / (18 · 1.83e-5)
    = 1.5e-3 / 3.294e-4
    ≈ 4.55 s   (only qualitative — the particle is far outside Stokes here,
                 Re ≈ 210, so τ_p is a reference scale not a true time constant.)

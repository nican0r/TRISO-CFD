# Day 1 

Create the repo skeleton from §2.2. Implement in src/gasprops.py a function gas_props(T_K, P_Pa, composition: dict) -> (rho, mu) using Cantera with gri30.yaml. Implement in src/correlations.py: reynolds, cd_schiller_naumann (with the Newton cap at Re > 1000), terminal_velocity (solve with scipy.optimize.brentq), relaxation_time. Add unit tests. Start analytics.ipynb with a table and plot of terminal velocity vs particle diameter (200–1200 µm) for densities {2500, 6000, 10800} kg/m³ in (a) air at 20 °C, (b) Ar at 20 °C, (c) Ar at 1400 °C, (d) Ar/H₂ 50/50 mol at 1400 °C. 

Done when: tests pass; the 20 °C air value matches manual calculation of Terminal velocity for a 500 µm, 6000 kg/m³ sphere in air at 20 °C; the plot is saved to results/fig_d1_terminal_velocity.png.
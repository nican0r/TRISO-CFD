# TRISO-CFD

CFD of a laboratory-scale conical spouted-bed coater of the type used for
TRISO fuel particle CVD, in MFiX. The goal, scope and cases are in
[PLAN.md](PLAN.md); the day-by-day steps are in
[implementation-steps/](implementation-steps/).

## Primary results

The primary results of this project are the two most recent c01–c09 sweeps
of the C26 benchmark (the ORNL benchtop spout-fluidized bed of Conry et al.
2026, Nucl. Eng. Des. 456:115009). They differ only in the gas inlet:

| | 2×2 inlet | Single-cell inlet |
|---|---|---|
| run suffix | `_apex_v232` | `_fullface_bore_mom_nointerp_v232` (sweep tag `fullface_nointerp_v232`) |
| gas source | mass-only point source over the 2×2 axis cells at the cone tip | point source in one axis cell at the foot of a one-cell-wide (2.1 mm) square bore, momentum `ps_v_g` = Q/dx² |
| geometry | closed-apex cone (`src/make_stl_c26_apex.py`) | closed-apex cone with the bore cut into its tip (`src/make_stl_c26_fullface.py`) |
| mesh | 24 × 104 × 24, uniform 2.1 mm (axis on a cell face) | 25 × 104 × 25, uniform 2.1 mm (axis at a cell centre) |
| interpolation | LINEAR_HAT | NONE (MFiX default) |
| finished | 2026-10-27 | 2026-10-29 |
| results note | [Control-configuration sweep results](notes/c26_conry2026_case.md#control-configuration-sweep-results-2026-10-27) | [Full-face 3 s sweep results](notes/c26_conry2026_case.md#full-face-3-s-sweep-results-mfix-default-interpolation-2026-10-29) |
| sweep comparison | [c26_compare_apex_v232.csv](results/c26_compare_apex_v232.csv) | [c26_compare_fullface_nointerp_v232.csv](results/c26_compare_fullface_nointerp_v232.csv) |
| sweep summary | [c26_sweep_summary_apex_v232.csv](results/c26_sweep_summary_apex_v232.csv) | [c26_sweep_summary_fullface_nointerp_v232.csv](results/c26_sweep_summary_fullface_nointerp_v232.csv) |
| sweep figures | [flow](results/fig_c26_sweep_flow_apex_v232.png) · [diameter](results/fig_c26_sweep_diameter_apex_v232.png) · [charge](results/fig_c26_sweep_charge_apex_v232.png) | [flow](results/fig_c26_sweep_flow_fullface_nointerp_v232.png) · [diameter](results/fig_c26_sweep_diameter_fullface_nointerp_v232.png) · [charge](results/fig_c26_sweep_charge_fullface_nointerp_v232.png) |
| PIVlab re-run (c01–c03) | [c26_pivlab_rerun.csv](results/c26_pivlab_rerun.csv) · [px50](results/c26_pivlab_rerun_px50.csv) · [maps](results/fig_c26_pivlab_maps.png) · [wall](results/fig_c26_pivlab_wall.png) · [time mean](results/fig_c26_pivlab_timemean.png); method: [notes/c26_pivlab_rerun.md](notes/c26_pivlab_rerun.md) | PIVlab scores are in the sweep comparison CSV |

Both sweeps use MFiX 23.2 CFD-CGDEM, argon at 298 K, YSZ (measured PSD),
statistical weight W = 10, TGS drag, 3 s of simulated time per case. Maps,
wall velocities and sweep quantities are time averages over t = 2.4–3.0 s,
after the startup transient.

**Status against C26 (no case is validated):**

- 2×2 inlet: c02 matches C26's Table 3 floor on the in-repo PIV (CC / SSI /
  MSE 0.914 / 0.883 / 0.0015). On PIVlab with C26's image settings it scores
  0.778 / 0.825 / 0.0048, below C26's CC floor of 0.817. c01 and c03 do not
  match. The wall velocity is 2.9–6.3× C26's model.
- Single-cell inlet: c01–c03 do not match (MSE fails on all three; c02 PIVlab
  0.829 / 0.732 / 0.0075). The maps have C26's shape but run 1.3–1.6× its
  model's speed. The wall velocity is 0.64–1.06× C26's model, the first
  layout to reach that level.

### Per-case data

Each case has a PIV table (CSV), a velocity-map grid (CSV) and a figure. For
the single-cell inlet, c03 and c08 are scored on the 90th-percentile map
bottom (`_p90`; see "Map bottom for the full-face runs" in the results note).
Case definitions are in [params/sweep_C26.yaml](params/sweep_C26.yaml).

| Case | Condition | 2×2 inlet | Single-cell inlet |
|---|---|---|---|
| c01 | YSZ, 10 slpm, 60 g (C26 case 1) | [table](results/c26_piv_c01_apex_v232.csv) · [map](results/c26_piv_map_c01_apex_v232.csv) · [figure](results/fig_c26_piv_c01_apex_v232.png) | [table](results/c26_piv_c01_fullface_bore_mom_nointerp_v232.csv) · [map](results/c26_piv_map_c01_fullface_bore_mom_nointerp_v232.csv) · [figure](results/fig_c26_piv_c01_fullface_bore_mom_nointerp_v232.png) |
| c02 | YSZ, 12 slpm, 60 g (C26 case 2; sweep baseline) | [table](results/c26_piv_c02_apex_v232.csv) · [map](results/c26_piv_map_c02_apex_v232.csv) · [figure](results/fig_c26_piv_c02_apex_v232.png) | [table](results/c26_piv_c02_fullface_bore_mom_nointerp_v232.csv) · [map](results/c26_piv_map_c02_fullface_bore_mom_nointerp_v232.csv) · [figure](results/fig_c26_piv_c02_fullface_bore_mom_nointerp_v232.png) |
| c03 | YSZ, 14 slpm, 60 g (C26 case 3) | [table](results/c26_piv_c03_apex_v232.csv) · [map](results/c26_piv_map_c03_apex_v232.csv) · [figure](results/fig_c26_piv_c03_apex_v232.png) | p90: [table](results/c26_piv_c03_fullface_bore_mom_nointerp_v232_p90.csv) · [map](results/c26_piv_map_c03_fullface_bore_mom_nointerp_v232_p90.csv) · [figure](results/fig_c26_piv_c03_fullface_bore_mom_nointerp_v232_p90.png); p95: [table](results/c26_piv_c03_fullface_bore_mom_nointerp_v232.csv) · [map](results/c26_piv_map_c03_fullface_bore_mom_nointerp_v232.csv) · [figure](results/fig_c26_piv_c03_fullface_bore_mom_nointerp_v232.png) |
| c04 | + buffer, 12 slpm, 60 g | [table](results/c26_piv_c04_apex_v232.csv) · [map](results/c26_piv_map_c04_apex_v232.csv) · [figure](results/fig_c26_piv_c04_apex_v232.png) | [table](results/c26_piv_c04_fullface_bore_mom_nointerp_v232.csv) · [map](results/c26_piv_map_c04_fullface_bore_mom_nointerp_v232.csv) · [figure](results/fig_c26_piv_c04_fullface_bore_mom_nointerp_v232.png) |
| c05 | + IPyC, 12 slpm, 60 g | [table](results/c26_piv_c05_apex_v232.csv) · [map](results/c26_piv_map_c05_apex_v232.csv) · [figure](results/fig_c26_piv_c05_apex_v232.png) | [table](results/c26_piv_c05_fullface_bore_mom_nointerp_v232.csv) · [map](results/c26_piv_map_c05_fullface_bore_mom_nointerp_v232.csv) · [figure](results/fig_c26_piv_c05_fullface_bore_mom_nointerp_v232.png) |
| c06 | + SiC, 12 slpm, 60 g | [table](results/c26_piv_c06_apex_v232.csv) · [map](results/c26_piv_map_c06_apex_v232.csv) · [figure](results/fig_c26_piv_c06_apex_v232.png) | [table](results/c26_piv_c06_fullface_bore_mom_nointerp_v232.csv) · [map](results/c26_piv_map_c06_fullface_bore_mom_nointerp_v232.csv) · [figure](results/fig_c26_piv_c06_fullface_bore_mom_nointerp_v232.png) |
| c07 | + OPyC, 12 slpm, 60 g | [table](results/c26_piv_c07_apex_v232.csv) · [map](results/c26_piv_map_c07_apex_v232.csv) · [figure](results/fig_c26_piv_c07_apex_v232.png) | [table](results/c26_piv_c07_fullface_bore_mom_nointerp_v232.csv) · [map](results/c26_piv_map_c07_fullface_bore_mom_nointerp_v232.csv) · [figure](results/fig_c26_piv_c07_fullface_bore_mom_nointerp_v232.png) |
| c08 | YSZ, 12 slpm, 40 g | [table](results/c26_piv_c08_apex_v232.csv) · [map](results/c26_piv_map_c08_apex_v232.csv) · [figure](results/fig_c26_piv_c08_apex_v232.png) | p90: [table](results/c26_piv_c08_fullface_bore_mom_nointerp_v232_p90.csv) · [map](results/c26_piv_map_c08_fullface_bore_mom_nointerp_v232_p90.csv) · [figure](results/fig_c26_piv_c08_fullface_bore_mom_nointerp_v232_p90.png); p95: [table](results/c26_piv_c08_fullface_bore_mom_nointerp_v232.csv) · [map](results/c26_piv_map_c08_fullface_bore_mom_nointerp_v232.csv) · [figure](results/fig_c26_piv_c08_fullface_bore_mom_nointerp_v232.png) |
| c09 | YSZ, 12 slpm, 80 g | [table](results/c26_piv_c09_apex_v232.csv) · [map](results/c26_piv_map_c09_apex_v232.csv) · [figure](results/fig_c26_piv_c09_apex_v232.png) | [table](results/c26_piv_c09_fullface_bore_mom_nointerp_v232.csv) · [map](results/c26_piv_map_c09_fullface_bore_mom_nointerp_v232.csv) · [figure](results/fig_c26_piv_c09_fullface_bore_mom_nointerp_v232.png) |

c04–c09 are this project's cases, not C26 Table 2 cases 4–9. C26 has no data
for them, so they are compared only between the two inlets.

Input decks: `cases/C26_3d/c0X_apex_v232/` and
`cases/C26_3d/c0X_fullface_bore_mom_nointerp_v232/`. Every run, with its
parameters, mesh, dt, wall-clock time and status, is in
[results/run_log.csv](results/run_log.csv).

### Other data

Results from every other run (earlier C26 source layouts, smokes, the
interpolation A/B, and the Day 1–9 fluidized-bed and Case V figures) are no
longer tracked in git; `.gitignore` keeps only the primary results above.
Those files remain in local working copies and in git history. Raw MFiX
output is never committed; the AWS runs are on
`s3://triso-cfd-mfix/results/<run>/`.

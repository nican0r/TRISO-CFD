# Day 13

Task: Case D parametric study on the ORNL/UTK geometry (same column, cone, orifice as Case V).

Sweep axes (each one a separate run, nominal Day-10 PIC 3D case as the baseline since the TFM path failed; PIC baseline is laptop-feasible for ~16 parametric runs at ~30-60 min each):
- Flow rate: U/U_ms ∈ {1.1, 1.3, 1.5, 1.7, 1.9}. S2's reported stable range is 1.0–1.9.
- Gas temperature: air at ambient (Case V reference), Ar at ambient, Ar at coating T (~1573 K), Ar/H₂ 50/50 at coating T. Gas properties from `src/gasprops.py` / Cantera gri30.
- Particle density: ZrO₂ 6050 kg/m³ (reference) and UO₂ 10 800 kg/m³ (real-kernel sensitivity).
- Drag model: SYAM_OBRIEN (reference, isothermal) vs GIDASPOW blended drag (S1 recommends GIDASPOW for non-isothermal runs).
- Restitution coefficient: for the PIC baseline, map `c_e` to `mppic_coeff_en1` (parcel-parcel damping). Sweep ∈ {0.80, 0.85, 0.90}. (If one DEM cross-check run per axis is desired, map to `des_en_input(1)` with the same values.) c_e is assumed in Case V (not reported in S1 or S2), so it is a required sensitivity axis.

Deliverables:
- Per-axis fan plot of ΔP, PSD peak frequency, and spout-core centerline v_g.
- `results/fig_d13_D_sweep_*.png` (one per axis).
- Short notebook section in `analytics.ipynb` summarizing which axes matter most.
- `notes/day13_summary.md`.

Done when: each sweep axis has at least 3 runs completed and reported; summary identifies the two or three most-sensitive axes with numeric evidence; all runs logged.

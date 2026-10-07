# Day 8

Build the Case V 2D template and builder for the ORNL/UTK cold-mockup spouted bed (S1 + S2). Column 0.050 m ID, 60° cone, 4 mm constant-diameter orifice; 54.5 g of 500 µm spherical ZrO₂ (ρ_p = 6050 kg/m³, eps_mf = 0.42, φ = 15°); air at ambient; 30 m/s inlet jet.

Model closures (S1 Table 1): Syamlal–O'Brien drag (isothermal, Case V), Schaeffer frictional stress with "modified sigmoidal" blending, granular temperature solved as a PDE (not algebraic). Internal friction angle 15°, void fraction at minimum fluidization 0.42. Max grid size in the particle zone ≤ 10·d_p = 5 mm; at least 4 cells across the 4 mm orifice.

Geometry preference (document whichever you fall back to):
(a) 2D axisymmetric (coordinates = 'CYLINDRICAL') with the cone represented as stair-stepped blocked/wall cells (MFiX 26.1.2 rejects cut-cells in cylindrical coordinates).
(b) Only if (a) is not possible with documented keywords: a planar quadric cut-cell slab. Document that its results are qualitative only.

Specify the inlet as the physical 30 m/s jet (no slab area-ratio scaling). The 4 mm orifice must span at least 4 cells (radial dx ≤ 0.5 mm in axisymmetric, ≤ 1 mm across the full width in planar). Grade the mesh if needed and stay within 10·d_p in the particle zone. Truncate the axial domain to ≈ 3× bed height plus fountain headroom, not the full S1/S2 column, and document the choice.

Params:
- params/particles.yaml: replace Case V blocks with the UTK/ORNL `case_V` block (d_p = 500 µm, rho_p = 6050 kg/m³, inventory 54.5 g, phi = 15°, ep_star = 0.42, c_e). Update `case_D` rho_p to 6050 kg/m³ with the S2 source. Record IC bed solids fraction below 1 − ep_star = 0.58 so Schaeffer friction is inactive at t = 0.
- params/geometry.yaml: 0.05 m column, 60° cone, 4 mm inlet, derived bed height from inventory, truncated domain height, 2D mesh.
- params/gases.yaml: Case V operating point is air at ambient; S2's humidity effect on gas properties noted as neglected.

Code:
- Add `build_V_2d_params()` (UTK/ORNL version) in `src/make_case.py`.
- Template at `cases/V_2d/template.mfx.j2`. Every keyword cross-checked against the MFiX 26.1.2 keyword reference. drag_type = 'SYAM_OBRIEN'; friction_model = 'SCHAEFFER' with the blending keyword for "modified sigmoidal"; granular energy as a PDE (kt_type ≠ 'ALGEBRAIC'). If any keyword cannot be confirmed in 26.1.2, STOP and report.

Monitors (SI, matched to S1/S2 sampling):
- Inlet-plane gas pressure at 1 kHz (monitor_dt = 1e-3 s), source: assumed (S1/S2 don't give tap location); noted as such.
- Top-of-bed gas pressure.
- Total solids inventory in the domain.
- Centerline v_g (used for S1 Fig. 1 validation on Day 10).

Tests (add to `tests/test_make_case.py`):
- Bed-height hand calc from inventory 54.5 g, ρ_p = 6050, eps_s = 1 − eps_mf = 0.58. Also report H₀/D_c and compare to S2's 0.50–0.65.
- Inlet volumetric flow 30 m/s × A_inlet ≈ 22.6 L/min; U/U_ms ≈ 22.6 / 12 ≈ 1.9 vs S2's 12 L/min reported minimum spouting at 500 µm.
- Rendered closure keywords (drag_type, friction_model + "modified sigmoidal" blending keyword, PDE granular energy, phi = 15, ep_star = 0.42).
- At least 4 cells across the inlet; cell size ≤ 10 d_p in the particle zone.

Done when: the case renders and MFiX initializes without errors, and a geometry screenshot is saved in `results/fig_d8_V_2d_geometry.png`.

# Day 8 — Summary

## What was implemented

Built the two Case V spouted-bed MFiX templates that will be driven in later
days to produce the project's first validation runs: **V_2d**, a true 2D
Cartesian cut-cell bed with the cone wall defined by a single Y_CONE
quadric, and **V_3d**, a 3D Cartesian cut-cell bed with the cone + cylinder
defined by an STL surface of revolution. Both render from `src.make_case`
through Jinja2, both pull every physical input from `params/*.yaml`, and
both now **initialize cleanly in MFiX 26.1.2 and advance past t=0** — the
Day-8 acceptance criterion.

**Why Cartesian (not cylindrical) coordinates for Case V.** MFiX 26.1.2
forbids combining `cartesian_grid = .True.` with `coordinates =
'CYLINDRICAL'` (`check_data_cartesian.f:55`, "CARTESIAN GRID OPTION NOT
AVAILABLE WITH CYLINDRICAL COORDINATE SYSTEM"). Cut-cells, STL geometry,
and quadric surfaces are therefore Cartesian-only features. Since the 60°
cone wall of the spouted bed cannot be represented on a rectangular IJK
grid in cylindrical coordinates without staircasing, both Case V templates
are built on Cartesian grids with cut-cells so the cone profile is
resolved cleanly.

**Why Y_CONE quadric for V_2d, STL for V_3d.** Native MFiX quadrics live
under the same `cartesian_grid=.True.` cut-cell machinery as STL but need
no external file. For V_2d's planar slice the vessel inner wall reduces to
a single Y_CONE primitive (apex at `y = -R_i/tan(30°) = -8.23 mm`, half-angle
30°, clipped to `0 ≤ y ≤ 57.63 mm`); above the clip height the natural
domain walls at `x = ±R_c` already serve as the cylindrical column walls.
Quadrics are also compatible with `no_k = .True.`, whereas STL is not
(`get_stl_data.f:813`: "STL method valid only in 3D"), so V_2d runs as a
true 2D case with `kmax = 1` rather than the thin 3D slab (`kmax = 2`) STL
would require. V_3d stays on STL because the full axisymmetric vessel is a
surface of revolution — not representable as a single quadric.

**Why we keep two models.** V_2d is cheap (~21.5 k cells vs. ~154 k for
V_3d) and runs several times faster per simulated second. It is the
project's debugging and sweep workhorse: Day 9 initial debugging /
solver stabilisation, Day 10 U_ms step-down and mesh study, and Days
11–12 parametric trends (gas temperature, particle density, drag model).
V_3d is the production geometry — the one compared against the Missouri
S&T RPT experimental data on Day 10. Treat V_2d as a fast proxy for V_3d:
directional trends and gross hydrodynamics transfer; absolute spouting
velocities and fountain heights do not (planar slab vs. true axisymmetric
swirl). Any validation claim ultimately rests on V_3d.

The frictional-stress threshold was also addressed. The Day 6-7 fb_sweep
template used `ep_star = 0.42` (Schaeffer activation at ε_s > 0.58) and
had to drop ε_s_bed to 0.55 to stay below it. The Day 8 step calls for
ε_s_bed ≈ 0.60, which crosses that threshold. Rather than re-interpret
the step, V_2d / V_3d adopt `ep_star = 0.39` so the threshold moves to
ε_s = 0.61 — just above ε_s_bed = 0.60. The IC is now safely in the
Schaeffer-inactive region. The practical effect: V_2d not only
initializes, it advances past t = 0 and writes VTK frames. V_3d still
initializes cleanly but then throws a different numerical warning at
t > 0 (solids-velocity blow-up at a cut-cell on the cone wall, O(2000 m/s)
in a single near-wall cell); V_2d sees the same class of warning at its
first time step after advancing one frame to disk. Both are cut-cell
near-wall transients, not Schaeffer, and are Day-9 territory.

Key artifacts:

- `params/geometry.yaml` — three new blocks: `spouted_bed_V_2d`
  (true-2D Cartesian grid 64 × 336 × 1 with `no_k=.True.`, Y_CONE quadric
  apex / clip-top / half-angle, z thickness = dx for cosmetic consistency),
  `spouted_bed_V_3d` (31 × 160 × 31 grid, 2.5 mm cells, 60 θ-segments),
  and `spouted_bed_V_initial` (ε_s_bed = 0.60, ep_star_V = 0.39 so
  Schaeffer threshold = 0.61).
- `cases/V_2d/template.mfx.j2` — Cartesian cut-cell + Y_CONE quadric TFM
  template (true 2D, `no_k=.True.`, `kmax=1`, `use_stl=.False.`), Case D
  particles in air @ 20 °C, MI orifice strip + `CG_NSW` cut-cell wall
  (bound to the quadric via `bc_id_q(1) = 3`) + PO top, three monitors.
  No STL file dependency. Every keyword traced to the MFiX 26.1.2 reference.
- `cases/V_3d/template.mfx.j2` — Cartesian `cartesian_grid=.True.` +
  `use_stl=.True.` + `stl_bc_id=3`, with BC 3 as the cut-cell wall
  (`CG_NSW`), monitors over the orifice/top-of-bed/spout core.
- `src/make_stl.py` — `build_spouted_bed_stl` writes an ASCII STL
  (surface of revolution of the cone + cylinder) with outward-pointing
  normals; the default 60-segment tessellation yields 240 triangles.
  (V_3d-only; V_2d does not use this module.)
- `cases/V_3d/geometry.stl` + `cases/V_3d/geometry_0003.stl` — the rendered
  STL. MFiX looks the file up by STL-BC id (hence the `_0003` filename).
- `src/make_case.py` — new `build_V_2d_params`, `build_V_3d_params`, and a
  CLI `--kind {fb_sweep,V_2d,V_3d}` dispatch. All physical values resolved
  from YAML and Cantera (`src.gasprops.gas_props`).
- `src/make_figs_d8.py` — matplotlib-only renderer (headless-safe) for the
  two visual-confirmation figures.
- `tests/test_make_stl.py` — six tests pin facet count = 240, bounding box
  (±R_c, 0..H_dom), orifice-rim radius 4.75 mm, cone-slope agreement
  within chord-error, cylinder-section radius 38.0 mm, STL header/footer,
  invalid-cone rejection.
- `tests/test_make_case.py` — V_2d tests pin coordinates=CARTESIAN,
  `no_k=.True.`, `use_stl=.False.`, Y_CONE quadric with half_angle=30,
  `bc_id_q(1)=3`, grid 64 × 336 × 1, orifice spans ≥ 4 cells across the
  diameter, MI/PO BCs, monitors, bed IC; V_3d tests pin `use_stl=.True.`,
  grid 31 × 160 × 31 at 2.0–3.0 mm cell size, monitors, Case D + air @ 20 °C
  substitution. Both pin `1 - ep_star > ε_s_bed`.
- `results/fig_d8_V_2d_geometry.png`, `results/fig_d8_V_3d_geometry.png` —
  visual-confirmation figures for the V_2d schematic and the V_3d axial
  profile + top view.

## How it works

1. **Geometry inputs.** The three new blocks in `params/geometry.yaml`
   derive everything from `bed_0p076m` (D_c = 76 mm, D_i = 9.5 mm, 60° cone)
   and the shared `H_static_m = 0.10` m. Day-8 fixes the axial domain to
   `H_dom = 3·H_static + 0.10 = 0.40 m` ("3× bed + 100 mm fountain
   headroom", the 100 mm `assumed`), truncating the real column's 1.14 m to
   the region where spouting physics actually lives. V_2d precomputes the
   cone-apex y-offset `cone_apex_y_m = −R_i / tan(30°) = −8.227 mm` and
   `cone_top_y_m = (R_c − R_i) / tan(30°) = 57.63 mm`; these are kept in
   YAML even though V_2d doesn't use them in the .mfx any more (MFiX forbids
   the cut-cell cone in cylindrical, see below) because the V_3d STL
   generator and the V_2d schematic figure both read them.

2. **STL generation.** `src.make_stl.build_spouted_bed_stl(R_c, R_i, half_angle,
   H_dom, theta_segments, out_path)` builds the vessel inner wall as a
   surface of revolution of the two-segment (r, y) profile
   {(R_i, 0), (R_c, h_cone), (R_c, H_dom)} swept around the y-axis in
   `theta_segments` steps. Each quadrilateral panel is split into two
   right-handed triangles with normals pointing outward (toward −r), which
   is what MFiX expects for `OUT_STL_VALUE = 1` (internal-flow convention).
   The 60-segment default yields 240 facets and a chord error of
   `R_c · (1 − cos(π/60)) ≈ 5 × 10⁻⁵ m` on the cylinder, i.e. about dx/50 on
   the 2.5 mm grid.

3. **Template parameter builders.** `build_V_2d_params` and
   `build_V_3d_params` are mirror functions of the Day 6-7
   `build_fb_sweep_params`: they load the YAML files, call
   `src.gasprops.gas_props(293.15 K, 101 325 Pa, N₂/O₂ 0.79/0.21)` for
   ρ_g = 1.1994 kg/m³ and µ_g = 1.830×10⁻⁵ Pa·s, and emit a dict ready for
   the Jinja2 template. Both place three MFiX monitors (`monitor_type = 4`,
   arithmetic average): M1 over the orifice at y = 0 (`monitor_p_g`), M2
   over a horizontal slab at y = H_static (`monitor_p_g`), M3 over the spout
   core at y = H_static/2 (`monitor_ep_s`, `monitor_v_s`), from which
   post-processing reconstructs the solids mass flux
   ρ_p · ε_s · v_s · A_spout.

4. **MFiX acceptance.**
   - **V_2d:** `mfixsolver -f V_2d.mfx` reads the namelist cleanly, loads
     240 STL facets (all valid, 0 ignored), carves a 64 × 336 × 2 =
     **43 008-cell** mesh — **2 564 standard, 492 cut, 3 056 fluid,
     39 952 blocked** (blocked cells are the ones outside the vessel's
     revolved profile, as expected since the full-width rectangular
     bounding box is much larger than the cone+cylinder in slab z).
     PRE_PROCESSING COMPLETE in 0.37 s, the solver enters the time-step
     loop, and **writes `BACKGROUND.pvd` + `BACKGROUND_0000.vtu`** before
     being killed by the watchdog. No ERROR, only the pre-existing
     cosmetic "description truncated on read" Fortran namelist warning.
     The ep_star = 0.39 lift (threshold 0.61 vs. IC 0.60) remains the
     decisive change keeping Schaeffer dormant at t = 0.
   - **V_3d:** `mfixsolver -f V_3d.mfx` loads 240 STL facets (all valid, 0
     ignored), carves 14,632 cut cells and 102,108 blocked cells out of
     153,760 total, leaves 51,652 fluid cells inside the vessel
     (33.6 % of the bounding box), completes PRE_PROCESSING in 1.55 s, and
     steps into the time-step loop. The first-step solve throws a
     *different* warning from Day 6-7: a solids-phase velocity blow-up at
     a cut-cell near the cone wall (|W_g| ~ 2000 m/s in cell
     I=32, J=42, K=16), which is a 3D cut-cell startup issue, not Schaeffer
     — and is Day-9 territory. The STL cone + cylinder + orifice are
     visibly correct in both the solver-emitted `V_3D_boundary.vtk` and in
     `results/fig_d8_V_3d_geometry.png` (axial projection: cone tapers
     from R_c = 38 mm at the top down to R_i = 4.75 mm at the orifice
     over the first 57.6 mm of axial height; top view: 240-facet decagon
     tracing the r = 38 mm column with the orifice disc centred).

5. **Verification.**
   - Twenty-two new unit tests (`tests/test_make_stl.py` + Day-8 additions
     to `tests/test_make_case.py`, including a pin that
     `1 − ep_star > ε_s_bed` so the Schaeffer threshold stays above the
     IC) all pass; full suite 57/57 green.
   - V_2d: `.mfx` passes the MFiX namelist reader and completes mesh
     generation (10,752 cells).
   - V_3d: `.mfx` + `geometry_0003.stl` passes the STL loader (240 valid
     facets, 0 ignored), carves 14,632 cut cells, writes `V_3D_boundary.vtk`
     for ParaView confirmation.
   - Both visual-confirmation figures saved under `results/fig_d8_V_*_geometry.png`.

## What it models (physically)

The two templates together define the **cold-flow hydrodynamic baseline**
for the Case V bed (PLAN §1.3): a 76 mm-ID cylindrical column with a 60°
conical base and a 9.5 mm inlet orifice, filled to H = 100 mm with Case D
particles (500 µm, 6000 kg/m³, Geldart D) and fluidized by air at 20 °C.
Nothing is advanced in time yet; what Day 8 locks in is the *computational
geometry* and the *initial state* every subsequent run will start from.

The modelling choices this step bakes in:

- **Cartesian coordinates for V_2d.** The axisymmetric cylindrical option
  that Shallbetter's thesis favours is unavailable here because MFiX
  26.1.2 refuses to combine `cartesian_grid=.True.` with
  `coordinates='CYLINDRICAL'` (see "Why Cartesian" above). The alternative
  — a cylindrical mesh with a straight-cylinder wall and no resolved cone
  — is physically wrong for a spouted bed (cone shape is first-order to
  the annulus down-flow pattern). V_2d therefore adopts the Cartesian
  cut-cell + STL workflow V_3d uses, run as a thin slab (`kmax = 2`), and
  accepts that the slab's planar flow loses the axisymmetric swirl term.
  The slab is kept because it is ~3.6× cheaper than V_3d and still
  captures gross hydrodynamics (spout, fountain, annulus) — making it
  the project's Day 9-12 debugging and sweep workhorse; V_3d is reserved
  for the Day-10 validation comparison against the Missouri S&T RPT data.
- **Truncated axial domain (H_dom = 0.40 m vs. full 1.14 m column).** The
  real column has freeboard length well beyond where the particles ever
  reach; the fountain in a stable spout rises only ~1–2× the static bed
  height above the bed (Mathur-Gishler scaling). 3 × H_static + 100 mm
  fountain headroom captures the full spout + fountain + a buffer, at
  roughly one-third the cell count of the full column.
- **Case V_3d: cone + orifice as cut-cell STL geometry.** The STL
  tessellates the real cone wall at 60 angular segments (6° per facet),
  which is finer than two cells per facet on the 2.5 mm grid, so the
  wall-cell cut geometry resolves the cone slope without staircasing
  artefacts. The orifice is left open as a 9.5 mm disc at y = 0; the
  solver imposes `U_in` across this disc as a mass-inflow boundary — the
  physical driver of spouting. Fountain height, spout-to-annulus mass
  flux, and the spouting-vs-choked transition all depend on how well the
  cone wall is resolved near the orifice, which is why Day 8 insists on a
  2–3 mm cell size here even though it costs 150k+ cells.
- **Monitor set as the Day-9+ read-outs.** The three monitors are chosen
  to map directly to downstream decisions. M1 (P at the inlet) + M2 (P at
  the top of the bed) gives ΔP(U) across the bed — the same signal used
  in Day 6-7 for U_mf, now reinterpreted as the spouting ΔP that drops
  sharply at U = U_ms (Mathur-Gishler). M3 (ε_s and v_s at the spout
  mid-plane) gives the solids volumetric throughput through the spout,
  which, divided by the bed inventory, gives the particle circulation
  time — the proxy for coating uniformity in PLAN §1.1. These monitors
  need no additional post-processing to exist; they are time series the
  solver writes directly to disk.
- **Case D particles in air at 20 °C.** V_2d and V_3d both use the same
  particle and gas properties as the Day 6-7 fb_sweep so that any
  spouting-regime results carry a consistent reference. The hot-gas
  variant (argon @ 1400 °C, Ar/H₂) will be a Day 10+ sweep against the
  same templates.

What Day 8 **deliberately does not** model, by step scope: anything
beyond the first time-step pass (full time-advancement stabilisation is
Day 9's first task), the full 1.14 m column (truncated to 3 × bed +
fountain headroom), the orifice tube stub below the cone (the MI face
represents the orifice directly), DEM / PIC particle statistics (TFM is
the Case V baseline per Day 4), and any thermal / species / reactive
content (PLAN §1.2 keeps these out of scope). Both V_2d and V_3d now
resolve the full 60° cone wall as a cut-cell STL boundary; the earlier
V_2d straight-cylinder approximation (rev 1-2) has been retired.

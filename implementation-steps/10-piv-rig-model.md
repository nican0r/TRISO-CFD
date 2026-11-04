# 10. PIV rig model: changes from the C26 model

Branch `piv-rig-model`. This adapts the C26 model to the benchtop rig built
to replicate the C26 experiment, using the measurements in `data/piv_rig/`
(IDs such as T1 and N1 below are rows of its Apparatus sheet, measured
2026-11-02).

**Status: built and dry-run only (1 ms). No smoke or production run exists,
so nothing here is validated against the experiment yet.**

The rig is added alongside C26, not in place of it. C26 decks and STLs still
render byte-identical (checked for three decks, including the production
single-cell deck's STL, md5 `59ddf7c1`).

```bash
.venv/bin/python -m src.run_C26_3d --rig piv_rig --case all
# writes cases/PIV_rig_3d/r01, r02, r03; decks RIG_r0X.mfx
```

## Values that changed

| Quantity | C26 model | PIV rig | Source |
|---|---|---|---|
| Column ID | 50.4 mm | **51.0 mm** | T1 |
| Orifice / throat diameter | 2.38 mm | **2.70 mm** | N1 |
| Cone angle | 60° | **60° exact** (sheet: 2 × 29.8°) | your instruction |
| Cone top | r = 25.2 mm at y = 41.6 mm | **r = 27.5 mm at y = 45.29 mm** | F8, 60° |
| Rim pocket and step | none | **vertical wall to y = 50.0 mm, then a 2 mm downward step to the tube** | F8, F11, T1; your choice |
| Outlet height | 210 mm | **354.57 mm** (tube top 354.80 mm, nearest whole row) | F11 + T3 |
| Cells across column | 24 (2.1 mm) | **21 (2.4286 mm)** | your choice |
| Box | 25 × 104 × 25 | **23 × 150 × 23** (50,799 fluid cells; C26 42,520) | derived |
| Source face vs orifice | 4.41 vs 4.45 mm² (−1 %) | **5.90 vs 5.73 mm² (+3 %)** | derived |
| Gas | Ar, 25 °C, 1.633 kg/m³ | **N₂, 25 °C, 1.145 kg/m³, 1.80e-5 Pa·s** | your instruction |
| slpm reference | 25 °C, 1 atm, Ar | **25 °C, 1 atm, N₂** | your choice |
| Particle | YSZ, 411.6 µm, σ 12.2 µm | **steel, 822 µm, σ 10.3 µm** (±3σ) | Particles sheet |
| Density | 6000 kg/m³ | **7585 kg/m³** | Particles sheet |
| Statistical weight W | 10 | **1** (C26 Table 2 steel) | your choice |
| Parcels | 27,317 estimated | **27,188 estimated**; 27,031 seeded (59.68 g) | 27,200 × 2.2059 mg = 60 g |
| Packed fraction (H_static only) | 0.57 | **0.60**, giving H_static = 31.2 mm | implied by the sheet's 13,184 mm³ |

Mesh note: at a 2.7 mm orifice, 22 or more cells do not fit the bore because
the cone's apex (y = −2.338 mm) falls below the bore floor. So 24 cells could
not carry over; 21 keeps one cell face about equal to the orifice.

## Cases

The inlet layout is unchanged: one walled source cell on the axis in the bore
row y ∈ [−dx, 0], with momentum, and gas leaving only through the top face.

| Case | Flow | ṁ | Orifice velocity | `ps_v_g` = Q/dx² |
|---|---|---|---|---|
| r01 | 24 slpm | 4.580e-4 kg/s | 69.9 m/s | 67.8 m/s |
| r02 | 26.5 slpm | 5.057e-4 kg/s | 77.1 m/s | 74.9 m/s |
| r03 | 30 slpm | 5.725e-4 kg/s | 87.3 m/s | 84.8 m/s |

All cases use 60 g.

## Carried over unchanged (not measured on the rig)

- C26 Table 1 contact values, whose steel column equals its YSZ column:
  kn 1000 N/m, kt factor 0.285, restitution 0.5, friction 0.1, damping
  factor 0.5. Steel-on-glass and steel-on-acrylic are not distinguished.
- TGS drag, NONE interpolation, gas-fraction floor 0.42.
- dt 1e-5 s (the CFL cap stays above it at all three flows), MFiX 23.2,
  3 s run, 300 frames at 2 ms, 4 ranks.
- Assumptions: pure N₂, 1 atm, normal PSD truncated at ±3σ, open outlet at
  the tube top, W = 1 run through MFiX's CGP path (identical to DEM physics
  at W = 1, and the same code path as the production decks).

## Measured but not modeled

- **F1, F12, N2:** the 0.65 mm plug-top ledge and the 1.5 mm straight glass
  bore. Both are well under dx. One cone runs from the orifice, as in the
  Particles sheet's own bed prediction. Above the ledge the real wall sits
  0.21 mm from the modeled line.
- **F4–F7:** superseded by 60° exact.
- **F9, F10:** lip outer diameter and lip height. The vertical rim band
  stands in for the lip.
- **F2, F3, N3, T2, T4:** outside surfaces and fit.

## Code changes

| File | Change |
|---|---|
| `params/geometry.yaml` | new `bed_piv_rig`, `spouted_bed_piv_rig_3d`, `spouted_bed_piv_rig_3d_dem` |
| `params/particles.yaml`, `gases.yaml` | new `case_piv_rig_ss`, `n2_25C`, `slpm_reference_piv_rig`; `short` label on the gas |
| `params/sweep_piv_rig.yaml` | new: r01–r03, W floor of 1 for this rig only |
| `src/make_case_c26.py` | `RIGS` table and `rig=` argument; rim pocket; `ic_top_on_grid_row` |
| `src/make_stl_c26_fullface.py` | optional `R_rim_m`, `y_rim_top_m`: the pocket and step |
| `src/run_C26_3d.py` | `--rig`; rig-specific output directory, run names and manifest labels |
| `cases/C26_3d/template_dem.mfx.j2` | comments only, by rig; no new keywords |
| `tests/test_make_case_piv_rig.py` | new: 28 tests, including yaml against the CSVs and hand values |
| `data/piv_rig/` | copies of the two measurement sheets |

## Checks done

- Full test suite: 311 passed.
- 1 ms dry run of r01 on the 23.2 image (`results/run_log.csv`): exit 0,
  12,496 facets read and 0 ignored, source volume exactly dx³, one fluid cell
  in the source row with v_y = 67.82 m/s and no sideways velocity. Rows up to
  the step reach r = 27.2 mm; rows above it reach r = 25.2 mm. No parcels lost
  over the 1 ms.
- **A bug the dry run caught:** the first attempt seeded only 25,011 parcels
  (55.2 g). MFiX snaps the IC top to a grid row and counts parcels from the
  snapped volume, and 16.44 rows rounded down to 16. The builder now puts the
  top on a row and solves `ic_ep_s` for the charge (0.3242). The C26 decks
  have the smaller opposite offset (27,765 seeded against 27,317, +1.6 %) and
  are left as they are.

## Open items

- **Rim height does not close.** A 60° cone plus the lip gives 47.9 mm against
  the measured 50.0 mm. The model puts the 2.1 mm in the vertical band. The
  sheet asks for a re-check of F11, F12 or the cone marks above ±1 mm. F5
  (2.01 mm) also sits inside the 2.26 mm lip, where the sheet asked for 4–5 mm.
- **Density** comes from a 0.01 g balance on 0.3 g total, and lands 4.4 %
  below nominal 304 SS. Resolution alone is a few percent.
- **Gas bypass.** The model sends all gas through the orifice. The sheet's
  radial plug clearance is 0.33 mm (aim 0.1–0.3) and the stem is 27 mm (fit
  note ≥ 33), so check that the plug seals.
- **Not adapted yet:** the PIV scoring in `src/post/` targets C26's digitized
  figures, so comparing against your rig's PIV needs your data and a new
  target definition. AWS submission of these run directories has not been
  exercised.
- **Next step (CLAUDE.md rule 11):** smoke r01 for 0.5 s and check it spouts
  (`--tstop 0.5 --run-dir-suffix _smoke05`). Check the first-vs-last parcel
  count at the new rim corners too, since 23.2 can delete parcels at cut-cell
  corners. Expect roughly C26's AWS cost plus about 20 % for the extra fluid
  cells; that is an estimate, not a measurement.

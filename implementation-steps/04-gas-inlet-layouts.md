# 4. Gas inlet layouts tested

C26 says only that the gas enters through "a point source centered with
the reactor and placed at the bottom of the cone", on cells sized to the
2.38 mm nozzle. Most of the C26 work went into finding out what that means
on a real MFiX grid, because the inlet turned out to control most of the
bed behavior.

## How an MFiX point source works

A point source (`ps_*` keywords) adds a gas mass flow ṁ to a region of
cells. With `ps_v_g` set it also adds momentum ṁ·v_ps. Two facts govern
every layout below:

1. **Continuity sets the jet, not the typed velocity.** If the source
   cells are walled on every side but the top, the gas must leave through
   the top faces at Q/A_open, whatever `ps_v_g` says. The typed velocity
   then only changes the pressure inside the source cell. A 1 ms test at
   1/64 of the velocity changed the exit velocity by 0.1 %; only the
   source-cell pressure changed, by ρV².
2. **The jet momentum entering the bed is set by the first open row's
   area, ρQ²/A.** A nozzle-sized opening carries the full nozzle momentum
   ṁU. Spreading the same gas over 2 × 2 cells (4 × the area) carries at
   most ¼ of it. If the source cell's side faces are open, gas also leaves
   sideways into the cone tip.

The quantities used to compare layouts: **up share** (the fraction of the
source gas leaving through the top), **lean** (offset of the spout from the
axis), **fountain height and reach**, the **wall velocity** in the annulus,
and the map scores (file 07).

## The layouts, in order

| # | Layout (run-dir suffix) | Geometry / grid | Result | Fate |
|---|---|---|---|---|
| 1 | Stub pipe + mass inlet (Day-10 recipe) | – | at R_t/dx = 0.57 the stub has no fluid cells | never built |
| 2 | Zero-width point at the axis-cell center (first smoke) | open apex, 27 cells | the center is a rounding tie; the source landed in a corner sliver cell (3.3e-10 m³) and drove a 150–170 m/s wall-hugging jet | bug; source regions are now written as cell faces |
| 3 | One cell + momentum (`c01`) | open apex, 27 cells | spout centered but overdriven: axis parcels 1.75 m/s, fountain 55 mm; wall 1.1e-2 m/s | reference only; its open side faces let the momentum entrain gas |
| 4 | One cell, mass-only (`_massonly`, production 2026-10-19) | open apex, 27 cells | ~⅔ of the gas left sideways (up share 0.34–0.38); spout leaned −x by 0.8–1.2 mm; maps nearly static | dropped |
| 5 | 2 × 2, mass-only (`_massonly_axisface`) | open apex, 26 cells | up share 1.0, spout centered, spout checks pass; wall velocity unchanged (1.8–3.4e-2 m/s, ~5× C26's model) | fixed the spout, not the wall |
| 6 | 2 × 2 + momentum (`_momentum24`) | open apex, 24 cells | identical to #5 except the source-cell pressure (fact 1) | source type closed as a candidate |
| 7 | 2 × 2, mass-only (`_apex`, `_apex_v232`) | **closed apex**, 24 cells | tip jet ~70 m/s (the tip cells are cut as tetrahedra); the jet picks one of the four axis cells by round-off (~1 mm lean, different cell on different machines). 3 s on 23.2: c02's map came closest to C26's of the layouts tried; PIVlab CC / SSI / MSE 0.778 / 0.825 / 0.0048 (CC and MSE miss their floors); wall 3.0× C26's model | **control inlet**, c01–c09 sweep (2026-10-26) |
| 8 | One-cell point, mass-only (`_point_v232`) | closed apex, 24 cells | the point snaps to the −x/−z cell, 1.5 mm off axis; up share 0.27; faster wall layer than #7 | dropped |
| 9 | One-cell point + momentum (`_point_mom_v232`) | closed apex, 24 cells | fountain height close to C26's model, but the spout leaned off axis, the lower map ran fast, and dt was cut in 23 % of steps | chosen by the pre-registered rule, overridden by the user; its sweep was stopped |
| 10 | 2 × 2 + momentum (`_apex_mom_v232`) | closed apex, 24 cells | worse than #7 (shorter fountain, faster wall layer); not explained | dropped |
| 11 | **Full-face bore, one cell + momentum** (`_fullface_bore_mom_v232` smoke, `_fullface_bore_mom_nointerp_v232` 3 s) | bore, 25 cells | up share 1.000, lean 0.03–0.28 mm; wall at C26's model level (0.64–1.06×); maps 1.3–1.4× C26's model speed | **single-cell inlet**, c01–c09 sweep (2026-10-27) |

Results for #7 and #11 are in the [README](../README.md) and file 08.

## Why two primary inlets

**The authors' answer (personal communication, 2026-10-27):** "the point
sources were defined as mass flow with velocity vector and point. Velocity
mag is calculated based on the local mesh characteristics." This rules out
a mass-only source and points to one cell with momentum, its speed typed by
hand as Q / (face area). MFiX does not compute it: the GUI default is 0.

That leaves a split between what C26 *describes* and what reproduces C26's
*result*:

- **2 × 2 closed-apex control (#7)** came closest to C26's c02 map of the
  layouts tried when the sweep inlet was chosen (2026-10-26). On PIVlab
  even its c02 misses the CC and MSE floors. It does not match at 10 or
  14 slpm either, its wall layer is 3–6× too fast, and it is not what the
  authors describe.
- **Single-cell bore (#11)** is the closest reading of the authors' reply
  that MFiX can represent with the source cell on the axis:
  - one standard cell, walled on its sides and bottom (file 03), so the
    gas leaves only through its whole top face (dx² = 4.41 mm², nozzle
    4.45 mm²);
  - `ps_v_g` = Q/dx² = 37.8 / 45.4 / 52.9 m/s at 10 / 12 / 14 slpm.

  On a face-centered grid a "point" lands 1.5 mm off axis (#8, #9), and on
  a center-centered grid a closed apex has no source cell. The bore solves
  both.

The two primary sweeps differ in the inlet, the grid and the
interpolation scheme (LINEAR_HAT for the control, NONE for the bore;
file 05). No full-face LINEAR_HAT 3 s run exists, so the interpolation
effect on the bore is known only from a 0.3–0.5 s matched-time comparison.
There it did not change the fountain (file 08).

## Still open about the source

- Whether C26's source cell was whole or cut by the cone wall. A cut
  cell's side faces leak momentum into the cone tip (our #9 had an up share
  of 0.58). The up share of any point-source deck is a free parameter until
  the authors answer.
- What velocity they typed at 10 slpm, and where the axis sat relative to
  their grid. Questions are drafted for the authors.
- The one-cell momentum point source with NONE interpolation (#9 with
  NONE) has not been run.
- `params/geometry.yaml` defaults stay `ps_momentum: false`,
  `ps_point: false`, `stl_builder: implicit`, pending a user decision. Each
  primary deck selects its layout explicitly (file 06).

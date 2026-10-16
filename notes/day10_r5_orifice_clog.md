# Day-10 R5 — orifice "clog" in the R4 AWS run: cause note (CLAUDE.md rule 9)

Written 2026-10-16, BEFORE any code change.

## Symptom (R4 AWS job bb113f48, t = 0 – 0.521 s)

- t ≈ 30–40 ms: inlet–top ΔP spikes to +19–21 kPa (bed lifts as a plug),
  then collapses; afterwards median dt drops ~30× (≈1e-4 → 4–9e-6 s) and
  stays there; cost rises from 15 h to 76–125 h wall per simulated second.
- 52 → 82 parcels sit at **y = −0.4 mm** (all at the same height) inside the
  stub, vy ≈ −0.2 m/s but not falling.
- Orifice (y ∈ [0,1) mm, r < 2 mm) mean v_g ≈ 100 m/s vs bc_v_g = 30 m/s;
  ε_g ≈ 0.56–0.62.
- ΔP(inlet − top) mean −906 Pa (wrong sign vs W/A ≈ +401 Pa), σ ≈ 4.6 kPa;
  63 % of PSD power > 100 Hz.

## What the cell-level fields show

Gas cells at y = −0.5 mm (just below the orifice plane), every frame t ≥ 0.05 s:

- 4 central stub cells (x, z = ±0.5 mm): **v_g = 0.0 exactly**, ε_g = 0.33 (clip).
- 4 ring cut-cells (r ≈ 1.4 mm): v_g = 110 – 500 m/s, ε_g = 0.33.

Above the plane (y = 0.5 – 4.5 mm) the jet is **annular and off-axis**
(100 – 250 m/s in the r ≈ 1.5 mm cells) with **reverse flow on the axis**
(−10 to −50 m/s in the central cells). That is not a free round jet.

## Cause: a spurious STL membrane across the orifice at y = 0

`src/make_stl.py::build_spouted_bed_stl_implicit` builds the vessel as the
union (min) of three capped primitives:

    sdf_stub = max(r − R_stub, −L_stub − y, y − 0)      # y ∈ [−L_stub, 0]
    sdf_cone = max(r − (R_i + y·tanθ), 0 − y, y − h_cone)  # y ∈ [0, h_cone]

Both primitives have a flat cap at exactly y = 0. The sampling grid
`ys = −L_stub − 2·dx_s + k·dx_s` lands **exactly** on y = 0.0 (verified:
`ys[k] == 0.0` for both dx_s = 0.5 mm and 1 mm). On that plane
`sdf_stub = sdf_cone = 0` for every r < R_i, so the union SDF is identically
zero over the whole orifice disc and marching cubes emits a surface there.

Evidence, from the STL actually used by R4 (`cases/V_3d/day10/geometry_0001.stl`):
in y ∈ [−0.25, 0) mm there are **144 facets with |n_y| > 0.7 covering
r = 0.24 – 1.84 mm**, i.e. a near-complete disc that leaves only a thin ring
gap at r ≈ 1.84 – 2.0 mm. Away from y = 0, the stub has only its sidewall
facets (r = 1.98 mm), as expected.

MFiX treats that disc as a CG_NSW wall. So:

1. **All the inlet gas squeezes through the thin ring gap** at the orifice.
   The local velocity is 3 – 15× bc_v_g, which gives the off-axis annular jet,
   the on-axis recirculation and the "orifice v ≈ 100 m/s" average.
2. **The central stub cells under the disc are a dead-end pocket** (v_g = 0).
   Parcels that fall through the ring gap or recirculate down the axis stall
   in that pocket and accumulate at y = −0.4 mm (disc underside). ε_g clips to
   0.333.
3. The 30 ms plug lift is the bed sealing over a jet that cannot form a
   central channel. Pressure builds under the bed and releases. After the
   collapse, the trapped parcels and the 200 – 500 m/s ring cut-cells force
   dt down by ~30×. They also produce the kHz pressure noise and the negative
   mean ΔP: the inlet monitor sits in the low-pressure ring jet / recirculation.

The same membrane is present in **every** implicit-SDF STL in the repo:
`cases/V_3d/day10_r3_smoke` (2 mm: 38 disc facets, r = 0.47 – 1.89 mm),
`cases/V_3d/day11` (DEM), `cases/D_sweep13/case_00` (sweep). Day-11 and
Day-13 results inherit it. The R1 – R3 "orifice velocity only 12 m/s on the
2 mm mesh" diagnosis may also be partly this artifact rather than resolution
alone. That is not re-tested here.

## Contributing (secondary) factor

`dt_max = 1e-3 s` is hard-coded in `src/make_case.py` (both V_2d and V_3d
builders). This violates CLAUDE.md rule 2 and spouted-bed recipe Rule 9
(`dt_max <= 0.5·dx/U_in = 1.7e-5 s`). The "Fix B" dt cap described in
`notes/day10_pathology_review.md` was never wired into the builder. Before the
collapse, dt grew to 2.9e-4 s (CFL ≈ 9 at 30 m/s; far higher in the ring
gap). This did not create the membrane, but it is a deviation that should be
moved into params and decided on deliberately.

## Proposed fix (to be implemented and smoke-tested next)

**Primary:** remove the y = 0 seam from the SDF union by overlapping the
primitives. Extend the stub primitive's top cap to y = +δ and the cone's
bottom cap to y = −δ, with δ = dx_sample. The stub radius (R_stub = R_i) is
≤ the cone radius for y > 0, so the overlap changes no wall location. It only
eliminates the zero-valued interior plane. Add a regression test asserting no
facet with |n_y| > 0.7 lies at 0 < r < R_i − dx_s anywhere except the inlet
and top caps.

**Secondary:** move `dt_max_s` (and `dt_s`) into `params/numerics.yaml` with
a source comment. Keep the numeric value unchanged for the first smoke, so the
membrane fix is tested in isolation (CLAUDE.md rule 11). Then decide on the
Rule-9 cap using the smoke's dt/NIT data.

Rejected or not retried (per `.claude/skills/mfix-spouted-bed-recipe`):
graded mesh, 2 mm mesh, widened stub, taper, widened gas column,
single-rectangle IC, TFM. None of them address a wall that should not exist.

## Update 2026-10-16 (after smoke R5a): the membrane was only half the story

R5a = R4 deck + membrane-free STL only (`cases/V_3d/day10_r5_smoke`, 5 ranks,
tstop 0.08 s, VTK every 10 ms). Through t = 30 ms the orifice behaves like a
real jet: mean v_g/bc_v_g = 0.84-0.87 at the orifice plane, peak 31-40 m/s,
ε_g 0.82-1.00, no parcels below y = 0, mid-bed centerline v_g 18-25 m/s.
**dt does not collapse** (2.1-3.2e-4 s through t = 49 ms vs R4's 4-9e-6 s).
The membrane was therefore the cause of the dt collapse, the 200-500 m/s
ring jet and the stagnant stub pocket.

**But the bed still slumps and seals the orifice at t = 30-40 ms.** dP rises to
+23 kPa (30-40 ms) and +26.5 kPa (40-50 ms), about 65× W/A. Mid-bed v_g falls
to 4.9 m/s. At t = 0.04 s there are 146 parcels below y = 0 and orifice
ε_g = 0.45.

Second mechanism (IC): R1-R4 seed the bed quadrants at ε_s = 0.57, only 0.01
below PIC close packing (1 − des_epg_clip = 0.58). The Harris-Crighton
solids-pressure model is near-singular there. The packed strips sit beside
an empty (ε_s = 0) central gas column and empty orifice. The resulting stress
gradient pushes parcels into the low-ε_s column and down into the orifice,
which seals. The fixed-flow MI then pressurizes the stub until the plug
lifts. NETL's shipped `tutorials/pic/spouted_bed_3d` uses the same quadrant
pattern but seeds it **loose, ic_ep_s = 0.2**, and lets it settle under
gravity onto a jet that is already flowing.

Fix 2 (R5b): seed the PIC quadrants at `spouted_bed_V_3d_pic.eps_s_ic_pic`
= 0.2 (NETL value) up to H_ic = 48.35 mm. That is the height that conserves
the 54.5 g inventory, extending past the 39.8 mm cone top into the column.
Monitors stay at the settled H_static = 32.2 mm. DEM template unchanged
(no Harris-Crighton singularity; contact model).

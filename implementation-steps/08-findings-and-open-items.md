# 8. Findings, dead ends and open items

**No C26 case is validated.** The numbers for both primary sweeps are in
the [README](../README.md). This file explains *why* the model behaves as
it does and what is left to try.

## Where the two primary inlets stand

| | 2 × 2 control (`_apex_v232`) | Single-cell bore (`_fullface_bore_mom_nointerp_v232`) |
|---|---|---|
| c02 map (12 slpm) | matches on the in-repo PIV (0.914 / 0.883 / 0.0015); PIVlab 0.778 / 0.825 / 0.0048, below the CC floor | no match (PIVlab 0.829 / 0.732 / 0.0075: SSI and MSE fail) |
| c01, c03 maps | no match (short fountains) | no match; C26's shape, 1.3–1.6× C26's model speed |
| Wall velocity | 2.9–6.3× C26's model | 0.64–1.06× C26's model, the first layout at that level |
| Spout | pulses, closes ~⅓ of the time | one-cell column, almost never closes |
| Coated beds at 12 slpm | c04 intermittent; c05–c07 do not spout | c04 internal jet; c05–c07 do not spout |

## What the runs showed

**1. The annulus is fluidized by a spout-closure cycle** (control and
every earlier layout). The "fluidized annulus" behind the fast wall layer
is a 10–12 Hz limit cycle with four steps:

1. Lower-cone wall-layer parcels slide down the 60° wall and plug the
   spout base.
2. The gas lifts the annulus; the support ratio reaches 1.6–2.9.
3. The plug leaves as a slug at ~0.44 m/s.
4. The spout reopens and the annulus voids (ε_g 0.69–0.74).

The jet's dynamic pressure (~80 Pa) is far below the annulus lateral
stress (280–830 Pa), so the gas cannot hold the spout wall open. The
gas-fraction floor and the interpolation scheme do not remove the cycle.
C26's own model pulses too, but at ⅓–⅕ of our amplitude with the same
frictions. An earlier hypothesis, a missing retarding force on fountain
parcels, was rejected: drag on them is ≤ 5 % of their weight.

**2. "Describes" vs "reproduces."** A one-cell momentum jet is what C26
describes (file 04) and the only kind that keeps the spout base clear. The
2 × 2 mass-only source is what reproduces C26's c02 map. The control sweep
was run on the second (2026-10-26) and the single-cell sweep on the first
(2026-10-27).

**3. The single-cell bore overdrives the spout, and the excess is
physical.** All of the inlet momentum goes up the axis (up share 1.0). The
jet never closes (closure duty 0.02–0.14 against 0.28–0.36 in the
controls). The spout stays one cell wide and carries few solids
(circulation 8.5–10.8 g/s against 20–32 g/s). The gas therefore throws a
small solids flux high: the fountain is 2.5–3× the experiment's height at
12 and 14 slpm, while the controls are 0.6–0.7×. The annulus stays packed
(ε_g 0.44), which is why the wall layer is right.

The measurement is not the cause: a parcel-truth map agrees with the PIV
within 6 % (file 07). Across the c02 source variants, the map mean over the
experiment ranks with how much inlet momentum survives above the bed:

| c02 source variant | Map mean / experiment |
|---|---|
| 2 × 2 + momentum, LINEAR_HAT | 0.88 |
| 2 × 2 + momentum, NONE (open apex, `_momentum24_nointerp`) | 0.95 |
| 2 × 2 mass-only, LINEAR_HAT (the control) | 1.00 |
| one-cell point, LINEAR_HAT (`_point_v232`) | 1.09 |
| one-cell point + momentum, LINEAR_HAT (`_point_mom_v232`) | 1.19 |
| full-face bore + momentum, NONE | 1.44 |

**4. Interpolation is not the fountain knob.** At matched times
(0.3–0.5 s) the bore with LINEAR_HAT and with NONE has the same fountain
(99th percentile 72 vs 74 mm), jet and parcel speeds. The remaining wrong
quantity is the momentum sink inside the bed: too few solids enter the
spout and too little gas crosses into the annulus. On this grid that is
governed by the drag closure at annulus packing, where TGS is 2.2–2.6× the
other closures (file 05).

**5. Coated beds do not spout at 12 slpm.**

- The Mathur–Gishler margin separates them: spouting cases sit at
  U/U_ms ≥ 0.56, non-spouting at ≤ 0.54.
- Flows that would give each coated stage c02's margin (predictions, not
  run): c04 13.6, c05 15, c06 17, c07 18.4 slpm. In an early smoke c07
  spouted at 18 slpm.
- Bed depth and particle size change together in this block, so the two
  effects cannot be separated.
- On the control, more charge gave a smaller fountain share and a slower
  wall (c08 → c02 → c09). Circulation stayed flat at 25–28 g/s, so turnover
  time scales with bed mass.

**6. The 27-cell mass-only spout leaned −x.** About ⅔ of the gas left the
source cell sideways, and the side outflows were unequal from the first
frame. The suspected cause is MFiX's unlimited x-interpolation factor
(`Theta_Ue` in `get_alpha.f`) in cut cells. That was never proven, and the
centred layouts made it moot.

## Dead ends on C26 (do not retry)

- **A point source.** A zero-width point at a cell centre: it snaps into a
  sliver cell. A typed `ps_v_g` sweep in a walled source cell: continuity
  fixes the exit velocity. A momentum source on the open-apex 2 × 2 layout:
  only the source-cell pressure changed. A one-cell point on the
  face-centred grid for production: 1.5 mm off axis, 5 mm map lean, dt
  cuts.
- **A closed apex on the odd grid without a bore**: the tip cell is
  blocked (`PS_VOLUME == ZERO`).
- **The flat-ledge bore**: parcels rest on the ledge, and it is slow.
- **Removing `des_epg_clip`**: Error 1100. Raising `TOL_SMALL_CELL` does
  not help.
- **`remove_rogue_particles = .F.`** against 23.2 parcel deletions: it
  becomes fatal.
- **The DMP partition or the STL** as the cause of a spout lean: both
  ruled out.
- **A start-up ramp** to stop the closure cycle: the cycle regenerates
  every period.
- **PIV.** The "empty = 0" PIV convention. Moving the map bottom to the
  best-scoring height (that is tuning). Using the source-cell ΔP monitor
  as a bed ΔP.
- **W < 10, or a grid finer than 2.1 mm**: out of scope (user decision).

## Open items

- **Drag-closure screen (running, submitted 2026-11-01).** Three 0.5 s c02
  bore smokes on AWS, each changing one keyword: Gidaspow (annulus drag
  0.39× TGS), BVK (0.45× annulus, 1.49× spout core), and TGS with a 0.50
  gas-fraction floor. Criteria are in `params/c26_acceptance.yaml`
  `drag_closure_smoke`.
  - *Supports* "annulus cross-flow is drag-limited" if fountain p99 is
    ≤ 50 mm and the spout gas share at 30 mm is ≤ 0.60.
  - A closure that brings the fountain to the experiment's 24–29 mm with a
    packed annulus is the candidate for a 3 s run.
  - TGS stays the production closure (C26 Table 1). A closure change is a
    sensitivity screen.
- **Targets for the next run at 12 slpm**, independent of the PIV mapping
  and written before any run:
  - fountain p99 24–29 mm above the annulus (C26 Fig. 3a);
  - launch velocity p90 ≈ 0.45 m/s;
  - axis net-down frames > 0.15;
  - wall ≤ 2× C26's model.
- **Questions for the C26 authors.** Was the source cell whole or cut by
  the cone? What velocity was typed at 10 slpm? Where did the axis sit on
  the grid?
- **Runs not yet done.** The one-cell momentum point source with NONE
  interpolation. The bore with LINEAR_HAT for 3 s, to separate source from
  interpolation. A friction screen (wall friction 0.3, then rolling
  friction 0.03, then particle friction 0.3), which would move the model
  toward the experiment but cannot explain the gap to C26's model.
- **Pipeline fixes.**
  - Move the inlet-pressure monitor one row above the source cell.
  - Make the PIVlab scorer skip c04–c09, which it scores against C26's
    12 slpm panel (those numbers are not reported).
  - Adopt a map-bottom rule that ignores airborne parcels, so c08's map is
    usable.
  - Check the circulation measure: on the bore runs, downward crossings are
    only 0.51–0.65 of upward ones.
- **Unexplained.** The 2 × 2 + momentum closed-apex run scored worse than
  the control. c06 on the earlier layouts showed a bed expansion (surface
  ~61 mm against 47.7 mm static).
- **Production defaults.** `params/geometry.yaml` (`ps_momentum`,
  `ps_point`, `stl_builder`) still hold the pre-2026-10-24 settings
  pending a user decision. Each primary deck sets its layout explicitly.

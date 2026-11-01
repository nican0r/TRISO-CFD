# 5. Solver and numerics

## MFiX release: 26.1.2 → 23.2

The project started on MFiX 26.1.2. C26 ran 23.2, so since 2026-10-24
every new run uses **23.2** (`params/numerics.yaml` `mfix_version`).

- **Deck changes.** The 23.2 keyword list (980 keywords against 1,159)
  was parsed from its namelists. Four 26.1.2-only keywords are dropped:
  `mu_g_model`, `bc_default_walls`, `ic_des_abort_on_failed_seeding`,
  `ic_des_tight_boundary_spacing`. 23.2 aborts on unknown keywords
  (Error 2000).
- **Physics check.** A source diff showed drag, CGP scaling, gas-fraction
  deposition, interpolation weights and point-source terms are identical in
  the two releases. Per-particle wall damping differs by < 0.1 % of the
  dissipation. Seeding and the mesh are identical.
- **Where it runs.** 23.2 runs in the AWS container image
  `triso-cfd/mfix:23.2`, which builds the DMP solver from source. Locally
  the same image runs under docker. `src/aws/submit.py` refuses a deck
  whose release does not match the image. Re-rendering an old case:
  `run_C26_3d --mfix-version 26.1.2`.
- **One real difference.** 23.2 deletes parcels at some cut-cell corners
  (file 03). Check the parcel count on every 23.2 run.

## Particle–mesh interpolation: LINEAR_HAT → NONE

The interpolation scheme spreads each parcel's volume and drag onto the
mesh. C26 does not report it.

- **Until 2026-10-27:** `LINEAR_HAT` with mean-field interpolation,
  carried over from the Day-10 spouted-bed recipe. Used by the control
  sweep.
- **Since 2026-10-27:** `NONE` (cell-centered deposition, MFiX's own
  default; user decision), written with `des_interp_on` and
  `des_interp_mean_fields` false. Used by the single-cell sweep.
  `run_C26_3d --interp-scheme LINEAR_HAT` re-renders older decks.
- **Evidence.** On the 2 × 2 momentum layout, NONE halved the wall
  velocity (0.026 → 0.013 m/s) and left the annulus fluidization and spout
  gas share unchanged. The pre-registered verdict was inconclusive.
- **Fountain.** On the bore, LINEAR_HAT and NONE give the same fountain,
  jet and parcel speeds at matched times (0.3–0.5 s). Interpolation acts on
  the annulus and wall layer, not on how much momentum the jet keeps.
- **Cost.** NONE needs more outer iterations per step and smaller steps:
  the single-cell 3 s runs took 1.5–2.5× as long as the control sweep.

## Gas-fraction floor: `des_epg_clip = 0.42`

`des_epg_clip` floors the **whole** gas-fraction field. It is not just a
drag limiter: `EP_G = 1 − min(1 − clip, ΣEP_S)` enters gas continuity and
momentum too. An `ep_star` in a DEM deck silently becomes the clip.

- **Why it is needed.** Sliver cut cells hold 1.0–1.8 % of a cell, and one
  parcel is bigger than that. Without a floor the solids fraction exceeds
  1 (Error 1100) at t = 0 with LINEAR_HAT and at 5 ms with SQUARE_DPVM.
  Depositing recorded parcels the way NONE does puts ~10 cells per frame
  at ε_s ≥ 1 in every bore and closed-apex frame. Raising
  `TOL_SMALL_CELL` would not help, because a third of those cells hold two
  or more parcel centers.
- **How much it binds.** With 0.42 the floor binds in ~900–1,250 cells per
  frame, mostly in the packed annulus. 0.42 is assumed (copied from
  Case V). C26 does not report a floor.
- **As a lever it is weak.** A 0.30 floor moved the wall velocity +15 %
  (inconclusive). A 0.50 floor (the densest state inside the TGS drag
  fit) is running as a comparator in the drag-closure smokes (file 08).

## Drag closure

TGS (C26 Table 1). Under coarse graining MFiX evaluates drag at the real
particle diameter and multiplies by the parcel volume; contacts have no CGP
term. Evaluating drag at the parcel diameter would raise U_mf 3.4–3.7×.

At the packed annulus (ε_g 0.44) TGS is 2.2–2.6× Wen–Yu, Gidaspow, BVK
and Koch–Hill, and 4.3× Syamlal–O'Brien. Most of that comes from its
5.81 φ/ε³ term, which is not an extrapolation artifact: inside its fitted
range (ε_g 0.50–0.60) the others are still 0.42–0.60× TGS. In the spout
core (ε_g 0.95) Wen–Yu and Gidaspow equal TGS, and BVK and Koch–Hill are
~1.45×. A drag-closure sensitivity screen is running (file 08). TGS stays
the production closure.

## Time step

dt_max = min(1e-5 s, 0.5·dx/U_nozzle). C26's 0.01 ms binds at every flow;
the CFL bound would bind only above ~28 slpm. The DEM sub-step is 2.13e-6 s
(1/50 of the collision time). A 1e-3 s cap, as in the Case V decks, must
never be copied (file 01).

## Run cost

| Run | Hardware | Wall time |
|---|---|---|
| 3 s, control sweep (`_apex_v232`, LINEAR_HAT) | AWS, 4 vCPU, 4 ranks | 11.9–16.7 h |
| 3 s, single-cell sweep (`_fullface_bore_mom_nointerp_v232`, NONE) | AWS, 4 vCPU | 12.1–29.9 h |
| 0.5 s smoke | AWS, 4 vCPU | ~2.5 h (LINEAR_HAT); expect 3–6 h with NONE |
| 1 ms dry run | laptop, docker, 4 ranks | ~20–40 s |

- **Cost driver.** DEM wall contact against STL facets. Coarsening the
  implicit STL 4× in facets made the DEM step 3.5× cheaper with identical
  cut cells.
- **Laptop vs AWS.** An AWS c7a core is 1.6–1.7× slower per step than the
  laptop's M4. Bind/MPI settings in the container were suspected but not
  tested.
- **Sizing.** 4 ranks on 4 vCPU and 6000 MB (8000 MB for c09) fits 8 jobs
  in the 32-vCPU quota. Memory peaks near 1.1 GB.
- **Timeouts.** Batch jobs carry no timeout and never retry. C26 needed
  167 h for 3 s at W = 10 (Table 5).

# Day-10 R5 — independent review of the loose-IC fix (R5b)

Written 2026-10-16 while R5b was running. Analysis only: no MFiX runs and no
deck edits. Question: is seeding the PIC bed loose (ε_s = 0.2, H_ic = 48.35 mm)
the right fix for R5a's 30–50 ms slump? Or is something more consistent with
NETL and the spouted-bed literature better?

## Verdict: partially right

Loose seeding removes the dominant **numerical** driver of the slump. It is the
right next test. But the justification "NETL seeds 0.2 around an empty central
column" misreads the NETL tutorial. Two pieces of R5a's behaviour are physical
rather than artifacts and should not be "fixed": part of the ΔP peak, and the
breakthrough that was under way at 80 ms.

## 1. NETL ground truth

Source: `/Users/nelsonpereira/mamba/envs/mfix-26.1.2/share/mfix/templates/tutorials/pic/spouted_bed_3d/`
(`spouted_bed_pic_3d.mfx`, `geometry.stl`).

- **NETL's tutorial is a draft-tube spouted bed.** Its STL contains a solid
  tube with inner r ≈ 7.3 mm (inward normals) and outer r ≈ 11 mm (outward
  normals), running from y = −150 to −52 mm (end caps at both heights). Below
  it is a ~35 mm entrainment gap above a 15 mm inlet tube (r = 7.5 mm,
  y = −200 to −185 mm). The cone half-angle is ≈ 26.5°; the column is
  D = 100 mm. The tube was identified by sampling STL facets; it is not
  mentioned in the deck.
- **ic_ep_s = 0.2 is confirmed** for the 4 strips (IC 2–5, y = −0.199 to
  0.0). The excluded centre is |x|, |z| < 7.5 mm, i.e. **the draft-tube
  bore**. In NETL that column is held open by a wall. Our vessel (ORNL/UTK
  cold mockup) has no draft tube, so the R1–R5 "Fix A" quadrant pattern
  copied a geometric feature that has no support in our geometry. Over its
  full 32 mm height, the empty 6 mm column is an unsupported vertical free
  face of granular material.
- **Every closure is identical to our deck** (R5a/R5b verified by grep):
  SYAM_OBRIEN with drag_c1/d1 = 0.8/2.65; des_epg_clip = 0.42;
  friction_model = 'SCHAEFFER'; psfac_fric_pic = 100; fric_exp_pic = 3;
  fric_non_sing_fac = 1e-7; mppic_coeff_en1 = en_wall = 0.85;
  et_wall = 1.0; velfac = 1.0; ic statwt 4.2; bc statwt 5; LINEAR_HAT;
  dt_max = 1e-3. The closures therefore do not explain the different
  behaviour.
- **What does differ:**

  | | NETL | Case V |
  |---|---|---|
  | Draft tube | yes | no |
  | D_c/D_i | 100/15 = 6.7 | 50/4 = 12.5 |
  | Cone half-angle | 26.5° | 30° |
  | d_p, ρ_p | 1.75 mm, 2000 | 0.5 mm, 6050 |
  | Jet velocity | 50 m/s | 30 m/s |
  | Cells across orifice | 3.75 | 4 |
  | dx/d_p | 2.3 | 2 |
  | Inflow ramp | none (impulsive) | none (impulsive) |
  | IC packing | 0.2 loose | **0.57 (R1–R5a)** |
  | tstop | 1.0 s | — |

## 2. Is R5a's +23–26.5 kPa transient physical, numerical, or both?

**Both (c).** The trigger is numerical; the size of the peak is plausible
physically.

- **Numerical trigger.** Harris–Crighton solids pressure
  (`des/pic/calc_ps_pic.f90:85-86`):

  P_s = psfac·ε_s³ / max(0.58 − ε_s, 1e-7·ε_g)

  At ε_s = 0.57 this gives P_s ≈ 1.85 kPa; at 0.2 it gives ≈ 2 Pa. Across the
  1 mm interface to the empty column (ε_s = 0), the gradient is ~1.9 MPa/m.
  The solids weight per volume is ρ_s ε_s g ≈ 34 kPa/m, so at t = 0 the
  lateral push is **~55× gravity**. That is not physics: a packed granular
  wall at rest carries only gravity-scale stresses. Observed in R5a over
  0–30 ms:
  - the annulus moves down at −0.075 to −0.095 m/s;
  - the mean particle height drops 23.7 → 21.4 mm;
  - core parcels go 0 → 89, then 602 at 60 ms.

  The column fills faster than a 30 m/s pipe jet can clear it. Seeding near
  the clip next to ε_s = 0 is the artifact. Loose seeding removes it directly.
- **Physical magnitude.**
  - The "65 × W/A" ratio is misleading for a cone: W/A uses the top area
    (401 Pa). The orifice area is ~1/106 of the top, so W/A_orifice ≈ 42.5 kPa.
  - Our own Ergun estimate for the full 22.6 L/min permeating a packed 60°
    cone radially from the orifice: virtual apex at r₀ = 4 mm, solid angle
    0.84 sr, ε = 0.43. It gives viscous ≈ 5 kPa + inertial ≈ 31 kPa ≈ 36 kPa.
    This is an order-of-magnitude estimate only, not validated.
  - R5a's peak (26.5 kPa = 0.62 × W/A_orifice) is therefore in the physically
    plausible range for incipient spouting from a packed conical bed.
  - The literature agrees qualitatively: the startup peak is well above the
    steady spouting ΔP and grows with bed height; hysteresis comes from the
    change in granular stress state after collapse
    (OSTI 1477173; Springer DOI 10.1007/978-1-4684-1045-7_25).
  - **Weak point:** I could not verify the numeric form of the Olazar /
    Saldarriaga ΔP_M correlation (J. Taiwan Inst. Chem. Eng. 80 (2017)
    678–685, DOI 10.1016/j.jtice.2017.09.001). No literature ratio is quoted
    here; verify it from the paper before any REPORT claim.
- **Experiments avoid this peak.** S2 reaches the operating point by raising
  the flow, then holds an established spout; the measured 17–25 Hz data are on
  the established-spout (hysteresis) branch. An impulsive 30 m/s start on a
  packed bed has to pass through ΔP_max once. That is expected, but it is only
  startup cost and must be excluded from the averaging window.

## 3. R5a after the slump: still evolving, heading toward breakthrough

From the VTK frames (10 ms) and monitor CSVs:

| t (ms) | mean particle y (mm) | annulus v_y (m/s) | parcels y > 35 mm | stub parcels | orifice ε_g (y = 2 mm) | ΔP (Pa, 5 ms mean) |
|---|---|---|---|---|---|---|
| 50 | 20.70 | +0.022 | 23 | 236 | 0.39 | 23 155 |
| 60 | 21.26 | +0.061 | 48 | 136 | 0.47 | 8 119 |
| 70 | 22.15 | +0.084 | 104 | 34 | 0.83 | 1 571 |
| 80 | 23.10 | +0.080 | 193 | 12 | 0.87 | 783 (≈ 2 × W/A) |

- The bed is lifting.
- The fountain is growing.
- ΔP falls monotonically towards spouting-order values.
- An internal jet cavity of ~10–12 mm forms above the orifice
  (v_g > 4 m/s up to y ≈ 10–12 mm).

This is the textbook incipient-spouting sequence. It is **not** a stable
non-spouting state. Two caveats:
- the trend could also be a whole-bed plug rising that later falls back;
- 0.08 s is too short to tell.

Also note that in R5a's 0–30 ms the "jet" kept 25–39 m/s **above the bed
surface** (y = 36 mm). That is pipe flow up the pre-cut column, not a spout.
The pre-slump state was the artificial one.

## 4. Critique of R5b (ε_s = 0.2 to 48.35 mm, 6 mm central gap kept)

- **Good:**
  - removes the ~55 × g interface push (P_s ≈ 2 Pa);
  - matches NETL's IC packing;
  - the bed settles onto an already-flowing jet, which approaches the
    established-spout branch from above as a draft-tube-free analogue of
    NETL's start;
  - closures unchanged.
- **The kept central gap** is harmless at 0.2, but still has no physical basis
  here (it is NETL's draft-tube bore). It does not need removing for R5b.
  Future decks should keep only a small orifice-footprint exclusion, which was
  Fix A's original purpose: no parcels seeded in or over the orifice cells.
- **Risks to watch:**
  1. **Settling time.** The top must fall ~16 mm; free fall alone takes
     ≈ 57 ms, and hindered settling is slower. Expect the bed to settle by
     t ≈ 60–100 ms. tstop = 0.1 s will show the end of settling at best, not
     stationary spouting.
  2. **Whole-bed fluidization of the loose lower cone** during settling. This
     would show as ΔP ≈ W/A with bubbling and no distinct annulus downflow.
  3. **Parcels raining into the 30 m/s jet** are thrown into the freeboard: a
     big transient fountain. This is benign if they return.
  4. **Settled packing** ends near 0.58 (the HC limit), same as R5a. That is
     expected, and a late packing-driven slump cannot be excluded until the
     run passes settling.
- **Acceptance window for R5b:**
  - t_settle = first time the mean particle height and the 95th-percentile
    surface change < 0.5 mm over 20 ms;
  - evaluate the smoke checks on [t_settle + 20 ms, tstop];
  - require tstop ≥ t_settle + 100 ms, i.e. **≈ 0.2 s**.
  - At t = 0.1 s R5b can only answer "did the bed settle without sealing the
    orifice / spiking ΔP?". It cannot answer "is spouting stationary?".

## 5. Ranked options

All keywords below were verified in the 26.1.2 sources:
`share/mfix/src/model/init_namelist.f` (line given) or the files named.

1. **Loose IC, as R5b.**
   - **Keep:** `IC_EP_S`, `IC_EP_G`, `IC_PIC_CONST_STATWT` (already in the
     deck).
   - **Basis:** NETL tutorial; removes the HC interface artifact.
   - **Cost:** +60–100 ms settling transient.
   - **Follow-up for production decks:** shrink the central exclusion to the
     orifice footprint near y ≈ 0 instead of full height.
2. **Add an inflow ramp only if R5b shows a startup shock.**
   - **Keywords:** `BC_V_g_ramp_t0/t1/v0/v1` (init_namelist.f:3477, 3483,
     3489, 3495). Applied every step for MASS_INFLOW in `set_bc1.f:101-104`
     and `673-685`, linear between t0 and t1, held at v1 afterwards
     (`set_bc1.f:723-737`). Mutually exclusive with `BC_JET_G*`
     (`check_boundary_conditions.f:372-379`). `MAX_INLET_VEL` is recomputed
     every step (`time_step.f:123`).
   - **Basis:** S2's experimental procedure (flow raised to the operating
     point).
   - **Caution:** a ramp on a **packed** bed reproduces the physical ΔP_max
     rather than avoiding it. It does not fix a slump caused by the IC
     gradient, and it lengthens startup by the ramp time. Low priority for
     production. Useful later as a check that the model's ΔP_max is
     consistent with the literature.
3. **"Established spout" IC.**
   - **Layout:** a dilute spout channel (ε_s ≈ 0.05–0.1, D_s ≈ 4.5 mm from
     McNab, D_s = 2.0 G^0.49 D_c^0.68 / ρ_b^0.41 with G = 0.226 kg/m²s,
     ρ_b = 3448 kg/m³). Annulus seeded below the clip (ε_s ≈ 0.5,
     P_s ≈ 156 Pa) with a stepped interface.
   - **Keywords:** the same IC keywords.
   - **Basis:** hysteresis / established-branch argument.
   - **Trade-off:** shortest possible transient if the guess is good, but
     contrived. Hold in reserve if R5b's settling proves too slow or too
     violent.
4. **Dense IC below the clip (0.50–0.55) with no central column, impulsive jet.**
   - **Effect:** removes the numerical push. But the run must then pass
     through the physical ΔP_max (~tens of kPa) and a breakthrough lasting
     ~50 ms, as R5a did. That costs more startup and gives worse dt than
     option 1.
5. **Not recommended: retune PIC stress parameters** (`PSFAC_FRIC_PIC`,
   `FRIC_EXP_PIC`, `FRIC_NON_SING_FAC`, `MPPIC_COEFF_EN1`; definitions in
   `des/des_init_namelist.f`; used in `calc_ps_pic.f90:85-86`).
   - These are NETL's validated defaults and act everywhere (packing,
     annulus flow).
   - There is no case-specific literature basis to change them.
   - Changing them to suppress a startup artifact would bias the
     production physics.

## 6. What R5b's outcome will and won't tell us

- **Will tell us:** whether loose seeding prevents the orifice seal and ΔP
  spike during startup, and whether dt stays O(1e-4) through settling.
- **Won't tell us:**
  - whether a stationary spout–fountain–annulus regime forms (it needs
    ~0.2 s);
  - PSD or circulation numbers;
  - whether R5a would have reached the same state on its own.
  If R5b is clean, that last point matters only for cost: R5a's route
  costs ~50 ms of extra transient.

## 7. Next smoke I would run

After R5b, on the same 5 cores (not in parallel):

1. **Extend R5b to tstop = 0.2 s** (RESTART_1, `init_namelist.f:113`), so the
   check window covers [t_settle + 20 ms, 0.2 s]. This is the decision run
   for AWS.
2. **Optional, cheap: restart R5a from its t = 0.08 s `.RES` to 0.15 s.**
   This settles whether R5a breaks through into spouting. It tells us
   whether the slump was only a long, physical startup rather than a dead
   end.
3. **If R5b settles badly** (whole-bed fluidization or a late seal), try
   option 3 (the established-spout IC) next. Try option 2 (the ramp) only
   if the startup shock itself is the problem.

## Sources

- NETL tutorial: `share/mfix/templates/tutorials/pic/spouted_bed_3d/spouted_bed_pic_3d.mfx`
  and `geometry.stl` (MFiX 26.1.2 install).
- MFiX source: `share/mfix/src/model/init_namelist.f`, `set_bc1.f`,
  `time_step.f`, `check_data/check_boundary_conditions.f`,
  `des/pic/calc_ps_pic.f90`.
- Saldarriaga et al. (Olazar group), "Correlations for calculating peak and
  spouting pressure drops in conical spouted beds of biomass", J. Taiwan Inst.
  Chem. Eng. 80 (2017) 678–685, DOI 10.1016/j.jtice.2017.09.001. Equation not
  verified: https://ekoizpen-zientifikoa.ehu.eus/documentos/5ee0ca532999527595f996d7?lang=en
- Analysis of maximum pressure drop for a flat-base spouted fluid bed (NETL),
  OSTI 1477173: https://www.osti.gov/biblio/1477173
- Spouted-bed startup peak / hysteresis chapter: https://link.springer.com/chapter/10.1007/978-1-4684-1045-7_25
- Mathur & Epstein (1974) *Spouted Beds* and Epstein & Grace (2011) *Spouted
  and Spout-Fluid Beds*: cited for background only. Not checked in this
  review.
- McNab (1972) spout-diameter correlation: quoted from memory. Verify before
  use.

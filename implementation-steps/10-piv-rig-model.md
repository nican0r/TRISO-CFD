# 10. PIV rig model: changes from the C26 model

Branch `piv-rig-model`. This adapts the C26 model to the benchtop rig built
to replicate the C26 experiment, using the measurements in `data/piv_rig/`
(IDs such as T1 and N1 below are rows of its Apparatus sheet, measured
2026-11-02).

**Status: r01, r02, r03 and the 20 slpm step-down were submitted to AWS on
2026-11-05. r03 failed at t = 0.13 s (see "r03 failure") and was resubmitted the
same day from its t = 0.12 s state with a larger pressure-solve sweep limit, as
`r03_leq200`. On 2026-11-06 r01 and r02 finished their 3 s and both scored
`SPOUTS`, r01 only marginally (see "3 s results"). The step-down scored
`in_between`, so it does not yet decide the drag question (see "Step-down
result"). `r03_leq200` is still running. Also on 2026-11-06 the flows turned out
to be air-scale rotameter readings, so r01–r03 ran 5.5–6.7 % below the rig's test
conditions; `q1`–`q3` hold the corrected flows and were submitted the same day
(see "Flow scale correction"). Nothing here is validated against the
experiment.**

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
| Wall below the rim | one 60° cone from the throat | **measured profile, table below** | N1, N2, F1, F12, F8, 60° |
| Cone top | r = 25.2 mm at y = 41.6 mm | **r = 27.5 mm at y = 47.76 mm** | F8, 60° |
| Rim pocket and step | none | **vertical wall to y = 50.0 mm, then a 2 mm downward step to the tube** | F8, F11, T1; your choice |
| Outlet height | 210 mm | **354.57 mm** (tube top 354.80 mm, nearest whole row) | F11 + T3 |
| Cells across column | 24 (2.1 mm) | **21 (2.4286 mm)** | your choice |
| Box | 25 × 104 × 25 | **23 × 150 × 23** (50,407 fluid cells; C26 42,520) | derived |
| Source face vs orifice | 4.41 vs 4.45 mm² (−1 %) | **5.90 vs 5.73 mm² (+3 %)** | derived |
| Gas | Ar, 25 °C, 1.633 kg/m³ | **N₂, 25 °C, 1.145 kg/m³, 1.80e-5 Pa·s** | your instruction |
| slpm reference | 25 °C, 1 atm, Ar | **25 °C, 1 atm, N₂** | your choice |
| Particle | YSZ, 411.6 µm, σ 12.2 µm | **steel, 822 µm, σ 10.3 µm** (±3σ) | Particles sheet |
| Density | 6000 kg/m³ | **7585 kg/m³** | Particles sheet |
| Statistical weight W | 10 | **1** (C26 Table 2 steel) | your choice |
| Parcels | 27,317 estimated | **27,188 estimated**; 27,028 seeded (59.64 g) | 27,200 × 2.2059 mg = 60 g |
| Packed fraction (H_static only) | 0.57 | **0.60**, giving H_static = 33.2 mm | implied by the sheet's 13,184 mm³ |

Mesh note: at a 2.7 mm orifice, 22 or more cells do not fit the bore because
the cone's apex (y = −2.338 mm) falls below the bore floor. So 24 cells could
not carry over; 21 keeps one cell face about equal to the orifice.

## The wall profile

| y (mm) | Wall | Set by |
|---|---|---|
| below 0 | one-cell square bore holding the source cell, as in the C26 single-cell inlet | dx |
| 0 to 11.79 | 60° recess line, from r = 1.35 mm to the glass radius, r = 8.155 mm | N1, 60°, F1 |
| 11.79 to 14.25 | vertical, r = 8.155 mm (the glass bore) | F1; plug top 10.65 + F12 |
| 14.25 to 47.76 | 60° glass cone, to r = 27.5 mm | F8, 60° |
| 47.76 to 50.00 | vertical, r = 27.5 mm (2.24 mm; the lip measures 2.26 mm) | F8, F11 |
| 50.00 | 2 mm flat step, facing down | F8, T1 |
| 50.00 to 354.57 | column, r = 25.5 mm, to the outlet | T1, F11 + T3 |

The stack closes against the measured rim height: 10.65 + 3.60 + 33.51 + 2.26
= 50.02 mm, against F11 = 50.00 mm. It only closes with the remeasured
F12 = 3.60 mm; the first reading of 1.50 mm left a 2.08 mm gap.

**The ledge left out.** At the plug top (y = 10.65 mm) the real wall has a flat,
upward-facing step 0.65 mm wide: 0.32 mm of plug wall top (N2 to N3), then the
0.33 mm clearance gap between plug and glass (N3 to F1). The model instead
continues the recess line up to the glass radius. That changes the bed volume
by 0.14 %, is about a quarter of a cell, and avoids a flat surface for parcels
to rest on, which slowed an earlier bore design badly enough to retire it.

## Cases

The inlet layout is unchanged: one walled source cell on the axis in the bore
row y ∈ [−dx, 0], with momentum, and gas leaving only through the top face.

| Case | Flow | ṁ | Orifice velocity | `ps_v_g` = Q/dx² |
|---|---|---|---|---|
| r01 | 24 slpm | 4.580e-4 kg/s | 69.9 m/s | 67.8 m/s |
| r02 | 26.5 slpm | 5.057e-4 kg/s | 77.1 m/s | 74.9 m/s |
| r03 | 30 slpm | 5.725e-4 kg/s | 87.3 m/s | 84.8 m/s |
| q1 | 25.40 slpm | 4.847e-4 kg/s | 73.9 m/s | 71.8 m/s |
| q2 | 28.19 slpm | 5.380e-4 kg/s | 82.1 m/s | 79.7 m/s |
| q3 | 32.16 slpm | 6.137e-4 kg/s | 93.6 m/s | 90.9 m/s |

All cases use 60 g. q1–q3 are the rig's three test conditions. r01–r03 are the
same three rotameter readings entered as N₂ flows before the scale was known
(see "Flow scale correction"); the rig never ran at them.

## Measured spouting limit

The rig's lowest spouting flow is a reading of **22 on the rotameter's air
scale, with N₂ flowing**, for the 60 g steel charge, found by lowering the flow
from an already spouting bed (your measurement, reported 2026-11-04; scale and
gas confirmed 2026-11-06). That is 22.8 LPM of N₂ at the scale's 20 °C
reference, or **23.2 slpm** at this model's 25 °C reference (see "Flow scale
correction"). Until 2026-11-06 this section read 22.4 slpm, from taking the test
gas to be air.

| | Flow, slpm N₂ | Multiple of the limit |
|---|---|---|
| Step-down | 20 | 0.86 |
| r01 / r02 / r03 | 24 / 26.5 / 30 | 1.03 / 1.14 / 1.29 |
| q1 / q2 / q3 (the rig's tests) | 25.40 / 28.19 / 32.16 | 1.09 / 1.21 / 1.39 |

This is a scoring target, not a model input. The model should hold a spout at
q1 and lose it not far below 23.2 slpm. Because the limit was measured on
falling flow, the matching model test steps down from a running r01 spout
instead of starting from rest.

The plan was to step down from r01's end at t = 3 s. To run alongside r01, it
starts instead from the r01 smoke's end at t = 0.5 s, the only spouting state
that exists yet. That spout is younger, so a collapse is a weaker result than a
hold. If the result is unclear, the same step can be rerun from r01's t = 3 s
state once r01 finishes.

```bash
# continue the r01 smoke at 20 slpm (0.86× the limit) for 2 s
.venv/bin/python -m src.run_C26_3d --rig piv_rig --case r01 \
    --restart-from cases/PIV_rig_3d/r01_smoke05/aws_out \
    --flow-slpm 20 --tstop 2.5 --run-dir-suffix _step20
.venv/bin/python src/post/piv_rig_score.py cases/PIV_rig_3d/r01_step20/aws_out --stage step_down_20
```

`params/piv_rig_acceptance.yaml` sets the pass and collapse thresholds, written
before any rig result existed. A 1 ms restart dry run at 20 slpm from the smoke's
state confirmed the restart takes the new flow from the deck (source velocity
56.52 m/s, was 67.82) and continues from t = 0.5 s with the smoke's 27,012 parcels.

The closure metric d20 needed one rescaling. C26 counts a cell as spout above
2 m/s of gas velocity, but all the rig's gas seeping through the packed annulus
at y = 20 mm would reach 2.3 m/s at 24 slpm, so a dead spout would score as open.
The rig uses 6 m/s, the same multiple of that seepage velocity as C26 at its
10 slpm.

## First smoke

`cases/PIV_rig_3d/r01_smoke05` ran 0.5 s on AWS (24 slpm, TGS, job `64e3a53e`,
3.52 h on 4 vCPU). Solver health passed: exit 0, no MFiX error, dt steady at
9e-6 s, no parcel below the orifice plane. The regime floor passed over
0.3–0.5 s: d20 0.39, dense top 7 mm.

**One check failed as written, and it was replaced before the longer runs.** The
inlet monitor minus bed top read −1483 Pa (−2.89 W/A) against a limit of above 0.
That monitor sits inside the source cell, which reads about ½ρV² below the row
above it on both models, so it does not read the bed:

| Quantity | Rig r01 | C26 c02 |
|---|---|---|
| Source cell minus outlet (0.3–0.5 s) | −1472 Pa | −766 Pa |
| Row above the source minus outlet (0.3–0.5 s) | +1079 Pa, 2.1 W/A | +832 Pa, 2.0 W/A |
| Jump across the source cell's top, as a share of ½ρV² (0.3–0.5 s) | 0.97 | 0.95 |
| Annulus minus bed top (0.3–0.5 s) | 680 Pa, 1.33 W/A | 575 Pa, 1.39 W/A |

The 3 s C26 c02 runs read the same sign on that monitor (TGS −1.75, BVK
−1.70 W/A), so the gate cannot tell a good run from a bad one. The acceptance
file's reason for it, that the monitor reads a positive jet head, was wrong.
You approved the replacement on 2026-11-05, after seeing this result: the cell
directly above the source minus the bed-top slab, positive and within 0.2–3 W/A.
The smoke reads +2.12 W/A on it over 0.1–0.5 s (C26 c02: +2.05).

The exit velocity check passes unchanged: 0.97 of `ps_v_g` in the cell above the
source, with a gas fraction of 0.92 (C26 c02: 0.97 and 0.96, 0.1–0.5 s). It was
first reported as a second failure at 0.27 because the check was run with a 1.8 mm
radius instead of the 1.35 mm orifice radius, which averaged four annulus cells
into the reading.

Lost parcels: 11 of 27,023 (0.041 %), 22 per second. All 11 were last seen at
cone-wall grid nodes and none near the outlet, which matches 23.2's cut-cell
deletion. At that rate, a 3 s run loses about 0.24 %, within the 0.5 % limit.

## Runs submitted

All four went to AWS Batch on 2026-11-05 in parallel, on 4 vCPU each.

| Case dir | Flow | Simulated time | Job ID | Scoring stage |
|---|---|---|---|---|
| `r01` | 24 slpm | 0–3 s | `e01d2eaf` | `run_r01_3s` |
| `r01_step20` | 20 slpm | 0.5–2.5 s, restart | `99543290` | `step_down_20` |
| `r02` | 26.5 slpm | 0–3 s | `2a8aa443` | `run_r02_3s` |
| `r03` | 30 slpm | 0–3 s | `e01813de`, failed at 0.13 s | not scored |
| `r03_leq200` | 30 slpm | 0.12–3 s, restart of `r03` with `leq_it(1) = 200` | `bdce6454` | `run_r03_3s` |
| `q1` | 25.40 slpm | 0–3 s | `e0a06ceb`, submitted 2026-11-06 | `run_q1_3s` |
| `q2` | 28.19 slpm | 0–3 s | `c62e4a6b`, submitted 2026-11-06 | `run_q2_3s` |
| `q3` | 32.16 slpm | 0–3 s | `cb5d35f0`, submitted 2026-11-06 | `run_q3_3s` |

At the smoke's rate, a 3 s run takes about 21 h and the step-down about 14 h. r02
and r03 had no 0.5 s smoke of their own: their decks differ from r01's only in
flow, and each passed a 1 ms local dry run. Their limits match r01's and were set
before any run above 24 slpm. The pressure check may run close to its 3 W/A limit
at r03, because the row above the source carries jet head that rises with flow.

## r03 failure

r03 (30 slpm) stopped with `DT < DT_MIN` at t = 0.1296 s of 3 s, after 74 min. It
has one gas frame past the startup response, so it cannot be scored.

- **What happened:** dt held near 9e-6 s until t ≈ 0.127 s. Then 58 failed
  iterations (`Run diverged/stalled`, NIT ≥ 15) cut it to 1.33e-7 s in about
  2.6 ms. The last three cuts, to 0.97e-7 s, had no failed iteration: every step
  converged at NIT 2 with the residual falling, and the adaptive dt rule walked
  through the 1e-7 s floor. The stall had already ended.
- **What the jet was doing:** cutting up through a plug on the axis (EP_G
  0.42–0.6 over y ≈ 18–30 mm, peak +2100 Pa, about 4 W/A) during startup. The inlet
  held: the cell above the source read 0.98 × `ps_v_g`, with no parcels below the
  orifice plane and no velocity-limit warning.
- **Replay:** a local restart from r03's t = 0.12 s state, with only output settings
  changed (`cases/PIV_rig_3d/r03/replay_t012`), reproduced the same 58 failed
  iterations at the same time, bottomed at 1.3e-7 s, was back above 5e-6 s by
  t = 0.1304 s and reached t = 0.135 s. It also had a brief second stall at
  t = 0.1344 s.
### Cause

The pressure solve was being cut off. The decks set `norm_g = 0` and leave
`leq_it(1)`, the most BiCGSTAB sweeps MFiX spends on one pressure-correction
solve, at its default of 20. MFiX's keyword note says to raise it with that
setting. A solve here needs about 18 sweeps on average, and when one hits the
limit MFiX uses the unfinished result without a warning.

- **What that does:** the leftover mass error holds the pressure residual at
  1–2e-3, just above the 1e-3 tolerance, until the stall detector fails the step
  at iteration 15, 20 or 25. It also shows as one-step pulses on the bed-top
  pressure monitor (−196 to +32 Pa in the first replay, against a steady 3.7 Pa).
- **Why cutting dt does not help:** in the first replay about 8 % of steps failed
  at every dt from 7e-6 to 1.5e-7 s. Each failure cuts dt by 10 % and six clean
  steps are needed to win one cut back, so above roughly that rate dt runs away
  through any floor. A lower `dt_min` would not have saved the run.
- **Not only r03:** by 7.6 h of wall time r01, r02 and the step-down had 61, 64
  and 79 such recoveries, the deepest at 3.9e-6, 3.5e-6 and 2.3e-6 s. The 3 s C26
  runs on `main` have 10–80 each, the deepest at 1.9e-6 s. r03 differs by one long
  episode. Whether its stronger jet made that episode longer is not established.

### Tests of the fix

All are local restarts from r03's t = 0.12 s state (`cases/PIV_rig_3d/r03/replay_*`).

| Replay | Sweep limit | Result |
|---|---|---|
| First replay, to 0.135 s | 20 | 58 + 2 failed iterations, dt to 1.3e-7 s |
| Stall detection off, to 0.135 s | 20 | steps of 21–50 iterations, two iteration-limit failures near 0.134 s |
| Raised limit, to 0.135 s | 200 | 0 failures, at most 5 iterations in 1,725 steps, bed-top pressure 3–4 Pa |
| Captured stalled state at 0.1342 s, 0.8 ms more | 20 | 3 failed iterations, up to 20 iterations per step |
| Same state, 0.8 ms more | 200 | 0 failures, at most 4 iterations |

Replays are not bit-reproducible: identical reruns drift apart at round-off level
within about 1 ms, so the stall appears on some trajectories and not others. That
is why the last two rows start from one captured state. Over 0.120–0.135 s the
raised-limit replay matches the default-limit ones on jet tip height, centerline
gas velocity (29.8 against 29.8 and 30.0 m/s) and inlet pressure. The cost is
about 5 % per step. Not tested: a full 3 s run, and the effect on the scored
quantities.

### Fix and restart

`leq_it(1) = 200` is now set for this rig only (`params/numerics.yaml`,
`spouted_bed_piv_rig_3d_numerics`); C26 decks still render byte-identical.
`r03_leq200` continues r03 from its last restart state at t = 0.12 s to 3 s:

```bash
.venv/bin/python -m src.run_C26_3d --rig piv_rig --case r03 \
    --restart-from cases/PIV_rig_3d/r03/aws_out --restart-unfinished --run-dir-suffix _leq200
```

Its deck differs from r03's only in the restart type and that line. Two edits
were made by hand in the run dir before submitting: the monitor CSVs were cut
back to t ≤ 0.12 s (nine rows from the abandoned 0.121–0.129 s were dropped), and
r03's three frames from before the restart (t = 0, 0.05, 0.10 s) were copied in so
the frame series is complete. A 0.3 ms local dry run passed (restart read at
0.12001 s, 27,019 parcels, source cell 84.775 m/s, no errors).

r01, r02 and the step-down were left running at the default limit, so r03 will
differ from them in this one numerics setting. They can still fail the way r03
did; if one does, it should be continued the same way from its last restart file.

## Step-down result (20 slpm)

`r01_step20` reached t = 2.5 s in 14.4 h with dt at 9e-6 s throughout (132
recovered steps, none near the floor). It ran at the default `leq_it(1) = 20`.

| Stage `step_down_20` | Value | Needed for `held` | Needed for `collapsed` |
|---|---|---|---|
| Dense top, 1.9–2.5 s | 3 mm (both halves 3 mm) | ≥ 8 mm | < 3 mm |
| d20, 1.5–2.5 s | 0.80 (16 of 20 frames) | < 0.5 | ≥ 0.9 |

**Outcome: `in_between`.** The stage leaves this undecided. Health passed over
0.6–2.5 s: the cell above the source ran at 0.98 × ps_v_g with EP_G 0.92 and sat
+2.12 W/A above the bed top; 0.063 % of parcels were lost; none entered the bore.
`check_smoke.py` passed finish, dt, the inlet velocity and noise. Its S6 fails
at −1.26 W/A because it reads the source cell, which the rig's pressure gate
replaced.

What the bed did (context only, from the 0.05 s particle frames; it does not
change the outcome):

- **Fountain:** 12 mm dense top and a 59 mm p99 at 24 slpm just before the step.
  It was down to 3 mm within 0.25 s of the cut and stayed at 3 mm in every
  0.25 s window to 2.5 s (p99 2–5 mm). The packed bed at t = 0 reads 1 mm, so
  the surface is still being stirred, but there is no fountain.
- **Jet at 20 mm:** closed in 50–100 % of the frames in each 0.2 s window from
  0.5 s on, with no trend. The fast cells carry 0.19 of the inlet gas on average.

So TGS lost the fountain below your measured limit, as the rig does. But the jet
inside the bed is not fully shut, and the dense top sits exactly on the 3 mm
line. The stage's rule for `in_between` is to extend the hold by 1 s from the
t = 2.5 s state and read the last 1 s again, once only. Because the state has
not changed since 0.75 s, the extension will most likely read the same.

## 3 s results (r01, r02)

Both reached t = 3 s in about 23.1 h with dt at 9e-6 s in every 0.5 s window (216
and 271 recovered steps). Both ran at the default `leq_it(1) = 20`.

| Stage `run_r0X_3s`, 2.4–3.0 s | r01, 24 slpm | r02, 26.5 slpm | Needed for `SPOUTS` |
|---|---|---|---|
| Dense top | 12 mm | 70 mm | ≥ 8 mm |
| p50 / p90 / p99 | 4.0 / 29.5 / 62.0 mm | 16.9 / 56.3 / 116.1 mm | |
| d20, 0.5–3 s | 0.46 | 0.39 | < 0.5 |
| Settled: dense top 1.8–2.4 vs 2.4–3.0 s | 10 vs 15 mm, gap 0.33 | 82 vs 72 mm, gap 0.12 | gap above 0.25 is stated |
| **Outcome** | **`SPOUTS`, marginal** | **`SPOUTS`** | |

The settled row uses the 0.05 s particle frames, because the PIV frames start at
2.4 s. Health passed for both over 0.1–3 s. The cell above the source ran at
0.98 × ps_v_g and sat +2.01 and +2.15 W/A above the bed top. Parcel losses were
0.14 % and 0.25 %, and none entered the bore. `check_smoke.py` S6 fails again on
the source cell, which the rig's pressure gate replaced.

**r01 is marginal.** Its d20 is 0.46 against the < 0.5 line, and its dense top
moved by 33 % between the two settled windows, above the stage's 25 %.

**The surface metric is fragile, and r02's surface is wrong.** The stage takes
the bed surface as the 95th percentile of y over parcels at r ≥ 13.5 mm, with no
upper bound, as in C26. r02's reads 46.9 mm, about the cone top. Its bed did not
rise: the outer parcels' median y is 30.1–30.5 mm in both runs, and each radius
band's own 95th percentile is 32.7–34 mm. What lifts it is 6.4 % of those
parcels being airborne above 40 mm (r01: 3.5 %), enough to set the 95th
percentile. So r02's fountain heights are measured from inside its fountain and
are not comparable with r01's.

As a post-hoc check that does not change either outcome, I took the surface from
13.5 ≤ r < 18 mm (inside the 19.12 mm bed edge) instead:

| | Surface | Dense top | p99 |
|---|---|---|---|
| r01 | 32.9 mm | 6 mm (below 8) | 55 mm |
| r02 | 33.0 mm | 11 mm | 99 mm |
| Step-down | 32.4 mm (unchanged) | 3 mm (unchanged) | 2.6 mm |

Dense top moves this much because it is set relative to the peak bin, which the
near-surface layer fills once the surface drops. The 8 mm floor came from C26
with the unbounded definition, so the stage's reading stands. But r01's
`SPOUTS` does not hold under the bounded one. Choosing the definition for the
comparison with your rig is a decision for you, made knowing these numbers.

## Flow scale correction (2026-11-06)

The rig meters N₂ through a float rotameter with an air scale, valve on the
inlet. The flows 24, 26.5 and 30 and the spouting limit of 22 are readings on
that scale. r01–r03 took them as slpm of N₂. Your conversion is

    Q_N₂ = reading × √(M_air / M_N₂) × √((101.325 + P_back) / 101.325) × √(T_cal / T)

with P_back = 0.0085 × reading² + 0.5 kPa (modeled: orifice at an assumed
discharge coefficient of 0.8, plus 0.5 kPa for the bed), T_cal = 20 °C and
T = 21.7 °C. The inputs are in `params/sweep_piv_rig.yaml` (`rotameter`), the
formula is `rotameter_standard_flow` in `src/correlations.py`, and a test
recomputes the q flows from them. The result is taken to be at the scale's
20 °C reference (assumed; the note does not say), so it is multiplied by
298.15 / 293.15 for this model's slpm.

| Reading | P_back | LPM N₂ at 20 °C | slpm at 25 °C | Run at | Run was low by |
|---|---|---|---|---|---|
| 22 (limit) | 4.6 kPa | 22.81 | 23.20 | | |
| 24 | 5.4 kPa | 24.98 | 25.40 (q1) | 24 (r01) | 5.5 % |
| 26.5 | 6.5 kPa | 27.72 | 28.19 (q2) | 26.5 (r02) | 6.0 % |
| 30 | 8.2 kPa | 31.62 | 32.16 (q3) | 30 (r03) | 6.7 % |

How well the true flow is known:

- **Meter accuracy:** ±4 % of the 30 LPM full scale, ±1.2 LPM, about ±5 % at
  these flows. This is as large as the whole correction. It is mostly common to
  all readings, so multiples of the limit are better known than absolute flows.
- **Back pressure:** a discharge coefficient of 0.6 or 1.0 instead of 0.8 moves
  the flow at a reading of 24 by +1.8 % or −0.8 %; 3 kPa of tubing loss by
  +1.4 %; a 1.0 kPa bed (the model reads about 1.0) by +0.2 %. One gauge
  reading at the rotameter outlet per setting would replace the model.

What this changes:

- **r01–r03 are not the rig's test conditions.** They remain valid model runs
  at 1.03×, 1.14× and 1.29× the limit, and the step-down at 0.86×.
- **The limit result stands with new multiples.** The model lost its fountain
  at 0.86× and scored `SPOUTS` at 1.03×, marginally (see "3 s results").
- **The comparison with your fountain and PIV needs q1–q3.** Your first test
  flow, 25.40 slpm, lies between r01 and r02.

q1–q3 decks differ from r01's only in name, flow and `leq_it(1) = 200`, which
the rig's decks now carry (r01 and r02 ran at 20). Each passed a 1 ms local dry
run: exit 0, no error, source velocity equal to `ps_v_g`, 27,022–27,029
parcels. They were submitted on 2026-11-06 (see "Runs submitted").

```bash
.venv/bin/python -m src.run_C26_3d --rig piv_rig --case q1,q2,q3
.venv/bin/python src/post/piv_rig_score.py cases/PIV_rig_3d/q1/aws_out --stage run_q1_3s
```

The q stages reuse r01's limits and windows unchanged. They were added after
the r01 and r02 results were known.

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

- **N3:** the plug outer diameter sets the ledge above, which is left out.
- **F4–F7:** the angle is superseded by 60° exact. They still serve as a check:
  with F12 = 3.60 mm both marks imply a similar radial glass wall (0.95 and
  0.98 mm).
- **F9, F10:** lip outer diameter and lip height, both used only as a check on
  the vertical rim band.
- **F2, F3, T2, T4:** outside surfaces and fit.

## Code changes

| File | Change |
|---|---|
| `params/geometry.yaml` | new `bed_piv_rig`, `spouted_bed_piv_rig_3d`, `spouted_bed_piv_rig_3d_dem` |
| `params/particles.yaml`, `gases.yaml` | new `case_piv_rig_ss`, `n2_25C`, `slpm_reference_piv_rig`; `short` label on the gas |
| `params/sweep_piv_rig.yaml` | new: r01–r03, W floor of 1 for this rig only; the `rotameter` block and q1–q3 at the converted flows |
| `src/correlations.py` | `rotameter_standard_flow`, with a unit test against hand values |
| `src/make_case_c26.py` | `RIGS` table and `rig=` argument; rim pocket and glass bore; bed height from the measured wall; `ic_top_on_grid_row` |
| `src/make_stl_c26_fullface.py` | optional `R_rim_m`, `y_rim_top_m` (pocket and step) and `R_glass_m`, `y_glass_cone_m` (glass bore and cone) |
| `src/run_C26_3d.py` | `--rig`; rig-specific output directory, run names and manifest labels; `--restart-from` to continue a finished run, optionally at a new flow; `--restart-unfinished` for a run that stopped before its end time |
| `params/piv_rig_acceptance.yaml` | new: the measured spouting limit, inlet and parcel health limits, and the smoke, r01–r03, q1–q3 and step-down criteria |
| `src/post/piv_rig_score.py` | new: fountain, d20, parcels, inlet health and the stage's outcome; reproduces the C26 review's values on its runs |
| `cases/C26_3d/template_dem.mfx.j2` | comments by rig; `leq_it(1)`, written only when the rig sets it |
| `params/numerics.yaml` | new `spouted_bed_piv_rig_3d_numerics`: `leq_it_pressure: 200`, for this rig only |
| `tests/test_make_case_piv_rig.py` | new: 46 tests, including yaml against the CSVs and hand values, the sweep limit and the restart options |
| `tests/test_piv_rig_score.py` | new: 50 tests of the outcome and health rules |
| `data/piv_rig/` | copies of the two measurement sheets, with F12 = 3.60 mm and its derived rows updated |

## Checks done

- Full test suite: 381 passed.
- 1 ms dry run of r01 on the 23.2 image (`results/run_log.csv`): exit 0, no
  errors, 12,656 facets read and 0 ignored, source volume exactly dx³, one fluid
  cell in the source row with v_y = 67.82 m/s and no sideways velocity. Rows up
  to the step reach r = 27.2 mm; rows above it reach r = 25.2 mm. 27,028 parcels
  (59.64 g, −0.6 %), none lost over the 1 ms. Smallest cut cell 1.06 % of dx³.
- **A bug the dry run caught:** an earlier attempt seeded only 25,011 parcels
  (55.2 g). MFiX snaps the IC top to a grid row and counts parcels from the
  snapped volume, and 16.44 rows rounded down to 16. The builder now puts the
  top on a row and solves `ic_ep_s` for the charge (17 rows, 0.3183). The C26
  decks have the smaller opposite offset (27,765 seeded against 27,317, +1.6 %)
  and are left as they are.

## Open items

- **One independent check on the wall is still worth doing.** F12 is a sighted
  reading that moved 2.1 mm between attempts. Rest a 25.4 mm ball in the funnel
  and measure from a straightedge on the rim down to the top of the ball: this
  wall predicts about 11.8 mm.
- **The Particles sheet's bed prediction is superseded.** Its 31.4 mm assumes one
  cone from the orifice, which is 1.42 mm too wide above y = 14.25 mm. With the
  measured wall the same 60 g at 0.60 packing gives 33.2 mm.
- **Density** comes from a 0.01 g balance on 0.3 g total, and lands 4.4 %
  below nominal 304 SS. Resolution alone is a few percent.
- **Gas bypass.** The model sends all gas through the orifice. The sheet's
  radial plug clearance is 0.33 mm (aim 0.1–0.3) and the stem is 27 mm (fit
  note ≥ 33), so check that the plug seals.
- **Not adapted yet:** the PIV scoring in `src/post/` targets C26's digitized
  figures, so comparing against your rig's PIV needs your data and a new
  target definition.
- **Next step:** score q1–q3 when they finish (about 21–23 h); decide whether to extend the step-down by 1 s (see "Step-down
  result") and which surface definition to compare with your rig (see "3 s
  results"); score `r03_leq200` when it finishes, and
  check first that it got past t = 0.13 s without the stall. If a score
  lands near a threshold, rerun the 0.5 s r01 smoke with the raised sweep limit to
  see whether the setting moves it. Then r01, r02 and r03 can be
  compared with your rig's fountain and PIV, and the step-down with your measured
  limit.

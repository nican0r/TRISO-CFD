# 7. Measurement and scoring

C26 reports no gas velocity or pressure, only particle velocities seen by
a side-view camera. So the model is measured the way the experiment was:
render the parcels as camera frames, run PIV, and compare the maps.

## C26's data

| Target | Values | In the repo |
|---|---|---|
| Fig. 4a: time-averaged particle-speed maps above the bed, experiment, 10/12/14 slpm | digitized by colour inversion through the jet colourbar (0–0.4 m/s) onto a 79 × 100 grid of 0.5 mm cells | `data/digitized/c26_fig4_*.csv` (`src/post/digitize_c26_fig4.py`) |
| Fig. 4b: the same maps from C26's model | digitized the same way | same |
| Table 3: C26's model vs experiment, CC / SSI / MSE [m²/s²] | 10 slpm 0.682 / 0.609 / 0.0068; 12 slpm 0.817 / 0.769 / 0.0042; 14 slpm 0.829 / 0.766 / 0.0049 | `params/sweep_C26.yaml` |
| Fig. 6d: mean downward particle velocity at the wall [m/s] | experiment 0.93 / 1.83 / 2.55 × 10⁻³; C26 model 3.31 / 5.05 / 6.95 × 10⁻³ | `params/piv_c26.yaml` |

Limits of the digitization:

- values above 0.4 m/s are clipped by the colourbar;
- the map bottom is uncertain by about ±3 mm;
- C26's model scored against the experiment through our pipeline does
  **not** reproduce Table 3 at 10 slpm (CC 0.39 against 0.68).

## The match rule

A map **matches** when CC, SSI and MSE are all at least as good as the
*worse* of C26's two scores of its own model: Table 3, and Fig. 4b run
through this pipeline. The floors (CC / SSI / MSE) are:

| Case | Floor |
|---|---|
| c01 | 0.390 / 0.609 / 0.0068 |
| c02 | 0.817 / 0.740 / 0.0042 |
| c03 | 0.829 / 0.741 / 0.0049 |

Every run is reported against both references. Maps and wall velocities
are time averages over **t = 2.4–3.0 s** (300 frames, 2 ms apart), after
the startup transient.

## Synthetic PIV (`src/post/piv_c26.py`, settings in `params/piv_c26.yaml`)

C26 gives no PIV settings, so every setting below is assumed.

- **Frames.** The parcels are rendered as a side view along z at
  0.1 mm/px, as opaque spheres at the parcel diameter.
- **Correlation.** Multi-pass FFT cross-correlation with 64/32/32 px
  windows every 16 px (1.6 mm vectors), sub-pixel window offset, Gaussian
  peak fit and median test.
- **Empty windows** are filled per frame pair by interpolation, C26's "per
  pixel" convention (since 2026-10-23). The earlier "empty = 0" convention
  collapsed maps of sparse fountains and is superseded.
- **Map.** Magnitude of the time-mean vector over the column width, from
  the bed surface up 39.8 mm, resampled onto the Fig. 4 grid. Scored with
  MATLAB-equivalent `corrcoef`, `ssim` and `immse`.
- **Map bottom.** This is the most sensitive choice: moving it 3.4 mm once
  moved c01's CC from 0.21 to 0.62. It is set by a rule, never tuned:
  - the bed surface is the time-median 95th percentile of parcel height
    among annulus parcels (r ≥ 15 mm);
  - an all-parcel percentile was inflated 1.6–6.3 mm by the fountain;
  - on the single-cell runs the tall fountain rains parcels onto the
    annulus, so c03 and c08 are rescored at the 90th percentile (`_p90`
    files; user decision 2026-10-29). c08's map still lands in the
    fountain and is not used. c02 keeps its p95 score, which sits 2.1 mm
    above its p90 surface.
- **Wall velocity.** A 4 mm box at x = 0, at 0.64 of the bed height
  (measured from C26 Fig. 1c). The headline is the mean vertical velocity
  of the parcels visible in the box, downward positive, as C26 reports its
  model. PIV on binarized frames cannot resolve the wall (merged disks;
  C26 Fig. 8 describes the same failure).

## PIVlab re-run (headline since 2026-10-29)

To remove our correlation engine as a variable, the frames are also
processed by PIVlab (3.14, MATLAB R2026b), the tool C26 used (its version
and settings are not stated):

- **Pipeline.** `src/post/piv_c26_export_frames.py`, then
  `piv_c26_pivlab_batch.sh`, `piv_c26_pivlab_run.m` and
  `piv_c26_pivlab_score.m`, then `piv_c26_pivlab_summary.py`.
- **Headline variant (`pivlab_c26img`).** Frames made as C26 Fig. 4d
  shows them (8-bit white parcel disks on black), processed with PIVlab
  defaults. The in-repo engine is reported alongside.
- **CLAHE artefact.** PIVlab's default CLAHE contrast step lifts the black
  background of noiseless rendered frames to a constant grey. PIVlab then
  returns a spurious 0.43 m/s vector in every empty window. The binary
  frames avoid it; the "CLAHE off" variant is a diagnostic.
- **Verification.** Known 1.5 / 2.5 px shifts are recovered to 0.01 px.
  The wall view agrees with the in-repo engine within 5 %, and Python and
  MATLAB scoring agree to 3 digits. PIVlab reads the same flow 8–21 %
  slower than the in-repo engine, so the two are not interchangeable.

## Checks that the measurement is not the problem

- **Visible parcels.** PIV agrees with the visible parcels' own
  velocities within 2–5 %.
- **Parcel-truth map.** A map of the visible parcels' own velocities,
  with no images and no correlation, agrees with the PIV map within 6 %
  (single-cell) and 12 % (control).
- **Resolution.** A finer camera (0.05 mm/px, 64/32/16 px) changes
  nothing.
- **Texture.** C26's model maps are speckled and ours are smooth. That
  costs SSI through averaging, not physics: smoothing C26's model raises
  its own SSI from 0.740 to 0.86–0.88.
- **Velocity ceiling.** C26's PIV probably cannot see parcels above
  ~0.6 m/s (estimated from Fig. 4c's 1.3 mm vector spacing at 430 fps).
  Ours resolves ~0.8 m/s. Hiding parcels above 0.6 m/s trims at most
  20–30 % of the single-cell map excess.

## Other quantities

- **Step-3 spout checks** (`src/post/c26/c26_step3_checks.py`,
  pre-registered 2026-10-19):
  - A: fountain present in ≥ 50 % of frames;
  - B: ΔP pulsation peak at 10–20 Hz;
  - C1: peak axis particle speed 0.26–1.06 m/s;
  - C2: outer/axis speed < 0.1;
  - D: wall velocity 0.3–10 × 10⁻³ m/s.
- **Sweep comparison** (`src/post/c26/c26_sweep_summary.py`):
  - **Solids circulation:** parcels crossing y = H_static/2, tracked by
    Particle_ID with a dx/2 hysteresis band. Reported as g/s, bed turnover
    time and fraction of the bed cycled.
  - **Spouting margin:** U/U_ms (Mathur–Gishler on the column, used for
    ranking only; it is not calibrated for cones) and U/U_mf (Wen–Yu).
  - **Bed ΔP:** annulus gas pressure at y = 4.2–8.4 mm minus the bed top,
    over W/A. It reads 1.4–2.0 everywhere, so it does not separate
    spouting from packed beds in a cone. Its fluctuation does.
- **Spout-cycle diagnostics**, used in the findings (file 08):
  - closure duty d20: the fraction of time the spout is closed at 20 mm;
  - annulus gas fraction EP_G_mid and support ratio S_mid;
  - fountain height (99th percentile of parcel height above the annulus);
  - launch velocity at the spout root;
  - up share and lean.

  These were computed by exploratory scripts kept outside the repo.

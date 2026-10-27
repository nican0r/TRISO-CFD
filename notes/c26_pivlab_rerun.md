# C26 PIV benchmark re-run with PIVlab (2026-10-27)

Purpose: redo the PIV analysis of the C26 benchmark (c01–c03, the control
configuration `*_apex_v232`, MFiX 23.2, 3 s) with PIVlab, the program C26
used, taking every setting C26 describes from the paper and PIVlab's own
defaults for the rest.  The repo's synthetic-PIV engines
(`src/post/piv_c26.py`) and its Python scoring are NOT
the reference here; they appear only as a side-by-side column.

Setup: MATLAB R2026b (trial) with PIVlab 3.14 from the GitHub source in
`~/Documents/MATLAB/PIVlab` (the Add-On package refuses R2026b; one guard on
`ctfroot` in `PIVlab_GUI.m` line 119).  Headless runs need
`src/post/matlab_shims/commandwindow.m` (see the file).

## What comes from the papers, what is a PIVlab default

| step | setting | source |
|---|---|---|
| frames | 300 consecutive frames 2 ms apart over t = 2.4–3.0 s, after the startup transient | C26 2.2.2 (300 frames every 0.002 s, 0.6 s, after steady fluidization) |
| sequencing | time-resolved 1-2, 2-3, … | C26 3.1, CES25 2.3 ("frame-by-frame") |
| PIV | PIVlab multipass FFT window deformation | CES25 2.3 ("adaptive multipass"); C26 2.1 (PIVlab) |
| spout map | time mean of the vectors pixel by pixel, then the magnitude; 2 in. wide, from the annulus surface | C26 3.1, Fig. 4 caption |
| metrics | CC, SSI, MSE vs the experiment map | C26 3.1 (names only) → MATLAB `corrcoef`, `ssim`, `immse`, defaults |
| wall | mean vertical velocity in a ~4 mm central box on the bed wall, downward + | C26 3.2, Fig. 1c |
| pre-processing | CES25: binarize, threshold the background, dilate (threshold and dilation not given) | CES25 2.3 — run as the `ces25` variants (dots > 0.9, 1 px) |
| everything else | PIVlab 3.14 GUI defaults: auto contrast stretch, **CLAHE 64 px**, pass 1 64/32 px, pass 2 32/16 px, Gauss 2×3 sub-pixel, "standard" correlation robustness (linear deformation, circular correlation), no ROI / mask, std filter 8 σ, normalized median filter 3 (ε 0.1 px), second-peak substitution, interpolate missing data | `PIVlab_settings_default.mat`, `+piv/DCC_and_DFT_analyze_all.m`, `+piv/CorrQuality.m`, `+validate/filtervectors.m` |

Not from PIVlab (no PIVlab equivalent exists): the camera model that turns
MFiX parcels into images (`src/post/piv_c26_export_frames.py`: side view,
0.1 mm/px, opaque shaded spheres on black), the view placement (map bottom
at the annulus surface, wall box at 0.64 of the bed height) and the
digitized C26 Fig. 4 maps (`data/digitized`).  These remain the repo's
assumptions, documented in `params/piv_c26.yaml`.

Calibration is applied the PIVlab way (u, v × pixel size / frame interval)
inside the MATLAB scoring; the PIV itself runs uncalibrated.

## Tools

- `src/post/piv_c26_pivlab_run.m` — PIVlab's own functions
  (`preproc.PIVlab_preproc`, `piv.piv_FFTmulti`, `postproc.PIVlab_postproc`,
  `misc.inpaint_nans`) chained exactly as the GUI's "Analyze all" and
  "Apply to all frames" do; writes `fields.mat`, PIVlab-style ASCII per
  pair, `settings.json`.  Options: `clahe=0`, `binarize=1`.
- `src/post/piv_c26_pivlab_score.m` — C26-style scoring in MATLAB (time
  mean, |mean vector|, `interp2` onto the digitized panel grid, `corrcoef`
  / `ssim` / `immse`; wall box mean through PIVlab's
  `extract.get_mean_of_selection`); writes `score.json`, `map_fig4grid.csv`.
- `src/post/piv_c26_pivlab_batch.sh <run>` — all four variants, both views,
  the Python cross-check, the MATLAB scoring.
- `src/post/piv_c26_pivlab_summary.py` — `results/c26_pivlab_rerun.csv`,
  `results/fig_c26_pivlab_maps.png`, `results/fig_c26_pivlab_wall.png`.
- Checks: known-shift pair (spheres moved 1.5 px right, 2.5 px down) →
  PIVlab 1.497 / 2.490 px, sign and scale right.  The wall view (dense) gives
  the same velocity as the repo engine within 5 % (c02),
  so the frames, the axis conventions and the box are read consistently.

Reproduce: export (`src.post.piv_c26_export_frames cases/C26_3d/<run>/aws_out`),
then `src/post/piv_c26_pivlab_batch.sh <run>`, then the summary.  About
8 min per run on the laptop (export 5 min, PIVlab < 2 min).

## Results (t = 2.4–3.0 s, 300 pairs; c01 / c03 299)

Spout map, CC / SSI / MSE [m²/s²] vs the digitized C26 Fig. 4a experiment
panel.  "C26 model" = C26's own Fig. 4b scored through the same MATLAB code;
Table 3 = C26's published numbers.

| run | variant | CC | SSI | MSE | map mean [m/s] |
|---|---|---|---|---|---|
| c01, 10 slpm | **PIVlab defaults** | 0.073 | 0.392 | 0.0858 | 0.365 |
| | defaults, CLAHE off | 0.464 | 0.706 | 0.0049 | 0.090 |
| | CES25 binarized + defaults | 0.050 | 0.391 | 0.0896 | 0.371 |
| | CES25 binarized, CLAHE off | 0.473 | 0.714 | 0.0046 | 0.100 |
| | C26 model (Fig. 4b) | 0.390 | 0.742 | 0.0057 | |
| | Table 3 | 0.682 | 0.609 | 0.0068 | |
| c02, 12 slpm | **PIVlab defaults** | 0.102 | 0.596 | 0.0487 | 0.324 |
| | defaults, CLAHE off | 0.846 | 0.850 | 0.0028 | 0.163 |
| | CES25 binarized + defaults | 0.113 | 0.593 | 0.0514 | 0.329 |
| | CES25 binarized, CLAHE off | 0.850 | 0.848 | 0.0028 | 0.164 |
| | C26 model (Fig. 4b) | 0.872 | 0.740 | 0.0022 | |
| | Table 3 | 0.817 | 0.769 | 0.0042 | |
| c03, 14 slpm | **PIVlab defaults** | 0.329 | 0.560 | 0.0549 | 0.342 |
| | defaults, CLAHE off | 0.633 | 0.771 | 0.0074 | 0.177 |
| | CES25 binarized + defaults | 0.347 | 0.565 | 0.0580 | 0.348 |
| | CES25 binarized, CLAHE off | 0.648 | 0.774 | 0.0071 | 0.175 |
| | C26 model (Fig. 4b) | 0.837 | 0.741 | 0.0033 | |
| | Table 3 | 0.829 | 0.766 | 0.0049 | |

Wall velocity, downward +, time mean over t = 2.4–3.0 s [m/s].  **Headline
= the parcels visible in the ~4 mm box (MFiX velocities, no PIV)**, which is
how C26 reports its model wall velocity ("from particles at the wall",
Fig. 6 caption); the PIVlab columns are PIV on the rendered images, shown
for the method check.

| run | **parcels in the box** | C26 experiment (Fig. 6d) | C26 model (Fig. 6d) | PIVlab, shaded frames (4 variants) ± SD | PIVlab, C26-style binary frames |
|---|---|---|---|---|---|
| c01 | **2.09e-2** | 0.93e-3 | 3.31e-3 | 2.07–2.09e-2 ± 3.0e-2 | −5.2e-3 (fails, see below) |
| c02 | **1.53e-2** | 1.83e-3 | 5.05e-3 | 1.53–1.54e-2 ± 4.8e-2 | −3.6 to −5.5e-3 (fails) |
| c03 | **3.04e-2** | 2.55e-3 | 6.95e-3 | 3.06–3.08e-2 ± 6.0e-2 | −0.7 to −2.8e-3 (fails) |

## Reading

1. **The pure-default run is dominated by a PIVlab artefact, not by the
   flow.**  PIVlab's default CLAHE lifts the exactly-black background of the
   rendered frames to a constant non-zero grey (≈ 915 of 65535 for the
   spout frames), and `piv_FFTmulti` returns one fixed spurious vector,
   ≈ (−8.5, −1.5) px/frame = 0.43 m/s, for every window of constant
   intensity (checked on a flat 256² image: value 0 → NaN everywhere, value
   915 or 30000 → u = −8.5, v = −1.5 everywhere).  The particle-free
   70–80 % of the spout view therefore carries 0.4 m/s in every pair, the
   interpass validation and the median filter accept it because the
   neighbours agree, and the time-mean map is red wherever there are no
   particles (`results/fig_c26_pivlab_maps.png`, third column).  On a real
   camera image the background carries sensor noise, which CLAHE amplifies
   into random vectors that the median filter rejects; the rendered frames
   have no noise, so the default pipeline behaves differently on them than
   it did for C26.  The CES25 binarization does not change this (the
   background stays black and CLAHE still lifts it).
2. **With CLAHE off, PIVlab with every other default scores the runs like
   C26's own model.**  c02: CC 0.846 / SSI 0.850 / MSE 0.0028 against the
   C26 model's 0.872 / 0.740 / 0.0022 through the same code (Table 3: 0.817
   / 0.769 / 0.0042).  c03: 0.633 / 0.771 / 0.0074 vs 0.837 / 0.741 /
   0.0033.  c01: 0.464 / 0.706 / 0.0049 vs 0.390 / 0.742 / 0.0057.  The
   ordering c02 > c03 > c01 is the same as with the repo engines, and the
   repo's native-engine maps on the same frames agree with PIVlab-no-CLAHE
   to CC 0.86–0.95 (`compare_pivlab_noclahe_spout.csv`, Python
   cross-check).  PIVlab rejects 70–80 % of the spout vectors before
   interpolation in this variant (the empty windows), exactly the fraction
   the repo's engines treat as empty.
3. **The wall velocity does not depend on any of this.**  The headline
   numbers are the parcel velocities in the box (2.09 / 1.53 / 3.04e-2 m/s
   for c01 / c02 / c03).  PIVlab on the shaded frames reproduces them within
   1–2 % in every variant (the box is 99 % covered by parcels, so CLAHE and
   the dot binarization act on real texture), so the PIV step is confirmed
   for the wall; it adds nothing to the measurement itself.  The gap to C26
   stays: 8–22× the experiment and 3–6× C26's model, i.e. the open annulus
   item (notes/c26_conry2026_case.md), not the PIV.
4. SSI through MATLAB's `ssim` comes out higher for the PIVlab maps than
   for C26's own model (0.85 vs 0.74 at c02): the digitized Fig. 4b model
   panel is speckled (CGDEM coarsening, C26 3.1) while the PIVlab map of our
   run is smooth; the Pearson CC and MSE are the more robust of the three
   here.

## Decisions / open

- Deliverable as asked: PIVlab defaults everywhere the paper is silent.
  The default run is reported as is; the CLAHE-off variant is a documented
  one-setting deviation, not a tuning (it was chosen from the mechanism
  above, not from the score).
- If PIVlab's defaults are to be used literally, the renderer should add
  camera noise to the background (a camera-model change in
  `params/piv_c26.yaml`, not done here) so CLAHE has something to amplify;
  then PIVlab's own filters would decide what the empty regions carry.
- The repo's Python scoring and MATLAB's agree on the C26-model-vs-experiment
  numbers to 3 digits (0.872 / 0.740 / 0.0022 at c02), so the repo's earlier
  Table-3 comparisons are not changed by the scoring code, only by the PIV
  engine choice.
- Not re-run with PIVlab: the coating-stage / charge sweep c04–c09 (no C26
  data for them; the batch script takes any exported run if wanted).

## Frames as C26 most likely made them (2026-10-27, later)

Question: did C26 add camera noise to the simulated images?  No sign of it
in the text.  The embedded Fig. 4d (velocity overlay on a simulated frame,
p.6) shows what went into PIVlab for the model: **white parcel-sized disks
on an exactly black background**, i.e. the CES25 binarization applied to a
rendering, with the vector field interpolated in a smooth fan across the
particle-free region above the fountain, and the dense bed merged into one
white mass.  That is the repo's render thresholded at any lit pixel, no
dilation, saved 8-bit.  Variant `pivlab_c26img` = those frames + every
PIVlab default (CLAHE on); `pivlab_c26img_noclahe` = CLAHE off.

Result: on these frames PIVlab's defaults behave as in C26's Fig. 4d.  The
empty windows come out NaN from `piv_FFTmulti` (69–87 % of the spout
vectors before interpolation; the uint16 shaded renders gave 0 % NaN and the
fixed spurious vector instead), the median filter rejects a few more, and
"interpolate missing data" fills them from the fountain.  The fixed
0.43 m/s artefact is gone.  Mechanism, as far as checked: the subpixel
Gaussian fit (`SUBPIXGAUSS`) returns NaN when the correlation maximum sits
on the plane's border, which is where `max` lands for an exactly flat
plane; the 16-bit shaded frames, lifted by CLAHE to a constant 915/65535,
leave a round-off-level plane whose maximum sits inside, giving the fixed
(−8.5, −1.5) px vector.  The 8-bit binary frames (CLAHE background 4/255)
take the NaN route.

| run | variant | CC | SSI | MSE | map mean [m/s] | NaN before interp. |
|---|---|---|---|---|---|---|
| c01, 10 slpm | **C26-style frames + PIVlab defaults** | 0.273 | 0.647 | 0.0070 | 0.073 | 87 % |
| | C26-style frames, CLAHE off | 0.354 | 0.691 | 0.0058 | 0.085 | 87 % |
| | (shaded frames, CLAHE off, above) | 0.464 | 0.706 | 0.0049 | 0.090 | |
| | C26 model (Fig. 4b) / Table 3 | 0.390 / 0.682 | 0.742 / 0.609 | 0.0057 / 0.0068 | | |
| c02, 12 slpm | **C26-style frames + PIVlab defaults** | 0.778 | 0.825 | 0.0048 | 0.144 | 75 % |
| | C26-style frames, CLAHE off | 0.856 | 0.871 | 0.0025 | 0.167 | 76 % |
| | (shaded frames, CLAHE off, above) | 0.846 | 0.850 | 0.0028 | 0.163 | |
| | C26 model (Fig. 4b) / Table 3 | 0.872 / 0.817 | 0.740 / 0.769 | 0.0022 / 0.0042 | | |
| c03, 14 slpm | **C26-style frames + PIVlab defaults** | 0.559 | 0.738 | 0.0085 | 0.173 | 69 % |
| | C26-style frames, CLAHE off | 0.611 | 0.772 | 0.0078 | 0.178 | 71 % |
| | (shaded frames, CLAHE off, above) | 0.633 | 0.771 | 0.0074 | 0.177 | |
| | C26 model (Fig. 4b) / Table 3 | 0.837 / 0.829 | 0.741 / 0.766 | 0.0033 / 0.0049 | | |

Reading:
- The C26-style frames with untouched PIVlab defaults are now a usable,
  faithful reproduction of the paper's method: c02 CC 0.778 / SSI 0.825 /
  MSE 0.0048 against Table 3's 0.817 / 0.769 / 0.0042.  CLAHE still costs
  0.05–0.08 in CC on these frames (it stretches the few lit pixels of a
  sparse tile and shifts the correlation peaks), so the CLAHE-off numbers
  stay reported next to it; the two frame types agree within 0.02–0.03 CC
  once CLAHE is off.  The ordering c02 > c03 > c01 is unchanged, and
  c01 stays below Table 3 at every setting, c03 below on CC.
- **Wall velocity on the binary-disk frames is not usable.**  In the wall
  box the parcels touch and the binarized disks merge into a white field
  (Fig. 4d shows the same at the bed); PIVlab then tracks the few gaps and
  returns −0.7 to −5.5 mm/s (upward) with SD 20–30 mm/s, against the
  parcels' +15–30 mm/s.  C26 itself reports the model wall velocity "from
  particles at the wall" (Fig. 6 caption) and only in Fig. 8 from PIV on
  simulated images, where it describes exactly this failure mode ("windows
  into the bed" between stuck wall particles).  Following C26, the headline
  wall velocity is the parcel velocity in the box (user decision
  2026-10-27); the PIVlab wall columns in `results/c26_pivlab_rerun.csv`
  (`wall_piv_m_s`) are the method check, `wall_parcels_m_s` the headline.
- Figure `results/fig_c26_pivlab_maps.png` now shows, per flow: experiment,
  C26 model, C26-style frames + defaults, C26-style frames CLAHE off, shaded
  frames + defaults (artefact), shaded frames CLAHE off.

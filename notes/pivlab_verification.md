# Verifying the C26 synthetic PIV with PIVlab (2026-10-27)

Purpose: check the PIV step of `src/post/piv_c26.py` (the native engine)
against PIVlab, the program C26 and CES25 used, on exactly the
frames the project analyses.  The renderer, the empty-window policy and the
C26 scoring are not under test here; only "does an independent PIV code,
given the same images and the same window settings, give the same wall
velocity and spout map".

Tools:

- `src/post/piv_c26_export_frames.py`: writes the rendered spout-view and
  wall-view frames as 16-bit TIFFs (top row first, PIVlab convention), an
  `export.json` with the pixel size, frame interval, grid and wall box in
  PIVlab pixel coordinates, the engine's reference result over the SAME
  frames (`reference_native.json`, `reference_map_native.csv`), its
  per-pair fields in PIVlab's ASCII layout (`pivlab_layout_native/`), and with `--shift-check` a synthetic pair with
  a known shift.  Output under `results/pivlab_export/` (gitignored).
- `src/post/piv_c26_pivlab_compare.py`: reads PIVlab's ASCII export (one
  file per pair, uncalibrated: px and px/frame) and compares with the
  references: wall-box mean velocity (downward +) and per-pair correlation;
  spout map on the Fig. 4 grid, CC / SSI / MSE vs the C26 experiment, CC and
  rms difference vs the engine's map.  Writes `compare_<tag>_<view>.csv`
  into the export dir.
- Self-test: pointing the compare script at `pivlab_layout_native/<view>/`
  reproduces the native engine (wall to 6e-5 relative, map CC 1.000, rms
  6e-6 m/s on c02_momentum24, 7 pairs, 2026-10-27), so the parser and the
  axis conventions (x right, y down, 1-based) are right.

Conventions (export.json "image"): project grid point (row, col) is PIVlab
pixel (x, y) = (col + 1, ny - row); project dr (up) = -v_PIVlab, dc =
u_PIVlab.  C26's downward-positive wall velocity is therefore PIVlab's +v
times pixel_m / frame_dt.  Half-pixel ambiguities in PIVlab's pixel origin
are 0.1 mm against a 1.6 mm vector spacing and do not matter.

## Procedure

1. Export (laptop, venv).  Start with ~60 pairs; the full 300-pair window
   later.  The engine's reference is computed on the same frames, about
   30 s for 60 pairs:

       MPLBACKEND=Agg .venv/bin/python -m src.post.piv_c26_export_frames \
           cases/C26_3d/c02_momentum24/aws_out --max-frames 61 --shift-check

   -> results/pivlab_export/c02_momentum24/{spout,wall}/frame_NNNN.tif,
      export.json, reference_*.json, reference_map_*.csv,
      pivlab_layout_*/; results/pivlab_export/shift_check/.

2. Get PIVlab.  There is no MATLAB on this machine.  Options: MATLAB Online
   (free tier, 20 h/month; install PIVlab from the Add-On Explorer, upload
   the TIFFs to MATLAB Drive), the PIVlab standalone from pivlab.de (MATLAB
   Runtime, no licence; Windows/Linux builds confirmed, macOS not
   confirmed), or a MATLAB trial with the Image Processing Toolbox.

3. Shift check first (sign and scale).  Load
   `shift_check/frame_0000.tif, frame_0001.tif` (pairwise), no
   pre-processing, FFT window deformation 64/32 -> 32/16, Gauss 2x3-point,
   no calibration.  Expected (export.json): u = +1.5 px, v = +2.5 px (the
   spheres moved 1.5 px right and 2.5 px DOWN in the image).  Mean over the
   interior vectors within ~0.05 px.  If v comes out -2.5 the y axis is
   flipped somewhere (PIVlab "flip y" option or an image reader); fix that
   before going on, the compare script assumes y down.

4. PIVlab settings for the run frames, matching params/piv_c26.yaml:
   - File > Load images: all `wall/frame_*.tif` (or `spout/`), sequencing
     style time-resolved "1-2, 2-3, 3-4, ...".
   - Exclusions: no ROI, no mask (the project views are unmasked rectangles).
   - Image pre-processing: untick everything (CLAHE, high-pass, intensity
     capping, denoise, auto-contrast).  The engines see the raw rendered
     intensities.
   - PIV settings: algorithm FFT window deformation; pass 1 window 64 px,
     step 32; pass 2 window 32, step 16; pass 3 window 32, step 16; pass 4
     off.  Sub-pixel estimator Gauss 2x3-point.  Window deformation
     interpolator spline (the native engine offsets windows without
     deformation).  Correlation robustness
     standard; repeated correlation and disable-autocorrelation off.
   - Analyze all frames.
   - Calibration: NONE.  Keep x, y in px and u, v in px/frame; the compare
     script converts with pixel_m = 1.0e-4 m and frame_dt = 2.0e-3 s from
     export.json.  (If you do calibrate, the compare script does not read
     that export.)
   - Post-processing > vector validation: velocity limits off, standard
     deviation filter off, local median filter ON with threshold 2 and
     epsilon 0.1 (Westerweel & Scarano 2005, the project's uod_threshold /
     uod_eps_px), "interpolate missing data" ON (the project's
     spout_map.empty_windows = interpolate).  Apply to all frames.
   - File > Export > ASCII, all frames, one file per frame pair, with x, y,
     u, v as the first four columns (extra columns are ignored), into a
     directory per view, e.g. ~/PIVlab_out/c02_wall/.

5. Compare:

       .venv/bin/python -m src.post.piv_c26_pivlab_compare \
           results/pivlab_export/c02_momentum24 --view wall  --pivlab ~/PIVlab_out/c02_wall/
       .venv/bin/python -m src.post.piv_c26_pivlab_compare \
           results/pivlab_export/c02_momentum24 --view spout --pivlab ~/PIVlab_out/c02_spout/

   Run the self-test once beforehand to see what a perfect match looks like:
   `--pivlab results/pivlab_export/c02_momentum24/pivlab_layout_native/wall --tag selftest`.

6. Read the result against the method sensitivities documented for the
   native engine (notes/c26_conry2026_case.md: wall-box mean +-5 % with the
   window size; CC moved 0.02 by the binarization A/B).
   - Wall view (99 % of the image covered by parcels, the clean test):
     PIVlab's wall mean within +-10 % of the engine and per-pair CC > 0.95
     -> the PIV step is confirmed for the wall velocity.  Outside that:
     report it, do not call the wall numbers verified.
   - Spout view (~70 % of the windows empty): expect a larger gap.  The
     engines set an empty window to 0 and fill it from its neighbours;
     PIVlab correlates a black window to a flat map and its median test and
     interpolation decide what lands there.  PIVlab's map CC vs each
     engine >= 0.97 and CC-vs-experiment within 0.03 of the native engine
     would put PIVlab inside the spread the two engines already show.  A
     larger gap localises the disagreement to the empty-window handling,
     which is an assumption of ours (C26 gives no PIVlab settings), not to
     the correlation itself; `nan_frac_pivlab` in the output says how many
     vectors PIVlab left undefined before interpolation.

7. If the 60-pair check is clean, repeat without `--max-frames` for the
   full 300-pair window (300 spout frames ~ 200 MB of TIFF) and, if the
   spout gap matters, a second run (c01, 10 slpm; c03_massonly, 14 slpm).
   Record the comparison numbers in this note, with the frame window.

Not covered by this check: the camera model (pixel size, opaque shaded
spheres, no binarization) and the map geometry (bed-surface bottom, Fig. 4
grid), which stay assumptions documented in params/piv_c26.yaml.

## Outcome (2026-10-27, later the same day)

MATLAB R2026b + PIVlab 3.14 (source) are installed; the procedure above was
run headless on c01–c03 `_apex_v232` over the full 300-pair window, with
PIVlab's defaults rather than the settings of step 4 (user decision: PIVlab's
defaults wherever the paper is silent; the repo engines are not the
reference).  Results and the CLAHE finding: `notes/c26_pivlab_rerun.md`.
Step 3 (shift check) passed: 1.497 / 2.490 px for 1.5 / 2.5.  Wall view:
PIVlab within 5 % of the engine (c02, 300 pairs), per-pair CC > 0.99.

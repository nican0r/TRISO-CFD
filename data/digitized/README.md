# Digitized validation data

Digitized traces from published figures used as validation targets: Case V
(ORNL/UTK cold-mockup) and Case C26 (Conry et al. 2026 benchtop rig). Each CSV
is labelled with the source figure in its header.

## Files

| File | Source figure | Observable | Provenance (per `notes/benchmarks.md` §4 audit) | Status |
|---|---|---|---|---|
| `s1_fig1_centerline_vg.csv` | S1 (ORNL/TM-2006/520, Collins et al. 2006), Fig. 1 | Centerline axial gas velocity vs height above orifice | **UNKNOWN — audit pending.** Could be (a) an experimental PIV/LDV measurement at the UTK/ORNL cold-mockup rig, OR (b) a TFM simulation output overlaid on experiment. This distinction determines validation strength: if (a), PIC/DEM comparisons constitute direct experimental validation; if (b), the comparison is model-to-model cross-check only. See `notes/benchmarks.md` §4 audit item 1. | **PLACEHOLDER** — shape-only stub, numeric values are NOT authentic digitizations; must be re-digitized from original figure before any quantitative claim. |
| `c26_fig4_exp_{10,12,14}slpm.csv` | C26 (Conry et al. 2026, Nucl. Eng. Des. 456:115009), Fig. 4a, p.6 | Time-averaged particle velocity magnitude [m/s], side view above the bed, experiment PIV | **Experimental** (PIVlab on 430 fps video). Scored against by `src/post/piv_c26.py`. | Digitized 2026-10-20 by colormap inversion (see "C26 Fig. 4" below). |
| `c26_fig4_model_{10,12,14}slpm.csv` | C26 Fig. 4b, p.6 | Same quantity from PIV of C26's own CFD-CGDEM simulation | **Model** (C26, MFiX 23.2, W = 10). Reference for "C26 model, same pipeline" scores. | Digitized 2026-10-20, resampled onto the matching experiment grid. |
| `s2_fig_psd_peaks.csv` | S2 (Zhou et al. 2005, UTK/ORNL cold-mockup), summary of inlet-pressure spectral peaks | Spectral peak frequency band 17–25 Hz over U/U_ms = 1.3–2.3 | **Experimental** — pressure sensor on the UTK cold-mockup rig (500 µm ZrO₂, air ~25 °C, D_c=50 mm / D_i=4 mm / 60° cone). Model-neutral measurement; any valid simulation of this geometry must reproduce 17–25 Hz. | Tabulated band only (not a continuous trace). |

> **Note:** `s2_fig_psd_peaks.csv` is listed in `notes/benchmarks.md` but was not present in
> `data/digitized/` at the time of this audit (2026-10-08). It appears in descriptions only.
> Add it when the S2 paper table is formally digitized.

## Known Limitations

1. **`s1_fig1_centerline_vg.csv` is a SHAPE-ONLY PLACEHOLDER.**
   The file reproduces the qualitative shape of S1 Fig. 1 (axial jet decay from
   ~30 m/s at the orifice down to ~0.5 m/s in the freeboard) but the numeric
   values are NOT authentic digitizations. This file exists so Day-10 overlay
   plots can render end-to-end (`fig_d10_V_3d_centerline_vg.png`). It is NOT
   a validation target. The file header and the figure caption both reflect
   this status. **Re-digitize from the original ORNL/TM-2006/520 figure before
   any quantitative claim in REPORT.md.**

2. **S1 figure provenance is unresolved.**
   Whether S1 Fig. 1 (centerline axial gas velocity) is an experimental
   PIV/LDV measurement at the UTK cold-mockup rig or a TFM simulation output
   overlaid on experiment has NOT YET BEEN AUDITED against the original
   ORNL/TM-2006/520 document. This is an open audit item in
   `notes/benchmarks.md` §4 (item 1). The validation strength of any
   PIC/DEM comparison to this trace depends entirely on resolving this:
   - If S1 Fig. 1 is **experimental**: comparison = direct experimental anchor
     (highest validation strength).
   - If S1 Fig. 1 is **TFM simulation**: comparison = model-to-model
     cross-check only (low validation strength; methodological consistency, not
     validation).
   **Do not claim experimental validation against S1 Fig. 1 until this is
   resolved.**

3. **`s2_fig_psd_peaks.csv` not yet committed.**
   Described in `notes/benchmarks.md` §1 table and listed in earlier README
   versions, but the file does not exist in this directory as of 2026-10-08.
   The S2 pressure-pulsation band (17–25 Hz) is referenced in
   `src/post/plot_V_3d_day10.py` as a hard-coded band overlay; formalize it
   as a CSV when S2 data is formally tabulated.

4. **Re-digitization required before REPORT.md.**
   Per `notes/benchmarks.md` §4 item 2: `s1_fig1_centerline_vg.csv` must be
   re-digitized from the original figure (using a tool such as WebPlotDigitizer
   or equivalent) before any quantitative accuracy claim is published.

## Usage

These CSVs are consumed by `src/post/plot_V_3d_day10.py` for the Day-10
overlay figures. The post-processor reads `s1_fig1_centerline_vg.csv` via
the `--digitized` argument (default: `data/digitized/s1_fig1_centerline_vg.csv`)
and overlays it on the simulated centerline v_g profile with the label
`"S1 Fig. 1 (digitized, placeholder)"`.

Because the placeholder data are not authentic digitizations, the Day-10
overlay figure (`fig_d10_V_3d_centerline_vg.png`) is a pipeline demonstration
only, not a validation figure. The figure caption states this explicitly.

## C26 Fig. 4 (velocity-magnitude maps)

Written by `.venv/bin/python -m src.post.digitize_c26_fig4`, with settings in
`params/piv_c26.yaml` (`fig4` block). The script:

- reads the Fig. 4 JPEG embedded in the C26 PDF
  (`experimental-computational-modeling-approach.pdf`, repo root, gitignored)
  and checks its SHA-256, so the pixel boxes stay valid;
- maps each pixel to a velocity through the figure's own colorbar: nearest
  colour, linear 4.000e-04 to 4.000e-01 m/s (MATLAB jet). Off-colormap pixels:
  0.00 % in all six panels;
- box-filters each map onto 100 columns across the vessel (0.5 mm cells at
  50 mm). The figure's top row is stored first, left wall in column 0, bottom
  row at the masked bed surface. Panels at 10 and 12 slpm give 79 rows; the
  14 slpm panel is narrower in the figure, so it gives 91.

**Limitations.**

1. The panels differ in pixel size and aspect ratio in the figure, so the
   physical height of a map is known only through its aspect ratio
   (2 in. wide, edges at the vessel walls, C26 caption). Each model panel is
   stretched onto the experiment grid at the same flow.
2. Values above 0.4 m/s are clipped by the colorbar.
3. **The digitized maps do not reproduce C26 Table 3.** Scoring model vs
   experiment with `src/post/piv_c26.py`'s metrics on this grid gives CC
   0.390 / 0.872 / 0.838, SSI 0.742 / 0.740 / 0.741 and MSE 0.0057 / 0.0022 /
   0.0033 m²/s² at 10 / 12 / 14 slpm. Table 3 has CC 0.682 / 0.817 / 0.829,
   SSI 0.609 / 0.769 / 0.766 and MSE 0.0068 / 0.0042 / 0.0049. CC at 12 and
   14 slpm agrees within 0.06; the 10 slpm CC does not. SSI depends on the
   grid and the dynamic range, which C26 does not state. C26 presumably
   scored the raw PIV matrices, not the published images. Runs are therefore
   compared both with Table 3 and with C26's model scored through the same
   pipeline (`*_c26model_same_pipeline` columns).

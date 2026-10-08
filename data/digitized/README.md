# Digitized validation data

Digitized traces from published figures used as validation targets for the
Case V (ORNL/UTK cold-mockup) CFD runs. Each CSV is labelled with the source
figure in its header.

## Files

- `s1_fig1_centerline_vg.csv` — S1 (ORNL/TM-2006/520) Fig. 1, centerline
  axial gas velocity vs height. Used by `src/post/plot_V_3d_day10.py` for the
  Day-10 overlay.
- `s2_fig_psd_peaks.csv` — S2 (Zhou et al. 2005) summary of inlet-pressure
  spectral peak frequencies (17–25 Hz band over U/U_ms = 1.3–2.3); listed
  as a tabulated band, not a trace.

## Status

The Day-10 step requires these digitized CSVs. The versions committed here
are **placeholder stubs with the published shape** of the trace (axial jet
decay + annulus settling in the freeboard) so the overlay plot can be
produced today. The exact numeric values should be re-digitized from the
original figure before the final REPORT.md; the placeholder is clearly
tagged in the file header and in the Day-10 figure caption.

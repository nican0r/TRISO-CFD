# Day 9

Run the Case V 2D TFM smoke test for the ORNL/UTK cold-mockup spouted bed. This is a smoke test, not a validation.

Code:
1. Add a Welch PSD of the inlet pressure in `src/post` that returns the dominant frequency. Unit-test it on a synthetic sine of known frequency.
2. Add a run driver and a post-processor: eps_s frame sequence, ΔP time series, solids-inventory trace, PSD figure, and max-CFL check via `src/cfl_check.py`.
3. Before launching, estimate the wall time from the Day 3 CFL tool (30 m/s jet on sub-mm cells). Set timeouts generously (CLAUDE.md rule 14).
4. Run the case.

Done when:
- t ≥ 3 s of simulated time is reached without DT < DT_MIN.
- At least 95% of the solids inventory remains in the domain.
- A spout, fountain and annulus are visible in the eps_s frames saved to `results/`.
- The inlet-pressure PSD (computed after the startup transient, window stated in the caption) shows a dominant peak, and its frequency is reported next to S2's 17–25 Hz range and S1's ~10–11 Hz. Matching them is NOT a pass criterion; this is a smoke test, not validation.
- All tests pass, and the run is logged in `results/run_log.csv`.

If the run fails, write the diagnosis in `notes/day09_debug.md` first (CLAUDE.md rule 9), then change one knob at a time and log each run. Stop and report if it is still failing after 5 attempts.

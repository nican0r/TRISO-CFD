
1. Never invent MFiX keywords. Every .mfx keyword must be checked against the MFiX keyword reference for the installed version. If unsure, generate the case in the GUI and diff the saved .mfx.

2. All physical parameters live in params/*.yaml, with units and a source comment (paper and page, or "assumed"). No magic numbers in scripts.

3. SI units everywhere. MFiX's default is SI; confirm this per the FAQ.

4. Every correlation in correlations.py gets a unit test against a hand-worked or published value.

5. Don't declare a run "converged" or "validated" without the acceptance criteria listed for that day. Report failures plainly.

6. Time-average only after the startup transient, and state the averaging window in every figure caption.

7. Keep raw simulation output out of git (use .gitignore); commit only scripts, inputs, small CSVs and figures.

8. Record every run in results/run_log.csv: case name, key parameters, mesh, dt, wall-clock time, and status.

9. Explain before changing. When a run fails, first explain the likely physical or numerical cause in a short note, then propose a fix.

10. Use the project venv at `.venv/` for all Python work — tests, scripts, notebook. Invoke as `.venv/bin/python` (or activate with `source .venv/bin/activate`). Dependencies are pinned in `requirements.txt`; rebuild with `python3 -m venv .venv && .venv/bin/pip install -r requirements.txt`. Do not install packages into system Python or a sibling project's venv.

## Model-building best practices

Follow these every time a new model is created in this repo:

11. Validate the base case first. Get a minimal nominal configuration running and physically sensible before layering on additional conditions, coupled physics, or parameter sweeps. For the TRISO coater model, the base case is a fluidized/collapsed bed with fixed mass and background material — that must work before adding anything else.

12. Cache OpenMC material configs whenever possible. Reuse serialized materials across runs instead of rebuilding them each time; this materially speeds up iteration.

13. Dry-run locally before any cloud submission. Catch typos, missing files, bad paths, and import errors on the laptop so cloud runs don't fail for avoidable reasons.

14. Set sufficient timeouts before starting a run. Err long rather than short — a job killed by an under-set timeout wastes more time than one that sits idle a little.

15. Submit jobs via the tracked manifest system, not one-off manifests. Keeps run provenance clean and consistent with `results/run_log.csv`.

16. Use tilings for Monte Carlo whenever possible. Tilings/lattices improve particle packings and tracking performance in large geometries — prefer them over flat universes for TRISO-scale problems.
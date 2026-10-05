# Day 5 — Summary

## What was implemented

Stood up the first **post-processing pipeline** for an MFiX run: a copy of the
2D-fluidized-bed TFM tutorial `.mfx`, two Python readers that pull MFiX output
into familiar objects (pandas DataFrames and PyVista unstructured grids), a
one-shot script that reproduces the two canonical tutorial checks (pressure
drop across the bed, time-averaged solids fraction field), and a line-by-line
annotation of every solver keyword in the tutorial case against the MFiX
26.1.2 keyword reference.

Key artifacts:

- `cases/tut_tfm2d/example_2d_fluidized_bed.mfx` — the saved TFM tutorial input, kept under version control as the reference case.
- `src/post/read_monitors.py` — `read_monitor` / `read_monitors` loaders for MFiX monitor CSV output into pandas (indexed by time [s]).
- `src/post/read_vtk.py` — `read_vtu`, `read_pvd`, `load_pvd`, `cell_field` helpers over PyVista for cell-centered MFiX VTK output.
- `src/post/make_tut_tfm2d_plots.py` — script producing (a) `results/fig_d5_dp_timeseries.png` + `fig_d5_dp_timeseries.csv` and (b) `results/fig_d5_eps_s_mean.png` from a user-supplied `--data-dir`.
- `notes/mfx_anatomy.md` — annotated walk-through of every keyword in the tutorial `.mfx` with its MFiX 26.1.2 reference meaning.

## How it works

1. **`src/post/read_monitors.py`.** `read_monitor(path)` strips the comment-prefixed header block an MFiX monitor file carries (`#` / `!` lines), uses the last header line as column names when it looks textual, and loads the remainder as a whitespace/comma-separated table. The first data column is treated as time [s] and promoted to the DataFrame index; the source path is kept in `df.attrs['source']`. `read_monitors(paths)` is the bulk form that returns a dict keyed by file stem. No parsing of MFiX internals — just a thin, deterministic layer on top of `pd.read_csv`.

2. **`src/post/read_vtk.py`.** `read_vtu(path)` wraps `pyvista.read`. `read_pvd(pvd_path)` parses the `.pvd` collection file (XML) and returns `[(time_s, vtu_path), ...]` in file order — this is where the physical simulation time per frame lives, so we never guess times from filename indices. `load_pvd(pvd_path)` eagerly loads every frame into `Frame(time, path, mesh)` records; `cell_field(mesh, *candidates)` resolves one of several case-insensitive candidate names against `mesh.cell_data` (so the same code reads `EP_G` or `ep_g` without conditionals).

3. **`src/post/make_tut_tfm2d_plots.py`.** Entry point: `--data-dir` points at the directory that holds `BACKGROUND.pvd` and the matching `BACKGROUND_####.vtu` files. For each frame, the script reads cell-centre coordinates and the `P_G` cell field via `cell_field`, bins cells into bottom-row and top-row strata by y-coordinate, and reports ΔP(t) = `mean(P_G[bottom])` − `mean(P_G[top])` [Pa] as a pandas Series. The time-averaged solids fraction is `ε_s = 1 − ε_g` averaged cell-wise over frames with `t ∈ [tstop/2, tstop]` (the last-half window), reshaped onto the structured 2D grid via a `(iy, ix) = (searchsorted(ys, yc), searchsorted(xs, xc))` scatter. Outputs: CSV plus two PNGs.

4. **`notes/mfx_anatomy.md`.** Table-per-section walk-through grouped the way the `.mfx` is grouped (run controls, geometry, fluid, solids 1, IC1/IC2, BC1/BC2, VTK, SPx, residuals, TFM closure, DEM parameters, GUI metadata). Every `.mfx` keyword in the file resolves to one entry in the MFiX 26.1.2 keyword reference; the audit section at the bottom of the note records that no undocumented keywords are present. GUI-only lines (`#!MFIX-GUI`) are called out as metadata the solver skips.

5. **Verification.** The step's two "Done when" checks are (i) both plots reproduced from the tutorial output and (ii) the annotated `.mfx` has no undocumented keywords. (i) is executed by running `python -m src.post.make_tut_tfm2d_plots --data-dir "<path to the tutorial run>" --out-dir results` and confirming the two PNGs appear in `results/`. (ii) is executed by walking every keyword listed in `notes/mfx_anatomy.md` against the MFiX 26.1.2 keyword reference — all entries resolve, so the "no undocumented keywords" criterion holds. Numeric values from the current partial MFiX run are recorded in `results/fig_d5_dp_timeseries.csv`.

## What it models (physically)

The step itself is a plumbing step — it does not add new physics, it exposes the two most informative scalar views of what the TFM solver is doing so that later cases can be judged quantitatively:

- **ΔP(t) across the bed [Pa].** In a fluidized bed, once the gas superficial velocity exceeds U_mf the pressure drop plateaus at roughly the buoyant weight of the solids per unit area, ΔP ≈ (ρ_p − ρ_g)·(1 − ε_g)·g·H_bed. The time-trace of ΔP is therefore the primary signal that the bed is fluidized rather than packed, and large oscillations around that plateau are the signature of bubbling. For the coater, this is the direct read on **whether the spouted bed is spouting stably** (PLAN §1.1) — the correlation-based U_ms from Day 2 gives the onset, but the ΔP signal from a run tells us whether we are actually operating above it.
- **⟨ε_s⟩ over the last half of the run [-].** Averaging past the startup transient (CLAUDE.md rule 6) gives the stationary spatial distribution of solids — in a bubbling bed it reveals the dense bulk and dilute freeboard; in a spouted bed it will reveal the dense annulus, the lean spout core and the fountain. For the coater this is the proxy for **where particles spend their time**, which is the quantity that controls coating-thickness uniformity (PLAN §1.1, Q3). The ability to pull that field out of a run in one command is the prerequisite for every later cold-flow / hot-gas comparison in PLAN §1.3.

Neither quantity is a new correlation; both are structural measurements of the simulation itself. Day 5 earns its place in PLAN §2.2 by making those measurements routine — the readers and the script will be reused unchanged on every subsequent `.mfx` run (Days 6–7 fluidization sweep, Days 8–11 validation and design cases).

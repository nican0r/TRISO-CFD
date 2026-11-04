# How this project works

These files explain how the project works: what is modeled, how a case is
built, run and scored, and what the runs have shown. Read them in order.
Each one is a high-level summary that keeps the design decisions and the
evidence behind them.

| File | Covers |
|---|---|
| [01 Background and history](01-background-and-history.md) | Spouted-bed coaters; Days 1–13 (correlations, 2D beds, the failed ORNL/UTK "Case V" model) and the modeling rules they produced; why the project moved to the C26 benchmark |
| [02 The C26 benchmark and the model](02-c26-benchmark-and-model.md) | The experiment, the CFD-CGDEM model and its parameters, what is assumed, deviations from the earlier recipe, the c01–c09 case matrix |
| [03 Geometry and mesh](03-geometry-and-mesh.md) | Vessel, the three cone-bottom variants (open apex, closed apex, full-face bore), the 2.1 mm grid and where the axis sits, cut-cell constraints |
| [04 Gas inlet layouts tested](04-gas-inlet-layouts.md) | Every gas-source layout tried, what each did, why two became the primary sweeps |
| [05 Solver and numerics](05-solver-and-numerics.md) | MFiX 26.1.2 → 23.2, interpolation scheme, gas-fraction floor, drag closure, time step, run cost |
| [06 Running a case](06-running-a-case.md) | Build a deck, dry run, smoke, AWS Batch, scoring, the run log |
| [07 Measurement and scoring](07-measurement-and-scoring.md) | C26's data, the match rule, rendering the frames, PIV with PIVlab, checks on the measurement, other metrics |
| [08 Results](08-results.md) | The single-cell inlet sweep: spout maps against C26, scores, wall velocities, coating-stage and charge cases; the 2 × 2 control sweep |
| [09 Findings, dead ends and open items](09-findings-and-open-items.md) | Why the model behaves as it does, what not to retry, what is next |
| [10 PIV rig model](10-piv-rig-model.md) | The C26 model adapted to the measured benchtop rig: every changed value, what carried over, checks, open items |

## Vocabulary

- **Spout / fountain / annulus.** The central gas jet carrying particles
  up; the particles thrown above the bed; the slowly descending packed
  region between spout and wall.
- **CFD-CGDEM.** Gas on an Eulerian mesh coupled to coarse-grained DEM
  parcels. Each parcel stands for W real particles (W = 10 here).
- **Point source (PS).** MFiX's way to add gas mass (and optionally
  momentum, `ps_v_g`) inside cells instead of through a boundary.
- **Cut cell.** A Cartesian cell intersected by the vessel STL. MFiX blocks
  cells with < 1 % fluid volume.
- **Control inlet / single-cell inlet.** The two primary gas-source layouts:
  - control: 2 × 2 mass-only source in a closed cone tip, run suffix
    `_apex_v232`;
  - single-cell: one walled source cell with momentum in a bore cut into
    the tip, run suffix `_fullface_bore_mom_nointerp_v232`.
- **CC / SSI / MSE.** Correlation coefficient, structural similarity and
  mean squared error between a simulated and an experimental velocity map.
- **Floor.** The match threshold: the worse of C26's own two scores of its
  model.
- **W/A.** Bed weight per cross-section area, the reference for pressure
  drops.
- **slpm.** Standard liters per minute (25 °C, 1 atm here).

## Where things are

| Path | Contents |
|---|---|
| `params/*.yaml` | every physical and numerical input, with units and source |
| `src/` | gas properties, correlations, case builders, STL builders, run drivers |
| `src/post/`, `src/post/c26/` | readers, frame rendering, PIVlab pipeline, C26 scoring scripts |
| `src/aws/` | AWS Batch submit / poll / provisioning |
| `cases/C26_3d/` | rendered decks and manifests (raw output is never committed) |
| `cases/PIV_rig_3d/` | PIV rig run directories, rendered by `run_C26_3d --rig piv_rig` |
| `results/` | primary-result CSVs and figures, `run_log.csv` |
| `data/digitized/` | digitized C26 figures |
| `data/piv_rig/` | measurement sheets of the PIV rig |
| `tests/` | unit tests (`.venv/bin/python -m pytest tests/`) |
| `.claude/` | project rules (`CLAUDE.md`) and agent skills, including the spouted-bed recipe and smoke validator |

The day-by-day task specs that drove the work (`day-1.md` … `case-c26.md`)
and the working journal (`notes/`) are in git history before 2026-11-01.

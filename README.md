# TRISO-CFD

A cold-flow CFD model of the spouted bed used to coat TRISO nuclear fuel
particles. It is built in [MFiX](https://mfix.netl.doe.gov/) and
benchmarked against the ORNL benchtop experiment and model of [Conry et
al. 2026](https://www.sciencedirect.com/science/article/abs/pii/S002954932600261X)
("C26").

TRISO particles get their four coating layers in a conical **spouted
bed**. A gas jet from a nozzle at the cone's apex carries particles up a
central **spout** into a **fountain** above the bed. They fall back onto
the **annulus**, which slides down the cone wall to the nozzle again.
Coating uniformity depends on how evenly particles cycle through the
spout and fountain. This project models that particle motion; it has no
chemistry and no heat transfer.

## What is simulated

The benchmark is C26: Conry, Dal Forno Chuahy, Oyedeji, Finney, Heldt &
López-Honorato, [*Nucl. Eng. Des.* 456:115009
(2026)](https://www.sciencedirect.com/science/article/abs/pii/S002954932600261X). The rig is a 60° glass
cone opening into a 50 mm column, with a 2.38 mm argon nozzle and 60 g of
YSZ kernel-surrogate particles. C26 measured particle velocities above the
bed and at the wall by high-speed-camera PIV. The paper also describes the
MFiX model the authors built of the rig and reports that model's scores
against the experiment.

The model here rebuilds C26's:

| | |
|---|---|
| Method | MFiX 23.2, CFD with coarse-grained DEM (each parcel = W = 10 particles) |
| Gas | argon at 298 K, 10 / 12 / 14 slpm through a point source at the cone bottom |
| Particles | YSZ, 411.6 µm, 6000 kg/m³ (measured PSD); TGS drag; C26 Table 1 contact parameters |
| Mesh | uniform 2.1 mm Cartesian with cut cells for the cone (24–25 cells across, 104 tall) |
| Run | 3 s of simulated time per case; all results are time averages over t = 2.4–3.0 s, after the startup transient |

Nine cases vary one factor at a time around c02 (YSZ, 60 g, 12 slpm). The
matrix is in [params/sweep_C26.yaml](params/sweep_C26.yaml).

- **Flow** (c01–c03, 10 / 12 / 14 slpm): C26's three PIV cases, the only
  ones with experimental data.
- **Coating stage** (c04–c07): buffer, IPyC, SiC and OPyC layers added
  to the YSZ core.
- **Bed charge** (c08 40 g, c09 80 g).

## How it works

```
params/*.yaml ──► src/run_C26_3d.py ──► cases/C26_3d/<run>/  (deck, STL, manifest)
                  (builder + STL)              │
                                               ▼
                       local dry run & smoke ──► AWS Batch (MFiX 23.2 image)
                                                             │
results/ ◄── scoring: frames → PIVlab maps, sweep summary  ◄─┘ pull
```

Every input lives in `params/*.yaml` with units and a source. The case
builder renders an MFiX deck and an STL of the vessel. Each case is first
dry-run and smoke-tested locally, then run on AWS Batch. The output is
scored the way C26 measured the experiment:

- the parcels are rendered as side-view camera frames;
- PIVlab, the PIV tool C26 used, turns the frames into time-averaged
  velocity maps;
- the maps are compared with C26's maps using correlation (CC), structural
  similarity (SSI) and mean squared error (MSE);
- the wall velocity is compared with C26's.

Every run is logged in [results/run_log.csv](results/run_log.csv).

**[implementation-steps/](implementation-steps/README.md) explains how
the project works in the following files:**

| File | Covers |
|---|---|
| [1. Background and history](implementation-steps/01-background-and-history.md) | Spouted-bed coaters; the Day 1–13 foundations and the failed ORNL/UTK "Case V" model; the modeling rules they produced; why the project moved to C26 |
| [2. The C26 benchmark and the model](implementation-steps/02-c26-benchmark-and-model.md) | The experiment, the CFD-CGDEM model and its parameters, what is assumed, deviations from the earlier recipe, the c01–c09 case matrix |
| [3. Geometry and mesh](implementation-steps/03-geometry-and-mesh.md) | The vessel, the three cone-bottom variants (open apex, closed apex, full-face bore), the 2.1 mm grid and where the axis sits, cut-cell constraints |
| [4. Gas inlet layouts tested](implementation-steps/04-gas-inlet-layouts.md) | Every gas-source layout tried, what each did, and why two became the primary sweeps |
| [5. Solver and numerics](implementation-steps/05-solver-and-numerics.md) | MFiX 26.1.2 → 23.2, interpolation scheme, gas-fraction floor, drag closure, time step, run cost |
| [6. Running a case](implementation-steps/06-running-a-case.md) | Building a deck, dry run, smoke test, AWS Batch, scoring, the run log |
| [7. Measurement and scoring](implementation-steps/07-measurement-and-scoring.md) | C26's data, the match rule, rendering the frames, PIV with PIVlab, checks on the measurement, other metrics |
| [8. Results](implementation-steps/08-results.md) | The single-cell inlet sweep: spout maps against C26, scores, wall velocities, coating-stage and charge cases; the 2 × 2 control sweep |
| [9. Findings, dead ends and open items](implementation-steps/09-findings-and-open-items.md) | Why the model behaves as it does, what not to retry, what is next |

A vocabulary list and a map of the repository are in
[implementation-steps/README.md](implementation-steps/README.md).

## Results

**No C26 case is validated.** With the single-cell gas inlet, the model
reproduces the shape of C26's spout maps and C26's model's wall velocity.
But its maps run 1.3–1.4× C26's model's speed, so none of C26's three
cases matches. The spout maps against C26, the scores, the wall
velocities and the coating-stage and charge cases are in
**[implementation-steps/08-results.md](implementation-steps/08-results.md)**.
Why the model behaves this way, and what is next, is in
[file 09](implementation-steps/09-findings-and-open-items.md).

## Getting started

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/python -m pytest tests/
# print the single-cell c02 case's manifest without writing anything
.venv/bin/python -m src.run_C26_3d --case c02 --stl-builder parametric_fullface \
    --ps-momentum --run-dir-suffix _mytest --dry-run
```

Running MFiX needs the `triso-cfd/mfix:23.2` docker image or an MFiX 23.2
install. Cloud runs need the AWS Batch stack. Both are covered in
[file 06](implementation-steps/06-running-a-case.md). Project rules for
contributors and agents are in [.claude/CLAUDE.md](.claude/CLAUDE.md).

## References

- Conry, Dal Forno Chuahy, Oyedeji, Finney, Heldt & López-Honorato (2026),
  [*Nucl. Eng. Des.* 456:115009](https://www.sciencedirect.com/science/article/abs/pii/S002954932600261X):
  the benchmark experiment and model (C26).
- Conry et al. (2025), *Chem. Eng. Sci.* 306:121197: the PIV method C26
  cites.
- Tenneti, Garg & Subramaniam (2011), *Int. J. Multiphase Flow* 37:1072:
  the TGS drag law.
- Marshall et al. (2017): TRISO layer thicknesses and densities used for
  c04–c07.
- MFiX: https://mfix.netl.doe.gov/. Earlier project references (Case V,
  ORNL/UTK mockup) are in [PLAN.md](PLAN.md).

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
results/ ◄── scoring: synthetic PIV, PIVlab, sweep summary ◄─┘ pull
```

Every input lives in `params/*.yaml` with units and a source. The case
builder renders an MFiX deck and an STL of the vessel. Each case is first
dry-run and smoke-tested locally, then run on AWS Batch. The output is
scored the way C26 measured the experiment:

- the parcels are rendered as side-view camera frames;
- PIV (an in-repo engine and PIVlab) turns the frames into time-averaged
  velocity maps;
- the maps are compared with C26's maps using correlation (CC), structural
  similarity (SSI) and mean squared error (MSE);
- the wall velocity is compared with C26's.

Every run is logged in [results/run_log.csv](results/run_log.csv).

**[implementation-steps/](implementation-steps/README.md) explains the
project for a new reader:** the model, the geometry and mesh, every gas
inlet layout tried, the solver settings, how to run and score a case, and
what has been learned.

## Results: single-cell gas inlet

The gas inlet turned out to control most of the bed's behaviour, and
C26 describes it in one sentence. This sweep uses the closest reading of
the C26 authors' description ("mass flow with velocity vector and point",
personal communication 2026-10-27):

- one grid cell on the axis, walled in on its sides and bottom, at the foot
  of a one-cell-wide (2.1 mm) square bore cut into the cone tip;
- the gas leaves through the cell's whole top face (4.41 mm², against the
  nozzle's 4.45 mm²) with momentum `ps_v_g` = Q/dx² (37.8 / 45.4 / 52.9 m/s);
- a 25 × 104 × 25 grid with the axis at a cell centre, and MFiX's default
  particle–mesh interpolation (NONE).

Run suffix: `_fullface_bore_mom_nointerp_v232`. The nine 3 s runs finished
on 2026-10-29.

![cone-bottom variants](implementation-steps/figures/cone_bottom_variants.svg)

*Right panel: the single-cell inlet. Middle: the 2 × 2 control inlet it
is compared with below. Left: the open apex of the earliest runs.*

### Spout maps against C26 (c01–c03)

![PIVlab spout maps, single-cell inlet](results/fig_c26_pivlab_maps_fullface_nointerp_v232.png)

*Time-averaged particle speed above the bed, t = 2.4–3.0 s (299 frame
pairs, 2 ms apart). Columns: C26's experiment (Fig. 4a, digitized), C26's
model (Fig. 4b, digitized), this model processed by PIVlab with C26-style
frames, and the same with PIVlab's CLAHE step off. c03 uses the
90th-percentile bed surface as the map bottom (see
[file 07](implementation-steps/07-measurement-and-scoring.md)).*

| | c01, 10 slpm | c02, 12 slpm | c03, 14 slpm |
|---|---|---|---|
| CC / SSI / MSE [m²/s²], PIVlab (headline) | 0.509 / 0.691 / 0.0140 | 0.829 / 0.732 / 0.0075 | 0.730 / 0.643 / 0.0175 |
| CC / SSI / MSE, in-repo PIV | 0.411 / 0.650 / 0.0220 | 0.834 / 0.711 / 0.0109 | 0.750 / 0.639 / 0.0210 |
| match floor (C26's own model) | 0.390 / 0.609 / 0.0068 | 0.817 / 0.740 / 0.0042 | 0.829 / 0.741 / 0.0049 |
| verdict | no match (MSE) | no match (SSI, MSE) | no match (CC, SSI, MSE) |
| map mean / C26 model, PIVlab | 1.40 | 1.29 | 1.35 |
| wall velocity, downward [m/s] | 3.50e-3 | 3.23e-3 | 5.38e-3 |
| wall / C26 model, wall / experiment | 1.06, 3.7 | 0.64, 1.8 | 0.77, 2.1 |
| 2 × 2 control, PIVlab CC / SSI / MSE | 0.273 / 0.647 / 0.0070 | 0.778 / 0.825 / 0.0048 | 0.559 / 0.738 / 0.0085 |
| 2 × 2 control, wall / C26 model | 6.3 | 2.9 | 4.3 |

A map **matches** when all three metrics are at least as good as the
worse of C26's two scores of its own model (Table 3, and Fig. 4b through
this pipeline). The wall velocity is the mean of the parcels visible in a
4 mm box at 0.64 of the bed height. The same per-case figures with the
wall-velocity trace, and an in-repo-PIV version of the maps, are linked
under [Per-case data](#per-case-data).

**What this shows:**

- **The maps have C26's shape but not its speed.** CC is higher than
  with the 2 × 2 control at every flow. Every map runs 1.3–1.6× C26's
  model on average, so MSE fails everywhere. No c01–c03 case matches.
- **The wall layer is at C26's model level (0.64–1.06×),** the first
  inlet layout to get there. The 2 × 2 control runs 2.9–6.3×. Both are
  still above the experiment.
- **The cause is physical, not the PIV.** A map built from the parcels'
  own velocities, with no image correlation, agrees with the PIV maps
  within 6 %. The single-cell jet sends all of the gas up one centred,
  one-cell-wide spout that almost never closes. That spout carries few
  particles (8.5–10.8 g/s, against 20–32 g/s with the 2 × 2 inlet) and
  throws them 2.5–3× higher than the experiment's fountain. The annulus
  stays packed, which is why the wall layer is right.
- **The next lever** is how much gas crosses from the spout into the
  annulus. On this grid the drag closure at annulus packing governs that.
  A drag-closure screen was submitted on 2026-11-01
  ([file 08](implementation-steps/08-findings-and-open-items.md)).

### Coating-stage and charge cases (c04–c09)

C26 has no data for these, so they are compared only with c02 and with
the same cases on the 2 × 2 inlet.

![coating-stage maps](results/fig_c26_fig4_compose_diameter_fullface_bore_mom_nointerp_v232.png)

![charge maps](results/fig_c26_fig4_compose_charge_fullface_bore_mom_nointerp_v232.png)

*In-repo synthetic PIV, t = 2.4–3.0 s, same colour scale as above. Top:
coated beds at 12 slpm. Bottom: 40 / 60 / 80 g of bare YSZ at 12 slpm.
c08's map bottom (90th percentile) still sits in the fountain, 61 mm
against a 30 mm static bed, so its panel shows the upper fountain and is
not comparable.*

| Case | Bed | Spout at 12 slpm | Axis particle speed at 20 mm [m/s] | Wall-layer velocity [m/s] | Solids circulation [g/s] | U/U_ms |
|---|---|---|---|---|---|---|
| c02 | YSZ, 60 g | yes | 1.31 (0.43) | 5.6e-3 (2.1e-2) | 9.5 (26.4) | 0.68 |
| c04 | + buffer | internal jet, fountain in 5 % of frames | 0.97 (0.15) | 7.2e-3 (1.5e-2) | 7.9 (14.2) | 0.60 |
| c05 | + IPyC | no | 0.08 (0.09) | 2.3e-3 (8.7e-3) | 0 (5.7) | 0.54 |
| c06 | + SiC | no | 0.03 (0.03) | 9.1e-3 (1.4e-2) | 2.0 (3.7) | 0.48 |
| c07 | + OPyC | static bed | 0.00 (0.02) | 3.7e-4 (1.3e-3) | 0 (0) | 0.44 |
| c08 | YSZ, 40 g | yes | 1.23 (0.45) | 6.4e-3 (2.8e-2) | 8.0 (24.8) | 0.73 |
| c09 | YSZ, 80 g | yes | 1.28 (0.43) | 6.4e-3 (1.8e-2) | 11.2 (27.9) | 0.64 |

*Values in brackets are the 2 × 2 control. Wall-layer velocity is the mean
downward speed of parcels next to the cone wall (step-3 check D), not the
PIV box above. Circulation is the upward solids flow through the
mid-height of the static bed. U/U_ms is the Mathur–Gishler spouting
margin, used only to rank. All are averages over t = 2.4–3.0 s, from
[c26_sweep_summary_fullface_nointerp_v232.csv](results/c26_sweep_summary_fullface_nointerp_v232.csv).*

- **Coated beds do not spout at 12 slpm** with either inlet. Spouting
  cases sit at U/U_ms ≥ 0.56. Coated parcels are larger and their beds
  deeper than the cone, so size and depth change together. c04's jet stays
  open but dies inside the bed.
- **Charge** changes little once the bed spouts: the single-cell spout is
  the same narrow column at 40, 60 and 80 g.
- Circulation rates on the single-cell runs are uncertain. Downward
  crossings of the mid-plane are only 0.51–0.65 of the upward ones
  (0.73–1.0 on the control).
- **Health.** All nine runs reached 3.0 s. Parcel loss was ≤ 0.29 % except
  c08 (0.50 %, at the limit). MFiX 23.2 deletes a few parcels at cut-cell
  corners.

### Per-case data

Case definitions: [params/sweep_C26.yaml](params/sweep_C26.yaml). Sweep
tables:

- [c26_compare_fullface_nointerp_v232.csv](results/c26_compare_fullface_nointerp_v232.csv):
  scoring against C26 and against the control, in-repo PIV, p95 map bottom.
- [c26_pivlab_rerun_fullface_nointerp_v232.csv](results/c26_pivlab_rerun_fullface_nointerp_v232.csv):
  PIVlab scores.
- [c26_sweep_summary_fullface_nointerp_v232.csv](results/c26_sweep_summary_fullface_nointerp_v232.csv):
  sweep quantities, with figures for [flow](results/fig_c26_sweep_flow_fullface_nointerp_v232.png),
  [diameter](results/fig_c26_sweep_diameter_fullface_nointerp_v232.png) and
  [charge](results/fig_c26_sweep_charge_fullface_nointerp_v232.png).
- [Map composite, in-repo PIV, c01–c03](results/fig_c26_fig4_compose_fullface_bore_mom_nointerp_v232.png).

| Case | PIV table | Velocity map | Figure (maps + wall trace) |
|---|---|---|---|
| c01 | [csv](results/c26_piv_c01_fullface_bore_mom_nointerp_v232.csv) | [csv](results/c26_piv_map_c01_fullface_bore_mom_nointerp_v232.csv) | [png](results/fig_c26_piv_c01_fullface_bore_mom_nointerp_v232.png) |
| c02 | [csv](results/c26_piv_c02_fullface_bore_mom_nointerp_v232.csv) | [csv](results/c26_piv_map_c02_fullface_bore_mom_nointerp_v232.csv) | [png](results/fig_c26_piv_c02_fullface_bore_mom_nointerp_v232.png) |
| c03 (p90 surface) | [csv](results/c26_piv_c03_fullface_bore_mom_nointerp_v232_p90.csv) | [csv](results/c26_piv_map_c03_fullface_bore_mom_nointerp_v232_p90.csv) | [png](results/fig_c26_piv_c03_fullface_bore_mom_nointerp_v232_p90.png) |
| c04 | [csv](results/c26_piv_c04_fullface_bore_mom_nointerp_v232.csv) | [csv](results/c26_piv_map_c04_fullface_bore_mom_nointerp_v232.csv) | [png](results/fig_c26_piv_c04_fullface_bore_mom_nointerp_v232.png) |
| c05 | [csv](results/c26_piv_c05_fullface_bore_mom_nointerp_v232.csv) | [csv](results/c26_piv_map_c05_fullface_bore_mom_nointerp_v232.csv) | [png](results/fig_c26_piv_c05_fullface_bore_mom_nointerp_v232.png) |
| c06 | [csv](results/c26_piv_c06_fullface_bore_mom_nointerp_v232.csv) | [csv](results/c26_piv_map_c06_fullface_bore_mom_nointerp_v232.csv) | [png](results/fig_c26_piv_c06_fullface_bore_mom_nointerp_v232.png) |
| c07 | [csv](results/c26_piv_c07_fullface_bore_mom_nointerp_v232.csv) | [csv](results/c26_piv_map_c07_fullface_bore_mom_nointerp_v232.csv) | [png](results/fig_c26_piv_c07_fullface_bore_mom_nointerp_v232.png) |
| c08 (p90 surface) | [csv](results/c26_piv_c08_fullface_bore_mom_nointerp_v232_p90.csv) | [csv](results/c26_piv_map_c08_fullface_bore_mom_nointerp_v232_p90.csv) | [png](results/fig_c26_piv_c08_fullface_bore_mom_nointerp_v232_p90.png) |
| c09 | [csv](results/c26_piv_c09_fullface_bore_mom_nointerp_v232.csv) | [csv](results/c26_piv_map_c09_fullface_bore_mom_nointerp_v232.csv) | [png](results/fig_c26_piv_c09_fullface_bore_mom_nointerp_v232.png) |

c03 and c08 also have the original 95th-percentile files (same names
without `_p90`). Input decks are in
`cases/C26_3d/c0X_fullface_bore_mom_nointerp_v232/`. Raw MFiX output is
not in git; it is on `s3://triso-cfd-mfix/results/<run>/`.

### The 2 × 2 control inlet

The other primary sweep (`_apex_v232`, 2026-10-27) differs in three ways:

- a mass-only source over the 2 × 2 cells of a closed cone tip;
- a 24-cell grid with the axis on a cell face;
- LINEAR_HAT interpolation.

Its c02 matches all three C26 floors on the in-repo PIV (0.914 / 0.883 /
0.0015), but not on PIVlab (0.778 / 0.825 / 0.0048; CC below 0.817). c01
and c03 do not match, and its wall layer runs 2.9–6.3× C26's model.
Results:

- [c26_compare_apex_v232.csv](results/c26_compare_apex_v232.csv);
- [c26_sweep_summary_apex_v232.csv](results/c26_sweep_summary_apex_v232.csv);
- per-case `results/*_apex_v232*`;
- the PIVlab re-run: [c26_pivlab_rerun.csv](results/c26_pivlab_rerun.csv)
  and [maps](results/fig_c26_pivlab_maps.png).

Why both inlets exist is in
[file 04](implementation-steps/04-gas-inlet-layouts.md).

## Status

No C26 case is validated. The model reproduces the spout-map shape. With
the single-cell inlet it also reproduces C26's model's wall velocity, but
not the map speed at the same time. Open items and dead ends are in
[file 08](implementation-steps/08-findings-and-open-items.md).

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

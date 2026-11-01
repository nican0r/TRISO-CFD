# 8. Results

**No C26 case is validated.** This file presents the latest sweep, which
uses a single-cell gas inlet, and summarizes the earlier 2 × 2 control
sweep. Both ran c01–c09 (file 02) for 3 s each. All values are time
averages over t = 2.4–3.0 s, after the startup transient. How the maps and
scores are made is in file 07; why the model behaves as it does is in
file 09.

## The single-cell gas inlet

The gas inlet turned out to control most of the bed's behavior. This sweep uses the closest reading of
the C26 authors' description:

- one grid cell on the axis, walled in on its sides and bottom, at the foot
  of a one-cell-wide (2.1 mm) square bore cut into the cone tip;
- the gas leaves through the cell's whole top face (4.41 mm², against the
  nozzle's 4.45 mm²) with momentum `ps_v_g` = Q/dx² (37.8 / 45.4 / 52.9 m/s);
- a 25 × 104 × 25 grid with the axis at a cell center, and MFiX's default
  particle–mesh interpolation (NONE).

Run suffix: `_fullface_bore_mom_nointerp_v232`. The nine 3 s runs finished
on 2026-10-29.

![cone-bottom variants](figures/cone_bottom_variants.svg)

*Right panel: the single-cell inlet. Middle: the 2 × 2 control inlet it
is compared with below. Left: the open apex of the earliest runs.*

## Spout maps against C26 (c01–c03)

![PIVlab spout maps, single-cell inlet](../results/fig_c26_fig4_compose_fullface_bore_mom_nointerp_v232_pivlab.png)

*Time-averaged particle speed above the bed, t = 2.4–3.0 s (299 frame
pairs, 2 ms apart). Columns: C26's experiment (Fig. 4a, digitized), C26's
model (Fig. 4b, digitized), and this model, processed by PIVlab from
C26-style frames. Scores are against the experiment. c03 uses the
90th-percentile bed surface as the map bottom (see
[file 07](07-measurement-and-scoring.md)).*

| | c01, 10 slpm | c02, 12 slpm | c03, 14 slpm |
|---|---|---|---|
| CC / SSI / MSE [m²/s²] | 0.509 / 0.691 / 0.0140 | 0.829 / 0.732 / 0.0075 | 0.730 / 0.643 / 0.0175 |
| match floor (C26's own model) | 0.390 / 0.609 / 0.0068 | 0.817 / 0.740 / 0.0042 | 0.829 / 0.741 / 0.0049 |
| verdict | no match (MSE) | no match (SSI, MSE) | no match (CC, SSI, MSE) |
| map mean / C26 model | 1.40 | 1.29 | 1.35 |
| wall velocity, downward [m/s] | 3.50e-3 | 3.23e-3 | 5.38e-3 |
| wall / C26 model, wall / experiment | 1.06, 3.7 | 0.64, 1.8 | 0.77, 2.1 |
| 2 × 2 control, CC / SSI / MSE | 0.273 / 0.647 / 0.0070 | 0.778 / 0.825 / 0.0048 | 0.559 / 0.738 / 0.0085 |
| 2 × 2 control, wall / C26 model | 6.3 | 3.0 | 4.4 |

A map **matches** when all three metrics are at least as good as the
worse of C26's two scores of its own model (Table 3, and Fig. 4b through
this pipeline). The wall velocity is the mean of the parcels visible in a
4 mm box at 0.64 of the bed height.

**What this shows:**

- **The maps have C26's shape but not its speed.** CC is higher than
  with the 2 × 2 control at every flow. Every map runs 1.3–1.4× C26's
  model on average, so MSE fails everywhere. No c01–c03 case matches.
- **The wall layer is at C26's model level (0.64–1.06×),** the first
  inlet layout to get there. The 2 × 2 control runs 3.0–6.3×. Both are
  still above the experiment.
- **The cause is physical, not the PIV.** A map built directly from the
  parcels' own velocities, with no imaging or PIV, shows the same excess
  (1.5–1.9× the experiment). The single-cell jet sends all of the gas up one centered,
  one-cell-wide spout that almost never closes. That spout carries few
  particles (8.5–10.8 g/s, against 20–32 g/s with the 2 × 2 inlet) and
  throws them 2.5–3× higher than the experiment's fountain. The annulus
  stays packed, which is why the wall layer is right.
- **The next lever** is how much gas crosses from the spout into the
  annulus. On this grid the drag closure at annulus packing governs that.
  A drag-closure screen was submitted on 2026-11-01
  ([file 09](09-findings-and-open-items.md)).

## Coating-stage and charge cases (c04–c09)

C26 has no data for these, so they are compared only with c02 and with
the same cases on the 2 × 2 inlet.

![coating-stage maps](../results/fig_c26_fig4_compose_diameter_fullface_bore_mom_nointerp_v232_pivlab.png)

![charge maps](../results/fig_c26_fig4_compose_charge_fullface_bore_mom_nointerp_v232_pivlab.png)

*PIVlab maps, t = 2.4–3.0 s, same color scale as above. Top:
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
[c26_sweep_summary_fullface_nointerp_v232.csv](../results/c26_sweep_summary_fullface_nointerp_v232.csv).*

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

## Data files

Case definitions are in [params/sweep_C26.yaml](../params/sweep_C26.yaml).

- [c26_pivlab_rerun_fullface_nointerp_v232.csv](../results/c26_pivlab_rerun_fullface_nointerp_v232.csv):
  c01–c03 map scores and wall velocities.
- [c26_sweep_summary_fullface_nointerp_v232.csv](../results/c26_sweep_summary_fullface_nointerp_v232.csv):
  sweep quantities for c01–c09, with figures for the
  [flow](../results/fig_c26_sweep_flow_fullface_nointerp_v232.png),
  [diameter](../results/fig_c26_sweep_diameter_fullface_nointerp_v232.png) and
  [charge](../results/fig_c26_sweep_charge_fullface_nointerp_v232.png) blocks.
- [fig_c26_pivlab_maps_fullface_nointerp_v232.png](../results/fig_c26_pivlab_maps_fullface_nointerp_v232.png):
  the c01–c03 maps above, plus a PIVlab variant with its CLAHE contrast
  step off (a diagnostic, [file 07](07-measurement-and-scoring.md)).
- Input decks: `cases/C26_3d/c0X_fullface_bore_mom_nointerp_v232/`. Raw
  MFiX output is not in git; it is on `s3://triso-cfd-mfix/results/<run>/`.

## The 2 × 2 control inlet

The other primary sweep (`_apex_v232`, 2026-10-27) differs in three ways:

- a mass-only source over the 2 × 2 cells of a closed cone tip;
- a 24-cell grid with the axis on a cell face;
- LINEAR_HAT interpolation.

None of its c01–c03 maps match C26 either. Its c02 scores 0.778 / 0.825 /
0.0048 (CC and MSE miss their floors), and its wall layer runs 3.0–6.3×
C26's model. Results:

- [c26_pivlab_rerun.csv](../results/c26_pivlab_rerun.csv) and
  [maps](../results/fig_c26_pivlab_maps.png);
- [c26_sweep_summary_apex_v232.csv](../results/c26_sweep_summary_apex_v232.csv).

Why both inlets exist is in
[file 04](04-gas-inlet-layouts.md).

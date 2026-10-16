# Day-10 R5 — spouted-bed startup: published precedents and decision

Written 2026-10-16. Companion to `notes/day10_r5_orifice_clog.md` (cause) and
`notes/day10_r5_ic_review.md` (independent IC review).

## Question

R1–R5b seed the bed as four quadrant strips around a 6 × 6 mm empty central
column (|x|, |z| < 3 mm) spanning the full bed height. Together with the
impulsive 30 m/s MI start, this is the "Fix A" pattern copied from NETL's
`pic/spouted_bed_3d` tutorial. Does published practice support it, or how do
others actually start a spouted-bed simulation?

## Precedents

| Source | Initial bed | Empty channel left in the bed? | Gas start | Startup discarded |
|---|---|---|---|---|
| Moliner et al. 2018, CEJ, doi 10.1016/j.cej.2018.11.164 (TFM vs CFD-DEM, spouted bed, glass beads) | TFM: solids "evenly patched all over the domain with ε_s = 0.65". DEM: particles "settled into the lower part of the reactor through an injection step in absence of air" | No | Impulsive, 26.68 m/s at t = 0 | First 1 s ("to achieve steady state"); statistics 1–2 s |
| Lethe CFD-DEM, arXiv 2211.01185, spouted-bed case (vs Yue et al. PIV) | Static packed bed, 200 mm, 175 800 particles | No. The "small channel at the bottom" is an inlet pipe added to the geometry (like our stub), not a gap in the bed | **Ramped** 0 → 20.8 m/s over 0.05 s, "to prevent an initial shock … due to the sudden introduction of a high velocity at t = 0" | Not stated in excerpt |
| Same paper, fluidized-bed case | DEM fill, then 1 s with no gas "to allow the particles to settle" | No | Gas started after settling | — |
| NETL MFiX-DEM spout-fluid study (WCPT8 2018; Link et al. 2005 benchmark) | "Particles automatically in bed generated based on initial specified height" | No | — | — |
| NETL tutorial `tutorials/pic/spouted_bed_3d` (MFiX 26.1.2 install) | Loose, ic_ep_s = 0.2 | Yes, but it is the inside of a **draft tube** in the STL (inner r ≈ 7.3 mm) | Impulsive | — |
| NETL tutorial `tutorials/dem/spouted_bed_pseudo_2d` ("Pseudo-2D DEM spouted bed with draft tube") | Loose, ic_ep_s = 0.25 | Yes, |z| < 7.5 mm, also the inside of a **draft tube** | Impulsive (bc_v_g = 25 m/s) | — |
| ORNL/TM-2006/520 (S1) | Not described | — | — | — |

## Findings

1. **No precedent for a free-standing empty channel.** Every study without a
   draft tube starts from a packed or settled bed with no channel. The only
   empty central columns we found are NETL tutorials where the column is the
   inside of a draft tube, held open by a wall. Our vessel has no draft tube,
   so the copied gap is an unsupported pre-opened spout. It also biases
   PLAN §1.1 Q1: through spouting hysteresis, a pre-opened channel can spout
   below the true U_ms. That would corrupt the Day-13 low-flow (U/U_ms ≈ 1.1)
   sweep cases.
2. **Settle, then start the gas.** The common sequence is to let particles
   settle into a packed bed without gas, then apply the flow. In MP-PIC,
   seeding directly at ε_s = 0.57 sits at the Harris–Crighton close-pack
   singularity (1 − des_epg_clip = 0.58; see `day10_r5_ic_review.md`). NETL's
   loose seeding (0.2) is the practical MP-PIC way to do the settle step.
3. **Ramp the jet.** There is a direct precedent at nearly the same d_p
   (500 µm) and jet speed (0 → 20.8 m/s over 50 ms) for avoiding the
   artificial shock of an impulsive start. The physical breakthrough ΔP peak
   is expected regardless. It is part of startup and is excluded from
   averaging.
4. **Discard a real startup period.** Moliner discards 1 s. Averaging windows
   must start after breakthrough has stabilized, judged from the data
   (ΔP, dt, inventory trends), not at a fixed t = 0.1 s.

Weak spots: no MP-PIC spouted-bed paper stating its startup procedure was
found, so the MP-PIC-specific choice (loose seed) rests on NETL tutorial
practice. Only the sources above were read in full text.

## Bonus: S1 Fig. 1 provenance resolved

The S1 Fig. 1 caption reads: "The solid blue line is from simulation and the
red line with points from experiment". The text adds that velocities were
"measured with a Prandtl static probe". Fig. 1 is therefore an **experimental**
validation target (with ORNL's MFiX result overlaid). `data/digitized/README.md`
audit item 1 can be closed, and the placeholder CSV should be replaced by a
real digitization from the report.

## Decision (2026-10-16)

Let R5b (loose seed, chimney retained, impulsive start) finish as an
isolated-variable data point. Next smoke **R5c**:

- no central chimney: the loose seed fills the whole cone/column cross-section
  (no orifice-footprint exclusion unless the smoke shows it is needed);
- ic_ep_s = 0.2 loose seed (NETL), inventory-conserving height;
- inlet held at a small v0 during a settling period (~0.1 s), then ramped
  linearly to 30 m/s over ~0.05 s via `BC_V_g_ramp_*` (keywords per
  `day10_r5_ic_review.md`; must be re-verified, including the inlet
  velocity-limit behaviour at small v0);
- tstop ≈ 0.3 s; evaluate with `check_smoke.py`, using a metric window that
  starts after breakthrough.

Sources:
- Moliner et al. 2018 accepted manuscript: https://aura.abdn.ac.uk/server/api/core/bitstreams/e4189098-ba9a-420d-9423-249226bae4c9/content
- arXiv 2211.01185: https://arxiv.org/pdf/2211.01185
- NETL WCPT8 presentation: https://mfix.netl.doe.gov/download/2018-aiche-presentations/WCPT8%20Spout-Fluidized-Beds-Presentation-Subhodeep.pdf
- ORNL/TM-2006/520: https://info.ornl.gov/sites/publications/files/Pub2422.pdf

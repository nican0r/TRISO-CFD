# Day 6-7 — Summary

## What was implemented

Built the project's first **parametric MFiX sweep**: a 2D flat-bottom bubbling
bed of Case D particles (500 µm, 6 000 kg/m³) in air at 20 °C, driven at
eight inlet superficial velocities spanning U/U_mf = {0.3, 0.5, 0.7, 0.9, 1.1,
1.3, 1.6, 2.0} on a single base mesh (38×150 cells). The sweep is the first
head-to-head between a correlation prediction (Wen-Yu U_mf for Case D in air
@ 20 °C = **0.4062 m/s**) and the TFM solver's own answer. Deliverable is the
fluidization curve `fig_d7_fluidization_curve.png` — simulated time-averaged
ΔP(U) overlaid with the Ergun packed-bed branch and the bed-weight-per-area
plateau, with U_mf,sim read from the crossover.

The 2× fine-mesh refinement check originally planned at U/U_mf = {0.7, 1.1,
1.6} was dropped on 2026-10-04: on the 76×300 grid the Schaeffer frictional
stress at the IC1/IC2 bed/freeboard step was severely under-resolved, driving
solids velocities to O(100 m/s) and the time step down to `dt_min = 1e-7`.
Mesh independence is no longer part of this study.

Key artifacts:

- `params/particles.yaml` — added the `fb_sweep` block (H_static_m, eps_s_bed,
  u_over_umf, t_startup_s, t_avg_s, dt_init_s, gas_case).
- `params/geometry.yaml` — added the `fb_sweep_2d` block (width 0.050 m,
  height 0.30 m, base-mesh counts).
- `cases/fb_sweep/template.mfx.j2` — Jinja2 template; every keyword resolves
  to the MFiX 26.1.2 reference entry audited in `notes/mfx_anatomy.md`.
- `src/make_case.py` — `render_case()` + `build_fb_sweep_params()`: reads the
  three YAML parameter files, calls `gasprops.gas_props` for ρ_g [kg/m³] and
  µ_g [Pa·s] of air at the gas case, substitutes into the template. Uses
  `StrictUndefined` so a missing variable fails loudly rather than writing a
  silently-broken `.mfx`.
- `src/run_fb_sweep.py` — enumerates the full plan (8 base runs) from the
  YAML, renders one `.mfx` + `manifest.json` per run under
  `cases/fb_sweep/u<ratio>_base/`, and emits `run_all.sh`. The runner can
  be configured for serial or N-way parallel execution via `--parallel N`
  (xargs -P scheduler); the production run is 6-way parallel on 10 physical
  cores.
- `src/post/fluidization_curve.py` — `time_averaged_dp(run_dir, t_avg)` plus
  `collect_sweep(...)` and `plot_fluidization_curve(...)`. Produces
  `results/fb_sweep_dp.csv` and `results/fig_d7_fluidization_curve.png`.
- `tests/test_make_case.py` — six tests pin the template rendering (bc_v_g,
  imax/jmax, Case D particle, air gas properties, bed IC bounds,
  StrictUndefined behaviour). All six pass.

## How it works

1. **Parameter wiring.** `build_fb_sweep_params(U_in, imax, jmax, tstop, run_name)`
   loads `params/particles.yaml`, `params/geometry.yaml`, `params/gases.yaml`
   and asks Cantera (`src.gasprops.gas_props`) for ρ_g [kg/m³], µ_g [Pa·s] at
   the fb_sweep gas case (`air_20C`: 293.15 K, 101325 Pa, N₂/O₂ = 0.79/0.21).
   Result: ρ_g = 1.1994 kg/m³, µ_g = 1.830×10⁻⁵ Pa·s. These, together with
   d_p = 500 µm and ρ_p = 6000 kg/m³ from `case_D`, are substituted into the
   template directly as `ro_g0`, `mu_g0`, `d_p0(1)`, `ro_s0(1)` (CLAUDE.md
   rules 2 and 3).

2. **Sweep sizing.** The ratios `u_over_umf = {0.3, …, 2.0}` are multiplied
   by U_mf,Wen-Yu = 0.4062 m/s (recomputed via
   `src.correlations.umf_wen_yu` so it never drifts from the Day-2 hand
   calc). This gives inlet velocities U ∈ [0.122, 0.812] m/s. Each run uses
   `tstop_s = t_startup_s + t_avg_s = 5.0 + 3.0 = 8.0 s` on the 38×150 base
   grid.

3. **Template.** Based directly on the Day-5 tutorial `.mfx`, so the keyword
   audit in `notes/mfx_anatomy.md` applies unchanged. Modifications are
   mechanical: switch the solid to Case D, switch the fluid to air at 20 °C,
   put the MI inlet across the whole bottom edge (flat-bottom, uniform
   inflow — "no spout yet"), place the packed-bed IC2 only over
   0 ≤ y ≤ H_static with ε_s = 0.55, and leave the freeboard IC1 above
   with ε_s = 0. All TFM closure keywords (`c_e`, `ep_star`, `phi`,
   `friction_model = SCHAEFFER`, etc.) match the tutorial.

4. **Choice of ε_s = 0.55 for the initial bed.** The first smoke run used
   ε_s = 0.60 (the Day-4 "loose random packing" convention) and the solver
   crashed at t ≈ 0 with velocities ~45 m/s at the bed/freeboard interface
   and the time step collapsing below `dt_min`. The physical cause is the
   `SCHAEFFER` frictional-stress term activating when ε_s > 1 − ep_star =
   0.58; at ε_s = 0.60 the frictional stress is divergent exactly across
   the IC1/IC2 step at y = H_static, which the solver cannot resolve.
   Dropping to ε_s = 0.55 puts the initial bed cleanly below the Schaeffer
   threshold while remaining a realistic Geldart-D loose packing. The
   bed-weight-per-area target scales with the chosen ε_s and becomes
   (ρ_p − ρ_g)·ε_s·g·H_static = (6000 − 1.2)·0.55·9.80665·0.10 ≈ **3234 Pa**.

5. **Post-processing.** For each run directory, `time_averaged_dp` reads
   `BACKGROUND.pvd`, pulls the P_g cell field from every frame, computes
   ΔP(t) = ⟨P_g⟩_{bottom row} − ⟨P_g⟩_{top row} [Pa] (same recipe as Day
   5), and averages over the trailing `t_avg_s = 3.0 s` window
   (CLAUDE.md rule 6). `collect_sweep` tabulates `(ratio, U, mesh, imax,
   jmax, dp_mean_Pa, n_frames, t_end, status)` to `results/fb_sweep_dp.csv`.
   `plot_fluidization_curve` overlays the simulated points on the Ergun
   dP/L × H_static curve (ε = 1 − 0.55 = 0.45) and the bed-weight plateau,
   and reads U_mf,sim by linear interpolation of the crossover between
   simulated points and the plateau line.

6. **Verification.**
    - *Template rendering:* six unit tests pass — see `tests/test_make_case.py`.
      bc_v_g(1), imax, jmax, d_p0(1), ro_s0(1) all substitute correctly;
      ρ_g from Cantera lands in [1.15, 1.25] kg/m³ and µ_g in
      [1.7, 1.9]×10⁻⁵ Pa·s; missing-key renders raise `UndefinedError`.
    - *Smoke run:* U/U_mf = 1.1 at tstop = 0.5 s on the base mesh completed
      cleanly (wall time 3 min 41 s). ΔP(t) was recovered via the Day-5
      post-processing pipeline; mean over the last 0.3 s = **2149 Pa**,
      against a plateau target of 3234 Pa — i.e. only ~0.5 s simulated time
      is nowhere near the statistically steady state the step calls for.
      This is a pipeline-correctness check, not a physics check.
    - *Acceptance criteria (plateau within 5 % of bed weight, U_mf,sim
      within ~25 % of Wen-Yu):* **pending** the full sweep. First launch
      was serial; measuring the actual wall time on this box gave ~4.5 min
      of wall time per simulated second (`u0p30_base`: 20 vtu frames =
      1.0 s simulated in ~4.5 min wall), i.e. ~36 min per 8 s base-mesh
      case → ~5 h total serial. The serial scheduler was killed and the
      sweep is run 6-way parallel via `xargs -P6` (nohup, master pid in
      `cases/fb_sweep/sweep_master.pid`). All 8 base cases finish well
      inside a workday. Once outputs land, running
      `python -m src.post.fluidization_curve` produces the figure and the
      two acceptance read-outs are computed automatically from the
      returned summary dict (`plateau_sim_pa`, `plateau_rel_err`,
      `umf_sim_m_s`).

    - *Fine-mesh refinement check dropped on 2026-10-04.* The planned
      76×300 runs at U/U_mf = {0.7, 1.1, 1.6} tripped a Schaeffer
      frictional-stress blow-up at the IC1/IC2 bed/freeboard step: solids
      velocities of O(100 m/s) appeared on the first frames and the time
      step collapsed toward `dt_min = 1e-7`, so the solver ran but did
      not advance in simulated time. The underlying instability is the
      same one that forced ε_s = 0.55 on the base mesh — on the fine
      mesh the step is more sharply resolved, so the margin is lost. No
      numerical-error bound on the base mesh is reported for this study.

## What it models (physically)

A **flat-bottom bubbling bed** with a single Case D phase and air inflow.
This is deliberately simpler than the conical spouted bed the project is
aimed at — the point is to isolate one well-studied bed behaviour (the
packed → fluidized transition) and demonstrate that the TFM solver, with
the Day-5 closure choices, reproduces it quantitatively for Case D
particles before any geometric complication is added.

The physics the sweep actually captures:

- **Packed-bed pressure drop (U < U_mf).** For U below minimum fluidization
  the bed is stationary and dP/L is set by the Ergun equation (viscous +
  inertial terms through the random-packed solids). Multiplied by
  H_static, this gives ΔP vs U rising linearly then quadratically — the
  overlay curve in the figure. Mismatch here flags the drag closure (we use
  Syamlal-O'Brien defaults with `drag_c1 = 0.8`, `drag_d1 = 2.65`) and the
  initial packing fraction.
- **Fluidized-bed plateau (U > U_mf).** Once the bed fluidizes, the gas
  must support the full submerged weight of the solids, so ΔP = (ρ_p − ρ_g)
  ·ε_s·g·H_static is a mass-conservation statement, independent of U.
  Agreement to 5 % is a direct check that the TFM is conserving solids
  mass and bed inventory as it bubbles. **For the coater, this is the
  calibration point for every later ΔP read-out**: the gas-blower sizing
  for Case V and Case D (PLAN §1.3) uses measured ΔP vs U as the signal
  that spouting is sustained, so a solver that can't match the simpler
  bubbling plateau to within 5 % cannot be trusted at the spout.
- **U_mf transition.** The crossover of the Ergun branch with the plateau
  defines U_mf. Matching Wen-Yu to ~25 % — the scatter between published
  correlations themselves — demonstrates that the TFM reaches fluidization
  at roughly the right gas velocity. For the coater this is the first-order
  check that the solver's drag model is in the right ballpark **for Case D
  particles specifically** (which Day 2 flagged as Geldart D, where
  spouting rather than smooth bubbling is the operating regime). A U_mf
  that is far off here would mean the design-case spouting velocity
  predictions from Day-8 onward cannot be trusted.
What the sweep **deliberately does not** model: bubble statistics (size,
rise velocity, frequency), freeboard entrainment, particle-size
distribution effects, 3D wall effects, and anything about the conical
spout. All of those are deferred to later days per PLAN §1.3. Here the
only read-out is a one-number-per-run ΔP̄, chosen because it is the
cleanest and most forgiving signature of correct TFM behaviour at this
stage.

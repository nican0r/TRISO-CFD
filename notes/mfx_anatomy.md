# MFiX .mfx anatomy — TFM 2D fluidized-bed tutorial

Line-by-line annotation of `cases/tut_tfm2d/example_2d_fluidized_bed.mfx` (MFiX
26.1.2). Every keyword below is a documented entry in the MFiX 26.1.2 keyword
reference (NETL "MFiX User Manual / Keyword Reference", release notes for
26.1.2); no undocumented keywords are present.

Grouping follows the sections in the file. Values are quoted verbatim; units
are those documented for the keyword (SI, since `units = 'SI'`).

## Run controls

| keyword | value | meaning |
|---|---|---|
| `description` | `'2D fluidized bed'` | Free-text label stored in the output. |
| `run_name` | `'example_2d_fluidized_bed'` | Base name for all output files. |
| `units` | `'SI'` | Selects SI unit system (default). |
| `run_type` | `'new'` | Fresh start; alternatives are `restart_1`, `restart_2`. |
| `tstop` | `5.0` | Simulated end time [s]. |
| `dt` | `0.001` | Initial time-step [s]. Adaptive from here. |
| `dt_min` | `1.0e-7` | Smallest allowed `dt` before giving up [s]. |
| `dt_max` | `0.01` | Largest allowed adaptive `dt` [s]. |
| `dt_fac` | `0.9` | Multiplier applied when solver reduces/raises `dt`. |
| `res_dt` | `0.01` | Period for writing the `.RES` restart file [s]. |
| `chk_batchq_end` | `.True.` | Honour batch-queue end signals and shut down cleanly. |
| `drag_c1` | `0.8` | First tunable in the Syamlal-O'Brien drag law. |
| `drag_d1` | `2.65` | Second tunable in the Syamlal-O'Brien drag law. |
| `energy_eq` | `.False.` | Energy equation off — isothermal. |
| `des_rigid_motion(1)` | `.False.` | DES rigid-body motion flag for solids phase 1 (unused here, TFM). |
| `momentum_x_eq(0)` | `.True.` | Solve gas-phase x-momentum. The `(0)` index selects the gas phase. |
| `momentum_y_eq(0)` | `.True.` | Solve gas-phase y-momentum. |
| `momentum_z_eq(0)` | `.True.` | Solve gas-phase z-momentum (benign in 2D with `no_k`). |
| `project_version` | `'6'` | GUI project format version (metadata only). |
| `species_eq(0)` | `.False.` | No species transport for gas. |
| `species_eq(1)` | `.False.` | No species transport for solids phase 1. |

## Physical parameters

| keyword | value | meaning |
|---|---|---|
| `gravity_x` | `0.0` | Gravity vector component [m/s²]. |
| `gravity_y` | `-9.80665` | Standard gravity downward [m/s²]. |
| `gravity_z` | `0.0` | Gravity vector component [m/s²]. |

## Cartesian cut-cell grid

| keyword | value | meaning |
|---|---|---|
| `cartesian_grid` | `.False.` | No cut-cell geometry; use the regular IJK grid. |
| `use_stl` | `.False.` | No STL geometry file imported. |

## Numeric controls

| keyword | value | meaning |
|---|---|---|
| `detect_stall` | `.True.` | Abort if `dt` keeps shrinking without progress. |
| `max_nit` | `50` | Max outer SIMPLE iterations per time-step. |
| `norm_g` | `0.0` | Fluid-phase pressure correction normalization (0 = use default/absolute). |

## Geometry

| keyword | value | meaning |
|---|---|---|
| `coordinates` | `'CARTESIAN'` | Cartesian IJK grid (not cylindrical). |
| `imax` | `80` | Cell count along x. |
| `jmax` | `180` | Cell count along y. |
| `x_max` | `0.1` | Domain extent in x [m]. |
| `x_min` | `0.0` | Domain origin in x [m]. |
| `y_max` | `0.3` | Domain extent in y [m]. |
| `y_min` | `0.0` | Domain origin in y [m]. |
| `z_max` | `1.0` | Thickness in z [m] (ignored under `no_k`). |
| `z_min` | `0` | z-origin [m] (ignored under `no_k`). |
| `no_k` | `.True.` | Collapse to 2D (xy plane); no k-direction solve. |

## Fluid phase

| keyword | value | meaning |
|---|---|---|
| `mu_g0` | `1.8e-5` | Constant gas viscosity [Pa·s] (`mu_g_model='CONSTANT'`). |
| `mu_g_model` | `'CONSTANT'` | Use `mu_g0`, not a temperature model. |
| `ro_g0` | `1.0` | Constant gas density [kg/m³] (incompressible). |

## Solids phase 1 — "glass beads"

| keyword | value | meaning |
|---|---|---|
| `mmax` | `1` | Number of solids phases. |
| `solids_model(1)` | `'TFM'` | Two-fluid Eulerian model for phase 1. |
| `d_p0(1)` | `2.0e-4` | Particle diameter [m] (200 µm). |
| `ro_s0(1)` | `2500.0` | Particle density [kg/m³] (glass). |
| `nmax_s(1)` | `0` | Species count in solids phase 1 (none). |
| `k_s0(1)` | `1.0` | Solids thermal conductivity [W/(m·K)] (unused — energy_eq off). |
| `ks_model(1)` | `'BAUER'` | Solids-conductivity model (unused here). |

## Initial condition 1 — "Background" (whole domain)

Keyword family `ic_*` with `(1)` selecting IC region 1.

| keyword | value | meaning |
|---|---|---|
| `ic_x_e(1)` | `0.1` | East extent of IC region [m]. |
| `ic_x_w(1)` | `0.0` | West extent [m]. |
| `ic_y_s(1)` | `0.0` | South extent [m]. |
| `ic_y_n(1)` | `0.3` | North extent [m]. |
| `ic_z_b(1)` | `0.0` | Bottom extent [m]. |
| `ic_z_t(1)` | `1.0` | Top extent [m]. |
| `ic_ep_g(1)` | `1.0` | Initial gas volume fraction (empty). |
| `ic_t_g(1)` | `293.15` | Initial gas temperature [K]. |
| `ic_u_g(1)` | `0.0` | Initial gas x-velocity [m/s]. |
| `ic_v_g(1)` | `0.0` | Initial gas y-velocity [m/s]. |
| `ic_w_g(1)` | `0.0` | Initial gas z-velocity [m/s]. |
| `ic_p_g(1)` | `101325.0` | Initial gas pressure [Pa]. |
| `ic_p_star(1)` | `0.0` | Initial solids-phase pressure [Pa]. |
| `ic_ep_s(1,1)` | `0.0` | Initial solids volume fraction (phase 1) in region 1. |
| `ic_t_s(1,1)` | `293.15` | Initial solids temperature [K]. |
| `ic_u_s(1,1)` | `0.0` | Initial solids x-velocity [m/s]. |
| `ic_v_s(1,1)` | `0.0` | Initial solids y-velocity [m/s]. |
| `ic_w_s(1,1)` | `0.0` | Initial solids z-velocity [m/s]. |

## Initial condition 2 — "bed" (packed region)

Same keyword family, index `(2)`.

| keyword | value | meaning |
|---|---|---|
| `ic_x_e(2)`, `ic_x_w(2)` | `0.1`, `0.0` | Region x extents [m]. |
| `ic_y_s(2)`, `ic_y_n(2)` | `0.0`, `0.3` | Region y extents [m]. |
| `ic_ep_g(2)` | `0.6` | Gas fraction = 0.6 (bed voidage). |
| `ic_t_g(2)`, `ic_t_s(2,1)` | `293.15` | Temperatures [K]. |
| `ic_u_g(2)…ic_w_g(2)` | `0.0` | Zero initial gas velocity components. |
| `ic_p_star(2)` | `0.0` | Initial solids pressure [Pa]. |
| `ic_ep_s(2,1)` | `0.4` | Solids fraction = 0.4 (so ε_g + ε_s = 1). |
| `ic_u_s(2,1)…ic_w_s(2,1)` | `0.0` | Zero initial solids velocity components. |

> Note: the GUI assigned the "bed" region to the whole domain (`y_s=0, y_n=0.3`); so the mapping `ic_regions = [[[1], ["Background"]], [[2], ["bed"]]]` leaves the bed covering the full column at `ε_s=0.4`. This is a quirk of this tutorial file; it does not represent a packed bed sitting at the bottom.

## Boundary condition 1 — "inlet" (bottom face, Mass Inflow)

| keyword | value | meaning |
|---|---|---|
| `bc_type(1)` | `'MI'` | Mass-inflow boundary. |
| `bc_x_e(1)`, `bc_x_w(1)` | `0.1`, `0.0` | x extents of inlet [m]. |
| `bc_y_s(1)`, `bc_y_n(1)` | `0.0`, `0.0` | y extent — zero-thickness strip at y=0. |
| `bc_ep_g(1)` | `1.0` | Pure gas at the inlet. |
| `bc_p_g(1)` | `101325` | Gas pressure at inlet [Pa]. |
| `bc_t_g(1)` | `293.15` | Gas temperature at inlet [K]. |
| `bc_u_g(1)` | `0.0` | Gas x-velocity at inlet [m/s]. |
| `bc_v_g(1)` | `0.25` | Gas y-velocity at inlet [m/s] — the superficial inflow. |
| `bc_w_g(1)` | `0.0` | Gas z-velocity at inlet [m/s]. |
| `bc_ep_s(1,1)` | `0.0` | No solids injected. |
| `bc_t_s(1,1)` | `293.15` | Solids temperature placeholder [K]. |
| `bc_u_s(1,1)…bc_w_s(1,1)` | `0.0` | Solids velocities at inlet [m/s]. |

## Boundary condition 2 — "outlet" (top face, Pressure Outflow)

| keyword | value | meaning |
|---|---|---|
| `bc_type(2)` | `'PO'` | Pressure outflow boundary. |
| `bc_x_e(2)`, `bc_x_w(2)` | `0.1`, `0.0` | x extents [m]. |
| `bc_y_s(2)`, `bc_y_n(2)` | `0.3`, `0.3` | y extent — zero-thickness strip at y=y_max. |
| `bc_p_g(2)` | `101325` | Prescribed outlet pressure [Pa]. |

## Walls

| keyword | value | meaning |
|---|---|---|
| `bc_default_walls` | `.False.` | Do not auto-populate the remaining faces with free-slip walls — only the explicit BCs above apply. |

## VTK output

| keyword | value | meaning |
|---|---|---|
| `write_vtk_files` | `.True.` | Enable VTK output. |
| `time_dependent_filename` | `.True.` | Append time index to each `.vtu` name. |
| `vtk_filebase(1)` | `'Background'` | Prefix for VTK region 1 output. |
| `vtk_x_e(1)…vtk_z_t(1)` *(implicit via `vtk_x_e(1)=0.1`, `vtk_x_w(1)=0.0`, `vtk_y_s(1)=0.0`, `vtk_y_n(1)=0.3`)* | | Spatial extent of VTK region 1 [m]. |
| `vtk_data(1)` | `'C'` | Cell-centered output (vs `'P'` for particle data). |
| `vtk_dt(1)` | `0.1` | VTK output interval [s]. |
| `vtk_nxs(1)`, `vtk_nys(1)`, `vtk_nzs(1)` | `0` | Number of slice planes in each direction (0 = full volume). |
| `vtk_ep_g(1)` | `.True.` | Write `EP_G` to VTK. |
| `vtk_p_g(1)` | `.True.` | Write `P_G` to VTK. |
| `vtk_vel_g(1)` | `.True.` | Write the gas velocity vector to VTK. |

## SPx output

`spx_dt(1..9) = 0.1` sets the write interval [s] for each SPx binary (SP1–SP9),
which hold the legacy single-phase/multi-phase variable groups (pressure,
density, velocities, volume fraction, temperature, species, scalars,
granular-energy, reaction rates — one group per SPx index per the keyword
reference).

## Residuals

`resid_string(1..6) = 'P0','U0','V0','P1','U1','V1'` requests residual
monitoring for gas pressure (P0), gas u- and v-velocity (U0, V0) and the
corresponding solids-phase-1 fields (P1, U1, V1). Each slot is a documented
residual tag.

## Two-fluid model (TFM) closure

| keyword | value | meaning |
|---|---|---|
| `c_e` | `0.95` | Particle–particle restitution coefficient. |
| `c_f` | `0.1` | Friction coefficient (TFM frictional closure). |
| `ep_star` | `0.42` | Packing limit for gas volume fraction (minimum ε_g). |
| `friction_model` | `'SCHAEFFER'` | Frictional-stress model. |
| `kt_type` | `'ALGEBRAIC'` | Granular-temperature closure (algebraic form). |
| `phi` | `30.0` | Internal angle of friction [deg]. |
| `phi_w` | `11.3` | Wall angle of friction [deg]. |
| `rdf_type` | `'CARNAHAN_STARLING'` | Radial-distribution function used by the kinetic theory. |

## Discrete element model (DEM) parameters

These are only consumed when a phase has `solids_model='DEM'`; here they are
retained defaults from the GUI (phase 1 is TFM, so they are inert).

| keyword | value | meaning |
|---|---|---|
| `des_en_input(1)` | `0.9` | Particle–particle restitution for phase 1. |
| `des_en_wall_input(1)` | `0.9` | Particle–wall restitution for phase 1. |
| `kn` | `1000` | Particle–particle normal spring stiffness [N/m]. |
| `kn_w` | `1000` | Particle–wall normal spring stiffness [N/m]. |
| `mew` | `0.1` | Particle–particle friction coefficient. |
| `mew_w` | `0.1` | Particle–wall friction coefficient. |

## GUI metadata (not solver keywords)

Every line beginning with `#!MFIX-GUI` is metadata the GUI writes next to the
solver keywords. The MFiX solver skips these lines; they let the GUI
round-trip region definitions, parameter aliases, and job-queue options. Not
part of the keyword reference, and not required for the solver to run.

## Audit

Every solver keyword above resolves to one entry in the MFiX 26.1.2 keyword
reference. No undocumented keywords are present in the file.

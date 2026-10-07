# 1. Project definition

## 1.1 Goal

Build, sanity-check and document a hydrodynamic CFD model of a laboratory-scale conical spouted-bed coater of the type used for TRISO fuel particle CVD. Use it to answer:

- At what gas flow does the bed spout stably, and how does that change from room temperature to coating temperature?
- What is the solids circulation pattern (spout → fountain → annulus) and roughly how long is one circulation cycle?
- How much time does a particle spend in the regions where deposition mainly happens (spout and fountain), and how uniform is that exposure across particles? This is a proxy for coating-thickness uniformity.


## 1.2 Scope

In scope: cold-flow hydrodynamics of the ORNL/UTK coater geometry (air at ambient, surrogate ZrO₂ particles), hot-gas-property hydrodynamics (argon/hydrogen at coating temperature on the same geometry), a passive tracer for gas residence time, and a small CFD-DEM run for particle-level statistics.

Out of scope for these two weeks: reactive chemistry (MTS → SiC, hydrocarbon → PyC), full heat transfer with wall heating, and turbulence-model studies. These are listed as next steps in the final write-up.

Why the Missouri S&T glass-bead case was dropped. The Day 8–9 attempts used the 0.076 m / 2.18 mm-bead column but TFM crashed at t = 0.072 s. The particle (2.18 mm) is too large relative to a mesh resolving the 9.5 mm orifice at ≥ 4 cells (cells ~ 1.1 d_p; orifice ~ 4.4 d_p). The continuum assumption that underpins TFM requires cells ≫ d_p (standard rule 10 d_p). The ORNL/UTK coater uses 500 µm ZrO₂ in a 50 mm / 4 mm column; the same mesh rule gives cells ≥ 5 mm comfortably above d_p, which brings TFM back into its domain of validity.

## 1.3 Two cases

Validation case (Case V): the ORNL/UTK cold mockup spouted bed used in the DOE AGR / NGNP TRISO coater program (S1, S2). A 0.05 m-ID straight column with a 60° conical base and a constant-diameter 0.004 m inlet orifice, 54.5 g of 500 µm spherical ZrO₂ (ρ_p = 6050 kg/m³), run with humidified compressed air modelled as dry air at ambient (humidity effect on gas properties neglected and noted). Operating point: 30 m/s inlet jet (≈ 22.6 L/min, ≈ 1.9 × U_ms). Validation targets are the centerline axial gas velocity (S1 Fig. 1) and the inlet-pressure pulsation frequency peak (S2 reports 17–25 Hz over U/U_ms = 1.3–2.3; S1 Fig. 2 reports ~10–11 Hz for the ORNL 2-inch coater near 600 K). Particle properties, inventory, and operating points come from S1 and S2, not guessed.

Design case (Case D): the same geometry run with (i) argon at ambient, (ii) argon at coating temperature (~1300–1600 °C), (iii) Ar/H₂ 50/50 mol at coating temperature, with the same 500 µm ZrO₂ surrogate; plus a UO₂ real-kernel sensitivity at 10 800 kg/m³. Also parametric: flow rate across U/U_ms ≈ 1.1–1.9, drag model, and particle-particle restitution coefficient (c_e is assumed in Case V, so Day 13 sweeps it).

## 1.4 Final deliverables (end of Day 14)

- A git repo with reproducible .mfx cases, run scripts, and post-processing scripts
- A Jupyter notebook analytics.ipynb with all hand calculations (terminal velocity, U_mf, U_ms, property tables)
- A validation figure set comparing simulation to S1/S2 published data
- A parametric study (flow rate, gas temperature, particle density, drag model, restitution coefficient)
- REPORT.md: 4–6 pages with figures, methods, validation, limitations and next steps, written to be a portfolio piece


# 2. Setup

## 2.1 Software

- MFiX, current release (26.1.2 installed). Install with the official conda-based instructions for your OS.
- ParaView for 3D visualization.
- Python env: numpy scipy pandas matplotlib jupyter pyvista meshio cantera. Cantera provides temperature-dependent gas viscosity and density; the gri30.yaml mechanism includes Ar, H₂ and N₂ with transport data.

## 2.2 Repo layout
coater-cfd/
├── REPORT.md               # final write-up (Day 14)
├── analytics.ipynb         # hand calcs & correlations
├── params/
│   ├── particles.yaml      # dp, rho_p, restitution, etc. per case
│   ├── gases.yaml          # T, P, composition → rho_g, mu_g (generated)
│   └── geometry.yaml       # D_c, D_i, cone angle, heights
├── src/
│   ├── correlations.py     # Ergun, Wen-Yu, Mathur-Gishler, drag, etc.
│   ├── gasprops.py         # Cantera wrapper
│   ├── make_case.py        # renders .mfx from templates + params
│   ├── post/               # readers & metrics for MFiX outputs
│   └── plotting.py
├── cases/
│   ├── tut_tfm2d/          # Day 5
│   ├── fb_sweep/           # Days 6–7
│   ├── V_*/                # validation runs
│   └── D_*/                # design runs
├── results/                # CSV + figures (small files only)
├── notes/                  # your daily notes & self-checks
└── tests/                  # unit tests for correlations


# 3. Resources

- S1. Pannala, S., Daw, C. S., Boyalakuntla, D., & Finney, C. E. A. (2006). "Process Modeling Phase I Summary Report for the Advanced Gas Reactor Fuel Development and Qualification Program." ORNL/TM-2006/520. https://info.ornl.gov/sites/publications/files/Pub2422.pdf. Model closures (Table 1), grid rules, validation data (Fig. 1 centerline velocity, Fig. 2 pulsation frequency).

- S2. Zhou, J., Bruns, D. D., Finney, C. E. A., Daw, C. S., Pannala, S., & McCollum, D. L. (2005). "Hydrodynamic Correlations with Experimental Results from Cold Mockup Spouted Beds for Nuclear Fuel Particle Coating." AIChE Annual Meeting, Cincinnati, paper 284f. Full text: https://skoge.folk.ntnu.no/prost/proceedings/aiche-2005/non-topical/Non%20topical/papers/284f.pdf. Source for the UTK cold mockup geometry, particles, inventory, pressure data (17–25 Hz spectral peaks sampled at 1000 Hz for 1 min).

- S3. Pannala, S., et al. (2007). "Simulating the Dynamics of Spouted-Bed Nuclear Fuel Coaters." Chemical Vapor Deposition 13(9):481–490. DOI 10.1002/cvde.200606562. Paywalled; use only if accessible.

- S4. Boyalakuntla, D., Pannala, S., Daw, C. S., Finney, C. E. A., et al. (2005). "Simulating the Hydrodynamics of Spouted Beds Using a Continuum Formulation." AIChE Annual Meeting, paper 209d. Extended abstract is free; the full paper is in the paid proceedings. Use only if accessible.

- Supporting: Wendt, D. S., Bewley, R. L., & Windes, W. E. (2007). "A Spouted Bed Reactor Monitoring System for Particulate Nuclear Fuel." INL/CON-07-12569. https://inldigitallibrary.inl.gov/sites/sti/sti/3874553.pdf. Independently confirms a 4 mm constant-diameter inlet on an INL coater (used to cross-check S2's 4 mm number against S2's abstract-level "0.04 cm" typo).

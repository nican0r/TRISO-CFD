# 1. Project definition

## 1.1 Goal

Build, sanity-check and document a hydrodynamic CFD model of a laboratory-scale conical spouted-bed coater of the type used for TRISO fuel particle CVD. Use it to answer:

- At what gas flow does the bed spout stably, and how does that change from room temperature to coating temperature?
- What is the solids circulation pattern (spout → fountain → annulus) and roughly how long is one circulation cycle?
- How much time does a particle spend in the regions where deposition mainly happens (spout and fountain), and how uniform is that exposure across particles? This is a proxy for coating-thickness uniformity.


## 1.2 Scope

In scope: cold-flow hydrodynamics (air/argon, surrogate particles), hot-gas-property hydrodynamics (argon/hydrogen properties at coating temperature), a passive tracer for gas residence time, and a small CFD-DEM run for particle-level statistics.

Out of scope for these two weeks: reactive chemistry (MTS → SiC, hydrocarbon → PyC), full heat transfer with wall heating, and turbulence-model studies. These are listed as next steps in the final write-up.

## 1.3 Two cases

Validation case (Case V): the Missouri S&T cold-flow spouted beds built for TRISO coater scale-up research (DOE-NERI project). The 0.076 m (3-inch) bed has a 1.14 m tall cylindrical section, a 60° conical base, and a 9.5 mm inlet orifice, run with dry compressed air. Particle properties, static bed heights and superficial gas velocities must be taken from the papers (see Resources R13–R15), not guessed.

Design case (Case D): the same geometry run with a TRISO-like surrogate: ZrO₂/YSZ spheres of ~500 µm diameter and ~6000 kg/m³ density (typical surrogate kernels). Run it first in argon at room temperature and then with gas properties at coating temperature (~1300–1600 °C). Real UO₂ kernels are ~10,800 kg/m³ and can be a sensitivity case.

## 1.4 Final deliverables (end of Day 14)

- A git repo with reproducible .mfx cases, run scripts, and post-processing scripts
- A Jupyter notebook analytics.ipynb with all hand calculations (terminal velocity, U_mf, U_ms, property tables)
- A validation figure set comparing simulation to correlations and published data
- A parametric study (flow rate, gas temperature, particle density, drag model)
- REPORT.md: 4–6 pages with figures, methods, validation, limitations and next steps, written to be a portfolio piece


# 2. Setup

## 2.1 Software

- MFiX, current release (25.x as of this writing). Install with the official conda-based instructions for your OS (R6). Windows works, but Linux or WSL2 is smoother for batch runs.
- ParaView for 3D visualization (R8).
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

- R13. "An advanced evaluation of spouted beds scale-up for coating TRISO nuclear fuel particles using Radioactive Particle Tracking (RPT)" (Missouri S&T, 2016). OSTI record with accepted manuscript: https://www.osti.gov/pages/biblio/1533774. It includes solids velocity and turbulence profiles intended as CFD benchmark data.
- R14. "An advanced evaluation of the mechanistic scale-up methodology of gas–solid spouted beds using radioactive particle tracking" (Missouri S&T, 2017). OSTI accepted manuscript: https://www.osti.gov/pages/biblio/1538701. This is the source of the 0.076 m / 0.152 m bed geometry.
- R15. Final management report, DOE-NERI DE-FC07-07ID14822, "Advancing the fundamental understanding and scale-up of TRISO fuel coaters via advanced measurement and computational techniques": https://www.osti.gov/servlets/purl/1054926
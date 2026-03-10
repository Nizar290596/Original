# mmcPbeFoam

`mmcPbeFoam` is an LES/FDF solver that extends the MMC framework with **Population Balance Equations (PBE)** for aerosol and nanoparticle dynamics. It enables simultaneous simulation of turbulent reacting flows and the evolution of a dispersed particulate size distribution (e.g. soot, nanoparticles, droplet populations).

## Algorithm Overview

Each time step performs the following steps:

1. **Density update** — solve continuity equation (`rhoEqn`)
2. **Momentum** — solve `UEqn`
3. **Diffusivity** — compute molecular and turbulent diffusivities
4. **Scalar transport** — solve mixture fraction (`XiEqn`) and enthalpy/species (`hYEqvE_Eqn`)
5. **Pressure correction** — PIMPLE outer correctors
6. **Post-PIMPLE turbulence fields** — compute normalised fields: kinematic viscosity (`nu`), turbulent viscosity (`nut`), dissipation rate (`epsilon`), turbulent kinetic energy (`k`) — all divided by density for PBE coupling
7. **Particle evolution** — advance aerosol MMC stochastic particle cloud (`basicAerosolReactingPopeCloud`)

## Key Physics / Models

| Component | Model |
|-----------|-------|
| Turbulence | Compressible LES |
| Reactive scalars | Aerosol stochastic particle cloud (`basicAerosolReactingPopeCloud`) |
| Particle dynamics | Population balance for size distribution evolution |
| Turbulence fields for PBE | `nu`, `nut`, `epsilon`, `k` (normalised by `rho`) |

## Additional Turbulence Fields

Unlike other MMC solvers, `mmcPbeFoam` writes normalised turbulence quantities after each PIMPLE loop for use by the PBE sub-model on the particles:

```
nu     = turbulence->mu() / rho
nut    = turbulence->mut() / rho
epsilon = turbulence->epsilon() / rho
k      = turbulence->k() / rho
```

## Usage

Use this solver when:

- The flow contains aerosol particles or nanoparticles whose size distribution evolves (nucleation, growth, coagulation)
- Soot formation and growth is of interest
- Nanoparticle synthesis in flame reactors is being studied

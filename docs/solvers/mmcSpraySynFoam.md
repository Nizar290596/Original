# mmcSpraySynFoam

`mmcSpraySynFoam` is an LES/FDF solver for **turbulent reacting dilute spray** combustion with comprehensive spray diagnostics. It extends `mmcDropletSprayFoam` with additional post-processing capabilities, including droplet size statistics and spray flux computations, and exposes detailed per-component timing information for performance analysis.

## Algorithm Overview

Each time step performs the following steps:

1. **Density update** — solve continuity equation (`rhoEqn`)
2. **Momentum** — solve `UEqn`
3. **Diffusivity** — compute molecular and turbulent diffusivities and effective viscosity
4. **Scalar transport** — solve mixture fraction (`XiEqn`) and enthalpy/species (`hYEqvE_Eqn`)
5. **Pressure correction** — PIMPLE outer correctors
6. **Turbulence fields** — compute `nu`, `nut`, `epsilon`, `k` (normalised by `rho`)
7. **Gas-phase particle evolution** — advance MMC aerosol stochastic particles
8. **Spray particle evolution** — advance liquid fuel spray cloud
9. **Spray diagnostics**:
   - Compute spray flux (`calcFlux`)
   - Compute mean droplet diameters: arithmetic mean D10 and Sauter mean D32 (`calcMeanDiameters`)

## Key Physics / Models

| Component | Model |
|-----------|-------|
| Turbulence | Compressible LES |
| Gas-phase scalars | Aerosol MMC stochastic particle cloud |
| Liquid phase | Spray thermo cloud with interpolation lookup tables |
| Spray diagnostics | Flux, D10, D32 droplet diameter statistics |

## Spray Diagnostics

Two additional post-processing fields are computed each time step:

- **Spray flux** — mass flux of liquid fuel across specified surfaces
- **D10** — arithmetic mean droplet diameter
- **D32** — Sauter mean diameter (surface-area-weighted mean; relevant for evaporation rate)

## Performance Profiling

This solver tracks and reports cumulative wall-clock time for each major component, making it useful for identifying performance bottlenecks in spray simulations.

## Usage

Use this solver when:

- Spray combustion simulations require detailed droplet size statistics
- Sauter mean diameter and spray flux need to be tracked for model validation
- Performance profiling of the spray–turbulence interaction is needed
- Synthetic fuel spray combustion is studied

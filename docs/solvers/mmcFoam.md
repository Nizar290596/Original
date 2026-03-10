# mmcFoam

`mmcFoam` is the base LES/FDF (Large Eddy Simulation / Filtered Density Function) solver for turbulent non-premixed reacting flows using the Multiple Mapping Conditioning (MMC) method. It couples an Eulerian compressible LES field with a Lagrangian stochastic particle cloud that carries the reactive scalars.

## Algorithm Overview

Each time step performs the following steps:

1. **Density update** — solve continuity equation (`rhoEqn`)
2. **Momentum** — solve `UEqn` for velocity
3. **Diffusivity** — compute molecular (`D = nu/Sc`) and turbulent (`Dt = nut/Sct`) diffusivities; effective diffusivity `DEff = D + Dt`
4. **Mixture fraction** — solve scalar transport equation (`XiEqn`)
5. **Enthalpy / species** — solve combined `hYEqvE_Eqn`
6. **Pressure correction** — PIMPLE outer correctors (`pEqn`)
7. **Turbulence update** — update LES turbulence model
8. **Particle evolution** — stochastic particle cloud is advanced (mixing, transport, reaction)

## Key Physics / Models

| Component | Model |
|-----------|-------|
| Turbulence | Compressible LES (`turbulentFluidThermoModel`) |
| Thermodynamics | `psiReactionThermo` |
| Reactive scalars | Stochastic Pope particle cloud (`basicReactingPopeCloud`) |
| Pressure–velocity | PIMPLE |
| Scalar transport | Multivariate differencing schemes |

## Combustion Mode

Non-premixed (diffusion flame). The mixture fraction `Z` is the primary conserved scalar; species and enthalpy are carried by the stochastic particles.

## Usage

This solver is the starting point for most non-premixed MMC simulations. Use it when:

- You have a non-premixed or partially-premixed flame
- Chemistry is handled by the particle cloud (no pre-tabulated tables needed)
- LES turbulence closure is desired

See `tutorials/oneStepFlameSheet` for a minimal working case.

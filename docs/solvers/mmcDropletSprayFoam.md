# mmcDropletSprayFoam

`mmcDropletSprayFoam` is an LES/FDF solver for **spray combustion** that couples the gas-phase MMC stochastic particle framework with a Lagrangian liquid-phase droplet spray cloud. It models dilute spray conditions where discrete fuel droplets evaporate and react in a turbulent gas flow.

## Algorithm Overview

Each time step performs the following steps:

1. **Density update** — solve continuity equation (`rhoEqn`)
2. **Momentum** — solve `UEqn`
3. **Diffusivity** — compute molecular and turbulent diffusivities
4. **Scalar transport** — solve mixture fraction (`XiEqn`) and enthalpy/species (`hYEqvE_Eqn`)
5. **Pressure correction** — PIMPLE outer correctors
6. **Turbulence update** — update LES model
7. **Gas-phase particle evolution** — advance MMC stochastic particles (`basicReactingPopeCloud`)
8. **Droplet evolution** — advance liquid fuel droplet cloud (`basicDropletSprayThermoCloud`), including evaporation, drag, heat transfer

## Key Physics / Models

| Component | Model |
|-----------|-------|
| Turbulence | Compressible LES |
| Gas-phase scalars | MMC stochastic particle cloud (`basicReactingPopeCloud`) |
| Liquid phase | Droplet spray cloud (`basicDropletSprayThermoCloud`) |
| Droplet physics | Evaporation, drag, heat/mass transfer |
| Pressure–velocity | PIMPLE |

## Combustion Mode

Non-premixed spray. Fuel vapour from droplet evaporation mixes with oxidiser in the gas phase and reacts via the MMC particle chemistry model.

## Two-Phase Particle Tracking

This solver maintains two separate Lagrangian clouds simultaneously:
- **Gas-phase MMC particles** — stochastic particles carrying reactive scalar PDFs
- **Fuel droplets** — physical droplets undergoing evaporation and Lagrangian tracking

Mass and energy exchange from droplet evaporation are coupled to the Eulerian gas-phase equations.

## Usage

Use this solver when:

- Fuel is injected as a liquid spray (e.g. gas turbines, diesel engines)
- Dilute spray regime applies (droplet volume fraction is small)
- LES turbulence resolution is desired

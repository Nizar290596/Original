# mmcDropletSprayDNSFoam

`mmcDropletSprayDNSFoam` is a **Direct Numerical Simulation (DNS)** solver for two-phase **spray combustion** that combines fully resolved DNS of the gas phase with Lagrangian tracking of liquid fuel droplets. Like `mmcDNSFoam`, it simultaneously computes the DNS solution and projects filtered fields onto a coarser mesh for MMC particle validation, but extends this capability to spray systems.

## Algorithm Overview

Each time step performs the following steps:

1. **Time synchronisation** — synchronise DNS and filter mesh clocks
2. **Density update** — solve continuity equation (`rhoEqn`)
3. **Diffusivity** — compute molecular diffusivity (`D = nu/Sc`)
4. **Momentum** — solve `UEqn` on DNS mesh
5. **Mixture fraction** — solve `XiEqn`
6. **Species** — solve full species transport equations (`YEqn`)
7. **Energy** — solve full energy equation (`EEqn`)
8. **Pressure correction** — `pEqn`
9. **Turbulence update**
10. **Filter mesh update** — project DNS fields onto coarser filter mesh
11. **Gas-phase particle evolution** — advance MMC stochastic particles (`basicReactingPopeCloud`)
12. **Droplet evolution** — advance liquid droplet cloud (`basicDropletSprayDNSThermoCloud`), including DNS-resolved evaporation and heat transfer
13. **Output** — write DNS and filtered fields; report execution times

## Key Physics / Models

| Component | Model |
|-----------|-------|
| Turbulence | None — all scales fully resolved (DNS) |
| Thermodynamics | `psiReactionThermo` |
| Chemistry | Full mechanism via `CombustionModel` |
| Gas-phase scalars | Full species + energy equations |
| Gas-phase particles | MMC stochastic particle cloud (`basicReactingPopeCloud`) |
| Liquid phase | DNS droplet cloud (`basicDropletSprayDNSThermoCloud`) |
| Filtering | DNS → filter mesh projection |

## Two-Phase DNS

This is the highest-fidelity solver in the suite, combining:

- **Full DNS** of the gas-phase turbulent reacting flow (no turbulence modelling)
- **DNS-resolved droplet dynamics** — evaporation, heat transfer, drag at fully resolved scales
- **MMC particle cloud** as a passive LES-scale surrogate for model validation

## Comparison with mmcDropletSprayFoam

| Feature | mmcDropletSprayFoam (LES) | mmcDropletSprayDNSFoam (DNS) |
|---------|---------------------------|------------------------------|
| Turbulence | LES subgrid model | Fully resolved |
| Droplet cloud | `basicDropletSprayThermoCloud` | `basicDropletSprayDNSThermoCloud` |
| Filtering | No | Yes (DNS → LES mesh) |
| Cost | Moderate | Very high |
| Purpose | Production simulations | Validation/reference |

## Usage

Use this solver when:

- A ground-truth DNS reference for spray combustion is required
- Validating LES spray-combustion MMC models against DNS-filtered data
- Studying the fundamental interaction between turbulence and droplet evaporation at resolved scales

# mmcFullChemPremixedFoam

`mmcFullChemPremixedFoam` is an LES/FDF solver for turbulent **premixed** reacting flows that integrates a **full detailed or reduced chemistry mechanism** directly, without reliance on pre-tabulated flamelet tables. It extends `mmcPremixedFoam` by replacing the flamelet lookup with an ODE-based chemistry solver evaluated on each stochastic particle.

## Algorithm Overview

Each time step performs the following steps:

1. **Density update** — solve continuity equation (`rhoEqn`)
2. **Momentum** — solve `UEqnTurbForcing`
3. **Diffusivity** — compute molecular and turbulent diffusivities
4. **ATF properties** — compute ATF model quantities (`calcATFproperties`)
5. **Efficiency constant** — compute local efficiency function (`calcEffConst`)
6. **Enthalpy / species** — solve directly without flamelet coupling (`hYEqvE_EqnNoCoupling`)
7. **Progress variable** — set progress variable field (`setProgressVar`)
8. **Pressure correction** — PIMPLE outer correctors
9. **Particle evolution** — advance stochastic particles; sample chemistry directly

## Key Physics / Models

| Component | Model |
|-----------|-------|
| Turbulence | Compressible LES |
| Thermodynamics | `psiReactionThermo` |
| Reactive scalars | Premixed stochastic particle cloud (`basicPremixedReactingPopeCloud`) |
| Chemistry | Full mechanism via `BasicChemistryModel` / `chemistrySolver` |
| Turbulence–flame | ATF model with efficiency function |

## Combustion Mode

Premixed with direct (online) kinetics. No flamelet table is needed; species and enthalpy evolve according to the full kinetic mechanism.

## Comparison with mmcPremixedFoam

| Feature | mmcPremixedFoam | mmcFullChemPremixedFoam |
|---------|-----------------|-------------------------|
| Chemistry | Flamelet table (offline) | Full ODE integration (online) |
| Cost | Lower | Higher |
| Accuracy | Table-limited | Full mechanism fidelity |
| Flamelet table required | Yes | No |

## Usage

Use this solver when:

- High accuracy kinetics are required for premixed flames
- The chemistry mechanism cannot be adequately tabulated
- Computational cost of online ODE integration is acceptable

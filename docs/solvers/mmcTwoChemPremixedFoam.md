# mmcTwoChemPremixedFoam

`mmcTwoChemPremixedFoam` is an LES/FDF solver for turbulent **premixed** reacting flows with **full chemistry** integration, designed to support two chemistry configurations or dual-timescale chemistry approaches. It shares the same algorithmic structure as `mmcFullChemPremixedFoam` while accommodating an alternative or secondary chemistry model.

## Algorithm Overview

Each time step performs the following steps:

1. **Density update** — solve continuity equation (`rhoEqn`)
2. **Momentum** — solve `UEqnTurbForcing`
3. **Diffusivity** — compute molecular and turbulent diffusivities
4. **ATF properties** — compute ATF model quantities (`calcATFproperties`)
5. **Efficiency constant** — compute local efficiency function (`calcEffConst`)
6. **Enthalpy / species** — solve without flamelet coupling (`hYEqvE_EqnNoCoupling`)
7. **Progress variable** — set progress variable field (`setProgressVar`)
8. **Pressure correction** — PIMPLE outer correctors
9. **Particle evolution** — advance stochastic particles; sample chemistry

## Key Physics / Models

| Component | Model |
|-----------|-------|
| Turbulence | Compressible LES |
| Thermodynamics | `psiReactionThermo` |
| Reactive scalars | Premixed stochastic particle cloud (`basicPremixedReactingPopeCloud`) |
| Chemistry | Full mechanism via `BasicChemistryModel` / `chemistrySolver` |
| Turbulence–flame | ATF model with efficiency function |

## Combustion Mode

Premixed with direct (online) kinetics. Supports configurations with two distinct chemistry mechanisms or timescales applied to the same flow field.

## Comparison with mmcFullChemPremixedFoam

Both solvers use the same time-step algorithm. The key distinction is support for a secondary or alternative chemistry configuration, enabling:

- Dual-mechanism approaches (e.g. fast and slow reactions)
- Comparative studies with two different reduced mechanisms
- Two-timescale modelling for separated fast/slow chemistry

## Usage

Use this solver when:

- Two chemistry models or mechanisms need to be active simultaneously
- Dual-timescale premixed combustion modelling is required
- You are extending `mmcFullChemPremixedFoam` with a secondary chemistry pathway

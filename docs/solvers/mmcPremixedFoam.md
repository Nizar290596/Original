# mmcPremixedFoam

`mmcPremixedFoam` is an LES/FDF solver for turbulent **premixed** reacting flows using the Multiple Mapping Conditioning (MMC) method. Chemistry is handled via pre-tabulated flamelet lookup tables, and turbulence–chemistry interaction is modelled through the Algebraic Turbulent Flame (ATF) model.

## Algorithm Overview

Each time step performs the following steps:

1. **Density update** — solve continuity equation (`rhoEqn`)
2. **Momentum** — solve `UEqnTurbForcing` (includes turbulent forcing term)
3. **Diffusivity** — compute molecular and turbulent diffusivities
4. **Flamelet closure** — compute progress variable source term (`flameletClosure_pvSource`)
5. **ATF properties** — compute ATF model quantities (`calcATFproperties`)
6. **Scalar transport** — solve progress variable and mixture fraction equations (`XiEqnATF`)
7. **Enthalpy / species update** — retrieve thermochemical state from flamelet table (`flamelet_hYEqvE_Z42`)
8. **Pressure correction** — PIMPLE outer correctors
9. **Particle evolution** — advance stochastic particles; sample statistics

## Key Physics / Models

| Component | Model |
|-----------|-------|
| Turbulence | Compressible LES (`turbulentFluidThermoModel`) |
| Thermodynamics | `psiReactionThermo` |
| Reactive scalars | Premixed stochastic particle cloud (`basicPremixedReactingPopeCloud`) |
| Chemistry | Pre-tabulated flamelet table (`interpolateFlameletTable`) |
| Turbulence–flame | Algebraic Turbulent Flame (ATF) model |
| Mixing time | Tabulated mixing time model |

## Combustion Mode

Premixed. A progress variable `c` tracks the flame front; mixture fraction `Z` accounts for mixture inhomogeneity.

## Control Settings

The flamelet table file and ATF model parameters are specified in `constant/cloudProperties`. The mixing time table is loaded from the location defined in `cloudProperties`.

## Usage

Use this solver when:

- The flame is premixed or stratified premixed
- Chemistry cost is critical and pre-tabulated flamelets are acceptable
- ATF closure is appropriate for the resolved flame–turbulence interaction

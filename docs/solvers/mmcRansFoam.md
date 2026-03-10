# mmcRansFoam

`mmcRansFoam` is a **RANS-PDF** (Reynolds-Averaged Navier–Stokes / Probability Density Function) solver for turbulent non-premixed reacting flows using the Multiple Mapping Conditioning (MMC) method. It is the RANS counterpart to the LES-based `mmcFoam`, using a k–ε turbulence closure and computing the turbulent time scale explicitly.

## Algorithm Overview

Each time step performs the following steps:

1. **Density update** — solve continuity equation (`rhoEqn`)
2. **Momentum** — solve `UEqn`
3. **Diffusivity** — compute molecular and turbulent diffusivities
4. **Turbulence time scale** — compute `tauTurb = k / epsilon`
5. **Scalar transport** — solve mixture fraction (`XiEqn`) and enthalpy/species (`hYEqvE_Eqn`)
6. **Pressure correction** — PIMPLE outer correctors
7. **Particle evolution** — advance MMC stochastic particles using the RANS turbulent time scale

## Key Physics / Models

| Component | Model |
|-----------|-------|
| Turbulence | RANS k–ε (`turbulentFluidThermoModel`) |
| Thermodynamics | `psiReactionThermo` |
| Reactive scalars | Stochastic Pope particle cloud (`basicReactingPopeCloud`) |
| Mixing time scale | `tauTurb = k / epsilon` |
| Pressure–velocity | PIMPLE |

## Combustion Mode

Non-premixed (diffusion flame), same as `mmcFoam`.

## Comparison with mmcFoam

| Feature | mmcFoam (LES) | mmcRansFoam (RANS) |
|---------|--------------|---------------------|
| Turbulence closure | LES subgrid model | k–ε |
| Scales resolved | All except subgrid | None (fully modelled) |
| Computational cost | Higher | Lower |
| Suitable for | Research, complex geometries | Industrial configurations |

## Usage

Use this solver when:

- Computational cost of LES is prohibitive
- RANS accuracy is sufficient for the target quantities
- Large-scale industrial geometries are simulated
- Steady-state or time-averaged statistics are the primary goal

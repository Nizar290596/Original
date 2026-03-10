# SPFoam

`SPFoam` (Shadow Position Foam) is an LES/FDF solver for turbulent reacting flows that implements an advanced turbulence–flame interaction closure based on a **lambda-function** approach. It computes the flame wrinkling factor `vb` dynamically from local flow conditions using explicit Laplacian-based strain rate estimation, providing a physics-based subgrid flame–turbulence model.

## Algorithm Overview

Each time step performs the following steps:

1. **Density update** — solve continuity equation (`rhoEqn`)
2. **Momentum** — solve `UEqn`
3. **Laplacian evaluation** — compute the Laplacian field of velocity (`evaluateLaplacian`)
4. **Velocity fluctuation from strain rate** — `uPrime = 2 * |delta^3 * curl(Laplacian(U))|`
5. **Velocity from TKE** — `uPrimek = sqrt(2/3 * k)`
6. **Lambda-function calculation**:
   - `RL = delta / delta_l0` — length scale ratio
   - `RV = uPrime / s_l0` — velocity ratio
   - Compute closure functions: `a_vb`, `fu`, `fdelta`, `Re_sgs`, `fre`, `gammafit`
7. **Wrinkling factor** — compute `vb` from composite lambda function
8. **Scalar transport** — solve mixture fraction and enthalpy/species with `vb`-modified closure
9. **Pressure correction** — PIMPLE outer correctors
10. **Particle evolution** — advance stochastic particles

## Key Physics / Models

| Component | Model |
|-----------|-------|
| Turbulence | Compressible LES |
| Reactive scalars | MMC stochastic Pope particle cloud |
| Flame–turbulence interaction | Lambda-function closure |
| Strain rate | Explicit Laplacian-based calculation |
| Wrinkling factor | Dynamic `vb` from velocity and length scale ratios |

## Lambda-Function Closure

The lambda-function closure accounts for both turbulent velocity effects and length scale effects on the subgrid flame surface density:

- **`fu`** — turbulent velocity contribution to wrinkling
- **`fdelta`** — flame thickness / length scale contribution
- **`fre`** — subgrid Reynolds number effect
- **`vb`** — resultant flame wrinkling factor applied in the scalar transport equation

## Flame Parameters

The closure requires laminar flame speed `s_l0` and laminar flame thickness `delta_l0`, which are hard-coded for the target configuration. These values must be set to match the fuel/oxidiser combination and thermodynamic conditions of the case.

## Relationship to SPFoamDRM22 and SPFoamDarmstadt

`SPFoam` is the generic version of the Shadow Position solver. `SPFoamDRM22` and `SPFoamDarmstadt` are pre-configured variants with flame parameters tuned to specific experimental burners (Aachen and Darmstadt swirl combustors, respectively).

## Usage

Use this solver when:

- A physics-based dynamic flame wrinkling closure is required
- The lambda-function model is appropriate for the target flame regime
- You are setting up a custom burner and will specify your own `s_l0` and `delta_l0`

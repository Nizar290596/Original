# mmcDNSFoam

`mmcDNSFoam` is a **Direct Numerical Simulation (DNS)** solver that couples full-resolution DNS with the MMC stochastic particle framework. It simultaneously computes the fully resolved DNS solution and extracts filtered fields on a coarser mesh, enabling the MMC particles to operate as LES-scale surrogates. This allows direct comparison between the particle-based FDF model and the exact DNS-filtered fields.

This solver extends the description in the original [mmcDNSFoam documentation](../mmcDNSFoam.md).

## Algorithm Overview

Each time step performs the following steps:

1. **Time synchronisation** — set the filter mesh time step equal to the DNS time step; advance both clocks
2. **Density update** — solve continuity equation (`rhoEqn`)
3. **Momentum** — solve `UEqn` on full DNS mesh
4. **Diffusivity** — compute molecular diffusivity (`D = nu/Sc`); no turbulent model
5. **Mixture fraction** — solve `XiEqn` on DNS mesh
6. **Species** — solve full species transport equations (`YEqn`)
7. **Energy** — solve full energy equation (`EEqn`)
8. **Pressure correction** — `pcEqn` (consistency mode) or `pEqn`
9. **Filter mesh update** — project DNS fields onto the coarser filter mesh
10. **Particle evolution** — advance MMC stochastic particles using filtered fields
11. **Output** — write both DNS and filtered fields

## Key Physics / Models

| Component | Model |
|-----------|-------|
| Turbulence | None — all scales fully resolved (DNS) |
| Thermodynamics | `psiReactionThermo` |
| Chemistry | Full mechanism via `CombustionModel` |
| Reactive scalars | Full species equations (`YEqn`) + energy (`EEqn`) |
| Particle cloud | MMC stochastic particles as passive LES surrogates |
| Filtering | DNS → filter mesh projection |

## Dual-Mesh Architecture

Two meshes are active simultaneously:

- **DNS mesh** — high-resolution mesh resolving all flow scales; full governing equations are solved here
- **Filter mesh** — coarser LES-equivalent mesh; DNS fields are filtered onto this mesh and used to drive the MMC particle cloud

## Control Settings

Filtering is controlled in `system/filterMesh/fvSolution`:

```
PDFMETHOD
{
    filterDNSFields true;   // Switch filtering on or off
}
```

When `filterDNSFields true`, the `cloudProperties` file and thermodynamic property files must be placed in `constant/filterMesh/`.

## Local Time Stepping

This solver supports Local Time Stepping (LTS) mode for non-reactive DNS simulations. It is enabled through the standard OpenFOAM LTS settings.

## Usage

Use this solver when:

- High-fidelity DNS reference data is needed to validate MMC LES particle models
- Filtered DNS fields are needed to assess the accuracy of the FDF closure
- A ground-truth comparison dataset for a turbulent reacting flow is required

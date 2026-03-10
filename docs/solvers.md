# mmcFoam Solvers Overview

This page provides an overview of all solvers available in the mmcFoam framework. mmcFoam implements the **Multiple Mapping Conditioning (MMC)** method — a stochastic Lagrangian particle approach for simulating turbulent reacting flows. Reactive scalars are carried by stochastic particles whose positions in composition space are conditioned on a reference field, enabling accurate representation of turbulence–chemistry interactions.

The **Shadow Position** (SP) concept is a special variant in which each stochastic particle carries a "shadow" position that tracks its trajectory in a reference space, improving localness of mixing in the MMC framework.

---

## Solver Index

| Solver | Type | Combustion Mode | Chemistry | Spray |
|--------|------|-----------------|-----------|-------|
| [mmcFoam](#mmcfoam) | LES/FDF | Non-premixed | Particle ODE | No |
| [mmcPremixedFoam](#mmcpremixedfoam) | LES/FDF | Premixed | Flamelet table | No |
| [mmcFullChemPremixedFoam](#mmcfullchempremixedfoam) | LES/FDF | Premixed | Full ODE | No |
| [mmcTwoChemPremixedFoam](#mmctwochempremixedfoam) | LES/FDF | Premixed | Dual full ODE | No |
| [mmcRansFoam](#mmcransfoam) | RANS/PDF | Non-premixed | Particle ODE | No |
| [mmcPbeFoam](#mmcpbefoam) | LES/FDF | Non-premixed | Particle ODE | No (aerosol) |
| [mmcDropletSprayFoam](#mmcdropletsprayfoam) | LES/FDF | Spray | Particle ODE | Yes |
| [mmcSpraySynFoam](#mmcspraysynfoam) | LES/FDF | Spray | Particle ODE | Yes |
| [SPFoam](#spfoam) | LES/FDF | Non-premixed | Particle ODE | No |
| [SPFoamDRM22](#spfoamdrm22) | LES/FDF | Non-premixed | DRM22 | No |
| [SPFoamDarmstadt](#spfoamdarmstadt) | LES/FDF | Non-premixed | Particle ODE | No |
| [mmcDNSFoam](#mmcdnsfoam) | DNS | Non-premixed | Full ODE | No |
| [mmcDropletSprayDNSFoam](#mmcdropletspraydnsfoam) | DNS | Spray | Full ODE | Yes |

---

## LES / FDF Solvers

### mmcFoam

The base MMC solver for turbulent **non-premixed** (diffusion) flames using LES. Couples an Eulerian compressible LES flow field with a Lagrangian stochastic Pope particle cloud that carries the reactive scalars (species, enthalpy) and performs mixing and reaction.

- **Turbulence:** LES
- **Thermodynamics:** `psiReactionThermo`
- **Particle cloud:** `basicReactingPopeCloud`
- **Pressure–velocity:** PIMPLE

[Full documentation →](./solvers/mmcFoam.md)

---

### mmcPremixedFoam

MMC solver for turbulent **premixed** flames. Chemistry is handled via pre-tabulated flamelet lookup tables. Turbulence–flame interaction uses the Algebraic Turbulent Flame (ATF) model, and a progress variable tracks the flame front.

- **Turbulence:** LES
- **Particle cloud:** `basicPremixedReactingPopeCloud`
- **Chemistry:** Flamelet table (`interpolateFlameletTable`)
- **Closure:** ATF model

[Full documentation →](./solvers/mmcPremixedFoam.md)

---

### mmcFullChemPremixedFoam

MMC solver for turbulent **premixed** flames with **full kinetics** solved online via ODE integration on each particle. No flamelet table is required. Uses the ATF model for turbulence–flame interaction.

- **Turbulence:** LES
- **Particle cloud:** `basicPremixedReactingPopeCloud`
- **Chemistry:** Full mechanism via `BasicChemistryModel`
- **Closure:** ATF model + efficiency function

[Full documentation →](./solvers/mmcFullChemPremixedFoam.md)

---

### mmcTwoChemPremixedFoam

Variant of `mmcFullChemPremixedFoam` supporting **two chemistry configurations** or dual-timescale chemistry for premixed combustion. Uses the same algorithm with an alternative or secondary chemistry model.

- **Turbulence:** LES
- **Chemistry:** Dual full mechanism
- **Closure:** ATF model + efficiency function

[Full documentation →](./solvers/mmcTwoChemPremixedFoam.md)

---

### mmcRansFoam

The **RANS** counterpart of `mmcFoam`. Uses a k–ε turbulence closure and computes the mixing time scale as `tau = k / epsilon`. Suitable for industrial-scale simulations where LES cost is prohibitive.

- **Turbulence:** RANS k–ε
- **Particle cloud:** `basicReactingPopeCloud`
- **Mixing time scale:** `k / epsilon`

[Full documentation →](./solvers/mmcRansFoam.md)

---

### mmcPbeFoam

Extends `mmcFoam` with **Population Balance Equations (PBE)** for tracking the size distribution of a dispersed particulate phase (aerosols, nanoparticles, soot). Exports normalised turbulence quantities (`nu`, `nut`, `epsilon`, `k`) for PBE coupling.

- **Turbulence:** LES
- **Particle cloud:** `basicAerosolReactingPopeCloud`
- **Additional physics:** Nucleation, growth, coagulation of dispersed particles

[Full documentation →](./solvers/mmcPbeFoam.md)

---

### mmcDropletSprayFoam

MMC solver for **spray combustion** with dilute liquid fuel droplets. Maintains two Lagrangian clouds simultaneously: gas-phase MMC stochastic particles and liquid-phase fuel droplets undergoing evaporation and heat transfer.

- **Turbulence:** LES
- **Gas-phase cloud:** `basicReactingPopeCloud`
- **Liquid cloud:** `basicDropletSprayThermoCloud`

[Full documentation →](./solvers/mmcDropletSprayFoam.md)

---

### mmcSpraySynFoam

Extended spray solver that adds comprehensive **spray diagnostics** (spray flux, D10, D32 Sauter mean diameter) and per-component performance timing. Designed for synthetic fuel spray combustion studies.

- **Turbulence:** LES
- **Gas-phase cloud:** Aerosol MMC particle cloud
- **Diagnostics:** `calcFlux`, `calcMeanDiameters` (D10, D32)

[Full documentation →](./solvers/mmcSpraySynFoam.md)

---

## Shadow Position (SP) Solvers

Shadow Position solvers implement an advanced turbulence–flame interaction closure based on a **lambda-function** approach. A local flame wrinkling factor `vb` is computed dynamically from the ratio of turbulent velocity and length scales to their laminar equivalents, providing a physics-based subgrid flame–turbulence model.

### SPFoam

The generic Shadow Position solver. Computes `vb` via the full lambda-function closure using an explicit Laplacian-based strain rate. Flame parameters (`s_l0`, `delta_l0`) must be specified by the user for the target configuration.

- **Closure:** Lambda-function (dynamic `vb`)
- **Strain rate:** Explicit Laplacian evaluation
- **Flame parameters:** User-specified

[Full documentation →](./solvers/SPFoam.md)

---

### SPFoamDRM22

Pre-configured variant of `SPFoam` for the **Aachen swirl combustor (SC-Aachen)** using the DRM22 reduced methane/air mechanism (`s_l0 = 0.409 m/s`, `delta_l0 = 3.62 × 10⁻⁴ m`).

[Full documentation →](./solvers/SPFoamDRM22.md)

---

### SPFoamDarmstadt

Pre-configured variant of `SPFoam` for the **Darmstadt swirl combustor (SC-Darmstadt)** with methane/air combustion (`s_l0 = 0.354 m/s`, `delta_l0 = 3.93 × 10⁻⁴ m`).

[Full documentation →](./solvers/SPFoamDarmstadt.md)

---

## DNS Solvers

DNS solvers fully resolve all turbulent scales and use no turbulence model. They serve as high-fidelity reference benchmarks for validating the LES/FDF MMC models. Both DNS solvers simultaneously compute the full DNS solution and project filtered fields onto a coarser LES-equivalent mesh for MMC particle validation.

### mmcDNSFoam

DNS solver for non-premixed turbulent reacting flows. Solves full species and energy equations on a DNS mesh, filters fields to a coarser mesh, and advances MMC stochastic particles using the filtered fields as surrogate LES inputs.

- **Turbulence:** None (all scales resolved)
- **Chemistry:** Full mechanism via `CombustionModel`
- **Dual mesh:** DNS mesh + filter mesh

[Full documentation →](./solvers/mmcDNSFoam.md)

---

### mmcDropletSprayDNSFoam

DNS solver for **two-phase spray combustion**. Extends `mmcDNSFoam` with a second Lagrangian cloud for DNS-resolved liquid fuel droplets. The highest-fidelity solver in the suite.

- **Turbulence:** None (all scales resolved)
- **Gas-phase cloud:** `basicReactingPopeCloud`
- **Droplet cloud:** `basicDropletSprayDNSThermoCloud`
- **Dual mesh:** DNS mesh + filter mesh

[Full documentation →](./solvers/mmcDropletSprayDNSFoam.md)

---

## Choosing a Solver

```
Is your flame premixed?
├── Yes
│   ├── Flamelet tables available → mmcPremixedFoam
│   ├── Full kinetics needed → mmcFullChemPremixedFoam
│   └── Two chemistry configs → mmcTwoChemPremixedFoam
└── No (non-premixed / spray)
    ├── Spray (liquid fuel)?
    │   ├── Production LES → mmcDropletSprayFoam / mmcSpraySynFoam
    │   └── DNS reference → mmcDropletSprayDNSFoam
    ├── Aerosol / nanoparticles → mmcPbeFoam
    ├── DNS reference → mmcDNSFoam
    ├── RANS (cost-limited) → mmcRansFoam
    ├── Lambda-function closure (custom burner) → SPFoam
    ├── SC-Aachen + DRM22 → SPFoamDRM22
    ├── SC-Darmstadt → SPFoamDarmstadt
    └── Default LES → mmcFoam
```

# SPFoamDarmstadt

`SPFoamDarmstadt` (Shadow Position Foam – Darmstadt) is a pre-configured variant of `SPFoam` with flame parameters tuned to the **Darmstadt swirl combustor** (SC-Darmstadt) for methane/air premixed combustion.

## Relationship to SPFoam

`SPFoamDarmstadt` uses the identical algorithm and lambda-function closure as `SPFoam`. The only difference is that the laminar flame parameters are fixed to values appropriate for the SC-Darmstadt burner operating conditions:

| Parameter | Value | Description |
|-----------|-------|-------------|
| `s_l0` | `3.54 × 10⁻¹ m/s` | Laminar flame speed |
| `delta_l0` | `3.93 × 10⁻⁴ m` | Laminar flame thickness |

These values are calibrated for methane/air combustion at the Darmstadt swirl combustor conditions.

## Algorithm Overview

Identical to `SPFoam`. See [SPFoam documentation](./SPFoam.md) for the full time-step algorithm description.

## Key Physics / Models

Same as `SPFoam`, targeting:

| Component | Detail |
|-----------|--------|
| Target burner | SC-Darmstadt swirl combustor |
| Fuel | Methane/air premixed |

## Comparison: Aachen vs Darmstadt Variants

| Parameter | SPFoamDRM22 (Aachen) | SPFoamDarmstadt |
|-----------|----------------------|-----------------|
| `s_l0` | 0.409 m/s | 0.354 m/s |
| `delta_l0` | 3.62 × 10⁻⁴ m | 3.93 × 10⁻⁴ m |
| Burner | SC-Aachen | SC-Darmstadt |

## Usage

Use this solver when:

- Simulating the Darmstadt swirl combustor (SC-Darmstadt) configuration
- Methane/air premixed combustion at Darmstadt operating conditions is used
- No manual tuning of `s_l0` and `delta_l0` is desired for this specific configuration

# SPFoamDRM22

`SPFoamDRM22` (Shadow Position Foam – DRM22) is a pre-configured variant of `SPFoam` with flame parameters tuned to the **Aachen swirl combustor** (SC-Aachen) using the **DRM22** reduced chemical mechanism (22-species Directed Relation Graph mechanism for methane/air combustion).

## Relationship to SPFoam

`SPFoamDRM22` uses the identical algorithm and lambda-function closure as `SPFoam`. The only difference is that the laminar flame parameters are fixed to values appropriate for the SC-Aachen burner operating conditions:

| Parameter | Value | Description |
|-----------|-------|-------------|
| `s_l0` | `4.09 × 10⁻¹ m/s` | Laminar flame speed |
| `delta_l0` | `3.62 × 10⁻⁴ m` | Laminar flame thickness |

These values are calibrated for the DRM22 methane/air mechanism at the Aachen swirl combustor conditions.

## Algorithm Overview

Identical to `SPFoam`. See [SPFoam documentation](./SPFoam.md) for the full time-step algorithm description.

## Key Physics / Models

Same as `SPFoam` with the addition of:

| Component | Model |
|-----------|-------|
| Chemistry mechanism | DRM22 (22-species reduced methane/air) |
| Target burner | SC-Aachen swirl combustor |

## Usage

Use this solver when:

- Simulating the Aachen swirl combustor (SC-Aachen) configuration
- The DRM22 reduced mechanism is used for methane/air combustion
- No manual tuning of `s_l0` and `delta_l0` is desired for this specific configuration

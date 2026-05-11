# Second-conditioning MMC mixing — implementation reference

Status: snapshot at branch `claude/identify-branch-YYXe0`, HEAD = `a47cb3f`.

This file documents the **two mixing stages** executed inside one solver
time step when `secondConditioning.enabled` is `true`, how the
flag-based partition steers the work between them, which scalars are
mixed where, and which dictionary entries actually matter.

It is intended as a self-contained reference: a new reader should be
able to read this file alone and understand the per-step Lagrangian
mixing logic.

---

## 1. The four-step picture per particle, per Δt

```
moveParticles.H per Δt
├─ 1. inflowBoundary().inflow()
├─ 2. mixing().Smix()                       ← MMCcurl, 3-D shadow space
├─ 3. updatePhiReaction(Δt)                 ← W(φ) on every particle
├─ 4. updateOUProcess(Δt)                   ← ω_OU advance + phiModified refresh
├─ 5. secondCondMixing().Smix()             ← secondCondMMCcurl, 4-D (φ°,sP)
├─ 6. solve()                               ← transport + per-particle chemistry
└─ 7. Sreact()                              ← load-balanced chemistry
```

Steps 3–5 are gated by `secondCondMixingEnabled()`. With
`secondConditioning.enabled false` they are skipped and the solver
reproduces the original Aachen behaviour exactly.

The reaction-source W(φ) and OU process are documented in §4 and §5.
The mixing stages are §2 and §3.

---

## 2. Stage 1 — `MMCcurl::Smix()`

### 2.1 Reference space

3-D, the first-conditioning shadow positions

    XiR = (sPx, sPy, sPz)

driven by their own SDE / advection elsewhere in the solver
(`mmcVariablesDefinitions`). MMCcurl does **not** touch `XiR`.

### 2.2 Particle list

Built by `mixParticleModel::buildParticleList()`. **Every** particle in
the cloud is added — no filter. For each particle a
`eulerianFieldData` record carries:

* `XiR[0..2]` = `sPx`, `sPy`, `sPz` (from the particle itself)
* `position` (used only for the `dx` diagnostic, not for splitting)
* `D`, `D_t`, `D_eff`, `ΔE`, `μ`, `v_b` interpolated at the particle
  location (used by the aISO timescale)
* `|∇XiR_i|²` per axis (used only by the C&K timescale variant)
* a Gaussian rand `Rand()` (used by some pair-rejection variants)

### 2.3 Pair selection

`mixParticleModel::findPairs()` calls `KkdTreeLikeSearch`. At each
internal node:

1. For each axis `i` compute the normalised spread
   `dis_i = (max_i − min_i) / Xii_[i]`.
2. Pick the axis with the largest `dis_i`.
3. Sort the slice along that axis, split at the median, recurse.

Recursion exits when slice size ≤ 2 (sometimes 3). Each terminal slice
becomes a pair (or triple) in `particlePairs_`.

`Xii_` for this stage has **3 entries** populated from
`MMCcurlCoeffs.Xim_i`:

* `Xii_[0] = sPx_m`
* `Xii_[1] = sPy_m`
* `Xii_[2] = sPz_m`

### 2.4 Per-pair mixing — the flag-based branch

`MMCcurl::mixpair` computes an aISO timescale

    τ_p = (1/v_b,p) · (ΔE_p)² / (CE · (D_p + D_t,p))     // similarly τ_q

with `τ_mix = harmonic_mean(τ_p, τ_q)` (or `min`, controlled by
`meanTimeScale`). Then

    α = 1 − exp(−Δt / τ_mix)

and calls `particleType::mixProperties(p, q, α)`. That dispatches up
the particle inheritance chain to `MixingPopeParticle::mixProperties`,
which has the flag-conditional branch (HEAD = `a47cb3f`):

| `(f_p, f_q)` | What gets mixed |
| --- | --- |
| `(0, 0)` | `ParticleType::mixProperties` runs → **Y, T, hA** mixed by parent chain (Reacting/Aerosol/Soot/Thermo). φ left alone. |
| `(1, 1)` | IEM relaxation on **φ** only: `φ ← φ + α·(φ_avg − φ)`. Y, T, hA left alone. |
| mixed   | **no-op**. |

So at stage 1, **flag is a hard partition** — same-flag pairs mix one
set of scalars, opposite-flag pairs do nothing.

`phiModified` is **not** recomputed here. Step 4
(`updateOUProcess`) refreshes it after the mixing → reaction sequence,
just before stage 2.

---

## 3. Stage 2 — `secondCondMMCcurl::Smix()`

### 3.1 Reference space

4-D

    XiR = (phiModified, sPx, sPy, sPz)

Particles' `XiR()` field stays 3-D (the shadow positions). The
`eulerianFieldData::XiR()` used by the k-d tree is **resized to 4 in
buildParticleList** with `phiModified` prepended.

### 3.2 Particle list

`secondCondMMCcurl::buildParticleList()` filters to
`secondCondFlag() == 1` only. Same Eulerian fields as stage 1 are
interpolated.

### 3.3 Pair selection

Same `KkdTreeLikeSearch` algorithm as stage 1 but with `Xii_` resized
to 4 in the `secondCondMMCcurl` constructor (commit `ecca6e7`) and
populated from `secondCondMMCcurlCoeffs.Xim_i`:

* `Xii_[0] = phiMod_m`
* `Xii_[1] = sPx_m`
* `Xii_[2] = sPy_m`
* `Xii_[3] = sPz_m`

**Parallelism:** local pairing only. No cross-rank pair exchange.

### 3.4 Per-pair mixing

`secondCondMMCcurl::mixpair` uses the same aISO timescale formula as
stage 1, then calls `mixSpeciesOnly(p, q, α)`. That routine relaxes:

* `hA` (enthalpy)
* `T` (temperature)
* `Y[i]` (species mass fractions)

toward the mass-weighted pair mean by extent `α`. **φ is not touched.**
`XiR`, `XiC`, `secondCondFlag` are not touched.

---

## 4. The reaction source `W(φ)` (step 3)

Applied to **every** particle by `MixingPopeCloud::updatePhiReaction`
(`MixingPopeCloud.C:310-334`):

    W(φ) = A_phi · (1 − φ) · exp[−Z_phi · (1 − φ)]
    φ ← clamp(φ + Δt · W(φ), 0, 1)

Forward Euler. Peak rate at `φ = 1 − 1/Z_phi` of magnitude
`A_phi / (e · Z_phi)`. Time-scale at peak `≈ e · Z_phi / A_phi`.

`phiModified` is set to `φ` for non-flagged particles inside the same
routine; flagged particles' `phiModified` is refreshed by step 4.

---

## 5. The OU process (step 4)

`MixingPopeCloud::updateOUProcess` (`MixingPopeCloud.C:338-362`).
For each particle with `secondCondFlag == 1`, runs the exact OU update

    ω_OU(t+Δt) = ω_OU · exp(−Δt/τ_OU) + √(1 − exp(−2Δt/τ_OU)) · ξ,
    ξ ~ N(0, 1)

then refreshes

    phiModified = φ · exp(β · ω_OU)

so the new value goes into stage 2's k-d tree.

The OU stationary distribution is N(0, 1), so non-flagged particles
keep ω_OU = 0 (set in the default constructor, also re-initialised
from N(0,1) at particle creation for flagged particles — commit
`70a589d`).

---

## 6. Per-particle behaviour table per Δt

For the dictionary in the user's case (`R = 1`, `enabled true`):

| `secondCondFlag` | Step 2 (MMCcurl) | Step 3 W(φ) | Step 4 OU | Step 5 (secondCondMMCcurl) | Step 6 chemistry |
| --- | --- | --- | --- | --- | --- |
| `1` (everyone, since R = 1) | **φ** mixed (3-D shadow, flagged↔flagged) | applied | `ω_OU` advances, `phiModified` refreshed | **Y, T, hA** mixed (4-D φ°+shadow, flagged↔flagged) | runs |
| `0` (none with R = 1) | — | — | — | — | — |

For the general `0 < R < 1` design:

| `secondCondFlag` | Step 2 mixes | Step 3 | Step 5 mixes | Chemistry |
| --- | --- | --- | --- | --- |
| `0` | Y, T, hA (3-D shadow, with other non-flagged) | W(φ) | not in pool | **skipped** |
| `1` | φ (3-D shadow, with other flagged) | W(φ) | Y, T, hA (4-D φ°+shadow, with other flagged) | runs |

Mixed-flag pairs in step 2 are **no-op** by design.

---

## 7. Dictionary entries that actually matter

### `cloudProperties.secondConditioning`

| key | effect | typical range / sanity |
| --- | --- | --- |
| `enabled` | turns on the whole second-conditioning machinery (steps 3–5, the flag-based partition in step 2, chemistry filter, OU process) | `true` to enable |
| `R` | per-particle Bernoulli probability of being flagged. **`R = 1` ⇒ every particle is flagged ⇒ MMCcurl never mixes Y, T, hA in step 2** (no `(0,0)` pairs ever). All Y/T/h mixing then happens only in step 5, on the 4-D reference space. | `0.1–0.5` for true second-conditioning; `1.0` collapses to "second-conditioning only" |
| `Tu`, `Tb` | bounds for `phi = (T − Tu)/(Tb − Tu)` at particle creation | actual unburnt / burnt temperatures of the case |
| `A_phi` | BML rate amplitude in `W(φ) = A·(1−φ)·exp[−Z·(1−φ)]` | for τ_chem ≈ τ_res / few, with `Z = 5–7`: `A ≈ 10 000–15 000` |
| `Z_phi` | BML activation parameter | typically `5–10`. **`Z = 20` is essentially "no ignition possible from φ = 0"** because `W(0) = A·exp(−Z)` is exponentially suppressed. |
| `beta` | OU shadow strength in `phiModified = φ · exp(β · ω_OU)` | `0.1–0.5`. Larger β = wider phi° distribution |
| `tauOU` | OU autocorrelation time | typically `O(Δt · 10–100)`; longer = stronger memory |

### `cloudProperties.subModels.MMCcurlCoeffs`

| key | effect |
| --- | --- |
| `r_i` | position-space normaliser; **not used by the k-d split** (k-d uses `XiR` only) — only consumed by older pair-rejection paths and `mixingSubVolumes` |
| `Xim_i.{sPx_m, sPy_m, sPz_m}` | 3-D normalisers, define which shadow axis wins each split (only ratios matter) |
| `CE` | mixing-rate coefficient in aISO timescale `τ = (1/v_b)·ΔE²/(CE·(D + D_t))`. **Larger CE = faster mixing.** Default `0.1` in code. `CE = 20` (current dict) → 200× faster mixing than default |
| `aISO` | `true` selects the aISO timescale; `false` selects the C&K timescale (`tauPi = CL·CE·β·(XiR_p − XiR_q)² / (|∇XiR|²·D_eff)`) |
| `meanTimeScale` | `true` = harmonic mean of `τ_p, τ_q`; `false` = `min(τ_p, τ_q)` |
| `pairingMethod` | `local` (no cross-rank pairs) or `global` |
| `fLowMixOnMaster`, `fHighMixOnMaster` | unused under aISO – legacy mixture-fraction filter |
| `localnessLimited`, `fullSort` | legacy |
| `C_E` *(note underscore)* | typo, **not read by the code**. Only `CE` is consumed |

### `cloudProperties.subModels.secondCondMMCcurlCoeffs`

Same keys as `MMCcurlCoeffs` plus:

| key | effect |
| --- | --- |
| `Xim_i.phiMod_m` | normaliser for the φ° axis in step 2's 4-D k-d split |
| `Xim_i.{sPx_m, sPy_m, sPz_m}` | normalisers for the three shadow axes in stage 2 |

---

## 8. Audit of the supplied `cloudProperties`

Issues, ordered by likely impact on "no burning":

### 8.1 `thermoPhysicalCouplingModel none`  ← *most likely the smoking gun*

With this entry the FV density and temperature fields **never receive
particle data**. Whatever the particles do — burn, mix, heat up — the
Eulerian side keeps its initial unburnt density. No thermal expansion,
no flow acceleration, no flame-driven flow response. The visual FV
flame can therefore never appear, even if particles are reacting
correctly.

Fix: set to `ParticleInCell` (already optimised for flagged-particle
filtering and super-cell bucketing in commits `5be08f6` and `f606c0a`):

```
thermoPhysicalCouplingModel     ParticleInCell;
```

### 8.2 `R 1` ← changes the meaning of the model

With every particle flagged, MMCcurl in step 2 only ever sees
flagged-flagged pairs, which under HEAD design mix only **φ**. So
**Y, T, hA never get mixed in step 2**. They are mixed exclusively in
step 5 (4-D conditioning).

That is not wrong — it's just a different model (every particle gets
4-D conditioning for species). But it means:

* `MMCcurlCoeffs.{sPx_m, sPy_m, sPz_m, CE, ...}` parameters now control
  only φ-mixing.
* The Aachen-validated 3-D shadow mixing of species no longer
  happens at all.
* All species/T mixing is at the resolution set by
  `secondCondMMCcurlCoeffs.Xim_i`.

If the intent is "every particle is treated by the refined 4-D model",
this is consistent — but recognise it's not the same as the dual-pool
design we implemented.

If the intent is dual-pool (refined treatment for a subset, standard
treatment for the rest), set `R = 0.25–0.5` instead.

### 8.3 `A_phi 1000, Z_phi 20` ← too slow / too stiff

* `τ_chem(peak) = e · Z / A = e · 20 / 1000 ≈ 54 ms`
* Residence time `τ_res = 0.18 m / 30 m/s ≈ 6 ms`

So the chemical timescale is **9× the residence time**. With
`Z = 20` the reaction is essentially zero except in a thin layer at
`φ ≳ 0.95` (`W(0) = 1000 · e⁻²⁰ ≈ 2·10⁻⁶ /s`). Flagged particles
entering at `φ = 0` cannot reach the W-active layer in the time they
spend in the domain.

Recommended:

```
A_phi   12000;
Z_phi   5;
```

Gives `W_peak ≈ 880 /s`, `τ_chem(peak) ≈ 1.1 ms`, and
`W(0) = 12000 · e⁻⁵ ≈ 80 /s` so cold particles can still ignite.

### 8.4 `pairingMethod local` for MMCcurl

With local-only pairing, particles can only mix with other particles
owned by the same MPI rank. For coarse decompositions this is fine; for
fine decompositions you can starve mixing across processor boundaries
and get "patchy" stripes of unmixed scalar. If you have many ranks,
consider:

```
pairingMethod   global;
```

(Second-stage `pairingMethod local` is **correct** and should stay —
the second-conditioning model is explicitly documented as local-only.)

### 8.5 `CE 20`

This is OK — it makes step-2 mixing 200× faster than the in-code
default `0.1`. Just be aware that step 5 is also using `CE 20`
(same coefficient in the same aISO formula), so the two stages are
balanced to each other.

If you want to slow step 2 to let chemistry breathe more, raise `CE`
in step 5 only (counter-intuitive: larger `CE` here ⇒ shorter `τ` ⇒
faster mixing, so to *slow* it you'd lower `CE` in step 5… or change
the formula. Not needed for now.)

### 8.6 `C_E 20` (with underscore)

This key is **not read** by the code. Either remove the line or
rename it to `CE` if you wanted it to take effect — but you already
have `CE 20` above, so this entry is dead.

### 8.7 `Xim_i` ratios

```
MMCcurlCoeffs.Xim_i           : sPx_m = sPy_m = sPz_m = 3.5e-5
secondCondMMCcurlCoeffs.Xim_i : phiMod_m = 3.5e-5, sP*_m = 3.5e-5
```

The second-stage ratios `phiMod_m : sP*_m = 1 : 1` are reasonable for
balanced 4-D conditioning. With your latest log this gives roughly
**27 % phi, 23 % xi_x, 23 % xi_y, 27 % xi_z**. Healthy.

If φ° mixing in stage 1 (since R=1) should be tighter, lower
`MMCcurlCoeffs.Xim_i.sPx_m, sPy_m, sPz_m`. But that's a tuning knob
for stage 1 only, not a "no-burn" fix.

### 8.8 Missing `f_m` / `mixture fraction` filter

`MMCcurlCoeffs.Xim_i` historically had an `f_m` entry; it's commented
out in your dict. With the current 4-D / 3-D `Xii_` setup `f_m` is
unused anyway. No action needed.

---

## 9. Recommended minimal patch to ignite

In order of impact:

1. **Switch on FV coupling.** Set `thermoPhysicalCouplingModel ParticleInCell;`
2. **Soften W(φ).** `A_phi 12000; Z_phi 5;`
3. **Lower R if you want a true partition.** `R 0.25` (or leave at 1
   if you intend the model to be "everyone gets 4-D").
4. Leave the `Xim_i` ratios alone — they look fine.
5. Leave `CE 20` alone.

After these four changes the flame should establish. If it still
doesn't, the next things to check are:

* The chemistry filter at `ReactingPopeParticle.C:63` — with `R = 1`
  every particle is flagged so this never blocks chemistry, but
  worth noting.
* Initial-temperature field — flagged particles inherit `phi` from
  `(T - Tu)/(Tb - Tu)`. If T is everywhere ≈ Tu at t = 0, every
  particle starts at φ = 0 and the only thing that breaks symmetry
  is `W(φ)`. With Z=5, A=12000 that's enough; with Z=20 it isn't.
* An explicit pilot region (force `phi = 1, T = Tb` in a small zone)
  if the chemistry is so stiff that even a softened W can't ignite
  from cold reactants alone.

---

## 10. Files that implement the above

| Function | File:line |
| --- | --- |
| `MixingPopeCloud::setParticleProperties` (flag assignment, ω_OU init, φ init) | `src/lagrangian/mmc/clouds/Templates/MixingPopeCloud/MixingPopeCloud.C:367-` |
| `MixingPopeCloud::updatePhiReaction` (W(φ) source) | `MixingPopeCloud.C:310-334` |
| `MixingPopeCloud::updateOUProcess` (OU step + φ° refresh) | `MixingPopeCloud.C:338-362` |
| `MMCcurl::mixpair` (stage 1 timescale + mixing dispatch) | `src/lagrangian/mmc/submodels/Mixing/MixingModel/ParticleInteractionModels/sparseParticleModel/MMCcurl/MMCcurl.C:81-289` |
| `MixingPopeParticle::mixProperties` (flag-conditional dispatch to Y/T/h or φ mixing) | `src/lagrangian/mmc/popeParticles/Templates/MixingPopeParticle/MixingPopeParticle.C:99-`*3 overloads* |
| `secondCondMMCcurl::buildParticleList` (flagged filter, 4-D record, k-d call) | `secondCondMMCcurl.C:70-275` |
| `secondCondMMCcurl::mixpair` and `mixSpeciesOnly` (stage 2 mixing) | `secondCondMMCcurl.C:278-376` |
| `mixParticleModel::KkdTreeLikeSearch` (the k-d split itself) | `mixParticleModel.C:862-982` |
| Diagnostic log / `postProcessing/secondCondMMCcurl.log` | `secondCondMMCcurl.C:160-275` |
| Bernoulli flag assignment | `MixingPopeCloud.C:393-394` |

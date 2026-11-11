# Homogenising baryon-fraction measurements

Two measurements can be compared only once five conventions match. In order of
how big the effect usually is:

| # | Convention | Typical size of the effect | Handled by |
|---|------------|----------------------------|------------|
| 1 | **Mass proxy**: hydrostatic vs. lensing vs. scaling relation | HSE bias b = 0.2 shifts M500 by **+37%** and f_gas by **−17%** | `apply_hse_bias(b)`, `filter(mass_proxy=…)` |
| 2 | **Overdensity / aperture**: r2500c, r500c, r200c, r200m, shells | M200c/M500c = 1.38–1.63; f_gas(<r200)/f_gas(<r500) = 1.13 (η=0.3) to 1.5 (η=1) | `to_delta(…, eta_gas=…, f_bcg=…)` |
| 3 | **Stellar IMF** | Salpeter → Chabrier: **×0.63** | `to_imf("chabrier")` |
| 4 | **Assumed cosmology** (mostly h) | h 0.70 → 0.677: M ×1.03, f_gas ×1.05, f_star ×1.03 | `to_cosmology("Planck18")` |
| 5 | **Content**: ICL included? satellites mass-limited? | up to ×1.5 on f_star in groups | metadata only (`stellar.includes_icl`) |

The rule in BaryCon is **no silent homogenisation**. Every conversion returns a new object
and appends a line to `.history`. The model-dependent conversions (2 and 1) have **no
default** for their key astrophysical parameter.

## 1. Cosmology at fixed observables

The observables are fluxes and surface brightness inside a fixed angle θ. At fixed redshift,
changing cosmology rescales the distances by a ≡ d_A^new/d_A^old (d_L scales the same way):

* hydrostatic mass M(<θ) ∝ d_A, and to first order the same holds for lensing and scaling-relation masses (exponent `mass.distance_exponent`, default 1);
* gas mass M_gas(<θ) ∝ d_A^{5/2} (emission measure ∝ n² V);
* stellar mass M★(<θ) ∝ d_L² (luminosity).

The overdensity radius must then be re-solved, because ρ_c(z) ∝ H(z)² also changes. With local
slopes M(<r) ∝ r^{α_M} (≈0.9 at r500 for NFW) and f(<r) ∝ r^η, the aperture ratio s = θ_new/θ_old
satisfies

    s^(3−α_M) = a^(p_M−3) · ρ_c^old(z) / ρ_c^new(z)

so that

    M_Δ → a^{p_M} s^{α_M} M_Δ ,   f_gas → a^{5/2−p_M} s^{η_gas} f_gas ,   f★ → a^{2−p_M} s^{η★} f★ .

For a change of h alone, s = 1 exactly, and the familiar scalings come out:
**M ∝ h⁻¹, f_gas ∝ h^{−3/2}, f★ ∝ h⁻¹**. 

Size of the effect (from `cosmology_factors`):

| change | z | mass | f_gas | f★ |
|---|---|---|---|---|
| (70, 0.3) → Planck18 | 0.1 | 1.033 | 1.051 | 1.034 |
| (70, 0.3) → Planck18 | 1.0 | 1.023 | 1.040 | 1.028 |
| Ω_m 0.30 → 0.35 | 0.1 | 0.993 | 0.993 | 0.996 |
| Ω_m 0.30 → 0.35 | 1.0 | 0.951 | 0.951 | 0.972 |

Ω_m matters only for high-z samples (SPT, eRASS1, Mantz+22), but the code handles it anyway.

**Units of the x-axis.** Masses are stored in physical M☉ for the authors' cosmology. With
`h_units=True` they are multiplied by the *target* h. That is valid only *after*
`to_cosmology`. Multiplying each paper's mass by its own h (as the legacy scripts did) is
equivalent for pure-h differences but is not the same operation.

## 2. Overdensity definitions and apertures

**Total masses** are converted under the NFW profile with a concentration–mass relation
(Duffy+08 built in; any `colossus` model, e.g. `diemer19`, via `cm_model=`). Where the paper
gives concentrations they can be passed with `c=`. The conversion depends on redshift and
mass:

| M500c [M☉] | z = 0 | z = 1 |
|---|---|---|
| 1e13 | M200c/M500c = 1.38 | 1.48 |
| 1e14 | 1.44 | 1.55 |
| 1e15 | 1.50 | 1.63 |

The legacy code used a single z = 0 table with no redshift dependence.

**Fractions** cannot be converted with the total-mass profile alone:

* **Gas** is less concentrated than dark matter, so f_gas(<r) rises with r. BaryCon assumes a
  local power law f_gas ∝ r^η, and **η must be given explicitly**. X-COP-like massive clusters
  have η ≈ 0.2–0.4 between r500 and r200. The r500 and r200 fits of Popesso+24 imply η ≈ 0.7
  at 1e15 M☉ and ≈ 1 at 1e14 M☉, because groups have very extended gas. Over r500 → r200 this
  gives a factor of 1.13 (η = 0.3) to 1.5 (η = 1). `eta_gas` also accepts a callable η(M).
* **Stars**: the BCG and its ICL are fully enclosed in both apertures, and satellites trace the
  total mass. Then f★,to = f★,from [f_BCG M_from/M_to + (1 − f_BCG)]. `fbcg_akino22(M)` gives
  f_BCG = 0.41, 0.17 and 0.07 at 1e13, 1e14 and 1e15 M☉, so f★(<r200)/f★(<r500) = 0.89, 0.95
  and 0.98.
* **Shells** (e.g. Mantz+22, 0.8–1.2 r2500c) are never converted.

The legacy stellar plot had a secondary right-hand axis mapping f★,500 → f★,200 computed at one
pivot mass (1e14). Because the mapping depends on mass, that axis was only correct at the pivot.
It has been removed in favour of an explicit `to_delta`.

## 3. Mass proxy

| proxy (`mass.proxy`) | entries | note |
|---|---|---|
| `hydrostatic` | Arnaud07, Vikhlinin06, Sun09, Ettori10, Zhang11, Landry13, Sanderson13, Mantz22, Wicker23 | biased low by (1 − b), b ≈ 0.1–0.3 |
| `weak_lensing` | Reyes12, Mulroy19, Akino22 | Mulroy uses posterior masses |
| `xray_scaling_relation` | Giodini09, Bulbul24, Kravtsov18 (Y_X), Gonzalez13 (T_X) | Giodini and Bulbul WL-calibrated; Y_X/T_X calibrated on HSE masses |
| `sz_scaling_relation` | Chiu18 | de Haan+16 ξ–M relation |
| `dynamical` | Khalil26 (to be confirmed) | σ_v-based? |
| `mixed` | Grandis24, Popesso24 | per-sample |

`apply_hse_bias(b)` moves both the mass *and the aperture*. For b = 0.2: M500 ×1.37,
r500 ×1.11, f_gas ×0.83, f★ ×0.80 (η_gas = 0.3).

**Correlated errors.** For hydrostatic and scaling-relation masses the aperture is defined by
the same mass, so the errors on M and f are correlated. Plotting them as independent bars
overstates the scatter. Fitting code should use the long-format table and model the
correlation.

## 4. IMF and stellar content

`IMF_TO_CHABRIER` uses Madau & Dickinson (2014): Salpeter ×0.63, Kroupa ×0.94, and
"diet Salpeter" ×0.90. Giodini+09 uses Salpeter, so its stellar fractions drop by 37%. Several
entries do not yet record their IMF (`needs_check`).

ICL: Gonzalez+13 and Kravtsov+18 include it; Chiu+18, Giodini+09 and Zhang+11 do not. Chiu+18
also integrates satellites only down to 1e10 M☉. These differences explain part of the factor
of ~2 offset between SPT and the low-z samples in `docs/figures/stellar_fractions_500c.png`.

## 5. Per-entry assumptions (first release)

| id | H0 | Ω_m | Δ | aperture | proxy | IMF | z | status |
|---|---|---|---|---|---|---|---|---|
| Arnaud2007 | 70 | 0.3 | 500c | sphere | HSE | – | 0.03–0.16 | unverified |
| Vikhlinin2006 | **71** | 0.3 | 500c, 2500c | sphere | HSE | – | 0.01–0.23 | unverified |
| Sun2009 | 73 | 0.24 | 500c | sphere | HSE | – | ≤0.12 | unverified |
| Giodini2009 | 72 | 0.258 | 500c | sphere | WL-cal. L_X | Salpeter | 0–1 | unverified (f_gas derived) |
| Ettori2010 | 70 | 0.3 | 500c | sphere | HSE | – | 0.09–0.31 | unverified |
| Zhang2011 | 70 | 0.3 | 500c | sphere | HSE | ? | 0.02–0.09 | unverified |
| Reyes2012 | 70? | 0.3? | 200c | sphere | WL | ? | ~0.1 | needs_check |
| Gonzalez2013 | 67.15? | 0.3175? | 500c | sphere | T_X SR | ? | 0.06–0.12 | needs_check |
| Landry2013 | 70? | 0.3? | 500c, 2500c | sphere | HSE | – | 0.15–0.3 | needs_check |
| Sanderson2013 | 70? | 0.3? | 500c, 2500c | sphere | HSE | – | ~0.1 | unverified |
| Chiu2018 | 68 | 0.304 | 500c | sphere | SZ SR | Chabrier | 0.28–1.29 | unverified |
| Kravtsov2018 | 70 | 0.27 | 500c | sphere | Y_X SR | Chabrier | 0.02–0.09 | needs_check |
| Mulroy2019 | 70 | 0.3 | 500c | sphere | WL | – | 0.15–0.31 | unverified |
| Akino2022 | 70? | 0.3? | 500c | sphere | WL | ? | 0–1 | needs_check |
| Mantz2022 | 70 | 0.3 | 2500c | **shell 0.8–1.2** | HSE | – | 0–1.2 | unverified |
| Wicker2023 | 70? | 0.3? | 500c | sphere | HSE | – | 0.06–0.55 | unverified |
| Bulbul2024 | 70? | 0.3? | 500c | sphere | WL-cal. SR | – | 0–1.3 | needs_check |
| Grandis2024 | per survey | | 500c | sphere | per row (SZ, WL-cal., WL) | ? | – | needs_check |
| Popesso2024 | 67.74 | 0.3089 | 500c, 200c | sphere | mixed | – | ≲0.2 | needs_check |
| Khalil2026 | 70? | 0.3? | 500c, 200c | sphere | dynamical? | – | 0.01–0.03 | needs_check |

**Bold** = differs from what the legacy scripts used. `?` = to be checked against the paper
(tracked as GitHub issues labelled `verify-metadata`).

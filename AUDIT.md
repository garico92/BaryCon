# Audit of the legacy codes and data (merged, 6 Nov 2026)

Each numerical fix is pinned by a regression test in `tests/test_barycon.py` and recorded in
the `caveats` or `changelog` of the entry. `docs/figures/legacy_vs_fixed.png` shows the size
of the largest ones.

## A. Errors that change the numbers

| Dataset | Problem → fix | Effect | Check here |
|---|---|---|---|
| Kravtsov+18 | The header was shifted by one column, so the code read 15.60 (A2142) as M500. M500 is the 11.96 ± 0.20 column, which agrees with R500 and z to <7%. The 15.60 column behaves like M_gas,500 in 1e13 M☉. The paper says M500 comes from the Y_X proxy. | f★ −18% to +30% (mass-dependent) | [arXiv:1401.7329 Table 1](https://arxiv.org/pdf/1401.7329), [text](https://ar5iv.arxiv.org/html/1401.7329) |
| Ettori+10 | NFW conversion of the Table 2 M200, which is the **dark** mass. Now M500 = 500 ρc(z) 4π/3 R500³, from R500 (Corrigendum) and z (Table 1). | M500 was 18% too low | [arXiv:1009.3266](https://arxiv.org/pdf/1009.3266), [ADS + Corrigendum](https://ui.adsabs.harvard.edu/abs/2010A&A...524A..68E) |
| Popesso+24 | The band was built at M = 1 M☉ on a fit with no pivot. It is now re-expressed at a 10¹⁴ M☉ pivot and sampled (σA assumed at the pivot). | Band about 2× too wide at 10¹⁵ | [doi:10.1051/0004-6361/202453256](https://doi.org/10.1051/0004-6361/202453256) |
| Reyes+12 | The columns are percentiles of log(M200c/M★) but were used as errors. | Error bars meaningless | [arXiv:1110.4107 Table 5](https://arxiv.org/pdf/1110.4107) |
| Vikhlinin+06 | h = 0.72 in the codes; the paper and the file use 0.71. | M −1.4%, f_gas −2.1% | [doi:10.1086/500288](https://doi.org/10.1086/500288) |
| Zhang+11 | f★ rescaled as h^−3/2 instead of h⁻¹. | −1.7% | `docs/homologation.md` §1 |
| Bulbul+24 | Masses are in units of 1e13 M☉ (header said 1e14). The "±" columns are bounds, not errors, so `M500_err` summed two bounds. h = 0.6715 was copied from Gonzalez. | ×10 in one script; wrong errors | [CDS ReadMe](https://cdsarc.cds.unistra.fr/ftp/J/A+A/685/A106/ReadMe) |
| Giodini+09 | f_gas is a fit to V06, A07 and S09, not a measurement; the stellar masses use a Salpeter IMF. | Double counting; f★ ×1.6 vs Chabrier | [arXiv:0904.0448](https://arxiv.org/pdf/0904.0448) |
| several | The common-h rescaling was inconsistent (Chiu, Mulroy, Akino, Grandis and Khalil not rescaled; f★ never rescaled). | 0.5–5% | `docs/homologation.md` §1 |

## B. Code and fitting

* `gas_frac_M_simple`: A and M_p are exactly degenerate, so the pivot is now fixed.
* `gas_frac_M_double`: it is a saturating function, so it is renamed `saturating` (alias kept). M_t is now fitted in log10; its error had been ~1e−16.
* Fits: optional intrinsic scatter and mass errors; numerical-Hessian errors; E(z) from the chosen cosmology. The fit refuses mixed apertures or shells and warns on mixed mass methods.
* The f★,500 → f★,200 axis was evaluated at a single pivot mass; removed in favour of `to_delta`.
* M500↔M200: the z = 0 concentration table is replaced by a z-dependent c–M relation (Duffy+08 or colossus).
* Popesso's mass was rescaled with the target h instead of the paper's h. The Akino band was built from corners and its f_gas not rescaled.

## C. Data files and metadata

* Mass methods are recorded: Kravtsov (Y_X), Gonzalez (T_X), Chiu (SZ; h = 0.68, Ωm = 0.304, Chabrier).
* Mantz+22 is a **shell** (0.8–1.2 r2500c): never converted or fitted with r500c data.
* Landry+13 (no header) and Sanderson+13 (unicode) are now usable. The Ettori method-1 file has 43 rows against 44.
* Names and redshifts were restored by Douspis (Sun, Vikhlinin, Gonzalez, Sanderson, Chiu, Kravtsov, Ettori). Ettori's names and redshifts were transcribed independently by both authors and are identical.
* Label years corrected: Arnaud 2007, Sun 2009, Wicker 2023, Ettori 2010 (+ 2011 Corrigendum).

## D. Still open (`needs_check`)

Kravtsov M_gas column (journal table) · Gonzalez cosmology and IMF · Bulbul catalogue cosmology,
z and flags · Khalil cosmology and mass method · cosmologies for Reyes, Akino, Landry, Sanderson
and Wicker · IMFs for Reyes, Akino and Zhang · Grandis homogenisation and title · Popesso
covariance · Chiu redshifts (Simbad vs SPT; `retrieved_redshifts.txt` missing).

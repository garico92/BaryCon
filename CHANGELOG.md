# Changelog

## Unreleased (2026-11-09)
- Plots: colour identifies the publication (fixed across figures), marker the mass-proxy family
  (hydrostatic, WL / WL-calibrated, X-ray or SZ scaling relation, dynamics).
- Grandis+24 stored as points (one per survey, at the pivot mass) with the local slope and a per-row
  mass proxy, instead of full bands; `Dataset.get` returns `slope*` and `proxy` columns when present.
- Popesso+24 reference updated to the published version (A&A 707, A362, 2026).

## 0.1.0 (unreleased, 2026-11-06)
- `BaryonFractions` dictionary interface re-implemented on the metadata core (`barycon.compat`);
  `baryon_fractions` import name kept as an alias.
- `barycon.fitting`: `power_law` (fixed pivot), `saturating`, `double_power_law_z` (aliases
  `gas_frac_M_simple/double/z`), optional intrinsic scatter and mass errors.
- Plots against redshift (`Collection.plot(..., x="z")`, `plot_gas_fractions_z`).
- New data: Khalil+26 (r500c and r200c); names and redshifts for Sun, Vikhlinin, Gonzalez,
  Sanderson, Chiu, Kravtsov; DOIs for all entries. Mass methods: Kravtsov (Y_X), Gonzalez (T_X).
- Raw text tables kept in `data/` (merged with M. Douspis's corrections); entries rebuilt with `tools/migrate_from_text.py`.

## 0.1.0 (unreleased) — data 2026.11.0
First packaged release, migrated from the `baryon_fractions` compilation.

### Data
- 23 entries / 19 papers in a machine-readable format (ECSV + YAML) with full metadata.
- Fixed: Kravtsov+18 M500 column; Ettori+10 total M500 (was the dark-matter mass); Popesso+24
  error band (no pivot); Reyes+12 percentiles; Vikhlinin+06 h = 0.71; Zhang+11 f★ h-scaling;
  Bulbul+24 units and bounds. See `AUDIT.md`.
- Added previously unusable entries: Landry+13, Sanderson+13 (r500 and r2500), Vikhlinin+06 r2500.
- Giodini+09 f_gas flagged as derived (excluded by default).

### Code
- `Dataset`, `Relation`, `Collection`; `to_cosmology`, `to_delta`, `to_imf`, `apply_hse_bias`;
  history logging; JSON-schema validation; CLI; plotting.

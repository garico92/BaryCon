# Contributing to BaryCon

New measurements, corrections and code are all welcome. The most valuable
contribution is **your own published measurement, entered by you**: nobody knows its
conventions better. Thanks for contributing to BaryCon!

## Adding a dataset

### Option A — no code (GitHub issue form)
Open *Issues → New → "New dataset"*. Paste the table (any format, or a link to
VizieR/journal), and answer the questions about cosmology, aperture and mass proxy.
A maintainer will create the entry and ask you to review it.

### Option B — pull request
1. Fork and clone; `pip install -e ".[test]"`.
2. Copy a template:
   * tabulated points or bins → `templates/points_template.ecsv`
   * a fitted relation f = A (M/M_piv)^B → `templates/relation_template.yaml`
3. Name the file `<FirstAuthor><Year>[_suffix].ecsv`, where the id equals the file name
   (e.g. `Eckert2019.ecsv`; use a suffix for a second aperture: `Eckert2019_200c`).
4. Fill in the metadata (see the checklist below) and the table.
5. `barycon validate src/barycon/data/observations/YourEntry.ecsv` and `pytest`.
6. Add the BibTeX to `src/barycon/data/references.bib` (key = id without suffix).
7. Open the PR. CI validates the schema; a maintainer reviews it.

### Conventions (please read)
* **Store values as published.** Do not rescale to another cosmology or h; record the authors' cosmology instead.
* Masses in physical M☉. ECSV columns may carry a unit (e.g. `1e14 solMass`); it is converted on load.
* Uncertainties: positive 1σ deltas in `<col>_lo` and `<col>_hi`. For other confidence levels, add extra columns and explain them in `description`.
* Column names: `M_<delta>` (e.g. `M_500c`), `f_gas`, `f_star`, `f_bcg`, `f_b`, optional `z`, `name`. Keep any other useful column (temperature, radius, gas mass…).
* One entry = one paper × one aperture. A second aperture gets a second file.
* `extracted: true` if numbers were read off a figure.

### Metadata checklist
- [ ] `cosmology`: H0, Om0 (and Ode0) exactly as assumed in the paper
- [ ] `aperture`: delta (`500c`, `200m`, …), region (`sphere`/`shell`), shell radii if any
- [ ] `mass.proxy`: hydrostatic / weak_lensing / *_scaling_relation / dynamical / … and `proxy_details`
- [ ] `stellar.imf`, `includes_icl` (stellar entries)
- [ ] `quantities`, plus `derived_quantities` for anything not directly measured
- [ ] `redshift` range; `selection`; `table_source` (table number, journal or arXiv version)
- [ ] `caveats`: anything a careful user must know
- [ ] `overlaps_with`: entries sharing objects or measurements (to avoid double counting)

## Status and verification
`unverified` (default) → `verified` once an original author or a second maintainer has
checked the entry against the paper (`verified_by`). `needs_check` flags a known open
question; describe it in `caveats` and open an issue.

## Code contributions
Follow PEP 8 (`black`, line length 110; `flake8`), add tests for new functionality, update
docstrings and README, and write clear commit messages. Development setup:

```bash
pip install -e ".[dev]"
pytest                       # tests
python -m build              # sdist + wheel (twine upload dist/* for PyPI)
```

Any conversion must (i) return a new object, (ii) append to
`.history` and (iii) have no hidden default for a model-dependent astrophysical parameter.

## Code of conduct
Be kind and constructive. We follow the [Contributor Covenant](https://www.contributor-covenant.org/version/2/1/code_of_conduct/).

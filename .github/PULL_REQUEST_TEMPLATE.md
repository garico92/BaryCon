## What does this PR do?
<!-- new dataset / correction / code -->

## For data PRs
- [ ] File name = `id` (`FirstAuthorYEAR[_suffix]`), values **as published** (no rescaling)
- [ ] `barycon validate <file>` passes and `pytest` passes
- [ ] cosmology, aperture (delta + region), mass proxy, redshift range filled in
- [ ] IMF / ICL recorded (stellar entries)
- [ ] `derived_quantities`, `overlaps_with`, `caveats` considered
- [ ] BibTeX added to `src/barycon/data/references.bib`
- [ ] I am an author of the measurement (→ can be `verified`) / I am not

## For code PRs
- [ ] conversions return new objects and append to `.history`
- [ ] no hidden defaults for model-dependent astrophysics
- [ ] tests added

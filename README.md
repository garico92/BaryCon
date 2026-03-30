
---

# Halo Gas and Stellar Fraction Datasets

This repository contains a collection of halo gas fractions, and in some cases, stellar fractions, compiled from the literature. These datasets are primarily derived from X-ray or SZ observations, with total halo masses estimated under the assumption of thermal equilibrium or calibrated with weak lensing measurements.

> **Disclaimer:** The tables and data are provided *as-is*, without any warranty. The idea is to have a fast tool to compare different observations. Users are responsible for verifying and interpreting the data appropriately. If you use any of these datasets in your research, please cite the original publication.

---

## Summary of Datasets

| Paper                 | Type          | Dataset specifics                                                                | Redshift Range | Halo Mass Range | Reference                                                                                                                 |
| --------------------- | ------------- | -------------------------------------------------------------------------------- | -------------- | --------------- | ------------------------------------------------------------------------------------------------------------------------- |
| Popesso et al. 2024   | Gas           | X-ray (eROSITA + CHEXMATE)                                                       | 0.0–0.5        | 10¹²–10¹⁵ M☉    | [arXiv:2411.16555](https://arxiv.org/pdf/2411.16555)                                                                      |
| Bulbul et al. 2024    | Gas           | X-ray (eROSITA)                                                                  | 0.0–1.3        | 10¹²–10¹⁵ M☉    | [arXiv:2402.08452](https://arxiv.org/pdf/2402.08452)                                                                      |
| Grandis et al. 2024   | Gas + Stellar | Optical + X-ray + SZ (compilation)                                               | 0.0–1.0        | 10¹³–10¹⁵ M☉    | [arXiv:2309.02920](https://arxiv.org/pdf/2309.02920)                                                                      |
| Wicker et al. 2023    | Gas           | SZ-selected (Planck) + X-ray (XMM-Newton)                                        | 0.0–0.5        | 10¹⁴–10¹⁵ M☉    | [arXiv:2204.12823](https://arxiv.org/pdf/2204.12823)                                                                      |
| Akino et al. 2022     | Gas + Stellar | X-ray (XXL/XMM-Newton) + Weak lensing (HSC) + Optical (SDSS photo)               | 0.0–1.0        | 10¹³–10¹⁵ M☉    | [arXiv:2111.10080](https://arxiv.org/pdf/2111.10080)                                                                      |
| Mantz et al. 2022     | Gas           | X-ray (Chandra) + WL (HSC/CFHT)                                                  | 0.0–1.2        | 10¹⁴–10¹⁵ M☉    | [arXiv:2111.09343](https://arxiv.org/pdf/2111.09343)                                                                      |
| Mulroy et al. 2019    | Gas           | X-ray (ROSAT/Chandra/XMM-Newton) + Weak lensing (HSC) + SZ (SZA/Planck) + NIR    | 0.1–0.3        | 10¹⁴–10¹⁵ M☉    | [arXiv:1901.11276](https://arxiv.org/pdf/1901.11276)                                                                      |
| Chiu et al. 2018      | Gas + Stellar | SZ-selected (SPT) + X-ray (Chandra) + Optical (DES) + NIR (WISE/Spitzer)         | 0.2–1.3        | 10¹⁴–10¹⁵ M☉    | [MNRAS](https://academic.oup.com/mnras/article/478/3/3072/4996803) / [arXiv:1711.00917](https://arxiv.org/pdf/1711.00917) |
| Kravtsov et al. 2018  | Stellar       | Optical (SDSS) + X-ray (Chandra + XMM-Newton)                                    | 0.0–1.0        | 10¹²–10¹⁵ M☉    | [arXiv:1401.7329](https://arxiv.org/pdf/1401.7329.pdf)                                                                    |
| Gonzalez et al. 2013  | Gas + Stellar | Optical + X-ray (XMM-Newton)                                                     | 0.0–0.3        | 10¹³–10¹⁵ M☉    | [arXiv:1309.3565](https://arxiv.org/pdf/1309.3565.pdf)                                                                    |
| Landry et al. 2013    | Gas           | X-ray (Chandra)                                                                  | 0.2–1.0        | 10¹³–10¹⁵ M☉    | [arXiv:1211.4626](https://arxiv.org/pdf/1211.4626.pdf)                                                                    |
| Sanderson et al. 2013 | Gas           | Optical + X-ray (XMM-Newton)                                                     | 0.0–0.2        | 10¹³–10¹⁵ M☉    | [arXiv:1212.1613](https://arxiv.org/pdf/1212.1613.pdf)                                                                    |
| Reyes et al. 2012     | Stellar       | Optical + Weak lensing (SDSS)                                                    | 0.0–0.3        | 10¹³–10¹⁵ M☉    | [arXiv:1110.4107](https://arxiv.org/pdf/1110.4107.pdf)                                                                    |
| Ettori et al. 2011    | Gas           | X-ray (XMM-Newton)                                                               | 0.1–0.3        | 10¹⁴–10¹⁵ M☉    | [arXiv:1009.3266](https://arxiv.org/pdf/1009.3266.pdf)                                                                    |
| Zhang et al. 2011     | Gas + Stellar | X-ray (XMM-Newton/Rosat) + Optical (SDSS photo)                                  | 0.0–0.1        | 10¹⁴–10¹⁵ M☉    | [arXiv:1109.0390](https://arxiv.org/pdf/1109.0390)                                                                        |
| Giodini et al. 2009   | Gas + Stellar | X-ray (XMM-Newton + Chandra) + Optical (COSMOS photo)                            | 0.1–1.0        | 10¹³–10¹⁵ M☉    | [arXiv:0904.0448](https://arxiv.org/pdf/0904.0448.pdf)                                                                    |
| Sun et al. 2008       | Gas           | X-ray (Chandra)                                                                  | 0.0–0.1        | 10¹³–10¹⁴ M☉    | [arXiv:0805.2320](https://arxiv.org/pdf/0805.2320.pdf)                                                                    |
| Arnaud et al. 2007    | Gas           | X-ray (XMM-Newton)                                                               | 0.0–0.5        | 10¹³–10¹⁵ M☉    | [A&A](https://www.aanda.org/articles/aa/pdf/2007/42/aa8541-07.pdf)                                                        |
| Vikhlinin et al. 2006 | Gas           | X-ray (Chandra)                                                                  | 0.0–0.3        | 10¹³–10¹⁵ M☉    | [Journal](https://iopscience.iop.org/article/10.1086/500288/pdf) / [arXiv](https://arxiv.org/pdf/astro-ph/0507092.pdf)    |



> Note: the redshift and halo mass ranges are approximate, based on the values reported in the papers and automatically compiled by AI.

---

## Usage

This rep
Each dataset is provided in ASCII files in the repository. Typical columns include:

* Halo mass estimates
* Gas fractions
* Stellar fractions
* Associated uncertainties

An example of how to load and plot the datasets is provided in **`plot_gas_fractions.py`**, **`plot_stellar_fractions.py`**, **`baryon_fractions.ipynb`**.

> **Warning:** The plotting routines focus on mass fractions as a function of the halo mass, but they do not take into account the different redshift of the observations. Both halo masses and baryonic fractions are not directly observable, and in their derivations by the different groups and datasets there are probably different biases and assumptions on astrophysics and cosmology not taken into account here.  

## Gas fractions 

![gas fractions](plots/gas_fractions.png)

## Stellar fractions 

![stellar fractions](plots/stellar_fractions.png)
---

## Citation

Please cite the original papers corresponding to the dataset(s) you use.

---

## License / Terms of Use

* The datasets in this repository are provided **as-is**, for research and educational purposes only.
* No warranty is given regarding accuracy, completeness, or suitability for any purpose.
* Users are responsible for verifying and interpreting the data reported.
* When using any of these datasets, please cite the original publication(s).
* Redistribution is allowed, provided this notice and citations are included.

---




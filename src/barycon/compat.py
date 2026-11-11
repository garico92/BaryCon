"""

Re-implementing BaryonFractions on top of the BaryCon entries so that existing notebooks keep
working, while the numbers come from the metadata core:

>>> from barycon import BaryonFractions
>>> bf = BaryonFractions()                       # Planck18 h, as before
>>> d = bf.gas_data["wicker_2022"]               # M500 [Msun/h], f_gas, f_gas_p, f_gas_m, f_gas_err, redshift
>>> popt, pcov = bf.fit(["wicker_2022", "arnaud_2007"], fit_type="gas_frac_M_double")
>>> fig, ax = bf.plot_gas_fractions_M(sources=["wicker_2022"], fit_sources=["akino_2022"])

Differences with the original package (all intentional, see ``AUDIT.md``):

* datasets are homogenised with the general cosmology rescaling (h *and* Omega_m),
  with the corrected Hubble exponents (f_star ∝ h^-1) and published h values;
* the Kravtsov+18, Ettori+10, Reyes+12 and Bulbul+24 inputs are corrected;
* ``f_gas_p``/``f_gas_m`` are always positive 1-sigma *deltas* in data dicts
  (Bulbul+24 used bounds), while in fit dicts they are the band edges as before;
* ``gas_frac_M_simple`` fixes the pivot (the 3-parameter version was degenerate);
* Giodini+09 gas fractions are kept for continuity but flagged as derived.
"""
from __future__ import annotations

import numpy as np
from astropy.cosmology import FlatLambdaCDM, Planck18
from scipy import optimize

from . import fitting as _fitting
from .collection import Collection, load
from .dataset import Dataset

# legacy key -> BaryCon id
GAS_KEYS = {
    "arnaud_2007": "Arnaud2007", "vikhlinin_2006": "Vikhlinin2006", "sun_2008": "Sun2009",
    "giodini_2009": "Giodini2009", "ettori_2011": "Ettori2010", "zhang_2011": "Zhang2011",
    "gonzalez_2013": "Gonzalez2013", "sanderson_2013": "Sanderson2013", "landry_2013": "Landry2013",
    "chiu_2018": "Chiu2018", "mulroy_2019": "Mulroy2019", "wicker_2022": "Wicker2023",
    "bulbul_2024": "Bulbul2024", "grandis_2024": "Grandis2024", "khalil_2026": "Khalil2026",
}
STELLAR_KEYS = {
    "giodini_2009": "Giodini2009", "zhang_2011": "Zhang2011", "reyes_2012": "Reyes2012",
    "gonzalez_2013": "Gonzalez2013", "chiu_2018": "Chiu2018", "kravtsov_2018": "Kravtsov2018",
    "grandis_2024": "Grandis2024",
}
GAS_FIT_KEYS = {"popesso_2024": "Popesso2024", "akino_2022": "Akino2022"}
STELLAR_FIT_KEYS = {"akino_2022": "Akino2022"}


def _ref(e):
    return f"arXiv:{e.meta.get('arxiv') or ''}|doi:{e.meta.get('doi') or ''}"


class BaryonFractions:
    """Dictionary interface (``gas_data``, ``stellar_data``, ``gas_fits``, ``stellar_fits``).

    Parameters
    ----------
    data_dir : path, optional
        Alternative data folder (BaryCon layout).
    hubble : float, optional
        Target h. Default: Planck18. Only h is changed with respect to Planck18.
    f_universal : float, optional
        Universal baryon fraction for plots and the saturating model. Default
        Planck18 Omega_b/Omega_m.
    cosmology : astropy cosmology, optional
        Full target cosmology (overrides ``hubble``).
    imf : str or None
        Convert stellar masses to this IMF (default None = as published, like the
        original package; use ``"chabrier"`` to homogenise).
    """

    FIT_FUNCTIONS = dict(_fitting.ALIASES)

    def __init__(self, data_dir=None, hubble=None, f_universal=None, cosmology=None, imf=None):
        if cosmology is None:
            cosmology = Planck18 if hubble is None else FlatLambdaCDM(
                H0=100 * float(hubble), Om0=Planck18.Om0, Ob0=Planck18.Ob0, Tcmb0=Planck18.Tcmb0)
        self.cosmology = cosmology
        self.hubble = float(cosmology.h)
        self.f_universal = float(cosmology.Ob0 / cosmology.Om0) if f_universal is None else float(f_universal)
        col = load(data_dir).to_cosmology(cosmology)
        if imf:
            col = col.to_imf(imf)
        self.collection = col
        self.data_dir = data_dir
        self.gas_data, self.stellar_data, self.gas_fits, self.stellar_fits = {}, {}, {}, {}
        for key, eid in GAS_KEYS.items():
            if eid in col:
                self.gas_data[key] = self._as_dict(col[eid], "f_gas")
        for key, eid in STELLAR_KEYS.items():
            if eid in col:
                self.stellar_data[key] = self._as_dict(col[eid], "f_star")
        for key, eid in GAS_FIT_KEYS.items():
            if eid in col:
                self.gas_fits[key] = self._fit_dict(col[eid], "f_gas")
        for key, eid in STELLAR_FIT_KEYS.items():
            if eid in col:
                self.stellar_fits[key] = self._fit_dict(col[eid], "f_star")

    # ------------------------------------------------------------- builders
    def entry(self, key):
        """The underlying BaryCon entry (metadata, history, conversions)."""
        return self.collection[GAS_KEYS.get(key) or STELLAR_KEYS.get(key) or GAS_FIT_KEYS[key]]

    def _as_dict(self, e, q):
        mkey = "M200" if e.delta == "200c" else "M500"
        if isinstance(e, Dataset):
            d = e.get(q, h_units=True)
            out = {mkey: d["M"], f"{mkey}_p": d["M_hi"], f"{mkey}_m": d["M_lo"],
                   f"{mkey}_err": 0.5 * (d["M_hi"] + d["M_lo"]),
                   q: d["f"], f"{q}_p": d["f_hi"], f"{q}_m": d["f_lo"], f"{q}_err": 0.5 * (d["f_hi"] + d["f_lo"])}
            if "z" in e.table.colnames:
                out["redshift"] = d["z"]
            if "name" in e.table.colnames:
                good = np.isfinite(np.asarray(e.table[e.meta["mass"]["column"]], float)) & \
                       np.isfinite(np.asarray(e.table[q], float))
                out["cluster_name"] = [str(n) for n in np.asarray(e.table["name"])[good]]
                if "mass_proxy" in e.table.colnames:  # compilations of surveys (Grandis+24)
                    out["survey"] = out["cluster_name"]
        else:  # a relation stored as components: points at the pivots
            comps = e.components_for(q)
            M = np.array([float(c["pivot_mass"]) for c in comps])
            f = np.array([e.evaluate(np.array([m]), q, i)[0] for i, (m, c) in enumerate(zip(M, comps))])
            err = np.array([0.5 * (c["norm_err_lo"] + c["norm_err_hi"]) for c in comps])
            out = {mkey: M * e.cosmology.h, f"{mkey}_err": None, f"{mkey}_p": None, f"{mkey}_m": None,
                   q: f, f"{q}_err": err, "survey": [c.get("label") for c in comps]}
        out["reference"] = _ref(e)
        out["status"] = e.meta["status"]
        out["derived"] = q in e.meta.get("derived_quantities", [])
        return out

    def _fit_dict(self, e, q):
        lo, hi = e.mass_range(q, 0, h_units=True)
        M = np.geomspace(lo, hi, 50)
        p16, p50, p84 = e.evaluate(M, q, 0, n_samples=2000, h_units=True)
        key = "M500_ih" if e.id.startswith("Popesso") else "M500"
        d = {key: M, q: p50, f"{q}_m": p16, f"{q}_p": p84, "reference": _ref(e)}
        if key != "M500":
            d["M500"] = M
        return d

    # -------------------------------------------------------------- fitting
    def _selection(self, data):
        if isinstance(data, str):
            data = [data]
        if isinstance(data, dict):
            raise TypeError("Pass dataset keys (dictionaries cannot carry the metadata needed to fit safely)")
        unknown = sorted(set(data) - set(self.gas_data))
        if unknown:
            raise ValueError(f"Unknown gas dataset(s): {', '.join(unknown)}. Available: {', '.join(sorted(self.gas_data))}")
        return Collection(self.collection[GAS_KEYS[k]] for k in data)

    def fit(self, data, fit_type="gas_frac_M_simple", p0=None, **kw):
        """Fit gas fractions of one or more datasets. Returns ``(popt, pcov)`` like the original.

        Extra keywords go to :func:`barycon.fitting.fit` (``scatter=True``,
        ``mass_errors=True``, ``pivot=...``). The full result is in ``self.last_fit``.
        """
        res = _fitting.fit(self._selection(data), model=fit_type, p0=p0, fb=self.f_universal, **kw)
        self.last_fit = res
        print(res.summary())
        return res.popt, res.cov

    def plot_best_fit(self, popt=None, fit_type="gas_frac_M_simple", data=None, ax=None, overplot=False,
                      **plot_kwargs):
        """Plot the last fit (or ``popt`` for ``fit_type``) with/over the selected data."""
        import matplotlib.pyplot as plt
        res = getattr(self, "last_fit", None)
        if popt is not None and (res is None or not np.allclose(popt, res.popt)):
            res = _fitting.fit(self._selection(data), model=fit_type, fb=self.f_universal)
            res.params = dict(zip(res.params, map(float, popt)))
        if ax is None:
            ax = plt.gca() if overplot else plt.subplots(figsize=(10, 8), constrained_layout=True)[1]
        if not overplot and data is not None:
            self._selection(data).plot("f_gas", ax=ax, universal=self.f_universal)
        else:
            ax.axhline(self.f_universal, color="k", alpha=0.7, label=r"$\Omega_{\rm b}/\Omega_{\rm m}$")
        res.plot(ax=ax, **plot_kwargs)
        ax.set_xscale("log")
        return ax.get_figure(), ax

    # ------------------------------------------------------------- plotting
    def _plot(self, mapping, fits_map, q, sources, fit_sources, x, figsize, save_path, overplot, **kw):
        import matplotlib.pyplot as plt
        sel = []
        avail = mapping if x == "M" else {k: v for k, v in mapping.items()}
        for k in _norm(sources, avail):
            sel.append(self.collection[mapping[k]])
        if x == "M":
            for k in _norm(fit_sources, fits_map):
                sel.append(self.collection[fits_map[k]])
        if overplot:
            ax = plt.gca()
        else:
            _, ax = plt.subplots(figsize=figsize, constrained_layout=True)
        Collection(sel).plot(q, ax=ax, x=x, universal=None if overplot or q != "f_gas" else self.f_universal, **kw)
        if save_path:
            ax.get_figure().savefig(save_path, dpi=300, bbox_inches="tight")
        return ax.get_figure(), ax

    def plot_gas_fractions_M(self, figsize=(10, 8), save_path=None, sources="all", fit_sources="all",
                             overplot=False):
        fig, ax = self._plot(GAS_KEYS, GAS_FIT_KEYS, "f_gas", sources, fit_sources, "M", figsize, save_path, overplot)
        ax.set_xlim(10 ** 12.5, 3e15)
        ax.set_ylim(0, 0.25)
        return fig, ax

    def plot_gas_fractions_z(self, figsize=(10, 8), save_path=None, sources="all", overplot=False):
        fig, ax = self._plot(GAS_KEYS, {}, "f_gas", sources, None, "z", figsize, save_path, overplot)
        ax.set_ylim(0, 0.25)
        return fig, ax

    def plot_stellar_fractions(self, figsize=(10, 8), save_path=None, sources="all", fit_sources="all"):
        return self._plot(STELLAR_KEYS, STELLAR_FIT_KEYS, "f_star", sources, fit_sources, "M", figsize,
                          save_path, False)

    # ------------------------------------------------------ legacy utilities
    @staticmethod
    def _virial_mass_converter(delta1, mass1, conc1, delta2, conc2):
        """NFW: M2 = (Delta2 c2^3)/(Delta1 c1^3) M1 (same scale radius)."""
        return delta2 * conc2**3 / (delta1 * conc1**3) * mass1

    @staticmethod
    def _virial_concentration_converter(delta1, conc1, delta2):
        def g(c2):
            return (conc1 / c2) ** 3 * (np.log(1 + c2) - c2 / (1 + c2)) / (
                np.log(1 + conc1) - conc1 / (1 + conc1)) - delta2 / delta1
        return optimize.fsolve(g, conc1)


def _norm(selection, available):
    if selection is None:
        return []
    if isinstance(selection, str):
        selection = list(available) if selection == "all" else [selection]
    unknown = set(selection) - set(available)
    if unknown:
        raise ValueError(f"Unknown key(s) {sorted(unknown)}; available: {sorted(available)}")
    return list(selection)

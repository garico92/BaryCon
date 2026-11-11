"""Core data classes: :class:`Dataset` (tabulated points) and :class:`Relation`.

Conventions
-----------
* Masses are stored **as published**, in physical ``Msun`` for the cosmology
  given in ``meta['cosmology']`` (no hidden ``h`` factors).
* Uncertainties are 1-sigma (68%) and stored as positive deltas in columns
  ``<col>_lo`` and ``<col>_hi``.
* Every transformation returns a *new* object and appends to ``.history``.
"""
from __future__ import annotations

import copy
import warnings
from pathlib import Path
from typing import Callable, Optional

import numpy as np
import yaml
from astropy import units as u
from astropy.table import Table

from . import conversions as cv
from .cosmology import cosmology_to_dict, get_cosmology
from .schema import validate_meta

__all__ = ["Entry", "Dataset", "Relation", "read_entry"]

_FRACTION_FACTOR_KEY = {"f_gas": "f_gas", "f_star": "f_star", "f_bcg": "f_star",
                        "f_sat": "f_star", "f_icl": "f_star", "f_b": None, "f_cold": None}


class Entry:
    """Common behaviour for all entries (metadata, citation, history)."""

    def __init__(self, meta: dict, path: Optional[Path] = None, validate: bool = True):
        if validate:
            validate_meta(meta)
        self.meta = meta
        self.path = path
        self.history: list[str] = []

    # ---- convenient accessors ------------------------------------------
    @property
    def id(self) -> str:
        return self.meta["id"]

    @property
    def cosmology(self):
        return get_cosmology(self.meta["cosmology"])

    @property
    def delta(self) -> str:
        return self.meta["mass"]["definition"]

    @property
    def aperture(self) -> dict:
        return self.meta["aperture"]

    @property
    def quantities(self) -> list:
        return list(self.meta["quantities"])

    @property
    def mass_proxy(self) -> str:
        return self.meta["mass"]["proxy"]

    @property
    def z_eff(self) -> float:
        r = self.meta["redshift"]
        if r.get("effective") is not None:
            return float(r["effective"])
        lo, hi = r.get("min"), r.get("max")
        if lo is None or hi is None:
            return 0.1
        return 0.5 * (float(lo) + float(hi))

    @property
    def label(self) -> str:
        return self.meta.get("label") or self.id

    def cite(self) -> str:
        """Short citation string (full BibTeX lives in ``references.bib``)."""
        m = self.meta
        s = m["reference"]
        if m.get("doi"):
            s += f", doi:{m['doi']}"
        if m.get("arxiv"):
            s += f", arXiv:{m['arxiv']}"
        return s

    def info(self) -> str:
        m = self.meta
        cz = m["cosmology"]
        lines = [
            f"{m['id']}  [{m['data_type']}, {m['kind']}, status={m['status']}]",
            f"  {m['reference']}",
            f"  quantities : {', '.join(m['quantities'])}"
            + (f"  (derived: {', '.join(m['derived_quantities'])})" if m.get("derived_quantities") else ""),
            f"  aperture   : {m['aperture']['region']} @ {m['aperture']['delta']}",
            f"  mass       : {m['mass']['definition']} from {m['mass']['proxy']}",
            f"  cosmology  : H0={cz['H0']}, Om0={cz['Om0']}",
            f"  redshift   : {m['redshift']['min']}–{m['redshift']['max']}",
        ]
        if m.get("stellar", {}).get("imf"):
            lines.append(f"  IMF        : {m['stellar']['imf']}")
        for c in m.get("caveats", []):
            lines.append(f"  ! {c}")
        if self.history:
            lines.append("  history:")
            lines += [f"    - {h}" for h in self.history]
        return "\n".join(lines)

    def __repr__(self) -> str:
        return f"<{type(self).__name__} {self.id} ({', '.join(self.quantities)} @ {self.aperture['delta']})>"

    def _copy(self):
        new = copy.deepcopy(self)
        return new

    def _set_cosmology(self, cosmo):
        d = cosmology_to_dict(get_cosmology(cosmo))
        self.meta["cosmology"] = {"H0": d["H0"], "Om0": d["Om0"], "Ode0": d["Ode0"],
                                  "note": f"converted to {d['name'] or 'user cosmology'}"}
        if d["name"]:
            self.meta["cosmology"]["name"] = d["name"]


# =============================================================================
# Tabulated points
# =============================================================================
class Dataset(Entry):
    """Tabulated measurements (one row per object or per bin)."""

    def __init__(self, table: Table, meta: dict, path=None, validate=True):
        super().__init__(meta, path, validate)
        self.table = table

    @classmethod
    def from_ecsv(cls, path, validate=True) -> "Dataset":
        t = Table.read(path, format="ascii.ecsv")
        meta = dict(t.meta)
        t.meta = {}
        # any column with a mass unit is converted to plain floats in Msun
        for c in t.colnames:
            unit = t[c].unit
            if unit is not None and unit.is_equivalent(u.Msun):
                t[c] = (np.asarray(t[c], dtype=float) * unit).to_value(u.Msun)
                t[c].unit = None
        return cls(t, meta, Path(path), validate)

    def write(self, path) -> None:
        t = self.table.copy()
        for c in t.colnames:
            if c.startswith("M_") and t[c].dtype.kind == "f":
                t[c] = np.asarray(t[c]) / 1e14
                t[c].unit = u.Unit("1e14 solMass")
        t.meta = self.meta
        t.write(path, format="ascii.ecsv", overwrite=True)

    def __len__(self):
        return len(self.table)

    # ---- data access ----------------------------------------------------
    def _col(self, name):
        if name in self.table.colnames:
            return np.asarray(self.table[name], dtype=float)
        return np.full(len(self.table), np.nan)

    def redshifts(self) -> np.ndarray:
        z = self._col("z")
        if np.all(np.isnan(z)):
            z = np.full(len(self.table), self.z_eff)
        return np.where(np.isnan(z), self.z_eff, z)

    def get(self, quantity: str, h_units: bool = False) -> dict:
        """Arrays for plotting/fitting.

        Returns a dict with ``M, M_lo, M_hi, f, f_lo, f_hi, z`` (``M`` in Msun,
        or Msun/h of the *current* cosmology if ``h_units``).
        """
        if quantity not in self.quantities:
            raise KeyError(f"{self.id} has no {quantity}; available: {self.quantities}")
        mc = self.meta["mass"]["column"]
        h = self.cosmology.h if h_units else 1.0
        out = {
            "M": self._col(mc) * h,
            "M_lo": self._col(mc + "_lo") * h,
            "M_hi": self._col(mc + "_hi") * h,
            "f": self._col(quantity),
            "f_lo": self._col(quantity + "_lo"),
            "f_hi": self._col(quantity + "_hi"),
            "z": self.redshifts(),
        }
        # optional per-row columns: local logarithmic slope of the fraction (e.g. a published
        # fit evaluated at its pivot) and a per-row mass proxy (compilations of several surveys)
        for key, col in (("slope", "slope_" + quantity), ("slope_lo", "slope_" + quantity + "_lo"),
                         ("slope_hi", "slope_" + quantity + "_hi")):
            if col in self.table.colnames:
                out[key] = self._col(col)
        for key, col in (("proxy", "mass_proxy"), ("proxy_details", "mass_proxy_details")):
            if col in self.table.colnames:
                out[key] = np.asarray(self.table[col]).astype(str)
        good = np.isfinite(out["M"]) & np.isfinite(out["f"]) & (out["M"] > 0)
        return {k: v[good] for k, v in out.items()}

    # ---- transformations ------------------------------------------------
    def _scale(self, col, factor):
        for c in (col, col + "_lo", col + "_hi"):
            if c in self.table.colnames:
                self.table[c] = np.asarray(self.table[c], dtype=float) * factor

    def to_cosmology(self, target="Planck18", **kw) -> "Dataset":
        """Rescale masses and fractions to ``target`` cosmology.

        Keyword arguments are forwarded to :func:`barycon.conversions.cosmology_factors`.
        The mass-distance exponent defaults to the value in metadata (1 if absent).
        """
        new = self._copy()
        kw.setdefault("mass_distance_exponent", self.meta["mass"].get("distance_exponent", 1.0))
        fac = cv.cosmology_factors(self.redshifts(), self.cosmology, target, **kw)
        new._scale(self.meta["mass"]["column"], fac["mass"])
        for q in self.quantities:
            key = _FRACTION_FACTOR_KEY.get(q)
            if key is None:
                continue
            new._scale(q, fac[key])
        if "f_b" in self.quantities:  # total = gas + stars, each with its own scaling
            if {"f_gas", "f_star"} <= set(self.quantities):
                fb_old = self._col("f_b")
                fb_new = new._col("f_gas") + new._col("f_star")
                new._scale("f_b", np.where(fb_old > 0, fb_new / fb_old, np.nan))
            else:
                warnings.warn(f"{self.id}: f_b not rescaled (components unknown)")
        old = self.meta["cosmology"]
        new._set_cosmology(target)
        new.history.append(
            f"cosmology H0={old['H0']},Om0={old['Om0']} -> H0={new.meta['cosmology']['H0']:.2f},"
            f"Om0={new.meta['cosmology']['Om0']:.4f}; median mass factor {np.median(fac['mass']):.4f}, "
            f"f_gas factor {np.median(fac['f_gas']):.4f}, f_star factor {np.median(fac['f_star']):.4f}"
        )
        return new

    def to_delta(
        self,
        delta_to: str,
        *,
        cm_model: str = "duffy08",
        eta_gas=None,
        f_bcg=None,
        masses_only: bool = False,
    ) -> "Dataset":
        """Convert to another overdensity definition (NFW + c–M relation).

        Fractions are model dependent: pass ``eta_gas`` (gas) and/or ``f_bcg``
        (stars, see :func:`barycon.conversions.fbcg_akino22`), or set
        ``masses_only=True`` to convert only the mass axis (the dataset is then
        flagged ``aperture_mismatch``).
        """
        if self.aperture["region"] != "sphere":
            raise ValueError(f"{self.id}: cannot convert a {self.aperture['region']} aperture")
        new = self._copy()
        mc = self.meta["mass"]["column"]
        M = self._col(mc)
        z = self.redshifts()
        Mt, rr = cv.convert_mass_definition(
            M, z, self.delta, delta_to, self.cosmology, cm_model=cm_model, return_radius_ratio=True
        )
        mfac = Mt / M
        new_col = "M_" + delta_to
        for suf in ("", "_lo", "_hi"):
            if mc + suf in new.table.colnames:
                new.table[new_col + suf] = np.asarray(new.table[mc + suf], dtype=float) * mfac
        new.meta["mass"] = dict(new.meta["mass"], column=new_col, definition=delta_to)
        notes = [f"mass {self.delta}->{delta_to} (NFW, {cm_model}); median M ratio {np.nanmedian(mfac):.3f}"]
        if not masses_only:
            kept, dropped = [], []
            for q in self.quantities:
                if q == "f_gas" and eta_gas is not None:
                    ff = cv.convert_gas_fraction_aperture(1.0, rr, eta_gas, M) * np.ones_like(M)
                    notes.append(f"f_gas aperture via r^eta (eta={eta_gas if not callable(eta_gas) else 'callable'})")
                elif q in ("f_star", "f_bcg", "f_sat", "f_icl") and f_bcg is not None:
                    ff = cv.convert_stellar_fraction_aperture(1.0, M, Mt, f_bcg) * np.ones_like(M)
                    notes.append(f"{q} aperture via BCG+satellites model")
                else:
                    dropped.append(q)
                    continue
                new._scale(q, ff)
                kept.append(q)
            if not kept:
                raise ValueError(f"{self.id}: no aperture model for {self.quantities}. Pass eta_gas "
                                 "(gas) and/or f_bcg (stars), or masses_only=True.")
            if dropped:
                notes.append(f"dropped {dropped} (no aperture model given)")
            new.meta["quantities"] = kept
            new.meta["derived_quantities"] = [q for q in self.meta.get("derived_quantities", []) if q in kept]
            new.meta["aperture"] = dict(new.meta["aperture"], delta=delta_to)
        else:
            new.meta.setdefault("caveats", []).append(
                f"aperture_mismatch: fractions still measured within r{self.delta}")
        slope_cols = [c for c in new.table.colnames if c.startswith("slope_")]
        if slope_cols and not masses_only:
            new.table.remove_columns(slope_cols)
            notes.append("local slopes dropped (not valid after the aperture change)")
        new.history.append("; ".join(notes))
        return new

    def to_imf(self, imf_to="chabrier") -> "Dataset":
        imf = (self.meta.get("stellar") or {}).get("imf")
        if imf is None:
            return self._copy()
        new = self._copy()
        fac = cv.imf_factor(imf, imf_to)
        for q in self.quantities:
            if q in ("f_star", "f_bcg", "f_sat", "f_icl"):
                new._scale(q, fac)
        new.meta["stellar"] = dict(new.meta["stellar"], imf=imf_to)
        new.history.append(f"IMF {imf}->{imf_to} (x{fac:.3f})")
        return new

    def apply_hse_bias(self, b: float, **kw) -> "Dataset":
        """Correct hydrostatic masses for a bias ``M_HSE = (1-b) M_true``."""
        if self.mass_proxy != "hydrostatic":
            warnings.warn(f"{self.id}: mass proxy is {self.mass_proxy}; bias not applied")
            return self._copy()
        new = self._copy()
        fac = cv.hse_bias_factors(b, **kw)
        new._scale(self.meta["mass"]["column"], fac["mass"])
        for q in self.quantities:
            key = _FRACTION_FACTOR_KEY.get(q)
            if key:
                new._scale(q, fac[key])
        new.history.append(f"hydrostatic bias b={b}: mass x{fac['mass']:.3f}, f_gas x{fac['f_gas']:.3f}")
        return new

    # ---- plotting -------------------------------------------------------
    def plot(self, quantity, ax=None, h_units=True, errorbars=True, **kw):
        from .plotting import plot_points
        return plot_points(self, quantity, ax=ax, h_units=h_units, errorbars=errorbars, **kw)


# =============================================================================
# Fitted relations
# =============================================================================
class Relation(Entry):
    r"""Power-law relations :math:`f = A\,(M/M_{\rm piv})^{B}` per component.

    Each component (dict) has keys ``quantity, norm, norm_err_lo, norm_err_hi,
    slope, slope_err_lo, slope_err_hi, pivot_mass, mass_range`` and optionally
    ``label``, ``corr`` (A–B correlation, default 0).
    """

    def __init__(self, components: list, meta: dict, path=None, validate=True):
        super().__init__(meta, path, validate)
        self.components = components
        # chain of transformations: each maps output mass -> (input mass, {q: multiplier})
        self._chain: list[Callable] = []

    @classmethod
    def from_yaml(cls, path, validate=True) -> "Relation":
        with open(path) as fh:
            d = yaml.safe_load(fh)
        return cls(d["components"], d["meta"], Path(path), validate)

    def write(self, path) -> None:
        with open(path, "w") as fh:
            yaml.safe_dump({"meta": self.meta, "components": self.components}, fh, sort_keys=False)

    def components_for(self, quantity):
        return [c for c in self.components if c["quantity"] == quantity]

    def _native(self, M, quantity):
        mult = np.ones_like(np.asarray(M, dtype=float))
        Min = np.asarray(M, dtype=float)
        for T in reversed(self._chain):
            Min, m = T(Min)
            mult = mult * m.get(quantity, 1.0)
        return Min, mult

    def evaluate(self, M, quantity, component=0, n_samples=0, seed=1, h_units=False):
        """Evaluate the relation on masses ``M`` (current cosmology/definition).

        With ``n_samples>0`` returns ``(p16, p50, p84)`` from Monte Carlo
        sampling of the (possibly asymmetric) parameter uncertainties.
        """
        comp = self.components_for(quantity)[component]
        M = np.asarray(M, dtype=float)
        if h_units:
            M = M / self.cosmology.h
        Mn, mult = self._native(M, quantity)
        x = Mn / float(comp["pivot_mass"])
        A, B = float(comp["norm"]), float(comp["slope"])
        if n_samples <= 0:
            return A * x**B * mult
        rng = np.random.default_rng(seed)
        lnA = _split_normal(rng, np.log(A), comp.get("norm_err_lo", 0) / A,
                            comp.get("norm_err_hi", 0) / A, n_samples)
        Bs = _split_normal(rng, B, comp.get("slope_err_lo", 0), comp.get("slope_err_hi", 0), n_samples)
        rho = float(comp.get("corr", 0.0))
        if rho:
            Bs = B + rho * (lnA - np.log(A)) * np.std(Bs) / np.std(lnA) + np.sqrt(1 - rho**2) * (Bs - B)
        samples = np.exp(lnA)[:, None] * x[None, :] ** Bs[:, None] * mult[None, :]
        return tuple(np.percentile(samples, [16, 50, 84], axis=0))

    def mass_range(self, quantity, component=0, h_units=False):
        """Validity range of the component, mapped to the current definition/cosmology."""
        comp = self.components_for(quantity)[component]
        lo, hi = (float(v) for v in comp["mass_range"])
        grid = np.geomspace(lo / 3, hi * 3, 400)
        Mn, _ = self._native(grid, quantity)
        ok = (Mn >= lo) & (Mn <= hi)
        out = grid[ok].min(), grid[ok].max()
        h = self.cosmology.h if h_units else 1.0
        return out[0] * h, out[1] * h

    # ---- transformations ------------------------------------------------
    def to_cosmology(self, target="Planck18", **kw) -> "Relation":
        new = self._copy()
        kw.setdefault("mass_distance_exponent", self.meta["mass"].get("distance_exponent", 1.0))
        fac = {k: float(np.asarray(v)[0]) for k, v in
               cv.cosmology_factors(self.z_eff, self.cosmology, target, **kw).items()}

        def T(M, fac=fac):
            return M / fac["mass"], {q: fac[k] for q, k in _FRACTION_FACTOR_KEY.items() if k}

        new._chain.append(T)
        old = self.meta["cosmology"]
        new._set_cosmology(target)
        new.history.append(
            f"cosmology H0={old['H0']} -> {new.meta['cosmology']['H0']:.2f} at z_eff={self.z_eff:.2f}: "
            f"mass x{fac['mass']:.4f}, f_gas x{fac['f_gas']:.4f}, f_star x{fac['f_star']:.4f}")
        return new

    def to_delta(self, delta_to, *, cm_model="duffy08", eta_gas=None, f_bcg=None, masses_only=False):
        new = self._copy()
        delta_from, cosmo, z = self.delta, self.cosmology, self.z_eff

        def T(Mout):
            Min, rr = cv.convert_mass_definition(Mout, z, delta_to, delta_from, cosmo,
                                                 cm_model=cm_model, return_radius_ratio=True)
            mult = {}
            if not masses_only:
                if eta_gas is not None:
                    mult["f_gas"] = cv.convert_gas_fraction_aperture(1.0, 1.0 / rr, eta_gas, Min)
                if f_bcg is not None:
                    for q in ("f_star", "f_bcg", "f_sat"):
                        mult[q] = cv.convert_stellar_fraction_aperture(1.0, Min, Mout, f_bcg)
            return Min, mult

        if not masses_only:
            need = {"f_gas": eta_gas} | {q: f_bcg for q in ("f_star", "f_bcg", "f_sat")}
            kept = [q for q in self.quantities if need.get(q) is not None]
            if not kept:
                raise ValueError(f"{self.id}: converting {self.quantities} needs eta_gas/f_bcg "
                                 "(or masses_only=True)")
            new.meta["quantities"] = kept
            new.components = [c for c in self.components if c["quantity"] in kept]
        new._chain.append(T)
        new.meta["mass"] = dict(new.meta["mass"], definition=delta_to)
        if not masses_only:
            new.meta["aperture"] = dict(new.meta["aperture"], delta=delta_to)
        new.history.append(f"definition {delta_from}->{delta_to} (NFW, {cm_model}, z_eff={z:.2f})")
        return new

    def to_imf(self, imf_to="chabrier"):
        imf = (self.meta.get("stellar") or {}).get("imf")
        new = self._copy()
        if imf is None:
            return new
        fac = cv.imf_factor(imf, imf_to)
        new._chain.append(lambda M, fac=fac: (M, {q: fac for q in ("f_star", "f_bcg", "f_sat")}))
        new.meta["stellar"] = dict(new.meta["stellar"], imf=imf_to)
        new.history.append(f"IMF {imf}->{imf_to} (x{fac:.3f})")
        return new

    def plot(self, quantity, ax=None, h_units=True, **kw):
        from .plotting import plot_relation
        return plot_relation(self, quantity, ax=ax, h_units=h_units, **kw)


def _split_normal(rng, mu, s_lo, s_hi, n):
    s_lo = s_lo if s_lo and np.isfinite(s_lo) else 0.0
    s_hi = s_hi if s_hi and np.isfinite(s_hi) else 0.0
    g = rng.standard_normal(n)
    return mu + np.where(g < 0, g * s_lo, g * s_hi)


def read_entry(path, validate=True):
    """Read an ``.ecsv`` (Dataset) or ``.yaml`` (Relation) entry."""
    path = Path(path)
    if path.suffix == ".ecsv":
        return Dataset.from_ecsv(path, validate)
    if path.suffix in (".yaml", ".yml"):
        return Relation.from_yaml(path, validate)
    raise ValueError(f"Unknown entry format: {path}")

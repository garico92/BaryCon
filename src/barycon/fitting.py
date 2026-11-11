r"""Fit analytic f(M[, z]) relations to selected entries.

Three models (from the Douspis ``BaryonFractions`` package; names in brackets
are accepted as aliases):

``power_law`` [``gas_frac_M_simple``]
    :math:`f = A\,(M/M_{\rm piv})^{\beta}` with a **fixed** pivot. The original
    version also fitted :math:`M_{\rm piv}`, which is exactly degenerate with
    :math:`A` (singular covariance).
``saturating`` [``gas_frac_M_double``]
    :math:`f = f_{\rm b}\,x^{\beta}/(1+x^{\beta})`, :math:`x=M/M_t`; tends to the
    universal fraction :math:`f_{\rm b}=\Omega_{\rm b}/\Omega_{\rm m}` at high mass.
    Free: :math:`M_t, \beta`.
``double_power_law_z`` [``gas_frac_M_z``]
    :math:`f = \tfrac{\epsilon}{2}\,(x^{\alpha}+x^{\beta})\,[E(z)/E(z_{\rm ref})]^{\gamma}`,
    :math:`x=M/M_t`, :math:`z_{\rm ref}=0.3`. Free: :math:`M_t,\alpha,\beta,\gamma,\epsilon`.
    :math:`E(z)` is taken from the cosmology the data were homogenised to.

Two estimators:

* ``scatter=False`` (default): weighted least squares (``scipy.optimize.curve_fit``),
  as in the original package;
* ``scatter=True``: Gaussian likelihood with a free *fractional* intrinsic scatter
  :math:`\sigma_{\rm int}`, :math:`\sigma^2 = \sigma_f^2 + (\sigma_{\rm int} f_{\rm model})^2`.
  This is usually needed when combining samples: their reduced :math:`\chi^2` is ≫1.

``mass_errors=True`` adds the propagated mass error,
:math:`(\partial f/\partial M)^2 \sigma_M^2`, to each point's variance (effective variance).

The fit refuses to mix different apertures or shells, and warns when the
selected entries use different mass proxies or cosmologies.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass, field

import numpy as np
from scipy import optimize

from .cosmology import get_cosmology

__all__ = ["MODELS", "ALIASES", "FitResult", "fit"]

Z_REF = 0.3


def _power_law(M, z, p, *, pivot, fb, cosmo):
    A, beta = p
    return A * (M / pivot) ** beta


def _saturating(M, z, p, *, pivot, fb, cosmo):
    Mt, beta = p
    x = (M / Mt) ** beta
    return fb * x / (1.0 + x)


def _double_power_law_z(M, z, p, *, pivot, fb, cosmo):
    Mt, alpha, beta, gamma, eps = p
    x = M / Mt
    return 0.5 * eps * (x**alpha + x**beta) * (cosmo.efunc(z) / cosmo.efunc(Z_REF)) ** gamma


MODELS = {
    "power_law": dict(func=_power_law, names=("A", "beta"), needs_z=False,
                      lower=(0.0, -np.inf), upper=(np.inf, np.inf)),
    "saturating": dict(func=_saturating, names=("M_t", "beta"), needs_z=False,
                       lower=(1e8, -np.inf), upper=(1e18, np.inf)),
    "double_power_law_z": dict(func=_double_power_law_z, names=("M_t", "alpha", "beta", "gamma", "eps"),
                               needs_z=True, lower=(1e8, -np.inf, -np.inf, -np.inf, 0.0),
                               upper=(1e18, np.inf, np.inf, np.inf, np.inf)),
}
ALIASES = {"gas_frac_M_simple": "power_law", "gas_frac_M_double": "saturating",
           "gas_frac_M_z": "double_power_law_z"}


@dataclass
class FitResult:
    """Outcome of :func:`fit`. ``params`` and ``errors`` are dicts keyed by parameter name."""

    model: str
    params: dict
    errors: dict
    cov: np.ndarray
    chi2: float
    dof: int
    n: int
    ids: list
    quantity: str
    h_units: bool
    pivot: float
    fb: float
    cosmology: object
    delta: str
    scatter: float | None = None
    notes: list = field(default_factory=list)

    @property
    def popt(self):
        return np.array(list(self.params.values()))

    def evaluate(self, M, z=None):
        """Model at masses ``M`` (same units as the fitted data) and redshift(s) ``z``."""
        M = np.asarray(M, dtype=float)
        z = np.full_like(M, Z_REF) if z is None else np.broadcast_to(np.asarray(z, float), M.shape)
        return MODELS[self.model]["func"](M, z, self.popt, pivot=self.pivot, fb=self.fb, cosmo=self.cosmology)

    def summary(self) -> str:
        s = [f"{self.model} fit of {self.quantity} (<r{self.delta}) to {self.n} points from "
             f"{len(self.ids)} entries: {', '.join(self.ids)}"]
        for k in self.params:
            s.append(f"  {k:6s} = {self.params[k]:.5g} ± {self.errors[k]:.2g}")
        if self.scatter is not None:
            s.append(f"  sigma_int (fractional) = {self.scatter:.3f}")
        s.append(f"  chi2/dof = {self.chi2:.1f}/{self.dof} = {self.chi2 / max(self.dof, 1):.2f}")
        if self.model == "power_law":
            s.append(f"  (pivot fixed at {self.pivot:.3g} {'Msun/h' if self.h_units else 'Msun'})")
        s += [f"  ! {n}" for n in self.notes]
        return "\n".join(s)

    def __repr__(self):
        return self.summary()

    def plot(self, ax=None, z=None, mass_range=None, **kw):
        """Overplot the best fit (solid over the data range, dotted when extrapolated)."""
        import matplotlib.pyplot as plt
        ax = ax or plt.gca()
        lo, hi = mass_range or self._mrange
        M = np.geomspace(lo, hi, 100)
        kw.setdefault("color", "k")
        kw.setdefault("lw", 2)
        line, = ax.plot(M, self.evaluate(M, z), label=kw.pop("label", f"{self.model} fit"), **kw)
        Mx = np.geomspace(5e12 * (0.7 if self.h_units else 1), 4e15, 100)
        ax.plot(Mx, self.evaluate(Mx, z), ls=":", color=line.get_color(), lw=1)
        return ax


def _num_hessian(fun, x, rel=1e-4):
    """Central finite-difference Hessian."""
    x = np.asarray(x, float)
    n = len(x)
    h = rel * np.maximum(np.abs(x), 1.0)
    H = np.empty((n, n))
    f0 = fun(x)
    for i in range(n):
        for j in range(i, n):
            if i == j:
                e = np.zeros(n); e[i] = h[i]
                H[i, i] = (fun(x + e) - 2 * f0 + fun(x - e)) / h[i] ** 2
            else:
                ei = np.zeros(n); ei[i] = h[i]
                ej = np.zeros(n); ej[j] = h[j]
                H[i, j] = H[j, i] = (fun(x + ei + ej) - fun(x + ei - ej) - fun(x - ei + ej)
                                     + fun(x - ei - ej)) / (4 * h[i] * h[j])
    return H


def _as_collection(data):
    from .collection import Collection
    from .dataset import Dataset
    if isinstance(data, Collection):
        return data
    if isinstance(data, Dataset):
        return Collection([data])
    return Collection(list(data))


def fit(data, model="power_law", quantity="f_gas", *, p0=None, pivot=None, h_units=True,
        scatter=False, mass_errors=False, fb=None, zmin=None, zmax=None) -> FitResult:
    """Fit ``model`` to the *points* of a Dataset or Collection.

    Homogenise first (``col.homogenize(...)``); the fit uses masses in Msun/h of
    the current cosmology when ``h_units`` (default, as in the original package).
    """
    from .dataset import Dataset
    model = ALIASES.get(model, model)
    if model not in MODELS:
        raise ValueError(f"Unknown model {model!r}; available: {list(MODELS)} (aliases {list(ALIASES)})")
    spec = MODELS[model]
    col = _as_collection(data)
    entries = [e for e in col if isinstance(e, Dataset) and quantity in e.quantities]
    if not entries:
        raise ValueError(f"No tabulated {quantity} in the selection (relations cannot be fitted)")
    deltas = {(e.aperture["delta"], e.aperture["region"]) for e in entries}
    if len(deltas) > 1:
        raise ValueError(f"Selection mixes apertures {sorted(deltas)}: convert with to_delta() or filter first")
    delta, region = deltas.pop()
    if region != "sphere":
        raise ValueError("Cannot fit shell measurements together with f(<r) models")
    notes = []
    proxies = sorted({e.mass_proxy for e in entries})
    if len(proxies) > 1:
        notes.append(f"mixes mass proxies {proxies}: consider apply_hse_bias() or fitting them separately")
    cosmos = {(round(e.meta['cosmology']['H0'], 2), round(e.meta['cosmology']['Om0'], 4)) for e in entries}
    if len(cosmos) > 1:
        warnings.warn(f"entries are in different cosmologies {sorted(cosmos)}: call homogenize() first")
        notes.append("entries NOT homogenised to a single cosmology")
    cosmo = entries[0].cosmology
    fb = float(fb) if fb is not None else float(cosmo.Ob0 / cosmo.Om0) if cosmo.Ob0 else 0.157

    M, f, s, z, sM = [], [], [], [], []
    for e in entries:
        d = e.get(quantity, h_units=h_units)
        sf = np.where(np.isfinite(d["f_lo"]) & np.isfinite(d["f_hi"]), 0.5 * (d["f_lo"] + d["f_hi"]), np.nan)
        sm = np.where(np.isfinite(d["M_lo"]) & np.isfinite(d["M_hi"]), 0.5 * (d["M_lo"] + d["M_hi"]), 0.0)
        M.append(d["M"]); f.append(d["f"]); s.append(sf); z.append(d["z"]); sM.append(sm)
    M, f, s, z, sM = map(np.concatenate, (M, f, s, z, sM))
    ok = np.isfinite(M) & np.isfinite(f) & np.isfinite(s) & (s > 0)
    if zmin is not None:
        ok &= z >= zmin
    if zmax is not None:
        ok &= z <= zmax
    if (~ok).sum():
        notes.append(f"{int((~ok).sum())} points without usable errors/redshift dropped")
    M, f, s, z, sM = M[ok], f[ok], s[ok], z[ok], sM[ok]
    pivot = float(pivot) if pivot is not None else float(10 ** np.round(np.median(np.log10(M)), 1))

    if p0 is None:
        p0 = {"power_law": [np.median(f), 0.3], "saturating": [np.median(M), 0.5],
              "double_power_law_z": [np.median(M), 0.3, 0.7, 0.0, 0.1]}[model]
    p0 = list(np.atleast_1d(p0).astype(float))
    if model == "power_law" and len(p0) == 3:  # old (Mp, beta, amp) signature
        pivot = p0[0]
        p0 = [p0[2], p0[1]]
        notes.append("3-parameter power law: pivot fixed to p0[0] (M_p and amplitude are degenerate)")

    func = spec["func"]

    def model_f(p):
        return func(M, z, p, pivot=pivot, fb=fb, cosmo=cosmo)

    def variance(p, sint=0.0):
        var = s**2 + (sint * model_f(p)) ** 2
        if mass_errors:
            eps = 1e-3
            dfdM = (func(M * (1 + eps), z, p, pivot=pivot, fb=fb, cosmo=cosmo) - model_f(p)) / (M * eps)
            var = var + (dfdM * sM) ** 2
        return var

    names = list(spec["names"])
    # internal parameters: log10 for mass-like parameters (conditioning), plus log(sigma_int)
    is_mass = np.array([n == "M_t" for n in names])

    def to_q(u):
        q = np.array(u, dtype=float)
        q[:len(names)][is_mass] = 10 ** q[:len(names)][is_mass]
        return q

    u0 = np.array(p0, dtype=float)
    u0[is_mass] = np.log10(u0[is_mass])
    lo = np.array(spec["lower"], float)
    hi = np.array(spec["upper"], float)
    lo[is_mass], hi[is_mass] = np.log10(lo[is_mass]), np.log10(hi[is_mass])

    if not scatter and not mass_errors:
        uopt, ucov = optimize.curve_fit(lambda _x, *u: model_f(to_q(np.array(u))), M, f, p0=u0, sigma=s,
                                        absolute_sigma=True, bounds=(lo, hi), maxfev=20000)
        sint = None
    else:
        def nll(u):
            p = to_q(u[:len(names)])
            sint_ = np.exp(u[-1]) if scatter else 0.0
            var = variance(p, sint_)
            return 0.5 * np.sum((f - model_f(p)) ** 2 / var + np.log(var))
        uu0 = np.concatenate([u0, [np.log(0.1)]]) if scatter else u0
        res = optimize.minimize(nll, uu0, method="Nelder-Mead",
                                options=dict(maxiter=40000, maxfev=40000, xatol=1e-9, fatol=1e-12))
        res = optimize.minimize(nll, res.x, method="Powell", options=dict(xtol=1e-10, ftol=1e-13))
        H = _num_hessian(nll, res.x)
        ucov_full = np.linalg.pinv(H)
        uopt, ucov = res.x[:len(names)], ucov_full[:len(names), :len(names)]
        sint = float(np.exp(res.x[-1])) if scatter else None
        if scatter:
            notes.append(f"sigma_int = {sint:.3f} +/- {sint * np.sqrt(max(ucov_full[-1, -1], 0)):.3f}")
    popt = to_q(uopt)
    J = np.diag(np.where(is_mass, popt * np.log(10), 1.0))  # dq/du
    pcov = J @ ucov @ J.T

    chi2 = float(np.sum((f - model_f(popt)) ** 2 / variance(popt, sint or 0.0)))
    err = np.sqrt(np.clip(np.diag(pcov), 0, np.inf))
    out = FitResult(model=model, params=dict(zip(names, map(float, popt))), errors=dict(zip(names, map(float, err))),
                    cov=pcov, chi2=chi2, dof=len(f) - len(names) - (1 if scatter else 0), n=len(f),
                    ids=[e.id for e in entries], quantity=quantity, h_units=h_units, pivot=pivot, fb=fb,
                    cosmology=get_cosmology(cosmo), delta=delta, scatter=sint, notes=notes)
    out._mrange = (M.min(), M.max())
    return out

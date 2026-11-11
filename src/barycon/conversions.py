r"""Homogenisation utilities: cosmology, overdensity definition, aperture, IMF.

Design rule: **nothing is converted silently**. Each function below returns the
converted values *and* every :class:`~barycon.dataset.Dataset` method that uses
them appends a human-readable line to ``dataset.history``.

Physics summary (see ``docs/homologation.md`` for the full derivation)
--------------------------------------------------------------------
Observations measure quantities inside a fixed *angular* aperture :math:`\theta`.
Changing cosmology at fixed redshift rescales the angular-diameter distance by
:math:`a \equiv d_A^{\rm new}(z)/d_A^{\rm old}(z)` (the luminosity distance
scales identically at fixed z), so that inside a fixed angle

* total mass (hydrostatic, and to first order lensing / scaling relations)
  :math:`M \propto d_A^{p_M}` with :math:`p_M = 1`;
* gas mass :math:`M_{\rm gas} \propto d_A^{5/2}` (emission measure);
* stellar mass :math:`M_\star \propto d_L^2 \propto d_A^2` (luminosity).

The overdensity radius is then re-solved, because :math:`\rho_c(z)\propto H^2(z)`
also changes. Writing :math:`M(<r) \propto r^{\alpha_M}` and
:math:`f(<r) \propto r^{\eta}` locally around :math:`r_\Delta`, the new aperture is
:math:`\theta_{\rm new} = s\,\theta_{\rm old}` with

.. math:: s^{3-\alpha_M} = a^{p_M-3}\,\rho_c^{\rm old}(z)/\rho_c^{\rm new}(z),

and :math:`M_\Delta \to a^{p_M} s^{\alpha_M} M_\Delta`,
:math:`f_{\rm gas} \to a^{5/2-p_M} s^{\eta_{\rm gas}} f_{\rm gas}`,
:math:`f_\star \to a^{2-p_M} s^{\eta_\star} f_\star`.

For a pure change of :math:`h` one finds :math:`s=1` exactly and recovers the
textbook :math:`M\propto h^{-1}`, :math:`f_{\rm gas}\propto h^{-3/2}`,
:math:`f_\star\propto h^{-1}`.
"""
from __future__ import annotations

import re
from typing import Callable, Union

import numpy as np
from astropy import units as u
from scipy.optimize import brentq

from .cosmology import get_cosmology

ArrayLike = Union[float, np.ndarray]
SlopeLike = Union[float, Callable[[np.ndarray], np.ndarray]]

__all__ = [
    "distance_ratio",
    "cosmology_factors",
    "hse_bias_factors",
    "parse_delta",
    "delta_threshold",
    "concentration",
    "convert_mass_definition",
    "radius_ratio",
    "convert_gas_fraction_aperture",
    "convert_stellar_fraction_aperture",
    "fbcg_akino22",
    "IMF_TO_CHABRIER",
    "imf_factor",
]

# --------------------------------------------------------------------------
# 1. Cosmology rescaling at fixed observables
# --------------------------------------------------------------------------

#: Default local log-slope of the total mass profile at r500c (NFW, c500~3).
ALPHA_M_DEFAULT = 0.9
#: Default local log-slope of the enclosed gas fraction profile near r500c.
ETA_GAS_DEFAULT = 0.3
#: Default local log-slope of the enclosed stellar fraction profile near r500c.
ETA_STAR_DEFAULT = 0.0


def distance_ratio(z: ArrayLike, cosmo_from, cosmo_to) -> np.ndarray:
    """:math:`d_A^{\\rm to}(z)/d_A^{\\rm from}(z)` (identical for :math:`d_L`)."""
    c1, c2 = get_cosmology(cosmo_from), get_cosmology(cosmo_to)
    z = np.atleast_1d(np.asarray(z, dtype=float))
    zz = np.where(z <= 0, 1e-4, z)  # ratio is finite as z -> 0 (c/H0 limit)
    return (c2.angular_diameter_distance(zz) / c1.angular_diameter_distance(zz)).decompose().value


def cosmology_factors(
    z: ArrayLike,
    cosmo_from,
    cosmo_to,
    *,
    mass_distance_exponent: float = 1.0,
    alpha_M: float = ALPHA_M_DEFAULT,
    eta_gas: float = ETA_GAS_DEFAULT,
    eta_star: float = ETA_STAR_DEFAULT,
) -> dict:
    """Multiplicative factors to move a measurement between cosmologies.

    Parameters
    ----------
    z : float or array
        Redshift(s) of the objects.
    cosmo_from, cosmo_to : cosmology specs (see :func:`barycon.get_cosmology`)
    mass_distance_exponent : float
        :math:`p_M` in :math:`M(<\\theta)\\propto d_A^{p_M}`. 1 for hydrostatic
        masses; ~1 is also used for weak-lensing and scaling-relation masses
        (an approximation, recorded in the dataset history).
    alpha_M, eta_gas, eta_star : float
        Local logarithmic slopes of :math:`M(<r)`, :math:`f_{\\rm gas}(<r)`,
        :math:`f_\\star(<r)` at the aperture. Only matter when :math:`\\Omega_m`
        or :math:`w` change; for pure-:math:`h` changes the result is exact.

    Returns
    -------
    dict with keys ``mass``, ``f_gas``, ``f_star``, ``M_gas``, ``M_star``,
    ``radius`` (physical) and ``s`` (angular aperture ratio).
    """
    c1, c2 = get_cosmology(cosmo_from), get_cosmology(cosmo_to)
    z = np.atleast_1d(np.asarray(z, dtype=float))
    a = distance_ratio(z, c1, c2)
    rho_ratio = (c1.efunc(z) * c1.H0.value) ** 2 / (c2.efunc(z) * c2.H0.value) ** 2  # rho_old/rho_new
    pM = mass_distance_exponent
    s = (a ** (pM - 3.0) * rho_ratio) ** (1.0 / (3.0 - alpha_M))
    return {
        "a": a,
        "s": s,
        "mass": a**pM * s**alpha_M,
        "radius": a * s,
        "M_gas": a**2.5 * s ** (alpha_M + eta_gas),
        "M_star": a**2.0 * s ** (alpha_M + eta_star),
        "f_gas": a ** (2.5 - pM) * s**eta_gas,
        "f_star": a ** (2.0 - pM) * s**eta_star,
    }


def hse_bias_factors(
    b: float,
    *,
    alpha_M: float = ALPHA_M_DEFAULT,
    eta_gas: float = ETA_GAS_DEFAULT,
    eta_star: float = ETA_STAR_DEFAULT,
) -> dict:
    """Factors to correct hydrostatic masses for a mass bias ``M_HSE=(1-b) M_true``.

    The true mass inside the *same* radius is larger by ``k = 1/(1-b)``; the
    overdensity radius therefore grows by ``s = k**(1/(3-alpha_M))``.
    """
    k = 1.0 / (1.0 - b)
    s = k ** (1.0 / (3.0 - alpha_M))
    return {
        "mass": k * s**alpha_M,
        "radius": s,
        "f_gas": s**eta_gas / k,
        "f_star": s**eta_star / k,
    }


# --------------------------------------------------------------------------
# 2. Overdensity definitions and NFW mass conversion
# --------------------------------------------------------------------------

_DELTA_RE = re.compile(r"^\s*(\d+(?:\.\d+)?)\s*([cm])\s*$")


def parse_delta(delta: str) -> tuple:
    """Parse ``'500c'``, ``'200m'``, ``'2500c'`` or ``'vir'``."""
    if str(delta).lower() == "vir":
        return ("vir", "c")
    m = _DELTA_RE.match(str(delta))
    if not m:
        raise ValueError(f"Cannot parse overdensity {delta!r}; use e.g. '500c', '200m', 'vir'")
    return float(m.group(1)), m.group(2)


def _rho_crit(z, cosmo) -> np.ndarray:
    return cosmo.critical_density(z).to(u.Msun / u.Mpc**3).value


def delta_threshold(delta: str, z: ArrayLike, cosmo) -> np.ndarray:
    """Density threshold :math:`\\Delta\\,\\rho_{\\rm ref}(z)` in Msun/Mpc^3."""
    cosmo = get_cosmology(cosmo)
    z = np.atleast_1d(np.asarray(z, dtype=float))
    val, ref = parse_delta(delta)
    rc = _rho_crit(z, cosmo)
    if val == "vir":  # Bryan & Norman (1998)
        x = cosmo.Om(z) - 1.0
        return (18 * np.pi**2 + 82 * x - 39 * x**2) * rc
    if ref == "c":
        return val * rc
    return val * cosmo.Om(z) * rc


def _m_nfw(x):
    return np.log1p(x) - x / (1.0 + x)


def concentration(M200c: ArrayLike, z: ArrayLike, cosmo, model: str = "duffy08") -> np.ndarray:
    """Concentration :math:`c_{200c}` for masses in Msun (physical).

    ``model='duffy08'`` uses Duffy et al. (2008), full sample, NFW, 200c.
    Any other string is passed to ``colossus.halo.concentration`` (optional
    dependency), e.g. ``'diemer19'`` or ``'ishiyama21'``.
    """
    cosmo = get_cosmology(cosmo)
    M = np.atleast_1d(np.asarray(M200c, dtype=float))
    z = np.broadcast_to(np.atleast_1d(np.asarray(z, dtype=float)), M.shape)
    if model == "duffy08":
        Mh = M * cosmo.h
        return 5.71 * (Mh / 2e12) ** (-0.084) * (1.0 + z) ** (-0.47)
    try:  # pragma: no cover - optional dependency
        from colossus.cosmology import cosmology as ccos
        from colossus.halo import concentration as cconc
    except ImportError as err:  # pragma: no cover
        raise ImportError("Install colossus to use concentration model %r" % model) from err
    ccos.setCosmology(
        "barycon_tmp",
        params=dict(flat=True, H0=cosmo.H0.value, Om0=cosmo.Om0, Ob0=0.049, sigma8=0.81, ns=0.96),
    )
    return np.array([cconc.concentration(m * cosmo.h, "200c", zz, model=model) for m, zz in zip(M, z)])


def _nfw_convert_one(M200c, c200, z, delta_to, cosmo):
    thr200 = delta_threshold("200c", z, cosmo)[0]
    r200 = (3 * M200c / (4 * np.pi * thr200)) ** (1 / 3)
    rs = r200 / c200
    thr = delta_threshold(delta_to, z, cosmo)[0]

    def g(x):
        return M200c * _m_nfw(x) / _m_nfw(c200) - thr * 4 * np.pi / 3 * (x * rs) ** 3

    x = brentq(g, 1e-3, 1e3)
    return M200c * _m_nfw(x) / _m_nfw(c200), x * rs


def convert_mass_definition(
    M: ArrayLike,
    z: ArrayLike,
    delta_from: str,
    delta_to: str,
    cosmo,
    *,
    cm_model: str = "duffy08",
    c: ArrayLike | None = None,
    return_radius_ratio: bool = False,
):
    """Convert halo masses between overdensity definitions assuming NFW.

    Parameters
    ----------
    M : masses in **physical Msun** for the given ``cosmo``.
    z : redshift(s).
    delta_from, delta_to : e.g. ``'500c'``, ``'200c'``, ``'200m'``, ``'vir'``.
    cosmo : cosmology used to compute the reference densities.
    cm_model : concentration-mass relation (see :func:`concentration`).
    c : optional concentration :math:`c_{200c}` per object (overrides cm_model),
        useful when the paper provides measured concentrations (e.g. Ettori+10).
    return_radius_ratio : also return :math:`r_{\\rm to}/r_{\\rm from}`.
    """
    cosmo = get_cosmology(cosmo)
    M = np.atleast_1d(np.asarray(M, dtype=float))
    z = np.broadcast_to(np.atleast_1d(np.asarray(z, dtype=float)), M.shape).copy()
    out = np.full_like(M, np.nan)
    rr = np.full_like(M, np.nan)
    cc = None if c is None else np.broadcast_to(np.atleast_1d(np.asarray(c, dtype=float)), M.shape)
    for i, (m, zz) in enumerate(zip(M, z)):
        if not np.isfinite(m) or m <= 0:
            continue

        # find M200c (and c200c) such that the 'from' mass matches m
        def mismatch(lnM200):
            M200 = np.exp(lnM200)
            ci = cc[i] if cc is not None else concentration(M200, zz, cosmo, cm_model)[0]
            mf, _ = _nfw_convert_one(M200, ci, zz, delta_from, cosmo)
            return np.log(mf) - np.log(m)

        lnM200 = brentq(mismatch, np.log(m) - 3, np.log(m) + 3)
        M200 = np.exp(lnM200)
        ci = cc[i] if cc is not None else concentration(M200, zz, cosmo, cm_model)[0]
        mt, rt = _nfw_convert_one(M200, ci, zz, delta_to, cosmo)
        _, rf = _nfw_convert_one(M200, ci, zz, delta_from, cosmo)
        out[i], rr[i] = mt, rt / rf
    return (out, rr) if return_radius_ratio else out


def radius_ratio(M_from, M_to, z, delta_from, delta_to, cosmo) -> np.ndarray:
    """:math:`r_{\\rm to}/r_{\\rm from}` given the masses in the two definitions."""
    thr_f = delta_threshold(delta_from, z, cosmo)
    thr_t = delta_threshold(delta_to, z, cosmo)
    return (np.asarray(M_to) / thr_t / (np.asarray(M_from) / thr_f)) ** (1 / 3)


# --------------------------------------------------------------------------
# 3. Fractions in a different aperture (model-dependent: explicit choices!)
# --------------------------------------------------------------------------


def _eval_slope(eta: SlopeLike, M):
    return eta(M) if callable(eta) else eta


def convert_gas_fraction_aperture(f, r_ratio, eta: SlopeLike, M=None):
    """Move :math:`f_{\\rm gas}(<r)` to a different radius assuming a local
    power law :math:`f_{\\rm gas}(<r)\\propto r^{\\eta}`.

    ``eta`` has **no default**: it is strongly mass dependent (≈0.2–0.4 for
    massive clusters between r500 and r200, X-COP; ≳0.7 for groups, see e.g.
    Popesso et al. 2024 whose r500 and r200 fits imply η≈0.7–1). Pass a float
    or a callable of the (from-aperture) mass.
    """
    return np.asarray(f) * np.asarray(r_ratio) ** _eval_slope(eta, M)


def fbcg_akino22(M500c_msun) -> np.ndarray:
    """Fraction of the stellar mass within r500c that sits in the BCG.

    Ratio of the Akino et al. (2022, Appendix E) BCG and total stellar mass
    regressions (h=0.7, M500c in Msun):
    ``exp(-1.07-0.68) * (M/1e14)**(0.39-0.76)``. Clipped to [0, 0.99].
    """
    M = np.asarray(M500c_msun, dtype=float)
    return np.clip(np.exp(-1.07 - 0.68) * (M / 1e14) ** (0.39 - 0.76), 0.0, 0.99)


def convert_stellar_fraction_aperture(f, M_from, M_to, f_bcg: SlopeLike):
    """Move :math:`f_\\star` between apertures with a BCG + satellites model.

    The BCG (incl. its ICL, if counted) is fully enclosed in both apertures;
    satellites trace the total mass profile. ``f_bcg`` is the fraction of the
    stellar mass inside the *from* aperture contributed by the BCG (float or
    callable of ``M_from``), e.g. :func:`fbcg_akino22`.
    """
    fb = _eval_slope(f_bcg, M_from)
    M_from, M_to = np.asarray(M_from), np.asarray(M_to)
    return np.asarray(f) * (fb * M_from / M_to + (1.0 - fb))


# --------------------------------------------------------------------------
# 4. IMF
# --------------------------------------------------------------------------

#: Multiplicative factors to convert stellar masses to a Chabrier (2003) IMF.
#: From Madau & Dickinson (2014): M(Chabrier)=0.63 M(Salpeter),
#: M(Kroupa)=0.67 M(Salpeter); 'diet Salpeter' = 0.7 Salpeter (Bell+03).
IMF_TO_CHABRIER = {
    "chabrier": 1.0,
    "salpeter": 0.63,
    "kroupa": 0.63 / 0.67,
    "diet_salpeter": 0.63 / 0.70,
}


def imf_factor(imf_from: str, imf_to: str = "chabrier") -> float:
    """Factor to multiply stellar masses/fractions computed with ``imf_from``."""
    f = IMF_TO_CHABRIER
    a, b = imf_from.lower().replace(" ", "_"), imf_to.lower().replace(" ", "_")
    if a not in f or b not in f:
        raise KeyError(f"Unknown IMF; known: {sorted(f)}")
    return f[a] / f[b]

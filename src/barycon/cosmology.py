"""Cosmology helpers.

Every dataset in BaryCon stores the cosmology *assumed by the original authors*
(``meta['cosmology']``). Conversions to a common target cosmology are done
explicitly by the user (see :mod:`barycon.conversions`), never silently.
"""
from __future__ import annotations

from typing import Any, Mapping

import astropy.cosmology as _ac
from astropy.cosmology import FLRW, FlatLambdaCDM, LambdaCDM

__all__ = ["get_cosmology", "cosmology_to_dict", "DEFAULT_TARGET"]

#: Default target cosmology for homogenisation (astropy's Planck18).
DEFAULT_TARGET = "Planck18"


def get_cosmology(spec: Any) -> FLRW:
    """Return an astropy cosmology from a flexible specification.

    Parameters
    ----------
    spec : astropy FLRW, str or mapping
        * an astropy cosmology instance (returned unchanged);
        * the name of an astropy built-in realisation, e.g. ``"Planck18"``,
          ``"WMAP9"``;
        * a mapping with at least ``H0`` and ``Om0`` (``Ode0`` optional;
          if missing a flat model is assumed). Extra keys are ignored.

    Examples
    --------
    >>> get_cosmology({"H0": 70, "Om0": 0.3}).h
    0.7
    """
    if isinstance(spec, FLRW):
        return spec
    if isinstance(spec, str):
        try:
            return getattr(_ac, spec)
        except AttributeError as err:  # pragma: no cover - message only
            raise ValueError(
                f"Unknown cosmology name {spec!r}. Use an astropy realisation "
                "(e.g. 'Planck18') or a dict with H0, Om0[, Ode0]."
            ) from err
    if isinstance(spec, Mapping):
        name = spec.get("name")
        if name and isinstance(getattr(_ac, str(name), None), FLRW):
            return getattr(_ac, name)  # exact realisation (keeps neutrinos/radiation)
        H0 = float(spec["H0"])
        Om0 = float(spec["Om0"])
        Ode0 = spec.get("Ode0")
        if Ode0 is None or abs(float(Ode0) + Om0 - 1.0) < 1e-6:
            return FlatLambdaCDM(H0=H0, Om0=Om0)
        return LambdaCDM(H0=H0, Om0=Om0, Ode0=float(Ode0))
    raise TypeError(f"Cannot interpret {spec!r} as a cosmology")


def cosmology_to_dict(cosmo: FLRW) -> dict:
    """Minimal serialisable description of an astropy cosmology."""
    return {
        "name": getattr(cosmo, "name", None),
        "H0": float(cosmo.H0.value),
        "Om0": float(cosmo.Om0),
        "Ode0": float(cosmo.Ode0),
    }

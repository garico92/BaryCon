"""Backward-compatible import name of the original Douspis package.

``from baryon_fractions import BaryonFractions`` keeps working; new code should use
``import barycon``.
"""
import warnings as _w

from barycon.compat import BaryonFractions  # noqa: F401

_w.warn("`baryon_fractions` is now `barycon`; please `from barycon import BaryonFractions`",
        DeprecationWarning, stacklevel=2)

"""BaryCon — a community collection of baryon-fraction constraints in haloes.

Quick start
-----------
>>> import barycon as bc
>>> col = bc.load()                                     # every bundled entry
>>> gas = col.filter(quantity="f_gas", delta="500c")    # gas fractions in r500c
>>> gas_p18 = gas.homogenize(cosmology="Planck18")      # explicit, logged conversion
>>> gas_p18["Vikhlinin2006"].history
>>> ax = gas_p18.plot("f_gas", universal=0.157)

Inspired by CoReCon (Garaldi 2023, JOSS 8, 5407).
"""
from ._version import __version__, __data_version__
from .collection import Collection, data_dir, load
from .conversions import (
    convert_mass_definition,
    cosmology_factors,
    fbcg_akino22,
    imf_factor,
)
from .cosmology import get_cosmology
from .dataset import Dataset, Relation, read_entry
from .schema import validate_meta
from .fitting import fit, FitResult
from .compat import BaryonFractions

__all__ = [
    "__version__", "__data_version__", "load", "Collection", "Dataset", "Relation",
    "read_entry", "data_dir", "get_cosmology", "cosmology_factors",
    "convert_mass_definition", "fbcg_akino22", "imf_factor", "validate_meta",
    "fit", "FitResult", "BaryonFractions",
]

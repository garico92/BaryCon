"""Metadata schema for BaryCon entries and a validator.

Every data file (``.ecsv`` for tabulated points, ``.yaml`` for fitted
relations) carries the same metadata block. The schema is deliberately
explicit about everything needed to *homogenise* measurements: assumed
cosmology, overdensity/aperture, how the total mass was obtained, IMF, etc.
"""
from __future__ import annotations

import jsonschema

QUANTITIES = ["f_gas", "f_star", "f_bcg", "f_sat", "f_icl", "f_b", "f_cold"]
MASS_PROXIES = [
    "hydrostatic",
    "weak_lensing",
    "xray_scaling_relation",
    "sz_scaling_relation",
    "optical_scaling_relation",
    "dynamical",
    "abundance_matching",
    "simulation",
    "mixed",
]

SCHEMA = {
    "$schema": "https://json-schema.org/draft/2020-12/schema",
    "title": "BaryCon entry metadata",
    "type": "object",
    "required": [
        "id", "reference", "url", "description", "kind", "data_type",
        "quantities", "aperture", "mass", "cosmology", "redshift", "extracted",
        "status",
    ],
    "properties": {
        "id": {"type": "string", "pattern": r"^[A-Za-z]+[0-9]{4}[A-Za-z0-9_]*$"},
        "reference": {"type": "string", "minLength": 5},
        "bibcode": {"type": ["string", "null"]},
        "doi": {"type": ["string", "null"]},
        "arxiv": {"type": ["string", "null"]},
        "url": {"type": "string"},
        "description": {"type": "string", "minLength": 10},
        "kind": {"enum": ["observation", "simulation", "emulator"]},
        "data_type": {"enum": ["points", "binned", "relation"]},
        "quantities": {"type": "array", "items": {"enum": QUANTITIES}, "minItems": 1},
        "derived_quantities": {"type": "array", "items": {"enum": QUANTITIES}},
        "aperture": {
            "type": "object",
            "required": ["delta", "region"],
            "properties": {
                "delta": {"type": "string", "pattern": r"^(\d+[cm]|vir)$"},
                "region": {"enum": ["sphere", "shell", "cylinder"]},
                "r_inner": {"type": ["number", "null"]},
                "r_outer": {"type": ["number", "null"]},
            },
        },
        "mass": {
            "type": "object",
            "required": ["column", "definition", "proxy", "unit"],
            "properties": {
                "column": {"type": "string"},
                "definition": {"type": "string", "pattern": r"^(\d+[cm]|vir)$"},
                "proxy": {"enum": MASS_PROXIES},
                "proxy_details": {"type": "string"},
                "unit": {"enum": ["Msun"]},
                "distance_exponent": {"type": "number"},
            },
        },
        "cosmology": {
            "type": "object",
            "required": ["H0", "Om0"],
            "properties": {
                "H0": {"type": "number", "minimum": 40, "maximum": 100},
                "Om0": {"type": "number", "minimum": 0.05, "maximum": 1.0},
                "Ode0": {"type": ["number", "null"]},
                "note": {"type": "string"},
            },
        },
        "stellar": {
            "type": "object",
            "properties": {
                "imf": {"enum": ["chabrier", "salpeter", "kroupa", "diet_salpeter", None]},
                "includes_bcg": {"type": ["boolean", "null"]},
                "includes_icl": {"type": ["boolean", "null"]},
            },
        },
        "redshift": {
            "type": "object",
            "required": ["min", "max"],
            "properties": {
                "min": {"type": ["number", "null"]},
                "max": {"type": ["number", "null"]},
                "effective": {"type": ["number", "null"]},
            },
        },
        "selection": {"type": "string"},
        "extracted": {"type": "boolean"},
        "table_source": {"type": "string"},
        "caveats": {"type": "array", "items": {"type": "string"}},
        "status": {"enum": ["verified", "unverified", "needs_check", "draft"]},
        "verified_by": {"type": ["string", "null"]},
        "contributors": {"type": "array", "items": {"type": "string"}},
        "added": {"type": "string"},
        "changelog": {"type": "array", "items": {"type": "string"}},
        "overlaps_with": {"type": "array", "items": {"type": "string"}},
    },
}


def validate_meta(meta: dict) -> None:
    """Raise :class:`jsonschema.ValidationError` if ``meta`` is invalid."""
    jsonschema.validate(meta, SCHEMA)

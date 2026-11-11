"""A searchable, homogenisable collection of BaryCon entries."""
from __future__ import annotations

import warnings
from importlib import resources
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
from astropy.table import Table, vstack

from .dataset import Dataset, Relation, read_entry

__all__ = ["Collection", "load", "data_dir"]


def data_dir() -> Path:
    """Directory holding the bundled data files."""
    return Path(resources.files("barycon") / "data")


class Collection:
    """Ordered mapping ``id -> entry`` with filtering and bulk conversions."""

    def __init__(self, entries: Iterable = ()):
        self._entries = {e.id: e for e in entries}

    # ---- mapping protocol ----------------------------------------------
    def __getitem__(self, key):
        return self._entries[key]

    def __iter__(self):
        return iter(self._entries.values())

    def __len__(self):
        return len(self._entries)

    def __contains__(self, key):
        return key in self._entries

    @property
    def ids(self):
        return list(self._entries)

    def __repr__(self):
        return f"<Collection of {len(self)} entries: {', '.join(self.ids)}>"

    # ---- discovery -------------------------------------------------------
    def summary(self) -> Table:
        """One row per entry with the information needed to homogenise."""
        rows = []
        for e in self:
            m = e.meta
            rows.append({
                "id": e.id,
                "type": m["data_type"],
                "quantities": ",".join(m["quantities"]),
                "aperture": f"{m['aperture']['region']}<{m['aperture']['delta']}",
                "mass_proxy": m["mass"]["proxy"],
                "H0": m["cosmology"]["H0"],
                "Om0": m["cosmology"]["Om0"],
                "zmin": m["redshift"]["min"] if m["redshift"]["min"] is not None else np.nan,
                "zmax": m["redshift"]["max"] if m["redshift"]["max"] is not None else np.nan,
                "N": len(e) if isinstance(e, Dataset) else 0,
                "imf": (m.get("stellar") or {}).get("imf") or "",
                "status": m["status"],
            })
        return Table(rows=rows)

    def filter(
        self,
        quantity: Optional[str] = None,
        delta: Optional[str] = None,
        region: Optional[str] = None,
        mass_proxy=None,
        kind=None,
        data_type=None,
        zmin: Optional[float] = None,
        zmax: Optional[float] = None,
        status=None,
        exclude_derived: bool = True,
        ids=None,
    ) -> "Collection":
        """Return a sub-collection.

        ``exclude_derived`` drops entries where ``quantity`` is only *derived*
        from scaling relations (e.g. Giodini+09 gas fractions).
        """
        def _as_set(x):
            return None if x is None else ({x} if isinstance(x, str) else set(x))

        mp, kd, dt, st, ii = map(_as_set, (mass_proxy, kind, data_type, status, ids))
        out = []
        for e in self:
            m = e.meta
            if quantity and quantity not in m["quantities"]:
                continue
            if quantity and exclude_derived and quantity in m.get("derived_quantities", []):
                continue
            if delta and m["aperture"]["delta"] != delta:
                continue
            if region and m["aperture"]["region"] != region:
                continue
            if mp and m["mass"]["proxy"] not in mp:
                continue
            if kd and m["kind"] not in kd:
                continue
            if dt and m["data_type"] not in dt:
                continue
            if st and m["status"] not in st:
                continue
            if ii and e.id not in ii:
                continue
            zlo, zhi = m["redshift"]["min"], m["redshift"]["max"]
            if zmax is not None and zlo is not None and zlo > zmax:
                continue
            if zmin is not None and zhi is not None and zhi < zmin:
                continue
            out.append(e)
        return Collection(out)

    # ---- bulk transformations -------------------------------------------
    def _apply(self, method, *args, skip_errors=True, **kw):
        out = []
        for e in self:
            try:
                out.append(getattr(e, method)(*args, **kw))
            except Exception as err:  # keep going, but tell the user
                if not skip_errors:
                    raise
                warnings.warn(f"{e.id}: {method} skipped ({err})")
        return Collection(out)

    def to_cosmology(self, target="Planck18", **kw) -> "Collection":
        return self._apply("to_cosmology", target, **kw)

    def to_delta(self, delta_to, **kw) -> "Collection":
        return self._apply("to_delta", delta_to, **kw)

    def to_imf(self, imf_to="chabrier") -> "Collection":
        return self._apply("to_imf", imf_to)

    def homogenize(self, cosmology="Planck18", delta=None, imf="chabrier", **kw) -> "Collection":
        """Convenience: cosmology -> IMF -> (optionally) overdensity definition."""
        c = self.to_cosmology(cosmology).to_imf(imf)
        if delta is not None:
            c = c.to_delta(delta, **kw)
        return c

    # ---- export -----------------------------------------------------------
    def to_table(self, quantity: str, h_units: bool = False) -> Table:
        """Long-format table of all *points* for ``quantity`` (relations skipped)."""
        tabs = []
        for e in self:
            if not isinstance(e, Dataset) or quantity not in e.quantities:
                continue
            d = e.get(quantity, h_units=h_units)
            n = len(d["M"])
            if n == 0:
                continue
            t = Table(d)
            t["id"] = [e.id] * n
            t["delta"] = [e.aperture["delta"]] * n
            t["mass_proxy"] = [e.mass_proxy] * n
            tabs.append(t)
        return vstack(tabs) if tabs else Table()

    def bibtex(self) -> str:
        """BibTeX entries for every entry in the collection."""
        bib = (data_dir() / "references.bib").read_text()
        blocks = {}
        for blk in bib.split("\n@")[0:]:
            blk = blk if blk.startswith("@") else "@" + blk
            key = blk.split("{", 1)[1].split(",", 1)[0].strip() if "{" in blk else None
            if key:
                blocks[key] = blk.strip()
        keys = list(dict.fromkeys(e.id.split("_")[0] for e in self))  # Foo2013_2500c -> Foo2013
        missing = [k for k in keys if k not in blocks]
        if missing:
            warnings.warn(f"no BibTeX for {missing}")
        return "\n\n".join(blocks[k] for k in keys if k in blocks)

    def plot(self, quantity, ax=None, h_units=True, x="M", **kw):
        """Plot ``quantity`` against mass (``x="M"``) or redshift (``x="z"``)."""
        from .plotting import plot_collection
        return plot_collection(self, quantity, ax=ax, h_units=h_units, x=x, **kw)

    def fit(self, model="power_law", quantity="f_gas", **kw):
        """Fit an analytic relation to the points; see :func:`barycon.fitting.fit`."""
        from .fitting import fit
        return fit(self, model=model, quantity=quantity, **kw)


def load(path=None, include_pending=False, validate=True) -> Collection:
    """Load all entries shipped with the package (or from ``path``)."""
    base = Path(path) if path else data_dir()
    dirs = [base / "observations"] + ([base / "pending"] if include_pending else [])
    files = sorted(p for d in dirs if d.exists() for p in d.iterdir()
                   if p.suffix in (".ecsv", ".yaml", ".yml"))
    return Collection(read_entry(f, validate=validate) for f in files)

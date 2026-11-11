"""Consistent plotting helpers (matplotlib).

Colour identifies the *publication* (the same colour in every figure, so a dataset
can be followed across gas, stellar and redshift plots); the marker shape encodes
*how the total mass was obtained* (hydrostatic, weak lensing, scaling relations,
dynamics). Published fits are drawn as a median line with a 68% band. Points that
carry a local logarithmic slope (e.g. Grandis+24, a fit evaluated at its pivot) get
a short segment showing that slope around the point.
"""
from __future__ import annotations

import re

import numpy as np

#: one colour per publication; fixed, so a dataset keeps its colour in every figure
DATASET_COLORS = {
    "Vikhlinin2006": "#2a78d6",
    "Arnaud2007": "#eb6834",
    "Sun2009": "#1baf7a",
    "Ettori2010": "#6250d6",
    "Sanderson2013": "#e34948",
    "Landry2013": "#eda100",
    "Zhang2011": "#008300",
    "Wicker2023": "#a58d7f",
    "Mulroy2019": "#e87ba4",
    "Chiu2018": "#0b3d91",
    "Gonzalez2013": "#8c564b",
    "Giodini2009": "#17a2b8",
    "Khalil2026": "#b8860b",
    "Bulbul2024": "#7f7f7f",
    "Grandis2024": "#1a1a1a",
    "Akino2022": "#c2185b",
    "Popesso2024": "#7cb342",
    "Kravtsov2018": "#00897b",
    "Reyes2012": "#ffa000",
    "Mantz2022": "#455a64",
}
_FALLBACK = ["#5c6bc0", "#8d6e63", "#26a69a", "#ab47bc", "#d4e157", "#ef5350", "#29b6f6", "#9ccc65"]

#: marker and legend label per mass-proxy family
PROXY_MARKERS = {
    "hydrostatic": ("o", "X-ray hydrostatic"),
    "wl": ("s", "weak lensing / WL-calibrated"),
    "xray": ("D", "X-ray scaling relation"),
    "sz": ("^", "SZ scaling relation"),
    "dynamical": ("v", "galaxy dynamics"),
    "other": ("P", "other / mixed"),
}

_YLABEL = {
    "f_gas": r"$f_{\rm gas}$",
    "f_star": r"$f_\star$",
    "f_bcg": r"$f_{\rm BCG}$",
    "f_b": r"$f_{\rm b}$",
}


def publication(entry) -> str:
    """Publication key of an entry: its id without aperture suffixes (``Landry2013_2500c`` -> ``Landry2013``)."""
    return entry.id.split("_")[0]


def short_label(entry) -> str:
    """Compact legend label, e.g. ``Vikhlinin2006`` -> ``Vikhlinin+06``."""
    m = re.match(r"([A-Za-z]+)(\d{4})", publication(entry))
    return f"{m.group(1)}+{m.group(2)[2:]}" if m else entry.id


def dataset_color(entry) -> str:
    key = publication(entry)
    if key in DATASET_COLORS:
        return DATASET_COLORS[key]
    return _FALLBACK[sum(map(ord, key)) % len(_FALLBACK)]


def proxy_family(proxy: str, details: str = "") -> str:
    """Map a schema mass proxy (plus its free-text details) onto a marker family."""
    det = (details or "").lower()
    if proxy == "hydrostatic":
        return "hydrostatic"
    if proxy == "weak_lensing" or re.search(r"\bwl\b|weak.lensing|lensing", det):
        return "wl"
    if proxy == "xray_scaling_relation":
        return "xray"
    if proxy == "sz_scaling_relation":
        return "sz"
    if proxy == "dynamical":
        return "dynamical"
    return "other"


def entry_family(entry) -> str:
    m = entry.meta["mass"]
    return proxy_family(m["proxy"], m.get("proxy_details", ""))


def _ax(ax):
    if ax is None:
        import matplotlib.pyplot as plt
        _, ax = plt.subplots(figsize=(6, 4.5), constrained_layout=True)
    return ax


def plot_points(e, quantity, ax=None, h_units=True, errorbars=True, max_points_errorbars=500, x="M",
                slope_dex=0.2, **kw):
    """Scatter one dataset; returns the axes. Marker per row follows the mass proxy."""
    ax = _ax(ax)
    d = e.get(quantity, h_units=h_units)
    if x == "z":  # fraction against redshift
        if "z" not in e.table.colnames:
            return ax
        d = dict(d, M=d["z"], M_lo=np.full_like(d["z"], np.nan), M_hi=np.full_like(d["z"], np.nan))
        d.pop("slope", None)
    color = kw.pop("color", dataset_color(e))
    kw.pop("label", None)
    kw.pop("marker", None)
    n = len(d["M"])
    if "proxy" in d:
        fam = np.array([proxy_family(p, q) for p, q in zip(d["proxy"], d.get("proxy_details", [""] * n))])
    else:
        fam = np.full(n, entry_family(e))
    ms = kw.pop("ms", 5)
    if n > max_points_errorbars:  # dense catalogues: faint cloud, no error bars
        for f in np.unique(fam):
            s = fam == f
            ax.plot(d["M"][s], d["f"][s], ls="", marker=PROXY_MARKERS[f][0], ms=2.2, mew=0, color=color,
                    alpha=0.12, zorder=1, rasterized=True)
        return ax
    for f in np.unique(fam):
        s = fam == f
        xerr = np.vstack([d["M_lo"][s], d["M_hi"][s]]) if errorbars else None
        yerr = np.vstack([d["f_lo"][s], d["f_hi"][s]]) if errorbars else None
        if xerr is not None and not np.any(np.isfinite(xerr)):
            xerr = None
        if yerr is not None and not np.any(np.isfinite(yerr)):
            yerr = None
        has_slope = "slope" in d and np.any(np.isfinite(d["slope"]))
        ax.errorbar(d["M"][s], d["f"][s], xerr=xerr, yerr=yerr, ls="", marker=PROXY_MARKERS[f][0],
                    ms=ms + 1.5 * has_slope, color=color, mec="white", mew=0.5, elinewidth=0.7, alpha=0.9,
                    zorder=5 if has_slope else 3, **kw)
    if x == "M" and "slope" in d and slope_dex:
        for M0, f0, b in zip(d["M"], d["f"], d["slope"]):
            if np.isfinite(b):
                Ms = M0 * 10 ** np.array([-slope_dex, slope_dex])
                ax.plot(Ms, f0 * (Ms / M0) ** b, color=color, lw=1.6, zorder=4, solid_capstyle="round")
    return ax


def plot_relation(e, quantity, ax=None, h_units=True, n_samples=2000, **kw):
    """Median line and 68% band of a published fit."""
    ax = _ax(ax)
    color = kw.pop("color", dataset_color(e))
    ls_cycle = ["-", "--", ":", "-."]
    comps = e.components_for(quantity)
    for i, comp in enumerate(comps):
        lo, hi = e.mass_range(quantity, i, h_units=h_units)
        M = np.geomspace(lo, hi, 60)
        p16, p50, p84 = e.evaluate(M, quantity, i, n_samples=n_samples, h_units=h_units)
        ax.plot(M, p50, color=color, lw=1.6, ls=kw.get("ls", ls_cycle[i % 4]), zorder=2)
        ax.fill_between(M, p16, p84, color=color, alpha=0.15, lw=0, zorder=0)
    return ax


def plot_collection(col, quantity, ax=None, h_units=True, universal=None, legend=True,
                    proxy_legend_loc="upper left", legend_loc="below", legend_ncol=4, x="M", **kw):
    """Plot every entry of a collection that has ``quantity``.

    Two legends: the datasets (colour) and the mass proxies (marker). With
    ``legend_loc="below"`` the dataset legend is placed under the axes, across the
    figure width.
    """
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    from .dataset import Dataset

    ax = _ax(ax)
    handles, families = {}, set()
    for e in col:
        if quantity not in e.quantities:
            continue
        if isinstance(e, Dataset):
            if x == "z" and "z" not in e.table.colnames:
                continue
            plot_points(e, quantity, ax=ax, h_units=h_units, x=x, **kw)
            d = e.get(quantity)
            if "proxy" in d:
                families |= {proxy_family(p, q) for p, q in zip(d["proxy"], d["proxy_details"])}
            else:
                families.add(entry_family(e))
            handles.setdefault(short_label(e), Patch(color=dataset_color(e), label=short_label(e)))
        elif x == "M":
            plot_relation(e, quantity, ax=ax, h_units=h_units, **kw)
            handles.setdefault(short_label(e) + " (fit)",
                               (Patch(color=dataset_color(e), alpha=0.25, lw=0),
                                Line2D([], [], color=dataset_color(e), lw=1.6)))
    if universal is not None:
        ax.axhline(universal, color="#52514e", lw=1, ls="--", zorder=1)
        ax.text(0.99, universal, r"$\Omega_{\rm b}/\Omega_{\rm m}$", transform=ax.get_yaxis_transform(),
                ha="right", va="bottom", color="#52514e")
    deltas = sorted({e.aperture["delta"] for e in col if quantity in e.quantities})
    dl = "/".join(deltas)
    if x == "z":
        ax.set_xlabel("redshift $z$")
    else:
        ax.set_xscale("log")
        unit = r"h^{-1}\,{\rm M}_\odot" if h_units else r"{\rm M}_\odot"
        ax.set_xlabel(rf"$M_{{\rm {dl}}}\ [{unit}]$")
    ax.set_ylabel(_YLABEL.get(quantity, quantity) + rf" $(<r_{{\rm {dl}}})$")
    ax.grid(alpha=0.25, lw=0.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    if legend:
        fam = [Line2D([], [], ls="", marker=PROXY_MARKERS[k][0], color="#52514e", mec="white", ms=6,
                      label=PROXY_MARKERS[k][1]) for k in PROXY_MARKERS if k in families]
        leg1 = ax.legend(handles=fam, loc=proxy_legend_loc, frameon=False, handletextpad=0.2,
                         borderaxespad=0.2, labelspacing=0.25)
        ax.add_artist(leg1)
        from matplotlib.legend_handler import HandlerTuple
        labels = sorted(handles)
        hs = [handles[k] for k in labels]
        common = dict(handles=hs, labels=labels, frameon=False, ncol=legend_ncol, columnspacing=1.0,
                      handlelength=1.2, handletextpad=0.4, labelspacing=0.3,
                      handler_map={tuple: HandlerTuple(ndivide=1, pad=0)})
        if legend_loc == "below":
            ax.figure.legend(loc="outside lower center", **common)
        else:
            ax.legend(loc=legend_loc, **common)
    return ax

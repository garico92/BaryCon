"""Reproduce the figures of the README.

    python examples/make_figures.py [LEGACY_DATA_DIR] [OUTDIR]

LEGACY_DATA_DIR (optional) points to the old ``baryon_fractions/data`` folder and
enables the "legacy vs fixed" comparison figure.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import barycon as bc
from barycon.conversions import fbcg_akino22
from barycon.plotting import DATASET_COLORS

OUT = Path(sys.argv[2] if len(sys.argv) > 2 else "docs/figures")
OUT.mkdir(parents=True, exist_ok=True)
LEGACY = Path(sys.argv[1]) if len(sys.argv) > 1 else None
plt.rcParams.update({"font.size": 11, "axes.labelsize": 12})

from astropy.cosmology import Planck18
FB = Planck18.Ob0 / (Planck18.Om0 - Planck18.Onu0)  # baryons / cold matter, as in the legacy plots

col = bc.load()

# Paper figures: full text width of the JOSS template (~13.6 cm) and 10 pt text, so that
# labels, ticks and legends print at the size of the body text.
TEXTWIDTH = 5.36  # inches
PAPER_RC = {"font.size": 10, "axes.labelsize": 10, "xtick.labelsize": 10, "ytick.labelsize": 10,
            "legend.fontsize": 10, "legend.title_fontsize": 10, "mathtext.fontset": "dejavusans"}


def save(fig, name):
    fig.savefig(OUT / name, dpi=300)
    if Path("paper").is_dir():
        fig.savefig(Path("paper") / name, dpi=300)  # figures of the JOSS paper


# ------------------------------------------------------------------ gas, r500c
gas = col.filter(quantity="f_gas", delta="500c", region="sphere").homogenize("Planck18")
with plt.rc_context(PAPER_RC):
    fig, ax = plt.subplots(figsize=(TEXTWIDTH, 5.6), layout="constrained")
    gas.plot("f_gas", ax=ax, universal=FB, legend_ncol=4)
    ax.set_xlim(10**12.5, 3e15)
    ax.set_ylim(0, 0.25)
    save(fig, "gas_fractions_500c.png")

# --------------------------------------------------------------- stars, r500c
stars = col.filter(quantity="f_star").homogenize("Planck18", imf="chabrier")
st500 = []
for e in stars:
    if e.aperture["delta"] == "500c":
        st500.append(e)
    else:  # Reyes2012: isolated central galaxies -> essentially all stars in the central galaxy
        st500.append(e.to_delta("500c", f_bcg=lambda M: 0.99 * np.ones_like(M)))
st500 = bc.Collection(st500)
with plt.rc_context(PAPER_RC):
    fig, ax = plt.subplots(figsize=(TEXTWIDTH, 4.6), layout="constrained")
    st500.plot("f_star", ax=ax, proxy_legend_loc="lower center", legend_ncol=4)
    ax.set_yscale("log")
    ax.set_xlim(5e10, 3e15)
    ax.set_ylim(2e-3, 0.1)
    save(fig, "stellar_fractions_500c.png")

# ------------------------------------------------------------ gas versus z
fig, ax = plt.subplots(figsize=(10, 7.5), constrained_layout=True)
gas.filter(data_type="points").plot("f_gas", ax=ax, x="z", universal=FB)
ax.set_ylim(0, 0.25)
fig.savefig(OUT / "gas_fractions_z.png", dpi=160)

# ----------------------------------------------------------- fitting example
hse = gas.filter(ids=["Vikhlinin2006", "Arnaud2007", "Sun2009", "Wicker2023", "Ettori2010"])
fit_ls = hse.fit("saturating")
fit_sc = hse.fit("saturating", scatter=True)
print(fit_ls.summary(), fit_sc.summary(), sep="\n")
fig, ax = plt.subplots(figsize=(10, 7.5), constrained_layout=True)
hse.plot("f_gas", ax=ax, universal=FB)
fit_ls.plot(ax=ax, color="#52514e", ls="--", label="saturating, least squares")
fit_sc.plot(ax=ax, color="#0b0b0b", label=rf"saturating + {100 * fit_sc.scatter:.0f}% intrinsic scatter")
ax.legend(fontsize=7, ncol=2, frameon=False, loc="lower right")
ax.set_xlim(10 ** 12.8, 3e15)
ax.set_ylim(0, 0.25)
ax.set_title("Fitting hydrostatic samples (Planck18): with and without intrinsic scatter", fontsize=11)
fig.savefig(OUT / "fit_example.png", dpi=160)

# --------------------------------------------- explicit aperture choice: r200c
fig, ax = plt.subplots(figsize=(10, 7.5), constrained_layout=True)
g200 = gas.filter(data_type="points").to_delta("200c", eta_gas=0.3)
g200 = bc.Collection(list(g200) + [col["Popesso2024_200c"].to_cosmology("Planck18"),
                                   col["Khalil2026_200c"].to_cosmology("Planck18")])
g200.plot("f_gas", ax=ax, universal=FB)
ax.set_xlim(10**12.5, 4e15)
ax.set_ylim(0, 0.25)
ax.set_title(r"Gas fractions moved to $r_{200c}$ assuming $f_{\rm gas}(<r)\propto r^{0.3}$ "
             "(a user choice, not a default)", fontsize=11)
fig.savefig(OUT / "gas_fractions_200c_eta0.3.png", dpi=160)

# --------------------------------------------------- legacy vs fixed (optional)
if LEGACY is not None:
    fig, axs = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    old_c, new_c = "#8a8984", "#2a78d6"

    # (a) Kravtsov: wrong M500 column
    d = np.loadtxt(LEGACY / "Kravtsov_2018.dat", skiprows=2)
    ms = (d[:, 7] + d[:, 9]) * 1e12
    ax = axs[0, 0]
    ax.plot(d[:, 3] * 1e14 * 0.7, ms / (d[:, 3] * 1e14), "o", color=old_c, ms=8, label="legacy (column 3 as M500)")
    ax.plot(d[:, 1] * 1e14 * 0.7, ms / (d[:, 1] * 1e14), "s", color=new_c, ms=8, label="fixed (column 1 = M500)")
    for i in range(len(d)):
        ax.annotate("", (d[i, 1] * 0.7e14, ms[i] / (d[i, 1] * 1e14)), (d[i, 3] * 0.7e14, ms[i] / (d[i, 3] * 1e14)),
                    arrowprops=dict(arrowstyle="->", color="#52514e", lw=0.8))
    ax.set(xscale="log", xlabel=r"$M_{500c}\ [h_{70}^{-1}{\rm M}_\odot]$", ylabel=r"$f_\star$",
           title="(a) Kravtsov+18: column shift")
    ax.legend(frameon=False, fontsize=9)

    # (b) Ettori: dark-matter M200 -> M500 vs R500-based total mass
    m2 = np.loadtxt(LEGACY / "Ettori_2011_masses_method2.dat", skiprows=2)
    from scipy.optimize import brentq

    def c500(c2):
        f = lambda c5: (c2 / c5) ** 3 * (np.log1p(c5) - c5 / (1 + c5)) / (np.log1p(c2) - c2 / (1 + c2)) - 2.5
        return brentq(f, 0.05, c2)
    c5 = np.array([c500(c) for c in m2[:, 3]])
    M_old = 2.5 * (c5 / m2[:, 3]) ** 3 * m2[:, 6] * 1e14
    et = col["Ettori2010"].table
    ax = axs[0, 1]
    ax.plot(M_old * 0.7, et["f_gas"], "o", color=old_c, ms=7, label="legacy (NFW from dark-matter $M_{200}$)")
    ax.plot(np.asarray(et["M_500c"]) * 0.7, et["f_gas"], "s", color=new_c, ms=7,
            label=r"fixed ($500\rho_c(z)\,4\pi r_{500}^3/3$)")
    ax.set(xscale="log", xlabel=r"$M_{500c}\ [h_{70}^{-1}{\rm M}_\odot]$", ylabel=r"$f_{\rm gas}$",
           title=f"(b) Ettori+10: masses were {100*(1-np.median(M_old/et['M_500c'])):.0f}% too low")
    ax.legend(frameon=False, fontsize=9)

    # (c) Popesso band
    M = np.geomspace(1e13, 1e15, 50)
    ax = axs[1, 0]
    lo_old = (2.23 - 0.18) * 1e-7 * M ** (0.39 - 0.02)
    hi_old = (2.23 + 0.18) * 1e-7 * M ** (0.39 + 0.02)
    ax.fill_between(M, lo_old, hi_old, color=old_c, alpha=0.3, lw=0, label="legacy band (no pivot)")
    p16, p50, p84 = col["Popesso2024"].evaluate(M, "f_gas", n_samples=4000)
    ax.fill_between(M, p16, p84, color=DATASET_COLORS["Popesso2024"], alpha=0.4, lw=0, label="fixed band (pivot $10^{14}$)")
    ax.plot(M, p50, color=DATASET_COLORS["Popesso2024"], lw=2)
    ax.set(xscale="log", xlabel=r"$M_{500c}\ [{\rm M}_\odot]$", ylabel=r"$f_{\rm gas}$",
           title="(c) Popesso+24: error band")
    ax.legend(frameon=False, fontsize=9)

    # (d) size of each fix
    fixes = {
        "Kravtsov18 $f_\\star$ (A2142)": 100 * (11.96 / 15.60 - 1) * -1,
        "Kravtsov18 $f_\\star$ (MKW4)": 100 * (0.56 / 0.46 - 1) * -1,
        "Ettori10 $M_{500}$": 100 * (1 / np.median(M_old / et["M_500c"]) - 1),
        "Vikhlinin06 $f_{\\rm gas}$ (h 0.72→0.71)": 100 * ((0.71 / 0.72) ** 1.5 - 1),
        "Zhang11 $f_\\star$ ($h^{-3/2}$→$h^{-1}$)": 100 * ((0.7 / 0.6766) ** 1.0 / (0.7 / 0.6766) ** 1.5 - 1),
        "Giodini09 $f_\\star$ (Salpeter→Chabrier)": 100 * (0.63 - 1),
    }
    ax = axs[1, 1]
    y = np.arange(len(fixes))
    ax.barh(y, list(fixes.values()), color="#2a78d6", height=0.5)
    ax.axvline(0, color="#52514e", lw=0.8)
    ax.set_yticks(y, list(fixes.keys()), fontsize=9)
    for yi, v in zip(y, fixes.values()):
        ax.text(v + (1 if v >= 0 else -1), yi, f"{v:+.1f}%", va="center", ha="left" if v >= 0 else "right",
                fontsize=9, color="#0b0b0b")
    ax.set(xlabel="change w.r.t. legacy value [%]", title="(d) size of each correction", xlim=(-50, 40))
    for a in axs.flat:
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
    fig.savefig(OUT / "legacy_vs_fixed.png", dpi=160)

print("figures written to", OUT)

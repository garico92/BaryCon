#!/usr/bin/env python
"""Build the BaryCon entries from the plain-text tables in ``data/``.

Usage
-----
    python tools/migrate_from_text.py data src/barycon/data

``data/`` holds the tables as transcribed from the papers (merged from the
Aricò and Douspis repositories, Nov 2026).

Every non-trivial choice is written in the entry metadata (``caveats``,
``changelog``) so that it is reviewable in the pull request that introduces it.
Numbers are *never* rescaled here: entries keep the units/cosmology of the
original paper. Homogenisation happens at run time, explicitly.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import yaml
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table

TODAY = "2026-11-06"
CONTRIB = ["G. Aricò", "M. Douspis"]


def _rows(path):
    return [l for l in open(path, encoding="utf-8") if l.strip() and not l.lstrip().startswith("#")]


def _meta(**kw):
    base = dict(kind="observation", extracted=False, status="unverified", contributors=CONTRIB,
                added=TODAY, caveats=[], changelog=[])
    base.update(kw)
    return base


def _write_ecsv(table, meta, out):
    from barycon.schema import validate_meta
    validate_meta(meta)
    t = table.copy()
    # human-readable files: masses stored in units of 1e14 Msun (unit recorded in the ECSV header)
    for c in t.colnames:
        if c.startswith("M_") and t[c].dtype.kind == "f":
            t[c] = np.round(np.asarray(t[c]) / 1e14, 6)
            t[c].unit = u.Unit("1e14 solMass")
        elif c.startswith("f_") and t[c].dtype.kind == "f":
            t[c] = np.round(np.asarray(t[c]), 6)
    t.meta = meta
    t.write(out / f"{meta['id']}.ecsv", format="ascii.ecsv", overwrite=True)
    print(f"  wrote {meta['id']}.ecsv ({len(t)} rows)")


def _write_yaml(components, meta, out):
    from barycon.schema import validate_meta
    validate_meta(meta)
    with open(out / f"{meta['id']}.yaml", "w") as fh:
        yaml.safe_dump({"meta": meta, "components": components}, fh, sort_keys=False, width=100)
    print(f"  wrote {meta['id']}.yaml ({len(components)} components)")


def _ap(delta="500c", region="sphere", r_inner=None, r_outer=1.0):
    return dict(delta=delta, region=region, r_inner=r_inner, r_outer=r_outer)


def _mass(column, definition, proxy, details="", pM=1.0):
    return dict(column=column, definition=definition, proxy=proxy, proxy_details=details,
                unit="Msun", distance_exponent=pM)


def _ratio_err(num, dnum, den, dden):
    """1-sigma error of num/den for uncorrelated errors."""
    r = num / den
    return r * np.sqrt((dnum / num) ** 2 + (dden / den) ** 2)


# ---------------------------------------------------------------------------
def arnaud2007(src, out):
    d = np.loadtxt(src / "Arnaud_2007.dat", skiprows=2, dtype=str)
    f = d[:, 1:].astype(float)
    t = Table({
        "name": d[:, 0], "z": f[:, 0],
        "M_500c": f[:, 3] * 1e14, "M_500c_hi": f[:, 4] * 1e14, "M_500c_lo": -f[:, 5] * 1e14,
        "M_gas_500c": f[:, 6] * 1e13,
        "f_gas": f[:, 12], "f_gas_hi": f[:, 13], "f_gas_lo": -f[:, 14],
    })
    _write_ecsv(t, _meta(
        id="Arnaud2007",
        reference="Arnaud, M., Pointecouteau, E., Pratt, G. W. 2007, A&A, 474, L37",
        bibcode="2007A&A...474L..37A", doi="10.1051/0004-6361:20078541", arxiv="0709.1561",
        url="https://www.aanda.org/articles/aa/pdf/2007/42/aa8541-07.pdf",
        description="10 nearby relaxed clusters observed with XMM-Newton; hydrostatic M500 from NFW fits "
                    "to the mass profiles (calibration of the M500-YX relation).",
        data_type="points", quantities=["f_gas"], aperture=_ap(),
        mass=_mass("M_500c", "500c", "hydrostatic", "XMM-Newton, NFW fit to HSE mass profile"),
        cosmology=dict(H0=70.0, Om0=0.3, Ode0=0.7),
        redshift=dict(min=0.03, max=0.16, effective=None),
        selection="relaxed, morphologically regular clusters",
        table_source="Table 1",
        changelog=["migrated from legacy Arnaud_2007.dat (was labelled 'Arnaud 2017' in legacy plots)"],
    ), out)


def bulbul2024(src, out):
    d = np.loadtxt(src / "Bulbul_2024.dat", skiprows=2)
    bad = np.any(d < 0, axis=1)
    d = d[~bad]
    zero = d[:, 3] <= 0
    d = d[~zero]
    M, Mlo_b, Mhi_b = d[:, 0] * 1e13, d[:, 1] * 1e13, d[:, 2] * 1e13
    f, flo_b, fhi_b = d[:, 3], d[:, 4], d[:, 5]
    t = Table({"M_500c": M, "M_500c_lo": M - Mlo_b, "M_500c_hi": Mhi_b - M,
               "f_gas": f, "f_gas_lo": f - flo_b, "f_gas_hi": fhi_b - f})
    _write_ecsv(t, _meta(
        id="Bulbul2024",
        reference="Bulbul, E., Liu, A., Kluge, M., et al. 2024, A&A, 685, A106",
        bibcode="2024A&A...685A.106B", doi="10.1051/0004-6361/202348264", arxiv="2402.08452",
        url="https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A+A/685/A106",
        description="eRASS1 (SRG/eROSITA) primary cluster catalogue, western Galactic hemisphere. "
                    "M500 from the count-rate/luminosity-mass relation calibrated with weak lensing; "
                    "gas masses from the deprojected surface brightness.",
        data_type="points", quantities=["f_gas"], aperture=_ap(),
        mass=_mass("M_500c", "500c", "xray_scaling_relation",
                   "count-rate-mass relation, WL calibrated (Ghirardini+24, Grandis+24, Kleinebreil+24)"),
        cosmology=dict(H0=70.0, Om0=0.3, Ode0=0.7,
                       note="VERIFY: catalogue fiducial cosmology; legacy script used h=0.6715 (copy-paste)"),
        redshift=dict(min=0.003, max=1.32, effective=0.2),
        selection="X-ray (extent) selected; strongly incomplete at low mass",
        table_source="CDS J/A+A/685/A106 (catalogue)",
        status="needs_check",
        caveats=[
            "file had no redshifts/names: re-extract from CDS with z, Lext and IN_ZVLIM flags",
            "legacy header said [1e14 Msun]; values are in 1e13 Msun (plot_stellar_fractions.py used 1e14)",
            "columns 2-3 and 5-6 were lower/upper *bounds*, converted here to errors",
            "f_gas and M500 both derive from the same X-ray data -> correlated errors",
            f"{int(bad.sum())} rows with -1 sentinels and {int(zero.sum())} rows with f_gas=0 removed",
        ],
    ), out)


def chiu2018(src, out):
    d = np.genfromtxt(src / "Chiu_2018.dat", dtype=None, names=True, encoding="utf-8", skip_header=17)
    M, dM = d["M500"] * 1e14, d["M500_err"] * 1e14
    Mg, dMg = d["MICM"] * 1e13, d["MICM_err"] * 1e13
    Ms, dMs = d["Mstar"] * 1e12, d["Mstar_err"] * 1e12
    fg, fs = Mg / M, Ms / M
    efg, efs = _ratio_err(Mg, dMg, M, dM), _ratio_err(Ms, dMs, M, dM)
    t = Table({"name": d["Name"], "z": d["z"], "R_500c_arcmin": d["R500"],
               "M_500c": M, "M_500c_lo": dM, "M_500c_hi": dM,
               "M_gas_500c": Mg, "M_star_500c": Ms,
               "f_gas": fg, "f_gas_lo": efg, "f_gas_hi": efg,
               "f_star": fs, "f_star_lo": efs, "f_star_hi": efs})
    _write_ecsv(t, _meta(
        id="Chiu2018",
        reference="Chiu, I., Mohr, J. J., McDonald, M., et al. 2018, MNRAS, 478, 3072",
        bibcode="2018MNRAS.478.3072C", doi="10.1093/mnras/sty1284", arxiv="1711.00917",
        url="https://academic.oup.com/mnras/article/478/3/3072/4996803",
        description="91 SPT-SZ clusters (0.2<z<1.25). M500 from the SZE observable (de Haan+16 relation, "
                    "Eddington/Malmquist corrected); ICM masses from Chandra within R500(SZ); stellar masses "
                    "from SED fits to DES + WISE/Spitzer photometry.",
        data_type="points", quantities=["f_gas", "f_star"], aperture=_ap(),
        mass=_mass("M_500c", "500c", "sz_scaling_relation", "SPT xi-M relation (de Haan et al. 2016)"),
        cosmology=dict(H0=68.0, Om0=0.304, Ode0=0.696),
        stellar=dict(imf="chabrier", includes_bcg=True, includes_icl=False),
        redshift=dict(min=float(np.min(d["z"])), max=float(np.max(d["z"])), effective=None),
        selection="SZE selected (xi>~5)",
        table_source="Table 2 (journal version); redshifts from Simbad (M. Douspis, Oct 2026)",
        caveats=["redshifts retrieved from Simbad, not from the SPT catalogue (Bleem+15): cross-check",
                 "fraction errors propagated assuming independent errors (they are correlated)"],
        changelog=["metadata (cosmology, IMF, units) checked against the paper 2026-10-31",
                   "redshifts added (Douspis) 2026-11"],
    ), out)

def ettori2010(src, out):
    e = np.loadtxt(src / "Ettori_2011.dat", skiprows=2)
    rows = _rows(src / "Ettori_2011_masses_method2.dat")
    names = [r.split()[0] for r in rows]
    z = np.array([float(r.split()[1]) for r in rows])
    assert len(e) == len(z) == 44
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    rhoc = cosmo.critical_density(z).to(u.Msun / u.Mpc**3).value
    r500, dr500 = e[:, 8] / 1e3, e[:, 9] / 1e3  # method 2, Mpc
    M500 = 500 * rhoc * 4 / 3 * np.pi * r500**3
    dM500 = 3 * M500 * dr500 / r500
    t = Table({"name": names, "z": z, "R_500c_kpc": e[:, 8],
               "M_500c": M500, "M_500c_lo": dM500, "M_500c_hi": dM500,
               "f_gas": e[:, 10], "f_gas_lo": e[:, 11], "f_gas_hi": e[:, 11],
               "f_gas_method1": e[:, 4], "f_gas_method1_err": e[:, 5]})
    _write_ecsv(t, _meta(
        id="Ettori2010",
        reference="Ettori, S., Gastaldello, F., Leccardi, A., et al. 2010, A&A, 524, A68 "
                  "(values from the Corrigendum, A&A 526, C1, 2011)",
        bibcode="2010A&A...524A..68E", doi="10.1051/0004-6361/201015271", arxiv="1009.3266",
        url="https://arxiv.org/abs/1009.3266",
        description="44 X-ray luminous clusters at 0.1<z<0.3 observed with XMM-Newton; hydrostatic NFW "
                    "reconstruction. f_gas within R500 from method 2 (method 1 kept as extra column).",
        data_type="points", quantities=["f_gas"], aperture=_ap(),
        mass=_mass("M_500c", "500c", "hydrostatic",
                   "M500 = 500 rho_c(z) 4pi/3 R500^3 with R500 (method 2) from Corrigendum Table 1"),
        cosmology=dict(H0=70.0, Om0=0.3, Ode0=0.7),
        redshift=dict(min=0.09, max=0.31, effective=None),
        table_source="Corrigendum (A&A 526, C1) Table 1; names and redshifts from Table 1 of the original paper",
        changelog=[
            "Legacy code converted M200 of Table 2 to M500 with NFW. Table 2 M200 is the *dark* mass "
            "(Mtot - Mgas), so M500,tot was underestimated by ~18%. Now computed from R500 and z.",
            "names/redshifts independently transcribed by Aricò and Douspis: identical",
        ],
    ), out)

def giodini2009(src, out):
    d = np.loadtxt(src / "Giodini_2009.dat", skiprows=2)
    t = Table({"M_500c": d[:, 0], "f_star": d[:, 1], "f_star_lo": d[:, 2], "f_star_hi": d[:, 2],
               "f_gas": d[:, 3], "f_gas_lo": d[:, 4], "f_gas_hi": d[:, 4]})
    _write_ecsv(t, _meta(
        id="Giodini2009",
        reference="Giodini, S., Pierini, D., Finoguenov, A., et al. 2009, ApJ, 703, 982",
        bibcode="2009ApJ...703..982G", doi="10.1088/0004-637X/703/1/982", arxiv="0904.0448",
        url="https://arxiv.org/abs/0904.0448",
        description="91 COSMOS X-ray groups (0.1<z<1) + 27 nearby clusters (Lin+03); mean stellar "
                    "fractions in five M500 bins. M200 from a WL-calibrated LX-M relation converted to "
                    "M500 with c=5.",
        data_type="binned", quantities=["f_star", "f_gas"], derived_quantities=["f_gas"], aperture=_ap(),
        mass=_mass("M_500c", "500c", "xray_scaling_relation", "LX-M200 calibrated with COSMOS WL (Leauthaud+10)"),
        cosmology=dict(H0=72.0, Om0=0.258, Ode0=0.742),
        stellar=dict(imf="salpeter", includes_bcg=True, includes_icl=False),
        redshift=dict(min=0.0, max=1.0, effective=0.3),
        table_source="Table 3 (bin centres) and text",
        caveats=["f_gas is NOT measured: it is a mean fit to Vikhlinin+06, Arnaud+07 and Sun+09 "
                 "(already in this collection) -> excluded by default from f_gas analyses",
                 "Salpeter IMF: multiply by ~0.63 for Chabrier (use .to_imf())"],
        changelog=["metadata (cosmology, IMF, f_gas origin) checked against the paper 2026-10-31"],
        overlaps_with=["Vikhlinin2006", "Arnaud2007", "Sun2009"],
    ), out)


def gonzalez2013(src, out):
    rows = [r.split() for r in _rows(src / "Gonzalez_2013.dat")]
    names = [r[0] for r in rows]
    f = np.array([[float(x) for x in r[1:]] for r in rows])  # z, fgas, e, fstar, e, fb, e, M500, eM, r500, er
    t = Table({"name": names, "z": f[:, 0],
               "M_500c": f[:, 7] * 1e14, "M_500c_lo": f[:, 8] * 1e14, "M_500c_hi": f[:, 8] * 1e14,
               "R_500c_Mpc": f[:, 9],
               "f_gas": f[:, 1], "f_gas_lo": f[:, 2], "f_gas_hi": f[:, 2],
               "f_star": f[:, 3], "f_star_lo": f[:, 4], "f_star_hi": f[:, 4],
               "f_b": f[:, 5], "f_b_lo": f[:, 6], "f_b_hi": f[:, 6]})
    _write_ecsv(t, _meta(
        id="Gonzalez2013",
        reference="Gonzalez, A. H., Sivanandam, S., Zabludoff, A. I., Zaritsky, D. 2013, ApJ, 778, 14",
        bibcode="2013ApJ...778...14G", doi="10.1088/0004-637X/778/1/14", arxiv="1309.3565",
        url="https://arxiv.org/abs/1309.3565",
        description="12 clusters/groups at z~0.1 with XMM gas masses and stellar masses including BCG+ICL. "
                    "Values from the Planck-cosmology tables (the paper also gives WMAP7 tables).",
        data_type="points", quantities=["f_gas", "f_star", "f_b"], aperture=_ap(),
        mass=_mass("M_500c", "500c", "xray_scaling_relation",
                   "M500 from the X-ray temperature in (0.15-0.5) R500 via the Vikhlinin+09 M-T relation "
                   "(as described by Kravtsov+18)"),
        cosmology=dict(H0=67.15, Om0=0.3175, Ode0=0.6825,
                       note="Planck 2013; VERIFY exact values used in Tables 6-7 (legacy h=0.6715)"),
        stellar=dict(imf=None, includes_bcg=True, includes_icl=True),
        redshift=dict(min=float(f[:, 0].min()), max=float(f[:, 0].max()), effective=None),
        table_source="Tables 6 and 7 (Planck cosmology); redshifts from Simbad/Sanderson+13 (Douspis, Oct 2026)",
        status="needs_check",
        caveats=["IMF not recorded yet: add before IMF homogenisation",
                 "stellar mass includes ICL (unlike most other entries)"],
        overlaps_with=["Sanderson2013"],
    ), out)

def grandis2024(src, out):
    """Table 1 of Grandis+24: one point per survey at its pivot mass, with the local slope."""
    d = np.genfromtxt(src / "Grandis_2024.dat", comments="#", dtype=None, encoding=None)
    proxy = {"SPT": ("sz_scaling_relation", "SPT SZ (Chiu+18)"),
             "eFEDS": ("xray_scaling_relation", "count-rate-mass relation, WL calibrated"),
             "HSCXXL": ("weak_lensing", "HSC weak lensing (Akino+22)"),
             "WtG": ("weak_lensing", "Weighing the Giants weak lensing")}
    surveys = list(dict.fromkeys(str(r[0]) for r in d))
    nan = np.full(len(surveys), np.nan)
    cols = {k: nan.copy() for k in ("M_500c", "f_gas", "f_gas_lo", "f_gas_hi", "f_star", "f_star_lo", "f_star_hi",
                                    "slope_f_gas", "slope_f_gas_lo", "slope_f_gas_hi",
                                    "slope_f_star", "slope_f_star_lo", "slope_f_star_hi")}
    for survey, comp, Mp, f, df, s, ds in d:
        i = surveys.index(str(survey))
        q = "f_gas" if comp == "ICM" else "f_star"
        cols["M_500c"][i] = float(Mp)
        cols[q][i] = float(f)
        if df > 0:
            cols[q + "_lo"][i] = cols[q + "_hi"][i] = float(df)
        cols["slope_" + q][i] = float(s)
        cols["slope_" + q + "_lo"][i] = cols["slope_" + q + "_hi"][i] = float(ds)
    t = Table({"name": surveys, "mass_proxy": [proxy[x][0] for x in surveys],
               "mass_proxy_details": [proxy[x][1] for x in surveys], **cols})
    _write_ecsv(t, _meta(
        id="Grandis2024",
        reference="Grandis, S., et al. 2024, MNRAS (compilation in their Table 1), arXiv:2309.02920",
        bibcode=None, doi="10.1093/mnras/stae259", arxiv="2309.02920", url="https://arxiv.org/abs/2309.02920",
        description="Compilation of gas and stellar fractions at the pivot mass of power-law fits from SPT "
                    "(Chiu+18), eFEDS, HSC-XXL (Akino+22) and WtG, with the local logarithmic slope.",
        data_type="points", quantities=["f_gas", "f_star"], aperture=_ap(),
        mass=_mass("M_500c", "500c", "mixed", "differs per survey: see the mass_proxy column"),
        cosmology=dict(H0=67.0, Om0=0.3, Ode0=0.7, note="VERIFY: legacy used h=0.67 for all rows"),
        stellar=dict(imf=None, includes_bcg=True, includes_icl=False),
        redshift=dict(min=0.0, max=1.3, effective=0.3),
        table_source="Table 1",
        status="needs_check",
        caveats=["each survey row was derived under its own cosmology: check whether Grandis+24 homogenised",
                 "duplicates information in Chiu2018 and Akino2022 -> do not fit jointly",
                 "one point per survey at the pivot mass; slope_* columns give the local logarithmic slope",
                 "WtG row has no error and no slope",
                 "IMF to be recorded"],
        overlaps_with=["Chiu2018", "Akino2022", "Mantz2022"],
    ), out)


def kravtsov2018(src, out):
    d = np.loadtxt(src / "Kravtsov_2018.dat", skiprows=2)  # col 0 = z (added by Douspis)
    names = ["A2142", "A2029", "A85", "A1795", "MKW3s", "A2052", "A1991", "MKW4", "RXJ1159+5531"]
    z = d[:, 0]
    d = d[:, 1:]
    M, dM = d[:, 1] * 1e14, d[:, 2] * 1e14
    Mbcg, dMbcg = d[:, 7] * 1e12, d[:, 8] * 1e12
    Msat, dMsat = d[:, 9] * 1e12, d[:, 10] * 1e12
    Ms = Mbcg + Msat
    dMs = np.hypot(dMbcg, dMsat)
    t = Table({"name": names, "z": z, "R_500c_kpc": d[:, 0],
               "M_500c": M, "M_500c_lo": dM, "M_500c_hi": dM,
               "M_gas_500c_tentative": d[:, 3] * 1e13,
               "M_star_bcg": Mbcg, "M_star_sat": Msat,
               "f_star": Ms / M, "f_star_lo": _ratio_err(Ms, dMs, M, dM), "f_star_hi": _ratio_err(Ms, dMs, M, dM),
               "f_bcg": Mbcg / M, "f_bcg_lo": _ratio_err(Mbcg, dMbcg, M, dM),
               "f_bcg_hi": _ratio_err(Mbcg, dMbcg, M, dM)})
    _write_ecsv(t, _meta(
        id="Kravtsov2018",
        reference="Kravtsov, A. V., Vikhlinin, A. A., Meshcheryakov, A. V. 2018, Astronomy Letters, 44, 8",
        bibcode="2018AstL...44....8K", doi="10.1134/S1063773717120015", arxiv="1401.7329",
        url="https://arxiv.org/abs/1401.7329",
        description="9 nearby clusters with deep photometry of BCG outskirts and Chandra total masses; "
                    "stellar mass in BCG (extended profile) and satellites within R500.",
        data_type="points", quantities=["f_star", "f_bcg"], aperture=_ap(),
        mass=_mass("M_500c", "500c", "xray_scaling_relation",
                   "Chandra Y_X proxy (Vikhlinin+09); A1991, MKW4, RXJ1159 from the same Chandra analysis"),
        cosmology=dict(H0=70.0, Om0=0.27, Ode0=0.73),
        stellar=dict(imf="chabrier", includes_bcg=True, includes_icl=True),
        redshift=dict(min=float(z.min()), max=float(z.max()), effective=None),
        table_source="Table 1; redshifts from NED/literature (Douspis, Oct 2026)",
        status="needs_check",
        caveats=[
            "LEGACY BUG FIXED: M500 is the 2nd numeric column (11.96 for A2142), consistent with R500 and z "
            "to <7% (Douspis' file header reaches the same conclusion). The 4th column (15.60) was used as "
            "M500 by both legacy codes; it is most likely M_gas,500 in 1e13 Msun (f_gas ~0.08-0.13). "
            "Legacy f_star biased by -18%..+30% depending on mass.",
            "M_gas_500c_tentative kept until checked against the journal table",
        ],
        changelog=["column mapping and cosmology/IMF checked against arXiv:1401.7329 2026-10-31",
                   "mass proxy = Y_X (paper text) 2026-11-06"],
    ), out)

def landry2013(src, out):
    names, rows = [], []
    cur = None
    for line in open(src / "Landry_2013.dat", encoding="utf-8"):
        m = re.search(r"underline\{([^}]*)\}", line)
        if m:
            cur = m.group(1)
            continue
        if line.strip().startswith("Vikh"):
            names.append(cur)
            rows.append([float(x) for x in line.split()[1:]])
    a = np.array(rows)
    # sanity: f_gas = M_gas/M_tot for both apertures (validates the inferred column meaning)
    for off in (0, 12):
        ratio = a[:, off + 3] * 1e13 / (a[:, off + 6] * 1e14)
        assert np.allclose(ratio, a[:, off + 9], rtol=0.06), "Landry columns not understood"
    for delta, off in (("500c", 12), ("2500c", 0)):
        t = Table({"name": names,
                   f"M_{delta}": a[:, off + 6] * 1e14, f"M_{delta}_hi": a[:, off + 7] * 1e14,
                   f"M_{delta}_lo": a[:, off + 8] * 1e14, f"M_gas_{delta}": a[:, off + 3] * 1e13,
                   "f_gas": a[:, off + 9], "f_gas_hi": a[:, off + 10], "f_gas_lo": a[:, off + 11]})
        _write_ecsv(t, _meta(
            id="Landry2013" if delta == "500c" else "Landry2013_2500c",
            reference="Landry, D., Bonamente, M., Giles, P., et al. 2013, MNRAS, 433, 2790",
            bibcode="2013MNRAS.433.2790L", doi=None, arxiv="1211.4626",
            url="https://arxiv.org/abs/1211.4626",
            description="Complete sample of X-ray luminous clusters at 0.15<z<0.3 observed with Chandra; "
                        "hydrostatic masses using the Vikhlinin+06 temperature model ('Vikh' rows; the "
                        "polytropic 'Poly' rows are not migrated).",
            data_type="points", quantities=["f_gas"], aperture=_ap(delta),
            mass=_mass(f"M_{delta}", delta, "hydrostatic", "Chandra, Vikhlinin+06 model"),
            cosmology=dict(H0=70.0, Om0=0.3, Ode0=0.7, note="VERIFY"),
            redshift=dict(min=0.15, max=0.3, effective=0.22),
            table_source="Table 11",
            status="needs_check",
            caveats=["legacy file had no column header: meaning inferred (Mgas[1e13]/Mtot[1e14] reproduces "
                     "f_gas to <6%); radius column units (arcsec?) unverified",
                     "not plotted by the legacy scripts"],
        ), out)


def mantz2022(src, out):
    d = np.genfromtxt(src / "Mantz_2022.dat", dtype=None, names=True, encoding="utf-8", skip_header=16)
    t = Table({"name": d["name"], "z": d["z"],
               "M_2500c": d["M2500"] * 1e14, "M_2500c_lo": d["M2500_err"] * 1e14,
               "M_2500c_hi": d["M2500_err"] * 1e14,
               "M_shell": d["Mtot"] * 1e14, "M_gas_shell": d["Mgas"] * 1e13,
               "f_gas": d["fgas"], "f_gas_lo": d["fgas_err"], "f_gas_hi": d["fgas_err"]})
    _write_ecsv(t, _meta(
        id="Mantz2022",
        reference="Mantz, A. B., Morris, R. G., Allen, S. W., et al. 2022, MNRAS, 510, 131",
        bibcode="2022MNRAS.510..131M", doi=None, arxiv="2111.09343",
        url="https://arxiv.org/abs/2111.09343",
        description="Massive, dynamically relaxed clusters (Chandra) used for f_gas cosmology. Gas fraction "
                    "measured in the SHELL 0.8-1.2 r2500c; mass column is M2500c (0-1 r2500c).",
        data_type="points", quantities=["f_gas"],
        aperture=_ap("2500c", "shell", 0.8, 1.2),
        mass=_mass("M_2500c", "2500c", "hydrostatic", "Chandra HSE"),
        cosmology=dict(H0=70.0, Om0=0.3, Ode0=0.7, note="reference cosmology; f_gas ∝ d_A^{3/2}"),
        redshift=dict(min=0.0, max=1.2, effective=None),
        table_source="Table 1",
        caveats=["shell measurement: NOT comparable with f_gas(<r500c); never convert with to_delta()",
                 "relaxed-cluster selection"],
    ), out)


def mulroy2019(src, out):
    d = np.genfromtxt(src / "Mulroy_2019.dat", dtype=None, names=True, encoding="utf-8", skip_header=26)
    M, Mp, Mm = d["M_post"] * 1e14, d["M_post_err_plus"] * 1e14, d["M_post_err_minus"] * 1e14
    Mg, Mgp, Mgm = d["Mgas"] * 1e14, d["Mgas_err_plus"] * 1e14, d["Mgas_err_minus"] * 1e14
    f = Mg / M
    t = Table({"name": d["name"], "z": d["z"],
               "M_500c": M, "M_500c_lo": Mm, "M_500c_hi": Mp,
               "M_500c_WL": d["M_wl"] * 1e14, "M_gas_500c": Mg,
               "f_gas": f, "f_gas_lo": f * np.hypot(Mgm / Mg, Mp / M), "f_gas_hi": f * np.hypot(Mgp / Mg, Mm / M)})
    _write_ecsv(t, _meta(
        id="Mulroy2019",
        reference="Mulroy, S. L., Farahi, A., Evrard, A. E., et al. 2019, MNRAS, 484, 60",
        bibcode="2019MNRAS.484...60M", doi="10.1093/mnras/sty3484", arxiv="1901.11276",
        url="https://arxiv.org/abs/1901.11276",
        description="LoCuSS: 41 X-ray luminous clusters at 0.15<z<0.3 with Subaru weak-lensing masses, "
                    "Chandra/XMM gas masses and NIR stellar masses.",
        data_type="points", quantities=["f_gas"], aperture=_ap(),
        mass=_mass("M_500c", "500c", "weak_lensing", "M_post (posterior WL masses); M_WL kept as column"),
        cosmology=dict(H0=70.0, Om0=0.3, Ode0=0.7),
        redshift=dict(min=0.15, max=0.31, effective=None),
        table_source="Tables 1-2",
        caveats=["uses posterior ('M_post') masses as in legacy code; consider M_WL for an independent mass",
                 "check the aperture of Mgas (R500 from WL?)"],
    ), out)


def reyes2012(src, out):
    d = np.loadtxt(src / "Reyes_2012.dat", skiprows=2)
    Ms = 10 ** d[:, 0]
    p2, p16, p50, p84, p97 = d[:, 1], d[:, 2], d[:, 3], d[:, 4], d[:, 5]
    M = Ms * 10**p50
    f = 10.0**-p50
    t = Table({"log_M_star": d[:, 0], "M_200c": M, "M_200c_lo": M - Ms * 10**p16, "M_200c_hi": Ms * 10**p84 - M,
               "f_star": f, "f_star_lo": f - 10.0**-p84, "f_star_hi": 10.0**-p16 - f,
               "f_star_2sig_lo": f - 10.0**-p97, "f_star_2sig_hi": 10.0**-p2 - f})
    _write_ecsv(t, _meta(
        id="Reyes2012",
        reference="Reyes, R., Mandelbaum, R., Gunn, J. E., et al. 2012, MNRAS, 425, 2610",
        bibcode="2012MNRAS.425.2610R", doi="10.1111/j.1365-2966.2012.21472.x", arxiv="1110.4107",
        url="https://arxiv.org/abs/1110.4107",
        description="Galaxy-galaxy lensing of SDSS disc galaxies: M200c/M_star as a function of stellar mass "
                    "(galaxy-scale haloes, 1e11-1e13 Msun).",
        data_type="binned", quantities=["f_star"], aperture=_ap("200c"),
        mass=_mass("M_200c", "200c", "weak_lensing", "stacked galaxy-galaxy lensing"),
        cosmology=dict(H0=70.0, Om0=0.3, Ode0=0.7, note="VERIFY"),
        stellar=dict(imf=None, includes_bcg=True, includes_icl=False),
        redshift=dict(min=0.0, max=0.2, effective=0.1),
        table_source="Table 5",
        status="needs_check",
        caveats=["LEGACY BUG FIXED: columns are percentiles (2.5,16,50,84,97.5) of log(M200c/M*); legacy "
                 "code used them as errors on M* and f_star",
                 "binned in stellar mass: mass and fraction errors are fully anti-correlated",
                 "single central galaxies: the cluster BCG+satellite aperture model does not apply",
                 "IMF to be recorded"],
    ), out)


def _pm(s):
    s = s.replace("−", "-")
    v, e = s.split("±")
    return float(v), float(e)


def sanderson2013(src, out):
    names, zs, vals = [], [], []
    for line in _rows(src / "Sanderson_2013.dat"):
        line = line.replace("−", "-")
        tok = line.split()
        names.append(tok[0])
        zs.append(float(tok[1]))
        parts = re.findall(r"(-?[\d.]+)\s*±\s*([\d.]+)", line)
        vals.append([float(x) for p in parts for x in p])
    a = np.array(vals)  # M2500,e, M500,e, r2500,e, r500,e, f2500,e, f500,e, alpha,e
    for delta, (im, iff) in (("500c", (2, 10)), ("2500c", (0, 8))):
        t = Table({"name": names, "z": zs, f"M_{delta}": a[:, im] * 1e14, f"M_{delta}_lo": a[:, im + 1] * 1e14,
                   f"M_{delta}_hi": a[:, im + 1] * 1e14,
                   "f_gas": a[:, iff], "f_gas_lo": a[:, iff + 1], "f_gas_hi": a[:, iff + 1]})
        _write_ecsv(t, _meta(
            id="Sanderson2013" if delta == "500c" else "Sanderson2013_2500c",
            reference="Sanderson, A. J. R., O'Sullivan, E., Ponman, T. J., et al. 2013, MNRAS, 429, 3288",
            bibcode="2013MNRAS.429.3288S", doi="10.1093/mnras/sts586", arxiv="1212.1613",
            url="https://arxiv.org/abs/1212.1613",
            description="Five low-mass clusters from the Gonzalez et al. sample with deep XMM-Newton data; "
                        "hydrostatic masses and gas fractions.",
            data_type="points", quantities=["f_gas"], aperture=_ap(delta),
            mass=_mass(f"M_{delta}", delta, "hydrostatic", "XMM-Newton HSE"),
            cosmology=dict(H0=70.0, Om0=0.3, Ode0=0.7, note="VERIFY"),
            redshift=dict(min=float(min(zs)), max=float(max(zs)), effective=None),
            table_source="Table 4; redshifts from Table 1 / Simbad (Douspis, Oct 2026)",
            overlaps_with=["Gonzalez2013"],
        ), out)

def sun2009(src, out):
    rows = [r.split() for r in _rows(src / "Sun_2008.dat")]
    names = [r[0] for r in rows]
    d = np.array([[float(x) for x in r[1:]] for r in rows])
    t = Table({"name": names, "z": d[:, 6],
               "M_500c": d[:, 0] * 1e13, "M_500c_hi": d[:, 1] * 1e13, "M_500c_lo": -d[:, 2] * 1e13,
               "f_gas": d[:, 3], "f_gas_hi": d[:, 4], "f_gas_lo": -d[:, 5]})
    _write_ecsv(t, _meta(
        id="Sun2009",
        reference="Sun, M., Voit, G. M., Donahue, M., et al. 2009, ApJ, 693, 1142",
        bibcode="2009ApJ...693.1142S", doi="10.1088/0004-637X/693/2/1142", arxiv="0805.2320",
        url="https://arxiv.org/abs/0805.2320",
        description="43 nearby galaxy groups with Chandra; hydrostatic masses (subsample with r500 coverage).",
        data_type="points", quantities=["f_gas"], aperture=_ap(),
        mass=_mass("M_500c", "500c", "hydrostatic", "Chandra HSE"),
        cosmology=dict(H0=73.0, Om0=0.24, Ode0=0.76),
        redshift=dict(min=float(d[:, 6].min()), max=float(d[:, 6].max()), effective=None),
        table_source="Table 3; names and redshifts restored by Douspis (Oct 2026)",
        caveats=["legacy id 'Sun 2008' (arXiv year); published 2009"],
    ), out)

def vikhlinin2006(src, out):
    rows = [r.split() for r in _rows(src / "Vikhlinin_2006.dat")]
    names = [r[0] for r in rows]
    d = np.array([[float(x) for x in r[1:]] for r in rows])  # 17 values + z
    z = d[:, 17]
    for delta, (iM, iF) in (("500c", (10, 14)), ("2500c", (8, 12))):
        t = Table({"name": names, "z": z, "R_500c_kpc": d[:, 0], "T_spec": d[:, 2], "c500": d[:, 6],
                   f"M_{delta}": d[:, iM] * 1e14, f"M_{delta}_lo": d[:, iM + 1] * 1e14,
                   f"M_{delta}_hi": d[:, iM + 1] * 1e14,
                   "f_gas": d[:, iF], "f_gas_lo": d[:, iF + 1], "f_gas_hi": d[:, iF + 1]})
        _write_ecsv(t, _meta(
            id="Vikhlinin2006" if delta == "500c" else "Vikhlinin2006_2500c",
            reference="Vikhlinin, A., Kravtsov, A., Forman, W., et al. 2006, ApJ, 640, 691",
            bibcode="2006ApJ...640..691V", doi="10.1086/500288", arxiv="astro-ph/0507092",
            url="https://iopscience.iop.org/article/10.1086/500288",
            description="13 relaxed low-z clusters with Chandra temperature profiles to ~r500; hydrostatic "
                        "masses (10 with r500 measurements).",
            data_type="points", quantities=["f_gas"], aperture=_ap(delta),
            mass=_mass(f"M_{delta}", delta, "hydrostatic", "Chandra HSE"),
            cosmology=dict(H0=71.0, Om0=0.3, Ode0=0.7),
            redshift=dict(min=float(z.min()), max=float(z.max()), effective=None),
            table_source="Table 2 of the journal version (differs from arXiv v1); names and z added by Douspis",
            changelog=["LEGACY BUG FIXED: both legacy codes used h=0.72 while the paper and the data file "
                       "state h=0.71 (1.4% in M, 2.1% in f_gas)"],
        ), out)

def wicker2023(src, out):
    d = np.loadtxt(src / "Wicker_2022.dat", skiprows=2)
    t = Table({"z": d[:, 0], "M_gas_500c": d[:, 1] * 1e14,
               "M_500c": d[:, 4] * 1e14, "M_500c_hi": d[:, 5] * 1e14, "M_500c_lo": d[:, 6] * 1e14,
               "f_gas": d[:, 7], "f_gas_hi": d[:, 8], "f_gas_lo": d[:, 9]})
    _write_ecsv(t, _meta(
        id="Wicker2023",
        reference="Wicker, R., Douspis, M., Salvati, L., Aghanim, N. 2023, A&A, 674, A48",
        bibcode="2023A&A...674A..48W", doi="10.1051/0004-6361/202243922", arxiv="2204.12823",
        url="https://arxiv.org/abs/2204.12823",
        description="Planck SZ-selected clusters with XMM-Newton hydrostatic masses and gas fractions, used "
                    "to constrain the hydrostatic mass bias.",
        data_type="points", quantities=["f_gas"], aperture=_ap(),
        mass=_mass("M_500c", "500c", "hydrostatic", "XMM-Newton HSE (literature compilation)"),
        cosmology=dict(H0=70.0, Om0=0.3, Ode0=0.7, note="VERIFY"),
        redshift=dict(min=float(d[:, 0].min()), max=float(d[:, 0].max()), effective=None),
        table_source="Table A.1",
        caveats=["legacy id 'Wicker 2022' (arXiv year); published 2023", "cluster names missing"],
    ), out)


def zhang2011(src, out):
    d = np.genfromtxt(src / "Zhang_2011.dat", dtype=None, names=True, encoding="utf-8", skip_header=1)
    t = Table({"name": d["name"], "z": d["z"],
               "M_500c": d["M500"] * 1e14, "M_500c_lo": d["M500_err"] * 1e14, "M_500c_hi": d["M500_err"] * 1e14,
               "M_500c_from_Mgas": d["M500_Mgas"] * 1e14,
               "f_gas": d["fgas"], "f_gas_lo": d["fgas_err"], "f_gas_hi": d["fgas_err"],
               "f_star": d["fstar"], "f_star_lo": d["fstar_err"], "f_star_hi": d["fstar_err"]})
    _write_ecsv(t, _meta(
        id="Zhang2011",
        reference="Zhang, Y.-Y., Laganá, T. F., Pierini, D., et al. 2011, A&A, 535, A78",
        bibcode="2011A&A...535A..78Z", doi="10.1051/0004-6361/201116803", arxiv="1109.0390",
        url="https://arxiv.org/abs/1109.0390",
        description="19 HIFLUGCS clusters with XMM-Newton hydrostatic masses and SDSS stellar masses.",
        data_type="points", quantities=["f_gas", "f_star"], aperture=_ap(),
        mass=_mass("M_500c", "500c", "hydrostatic", "XMM-Newton HSE"),
        cosmology=dict(H0=70.0, Om0=0.3, Ode0=0.7),
        stellar=dict(imf=None, includes_bcg=True, includes_icl=False),
        redshift=dict(min=float(d["z"].min()), max=float(d["z"].max()), effective=None),
        table_source="Tables 1-2",
        caveats=["LEGACY BUG FIXED: f_star was rescaled as h^-3/2 (gas scaling); stellar fractions scale "
                 "as h^-1", "IMF to be recorded"],
    ), out)


def khalil2026(src, out):
    d = np.genfromtxt(src / "Khalil_eRositaGroups_2608.17735.dat", delimiter=",", names=True,
                      dtype=None, encoding="utf-8")
    for delta in ("500c", "200c"):
        k = delta[:3]
        M, dM = d[f"M{k}c_1e13_Msun"] * 1e13, d[f"M{k}c_err_1e13_Msun"] * 1e13
        Mg, dMg = d[f"Mgas_{k}_1e12_Msun"] * 1e12, d[f"Mgas_{k}_err_1e12_Msun"] * 1e12
        ok = np.isfinite(M) & np.isfinite(Mg) & np.isfinite(dMg)
        fg = Mg / M
        efg = _ratio_err(Mg, dMg, M, dM)
        t = Table({"name": [str(x) for x in d["Group_ID"][ok]], "z": d["z"][ok],
                   f"M_{delta}": M[ok], f"M_{delta}_lo": dM[ok], f"M_{delta}_hi": dM[ok],
                   f"M_gas_{delta}": Mg[ok], "sigma_v": d["sigma_v_km_s"][ok], "N_gal": d["Ngal"][ok],
                   "f_gas": fg[ok], "f_gas_lo": efg[ok], "f_gas_hi": efg[ok]})
        _write_ecsv(t, _meta(
            id="Khalil2026" if delta == "500c" else "Khalil2026_200c",
            reference="Khalil, H., Finoguenov, A., Eckert, D., et al. 2026, arXiv:2608.17735",
            bibcode=None, doi="10.48550/arXiv.2608.17735", arxiv="2608.17735",
            url="https://arxiv.org/abs/2608.17735",
            description="25 nearby galaxy groups from the 2MRS optical catalogue with eROSITA (eRASS1) surface "
                        "brightness profiles; hot gas mass and gas fraction out to R200.",
            data_type="points", quantities=["f_gas"], aperture=_ap(delta),
            mass=_mass(f"M_{delta}", delta, "dynamical",
                       "VERIFY: catalogue lists sigma_v and N_gal; mass method not stated in the abstract"),
            cosmology=dict(H0=70.0, Om0=0.3, Ode0=0.7, note="VERIFY (h=0.7 assumed in the Douspis loader)"),
            redshift=dict(min=float(np.nanmin(d["z"][ok])), max=float(np.nanmax(d["z"][ok])), effective=None),
            selection="optically selected (2MRS) groups detected by eROSITA; complete sample",
            table_source="catalogue table (transcribed by M. Douspis, Oct 2026)",
            status="needs_check",
            contributors=["M. Douspis"],
            caveats=["preprint: check against the published version",
                     "f_gas errors propagated assuming independent M and M_gas errors",
                     f"{int((~ok).sum())} groups without a gas-mass measurement at r{delta} dropped"],
        ), out)


def akino2022(out):
    def comp(q, a, a_hi, a_lo, b, b_hi, b_lo, label):
        A = np.exp(a) * 1e12 / 1e14
        return dict(quantity=q, label=label, norm=float(A),
                    norm_err_lo=float(A - np.exp(a - a_lo) * 1e-2), norm_err_hi=float(np.exp(a + a_hi) * 1e-2 - A),
                    slope=float(b - 1), slope_err_lo=float(b_lo), slope_err_hi=float(b_hi),
                    pivot_mass=1e14, mass_range=[1e13, 1e15], corr=0.0)
    comps = [comp("f_gas", 2.02, 0.08, 0.10, 1.28, 0.12, 0.09, "gas"),
             comp("f_star", 0.68, 0.15, 0.07, 0.76, 0.12, 0.08, "stars"),
             comp("f_bcg", -1.07, 0.18, 0.08, 0.39, 0.13, 0.08, "BCG")]
    _write_yaml(comps, _meta(
        id="Akino2022",
        reference="Akino, D., Eckert, D., Okabe, N., et al. 2022, PASJ, 74, 175",
        bibcode="2022PASJ...74..175A", doi="10.1093/pasj/psab115", arxiv="2111.10080",
        url="https://arxiv.org/abs/2111.10080",
        description="HSC-XXL: 136 X-ray selected groups/clusters with HSC weak-lensing masses; regressions "
                    "ln(Y/1e12 Msun) = alpha + beta ln(M500/1e14 Msun) for gas, stellar and BCG mass, "
                    "re-expressed here as fractions f = norm (M/1e14)^(beta-1).",
        data_type="relation", quantities=["f_gas", "f_star", "f_bcg"], aperture=_ap(),
        mass=_mass("M_500c", "500c", "weak_lensing", "HSC WL, joint population model"),
        cosmology=dict(H0=70.0, Om0=0.3, Ode0=0.7, note="VERIFY"),
        stellar=dict(imf=None, includes_bcg=True, includes_icl=False),
        redshift=dict(min=0.0, max=1.0, effective=0.3),
        table_source="Appendix E",
        status="needs_check",
        caveats=["asymmetric parameter errors treated as independent (covariance not available)",
                 "redshift-evolution term, if any, ignored", "IMF to be recorded"],
    ), out)


def popesso2024(out):
    for delta, (A0, dA0, B, dB, rng) in {
        "500c": (2.23e-7, 0.18e-7, 0.39, 0.02, [1e12, 1e15]),
        "200c": (2.09e-6, 0.14e-6, 0.33, 0.02, [1.5e12, 1.5e15]),
    }.items():
        piv = 1e14
        A = A0 * piv**B
        comps = [dict(quantity="f_gas", label=f"r{delta}", norm=float(A), norm_err_lo=float(A * dA0 / A0),
                      norm_err_hi=float(A * dA0 / A0), slope=B, slope_err_lo=dB, slope_err_hi=dB,
                      pivot_mass=piv, mass_range=rng, corr=0.0)]
        _write_yaml(comps, _meta(
            id="Popesso2024" if delta == "500c" else "Popesso2024_200c",
            reference="Popesso, P., Biviano, A., Marini, I., et al., A&A, doi:10.1051/0004-6361/202453256 (arXiv:2411.16555, 2024)",
            bibcode=None, doi="10.1051/0004-6361/202453256", arxiv="2411.16555", url="https://arxiv.org/abs/2411.16555",
            description="eROSITA (eFEDS) stacking of optically selected GAMA groups + eFEDS and CHEX-MATE "
                        "clusters: f_gas = A (M/Msun)^B, re-expressed at a 1e14 Msun pivot.",
            data_type="relation", quantities=["f_gas"], aperture=_ap(delta),
            mass=_mass(f"M_{delta}", delta, "mixed",
                       "GAMA group luminosity (groups), LX-M (eFEDS), CHEX-MATE HSE; NFW c-M (Dutton & Maccio 14)"),
            cosmology=dict(H0=67.74, Om0=0.3089, Ode0=0.6911),
            redshift=dict(min=0.0, max=0.2, effective=0.1),
            table_source="Eqs. 4-5",
            status="needs_check",
            caveats=["LEGACY BUG FIXED: published fit has no pivot (M in Msun); legacy band combined "
                     "(A±dA)(M)^(B±dB) at M=1 Msun, inflating the band by ~2x at 1e15. Here the "
                     "fractional normalisation error is ASSUMED to apply at the 1e14 pivot -> ask authors "
                     "for the covariance",
                     "mass range approximate; redshift range to verify"],
            changelog=["cosmology and relations checked against arXiv:2411.16555v2 2026-10-31"],
        ), out)


def main(src, dst):
    src, dst = Path(src), Path(dst)
    obs = dst / "observations"
    obs.mkdir(parents=True, exist_ok=True)
    for fn in (arnaud2007, bulbul2024, chiu2018, ettori2010, giodini2009, gonzalez2013, grandis2024,
               kravtsov2018, landry2013, mantz2022, mulroy2019, reyes2012, sanderson2013, sun2009,
               vikhlinin2006, wicker2023, zhang2011, khalil2026):
        fn(src, obs)
    akino2022(obs)
    popesso2024(obs)


if __name__ == "__main__":
    main(*sys.argv[1:3])

import numpy as np
import pytest

import barycon as bc
from barycon import conversions as cv
from barycon.dataset import Dataset, Relation


@pytest.fixture(scope="module")
def col():
    return bc.load()


# --------------------------------------------------------------------- data
def test_all_entries_valid_and_unique(col):
    ids = [e.id for e in col]
    assert len(ids) == len(set(ids))
    for e in col:
        bc.validate_meta(e.meta)
        assert e.path.stem == e.id, "file name must equal the entry id"


def test_every_entry_has_required_columns(col):
    for e in col:
        if isinstance(e, Dataset):
            mc = e.meta["mass"]["column"]
            assert mc in e.table.colnames, (e.id, mc)
            for q in e.quantities:
                assert q in e.table.colnames, (e.id, q)
                f = np.asarray(e.table[q], float)
                assert np.nanmax(f) < 1.0 and np.nanmin(f) > 0, (e.id, q)
            M = np.asarray(e.table[mc], float)
            assert 1e9 < np.nanmedian(M) < 1e16, (e.id, "mass units?")
        else:
            for c in e.components:
                assert c["quantity"] in e.quantities


def test_bibtex_for_every_entry(col):
    bib = col.bibtex()
    for e in col:
        assert e.id.split("_")[0] in bib, f"missing BibTeX for {e.id}"


# -------------------------------------------------------------- conversions
def test_pure_h_change_is_exact():
    f = cv.cosmology_factors(0.3, {"H0": 70, "Om0": 0.3}, {"H0": 67.66, "Om0": 0.3})
    r = 70 / 67.66
    assert np.allclose(f["s"], 1.0)
    assert np.allclose(f["mass"], r)          # M ∝ h^-1
    assert np.allclose(f["f_gas"], r**1.5)    # f_gas ∝ h^-3/2
    assert np.allclose(f["f_star"], r)        # f_star ∝ h^-1


def test_omega_m_change_small_at_low_z():
    f = cv.cosmology_factors(0.05, {"H0": 70, "Om0": 0.3}, {"H0": 70, "Om0": 0.25})
    assert abs(f["mass"][0] - 1) < 0.01 and abs(f["f_gas"][0] - 1) < 0.01


def test_nfw_round_trip_and_known_ratio():
    M = np.array([1e13, 1e14, 1e15])
    m500 = cv.convert_mass_definition(M, 0.0, "200c", "500c", "Planck18")
    back = cv.convert_mass_definition(m500, 0.0, "500c", "200c", "Planck18")
    assert np.allclose(back, M, rtol=1e-4)
    # c200=5 -> c500≈3.31, M500/M200 = 2.5 (c500/c200)^3 ≈ 0.722
    r = cv.convert_mass_definition(1e14, 0.0, "200c", "500c", "Planck18", c=5.0)[0] / 1e14
    assert abs(r - 0.722) < 0.003


def test_hse_bias_identity():
    f = cv.hse_bias_factors(0.0)
    assert np.isclose(f["mass"], 1) and np.isclose(f["f_gas"], 1)


def test_imf_factor():
    assert np.isclose(cv.imf_factor("salpeter", "chabrier"), 0.63)
    assert np.isclose(cv.imf_factor("chabrier", "chabrier"), 1.0)


# ---------------------------------------------------------------- datasets
def test_cosmology_round_trip(col):
    d = col["Vikhlinin2006"]
    back = d.to_cosmology("Planck18").to_cosmology(d.meta["cosmology"])
    assert np.allclose(back.table["f_gas"], d.table["f_gas"], rtol=1e-6)
    assert np.allclose(back.table["M_500c"], d.table["M_500c"], rtol=1e-6)
    assert len(d.to_cosmology("Planck18").history) == 1


def test_aperture_conversion_requires_explicit_model(col):
    d = col["Arnaud2007"]
    with pytest.raises(ValueError):
        d.to_delta("200c")
    m = d.to_delta("200c", masses_only=True)
    assert np.all(m.table["M_200c"] > d.table["M_500c"])
    g = d.to_delta("200c", eta_gas=0.3)
    assert np.all(g.table["f_gas"] > d.table["f_gas"])


def test_shell_cannot_be_converted(col):
    with pytest.raises(ValueError):
        col["Mantz2022"].to_delta("500c", eta_gas=0.3)


def test_filter_excludes_derived(col):
    gas = col.filter(quantity="f_gas", delta="500c")
    assert "Giodini2009" not in gas and "Vikhlinin2006" in gas
    assert "Giodini2009" in col.filter(quantity="f_gas", exclude_derived=False)


def test_to_table(col):
    t = col.filter(quantity="f_gas", delta="500c", data_type="points").to_table("f_gas", h_units=True)
    assert {"M", "f", "id"} <= set(t.colnames) and len(t) > 100


# --------------------------------------------------------------- relations
def test_popesso_matches_published_formula(col):
    r = col["Popesso2024"]
    M = np.array([1e13, 3e14])
    assert np.allclose(r.evaluate(M, "f_gas"), 2.23e-7 * M**0.39)
    lo, mid, hi = r.evaluate(np.array([1e15]), "f_gas", n_samples=4000)
    assert (hi - lo)[0] / mid[0] < 0.3  # legacy band was ~ x2


def test_relation_cosmology_shift(col):
    r = col["Akino2022"]
    r2 = r.to_cosmology({"H0": 67.66, "Om0": 0.3})  # pure-h change: exact scaling expected
    M = np.array([1e14])
    # f_gas at the *same physical* object must scale as (h_old/h_new)^1.5
    f_old = r.evaluate(M, "f_gas")
    f_new = r2.evaluate(M * 0.7 / 0.6766, "f_gas")
    assert np.allclose(f_new / f_old, (0.7 / 0.6766) ** 1.5, rtol=2e-3)


# ---------------------------------------------- regression tests (legacy bugs)
def test_kravtsov_uses_true_m500(col):
    t = col["Kravtsov2018"].table
    a2142 = t[t["name"] == "A2142"][0]
    assert np.isclose(a2142["M_500c"], 11.96e14)


def test_ettori_total_mass_from_r500(col):
    t = col["Ettori2010"].table
    assert np.all(np.isfinite(t["z"])) and len(t) == 44


def test_vikhlinin_h071(col):
    assert col["Vikhlinin2006"].meta["cosmology"]["H0"] == 71.0


def test_bulbul_units(col):
    M = np.asarray(col["Bulbul2024"].table["M_500c"])
    assert 5e13 < np.median(M) < 5e14


# ----------------------------------------------------------------- fitting
def test_power_law_recovers_truth():
    from astropy.table import Table
    rng = np.random.default_rng(3)
    M = 10 ** rng.uniform(13.5, 15, 200)
    f = 0.12 * (M / 3e14) ** 0.25
    err = 0.05 * f
    t = Table({"M_500c": M, "f_gas": f * (1 + 0.05 * rng.standard_normal(M.size)), "f_gas_lo": err, "f_gas_hi": err})
    meta = dict(bc.load()["Arnaud2007"].meta, id="Mock2026")
    d = Dataset(t, meta)
    r = bc.fit(d, "power_law", pivot=3e14, h_units=False)
    assert abs(r.params["A"] - 0.12) < 3 * r.errors["A"] + 1e-3
    assert abs(r.params["beta"] - 0.25) < 3 * r.errors["beta"] + 1e-3
    assert 0.6 < r.chi2 / r.dof < 1.5


def test_fit_aliases_and_old_signature(col):
    sel = col.filter(ids=["Wicker2023", "Arnaud2007"]).homogenize("Planck18")
    r = sel.fit("gas_frac_M_simple", p0=[1e14, 0.3, 0.1])
    assert r.model == "power_law" and r.pivot == 1e14
    r = sel.fit("gas_frac_M_double")
    assert r.model == "saturating" and r.errors["M_t"] > 0.01 * r.params["M_t"]
    r = sel.fit("gas_frac_M_z")
    assert len(r.params) == 5


def test_fit_with_scatter_has_sane_chi2(col):
    sel = col.filter(ids=["Wicker2023", "Arnaud2007", "Sun2009", "Vikhlinin2006"]).homogenize("Planck18")
    r = sel.fit("power_law", scatter=True)
    assert 0.7 < r.chi2 / r.dof < 1.3 and 0.05 < r.scatter < 0.4


def test_fit_refuses_mixed_apertures(col):
    with pytest.raises(ValueError):
        col.filter(ids=["Vikhlinin2006", "Vikhlinin2006_2500c"]).fit()
    with pytest.raises(ValueError):
        col.filter(ids=["Mantz2022"]).fit()


# ------------------------------------------------------- compatibility API
def test_baryonfractions_compat():
    bf = bc.BaryonFractions()
    for k in ("arnaud_2007", "wicker_2022", "sun_2008", "vikhlinin_2006", "ettori_2011", "khalil_2026"):
        d = bf.gas_data[k]
        assert {"M500", "f_gas", "f_gas_err", "reference"} <= set(d)
        assert np.all(d["M500"] > 1e12) and np.all(d["f_gas"] < 0.3)
    assert "M200" in bf.stellar_data["reyes_2012"]
    assert {"M500_ih", "f_gas_m", "f_gas_p"} <= set(bf.gas_fits["popesso_2024"])
    popt, pcov = bf.fit(["wicker_2022", "arnaud_2007"], fit_type="gas_frac_M_double")
    assert popt.shape == (2,) and pcov.shape == (2, 2)
    # Vikhlinin: corrected h (0.71) and Planck18 rescaling, f_gas ∝ h^-3/2
    v = bc.load()["Vikhlinin2006"]
    expect = np.asarray(v.table["f_gas"])[0] * (0.71 / bf.hubble) ** 1.5
    assert np.isclose(bf.gas_data["vikhlinin_2006"]["f_gas"][0], expect, rtol=5e-3)


def test_khalil_entries(col):
    for eid in ("Khalil2026", "Khalil2026_200c"):
        t = col[eid].table
        assert np.all(np.isfinite(t["z"])) and np.median(t["f_gas"]) < 0.1
    assert col["Khalil2026_200c"].delta == "200c"


def test_grandis_points_with_local_slopes():
    """Grandis+24 is stored as one point per survey at its pivot, with the local slope."""
    from barycon import Dataset
    e = bc.load()["Grandis2024"]
    assert isinstance(e, Dataset)
    d = e.get("f_star")
    assert len(d["M"]) == 2 and np.allclose(d["slope"], -0.2)
    assert set(e.get("f_gas")["proxy"]) == {"sz_scaling_relation", "xray_scaling_relation", "weak_lensing"}
    moved = e.to_delta("200c", eta_gas=0.3, f_bcg=lambda M: 0.5 * np.ones_like(M))
    assert not any(c.startswith("slope_") for c in moved.table.colnames)


def test_plot_collection_runs():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from barycon.plotting import dataset_color, short_label
    col = bc.load().filter(quantity="f_gas", delta="500c", region="sphere")
    fig, ax = plt.subplots(layout="constrained")
    col.plot("f_gas", ax=ax, universal=0.157)
    assert short_label(bc.load()["Landry2013_2500c"]) == "Landry+13"
    assert dataset_color(bc.load()["Landry2013_2500c"]) == dataset_color(bc.load()["Landry2013"])
    plt.close(fig)

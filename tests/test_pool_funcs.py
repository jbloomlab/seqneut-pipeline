"""Tests for ``scripts/pool_funcs.py``."""

import math

import pandas as pd
import pytest
from pool_funcs import (
    DEFAULT_LINEARITY_SLOPE_TOL,
    DEFAULT_LINEARITY_WINDOW_DILUTIONS,
    IDEAL_SLOPE,
    analysis_wells,
    barcode_fractions,
    diverging_stocks,
    fit_dilution_windows,
    log_log_slope,
    member_barcodes,
    pool_samples,
    pool_units,
    unit_shares,
)


def series(ratios, replicates=2):
    """Wells over a dilution series, `ratios` mapping each dilution to its viral:standard.

    The fit is on the neutralization-standard odds, the reciprocal of that ratio, which is
    what the report plots; a series is described by the ratio here because that is the way
    round the experiment is set up.

    Every replicate of a dilution is given a slightly different ratio, so that a fit over
    them has a standard error rather than the exactly zero one that identical replicates
    would give, as real wells never agree exactly.

    """
    rows = []
    for dilution, ratio in ratios.items():
        for replicate in range(replicates):
            rows.append(
                {
                    "plate": "p1",
                    "well": f"{'ABCDEFGH'[replicate]}{int(math.log2(dilution))}",
                    "dilution_factor": dilution,
                    "neut_standard_odds": 1 / (ratio * (1 + 0.02 * replicate)),
                }
            )
    return pd.DataFrame(rows)


def ideal_series(dilutions, titer=1000.0):
    """A series responding exactly to dilution, so every window of it is strictly linear."""
    return series({d: titer / d for d in dilutions})


def saturated_then_ideal(saturated, linear, titer=1000.0):
    """A series flat over `saturated` and then responding over `linear`.

    The flat part is what a pool too concentrated to respond to dilution looks like: its
    counts barely fall as it is diluted.

    """
    ratios = {d: titer / linear[0] for d in saturated}
    ratios.update({d: titer / d for d in linear})
    return series(dict(sorted(ratios.items())))


def library(barcodes):
    """A viral library, `barcodes` mapping each barcode to its (strain, stock)."""
    return pd.DataFrame(
        [
            {"barcode": barcode, "strain": strain, "stock": stock}
            for barcode, (strain, stock) in barcodes.items()
        ]
    )


# --- fitting the dilution series -------------------------------------------------------


def test_fit_dilution_windows_over_an_ideal_series():
    """Every window of a series responding exactly to dilution is strictly linear."""
    windows = fit_dilution_windows(ideal_series([2, 4, 8, 16, 32]))
    assert [w["window"] for w in windows] == [[2, 4, 8], [4, 8, 16], [8, 16, 32]]
    # Python numbers, so that a window written to a CSV reads as `128` and not as
    # `np.int64(128)`
    assert all(type(d) is int for w in windows for d in w["window"])
    for window in windows:
        assert window["slope"] == pytest.approx(IDEAL_SLOPE, abs=0.01)
        assert window["within_tol"]


def test_fit_dilution_windows_takes_the_window_length():
    """`linearity_window_dilutions` sets how many dilutions a window spans."""
    windows = fit_dilution_windows(ideal_series([2, 4, 8, 16, 32]), window_dilutions=5)
    assert [w["window"] for w in windows] == [[2, 4, 8, 16, 32]]


def test_fit_dilution_windows_marks_a_saturated_run_nonlinear():
    """A window over the flat end of a series is not within tolerance of the ideal slope."""
    windows = fit_dilution_windows(saturated_then_ideal([2, 4, 8], [16, 32, 64]))
    flat = next(w for w in windows if w["window"] == [2, 4, 8])
    assert flat["slope"] == pytest.approx(0, abs=0.05)
    assert not flat["within_tol"]


def test_fit_dilution_windows_rejects_too_few_dilutions():
    """A series shorter than one window has no window to fit."""
    assert fit_dilution_windows(ideal_series([2, 4])) == []


@pytest.mark.parametrize("window_dilutions", [1, 2, 4, 0, -1])
def test_fit_dilution_windows_rejects_an_unusable_window(window_dilutions):
    """A window needs at least three dilutions to fit a slope, and an odd number for a middle."""
    with pytest.raises(ValueError, match="linearity_window_dilutions"):
        fit_dilution_windows(
            ideal_series([2, 4, 8, 16, 32]), window_dilutions=window_dilutions
        )


# --- choosing the analysis wells ---------------------------------------------------------


def test_log_log_slope_fits_a_line_through_its_points():
    """The slope and intercept are those of the line the points lie on."""
    # y = 3 * x**2 on log-log axes: slope 2, intercept log10(3)
    xs = [1, 10, 100]
    slope, slope_se, intercept = log_log_slope(xs, [3 * x**2 for x in xs])
    assert slope == pytest.approx(2)
    assert intercept == pytest.approx(math.log10(3))
    # the points are exactly on the line, so the fit has no residual error
    assert slope_se == pytest.approx(0)


def test_fit_dilution_windows_places_each_window_s_own_line():
    """A window's slope and intercept put its fitted line through its own wells.

    The report draws each window's fit from these, so a window whose line is placed wrong
    would be drawn through nothing.

    """
    wells = series({2: 500.0, 4: 250.0, 8: 125.0, 16: 62.5})
    for window in fit_dilution_windows(wells):
        in_window = wells[wells["dilution_factor"].isin(window["window"])]
        predicted = [
            10 ** (window["intercept"] + window["slope"] * math.log10(d))
            for d in in_window["dilution_factor"]
        ]
        # the replicates of a dilution straddle the line, so the fit runs through their
        # mean rather than through either of them
        assert sum(
            math.log10(o) - math.log10(p)
            for (o, p) in zip(in_window["neut_standard_odds"], predicted, strict=True)
        ) == pytest.approx(0, abs=1e-9)


def test_fit_dilution_windows_closes_over_a_missing_dilution():
    """A window is consecutive among the dilutions given, not among the plate's.

    A dilution whose every well failed QC does not reach the fit, and the windows then
    close over the gap: a window can cover a wider span of the series than its length
    suggests. `analyze_pool.py` warns where this happens, as the fit cannot tell a
    dilution that was never measured from one that was dropped.

    """
    # 8 is absent, as it would be with both of its wells failing QC
    gapped = series({2: 500.0, 4: 250.0, 16: 62.5, 32: 31.25})
    windows = fit_dilution_windows(gapped)
    assert [w["window"] for w in windows] == [[2, 4, 16], [4, 16, 32]]
    # the first window spans 2- to 16-fold, three dilutions covering a threefold step
    assert windows[0]["window"][-1] / windows[0]["window"][0] == 8


def test_analysis_wells_takes_the_middle_of_the_first_strict_run():
    """The composition is read from the middle dilution of the first strictly linear window."""
    wells = series_wells = saturated_then_ideal([2, 4, 8], [16, 32, 64, 128])
    chosen = analysis_wells(fit_dilution_windows(wells), series_wells)
    assert chosen.dilution == 32
    assert chosen.titer_window == [16, 32, 64]
    assert chosen.strict
    assert not chosen.wells_given
    assert sorted(chosen.wells) == sorted(
        ("p1", w)
        for w in series_wells.loc[series_wells["dilution_factor"] == 32, "well"]
    )


def test_analysis_wells_chooses_nothing_where_only_the_noise_would_admit_a_run():
    """A pool whose every window misses the tolerance is not read at all.

    Modeled on the most saturated window of a real equal-volume pool, whose slope of
    -0.416 on a standard error of 0.106 is far off the ideal but not by much more than its
    own noise. Such a pool does not clearly respond to being diluted, so no wells are
    chosen and the report says so, unless `linear_range_wells` names the wells.

    """
    noisy = pd.DataFrame(
        [
            {
                "plate": "p1",
                "well": w,
                "dilution_factor": d,
                "neut_standard_odds": 1 / r,
            }
            for w, d, r in [
                ("A1", 2, 500.0),
                ("B1", 2, 300.0),
                ("A2", 4, 460.0),
                ("B2", 4, 250.0),
                ("A3", 8, 330.0),
                ("B3", 8, 230.0),
            ]
        ]
    )
    windows = fit_dilution_windows(noisy)
    assert not any(w["within_tol"] for w in windows)
    chosen = analysis_wells(windows, noisy)
    assert chosen.wells == []
    assert chosen.dilution is None
    # the windows are still fit and returned, the report tabulating them with its warning
    assert [w["window"] for w in windows] == [[2, 4, 8]]


def test_analysis_wells_passes_over_a_saturated_run_for_a_strict_one():
    """The first window within tolerance wins, however far out in the series it lies.

    The saturated end of a healthy series misses the tolerance however noisy it is, so the
    range is placed where the pool actually responds rather than at its most saturated.

    """
    saturated_then_linear = pd.DataFrame(
        [
            {
                "plate": "p1",
                "well": w,
                "dilution_factor": d,
                "neut_standard_odds": 1 / r,
            }
            for w, d, r in [
                # a saturated end, scattered enough that its slope is poorly determined
                ("A1", 2, 500.0),
                ("B1", 2, 300.0),
                ("A2", 4, 460.0),
                ("B2", 4, 250.0),
                ("A3", 8, 330.0),
                ("B3", 8, 230.0),
                # responding to dilution from here, and tightly enough to be strict
                ("A4", 16, 165.0),
                ("B4", 16, 161.7),
                ("A5", 32, 82.5),
                ("B5", 32, 80.85),
                ("A6", 64, 41.25),
                ("B6", 64, 40.43),
            ]
        ]
    )
    windows = fit_dilution_windows(saturated_then_linear)
    assert [w["within_tol"] for w in windows] == [False, False, True, True]
    chosen = analysis_wells(windows, saturated_then_linear)
    assert chosen.strict
    # the middle of the first strict window, not of the saturated end at 4-fold
    assert chosen.dilution == 16


def test_analysis_wells_chooses_nothing_when_no_run_qualifies():
    """A pool that never responds to dilution is read from no well at all."""
    flat = series({2: 100.0, 4: 100.0, 8: 100.0, 16: 100.0})
    chosen = analysis_wells(fit_dilution_windows(flat), flat)
    assert chosen.wells == []
    assert chosen.dilution is None
    assert chosen.titer_window is None
    assert not chosen.strict
    assert not chosen.wells_given


def test_analysis_wells_chooses_nothing_when_the_series_is_too_short():
    """Too few dilutions to fit a window, with no wells given, chooses no wells."""
    short = ideal_series([2, 4])
    chosen = analysis_wells(fit_dilution_windows(short), short)
    assert chosen.wells == []
    assert chosen.dilution is None


def test_analysis_wells_given_takes_the_run_they_center():
    """Given wells, the titer comes from the window whose middle dilution is theirs.

    A dilution lies in up to three overlapping windows, whose titers differ, so the window is
    named rather than left to whichever is found first.

    """
    wells = ideal_series([2, 4, 8, 16, 32])
    chosen = analysis_wells(
        fit_dilution_windows(wells), wells, given=["p1:A3", "p1:B3"]
    )
    assert chosen.dilution == 8
    assert chosen.titer_window == [4, 8, 16]
    assert chosen.wells_given
    assert sorted(chosen.wells) == [("p1", "A3"), ("p1", "B3")]


def test_analysis_wells_given_at_an_uncentered_dilution_has_no_titer():
    """Wells no strict window centers are still read, but carry no titer."""
    wells = ideal_series([2, 4, 8, 16, 32])
    chosen = analysis_wells(
        fit_dilution_windows(wells), wells, given=["p1:A1", "p1:B1"]
    )
    assert chosen.dilution == 2
    assert chosen.titer_window is None
    assert sorted(chosen.wells) == [("p1", "A1"), ("p1", "B1")]


def test_analysis_wells_given_without_a_strict_run_has_no_titer():
    """Only a strict window gives a titer, on the explicit path as on the automatic one."""
    noisy = series({2: 500.0, 4: 300.0, 8: 250.0}, replicates=3)
    windows = fit_dilution_windows(noisy)
    assert not any(w["within_tol"] for w in windows)
    chosen = analysis_wells(
        windows,
        noisy,
        given=[f"p1:{w}" for w in noisy["well"][noisy["dilution_factor"] == 4]],
    )
    assert chosen.titer_window is None


def test_analysis_wells_given_with_too_short_a_series():
    """Wells given on a series too short to fit are read, without a titer."""
    short = ideal_series([2, 4])
    chosen = analysis_wells(
        fit_dilution_windows(short), short, given=["p1:A1", "p1:B1"]
    )
    assert chosen.dilution == 2
    assert chosen.titer_window is None


@pytest.mark.parametrize(
    "given, match",
    [
        (["p1:A3", "p1:A4"], "one dilution"),
        (["p1:A3", "p1:A3"], "more than once"),
        ([], "at least one well"),
        (["p1:Z9"], "not wells of"),
        (["A3"], "named as 'plate:well'"),
        (["other:A3"], "not wells of"),
    ],
)
def test_analysis_wells_rejects_bad_given_wells(given, match):
    """Bad given wells raise: two dilutions, duplicated, none, unknown, or unqualified."""
    wells = ideal_series([2, 4, 8, 16, 32])
    with pytest.raises(ValueError, match=match):
        analysis_wells(fit_dilution_windows(wells), wells, given=given)


def test_analysis_wells_rejects_given_wells_that_all_fail_qc():
    """Wells that failed QC are not in the frame fit, so none of them can be read."""
    wells = ideal_series([2, 4, 8, 16, 32])
    passing = wells[wells["dilution_factor"] != 8]
    with pytest.raises(ValueError, match="not wells of"):
        analysis_wells(fit_dilution_windows(passing), passing, given=["p1:A3", "p1:B3"])


def replicate_plates(dilutions, plates=("p1", "p2"), titer=1000.0):
    """One ideal dilution series measured on each of `plates`, as replicates of it."""
    return pd.concat(
        [ideal_series(dilutions, titer=titer).assign(plate=plate) for plate in plates],
        ignore_index=True,
    )


def test_analysis_wells_spans_the_plates_of_a_pool():
    """The chosen dilution's wells are taken from every plate, the plates being replicates."""
    wells = replicate_plates([2, 4, 8, 16, 32])
    chosen = analysis_wells(fit_dilution_windows(wells), wells)
    # every window of an ideal series is strict, so the first one is taken, as on one plate
    assert chosen.dilution == 4
    # the same well of each plate, which only the pair distinguishes
    assert sorted(chosen.wells) == [
        ("p1", "A2"),
        ("p1", "B2"),
        ("p2", "A2"),
        ("p2", "B2"),
    ]


def test_analysis_wells_given_picks_one_plate_of_a_pool():
    """A given well names its plate, so the composition can be read from one of them."""
    wells = replicate_plates([2, 4, 8, 16, 32])
    chosen = analysis_wells(fit_dilution_windows(wells), wells, given=["p2:A3"])
    assert chosen.dilution == 8
    assert chosen.wells == [("p2", "A3")]
    # the titer still comes from the window the dilution centers, fit over both plates
    assert chosen.titer_window == [4, 8, 16]
    assert chosen.wells_given


def test_analysis_wells_tells_the_plates_of_a_well_apart():
    """A well of one plate is not taken for the same well of another."""
    wells = replicate_plates([2, 4, 8, 16, 32])
    with pytest.raises(ValueError, match="not wells of"):
        analysis_wells(fit_dilution_windows(wells), wells, given=["p3:A3"])


def test_fit_dilution_windows_pools_the_replicate_plates():
    """A plate's wells add points to the fit at their dilution rather than a series of their own."""
    one = fit_dilution_windows(ideal_series([2, 4, 8, 16, 32]))
    both = fit_dilution_windows(replicate_plates([2, 4, 8, 16, 32]))
    # the same windows, as the dilutions are the same, and each fit over twice the wells
    assert [w["window"] for w in one] == [w["window"] for w in both]
    assert all(w["within_tol"] for w in both)


# --- which wells of the plate are the pool's ---------------------------------------------


def plate_samples():
    """A plate holding two pools, each diluted down its own series."""
    return pd.DataFrame(
        {
            "well": ["A1", "B1", "A2", "B2"],
            "sample": ["pool_x", "pool_x", "pool_y", "pool_y"],
            "dilution_factor": [2, 4, 2, 4],
        }
    )


def test_pool_samples_defaults_to_the_whole_plate():
    """No `pool_wells` is the case of one pool per plate."""
    samples = plate_samples()
    assert pool_samples(samples, {}).equals(samples)


def test_pool_samples_selects_one_pool_of_several():
    """A plate of several pools is read one pool at a time."""
    chosen = pool_samples(plate_samples(), {"sample": "pool_y"})
    assert chosen["well"].tolist() == ["A2", "B2"]


def test_pool_samples_compares_as_strings():
    """A numeric column read as a number still matches the value given in the config."""
    samples = plate_samples().assign(subpool=[1, 1, 2, 2])
    assert pool_samples(samples, {"subpool": 2})["well"].tolist() == ["A2", "B2"]


def test_pool_samples_requires_every_selection_to_hold():
    """Several keys select the wells matching all of them."""
    samples = plate_samples()
    chosen = pool_samples(samples, {"sample": "pool_x", "dilution_factor": 4})
    assert chosen["well"].tolist() == ["B1"]


def test_pool_samples_rejects_a_column_the_plate_lacks():
    """Selecting on a column that is not there is a configuration error, not no wells."""
    with pytest.raises(ValueError, match="not a column"):
        pool_samples(plate_samples(), {"subpool": "a"})


def test_pool_samples_rejects_a_selection_matching_no_well():
    """A value that matches nothing would otherwise read as a pool with no wells."""
    with pytest.raises(ValueError, match="matches no well"):
        pool_samples(plate_samples(), {"sample": "pool_z"})


# --- members of the pool -----------------------------------------------------------------


def test_member_barcodes_defaults_to_the_whole_library():
    """No `pool_members` means the pool was made from everything in the library."""
    lib = library({"AA": ("s1", "s1_a"), "CC": ("s2", "s2_a")})
    assert sorted(member_barcodes(lib, {})) == ["AA", "CC"]


def test_member_barcodes_selects_on_a_column():
    """`pool_members` selects the strains the pool was made from by a library column."""
    lib = library({"AA": ("s1", "s1_a"), "CC": ("s2", "s2_a")}).assign(
        subtype=["H1", "H3"]
    )
    assert member_barcodes(lib, {"subtype": "H1"}) == ["AA"]


def test_member_barcodes_compares_as_strings():
    """A value typed differently by YAML than by the CSV still matches."""
    lib = library({"AA": ("s1", "s1_a"), "CC": ("s2", "s2_a")}).assign(group=[1, 2])
    assert member_barcodes(lib, {"group": "1"}) == ["AA"]
    assert member_barcodes(lib, {"group": 1}) == ["AA"]


def test_member_barcodes_rejects_a_partly_selected_strain():
    """A strain only some of whose barcodes are selected has no meaningful even share."""
    lib = library({"AA": ("s1", "s1_a"), "CC": ("s1", "s1_b")}).assign(tube=["x", "y"])
    with pytest.raises(ValueError, match="only some"):
        member_barcodes(lib, {"tube": "x"})


@pytest.mark.parametrize(
    "members, match",
    [({"nope": "x"}, "not a column"), ({"strain": "absent"}, "selects no")],
)
def test_member_barcodes_rejects_bad_selections(members, match):
    """An unknown column, or one selecting nothing, raises."""
    lib = library({"AA": ("s1", "s1_a")})
    with pytest.raises(ValueError, match=match):
        member_barcodes(lib, members)


@pytest.mark.parametrize(
    "pool_by, units", [("barcode", 3), ("strain", 2), ("stock", 2)]
)
def test_pool_units_counts_the_units_of_each_kind(pool_by, units):
    """The even share is one over the number of `pool_by` units among the members."""
    lib = library({"AA": ("s1", "s1_a"), "CC": ("s1", "s1_a"), "GG": ("s2", "s2_a")})
    assert pool_units(lib, pool_by) == units


def test_pool_units_rejects_an_unknown_pool_by():
    """`pool_by` is one of the three units a pool can be balanced by."""
    with pytest.raises(ValueError, match="pool_by"):
        pool_units(library({"AA": ("s1", "s1_a")}), "well")


@pytest.mark.parametrize(
    "lib, match",
    [
        (pd.DataFrame({"barcode": ["AA"], "strain": ["s1"]}), "'stock' column"),
        (
            pd.DataFrame({"barcode": ["AA"], "strain": ["s1"], "stock": [""]}),
            "no stock",
        ),
        (
            pd.DataFrame(
                {
                    "barcode": ["AA", "CC"],
                    "strain": ["s1", "s2"],
                    "stock": ["shared", "shared"],
                }
            ),
            "more than one strain",
        ),
    ],
)
def test_pool_units_by_stock_requires_usable_stocks(lib, match):
    """Analyzing by stock needs a stock per barcode, and each stock within one strain."""
    with pytest.raises(ValueError, match=match):
        pool_units(lib, "stock")


# --- shares of the pool ------------------------------------------------------------------


def test_barcode_fractions_averages_over_the_analysis_wells():
    """A barcode's fraction is the mean over the wells of its share of each well."""
    counts = pd.DataFrame(
        [
            {"plate": "p1", "well": "A1", "barcode": "AA", "count": 10},
            {"plate": "p1", "well": "A1", "barcode": "CC", "count": 90},
            {"plate": "p1", "well": "A2", "barcode": "AA", "count": 30},
            {"plate": "p1", "well": "A2", "barcode": "CC", "count": 70},
        ]
    )
    fractions = barcode_fractions(counts, ["AA", "CC"], [("p1", "A1"), ("p1", "A2")])
    assert dict(zip(fractions["barcode"], fractions["mean_fraction"])) == {
        "AA": pytest.approx(0.2),
        "CC": pytest.approx(0.8),
    }
    assert fractions["mean_fraction"].sum() == pytest.approx(1.0)


def test_barcode_fractions_keeps_a_member_with_no_reads():
    """A member that did not grow is reported at zero rather than left out."""
    counts = pd.DataFrame([{"plate": "p1", "well": "A1", "barcode": "AA", "count": 10}])
    fractions = barcode_fractions(counts, ["AA", "CC"], [("p1", "A1")])
    assert dict(zip(fractions["barcode"], fractions["mean_fraction"])) == {
        "AA": pytest.approx(1.0),
        "CC": 0.0,
    }


def test_barcode_fractions_excludes_non_members_from_the_denominator():
    """Fractions are over member reads only, so a non-member cannot dilute a share."""
    counts = pd.DataFrame(
        [
            {"plate": "p1", "well": "A1", "barcode": "AA", "count": 25},
            {"plate": "p1", "well": "A1", "barcode": "CC", "count": 25},
            {"plate": "p1", "well": "A1", "barcode": "GG", "count": 50},
        ]
    )
    fractions = barcode_fractions(counts, ["AA", "CC"], [("p1", "A1")])
    assert fractions["mean_fraction"].sum() == pytest.approx(1.0)
    assert set(fractions["barcode"]) == {"AA", "CC"}


def test_barcode_fractions_rejects_a_well_with_no_member_reads():
    """A well passing QC on non-member reads alone measures nothing about the pool."""
    counts = pd.DataFrame([{"plate": "p1", "well": "A1", "barcode": "GG", "count": 50}])
    with pytest.raises(ValueError, match="no reads of any member"):
        barcode_fractions(counts, ["AA"], [("p1", "A1")])


def test_unit_shares_sums_a_unit_over_its_barcodes():
    """A unit's measured share is the sum of the fractions of its own barcodes."""
    lib = library({"AA": ("s1", "s1_a"), "CC": ("s1", "s1_b"), "GG": ("s2", "s2_a")})
    fractions = pd.DataFrame(
        {"barcode": ["AA", "CC", "GG"], "mean_fraction": [0.5, 0.2, 0.3]}
    )
    shares = unit_shares(fractions, lib, "stock")
    assert dict(zip(shares["unit"], shares["measured_share"])) == {
        "s1_a": pytest.approx(0.5),
        "s1_b": pytest.approx(0.2),
        "s2_a": pytest.approx(0.3),
    }


def test_unit_shares_divides_a_strain_among_its_stocks():
    """The divided share splits a strain's share equally among its stocks.

    This is what the existing project scripts compute and what a project balancing per
    strain pipettes, and it differs from the measured share whenever a strain's stocks
    grew differently: here 0.5 and 0.2 measured against 0.35 divided.

    """
    lib = library({"AA": ("s1", "s1_a"), "CC": ("s1", "s1_b"), "GG": ("s2", "s2_a")})
    fractions = pd.DataFrame(
        {"barcode": ["AA", "CC", "GG"], "mean_fraction": [0.5, 0.2, 0.3]}
    )
    shares = unit_shares(fractions, lib, "stock", equal_volume=True)
    assert dict(zip(shares["unit"], shares["divided_share"])) == {
        "s1_a": pytest.approx(0.35),
        "s1_b": pytest.approx(0.35),
        "s2_a": pytest.approx(0.3),
    }


def test_unit_shares_omits_the_divided_share_without_equal_volume():
    """Dividing a strain's share among its stocks assumes they went in equally."""
    lib = library({"AA": ("s1", "s1_a"), "CC": ("s1", "s1_b")})
    fractions = pd.DataFrame({"barcode": ["AA", "CC"], "mean_fraction": [0.8, 0.2]})
    shares = unit_shares(fractions, lib, "stock", equal_volume=False)
    assert "divided_share" not in shares.columns
    assert dict(zip(shares["unit"], shares["measured_share"])) == {
        "s1_a": pytest.approx(0.8),
        "s1_b": pytest.approx(0.2),
    }


@pytest.mark.parametrize("pool_by", ["barcode", "strain", "stock"])
def test_unit_shares_always_names_the_strain(pool_by):
    """Every `pool_by` gives a `strain` column, which the report and the CSV look up by.

    Grouping by the unit would otherwise drop it under `pool_by: strain`, where the unit is
    the strain, and the shares could then not be joined to the per-strain annotations.

    """
    lib = library({"AA": ("s1", "s1_a"), "CC": ("s1", "s1_b"), "GG": ("s2", "s2_a")})
    fractions = pd.DataFrame(
        {"barcode": ["AA", "CC", "GG"], "mean_fraction": [0.5, 0.2, 0.3]}
    )
    shares = unit_shares(fractions, lib, pool_by, equal_volume=True)
    assert {"unit", "strain", "measured_share"} <= set(shares.columns)
    assert set(shares["strain"]) == {"s1", "s2"}
    assert shares["measured_share"].sum() == pytest.approx(1.0)


def test_unit_shares_by_strain_has_equal_measured_and_divided_shares():
    """A strain is its own unit, so there is nothing to divide."""
    lib = library({"AA": ("s1", "s1_a"), "CC": ("s1", "s1_b")})
    fractions = pd.DataFrame({"barcode": ["AA", "CC"], "mean_fraction": [0.8, 0.2]})
    shares = unit_shares(fractions, lib, "strain", equal_volume=True)
    assert shares["measured_share"].tolist() == pytest.approx([1.0])
    assert shares["divided_share"].tolist() == pytest.approx([1.0])


# --- stocks of one strain that disagree ---------------------------------------------------


def test_diverging_stocks_flags_a_strain_whose_stocks_disagree():
    """Stocks of one strain differing by more than the factor are reported.

    The numbers are those of a real strain whose two tubes measure 4.6-fold apart while
    being pipetted equally, which is what makes the divergence worth surfacing.

    """
    lib = library({"AA": ("darwin", "darwin_1"), "CC": ("darwin", "darwin_2")})
    shares = pd.DataFrame(
        {
            "unit": ["darwin_1", "darwin_2"],
            "strain": ["darwin", "darwin"],
            "measured_share": [0.024637, 0.0053222],
        }
    )
    flagged = diverging_stocks(shares, lib, factor=3)
    assert flagged["strain"].tolist() == ["darwin"]
    assert flagged["ratio"].tolist() == pytest.approx([4.6291], rel=1e-3)


def test_diverging_stocks_leaves_a_strain_within_the_factor():
    """Ordinary spread between tubes is not reported."""
    lib = library({"AA": ("s1", "s1_a"), "CC": ("s1", "s1_b")})
    shares = pd.DataFrame(
        {
            "unit": ["s1_a", "s1_b"],
            "strain": ["s1", "s1"],
            "measured_share": [0.012, 0.010],
        }
    )
    assert len(diverging_stocks(shares, lib, factor=3)) == 0


def test_diverging_stocks_ignores_a_single_stock_strain():
    """A strain of one stock has nothing to disagree with."""
    lib = library({"AA": ("s1", "s1_a")})
    shares = pd.DataFrame({"unit": ["s1_a"], "strain": ["s1"], "measured_share": [0.5]})
    assert len(diverging_stocks(shares, lib, factor=3)) == 0


def test_diverging_stocks_reports_a_stock_with_no_reads():
    """A stock with no reads beside one with reads is named, with no ratio against zero."""
    lib = library({"AA": ("s1", "s1_a"), "CC": ("s1", "s1_b")})
    shares = pd.DataFrame(
        {
            "unit": ["s1_a", "s1_b"],
            "strain": ["s1", "s1"],
            "measured_share": [0.5, 0.0],
        }
    )
    flagged = diverging_stocks(shares, lib, factor=3)
    assert flagged["strain"].tolist() == ["s1"]
    assert pd.isna(flagged["ratio"]).all()


def test_diverging_stocks_ignores_a_strain_with_no_reads_at_all():
    """A strain whose every stock is empty is a low-titer result, not a disagreement."""
    lib = library({"AA": ("s1", "s1_a"), "CC": ("s1", "s1_b")})
    shares = pd.DataFrame(
        {
            "unit": ["s1_a", "s1_b"],
            "strain": ["s1", "s1"],
            "measured_share": [0.0, 0.0],
        }
    )
    assert len(diverging_stocks(shares, lib, factor=3)) == 0


def test_defaults_are_the_documented_ones():
    """The defaults the README states, applied in one place."""
    assert DEFAULT_LINEARITY_WINDOW_DILUTIONS == 3
    assert DEFAULT_LINEARITY_SLOPE_TOL == 0.3
    assert IDEAL_SLOPE == 1.0

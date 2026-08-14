"""Tests for ``scripts/seqneut_funcs.py``."""

import numpy
import pandas as pd
import pytest
from seqneut_funcs import (
    get_median_bound,
    narrow_for_altair,
    padded_log_domain,
    pearson_r_log10,
    round_sig,
    viruses_in_plot_order,
)


def test_round_sig():
    """Rounds to four significant figures, at any magnitude."""
    assert round_sig(210.35000000000002) == 210.4
    assert round_sig(0.000123456) == 0.0001235
    assert round_sig(123456.0) == 123500.0
    assert round_sig(0.0) == 0.0
    assert isinstance(round_sig(1), float)


class TestGetMedianBound:
    """The bound reported for a median titer."""

    def test_odd_count_takes_the_middle(self):
        """With an odd count the median is one value, so it keeps that value's bound."""
        assert get_median_bound(["lower", "interpolated", "upper"]) == "interpolated"

    def test_even_count_with_matching_bounds(self):
        """Two straddling values that agree give that bound."""
        assert get_median_bound(["lower", "lower", "lower", "lower"]) == "lower"

    def test_even_count_prefers_the_real_bound_over_interpolated(self):
        """Where one straddling value is censored, the median is reported as censored."""
        assert get_median_bound(["interpolated", "upper"]) == "upper"
        assert get_median_bound(["lower", "interpolated"]) == "lower"

    def test_even_count_with_opposing_bounds_is_inconsistent(self):
        """A median straddled by an upper and a lower bound is not a bound at all."""
        assert get_median_bound(["lower", "upper"]) == "inconsistent"

    def test_accepts_any_sequence(self):
        """Called on a `pandas` group, so it must not require a list."""
        assert get_median_bound(pd.Series(["upper", "upper"])) == "upper"


class TestPearsonRLog10:
    """Correlation of log-transformed titers."""

    def test_perfect_correlation(self):
        """A plate against an exact duplicate of itself correlates perfectly."""
        df = pd.DataFrame({"titer_x": [1.0, 10.0, 100.0]})
        df["titer_y"] = df["titer_x"]
        assert pearson_r_log10(df) == pytest.approx(1.0)

    def test_is_computed_on_the_logs(self):
        """Titers span orders of magnitude, so the correlation is of their logs.

        `titer_y` is the square of `titer_x`, which is a perfect straight line in log
        space but not in linear space, so the two disagree.
        """
        df = pd.DataFrame({"titer_x": [1.0, 10.0, 100.0], "titer_y": [1.0, 100.0, 1e4]})
        assert pearson_r_log10(df) == pytest.approx(1.0)
        assert numpy.corrcoef(df["titer_x"], df["titer_y"])[0, 1] != pytest.approx(1.0)

    @pytest.mark.parametrize("n", [0, 1, 2])
    def test_too_few_points_is_none(self, n):
        """Fewer than three points is not a correlation worth reporting."""
        df = pd.DataFrame({"titer_x": [10.0] * n, "titer_y": [20.0] * n})
        assert pearson_r_log10(df) is None

    def test_returns_a_plain_float(self):
        """The value is written to a CSV, so it must not be a numpy scalar."""
        df = pd.DataFrame({"titer_x": [1.0, 10.0, 100.0], "titer_y": [2.0, 8.0, 90.0]})
        assert type(pearson_r_log10(df)) is float

    def test_named_columns(self):
        """The columns can be named, as the two callers use different frames."""
        df = pd.DataFrame({"a": [1.0, 10.0, 100.0], "b": [1.0, 10.0, 100.0]})
        assert pearson_r_log10(df, x_col="a", y_col="b") == pytest.approx(1.0)


class TestPaddedLogDomain:
    """The shared domain of the two axes of a titer-versus-titer plot."""

    def test_pads_proportionally_to_the_log_range(self):
        """Two decades of titers are padded by a twentieth of that range at each end."""
        s = pd.Series([1.0, 100.0])
        assert padded_log_domain(s, s) == [0.7943, 125.9]

    def test_spans_both_axes(self):
        """The domain is shared, so it must cover the titers on either axis."""
        domain = padded_log_domain(pd.Series([1.0, 5.0]), pd.Series([0.5, 100.0]))
        assert domain[0] < 0.5
        assert domain[1] > 100.0

    def test_floor_keeps_a_flat_range_from_collapsing(self):
        """A plate against an exact duplicate of itself has no spread to pad."""
        s = pd.Series([10.0, 10.0])
        assert padded_log_domain(s, s) == [9.091, 11.0]

    def test_returns_plain_floats(self):
        """The domain is embedded in a chart specification, so it must be JSON-able."""
        s = pd.Series([1.0, 100.0])
        assert all(type(x) is float for x in padded_log_domain(s, s))


# a configured plot order, deliberately not alphabetical
ORDER = ["A/Texas/1/2024", "A/Bangkok/2/2023", "A/Bhutan/3/2022"]


class TestVirusesInPlotOrder:
    """The order the viruses are plotted in."""

    def test_uses_the_configured_order(self):
        """The configured order is not alphabetical, and is what is followed."""
        assert viruses_in_plot_order(ORDER[::-1], ORDER) == ORDER

    def test_restricts_to_the_viruses_present(self):
        """The order names every virus in the project, a subset of any one chart's."""
        assert viruses_in_plot_order(["A/Bhutan/3/2022"], ORDER) == ["A/Bhutan/3/2022"]

    def test_alphabetical_when_unconfigured(self):
        """`viral_strain_plot_order` is optional, and unset means alphabetical."""
        assert viruses_in_plot_order(ORDER, None) == sorted(ORDER)

    def test_missing_virus_raises_and_is_named(self):
        """A virus absent from the order would be dropped from it, so it must fail."""
        with pytest.raises(ValueError, match="A/Perth/4/2021"):
            viruses_in_plot_order([*ORDER, "A/Perth/4/2021"], ORDER)

    def test_accepts_a_column_and_dedupes(self):
        """Both callers pass a `virus` column, which repeats a virus across rows."""
        column = pd.Series(["A/Bangkok/2/2023"] * 3 + ["A/Texas/1/2024"])
        assert viruses_in_plot_order(column, ORDER) == [
            "A/Texas/1/2024",
            "A/Bangkok/2/2023",
        ]

    @pytest.mark.parametrize("order", [ORDER, None])
    def test_each_group_keeps_all_of_its_own_viruses(self, order):
        """Groups need not share a virus set, and neither may truncate the other.

        Ordering one group's viruses using another group's would silently drop those the
        other group lacks, leaving them off the chart's axis.
        """
        group_a = ["A/Texas/1/2024", "A/Bangkok/2/2023"]
        group_b = ["A/Bangkok/2/2023", "A/Bhutan/3/2022"]
        assert set(viruses_in_plot_order(group_a, order)) == set(group_a)
        assert set(viruses_in_plot_order(group_b, order)) == set(group_b)


class TestNarrowForAltair:
    """Narrowing a frame before it is embedded in a chart."""

    def test_drops_requested_columns(self):
        df = pd.DataFrame({"keep": [1.0], "drop_me": ["x"]})
        assert list(narrow_for_altair(df, drop=["drop_me"])) == ["keep"]

    def test_rounds_floats_only(self):
        """Floats are rounded; ints and strings are untouched."""
        df = pd.DataFrame(
            {"f": [210.35000000000002], "i": [1234567], "s": ["1234.5678"]}
        )
        out = narrow_for_altair(df)
        assert out["f"].tolist() == [210.4]
        assert out["i"].tolist() == [1234567]
        assert out["s"].tolist() == ["1234.5678"]

    def test_does_not_mutate_the_input(self):
        """The caller keeps using the unrounded frame to write its CSVs."""
        df = pd.DataFrame({"f": [210.35000000000002]})
        narrow_for_altair(df)
        assert df["f"].tolist() == [210.35000000000002]

    def test_shrinks_the_embedded_json(self):
        """The point of the rounding is the size of the JSON `altair` embeds."""
        df = pd.DataFrame({"f": [210.35000000000002] * 10})
        assert len(narrow_for_altair(df).to_json()) < len(df.to_json())

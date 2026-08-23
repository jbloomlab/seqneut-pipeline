"""Tests for ``scripts/viral_library_validation.py``.

Every check exists to stop a bad library, so each is exercised by an input that breaks it
and by nothing else, and `TestEveryCheckFails` fails if a check is ever added without such
an input.

"""

import pandas as pd
import pytest
from viral_library_validation import (
    LibraryError,
    SpecError,
    describe_column,
    library_checks,
    numeric_lines,
    report,
    summary_lines,
    validate_spec,
)

# a spec exercising every check, with short barcodes and sequences so that a test frame
# can be read at a glance
SPEC = {
    "columns": {
        "strain": {"nulls": False, "matches": [r"^\S+$"]},
        "barcode": {"nulls": False, "unique": "barcode", "matches": ["^[ACGT]{4}$"]},
        "plasmid_id": {"nulls": True, "unique": "barcode"},
        "subtype": {"nulls": False, "values": ["H1N1", "H3N2"]},
        "strain_type": {"nulls": False, "values": ["vaccine", "circulating"]},
        "vaccine_type": {
            "nulls": True,
            "values": ["egg"],
            "non_null_when": {"strain_type": "vaccine"},
        },
        "haplotype": {"nulls": False, "unique": "strain", "unique_within": "subtype"},
        "collection_date": {"nulls": False, "numeric": True},
        # no `unique: strain`, which translation would make redundant: an nt sequence
        # naming two strains forces its protein to name them too
        "nt_seq": {
            "nulls": False,
            "matches": ["^([ACGT]{3})+$"],
            "translates_to": "prot_seq",
        },
        "prot_seq": {
            "nulls": False,
            "unique": "strain",
            "matches": ["^[ACDEFGHIKLMNPQRSTVWY]+$"],
            "length": {"min": 2, "max": 3},
            "by": "subtype",
            "groups": {"H1N1": {"matches": ["^M"]}, "H3N2": {"matches": ["K$"]}},
        },
    }
}


@pytest.fixture
def library():
    """Three strains, the first carrying two barcodes, passing every check in `SPEC`."""
    return pd.DataFrame(
        {
            "barcode": ["ACGT", "ACGA", "TTTT", "GGGG"],
            "strain": ["A/X/1", "A/X/1", "A/Y/2", "A/Z/3"],
            "plasmid_id": ["p1", "p2", None, "p4"],
            "subtype": ["H1N1", "H1N1", "H3N2", "H1N1"],
            "strain_type": ["vaccine", "vaccine", "circulating", "circulating"],
            "vaccine_type": ["egg", "egg", None, None],
            "haplotype": ["A", "A", "A", "B"],
            "collection_date": ["2020.5", "2020.5", "2021.25", "2019.0"],
            "nt_seq": ["ATGAAA", "ATGAAA", "AAAAAA", "ATGTGC"],
            "prot_seq": ["MK", "MK", "KK", "MC"],
        }
    )


def _set(df, row, **values):
    """Copy of `df` with `values` assigned to the rows indexed by `row`."""
    df = df.copy()
    for col, value in values.items():
        df.loc[row, col] = value
    return df


# one input per check `id`, breaking that check and only that check. Where a column has
# to change on every barcode of a strain to leave the other checks alone, it does.
BREAKERS = {
    "nulls": lambda df: _set(df, 2, haplotype=None),
    "values:subtype": lambda df: _set(df, 2, subtype="H5N1"),
    "values:strain_type": lambda df: _set(df, 3, strain_type="unknown"),
    "values:vaccine_type": lambda df: _set(df, [0, 1], vaccine_type="cell"),
    "non_null_when:vaccine_type": lambda df: _set(df, 2, vaccine_type="egg"),
    "numeric:collection_date": lambda df: _set(df, 3, collection_date="unknown"),
    "matches:strain:^\\S+$": lambda df: _set(df, [0, 1], strain="A/X 1"),
    "matches:barcode:^[ACGT]{4}$": lambda df: _set(df, 0, barcode="acgt"),
    # lowercase still translates, so this fails the format check alone
    "matches:nt_seq:^([ACGT]{3})+$": lambda df: _set(df, 3, nt_seq="atgtgc"),
    "matches:prot_seq:^[ACDEFGHIKLMNPQRSTVWY]+$": lambda df: _set(
        df, 3, nt_seq="ATGTAA", prot_seq="M*"
    ),
    "length:prot_seq": lambda df: _set(df, 3, nt_seq="ATGTGCATGTGC", prot_seq="MCMC"),
    "matches:prot_seq:H1N1:^M": lambda df: _set(df, 3, nt_seq="TGCATG", prot_seq="CM"),
    "matches:prot_seq:H3N2:K$": lambda df: _set(df, 2, nt_seq="AAAATG", prot_seq="KM"),
    "unique_barcode:barcode": lambda df: _set(df, 1, barcode="ACGT"),
    "unique_barcode:plasmid_id": lambda df: _set(df, 3, plasmid_id="p1"),
    "unique_strain:haplotype": lambda df: _set(df, 3, haplotype="A"),
    # a synonymous change, so the two strains share a protein but not a coding sequence
    "unique_strain:prot_seq": lambda df: _set(df, 3, nt_seq="ATGAAG", prot_seq="MK"),
    "constant_per_strain": lambda df: _set(df, 1, collection_date="2021.0"),
    "translates_to:nt_seq": lambda df: _set(df, 3, prot_seq="MM"),
}


def validate(df, spec=SPEC, tmp_path=None, name="library.csv"):
    """Write `df` to a CSV and validate it, returning the report."""
    csv = tmp_path / name
    df.to_csv(csv, index=False)
    return report(csv, spec)


class TestValidLibrary:
    """A library that meets its spec."""

    def test_passes(self, library, tmp_path):
        """Every check passes, so nothing is raised."""
        text = validate(library, tmp_path=tmp_path)
        assert "FAIL" not in text

    def test_report_holds_every_section(self, library, tmp_path):
        """The report a project tracks says what the library is and what was checked."""
        text = validate(library, tmp_path=tmp_path)
        for section in ["schema", "summary of", "columns", "patterns", "uniqueness"]:
            assert section in text

    def test_report_gives_the_range_of_a_numeric_column(self, library, tmp_path):
        """A number is summarized rather than merely checked, over strains not barcodes."""
        assert "`collection_date` over 3 strains:  min 2019.00" in validate(
            library, tmp_path=tmp_path
        )

    def test_a_column_the_spec_omits_is_named_rather_than_ignored(
        self, library, tmp_path
    ):
        """An unchecked column is visible in the report instead of passing unnoticed."""
        text = validate(library.assign(note=["a", "b", "c", "d"]), tmp_path=tmp_path)
        assert (
            "columns not listed in `viral_library_validations: columns`: ['note']"
            in text
        )
        assert "NOT VALIDATED" in text


class TestEveryCheckFails:
    """Each configured check fires on an input that breaks it, and on nothing else."""

    @pytest.mark.parametrize("check_id", sorted(BREAKERS))
    def test_the_right_check_fails(self, check_id, library, tmp_path):
        """Asserting which check fired, as passing on an earlier one proves nothing."""
        with pytest.raises(LibraryError) as excinfo:
            validate(BREAKERS[check_id](library), tmp_path=tmp_path)
        assert [c.id for c in excinfo.value.failures] == [check_id]

    def test_no_check_goes_untested(self, library):
        """A check with no breaking input is a check that might never fire."""
        assert {c.id for c in library_checks(library, SPEC)} == set(BREAKERS)

    def test_a_missing_column_is_reported_before_anything_else(self, library, tmp_path):
        """Nothing can be checked against a column that is not there."""
        with pytest.raises(LibraryError, match="missing listed columns"):
            validate(library.drop(columns="haplotype"), tmp_path=tmp_path)


class TestNotCrashing:
    """Inputs that must be reported as failures rather than raised out of `report`."""

    def test_a_sequence_that_cannot_be_translated_is_reported(self, library, tmp_path):
        """`ZZZ` is three letters but not a codon, so `Bio` refuses to translate it."""
        df = _set(library, 3, nt_seq="ZZZATG")
        with pytest.raises(LibraryError) as excinfo:
            validate(df, tmp_path=tmp_path)
        assert "translates_to:nt_seq" in {c.id for c in excinfo.value.failures}

    def test_strain_may_be_declared_distinct_per_barcode(self, library, tmp_path):
        """A library with one barcode per strain can say so, which `strain` groups by."""
        spec = {"columns": {**SPEC["columns"]}}
        spec["columns"]["strain"] = {**spec["columns"]["strain"], "unique": "barcode"}
        # the fixture's first strain carries two barcodes, so this is what it catches
        with pytest.raises(LibraryError) as excinfo:
            validate(library, spec=spec, tmp_path=tmp_path)
        assert "unique_barcode:strain" in {c.id for c in excinfo.value.failures}


class TestFailureReporting:
    """What a failing run leaves behind."""

    def test_the_report_holds_the_checks_that_passed(self, library, tmp_path):
        """The wrapper writes this to the log, so a failure says what it got through."""
        with pytest.raises(LibraryError) as excinfo:
            validate(BREAKERS["unique_barcode:barcode"](library), tmp_path=tmp_path)
        assert "=== summary of" in excinfo.value.report
        assert " OK " in excinfo.value.report

    def test_every_failure_is_reported_not_just_the_first(self, library, tmp_path):
        """A library with several faults is fixed in one pass rather than in several."""
        df = _set(library, 3, plasmid_id="p1", collection_date="unknown")
        with pytest.raises(LibraryError) as excinfo:
            validate(df, tmp_path=tmp_path)
        assert {c.id for c in excinfo.value.failures} == {
            "unique_barcode:plasmid_id",
            "numeric:collection_date",
        }

    def test_a_failure_names_the_offending_strains_and_values(self, library, tmp_path):
        """The message has to be enough to find the rows in the CSV."""
        with pytest.raises(LibraryError) as excinfo:
            validate(
                BREAKERS["matches:barcode:^[ACGT]{4}$"](library), tmp_path=tmp_path
            )
        assert "A/X/1: 'acgt'" in excinfo.value.report


class TestValidateSpec:
    """A malformed spec is a configuration error, never a problem with the data."""

    @pytest.mark.parametrize(
        ("bad", "message"),
        [
            ({"columns": {"strain": {"nulls": False}}}, "is missing \\['barcode'\\]"),
            ({"columns": {}}, "non-empty mapping"),
            ({"columns": {}, "barcode_length": 16}, "exactly the key `columns`"),
            ({"columns": {"strain": {}, "barcode": {}}}, "must declare `nulls`"),
            (
                {
                    "columns": {
                        "strain": {"nulls": False, "unqiue": "strain"},
                        "barcode": {"nulls": False},
                    }
                },
                "unknown keys",
            ),
            (
                {
                    "columns": {
                        "strain": {"nulls": False, "unique": "strain"},
                        "barcode": {"nulls": False},
                    }
                },
                "cannot have `unique: strain`",
            ),
            (
                {
                    "columns": {
                        "strain": {"nulls": False},
                        "barcode": {"nulls": False, "unique": "row"},
                    }
                },
                "expected one of",
            ),
            (
                {
                    "columns": {
                        "strain": {"nulls": False},
                        "barcode": {"nulls": False, "unique_within": "strain"},
                    }
                },
                "without `unique: strain`",
            ),
            (
                {
                    "columns": {
                        "strain": {"nulls": False, "values": [True, False]},
                        "barcode": {"nulls": False},
                    }
                },
                "must all be strings",
            ),
            (
                {
                    "columns": {
                        "strain": {"nulls": False, "matches": ["["]},
                        "barcode": {"nulls": False},
                    }
                },
                "invalid regex",
            ),
            (
                {
                    "columns": {
                        "strain": {"nulls": False, "length": {"min": 3}},
                        "barcode": {"nulls": False},
                    }
                },
                "exactly `min` and `max`",
            ),
            (
                {
                    "columns": {
                        "strain": {"nulls": False, "translates_to": "nope"},
                        "barcode": {"nulls": False},
                    }
                },
                "not another column listed",
            ),
            (
                {
                    "columns": {
                        "strain": {"nulls": False, "by": "barcode", "groups": {}},
                        "barcode": {"nulls": False},
                    }
                },
                "must declare the `values`",
            ),
            (
                {
                    "columns": {
                        "strain": {"nulls": False, "groups": {}},
                        "barcode": {"nulls": False},
                    }
                },
                "must have `by` and `groups` together",
            ),
            (
                {
                    "columns": {
                        "strain": {
                            "nulls": False,
                            "non_null_when": {"barcode": "x", "strain": "y"},
                        },
                        "barcode": {"nulls": False},
                    }
                },
                "must name exactly one column",
            ),
        ],
    )
    def test_rejects(self, bad, message):
        """Each malformed spec is rejected with a message naming what is wrong."""
        with pytest.raises(SpecError, match=message):
            validate_spec(bad)

    @pytest.mark.parametrize("key", ["by", "unique_within"])
    def test_a_grouping_column_may_not_be_null(self, key):
        """A row with no group has no group's checks to meet, so it must not exist."""
        spec = {
            "columns": {
                "strain": {"nulls": False},
                "barcode": {"nulls": False},
                "subtype": {"nulls": True, "values": ["H1N1"]},
                "other": (
                    {
                        "nulls": False,
                        "by": "subtype",
                        "groups": {"H1N1": {"matches": ["^M"]}},
                    }
                    if key == "by"
                    else {
                        "nulls": False,
                        "unique": "strain",
                        "unique_within": "subtype",
                    }
                ),
            }
        }
        with pytest.raises(SpecError, match="must be declared `nulls: false`"):
            validate_spec(spec)

    def test_groups_must_cover_exactly_the_values_they_key_on(self):
        """A group for a value the column cannot hold would never be applied."""
        spec = {
            "columns": {
                "strain": {"nulls": False},
                "barcode": {"nulls": False},
                "subtype": {"nulls": False, "values": ["H1N1"]},
                "prot": {
                    "nulls": False,
                    "by": "subtype",
                    "groups": {
                        "H1N1": {"matches": ["^M"]},
                        "H3N2": {"matches": ["^Q"]},
                    },
                },
            }
        }
        with pytest.raises(SpecError, match="every value needs an entry"):
            validate_spec(spec)

    def test_the_valid_spec_is_accepted(self):
        """The spec the rest of these tests use is itself valid."""
        validate_spec(SPEC)

    def test_a_bad_spec_is_caught_before_the_library_is_read(self, tmp_path):
        """`SpecError` rather than a `FileNotFoundError` from a CSV never opened."""
        with pytest.raises(SpecError):
            report(tmp_path / "does_not_exist.csv", {"columns": {"strain": {}}})


class TestNumericLines:
    """The range reported for a column declared `numeric`."""

    def test_says_so_when_nothing_parses(self, library):
        """Otherwise it reads as a second, garbled failure beside the real one."""
        df = _set(library, [0, 1, 2, 3], collection_date="unknown")
        (line,) = numeric_lines(df, SPEC)
        assert "nan" not in line
        assert "no numeric values" in line


class TestDescribeColumn:
    """The `checks` cell of the report's per-column table."""

    def test_names_each_configured_check(self):
        """A reader sees what a column was held to without opening `config.yml`."""
        assert describe_column("barcode", SPEC) == (
            "non-null; distinct per row; matches ^[ACGT]{4}$"
        )
        assert describe_column("haplotype", SPEC) == (
            "non-null; names one strain within its subtype"
        )

    def test_says_so_when_a_column_is_not_checked(self):
        """The point of the table is that an unchecked column cannot hide in it."""
        assert describe_column("anything_else", SPEC) == "NOT VALIDATED"


class TestSummaryLines:
    """The breakdown of the library that heads the report."""

    def test_breaks_down_every_categorical_column(self, library):
        """Driven by `values`, so the summary names no column of its own."""
        text = "\n".join(summary_lines(library, SPEC, "library.csv"))
        assert "4 barcodes for 3 strains" in text
        for column in ["subtype", "strain_type", "vaccine_type"]:
            assert f"  {column}" in text

    def test_counts_strains_and_barcodes_separately(self, library):
        """A strain carrying several barcodes must not be counted several times."""
        text = "\n".join(summary_lines(library, SPEC, "library.csv"))
        assert "    H1N1         2         3" in text

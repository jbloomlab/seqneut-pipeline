"""Tests for ``scripts/miscellaneous_plate_funcs.py``."""

import pandas as pd
import pytest
from miscellaneous_plate_funcs import (
    ANOTHER_PLATE,
    NOT_ASSAYED,
    PlateReport,
    all_origins,
    check_discarded_barcodes,
    check_fates_add_up,
    min_hamming,
    plate_and_well,
    read_per_well_csvs,
    reads_per_well,
    split_discarded,
    table_html,
    trace_origin,
    well_position,
    wells_apart,
)


def write_well_csvs(tmp_path, suffix, frames):
    """Write `frames`, keyed by (plate, well), where the counting rule would put them."""
    paths = []
    for (plate, well), frame in frames.items():
        (tmp_path / plate).mkdir(exist_ok=True)
        path = tmp_path / plate / f"{well}{suffix}.csv"
        frame.to_csv(path, index=False)
        paths.append(str(path))
    return paths


@pytest.mark.parametrize(
    "well_1, well_2, apart",
    [("A1", "A1", 0), ("B2", "C3", 1), ("B2", "B4", 2), ("A1", "H12", 11)],
)
def test_wells_apart_counts_diagonals_as_one(well_1, well_2, apart):
    """A step onto any of the eight surrounding wells is one, diagonals included."""
    assert wells_apart(well_1, well_2) == apart


@pytest.mark.parametrize("well", ["I1", "A13", "A0", "1A", "B"])
def test_well_position_rejects_wells_off_the_plate(well):
    """A well that is not on a 96-well plate raises rather than being placed anywhere."""
    with pytest.raises(ValueError):
        well_position(well)


def test_plate_and_well_from_path():
    """The plate is the directory and the well the file name less its suffix."""
    path = "results/miscellaneous_plates/plate_1/B6_counts.csv"
    assert plate_and_well(path, "_counts") == ("plate_1", "B6")


def test_read_per_well_csvs_labels_rows_by_plate_and_well(tmp_path):
    """Rows of every file are concatenated and labeled by the plate and well read."""
    frames = {
        (plate, well): pd.DataFrame({"barcode": ["AAAA"], "count": [n]})
        for n, (plate, well) in enumerate(
            [("p1", "A1"), ("p1", "A2"), ("p2", "A1"), ("p2", "A2")]
        )
    }
    paths = write_well_csvs(tmp_path, "_counts", frames)
    df = read_per_well_csvs(paths, "_counts", ["barcode", "count"], frames)
    assert set(zip(df["plate"], df["well"], strict=True)) == set(frames)
    assert df.set_index(["plate", "well"]).loc[("p2", "A1"), "count"] == 2


def test_read_per_well_csvs_accepts_an_empty_file(tmp_path):
    """A well with nothing to list has an empty CSV, which is not a missing well."""
    frames = {
        ("p1", "A1"): pd.DataFrame({"barcode": ["AAAA"], "count": [3]}),
        ("p1", "A2"): pd.DataFrame({"barcode": [], "count": []}),
    }
    paths = write_well_csvs(tmp_path, "_invalid", frames)
    df = read_per_well_csvs(paths, "_invalid", ["barcode", "count"], frames)
    assert len(df) == 1


def test_read_per_well_csvs_with_every_file_empty(tmp_path):
    """Wells with nothing to list at all give an empty frame with the columns read."""
    frames = {
        ("p1", "A1"): pd.DataFrame({"barcode": [], "count": []}),
        ("p1", "A2"): pd.DataFrame({"barcode": [], "count": []}),
    }
    paths = write_well_csvs(tmp_path, "_invalid", frames)
    df = read_per_well_csvs(paths, "_invalid", ["barcode", "count"], frames)
    assert len(df) == 0
    assert {"barcode", "count", "plate", "well"} <= set(df.columns)


def test_read_per_well_csvs_rejects_a_well_not_listed(tmp_path):
    """A file for a well the samples CSVs do not list raises."""
    frames = {
        ("p1", "A1"): pd.DataFrame({"barcode": ["AAAA"], "count": [3]}),
        ("p1", "A2"): pd.DataFrame({"barcode": ["AAAA"], "count": [1]}),
    }
    paths = write_well_csvs(tmp_path, "_counts", frames)
    with pytest.raises(ValueError, match="do not list"):
        read_per_well_csvs(paths, "_counts", ["barcode", "count"], [("p1", "A1")])


def test_read_per_well_csvs_requires_exactly_the_wells(tmp_path):
    """A well listed but not read, or read but not listed, raises."""
    frames = {("p1", "A1"): pd.DataFrame({"barcode": ["AAAA"], "count": [3]})}
    paths = write_well_csvs(tmp_path, "_counts", frames)
    with pytest.raises(ValueError, match="do not list"):
        read_per_well_csvs(
            paths, "_counts", ["barcode", "count"], [("p1", "A1"), ("p1", "A2")]
        )


def test_read_per_well_csvs_requires_columns_per_file(tmp_path):
    """A column missing from one file raises, rather than reading as nulls."""
    frames = {
        ("p1", "A1"): pd.DataFrame({"barcode": ["AAAA"], "count": [3]}),
        ("p1", "A2"): pd.DataFrame({"barcode": ["AAAA"]}),
    }
    paths = write_well_csvs(tmp_path, "_counts", frames)
    with pytest.raises(ValueError, match="lacks the columns"):
        read_per_well_csvs(paths, "_counts", ["barcode", "count"], frames)


def test_min_hamming():
    """The distance is to the nearest target."""
    assert min_hamming(["AAAA", "AATT", "CCCC"], ["AAAA", "AATA"]).tolist() == [0, 1, 4]
    assert min_hamming([], ["AAAA"]).tolist() == []


def test_min_hamming_rejects_mismatched_lengths():
    """Barcodes of another length than the targets have no Hamming distance to them."""
    with pytest.raises(ValueError):
        min_hamming(["AAAA"], ["AAA"])


def test_reads_per_well_counts_a_well_without_rows_as_zero():
    """Counts are summed per well, in the order given, a well without rows reading zero."""
    frame = pd.DataFrame(
        {"plate": ["p1", "p1", "p2"], "well": ["A1", "A1", "B2"], "count": [2, 3, 4]}
    )
    index = pd.MultiIndex.from_tuples(
        [("p2", "B2"), ("p1", "C3"), ("p1", "A1")], names=["plate", "well"]
    )
    assert reads_per_well(frame, index).tolist() == [4, 0, 5]


def test_split_discarded_takes_exact_matches_only():
    """Only a discarded barcode itself is set aside; its misreads are kept."""
    invalid = pd.DataFrame(
        {
            "barcode": ["GGTCCATC", "GGTCCATG", "GGTCCAAG", "TTTTTTTT"],
            "count": [10, 3, 2, 5],
        }
    )
    kept, discarded = split_discarded(invalid, ["GGTCCATC"])
    assert discarded["barcode"].tolist() == ["GGTCCATC"]
    assert kept["barcode"].tolist() == ["GGTCCATG", "GGTCCAAG", "TTTTTTTT"]


def test_split_discarded_with_nothing_to_discard():
    """No discarded barcodes, or no unmatched reads, discards nothing."""
    invalid = pd.DataFrame({"barcode": ["AAAA"], "count": [1]})
    kept, discarded = split_discarded(invalid, [])
    assert len(kept) == 1 and len(discarded) == 0
    kept, discarded = split_discarded(invalid.iloc[:0], ["AAAA"])
    assert len(kept) == 0 and len(discarded) == 0


def test_check_fates_add_up():
    """Tallies equal to the fates pass, a fate absent from a well counting as zero."""
    fates = pd.DataFrame(
        {
            "plate": ["p1", "p1", "p1"],
            "well": ["A1", "A1", "A2"],
            "fate": ["valid barcode", "invalid barcode", "valid barcode"],
            "count": [10, 2, 7],
        }
    )
    index = pd.MultiIndex.from_tuples([("p1", "A1"), ("p1", "A2")])
    check_fates_add_up(
        fates,
        {
            "valid barcode": pd.Series([10, 7], index=index),
            "invalid barcode": pd.Series([2, 0], index=index),
        },
    )
    with pytest.raises(ValueError, match=r"invalid barcode.*'p1', 'A2'"):
        check_fates_add_up(fates, {"invalid barcode": pd.Series([2, 1], index=index)})


def test_report_puts_warnings_at_the_reserved_slot(tmp_path):
    """Warnings added last are written where the slot was reserved, above what follows."""
    report = PlateReport(title="t")
    report.md("## First")
    report.reserve_warnings()
    report.md("## Second")
    report.warning("off the plate", pd.DataFrame({"plate": ["p1"], "well": ["Z9"]}))
    report.write(out := tmp_path / "r.html")
    html = out.read_text()
    assert html.index(">First</h2>") < html.index(">Warnings</h3>")
    assert html.index("<b>1</b> well off the plate") < html.index("Z9")
    assert html.index("Z9") < html.index(">Second</h2>")


def test_report_without_warnings_has_no_warnings_section(tmp_path):
    """The reserved slot stays empty when nothing was warned about."""
    report = PlateReport(title="t")
    report.reserve_warnings()
    report.md("## Only")
    report.write(out := tmp_path / "r.html")
    assert "Warnings" not in out.read_text()


def test_report_warnings_need_a_slot(tmp_path):
    """A warning with nowhere to go raises rather than being dropped."""
    report = PlateReport(title="t")
    report.warning("lost", pd.DataFrame({"well": ["A1"]}))
    with pytest.raises(ValueError, match="no slot"):
        report.write(tmp_path / "r.html")


def test_report_table_formats_floats_and_marks_rows(tmp_path):
    """Floats are written to four significant digits, and marked rows are colored."""
    report = PlateReport(title="t")
    df = pd.DataFrame({"strain": ["a<b", "c"], "frac": [0.123456, 2.0]})
    report.table(df, marked=pd.Series([False, True]))
    report.write(out := tmp_path / "r.html")
    html = out.read_text()
    assert "0.1235" in html
    assert "a&lt;b" in html
    assert '<span style="color: #cc00cc">c</span>' in html
    assert '<span style="color: #cc00cc">a' not in html


@pytest.mark.parametrize(
    "homes, expected",
    [
        ([("p1", "B2")], ("1 well away", "B2", 1, "B2, 1 well away")),
        ([("p1", "B5")], ("2 wells away", "B5", 2, "B5, 2 wells away")),
        ([("p1", "F9")], ("over 2 wells away", "F9", 6, "F9, over 2 wells away")),
        # a home on this plate wins over one on another
        ([("p2", "A1"), ("p1", "B2")], ("1 well away", "B2", 1, "B2, 1 well away")),
        ([("p2", "A1")], (ANOTHER_PLATE, "A1", None, "A1 of p2")),
        ([], (NOT_ASSAYED, None, None, NOT_ASSAYED)),
    ],
)
def test_trace_origin(homes, expected):
    """Material is traced to its home well on this plate, another plate, or nowhere."""
    assert trace_origin("p1", "B3", homes) == expected


def test_trace_origin_takes_the_nearest_home_on_a_plate():
    """Of several homes on this plate, the nearest is the source, the first of a tie."""
    nearest = ("1 well away", "B4", 1, "B4, 1 well away")
    assert trace_origin("p1", "B3", [("p1", "F9"), ("p1", "B4")]) == nearest
    assert trace_origin("p1", "B3", [("p1", "B4"), ("p1", "B2")]) == nearest


def test_all_origins_covers_every_traced_origin():
    """Every origin `trace_origin` gives is in `all_origins`, nearest first."""
    origins = all_origins()
    traced = {
        trace_origin("p1", "A1", homes)[0]
        for homes in [[("p1", w)] for w in ["A2", "A3", "A12"]] + [[("p2", "A1")], []]
    }
    assert traced == set(origins)
    assert origins.index("1 well away") < origins.index(NOT_ASSAYED)


COUNTED = ["AAAAAAAA", "CCCCCCCC"]


def test_check_discarded_barcodes_upper_cases():
    """A discarded barcode passes, upper-cased."""
    assert check_discarded_barcodes({"gggggggg": " a reason "}, COUNTED) == {
        "GGGGGGGG": "a reason"
    }
    assert check_discarded_barcodes({}, COUNTED) == {}


@pytest.mark.parametrize("barcode", ["AAAAAGGG", "AAAAAAGG", "AAAAAAAG"])
def test_check_discarded_barcodes_any_distance(barcode):
    """A discarded barcode can be any distance from a counted one, even 1 nt."""
    assert list(check_discarded_barcodes({barcode: "r"}, COUNTED)) == [barcode]


@pytest.mark.parametrize("barcode", ["AAAAAAAA", "cccccccc"])
def test_check_discarded_barcodes_rejects_a_counted_barcode(barcode):
    """A discarded barcode that is itself counted raises."""
    with pytest.raises(ValueError, match="counted rather than discarded"):
        check_discarded_barcodes({barcode: "r"}, COUNTED)


@pytest.mark.parametrize("barcode", ["GGGG", "GGGGGGGN", 12345678])
def test_check_discarded_barcodes_rejects_non_barcodes(barcode):
    """A wrong length, a letter other than ACGT, or a non-string raises."""
    with pytest.raises(ValueError):
        check_discarded_barcodes({barcode: "r"}, COUNTED)


@pytest.mark.parametrize("reason", [None, "", "  ", 3])
def test_check_discarded_barcodes_needs_a_reason(reason):
    """Every discarded barcode needs a reason, given as text."""
    with pytest.raises(ValueError, match="needs the reason"):
        check_discarded_barcodes({"GGGGGGGG": reason}, COUNTED)


def test_check_discarded_barcodes_rejects_a_list_and_repeats():
    """A bare list of barcodes, or one barcode given twice, raises."""
    with pytest.raises(TypeError, match="must map each barcode"):
        check_discarded_barcodes(["GGGGGGGG"], COUNTED)
    with pytest.raises(ValueError, match="more than once"):
        check_discarded_barcodes({"GGGGGGGG": "r", "gggggggg": "r"}, COUNTED)


def test_report_stat_tiles(tmp_path):
    """Each tile holds its value over its label, escaped, with the caption after."""
    report = PlateReport(title="t")
    report.stat_tiles({"plates": 2, "a<b": "12 / 36"}, caption="library: x.csv")
    report.write(out := tmp_path / "r.html")
    html = out.read_text()
    assert html.index(">2</div>") < html.index(">plates</div>")
    assert "a&lt;b" in html
    assert html.index("12 / 36") < html.index("library: x.csv")


def test_report_sidebar(tmp_path):
    """Sections get ids, linked from a sidebar that leaves out what is headed within one."""
    report = PlateReport(title="t")
    report.md("## First &amp; one\n\n### Warnings\n\n## Second\n\n## Second")
    report.write(out := tmp_path / "r.html")
    html = out.read_text()
    nav = html[html.index('<nav class="sidebar">') : html.index("</nav>")]
    for anchor in ["first-one", "second", "second-2"]:
        assert f'id="{anchor}"' in html
        assert f'href="#{anchor}"' in nav
    assert nav.index("first-one") < nav.index('"#second"')
    assert "<h3>Warnings</h3>" in html
    assert "Warnings" not in nav
    assert 'href="index.html"' in nav
    # writing again gives the same page, rather than ids on ids
    report.write(again := tmp_path / "r2.html")
    assert again.read_text() == html


def test_report_notes(tmp_path):
    """What is appended within a notes block sits inside a closed details element."""
    report = PlateReport(title="t")
    report.md("before")
    with report.notes("How to read <this>"):
        report.md("hidden one")
        report.md("hidden two")
    with report.notes("more"):
        report.md("hidden three")
    report.md("after")
    report.write(out := tmp_path / "r.html")
    html = out.read_text()
    assert "How to read &lt;this&gt;" in html
    assert html.count("details.notes {") == 1
    assert html.count('<details class="notes">') == 2
    assert '<details class="notes" open' not in html
    first_open = html.index('<details class="notes">')
    first_close = html.index("</details>", first_open)
    assert html.index("before") < first_open
    assert (
        first_open < html.index("hidden one") < html.index("hidden two") < first_close
    )
    assert html.rindex("</details>") < html.index("after")


def test_table_html_group():
    """Grouped columns sit under one header, and rows agreeing elsewhere share cells."""
    df = pd.DataFrame(
        {
            "well": ["B5", "B5", "C1"],
            "reads": [10, 10, 5],
            "barcode": ["AAA", "CCC", "GGG"],
            "fraction": [0.123456, 0.2, 0.3],
        }
    )
    out = table_html(df, group=("top <barcodes>", ["barcode", "fraction"]))
    assert '<th colspan="2" style="text-align: center">top &lt;barcodes&gt;</th>' in out
    assert out.count('<th rowspan="2">') == 2
    assert out.count('<td rowspan="2" style="vertical-align: top">') == 2
    assert out.count('<td rowspan="1" style="vertical-align: top">') == 2
    assert out.index("AAA") < out.index("CCC") < out.index("C1") < out.index("GGG")
    assert "0.1235" in out
    with pytest.raises(ValueError, match="cannot group"):
        table_html(df, group=("x", ["missing"]))
    with pytest.raises(ValueError, match="cannot show its index"):
        table_html(df, index=True, group=("x", ["barcode"]))


def test_table_html_group_relabeled():
    """A grouped column can be shown under a sub-header naming a column outside it."""
    df = pd.DataFrame({"well": ["B5", "B5"], "source": ["B4", "C5"]})
    out = table_html(df, group=("came from", {"source": "well"}))
    assert out.count("<th") - out.count("<thead") == 3
    assert '<th><div class="head">well</div></th>' in out
    assert "source" not in out


def test_report_section_warnings(tmp_path):
    """Section warnings are written where they are added, with their note, not in the slot."""
    report = PlateReport(title="t")
    report.reserve_warnings()
    report.md("## Titers")
    report.section_warnings([])
    report.section_warnings(
        [
            {
                "text": "with zero titer",
                "details": pd.DataFrame({"well": ["A1", "B2"]}),
                "note": "absent from <log>",
            }
        ]
    )
    report.md("## After")
    report.write(out := tmp_path / "r.html")
    html = out.read_text()
    assert html.count(">Warnings</h3>") == 1
    assert html.index(">Titers</h2>") < html.index(">Warnings</h3>")
    assert html.index("<b>2</b> wells with zero titer") < html.index(
        "absent from &lt;log&gt;"
    )
    assert (
        html.index("absent from &lt;log&gt;")
        < html.index("A1")
        < html.index(">After</h2>")
    )
    assert len(report.warnings) == 1


def test_table_html_headers_wrap(tmp_path):
    """Headers are escaped and wrapped in the element that caps their width."""
    out = table_html(pd.DataFrame({"non-neut <reads>": [3]}))
    assert '<div class="head">non&#8209;neut &lt;reads&gt;</div>' in out

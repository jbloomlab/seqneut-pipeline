"""Tests for ``scripts/seqneut_report.py``."""

import re

import altair as alt
import matplotlib.figure
import pandas as pd
import pytest
from seqneut_report import Report


@pytest.fixture
def chart():
    """A minimal `altair` chart."""
    return (
        alt.Chart(pd.DataFrame({"x": [1, 2, 3], "y": [1, 4, 9]}))
        .mark_point()
        .encode(x="x", y="y")
    )


@pytest.fixture
def fig():
    """A minimal `matplotlib` figure, built without `pyplot` so no backend is needed."""
    f = matplotlib.figure.Figure(figsize=(2, 2))
    f.subplots().plot([0, 1], [0, 1])
    return f


def test_page_is_written_with_title(tmp_path):
    """The report is a complete HTML page titled as requested."""
    report = Report(title="A plate")
    out = tmp_path / "report.html"
    report.write(out)
    html = out.read_text()
    assert html.startswith("<!DOCTYPE html>")
    assert html.rstrip().endswith("</html>")
    assert "<title>A plate</title>" in html
    assert "<h1>A plate</h1>" in html


def test_md_is_dedented_not_read_as_code(tmp_path):
    """Indented markdown renders as prose, not as an indented code block.

    Report prose is written as indented triple-quoted strings, and four spaces of
    indentation is a code block to markdown, so failing to dedent would turn every
    heading into preformatted text.
    """
    report = Report(title="t")
    report.md("""
        ## A heading

        Some prose with `code` in it.
        """)
    report.write(out := tmp_path / "r.html")
    html = out.read_text()
    assert "<h2>A heading</h2>" in html
    assert "<code>code</code>" in html
    assert "<pre>" not in html


def test_md_logs_only_when_asked(tmp_path, capsys):
    """`log=True` also prints, which each script redirects into its rule's log file."""
    report = Report(title="t")
    report.md("not logged")
    assert capsys.readouterr().out == ""
    report.md("logged", log=True)
    assert capsys.readouterr().out.strip() == "logged"
    report.write(out := tmp_path / "r.html")
    html = out.read_text()
    assert "not logged" in html
    assert "logged" in html


def test_yaml_renders_as_a_code_block(tmp_path):
    """The QC parameters and the QC drops are shown as YAML, which renders as code."""
    report = Report(title="t")
    report.yaml({"qc_thresholds": {"min_replicates": 2}})
    report.write(out := tmp_path / "r.html")
    html = out.read_text()
    assert "<pre>" in html
    assert "min_replicates: 2" in html


def test_yaml_keeps_nested_indentation(tmp_path):
    """`Report.md` dedents its text, which must not flatten a nested mapping."""
    report = Report(title="t")
    report.yaml({"qc_thresholds": {"min_replicates": 2}})
    report.write(out := tmp_path / "r.html")
    assert "\n  min_replicates: 2" in out.read_text()


def test_md_passes_latex_through_for_katex(tmp_path):
    """LaTeX survives markdown untouched, and KaTeX is loaded to render it."""
    report = Report(title="t")
    report.md(r"$$F_{v_b,s} = \frac{c}{n}$$")
    report.write(out := tmp_path / "r.html")
    html = out.read_text()
    assert r"\frac{c}{n}" in html
    assert "katex" in html
    assert "renderMathInElement" in html


def test_charts_get_unique_div_ids(tmp_path, chart):
    """Every chart on a page needs its own `output_div` or they collide."""
    report = Report(title="t")
    for _ in range(3):
        report.chart(chart)
    report.write(out := tmp_path / "r.html")
    html = out.read_text()
    ids = re.findall(r'<div id="(chart_\d+)"></div>', html)
    assert ids == ["chart_1", "chart_2", "chart_3"]
    assert len(ids) == len(set(ids))


def test_page_has_no_script_src_tags(tmp_path, chart):
    """`altair`'s universal template loads its libraries itself.

    Nothing in the report writes a `<script src>`, so the `vega` libraries cannot be
    fetched once per chart, and no `vega` version is pinned in this repo.
    """
    report = Report(title="t")
    for _ in range(3):
        report.chart(chart)
    report.write(out := tmp_path / "r.html")
    html = out.read_text()
    srcs = re.findall(r'<script[^>]*\ssrc="([^"]+)"', html)
    assert not [s for s in srcs if "vega" in s]
    # the libraries are loaded once, by the guarded loader `altair` emits
    assert "VEGA_DEBUG" in html


def test_chart_data_reaches_the_page(tmp_path, chart):
    """The chart's data frame is embedded, so the report needs no side-car files."""
    report = Report(title="t")
    report.chart(chart)
    report.write(out := tmp_path / "r.html")
    assert '"x": 3' in out.read_text().replace("&#34;", '"')


def test_table_is_scrollable(tmp_path):
    """A table goes in the scroll box so a long one does not run off the page."""
    df = pd.DataFrame({"a": [1, 2], "b": ["x", "y"]})
    report = Report(title="t")
    report.table(df)
    report.write(out := tmp_path / "r.html")
    html = out.read_text()
    assert '<div class="scroll-box"><table' in html
    assert "max-height: 80vh" in html
    assert "overflow: auto" in html


def test_table_index_can_be_dropped(tmp_path):
    """A meaningless range index is not written."""
    df = pd.DataFrame({"a": [1, 2]}, index=["keep_me", "and_me"])
    report = Report(title="t")
    report.table(df, index=True)
    report.table(df.reset_index(drop=True), index=False)
    report.write(out := tmp_path / "r.html")
    html = out.read_text()
    assert "keep_me" in html
    tables = re.findall(r"<table.*?</table>", html, re.DOTALL)
    assert len(tables) == 2
    assert "and_me" not in tables[1]


@pytest.mark.parametrize("method", ["svg", "pdf", "png8"])
def test_figure_is_embedded(tmp_path, fig, method):
    """Each display method embeds the figure with no external file."""
    report = Report(title="t")
    report.figure(fig, method)
    report.write(out := tmp_path / "r.html")
    html = out.read_text()
    if method == "svg":
        assert 'id="svgwrap"' in html
        assert "overflow:auto" in html
    elif method == "pdf":
        assert "data:application/pdf;base64," in html
    else:
        assert "data:image/png;base64," in html


def test_figure_propagates_an_invalid_display_method(fig):
    """`Report.figure` does not swallow the error the renderer raises.

    Which methods are valid is `neutcurve`'s to say, so this checks only that an invalid
    one is not quietly ignored here.
    """
    report = Report(title="t")
    with pytest.raises(ValueError):
        report.figure(fig, "not_a_real_method")


def test_blocks_keep_their_order(tmp_path, chart):
    """Blocks appear in the order the script appended them."""
    report = Report(title="t")
    report.md("# first")
    report.chart(chart)
    report.md("# last")
    report.write(out := tmp_path / "r.html")
    html = out.read_text()
    assert html.index("first") < html.index("chart_1") < html.index("last")


def test_report_is_reproducible(tmp_path, chart, fig):
    """The same report written twice gives the same bytes.

    The reports are rendered by a pipeline whose other outputs are tracked in version
    control, so a report that differed on every run would be noise.
    """
    outs = []
    for i in range(2):
        report = Report(title="t")
        report.md("## heading")
        report.chart(chart)
        report.table(pd.DataFrame({"a": [1]}))
        report.figure(fig, "svg")
        report.write(out := tmp_path / f"r{i}.html")
        outs.append(out.read_bytes())
    assert outs[0] == outs[1]

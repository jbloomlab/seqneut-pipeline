"""Tests for ``scripts/miscellaneous_plate_charts.py``."""

import altair as alt
import pandas as pd
import pytest
from miscellaneous_plate_charts import (
    LABEL_LIMIT,
    RETAINED_COLOR,
    WellPanels,
    carried_forward_stroke,
    category_fill,
    lookup_tooltip,
    stroke_width,
    well_hover,
)


def test_well_panels_rejects_another_axis():
    """The wells run along `x` or `y` and nothing else."""
    with pytest.raises(ValueError, match="'strain_axis' must be 'x' or 'y'"):
        WellPanels("z")


@pytest.mark.parametrize(
    "strain_axis, value_channel, concat, labeling_panel, angle",
    [("y", alt.X, alt.hconcat, 0, 0), ("x", alt.Y, alt.vconcat, 2, -90)],
)
def test_well_panels_layout(strain_axis, value_channel, concat, labeling_panel, angle):
    """Side by side along `y`, labeled on the first; stacked along `x`, on the last."""
    panels = WellPanels(strain_axis)
    assert panels.value_channel is value_channel
    assert panels.concat is concat
    assert [panels.labels_wells(i, 3) for i in range(3)] == [
        i == labeling_panel for i in range(3)
    ]
    axis = panels.well_axis(["a", "b"], labeled=True).to_dict()
    assert axis["sort"] == ["a", "b"]
    assert axis["axis"] == {"labelAngle": angle, "labelLimit": LABEL_LIMIT}
    assert panels.well_axis(["a"], labeled=False).to_dict()["axis"] is None
    across = "width" if strain_axis == "x" else "height"
    assert panels.panel_size(220)[across].to_dict() == {"step": 10}
    assert panels.bar_band(0.85) == {across: {"band": 0.85}}


def test_well_hover_and_stroke():
    """The hover is declared once and outlines the hovered well heavier."""
    hover = well_hover()
    chart = (
        alt.Chart(pd.DataFrame({"unit_label": ["a"], "n": [1]}))
        .mark_bar()
        .encode(
            alt.X("n:Q"),
            alt.Y("unit_label:N"),
            stroke=carried_forward_stroke(False, "top"),
            strokeWidth=stroke_width(hover, 2.5, 0.5),
        )
        .add_params(hover)
    )
    spec = chart.to_dict()
    assert [p["name"] for p in spec["params"]] == ["well_hover"]
    width = spec["encoding"]["strokeWidth"]
    assert width["condition"]["value"] == 2.5 and width["value"] == 0.5
    assert width["condition"]["empty"] is False
    assert spec["encoding"]["stroke"] == {"value": RETAINED_COLOR}


def test_carried_forward_stroke_with_final_library():
    """With a final library the outline is colored by whether a well was carried forward."""
    stroke = carried_forward_stroke(True, "right").to_dict()
    assert stroke["field"] == "carried_forward"
    assert stroke["scale"]["domain"] == ["carried forward", "dropped"]
    assert stroke["legend"]["orient"] == "right"


def test_lookup_tooltip_names_the_type():
    """A looked-up field's type comes from the lookup frame."""
    frame = pd.DataFrame({"n": [1, 2], "name": ["a", "b"]})
    assert lookup_tooltip(frame, "n").to_dict()["type"] == "quantitative"
    assert lookup_tooltip(frame, "name").to_dict()["type"] == "nominal"


def test_category_fill_from_a_scheme():
    """A scheme colors the categories and the sort order is the order given."""
    fill = category_fill("c", "a title", ["x", "y"], "right", scheme="dark2").to_dict()
    assert fill["field"] == "c"
    assert fill["title"] == "a title"
    assert fill["sort"] == ["x", "y"]
    assert fill["scale"] == {"scheme": "dark2"}
    assert fill["legend"] == {
        "labelLimit": LABEL_LIMIT,
        "orient": "right",
        "columns": 1,
    }


def test_category_fill_from_explicit_colors():
    """Explicit colors are mapped to the categories in order."""
    fill = category_fill("c", "t", ["x", "y"], "bottom", range_=["#000", "#fff"])
    assert fill.to_dict()["scale"] == {"domain": ["x", "y"], "range": ["#000", "#fff"]}


def test_category_fill_needs_exactly_one_of_scheme_and_range():
    """Giving both colorings, or neither, is an error rather than a silent precedence."""
    with pytest.raises(ValueError, match="exactly one"):
        category_fill("c", "t", ["x"], "right")
    with pytest.raises(ValueError, match="exactly one"):
        category_fill("c", "t", ["x"], "right", scheme="dark2", range_=["#000"])

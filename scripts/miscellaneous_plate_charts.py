"""Charts shared by the analyses of miscellaneous plates.

Every chart of these reports puts the wells analyzed on one shared axis, a step per well,
keyed on a ``unit_label`` field that names each well once, and draws its values in panels
beside or above one another along that axis. Hovering a well outlines it in every panel, and
where a final library is given each well's outline says whether it was carried forward.

"""

import altair as alt
import pandas

# pixels allowed for axis labels, enough for the longest strain name
LABEL_LIMIT = 300

# pixels per well along the shared axis
WELL_STEP = 10

# Width of the outline around a well's bars and points, ordinarily and when hovered. Nonzero
# either way, the outline being what carries the carried-forward mark. Thin around a bar,
# where it would otherwise crowd the divisions of a stacked one, and heavier around a point.
STROKE_WIDTH = 0.5
STROKE_WIDTH_HOVERED = 2.5
POINT_STROKE_WIDTH = 1.5
POINT_STROKE_WIDTH_HOVERED = 3

# How a unit is marked by a final library: by the color of the outline around its bars and
# points, and of its text in the tables. The dropped color has to read as body text as well
# as a thin outline, and to stay clear of the fill palettes the panels use.
CARRIED_FORWARD_LABELS = {True: "carried forward", False: "dropped"}
RETAINED_COLOR = "black"
DROPPED_COLOR = "#cc00cc"

# The color of each fate of a read, and the order the legends and notes list them in: the
# read that was counted first, then the ways reading one can fail. Fixed per fate so that a
# fate reads the same across reports and whether or not the other fates occur, a color taken
# from a `vega` scheme being assigned by a category's position among those present and so
# shifting as soon as one fate has no reads. The colors are from `tableau10`.
FATE_COLORS = {
    "valid barcode": "#4c78a8",
    "unparseable barcode": "#f58518",
    "read too short": "#eeca3b",
    "low quality barcode": "#e45756",
    "invalid outer flank": "#72b7b2",
    "invalid barcode": "#54a24b",
    "failed chastity filter": "#b279a2",
    "discarded barcode": "#ff9da6",
}


class WellPanels:
    """How the panels of a chart of wells are laid out.

    Along `y` each well is a row of horizontal bars and the panels of a chart sit side by
    side; along `x` each well is a column of vertical bars and the panels are stacked, so
    that the axis they share stays shared.

    Args:
        `strain_axis` (str)
            The axis the wells run along, ``"x"`` or ``"y"``.

    """

    def __init__(self, strain_axis):
        if strain_axis not in {"x", "y"}:
            raise ValueError(f"'strain_axis' must be 'x' or 'y', not {strain_axis!r}")
        self.strain_axis = strain_axis
        self.vertical = strain_axis == "x"
        # the channel a panel draws its values along
        self.value_channel = alt.Y if self.vertical else alt.X
        self.concat = alt.vconcat if self.vertical else alt.hconcat
        # how the panels of a chart are referred to in the text, first to last
        self.positions = (
            ["top", "middle", "bottom"]
            if self.vertical
            else ["left", "middle", "right"]
        )
        self.legend_orient = "right" if self.vertical else "top"

    def well_axis(self, order, labeled):
        """The shared axis of wells, labeled only on the panel that labels it.

        The panel that labels it is the first when the panels sit side by side and the last
        when they are stacked, so that the labels are beside or below every panel. Along `x`
        the labels are turned upright, as strain names do not fit across a well's width.

        """
        return (alt.X if self.vertical else alt.Y)(
            "unit_label:N",
            title=None,
            sort=order,
            axis=(
                alt.Axis(labelLimit=LABEL_LIMIT, labelAngle=-90 if self.vertical else 0)
                if labeled
                else None
            ),
        )

    def labels_wells(self, i_panel, n_panels):
        """Whether panel `i_panel` of `n_panels` is the one that labels the wells' axis."""
        return i_panel == (n_panels - 1 if self.vertical else 0)

    def panel_size(self, length):
        """A panel's size: a step per well, and `length` pixels for the values."""
        if self.vertical:
            return {"width": alt.Step(WELL_STEP), "height": length}
        return {"height": alt.Step(WELL_STEP), "width": length}

    def bar_band(self, fraction):
        """A bar's thickness across the wells' axis, as a fraction of a well's step."""
        return {("width" if self.vertical else "height"): {"band": fraction}}

    def configure(self, chart):
        """No grid, and the value axis' labels thinned out where they would collide.

        Only the value axis: thinning the wells' axis would drop wells' labels.

        """
        chart = chart.configure_axis(grid=False)
        if self.vertical:
            return chart.configure_axisY(labelOverlap="greedy")
        return chart.configure_axisX(labelOverlap="greedy")


def well_hover():
    """The selection of the well under the pointer, to outline it in every panel.

    Declare it on the concatenated chart, so it is in scope for all of its panels. It needs
    `empty=False` or every well matches when nothing is hovered, and is keyed on the axis
    label, which is one value per well.

    """
    return alt.selection_point(
        name="well_hover",
        fields=["unit_label"],
        on="pointerover",
        clear="pointerout",
        empty=False,
    )


def stroke_width(hover, hovered, ordinary):
    """Width of a well's outline, `hovered` on the well `hover` selects."""
    return alt.condition(hover, alt.value(hovered), alt.value(ordinary))


def carried_forward_stroke(has_final_library, legend_orient):
    """A well's outline color: whether it was carried forward, or black with no final library.

    One shared encoding, so one legend for the chart. It reads a ``carried_forward`` field
    labeled as in `CARRIED_FORWARD_LABELS`.

    """
    if not has_final_library:
        return alt.value(RETAINED_COLOR)
    return alt.Stroke(
        "carried_forward:N",
        title="in the library carried forward",
        scale=alt.Scale(
            domain=list(CARRIED_FORWARD_LABELS.values()),
            range=[RETAINED_COLOR, DROPPED_COLOR],
        ),
        legend=alt.Legend(orient=legend_orient, columns=1),
    )


def category_fill(field, title, categories, legend_orient, scheme=None, range_=None):
    """A stacked bar's fill, colored by `field` and drawn in the order `categories` gives.

    Args:
        `field` (str)
            The field to color by.
        `title` (str)
            The legend's title.
        `categories` (list)
            The category labels, in the order they stack and are listed.
        `legend_orient` (str)
            Which side the legend goes on.
        `scheme` (str or None)
            A `vega` color scheme to draw from.
        `range_` (list or None)
            Explicit colors, one per category, instead of a scheme.

    """
    if (scheme is None) == (range_ is None):
        raise ValueError("specify exactly one of `scheme` and `range_`")
    scale = (
        alt.Scale(scheme=scheme)
        if range_ is None
        else alt.Scale(domain=categories, range=range_)
    )
    return alt.Fill(
        f"{field}:N",
        title=title,
        sort=categories,
        scale=scale,
        legend=alt.Legend(labelLimit=LABEL_LIMIT, orient=legend_orient, columns=1),
    )


def order_fates(fates):
    """How a chart of read fates draws them, so a fate always looks the same.

    A bar stacks in the order its fates are colored in, which the label alone cannot convey
    to Vega, so the position in that order is returned to sort on.

    Args:
        `fates` (pandas.Series)
            The `fate` column of what the chart draws, not necessarily holding every fate.

    Returns:
        The 2-tuple `(fate_order, stack_order)`: the fates drawn in the order `FATE_COLORS`
        gives them, for the legend and the notes, and each row's position in that order.

    """
    unknown = sorted(set(fates) - set(FATE_COLORS))
    if unknown:
        raise ValueError(f"no color in `FATE_COLORS` for fate(s) {unknown}")
    fate_order = [fate for fate in FATE_COLORS if fate in set(fates)]
    return fate_order, fates.map({fate: i for (i, fate) in enumerate(fate_order)})


def fate_colors(fate_order):
    """Scale coloring each fate of `fate_order`, as `order_fates` returns it."""
    return alt.Scale(
        domain=fate_order, range=[FATE_COLORS[fate] for fate in fate_order]
    )


def fate_fill(fate_order, legend_orient):
    """Fill of a stacked bar of read fates, `fate_order` as `order_fates` returns it."""
    return category_fill(
        "fate",
        "fate of the read",
        fate_order,
        legend_orient,
        range_=[FATE_COLORS[fate] for fate in fate_order],
    )


def lookup_tooltip(lookup_frame, column):
    """Tooltip of a field a lookup adds from `lookup_frame`, as Altair cannot infer its type."""
    is_numeric = pandas.api.types.is_numeric_dtype(lookup_frame[column])
    return alt.Tooltip(column, type="quantitative" if is_numeric else "nominal")

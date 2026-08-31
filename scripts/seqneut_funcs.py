"""Functions shared by more than one of the pipeline's analysis scripts."""

import sys

import numpy
import pandas
from ruamel.yaml import YAML

# The QC-drop keys are names joined by spaces, so `ruamel` wraps them at its default
# width of 80, and a mapping key wrapped mid-string emits YAML that no longer parses.
# `width = None` would not help: `ruamel` falls back to 80 for a falsy width.
YAML_WIDTH = sys.maxsize

# Floats that came from arithmetic are rounded to this many significant figures before
# being embedded in a chart. A median over an even number of values otherwise runs to 18
# characters apiece, which is far finer than a pixel and finer than the chart tooltips.
SIG_FIGS_FOR_CHARTS = 4


def yaml_rt():
    """A round-trip `ruamel` YAML that never wraps a line."""
    yaml = YAML(typ="rt")
    yaml.width = YAML_WIDTH
    return yaml


def round_sig(x):
    """Round `x` to `SIG_FIGS_FOR_CHARTS` significant figures."""
    return float(f"{x:.{SIG_FIGS_FOR_CHARTS}g}")


def narrow_for_altair(df, drop=()):
    """Narrow a data frame that is about to be embedded in an `altair` chart.

    `altair` embeds the whole frame as inline JSON, so anything superfluous is paid for on
    every row. Drops `drop`, which should be the columns that are the same on every row:
    put them back with `transform_calculate` if the chart needs them, so that the value
    appears once in the chart rather than once per row. Also rounds the floats with
    `round_sig`. The frames written to the output CSVs are not passed through here, so
    this affects only the charts.

    """
    df = df.drop(columns=list(drop))
    return df.assign(
        **{
            col: df[col].map(round_sig)
            for col in df.columns
            if pandas.api.types.is_float_dtype(df[col])
        }
    )


def padded_log_domain(x, y):
    """Domain of a log axis spanning both `x` and `y`, padded away from the edges.

    Args:
        `x` and `y` (pandas.Series)
            Values plotted on the two axes, which share this domain so that the y = x
            line is a true diagonal.

    Returns:
        The rounded `[min, max]` domain.

    The padding is a fraction of the plotted range rather than a fixed factor, so that
    it does not swamp axes whose values span a narrow range, and the floor keeps values
    with almost no spread (such as a plate compared with an exact duplicate of itself)
    from collapsing. The end points are rounded for the same reason as the plotted
    values, and by so much less than the padding that the rounding cannot bring a point
    to the edge.

    """
    lo = min(x.min(), y.min())
    hi = max(x.max(), y.max())
    pad = max(0.05 * (numpy.log10(hi) - numpy.log10(lo)), numpy.log10(1.1))
    return [
        round_sig(10 ** (numpy.log10(lo) - pad)),
        round_sig(10 ** (numpy.log10(hi) + pad)),
    ]


def get_median_bound(s):
    """Get the bound on titer when taking median."""
    s = list(s)
    if len(s) % 2:
        return s[len(s) // 2]
    else:
        bounds = s[len(s) // 2 - 1 : len(s) // 2 + 1]
        assert len(bounds) == 2
        if len(set(bounds)) == 1:
            return bounds[0]
        elif "interpolated" in bounds:
            return next(b for b in bounds if b != "interpolated")
        else:
            return "inconsistent"


def viruses_in_plot_order(viruses, viral_strain_plot_order):
    """The unique `viruses`, ordered as configured or alphabetically.

    Args:
        `viruses` (iterable)
            Viruses to order, typically the `virus` column of a titer frame.
        `viral_strain_plot_order` (list or None)
            The configured order, or `None` to order alphabetically.

    Returns:
        List of the unique viruses in the order to plot them.

    Raises a `ValueError` naming any of `viruses` that `viral_strain_plot_order` lacks,
    as such a virus would be dropped from the order rather than positioned in it.

    """
    viruses = set(viruses)
    if viral_strain_plot_order is None:
        return sorted(viruses)
    if missing := viruses - set(viral_strain_plot_order):
        raise ValueError(
            f"`viral_strain_plot_order` lacks viruses with titers: {sorted(missing)}"
        )
    return [v for v in viral_strain_plot_order if v in viruses]


def pearson_r_log10(df, x_col="titer_x", y_col="titer_y"):
    """Pearson R of the log-transformed titers, or `None` if too few points.

    Fewer than three points always give a correlation of one or of nothing at all, so
    they are reported as not computed rather than as a correlation.

    """
    if len(df) < 3:
        return None
    return float(numpy.corrcoef(numpy.log10(df[x_col]), numpy.log10(df[y_col]))[0, 1])

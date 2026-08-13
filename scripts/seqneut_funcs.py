"""Functions shared by more than one of the pipeline's analysis scripts."""

import numpy

# Floats that came from arithmetic are rounded to this many significant figures before
# being embedded in a chart. A median over an even number of values otherwise runs to 18
# characters apiece, which is far finer than a pixel and finer than the chart tooltips.
SIG_FIGS_FOR_CHARTS = 4


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
        **{col: df[col].map(round_sig) for col in df.columns if df[col].dtype == float}
    )


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


def pearson_r_log10(df, x_col="titer_x", y_col="titer_y"):
    """Pearson R of the log-transformed titers, or `None` if too few points.

    Fewer than three points always give a correlation of one or of nothing at all, so
    they are reported as not computed rather than as a correlation.

    """
    if len(df) < 3:
        return None
    return float(numpy.corrcoef(numpy.log10(df[x_col]), numpy.log10(df[y_col]))[0, 1])

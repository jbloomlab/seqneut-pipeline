"""Correlate the titers on a plate with those on the other plates of its group."""

import sys

import altair as alt
import numpy
import pandas as pd
from seqneut_funcs import get_median_bound, pearson_r_log10, round_sig
from seqneut_report import Report

# `noqa: SIM115` as this log file must stay open for the life of the script
sys.stderr = sys.stdout = open(snakemake.log[0], "w")  # noqa: SIM115

_ = alt.data_transformers.disable_max_rows()

per_rep_titers_csvs = snakemake.input.per_rep_titers
corrs_csv = snakemake.output.corrs_csv
plate_date = snakemake.params.plate_date
comparator_dates = snakemake.params.comparator_dates
dilution_factor_or_concentration = snakemake.params.dilution_factor_or_concentration
concentration_units = snakemake.params.concentration_units
group = snakemake.wildcards.group
plate = snakemake.wildcards.plate

report = Report(
    title=f"Correlate the titers on {plate} with those on the other plates of {group}"
)

comparators_md = "\n".join(
    f" - `{p}` (dated {d})" for (p, d) in comparator_dates.items()
)
report.md(
    f"Correlating the titers on `{plate}` (dated {plate_date}) of `{group}` "
    f"with those on the following {len(comparator_dates)} plate(s), which "
    f"share at least one serum with it:\n\n{comparators_md}"
)

report.md("""
    ## Get the per-replicate titers for the sera on this plate
    Read the per-replicate titers for every serum measured on this plate; each of those
    files holds that serum's titers from all of the plates it was measured on, which is
    what the correlations below need.
    Each replicate has already passed the QC applied when its plate was processed, and
    records its plate and whether its virus was dropped by the per-serum QC on
    replicate-to-replicate variation:
    """)

if per_rep_titers_csvs:
    per_rep_titers = pd.concat(
        [pd.read_csv(f) for f in per_rep_titers_csvs], ignore_index=True
    )
else:
    # every serum on this plate lost all of its titers to the per-plate QC; the
    # columns are still specified so that the code below can be run unchanged
    per_rep_titers = pd.DataFrame(
        columns=[
            "group",
            "serum",
            "virus",
            "plate",
            "replicate",
            "titer",
            "titer_bound",
            "titer_as",
            "titer_units",
            "dropped_by_qc",
        ]
    )

assert (per_rep_titers["group"] == group).all(), (
    f"not all of the input titers are for `{group=}`:\n"
    f"{per_rep_titers[per_rep_titers['group'] != group]}"
)

# a serum has the same 'titer_as' and 'titer_units' on every plate, as these are
# configured per serum rather than per plate; the titers are therefore comparable
# between plates even if they differ between sera in the group
titer_types = per_rep_titers.groupby("serum")[["titer_as", "titer_units"]].nunique()
assert (
    (titer_types == 1).all().all()
), f"serum with more than one type of titer:\n{titer_types[titer_types.gt(1).any(axis=1)]}"

# a serum measured on this plate can only also appear on a plate that shares a serum
# with this one, which is exactly what makes a plate a comparator
if extra_plates := set(per_rep_titers["plate"]) - {plate, *comparator_dates}:
    raise ValueError(
        f"titers from plates that are neither `{plate}` nor one of its "
        f"comparators: {extra_plates}"
    )

report.md(
    f"Read {len(per_rep_titers)} per-replicate titers for "
    f"{per_rep_titers['serum'].nunique()} sera measured on `{plate}`"
)

report.md("""
    ## Get the titers for each strain on each plate
    A plate can have several barcodes for a strain, so the titer for a strain against a
    serum on a plate is the median over that plate's barcodes.
    We use the same handling of the bound on a median titer as when medians are taken
    across all replicates of a serum:
    """)

# the per-serum QC drops a virus rather than an individual replicate, so its outcome
# is a property of a serum-virus pair and can be carried through the median below
qc_outcomes = per_rep_titers.groupby(["serum", "virus"])["dropped_by_qc"].nunique()
assert (
    qc_outcomes == 1
).all(), f"serum-virus with inconsistent QC outcomes:\n{qc_outcomes[qc_outcomes != 1]}"

plate_titers = (
    per_rep_titers.sort_values("titer")  # for getting median bound
    .groupby(["serum", "virus", "plate"], as_index=False)
    .aggregate(
        titer=pd.NamedAgg("titer", "median"),
        n_barcodes=pd.NamedAgg("replicate", "count"),
        titer_bound=pd.NamedAgg("titer_bound", get_median_bound),
        dropped_by_qc=pd.NamedAgg("dropped_by_qc", "first"),
    )
)
report.md(
    f"Got {len(plate_titers)} serum-strain titers across the plates"
    + (
        f", {plate_titers['n_barcodes'].max()} barcodes at most per titer"
        if len(plate_titers)
        else ""
    )
)

report.md("""
    ## Correlate the titers on this plate with those on each other plate
    For each of the other plates, get the titers for every serum measured on both it and
    this plate, and the Pearson correlation of the log titers for each of those sera.
    These per-serum correlations are the ones written to the CSV file at the bottom of
    this report; those annotating the plots are computed from the plotted points instead,
    so that they follow the dropdowns:
    """)

# This plate's titers get the '_y' suffix and the other plate's the '_x' suffix, as
# this plate is on the y axis of every panel below. 'dropped_by_qc' is a merge key
# rather than a merged column, because the per-serum QC drops a virus for a serum on
# every plate at once and so it is the same on both plates.
merge_keys = ["serum", "virus", "dropped_by_qc"]
pair_titers = (
    plate_titers[plate_titers["plate"] == plate]
    .drop(columns="plate")
    .merge(
        plate_titers[plate_titers["plate"].isin(comparator_dates)].rename(
            columns={"plate": "other_plate"}
        ),
        on=merge_keys,
        suffixes=("_y", "_x"),
        validate="one_to_many",
    )
    .assign(
        censored=lambda x: (x["titer_bound_y"] != "interpolated")
        | (x["titer_bound_x"] != "interpolated"),
    )
)

comparator_info = []
for other in comparator_dates:  # config order, which orders the panels
    other_df = pair_titers[pair_titers["other_plate"] == other]
    if not len(other_df):
        continue
    comparator_info.append(
        {
            "other_plate": other,
            "other_date": comparator_dates[other],
            "n": len(other_df),
            # The correlations shown on the plots are not computed here. They are
            # computed by `altair` transforms from whatever the dropdowns are
            # currently showing, so a number computed here would go stale as soon as
            # a dropdown was changed. The per-serum correlations below are computed
            # here because they go to a CSV file rather than onto a plot.
            # Each panel gets just its own plate's titers, rather than every panel
            # getting the whole frame and filtering to its plate client-side. That
            # keeps the plotted data frame narrow: the other plate is named in the
            # panel title, so an `other_plate` column would be the same string
            # repeated on every row. It costs no extra rows, since `altair`
            # consolidates every frame into the chart's `datasets` and the rows are
            # partitioned among the panels either way. This plate is likewise named
            # in the title of the plot rather than in a column, and the QC status is
            # kept as a boolean that is expanded to a label when plotting.
            "titers": (
                other_df[
                    [
                        "serum",
                        "virus",
                        "titer_x",
                        "titer_y",
                        "n_barcodes_x",
                        "n_barcodes_y",
                        "censored",
                        "dropped_by_qc",
                    ]
                ]
                .assign(
                    titer_x=lambda x: x["titer_x"].map(round_sig),
                    titer_y=lambda x: x["titer_y"].map(round_sig),
                )
                .reset_index(drop=True)
            ),
            # the same correlations for each serum on its own, for the CSV file
            "per_serum": [
                {
                    "serum": serum,
                    "n_strains": len(serum_df),
                    "pearson_r": pearson_r_log10(serum_df),
                    "n_strains_retained": int((~serum_df["dropped_by_qc"]).sum()),
                    "pearson_r_retained": pearson_r_log10(
                        serum_df[~serum_df["dropped_by_qc"]]
                    ),
                }
                for (serum, serum_df) in other_df.groupby("serum", sort=True)
            ],
        }
    )

if comparator_info:
    # Every panel shares one domain, used for both of its axes. Sharing it between
    # the axes makes the y = x line a true diagonal, and sharing it between the
    # panels puts all of the plates on the same scale. The domain is the range of the
    # titers on this plate or any of the others, padded by a fraction of that range
    # on a log scale so that points do not sit on the plot edge. The padding is
    # proportional rather than a fixed factor so that it does not swamp a plate whose
    # titers span a narrow range, and the floor keeps a plate with almost no spread
    # (such as one compared with an exact duplicate of itself) from collapsing.
    lo = min(pair_titers["titer_x"].min(), pair_titers["titer_y"].min())
    hi = max(pair_titers["titer_x"].max(), pair_titers["titer_y"].max())
    pad = max(0.05 * (numpy.log10(hi) - numpy.log10(lo)), numpy.log10(1.1))
    # rounded for the same reason as the titers, and by so much less than the padding
    # that it cannot bring a point to the edge of a panel
    plot_domain = [
        round_sig(10 ** (numpy.log10(lo) - pad)),
        round_sig(10 ** (numpy.log10(hi) + pad)),
    ]
    # end points of the y = x line of each panel, which are the ends of its axes
    # rather than of its data so that the line spans the whole panel
    diagonal = pd.DataFrame({"titer": plot_domain})
    plot_sera = sorted(set(pair_titers["serum"]))
    report.md(
        f"Comparing `{plate}` with {len(comparator_info)} other plate(s), "
        f"covering {sum(info['n'] for info in comparator_info)} "
        "plate-serum-strain titer comparison(s)"
    )
else:
    plot_domain = None
    diagonal = None
    plot_sera = []
    report.md(
        f"No other plate has a serum-strain titer that is also measured on "
        f"`{plate}`, so there is no plate-to-plate correlation to show. Note "
        "that the plates listed above share sera with this one in the "
        "configuration, so this means the QC has dropped all of the titers they "
        "had in common."
    )

uncompared_sera = sorted(set(plate_titers["serum"]) - set(pair_titers["serum"]))
if uncompared_sera:
    report.md(
        f"The following {len(uncompared_sera)} sera on `{plate}` have no titers "
        f"in common with another plate and so are not shown: {uncompared_sera}"
    )

report.md("""
    ## Plot the plate-to-plate correlations
    Each point is one strain against one serum, at its titer on this plate (y axis)
    versus on the other plate (x axis), with the dashed `y = x` line for reference.
    Each panel is annotated with the number of titers, the number of sera, and the
    Pearson correlation of the log titers.

    Points are colored by whether they are retained in the final titers or were dropped
    by the per-serum QC on replicate-to-replicate variation.
    That QC drops a strain for a serum on all of its plates at once, so the color is the
    same in every panel, and a dropped point is not necessarily one that this panel's two
    plates disagreed about.
    Whether they disagreed is the distance from the `y = x` line, which the tooltip gives
    as the fold difference between the two plates.

    Two dropdowns below the plot restrict what is shown: `show strain-serum pairs by QC`
    to only the retained or only the dropped points, and `show serum` to a single serum.
    The annotated numbers are computed from whatever is currently shown, so the
    correlation over just the retained titers, or for one serum, comes from setting a
    dropdown.
    Mousing over a point highlights that strain-serum combination in every panel.

    The per-serum correlations are also in the CSV file at the bottom of this report:
    """)

if comparator_info:
    titer_label = (
        "titer"
        if dilution_factor_or_concentration == "dilution_factor"
        else f"titer ({concentration_units})"
    )

    qc_statuses = ["retained in final titers", "dropped by per-serum QC"]
    # The QC status is carried in the data as the boolean `dropped_by_qc` and expanded
    # to these labels here, so that the labels are not repeated on every row of the
    # data embedded in the chart. It is a property of a strain-serum combination, not
    # of a single titer or of a strain on its own: the per-serum QC drops a strain for
    # a serum on all of its plates at once, and a strain can be dropped for one serum
    # and retained for another. Each point is one strain-serum combination, so the
    # status is per point, and it is the same in every panel.
    qc_calculate = f"datum.dropped_by_qc ? '{qc_statuses[1]}' : '{qc_statuses[0]}'"

    # How much this panel's two plates disagree about a point, which is what the QC
    # status above does not say: the QC drops a strain-serum for varying too much over
    # all of its replicates, so a point can be dropped because of a plate other than
    # the one in the panel, and then sit right on the diagonal here. Computed as a
    # `vega` expression rather than a column so that it costs nothing in the data
    # embedded in the chart.
    fold_diff_calculate = (
        "datum.titer_y > datum.titer_x"
        " ? datum.titer_y / datum.titer_x"
        " : datum.titer_x / datum.titer_y"
    )

    # Dropdowns are used to select what to show rather than binding a selection to a
    # legend. Binding to the legend does not work here: the legend of a concatenated
    # chart belongs to whichever subplot draws it, so the binding never reaches the
    # other subplots. An input binding does work, because the param is declared on the
    # concatenated chart itself (see `add_params` below) and so is in scope for every
    # subplot. Selecting "all" in a dropdown shows everything.
    qc_selection = alt.selection_point(
        fields=["qc_status"],
        bind=alt.binding_select(
            options=[None, *qc_statuses],
            labels=["all", *qc_statuses],
            name="show strain-serum pairs by QC",
        ),
        name="qc_status_selection",
    )
    serum_selection = alt.selection_point(
        fields=["serum"],
        bind=alt.binding_select(
            options=[None, *plot_sera],
            labels=["all", *plot_sera],
            name="show serum",
        ),
        name="serum_selection",
    )

    # Mousing over a point highlights it in every panel, not just the panel being
    # moused over, because this param is also declared on the concatenated chart. The
    # selection is on the strain and the serum together rather than on the strain
    # alone, so that at most one point per panel is highlighted; selecting on `virus`
    # alone would highlight that strain for every serum.
    point_hover = alt.selection_point(
        fields=["serum", "virus"],
        on="mouseover",
        empty=False,
        name="strain_serum_hover",
    )

    # The lines of the label annotating each panel, in the order they are drawn, and
    # computed by `altair` transforms rather than in `pandas` so that they describe just
    # the titers the dropdowns are currently showing. Each is the last step of the
    # pipeline set up in `stats_base` below, which filters exactly as the points do and
    # then aggregates what is left, so the numbers cannot disagree with the points on
    # the panel. The Pearson R is of the log-transformed titers, as it is in the CSV
    # file. `sxx` or `syy` of zero means every remaining titer on one of the plates is
    # identical, which would make R a division by zero rather than a number.
    labels = [
        (
            "datum.n + (datum.n == 1 ? ' titer, ' : ' titers, ')"
            " + datum.n_sera + (datum.n_sera == 1 ? ' serum' : ' sera')"
        ),
        (
            "datum.n < 3 ? 'R = not computed (< 3 titers)' : "
            "(datum.sxx <= 0 || datum.syy <= 0"
            " ? 'R = not computed (no spread in titers)'"
            " : 'R = ' + format(datum.sxy / sqrt(datum.sxx * datum.syy), '.2f'))"
        ),
    ]

    scale = alt.Scale(type="log", nice=False, domain=plot_domain)
    n_columns = min(4, len(comparator_info))

    charts = []
    for i, info in enumerate(comparator_info):
        other = info["other_plate"]
        # Every panel has the same y axis, this plate, so it is drawn only on the
        # leftmost panel of each row. Each panel has its own x axis, the plate it
        # compares with, so that one is titled on every panel.
        leftmost = (i % n_columns) == 0
        x_title = f"{titer_label} on {other}"
        y_title = f"{titer_label} on {plate}" if leftmost else None
        y_axis = alt.Axis(labels=leftmost, ticks=leftmost)
        # Shared by the points and by the labels annotating them, so that both show
        # the same titers: the QC status has to be calculated before the dropdown that
        # selects on it can filter, and the labels have to be filtered exactly as the
        # points are. `altair` copies a chart when a method is chained onto it, so
        # this can be built on without the layers affecting each other.
        base = (
            alt.Chart(info["titers"])
            .transform_calculate(qc_status=qc_calculate)
            .transform_filter(qc_selection)
            .transform_filter(serum_selection)
        )
        points = (
            base.transform_calculate(fold_difference=fold_diff_calculate)
            .encode(
                alt.X("titer_x:Q", title=x_title, scale=scale),
                alt.Y("titer_y:Q", title=y_title, scale=scale, axis=y_axis),
                alt.Color(
                    "qc_status:N",
                    scale=alt.Scale(
                        domain=qc_statuses,
                        range=["MediumBlue", "OrangeRed"],
                    ),
                    legend=(
                        alt.Legend(
                            title=[
                                "strain-serum QC status",
                                "(use the 'show strain-serum",
                                "pairs by QC' dropdown below",
                                "the plot to show only one",
                                "of these)",
                            ],
                            titleLimit=400,
                            symbolLimit=0,
                        )
                        if i == 0
                        else None
                    ),
                ),
                alt.Shape(
                    "censored:N",
                    title=["is either median titer", "censored at a bound?"],
                    legend=alt.Legend(titleLimit=400) if i == 0 else None,
                ),
                tooltip=[
                    alt.Tooltip("serum:N", title="serum"),
                    alt.Tooltip("virus:N", title="strain"),
                    alt.Tooltip(
                        "titer_y:Q",
                        title=f"{titer_label} on {plate}",
                        format=".3g",
                    ),
                    alt.Tooltip(
                        "titer_x:Q",
                        title=f"{titer_label} on {other}",
                        format=".3g",
                    ),
                    alt.Tooltip(
                        "fold_difference:Q",
                        title="fold difference between these plates",
                        format=".2f",
                    ),
                    alt.Tooltip("n_barcodes_y:Q", title=f"barcodes on {plate}"),
                    alt.Tooltip("n_barcodes_x:Q", title=f"barcodes on {other}"),
                    alt.Tooltip("qc_status:N", title="QC status"),
                ],
                # the moused-over strain gets a larger point with a thick red stroke
                size=alt.condition(point_hover, alt.value(180), alt.value(60)),
                stroke=alt.condition(point_hover, alt.value("red"), alt.value("black")),
                strokeWidth=alt.condition(point_hover, alt.value(3), alt.value(0.5)),
                fillOpacity=alt.condition(point_hover, alt.value(1), alt.value(0.5)),
            )
            .mark_point(filled=True, strokeOpacity=1)
        )
        # y = x, drawn from the ends of the axes rather than from the titers so that
        # the line spans the whole panel. This layer repeats the axis encodings of the
        # points layer rather than leaving them unset, because the layers share their
        # axes and a title of `None` here blanks the shared title.
        diagonal_layer = (
            alt.Chart(diagonal)
            .encode(
                alt.X("titer:Q", title=x_title, scale=scale),
                alt.Y("titer:Q", title=y_title, scale=scale, axis=y_axis),
            )
            .mark_line(color="gray", strokeDash=[4, 4], strokeWidth=1)
        )
        # The Pearson R of the log titers and the counts of what is being correlated,
        # aggregated from the titers left by the filters in `base` and so recomputed
        # by `vega` whenever a dropdown changes what is shown. They are drawn in the
        # panel rather than in its title because a title is a fixed string in the chart
        # specification: it can refer to a param, but not to a value aggregated from
        # the data. Each line is its own layer, since a text mark draws one string and
        # the aggregate of a layer cannot be shared with another layer.
        stats_base = (
            base.transform_calculate(
                log_x="log(datum.titer_x) / LN10",
                log_y="log(datum.titer_y) / LN10",
            )
            .transform_joinaggregate(mean_x="mean(log_x)", mean_y="mean(log_y)")
            .transform_calculate(
                dxdy="(datum.log_x - datum.mean_x) * (datum.log_y - datum.mean_y)",
                dx2="pow(datum.log_x - datum.mean_x, 2)",
                dy2="pow(datum.log_y - datum.mean_y, 2)",
            )
            .transform_aggregate(
                sxy="sum(dxdy)",
                sxx="sum(dx2)",
                syy="sum(dy2)",
                n="count()",
                n_sera="distinct(serum)",
            )
        )
        stats = [
            stats_base.transform_calculate(label=label)
            .encode(
                x=alt.value(4),
                y=alt.value(4 + 12 * line),
                text=alt.Text("label:N"),
            )
            .mark_text(align="left", baseline="top", color="black", fontSize=10)
            for (line, label) in enumerate(labels)
        ]
        charts.append(
            alt.layer(diagonal_layer, points, *stats).properties(
                width=200,
                height=200,
                title=alt.Title(
                    f"{other} ({info['other_date']})",
                    fontSize=13,
                ),
            )
        )
    report.chart(
        alt.concat(*charts, columns=n_columns)
        # all of the params are declared here rather than on the panels so that they
        # are in scope for all of them
        .add_params(qc_selection, serum_selection, point_hover).properties(
            title=alt.Title(
                f"{plate} (dated {plate_date}) versus the other plates of {group}",
                subtitle=[
                    (
                        "each point is a viral strain for one serum; titers are "
                        "medians over each plate's barcodes"
                    ),
                    (
                        f"the y axis of every panel is the {titer_label} on "
                        f"{plate}, and the x axis that on the plate titling the panel"
                    ),
                ],
                fontSize=15,
                anchor="middle",
            )
        )
        # `titleFontWeight` is set because axis titles are bold by default, and here
        # only the overall and panel titles should be bold. `labelOverlap` is set here
        # rather than on each axis so that the chart spec carries it once instead of
        # once per axis of every panel.
        .configure_axis(
            grid=False,
            labelOverlap=True,
            titleFontSize=13,
            titleFontWeight="normal",
        )
        # `limit=0` means no limit, so that neither the titles nor the subtitles of
        # the plot or its panels are truncated. It is set in the chart's config rather
        # than on each title so that the spec carries it once.
        .configure_title(limit=0)
    )

report.md("""
    ## Write the correlations to a CSV
    Write the Pearson correlations for each other plate and serum to a CSV file.
    Note that each pair of plates appears in both plates' files, with the two plates
    swapped:
    """)

# the columns are specified so that the file still has a header if there are no
# correlations
corrs = pd.DataFrame(
    [
        {
            "plate": plate,
            "date": plate_date,
            "other_plate": info["other_plate"],
            "other_date": info["other_date"],
            **serum_corrs,
        }
        for info in comparator_info
        for serum_corrs in info["per_serum"]
    ],
    columns=[
        "plate",
        "date",
        "other_plate",
        "other_date",
        "serum",
        "n_strains",
        "pearson_r",
        "n_strains_retained",
        "pearson_r_retained",
    ],
)
report.md(f"Writing the correlations to `{corrs_csv}`")
print(f"Writing the correlations to {corrs_csv}")
corrs.to_csv(corrs_csv, index=False, float_format="%.4g")
report.table(corrs, index=False)

report.write(snakemake.output.html)
print(f"Wrote the report to {snakemake.output.html}")

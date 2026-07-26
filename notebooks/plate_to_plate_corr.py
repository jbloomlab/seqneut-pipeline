# /// script
# [tool.marimo.runtime]
# auto_instantiate = false
# output_max_bytes = 50_000_000
# ///

import marimo

__generated_with = "0.17.6"
app = marimo.App(width="full")


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Plate-to-plate correlations of titers for a plate
    Compare the titers measured on this plate with those measured for the same sera on the other plates of its group, to assess plate-to-plate reproducibility.
    """)
    return


@app.cell
def _():
    import altair as alt
    import marimo as mo
    import numpy
    import pandas as pd

    _ = alt.data_transformers.disable_max_rows()
    return alt, mo, numpy, pd


@app.cell
def _():
    # Load context from pickled file.
    #
    # This cell supports multiple ways to provide context:
    # 1. Via command-line: marimo export html notebook.py -- --context-pickle path/to/context.pickle
    # 2. Via saved pickle: Manually save a context pickle to test_example/results/context_dev.pickle
    # 3. Stub context: If no pickle available, creates minimal empty context for exploration
    #
    # For interactive development with `marimo edit`, you can:
    # - Run the pipeline once to generate a real context pickle, then copy it to context_dev.pickle
    # - Or work with the stub context (downstream cells will show warnings/empty data)

    import argparse
    import os
    import pathlib
    import pickle
    import sys

    # Check if context-pickle argument is provided (run by driver script)
    from_cmdline = "--context-pickle" in sys.argv

    if from_cmdline:
        # Running via driver script - parse args
        print("Loading context from command-line argument")
        p = argparse.ArgumentParser()
        p.add_argument("--context-pickle", required=True)
        args = p.parse_args()
        context_pickle_path = pathlib.Path(args.context_pickle)
    else:
        # Running in marimo edit - try to use development pickle
        print("Running in marimo edit mode")
        # set `context_pickle_path` to valid option if using edit mode
        context_pickle_path = None
        # context_pickle_path = pathlib.Path("test_example/results/plate_to_plate_corrs/plate_to_plate_corr_serum_plate2_context.pickle")

    # Load context if pickle path exists and is valid
    if context_pickle_path and context_pickle_path.exists():
        print(f"Reading context from {context_pickle_path}")
        with open(context_pickle_path, "rb") as f_context:
            context = pickle.load(f_context)

        # Handle working directory
        context_workdir = context["workdir"]
        current_workdir = os.getcwd()

        if from_cmdline:
            # Running via snakemake - verify workdir matches
            if context_workdir != current_workdir:
                raise RuntimeError(
                    f"Context workdir mismatch!\n"
                    f"  Context was created in: {context_workdir}\n"
                    f"  Currently running in:   {current_workdir}\n"
                    f"This should not happen when running via Snakemake."
                )
            print(f"Verified working directory: {current_workdir}")
        else:
            # Running in marimo edit - change to context workdir
            if context_workdir and context_workdir != current_workdir:
                print(f"Changing directory from {current_workdir} to {context_workdir}")
                os.chdir(context_workdir)
            elif context_workdir:
                print(f"Already in correct working directory: {context_workdir}")
    else:
        # Create a minimal stub context for interactive development
        print("Creating minimal stub context that you need to complete")
        context = {
            "input": {},
            "output": {},
            "params": {},
            "wildcards": {},
            "threads": 1,
            "resources": {},
        }
    return (context,)


@app.cell
def _(context, mo):
    # Extract variables from context - raises KeyError if required keys missing
    per_rep_titers_csvs = context["input"]["per_rep_titers"]
    corrs_csv = context["output"]["corrs_csv"]
    plate_date = context["params"]["plate_date"]
    comparator_dates = context["params"]["comparator_dates"]
    dilution_factor_or_concentration = context["params"][
        "dilution_factor_or_concentration"
    ]
    concentration_units = context["params"]["concentration_units"]
    group = context["wildcards"]["group"]
    plate = context["wildcards"]["plate"]

    # Show informative message about context mode
    if not context["input"]:
        mo.output.append(
            mo.callout(
                mo.md(
                    "**⚠️ Running in interactive mode with stub context**\n\n"
                    "To run with real data:\n"
                    "1. Run the pipeline to generate a context pickle\n"
                    "2. Copy it to `test_example/results/context_dev.pickle`\n"
                    "3. Or run: `marimo export html notebook.py -- --context-pickle path/to/context.pickle`"
                ),
                kind="warn",
            )
        )

    _comparators_md = "\n".join(
        f" - `{p}` (dated {d})" for (p, d) in comparator_dates.items()
    )
    mo.output.append(
        mo.md(
            f"Correlating the titers on `{plate}` (dated {plate_date}) of `{group}` "
            f"with those on the following {len(comparator_dates)} plate(s), which "
            f"share at least one serum with it:\n\n{_comparators_md}"
        )
    )
    return (
        comparator_dates,
        concentration_units,
        corrs_csv,
        dilution_factor_or_concentration,
        group,
        per_rep_titers_csvs,
        plate,
        plate_date,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Get the per-replicate titers for the sera on this plate
    Read the per-replicate titers for every serum measured on this plate; each of those
    files holds that serum's titers from all of the plates it was measured on, which is
    what the correlations below need.
    Each replicate has already passed the QC applied when its plate was processed, and
    records its plate and whether its virus was dropped by the per-serum QC on
    replicate-to-replicate variation:
    """)
    return


@app.cell
def _(comparator_dates, group, mo, pd, per_rep_titers_csvs, plate):
    if per_rep_titers_csvs:
        per_rep_titers = pd.concat(
            [pd.read_csv(f) for f in per_rep_titers_csvs], ignore_index=True
        )
    else:
        # every serum on this plate lost all of its titers to the per-plate QC; the
        # columns are still specified so that the cells below can be run unchanged
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
    _titer_types = per_rep_titers.groupby("serum")[
        ["titer_as", "titer_units"]
    ].nunique()
    assert (
        (_titer_types == 1).all().all()
    ), f"serum with more than one type of titer:\n{_titer_types[_titer_types.gt(1).any(axis=1)]}"

    # a serum measured on this plate can only also appear on a plate that shares a serum
    # with this one, which is exactly what makes a plate a comparator
    if _extra_plates := set(per_rep_titers["plate"]) - {plate, *comparator_dates}:
        raise ValueError(
            f"titers from plates that are neither `{plate}` nor one of its "
            f"comparators: {_extra_plates}"
        )

    mo.output.append(
        mo.md(
            f"Read {len(per_rep_titers)} per-replicate titers for "
            f"{per_rep_titers['serum'].nunique()} sera measured on `{plate}`"
        )
    )
    return (per_rep_titers,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Get the titers for each strain on each plate
    A plate can have several barcodes for a strain, so the titer for a strain against a
    serum on a plate is the median over that plate's barcodes.
    We use the same handling of the bound on a median titer as when medians are taken
    across all replicates of a serum:
    """)
    return


@app.cell
def _():
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

    return (get_median_bound,)


@app.cell
def _(get_median_bound, mo, pd, per_rep_titers):
    # the per-serum QC drops a virus rather than an individual replicate, so its outcome
    # is a property of a serum-virus pair and can be carried through the median below
    _qc_outcomes = per_rep_titers.groupby(["serum", "virus"])["dropped_by_qc"].nunique()
    assert (
        _qc_outcomes == 1
    ).all(), (
        f"serum-virus with inconsistent QC outcomes:\n{_qc_outcomes[_qc_outcomes != 1]}"
    )

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
    mo.output.append(
        mo.md(
            f"Got {len(plate_titers)} serum-strain titers across the plates"
            + (
                f", {plate_titers['n_barcodes'].max()} barcodes at most per titer"
                if len(plate_titers)
                else ""
            )
        )
    )
    return (plate_titers,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Correlate the titers on this plate with those on each other plate
    For each of the other plates, get the titers for every serum measured on both it and
    this plate, and the Pearson correlation of the log titers for each of those sera.
    These per-serum correlations are the ones written to the CSV file at the bottom of this
    notebook; those annotating the plots are computed from the plotted points instead, so
    that they follow the dropdowns:
    """)
    return


@app.cell
def _(comparator_dates, mo, numpy, pd, plate, plate_titers):
    def _pearson_r(df):
        """Pearson R of the log-transformed titers, or `None` if too few points."""
        if len(df) < 3:
            return None
        return float(
            numpy.corrcoef(numpy.log10(df["titer_x"]), numpy.log10(df["titer_y"]))[0, 1]
        )

    def _round_sig(x):
        """Round to four significant figures.

        The titers are medians, so an even number of barcodes gives values such as
        210.35000000000002 that cost 18 characters apiece in the JSON embedded in the
        chart. Four significant figures is far finer than a pixel on the log axes and
        finer than the `.3g` of the tooltips, and the correlations are computed from the
        unrounded titers, so this only removes wasted bytes.
        """
        return float(f"{x:.4g}")

    # This plate's titers get the '_y' suffix and the other plate's the '_x' suffix, as
    # this plate is on the y axis of every panel below. 'dropped_by_qc' is a merge key
    # rather than a merged column, because the per-serum QC drops a virus for a serum on
    # every plate at once and so it is the same on both plates.
    _merge_keys = ["serum", "virus", "dropped_by_qc"]
    pair_titers = (
        plate_titers[plate_titers["plate"] == plate]
        .drop(columns="plate")
        .merge(
            plate_titers[plate_titers["plate"].isin(comparator_dates)].rename(
                columns={"plate": "other_plate"}
            ),
            on=_merge_keys,
            suffixes=("_y", "_x"),
            validate="one_to_many",
        )
        .assign(
            censored=lambda x: (x["titer_bound_y"] != "interpolated")
            | (x["titer_bound_x"] != "interpolated"),
        )
    )

    comparator_info = []
    for _other in comparator_dates:  # config order, which orders the panels
        _other_df = pair_titers[pair_titers["other_plate"] == _other]
        if not len(_other_df):
            continue
        comparator_info.append(
            {
                "other_plate": _other,
                "other_date": comparator_dates[_other],
                "n": len(_other_df),
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
                    _other_df[
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
                        titer_x=lambda x: x["titer_x"].map(_round_sig),
                        titer_y=lambda x: x["titer_y"].map(_round_sig),
                    )
                    .reset_index(drop=True)
                ),
                # the same correlations for each serum on its own, for the CSV file
                "per_serum": [
                    {
                        "serum": _serum,
                        "n_strains": len(_serum_df),
                        "pearson_r": _pearson_r(_serum_df),
                        "n_strains_retained": int((~_serum_df["dropped_by_qc"]).sum()),
                        "pearson_r_retained": _pearson_r(
                            _serum_df[~_serum_df["dropped_by_qc"]]
                        ),
                    }
                    for (_serum, _serum_df) in _other_df.groupby("serum", sort=True)
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
        _lo = min(pair_titers["titer_x"].min(), pair_titers["titer_y"].min())
        _hi = max(pair_titers["titer_x"].max(), pair_titers["titer_y"].max())
        _pad = max(0.05 * (numpy.log10(_hi) - numpy.log10(_lo)), numpy.log10(1.1))
        # rounded for the same reason as the titers, and by so much less than the padding
        # that it cannot bring a point to the edge of a panel
        plot_domain = [
            _round_sig(10 ** (numpy.log10(_lo) - _pad)),
            _round_sig(10 ** (numpy.log10(_hi) + _pad)),
        ]
        # end points of the y = x line of each panel, which are the ends of its axes
        # rather than of its data so that the line spans the whole panel
        diagonal = pd.DataFrame({"titer": plot_domain})
        plot_sera = sorted(set(pair_titers["serum"]))
        mo.output.append(
            mo.md(
                f"Comparing `{plate}` with {len(comparator_info)} other plate(s), "
                f"covering {sum(info['n'] for info in comparator_info)} "
                "plate-serum-strain titer comparison(s)"
            )
        )
    else:
        plot_domain = None
        diagonal = None
        plot_sera = []
        mo.output.append(
            mo.md(
                f"No other plate has a serum-strain titer that is also measured on "
                f"`{plate}`, so there is no plate-to-plate correlation to show. Note "
                "that the plates listed above share sera with this one in the "
                "configuration, so this means the QC has dropped all of the titers they "
                "had in common."
            )
        )

    _uncompared_sera = sorted(
        set(plate_titers["serum"]) - set(pair_titers["serum"]),
    )
    if _uncompared_sera:
        mo.output.append(
            mo.md(
                f"The following {len(_uncompared_sera)} sera on `{plate}` have no titers "
                f"in common with another plate and so are not shown: {_uncompared_sera}"
            )
        )
    return comparator_info, diagonal, plot_domain, plot_sera


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Plot the plate-to-plate correlations
    Each point is one strain against one serum, at its titer on this plate (y axis) versus
    on the other plate (x axis), with the dashed `y = x` line for reference.
    Each panel is annotated with the number of titers, the number of sera, and the Pearson
    correlation of the log titers.

    Points are colored by whether they are retained in the final titers or were dropped by
    the per-serum QC on replicate-to-replicate variation.
    That QC drops a strain for a serum on all of its plates at once, so the color is the
    same in every panel, and a dropped point is not necessarily one that this panel's two
    plates disagreed about.
    Whether they disagreed is the distance from the `y = x` line, which the tooltip gives
    as the fold difference between the two plates.

    Two dropdowns below the plot restrict what is shown: `show strain-serum pairs by QC`
    to only the retained or only the dropped points, and `show serum` to a single serum.
    The annotated numbers are computed from whatever is currently shown, so the correlation
    over just the retained titers, or for one serum, comes from setting a dropdown.
    Mousing over a point highlights that strain-serum combination in every panel.

    The per-serum correlations are also in the CSV file at the bottom of this notebook:
    """)
    return


@app.cell
def _(
    alt,
    comparator_info,
    concentration_units,
    diagonal,
    dilution_factor_or_concentration,
    group,
    mo,
    plate,
    plate_date,
    plot_domain,
    plot_sera,
):
    if comparator_info:
        _titer_label = (
            "titer"
            if dilution_factor_or_concentration == "dilution_factor"
            else f"titer ({concentration_units})"
        )

        _qc_statuses = ["retained in final titers", "dropped by per-serum QC"]
        # The QC status is carried in the data as the boolean `dropped_by_qc` and expanded
        # to these labels here, so that the labels are not repeated on every row of the
        # data embedded in the chart. It is a property of a strain-serum combination, not
        # of a single titer or of a strain on its own: the per-serum QC drops a strain for
        # a serum on all of its plates at once, and a strain can be dropped for one serum
        # and retained for another. Each point is one strain-serum combination, so the
        # status is per point, and it is the same in every panel.
        _qc_calculate = (
            f"datum.dropped_by_qc ? '{_qc_statuses[1]}' : '{_qc_statuses[0]}'"
        )

        # How much this panel's two plates disagree about a point, which is what the QC
        # status above does not say: the QC drops a strain-serum for varying too much over
        # all of its replicates, so a point can be dropped because of a plate other than
        # the one in the panel, and then sit right on the diagonal here. Computed as a
        # `vega` expression rather than a column so that it costs nothing in the data
        # embedded in the chart.
        _fold_diff_calculate = (
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
        _qc_selection = alt.selection_point(
            fields=["qc_status"],
            bind=alt.binding_select(
                options=[None, *_qc_statuses],
                labels=["all", *_qc_statuses],
                name="show strain-serum pairs by QC",
            ),
            name="qc_status_selection",
        )
        _serum_selection = alt.selection_point(
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
        _point_hover = alt.selection_point(
            fields=["serum", "virus"],
            on="mouseover",
            empty=False,
            name="strain_serum_hover",
        )

        # The lines of the label annotating each panel, in the order they are drawn, and
        # computed by `altair` transforms rather than in `pandas` so that they describe just
        # the titers the dropdowns are currently showing. Each is the last step of the
        # pipeline set up in `_stats_base` below, which filters exactly as the points do and
        # then aggregates what is left, so the numbers cannot disagree with the points on
        # the panel. The Pearson R is of the log-transformed titers, as it is in the CSV
        # file. `sxx` or `syy` of zero means every remaining titer on one of the plates is
        # identical, which would make R a division by zero rather than a number.
        _labels = [
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

        _scale = alt.Scale(type="log", nice=False, domain=plot_domain)
        _n_columns = min(4, len(comparator_info))

        _charts = []
        for _i, _info in enumerate(comparator_info):
            _other = _info["other_plate"]
            # Every panel has the same y axis, this plate, so it is drawn only on the
            # leftmost panel of each row. Each panel has its own x axis, the plate it
            # compares with, so that one is titled on every panel.
            _leftmost = (_i % _n_columns) == 0
            _x_title = f"{_titer_label} on {_other}"
            _y_title = f"{_titer_label} on {plate}" if _leftmost else None
            _y_axis = alt.Axis(labels=_leftmost, ticks=_leftmost)
            # Shared by the points and by the labels annotating them, so that both show
            # the same titers: the QC status has to be calculated before the dropdown that
            # selects on it can filter, and the labels have to be filtered exactly as the
            # points are. `altair` copies a chart when a method is chained onto it, so
            # this can be built on without the layers affecting each other.
            _base = (
                alt.Chart(_info["titers"])
                .transform_calculate(qc_status=_qc_calculate)
                .transform_filter(_qc_selection)
                .transform_filter(_serum_selection)
            )
            _points = (
                _base.transform_calculate(fold_difference=_fold_diff_calculate)
                .encode(
                    alt.X("titer_x:Q", title=_x_title, scale=_scale),
                    alt.Y("titer_y:Q", title=_y_title, scale=_scale, axis=_y_axis),
                    alt.Color(
                        "qc_status:N",
                        scale=alt.Scale(
                            domain=_qc_statuses,
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
                            if _i == 0
                            else None
                        ),
                    ),
                    alt.Shape(
                        "censored:N",
                        title=["is either median titer", "censored at a bound?"],
                        legend=alt.Legend(titleLimit=400) if _i == 0 else None,
                    ),
                    tooltip=[
                        alt.Tooltip("serum:N", title="serum"),
                        alt.Tooltip("virus:N", title="strain"),
                        alt.Tooltip(
                            "titer_y:Q",
                            title=f"{_titer_label} on {plate}",
                            format=".3g",
                        ),
                        alt.Tooltip(
                            "titer_x:Q",
                            title=f"{_titer_label} on {_other}",
                            format=".3g",
                        ),
                        alt.Tooltip(
                            "fold_difference:Q",
                            title="fold difference between these plates",
                            format=".2f",
                        ),
                        alt.Tooltip("n_barcodes_y:Q", title=f"barcodes on {plate}"),
                        alt.Tooltip("n_barcodes_x:Q", title=f"barcodes on {_other}"),
                        alt.Tooltip("qc_status:N", title="QC status"),
                    ],
                    # the moused-over strain gets a larger point with a thick red stroke
                    size=alt.condition(_point_hover, alt.value(180), alt.value(60)),
                    stroke=alt.condition(
                        _point_hover, alt.value("red"), alt.value("black")
                    ),
                    strokeWidth=alt.condition(
                        _point_hover, alt.value(3), alt.value(0.5)
                    ),
                    fillOpacity=alt.condition(
                        _point_hover, alt.value(1), alt.value(0.5)
                    ),
                )
                .mark_point(filled=True, strokeOpacity=1)
            )
            # y = x, drawn from the ends of the axes rather than from the titers so that
            # the line spans the whole panel. This layer repeats the axis encodings of the
            # points layer rather than leaving them unset, because the layers share their
            # axes and a title of `None` here blanks the shared title.
            _diagonal = (
                alt.Chart(diagonal)
                .encode(
                    alt.X("titer:Q", title=_x_title, scale=_scale),
                    alt.Y("titer:Q", title=_y_title, scale=_scale, axis=_y_axis),
                )
                .mark_line(color="gray", strokeDash=[4, 4], strokeWidth=1)
            )
            # The Pearson R of the log titers and the counts of what is being correlated,
            # aggregated from the titers left by the filters in `_base` and so recomputed
            # by `vega` whenever a dropdown changes what is shown. They are drawn in the
            # panel rather than in its title because a title is a fixed string in the chart
            # specification: it can refer to a param, but not to a value aggregated from
            # the data. Each line is its own layer, since a text mark draws one string and
            # the aggregate of a layer cannot be shared with another layer.
            _stats_base = (
                _base.transform_calculate(
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
            _stats = [
                _stats_base.transform_calculate(label=_label)
                .encode(
                    x=alt.value(4),
                    y=alt.value(4 + 12 * _line),
                    text=alt.Text("label:N"),
                )
                .mark_text(align="left", baseline="top", color="black", fontSize=10)
                for (_line, _label) in enumerate(_labels)
            ]
            _charts.append(
                alt.layer(_diagonal, _points, *_stats).properties(
                    width=200,
                    height=200,
                    title=alt.Title(
                        f"{_other} ({_info['other_date']})",
                        fontSize=13,
                    ),
                )
            )
        mo.output.append(
            alt.concat(*_charts, columns=_n_columns)
            # all of the params are declared here rather than on the panels so that they
            # are in scope for all of them
            .add_params(_qc_selection, _serum_selection, _point_hover).properties(
                title=alt.Title(
                    f"{plate} (dated {plate_date}) versus the other plates of {group}",
                    subtitle=[
                        (
                            "each point is a viral strain for one serum; titers are "
                            "medians over each plate's barcodes"
                        ),
                        (
                            f"the y axis of every panel is the {_titer_label} on "
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
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Write the correlations to a CSV
    Write the Pearson correlations for each other plate and serum to a CSV file.
    Note that each pair of plates appears in both plates' files, with the two plates
    swapped:
    """)
    return


@app.cell
def _(comparator_info, corrs_csv, mo, pd, plate, plate_date):
    # the columns are specified so that the file still has a header if there are no
    # correlations
    _corrs = pd.DataFrame(
        [
            {
                "plate": plate,
                "date": plate_date,
                "other_plate": _info["other_plate"],
                "other_date": _info["other_date"],
                **_serum_corrs,
            }
            for _info in comparator_info
            for _serum_corrs in _info["per_serum"]
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
    mo.output.append(mo.md(f"Writing the correlations to `{corrs_csv}`"))
    _corrs.to_csv(corrs_csv, index=False, float_format="%.4g")
    mo.output.append(_corrs)
    return


if __name__ == "__main__":
    app.run()

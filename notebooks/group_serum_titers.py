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
    # Titers for a serum in a group
    Analyze titers for a serum assigned to a group, aggregating replicates which may be across multiple plates.
    """)
    return


@app.cell
def _():
    import io
    import itertools
    import pickle
    import sys

    import altair as alt
    import marimo as mo
    import matplotlib
    import matplotlib.pyplot as plt
    import neutcurve
    import numpy
    import pandas as pd
    from neutcurve.marimo_utils import display_fig_marimo
    from ruamel import yaml

    _ = alt.data_transformers.disable_max_rows()

    # faster plotting of neut curves
    matplotlib.style.use("fast")
    return (
        alt,
        display_fig_marimo,
        io,
        itertools,
        mo,
        neutcurve,
        numpy,
        pd,
        pickle,
        plt,
        sys,
        yaml,
    )


@app.cell
def _():
    def narrow_for_altair(df, drop=()):
        """Narrow a data frame that is about to be embedded in an `altair` chart.

        `altair` embeds the whole frame as inline JSON, so anything superfluous is paid
        for on every row. Drops `drop`, which should be the columns that are the same on
        every row: put them back with `transform_calculate` if the chart needs them, so
        that the value appears once in the chart rather than once per row. Also rounds
        the floats to four significant figures, since the values from the curve fits and
        from the medians of them otherwise run to 18 characters apiece, which is far
        finer than a pixel and finer than the formats of the tooltips. The frames written
        to the output CSVs are not passed through here, so this affects only the charts.
        """
        df = df.drop(columns=list(drop))
        return df.assign(
            **{
                col: df[col].map(lambda x: float(f"{x:.4g}"))
                for col in df.columns
                if df[col].dtype == float
            }
        )

    return (narrow_for_altair,)


@app.cell
def _(pickle, sys):
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
        # set `context_pickle_path` to valid option if using edit more
        context_pickle_path = None
        # context_pickle_path = pathlib.Path("test_example/results/sera/serum_M099d0/serum_M099d0_titers_context.pickle")

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
    pickle_fits = context["input"]["pickles"]
    fits_csvs = context["input"]["fits_csvs"]
    per_rep_titers_csv = context["output"]["per_rep_titers"]
    titers_csv = context["output"]["titers"]
    curves_pdf = context["output"]["curves_pdf"]
    output_pickle = context["output"]["pickle"]
    qc_drops_file = context["output"]["qc_drops"]
    viral_strain_plot_order = context["params"]["viral_strain_plot_order"]
    serum_titer_as = context["params"]["serum_titer_as"]
    qc_thresholds = context["params"]["qc_thresholds"]
    curve_display_method = context["params"]["curve_display_method"]
    dilution_factor_or_concentration = context["params"][
        "dilution_factor_or_concentration"
    ]
    concentration_units = context["params"]["concentration_units"]
    serum = context["wildcards"]["serum"]
    group = context["wildcards"]["group"]

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

    mo.output.append(mo.md(f"Processing `{group=}`, `{serum=}`"))
    return (
        concentration_units,
        curve_display_method,
        curves_pdf,
        dilution_factor_or_concentration,
        fits_csvs,
        group,
        output_pickle,
        per_rep_titers_csv,
        pickle_fits,
        qc_drops_file,
        qc_thresholds,
        serum,
        serum_titer_as,
        titers_csv,
        viral_strain_plot_order,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Get all titers for this plate
    Combine all the pickled `neutcurve.CurveFits` from plates for this serum into a single `neutcurve.CurveFits`:
    """)
    return


@app.cell
def _(group, mo, neutcurve, pickle, pickle_fits, serum):
    mo.output.append(
        mo.md(
            f"Combining the curve fits for `group={group!r}`, `serum={serum!r}` from `pickle_fits={pickle_fits!r}`"
        )
    )
    fits_to_combine = []
    for fname in pickle_fits:
        with open(fname, "rb") as f_combine:
            fits_to_combine.append(pickle.load(f_combine))
    fits_noqc = neutcurve.CurveFits.combineCurveFits(fits_to_combine, sera=[serum])
    return (fits_noqc,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Get the plate that each replicate was measured on.
    Each plate records this in its `curvefits.csv`, so it does not have to be parsed out
    of the replicate name (a `plate_barcode`), which would be ambiguous when one plate
    name is a prefix of another:
    """)
    return


@app.cell
def _(fits_csvs, mo, pd, serum):
    _plate_reps = pd.concat(
        [pd.read_csv(f)[["plate", "replicate"]] for f in fits_csvs], ignore_index=True
    ).drop_duplicates()
    assert (
        len(_plate_reps) == _plate_reps["replicate"].nunique()
    ), f"a replicate is assigned to more than one plate:\n{_plate_reps}"
    plate_of_replicate = _plate_reps.set_index("replicate")["plate"].to_dict()

    # plates in the order they are listed for this serum, used to order the plate pairs
    serum_plates = list(dict.fromkeys(_plate_reps["plate"]))
    mo.output.append(
        mo.md(f"`{serum}` was measured on {len(serum_plates)} plate(s): {serum_plates}")
    )
    return plate_of_replicate, serum_plates


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Indicate how we are calculating the titer:
    """)
    return


@app.cell
def _(dilution_factor_or_concentration, mo, serum_titer_as):
    mo.output.append(mo.md(f"Calculating with `serum_titer_as={serum_titer_as!r}`"))
    if dilution_factor_or_concentration == "dilution_factor":
        assert serum_titer_as in {"nt50", "midpoint"}, serum_titer_as
    else:
        assert serum_titer_as in {"ic50", "midpoint"}, serum_titer_as
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Get all the per-replicate fit params with the titers.
    We also convert the IC50 to NT50, and take inverse of midpoint to get it on same scale as NT50s:
    """)
    return


@app.cell
def _(
    concentration_units,
    dilution_factor_or_concentration,
    fits_noqc,
    group,
    mo,
    plate_of_replicate,
    serum,
    serum_titer_as,
    viral_strain_plot_order,
):
    _fit_params = fits_noqc.fitParams(average_only=False, no_average=True).assign(
        group=group,
        plate=lambda x: x["replicate"].map(plate_of_replicate),
    )
    assert _fit_params["plate"].notnull().all(), (
        "could not assign a plate to every replicate:\n"
        f"{_fit_params[_fit_params['plate'].isnull()]}"
    )
    if dilution_factor_or_concentration == "dilution_factor":
        # titer is a reciprocal serum dilution: convert IC50 to NT50 and take the
        # inverse of the midpoint to put it on the same scale. Taking the reciprocal
        # flips the direction of a bound, so map lower <-> upper.
        per_rep_titers = _fit_params.assign(
            nt50=lambda x: 1 / x["ic50"],
            midpoint=lambda x: 1 / x["midpoint_bound"],
            titer=lambda x: (
                x["midpoint"] if serum_titer_as == "midpoint" else x["nt50"]
            ),
            titer_bound=lambda x: (
                x["midpoint_bound_type"]
                if serum_titer_as == "midpoint"
                else x["ic50_bound"]
            ).map({"lower": "upper", "upper": "lower", "interpolated": "interpolated"}),
            titer_as=serum_titer_as,
            titer_units="reciprocal_dilution",
        )[
            [
                "group",
                "serum",
                "virus",
                "plate",
                "replicate",
                "titer",
                "titer_bound",
                "titer_as",
                "titer_units",
                "nt50",
                "midpoint",
                "top",
                "bottom",
                "slope",
            ]
        ]
    else:
        # titer is the fitted concentration directly (e.g. IC50 in `concentration_units`);
        # no reciprocal is taken, so bounds are reported as-is (no lower <-> upper swap).
        per_rep_titers = _fit_params.assign(
            midpoint=lambda x: x["midpoint_bound"],
            titer=lambda x: (
                x["midpoint"] if serum_titer_as == "midpoint" else x["ic50"]
            ),
            titer_bound=lambda x: (
                x["midpoint_bound_type"]
                if serum_titer_as == "midpoint"
                else x["ic50_bound"]
            ),
            titer_as=serum_titer_as,
            titer_units=concentration_units,
        )[
            [
                "group",
                "serum",
                "virus",
                "plate",
                "replicate",
                "titer",
                "titer_bound",
                "titer_as",
                "titer_units",
                "ic50",
                "midpoint",
                "top",
                "bottom",
                "slope",
            ]
        ]
    assert per_rep_titers.notnull().all().all()

    if len(
        invalid_titer_as := per_rep_titers.query(
            "(titer_as in ['nt50', 'ic50']) and top <= 0.5"
        )
    ):
        raise ValueError(
            "There are titers computed as nt50/ic50 when curve top <= 0.5:\n"
            f"{invalid_titer_as}"
        )
    assert len(per_rep_titers) == per_rep_titers["replicate"].nunique()

    # get viruses in the order to plot them
    viruses = sorted(per_rep_titers["virus"].unique())
    if viral_strain_plot_order is not None:
        if not set(viruses).issubset(viral_strain_plot_order):
            raise ValueError(
                "`viral_strain_plot_order` lacks some viruses with titers:\n"
                + str(set(viruses) - set(viral_strain_plot_order))
            )
        viruses = [v for v in viral_strain_plot_order if v in viruses]
    mo.output.append(
        mo.md(f"`{serum}` has titers for a total of {len(viruses)} viruses")
    )
    return per_rep_titers, viruses


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Correlate NT50s with midpoints of curves
    Plot the correlation of the NT50s with the midpoint (this is an interactive plot, mouse over points for details).
    This plot can help you determine if you made the correct choice of `serum_titer_as` when choosing to use the midpoint or NT50 for the titer.
    For titers where they are well correlated it should not matter which you chose.
    But if there are titers far from the correlation line, you should look at those measurements and curves to make sure you made the correct choice of calculating the titer as the NT50 versus midpoint:
    """)
    return


@app.cell
def _(
    alt,
    dilution_factor_or_concentration,
    group,
    mo,
    narrow_for_altair,
    per_rep_titers,
    serum,
):
    _virus_selection_1 = alt.selection_point(
        fields=["virus"], on="mouseover", empty=False
    )
    # the per-replicate potency column is 'nt50' (dilution) or 'ic50' (concentration)
    _potency_col = (
        "nt50" if dilution_factor_or_concentration == "dilution_factor" else "ic50"
    )
    # `group`, `serum`, and `titer_as` are dropped because they are not in the tooltips
    # below, and `titer_units` because it is the same on every row and so is put back as a
    # constant by the `transform_calculate`
    _titer_units = per_rep_titers["titer_units"].unique()
    assert len(_titer_units) == 1, _titer_units
    _chart_df = narrow_for_altair(
        per_rep_titers, drop=["group", "serum", "titer_as", "titer_units"]
    )
    midpoint_vs_nt50_chart = (
        alt.Chart(_chart_df)
        .add_params(_virus_selection_1)
        .transform_calculate(titer_units=f"'{_titer_units[0]}'")
        .encode(
            alt.X(_potency_col, scale=alt.Scale(type="log", nice=False, padding=8)),
            alt.Y("midpoint", scale=alt.Scale(type="log", nice=False, padding=8)),
            alt.Color("titer_bound"),
            strokeWidth=alt.condition(_virus_selection_1, alt.value(3), alt.value(0)),
            size=alt.condition(_virus_selection_1, alt.value(100), alt.value(60)),
            tooltip=[
                alt.Tooltip(c, format=".2g") if _chart_df[c].dtype == float else c
                for c in _chart_df.columns
            ]
            # a calculated field has no `pandas` dtype, so its type is given explicitly
            + [alt.Tooltip("titer_units:N")],
        )
        .mark_circle(stroke="black", fillOpacity=0.45, color="black")
        .properties(
            width=350,
            height=350,
            title=f"{_potency_col.upper()} versus midpoint for {group} {serum}",
        )
        .configure_axis(grid=False)
    )
    mo.output.append(midpoint_vs_nt50_chart)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Plot median titers and determine if they pass QC
    Get the median titers for each virus across replicates, then add these median titers to the per-replicate titers and calculate the fold-change in titer between each replicate and its median.
    Finally, for each virus indicate whether it passes the QC:
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
def _(get_median_bound, mo, numpy, pd, per_rep_titers, qc_thresholds, viruses):
    mo.output.append(mo.md(f"Using the following `qc_thresholds={qc_thresholds!r}`"))

    median_titers_noqc = (
        per_rep_titers.sort_values("titer")  # for getting median bound
        .groupby(["group", "serum", "virus", "titer_as", "titer_units"], as_index=False)
        .aggregate(
            titer=pd.NamedAgg("titer", "median"),
            n_replicates=pd.NamedAgg("replicate", "count"),
            titer_sem=pd.NamedAgg("titer", "sem"),
            titer_bound=pd.NamedAgg("titer_bound", get_median_bound),
        )
    )

    per_rep_titers_w_fc = (
        per_rep_titers.merge(
            median_titers_noqc[["group", "serum", "virus", "titer"]].rename(
                columns={"titer": "median_titer"}
            ),
            validate="many_to_one",
            on=["group", "serum", "virus"],
        )
        .assign(
            fc_from_median=lambda x: numpy.where(
                x["titer"] > x["median_titer"],
                x["titer"] / x["median_titer"],
                x["median_titer"] / x["titer"],
            ),
        )
        .drop(columns=["group", "serum", "titer_as", "median_titer"])
    )

    median_titers_noqc = median_titers_noqc.merge(
        per_rep_titers_w_fc.groupby("virus", as_index=False).aggregate(
            max_fc_from_median=pd.NamedAgg("fc_from_median", "max")
        ),
        on="virus",
        validate="one_to_one",
    ).assign(
        fails_min_reps=lambda x: x["n_replicates"] < qc_thresholds["min_replicates"],
        fails_max_fc=lambda x: (
            x["max_fc_from_median"] >= qc_thresholds["max_fold_change_from_median"]
        ),
        fails_qc=lambda x: x["fails_min_reps"] | x["fails_max_fc"],
        fails_qc_reason=lambda x: (
            x.apply(
                lambda r: ", ".join(
                    (["min_replicates"] if r["fails_min_reps"] else [])
                    + (["max_fold_change_from_median"] if r["fails_max_fc"] else [])
                ),
                axis=1,
            )
        ),
    )

    # get viruses failing QC in order to plot
    viruses_failing_qc = (
        median_titers_noqc.query("fails_qc")
        .set_index("virus")["fails_qc_reason"]
        .to_dict()
    )
    viruses_failing_qc = {
        v: viruses_failing_qc[v] for v in viruses if v in viruses_failing_qc
    }

    median_titers_noqc = median_titers_noqc.drop(
        columns=["fails_min_reps", "fails_max_fc", "fails_qc_reason"]
    )

    per_rep_titers_w_fc = per_rep_titers_w_fc.merge(
        median_titers_noqc[["virus", "fails_qc"]],
        on="virus",
        validate="many_to_one",
    )
    return median_titers_noqc, per_rep_titers_w_fc, viruses_failing_qc


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Now plot the per-replicate and median titers, indicating any viruses that failed QC.
    Note that potentially some of these titers may still be retained if the viruses in question are specified in `viruses_ignore_qc` of `qc_thresholds`.
    """)
    return


@app.cell
def _(
    alt,
    concentration_units,
    dilution_factor_or_concentration,
    group,
    median_titers_noqc,
    mo,
    narrow_for_altair,
    per_rep_titers_w_fc,
    qc_thresholds,
    serum,
    viruses,
):
    _virus_selection_2 = alt.selection_point(
        fields=["virus"], on="mouseover", empty=False
    )
    _titer_axis_title = (
        "titer"
        if dilution_factor_or_concentration == "dilution_factor"
        else f"titer ({concentration_units})"
    )
    # the columns dropped here are the same on every row, and are put back as constants by
    # the `transform_calculate` of each layer so that the tooltips still show them
    _constants = {
        "group": group,
        "serum": serum,
        "titer_as": median_titers_noqc["titer_as"].unique()[0],
        "titer_units": median_titers_noqc["titer_units"].unique()[0],
    }
    _calculate = {col: f"'{val}'" for col, val in _constants.items()}
    _per_rep_df = narrow_for_altair(per_rep_titers_w_fc, drop=["titer_units"])
    _median_df = narrow_for_altair(median_titers_noqc, drop=list(_constants))
    per_rep_chart = (
        alt.Chart(_per_rep_df)
        .transform_calculate(titer_units=_calculate["titer_units"])
        .encode(
            alt.X(
                "titer",
                title=_titer_axis_title,
                scale=alt.Scale(nice=False, padding=5, type="log"),
            ),
            alt.Y("virus", sort=viruses),
            alt.Fill(
                "fails_qc",
                title=f"fails qc_thresholds['min_replicates']={qc_thresholds['min_replicates']!r}, qc_thresholds['max_fold_change_from_median']={qc_thresholds['max_fold_change_from_median']!r}",
                legend=alt.Legend(titleLimit=500),
            ),
            alt.Shape("titer_bound"),
            strokeWidth=alt.condition(_virus_selection_2, alt.value(2), alt.value(0)),
            tooltip=[
                (alt.Tooltip(c, format=".3g") if _per_rep_df[c].dtype == float else c)
                for c in _per_rep_df
            ]
            # a calculated field has no `pandas` dtype, so its type is given explicitly
            + [alt.Tooltip("titer_units:N")],
        )
        .mark_point(
            size=35, filled=True, fillOpacity=0.5, strokeOpacity=1, stroke="black"
        )
    )
    median_chart = (
        alt.Chart(_median_df)
        .transform_calculate(**_calculate)
        .encode(
            alt.X(
                "titer",
                title=_titer_axis_title,
                scale=alt.Scale(nice=False, padding=5, type="log"),
            ),
            alt.Y("virus", sort=viruses),
            alt.Fill("fails_qc"),
            alt.Shape("titer_bound"),
            strokeWidth=alt.condition(_virus_selection_2, alt.value(2), alt.value(0.5)),
            tooltip=[
                (alt.Tooltip(c, format=".3g") if _median_df[c].dtype == float else c)
                for c in _median_df
            ]
            + [alt.Tooltip(f"{c}:N") for c in _calculate],
        )
        .mark_point(
            size=75, filled=True, fillOpacity=0.9, strokeOpacity=1, stroke="black"
        )
    )
    titer_chart = (
        (per_rep_chart + median_chart)
        .add_params(_virus_selection_2)
        .properties(
            height=alt.Step(11),
            width=250,
            title=f"{group} {serum} median (large points) and per-replicate (small points) titers",
        )
        .configure_axis(grid=False)
    )
    mo.output.append(titer_chart)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Plot individual curves for any viruses failing QC
    Plot individual curves for viruses failing QC.
    Note that potentially some of these titers may still be retained if the viruses in question are specified in `viruses_ignore_qc` of `qc_thresholds`.
    """)
    return


@app.cell
def _(
    curve_display_method,
    display_fig_marimo,
    fits_noqc,
    group,
    mo,
    plt,
    serum,
    viruses_failing_qc,
):
    mo.output.append(
        mo.md(
            f"Neutralization curves for the {len(viruses_failing_qc)} viruses failing QC:"
        )
    )
    if len(viruses_failing_qc):
        _fig_failing, _ = fits_noqc.plotReplicates(
            viruses=viruses_failing_qc,
            attempt_shared_legend=False,
            legendfontsize=8,
            titlesize=12,
            ncol=4,
            heightscale=1.2,
            widthscale=1.2,
            subplot_titles="{virus}",
            draw_in_bounds=True,
        )
        _ = _fig_failing.suptitle(
            f"neutralization curves for viruses failing QC for {group} {serum}",
            y=1,
            fontsize=18,
            fontweight="bold",
        )
        _fig_failing.tight_layout()
        mo.output.append(
            display_fig_marimo(_fig_failing, display_method=curve_display_method)
        )
        plt.close(_fig_failing)
    else:
        mo.output.append(mo.md("No curves fail QC"))
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Get the viruses to drop for QC failures
    Drop any viruses that fail QC and are not specified in `viruses_ignore_qc` of `qc_thresholds`.
    """)
    return


@app.cell
def _(io, mo, qc_drops_file, qc_thresholds, viruses_failing_qc, yaml):
    viruses_to_drop = {
        v: reason
        for v, reason in viruses_failing_qc.items()
        if v not in qc_thresholds["viruses_ignore_qc"]
    }
    mo.output.append(mo.md(f"Dropping {len(viruses_to_drop)} viruses for failing QC:"))
    yaml_buffer_viruses_drop = io.StringIO()
    yaml.YAML(typ="rt").dump(viruses_to_drop, stream=yaml_buffer_viruses_drop)
    mo.output.append(mo.md(f"```yaml\n{yaml_buffer_viruses_drop.getvalue()}```"))
    if nkept := (len(viruses_failing_qc) - len(viruses_to_drop)):
        kept_viruses = {
            v: reason
            for v, reason in viruses_failing_qc.items()
            if v in qc_thresholds["viruses_ignore_qc"]
        }
        mo.output.append(
            mo.md(
                f"Retaining {nkept} viruses that fail QC because they are in `viruses_ignore_qc`: "
                f"{kept_viruses}"
            )
        )
    mo.output.append(mo.md(f"Writing QC drops to `{qc_drops_file}`"))
    with open(qc_drops_file, "w") as f_qc_drops:
        yaml.YAML(typ="rt").dump(viruses_to_drop, f_qc_drops)
    return (viruses_to_drop,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Write the individual per-replicate titers to a file.
    These are all of the replicates that passed the QC applied when their plate was
    processed, and so are not filtered by the per-serum QC just applied above.
    Instead, the `dropped_by_qc` column indicates whether a replicate's virus was dropped
    by that per-serum QC, making this file a superset of the titers written below.
    Note that the per-serum QC drops a virus rather than an individual replicate, so this
    column has the same value for all replicates of a virus:
    """)
    return


@app.cell
def _(mo, per_rep_titers, per_rep_titers_csv, viruses_to_drop):
    mo.output.append(
        mo.md(f"Writing per-replicate titers to `{per_rep_titers_csv}`"),
    )
    per_rep_titers.assign(
        dropped_by_qc=lambda x: x["virus"].isin(set(viruses_to_drop))
    ).to_csv(per_rep_titers_csv, index=False, float_format="%.4g")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Plate-to-plate correlation of titers
    If this serum was measured on more than one plate, compare the titers between each
    pair of plates to assess plate-to-plate reproducibility.
    The titer for a strain on a plate is the median over all of that plate's barcodes,
    all of which have already passed the per-plate QC applied when the plate was
    processed.

    All strains are shown, colored by whether they are retained in the final titers or
    were dropped just above by the per-serum QC on replicate-to-replicate variation.
    **Use the `show strains` dropdown at the bottom of the plots to show only the retained
    or only the dropped strains**, and set it back to `all` to show all of them.
    Note that a strain listed in `viruses_ignore_qc` of `qc_thresholds` is colored as
    retained even if it fails the QC, because it is kept in the final titers; mouse over
    a point for the details.
    The Pearson correlation is reported both over all strains and over just the retained
    strains, and the dashed line is `y = x`.
    """)
    return


@app.cell
def _(
    get_median_bound,
    itertools,
    mo,
    narrow_for_altair,
    numpy,
    pd,
    per_rep_titers,
    serum,
    serum_plates,
    viruses,
    viruses_failing_qc,
    viruses_to_drop,
):
    # the titer for a strain on a plate is the median over that plate's barcodes
    plate_titers = (
        per_rep_titers.sort_values("titer")  # for getting median bound
        .groupby(["virus", "plate"], as_index=False)
        .aggregate(
            titer=pd.NamedAgg("titer", "median"),
            n_barcodes=pd.NamedAgg("replicate", "count"),
            titer_bound=pd.NamedAgg("titer_bound", get_median_bound),
        )
    )

    def _pearson_r(df):
        """Pearson R of the log-transformed titers, or `None` if too few strains."""
        if len(df) < 3:
            return None
        return numpy.corrcoef(numpy.log10(df["titer_x"]), numpy.log10(df["titer_y"]))[
            0, 1
        ]

    _dropped = set(viruses_to_drop)
    _pair_viruses = set()
    plate_pair_info = []
    for _plate_x, _plate_y in itertools.combinations(serum_plates, 2):
        _pair_df = (
            plate_titers[plate_titers["plate"] == _plate_x]
            .drop(columns="plate")
            .merge(
                plate_titers[plate_titers["plate"] == _plate_y].drop(columns="plate"),
                on="virus",
                suffixes=("_x", "_y"),
                validate="one_to_one",
            )
            .assign(
                censored=lambda x: (x["titer_bound_x"] != "interpolated")
                | (x["titer_bound_y"] != "interpolated"),
            )
        )
        if not len(_pair_df):
            continue
        _kept_df = _pair_df[~_pair_df["virus"].isin(_dropped)]
        # both axes get the same domain so that the y = x line is a true diagonal. The
        # domain is the range of the titers on either plate, padded by a fraction of that
        # range on a log scale so that points do not sit on the plot edge. The padding is
        # proportional rather than a fixed factor so that it does not swamp a pair whose
        # titers span a narrow range, and the floor keeps a pair with almost no spread
        # (such as a plate compared with an exact duplicate of itself) from collapsing.
        _lo = min(_pair_df["titer_x"].min(), _pair_df["titer_y"].min())
        _hi = max(_pair_df["titer_x"].max(), _pair_df["titer_y"].max())
        _pad = max(0.05 * (numpy.log10(_hi) - numpy.log10(_lo)), numpy.log10(1.1))
        plate_pair_info.append(
            {
                "pair": f"{_plate_x} vs {_plate_y}",
                "plate_x": _plate_x,
                "plate_y": _plate_y,
                "n": len(_pair_df),
                "r": _pearson_r(_pair_df),
                "n_kept": len(_kept_df),
                "r_kept": _pearson_r(_kept_df),
                # rounded for the same reason as the titers in `narrow_for_altair`, and by
                # so much less than the padding that it cannot bring a point to the edge
                "domain": [
                    float(f"{10 ** (numpy.log10(_lo) - _pad):.4g}"),
                    float(f"{10 ** (numpy.log10(_hi) + _pad):.4g}"),
                ],
                # keep the plotted data frame narrow: each subplot gets just its own pair
                # of plates, so no column names the pair. That costs no extra rows, since
                # `altair` consolidates every frame into the chart's `datasets` and the
                # rows are partitioned among the subplots either way. The plates are named
                # in the subplot's axis titles and tooltips rather than in a column, and
                # the QC status of a strain is determined by `virus` so it goes in the
                # `virus_qc_lookup` below.
                "titers": narrow_for_altair(
                    _pair_df[
                        [
                            "virus",
                            "titer_x",
                            "titer_y",
                            "n_barcodes_x",
                            "n_barcodes_y",
                            "censored",
                        ]
                    ].reset_index(drop=True)
                ),
            }
        )
        _pair_viruses.update(_pair_df["virus"])

    if plate_pair_info:
        # end points of each subplot's y = x line, which are the ends of its axes rather
        # than of its data, so that the line spans the whole subplot
        for _info in plate_pair_info:
            _info["diagonal"] = pd.DataFrame({"titer": _info["domain"]})

        def _qc_reason(v):
            """Explain the QC status of a strain."""
            if v in viruses_to_drop:
                return f"fails {viruses_to_drop[v]}"
            elif v in viruses_failing_qc:
                return f"fails {viruses_failing_qc[v]} but retained via `viruses_ignore_qc`"
            else:
                return "passes QC"

        virus_qc_lookup = pd.DataFrame(
            {"virus": [v for v in viruses if v in _pair_viruses]}
        ).assign(
            qc_status=lambda x: numpy.where(
                x["virus"].isin(_dropped),
                "dropped by per-serum QC",
                "retained in final titers",
            ),
            qc_reason=lambda x: x["virus"].map(_qc_reason),
        )
        assert len(virus_qc_lookup) == virus_qc_lookup["virus"].nunique()

        mo.output.append(
            mo.md(
                f"Comparing titers between {len(plate_pair_info)} pair(s) of plates for "
                f"`{serum}`"
            )
        )
    else:
        virus_qc_lookup = None
        if len(serum_plates) < 2:
            mo.output.append(
                mo.md(
                    f"`{serum}` was measured on only one plate "
                    f"(`{serum_plates[0]}`), so there is no plate-to-plate correlation "
                    "to show."
                )
            )
        else:
            mo.output.append(
                mo.md(
                    "No pair of plates has any strain measured on both plates, so there "
                    "is no plate-to-plate correlation to show."
                )
            )
    return plate_pair_info, virus_qc_lookup


@app.cell
def _(
    alt,
    concentration_units,
    dilution_factor_or_concentration,
    group,
    mo,
    plate_pair_info,
    serum,
    virus_qc_lookup,
):
    if plate_pair_info:
        _titer_label = (
            "titer"
            if dilution_factor_or_concentration == "dilution_factor"
            else f"titer ({concentration_units})"
        )

        # A dropdown is used to select which strains to show rather than binding the
        # selection to the QC-status legend. Binding to the legend does not work here: the
        # legend of a concatenated chart belongs to whichever subplot draws it, so the
        # binding never reaches the other subplots. An input binding does work, because the
        # param is declared on the concatenated chart itself (see `add_params` below) and
        # so is in scope for every subplot. Selecting "all" shows every strain.
        _qc_statuses = ["retained in final titers", "dropped by per-serum QC"]
        _qc_selection = alt.selection_point(
            fields=["qc_status"],
            bind=alt.binding_select(
                options=[None, *_qc_statuses],
                labels=["all", *_qc_statuses],
                name="show strains",
            ),
            name="qc_status_selection",
        )

        # mousing over a strain highlights it on every subplot, not just the one being
        # moused over, because this param is also declared on the concatenated chart
        _virus_hover = alt.selection_point(
            fields=["virus"], on="mouseover", empty=False, name="virus_hover"
        )

        def _scale(info):
            """Log scale shared by both axes of a subplot."""
            return alt.Scale(type="log", nice=False, domain=info["domain"])

        def _fmt_r(r):
            """Format a Pearson R that is `None` when there were too few strains."""
            return "not computed (< 3 strains)" if r is None else f"{r:.2f}"

        _charts = []
        for _i, _info in enumerate(plate_pair_info):
            _base = alt.Chart(_info["titers"]).transform_lookup(
                lookup="virus",
                from_=alt.LookupData(
                    virus_qc_lookup,
                    key="virus",
                    fields=["qc_status", "qc_reason"],
                ),
            )
            _x_title = f"{_titer_label} on {_info['plate_x']}"
            _y_title = f"{_titer_label} on {_info['plate_y']}"
            _points = (
                _base.transform_filter(_qc_selection)
                .encode(
                    alt.X("titer_x:Q", title=_x_title, scale=_scale(_info)),
                    alt.Y("titer_y:Q", title=_y_title, scale=_scale(_info)),
                    alt.Color(
                        "qc_status:N",
                        scale=alt.Scale(
                            domain=_qc_statuses,
                            range=["MediumBlue", "OrangeRed"],
                        ),
                        legend=(
                            alt.Legend(
                                title=[
                                    "strain QC status",
                                    "(use the 'show strains'",
                                    "dropdown below the plots",
                                    "to show only one of these)",
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
                        alt.Tooltip("virus:N", title="strain"),
                        alt.Tooltip(
                            "titer_x:Q",
                            title=f"{_titer_label} on {_info['plate_x']}",
                            format=".3g",
                        ),
                        alt.Tooltip(
                            "titer_y:Q",
                            title=f"{_titer_label} on {_info['plate_y']}",
                            format=".3g",
                        ),
                        alt.Tooltip(
                            "n_barcodes_x:Q", title=f"barcodes on {_info['plate_x']}"
                        ),
                        alt.Tooltip(
                            "n_barcodes_y:Q", title=f"barcodes on {_info['plate_y']}"
                        ),
                        alt.Tooltip("qc_status:N", title="QC status"),
                        alt.Tooltip("qc_reason:N", title="QC detail"),
                    ],
                    # the moused-over strain gets a larger point with a thick red stroke
                    size=alt.condition(_virus_hover, alt.value(180), alt.value(60)),
                    stroke=alt.condition(
                        _virus_hover, alt.value("red"), alt.value("black")
                    ),
                    strokeWidth=alt.condition(
                        _virus_hover, alt.value(3), alt.value(0.5)
                    ),
                    fillOpacity=alt.condition(
                        _virus_hover, alt.value(1), alt.value(0.5)
                    ),
                )
                .mark_point(filled=True, strokeOpacity=1)
            )
            # y = x, drawn from the ends of the axes rather than from the titers so that
            # the line spans the whole subplot. This layer repeats the axis titles of the
            # points layer rather than leaving them unset, because the layers share their
            # axes and a title of `None` here blanks the shared title.
            _diagonal = (
                alt.Chart(_info["diagonal"])
                .encode(
                    alt.X("titer:Q", title=_x_title, scale=_scale(_info)),
                    alt.Y("titer:Q", title=_y_title, scale=_scale(_info)),
                )
                .mark_line(color="gray", strokeDash=[4, 4], strokeWidth=1)
            )
            _charts.append(
                (_diagonal + _points).properties(
                    width=200,
                    height=200,
                    title=alt.Title(
                        _info["pair"],
                        subtitle=[
                            f"R = {_fmt_r(_info['r'])} (all {_info['n']} strains)",
                            (
                                f"R = {_fmt_r(_info['r_kept'])} "
                                f"({_info['n_kept']} retained strains)"
                            ),
                        ],
                        fontSize=14,
                        subtitleFontSize=10,
                    ),
                )
            )
        plate_correlation_chart = (
            alt.concat(*_charts, columns=min(4, len(_charts)))
            # both params are declared here rather than on the subplots so that they are
            # in scope for all of them
            .add_params(_qc_selection, _virus_hover).properties(
                title=alt.Title(
                    f"plate-to-plate correlation for {group} {serum}",
                    subtitle=(
                        "each point is a viral strain; titers are medians over each "
                        "plate's barcodes"
                    ),
                    fontSize=15,
                    anchor="middle",
                )
            )
            # `titleFontWeight` is set because axis titles are bold by default, and here
            # only the overall and subplot titles should be bold. `labelOverlap` drops
            # tick labels rather than letting them overlap, and is set here rather than on
            # each axis so that the chart spec carries it once instead of once per axis.
            .configure_axis(
                grid=False,
                labelOverlap=True,
                titleFontSize=13,
                titleFontWeight="normal",
            )
            # `limit=0` means no limit, so that neither the titles nor the subtitles of the
            # plot or its subplots are truncated. It is set in the chart's config rather
            # than on each title so that the spec carries it once.
            .configure_title(limit=0)
        )
        mo.output.append(plate_correlation_chart)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Get and plot the neutralization curves for all retained viruses
    First, get the `CurveFits` for just those retained viruses (dropping ones that fail QC), and plot:
    """)
    return


@app.cell
def _(
    curve_display_method,
    curves_pdf,
    display_fig_marimo,
    fits_noqc,
    group,
    mo,
    neutcurve,
    plt,
    serum,
    viruses,
    viruses_to_drop,
):
    fits_qc = neutcurve.CurveFits.combineCurveFits(
        [fits_noqc],
        viruses=[v for v in viruses if v not in viruses_to_drop],
    )
    assert len(viruses) == len(fits_qc.viruses[serum]) + len(viruses_to_drop)
    _fig_retained, _ = fits_qc.plotReplicates(
        attempt_shared_legend=False,
        legendfontsize=8,
        titlesize=12,
        ncol=4,
        heightscale=1.2,
        widthscale=1.2,
        subplot_titles="{virus}",
        viruses=[v for v in viruses if v not in viruses_to_drop],
        draw_in_bounds=True,
    )
    _ = _fig_retained.suptitle(
        f"neutralization curves for retained viruses for {group} {serum}",
        y=1,
        fontsize=18,
        fontweight="bold",
    )
    _fig_retained.tight_layout()
    mo.output.append(
        display_fig_marimo(_fig_retained, display_method=curve_display_method)
    )
    mo.output.append(mo.md(f"Saving to plot of curves to `{curves_pdf}`"))
    _fig_retained.savefig(curves_pdf)
    plt.close(_fig_retained)
    return (fits_qc,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Save the `CurveFits` to a pickle file:
    """)
    return


@app.cell
def _(fits_qc, mo, output_pickle, pickle):
    with open(output_pickle, "wb") as f_out_pickle:
        pickle.dump(fits_qc, f_out_pickle)

    mo.output.append(mo.md(f"Writing curve fits to {output_pickle}"))
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Write the titers (excluding QC dropped viruses) to a CSV:
    """)
    return


@app.cell
def _(median_titers_noqc, mo, titers_csv):
    mo.output.append(mo.md(f"Writing titers to `{titers_csv}`"))

    (
        median_titers_noqc.query("virus not in @viruses_to_drop")[
            [
                "group",
                "serum",
                "virus",
                "titer",
                "titer_bound",
                "titer_sem",
                "n_replicates",
                "titer_as",
                "titer_units",
            ]
        ].to_csv(titers_csv, index=False, float_format="%.4g")
    )
    return


if __name__ == "__main__":
    app.run()

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
    serum,
    serum_titer_as,
    viral_strain_plot_order,
):
    _fit_params = fits_noqc.fitParams(average_only=False, no_average=True).assign(
        group=group,
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
def _(alt, dilution_factor_or_concentration, group, mo, per_rep_titers, serum):
    _virus_selection_1 = alt.selection_point(
        fields=["virus"], on="mouseover", empty=False
    )
    # the per-replicate potency column is 'nt50' (dilution) or 'ic50' (concentration)
    _potency_col = (
        "nt50" if dilution_factor_or_concentration == "dilution_factor" else "ic50"
    )
    midpoint_vs_nt50_chart = (
        alt.Chart(per_rep_titers)
        .add_params(_virus_selection_1)
        .encode(
            alt.X(_potency_col, scale=alt.Scale(type="log", nice=False, padding=8)),
            alt.Y("midpoint", scale=alt.Scale(type="log", nice=False, padding=8)),
            alt.Color("titer_bound"),
            strokeWidth=alt.condition(_virus_selection_1, alt.value(3), alt.value(0)),
            size=alt.condition(_virus_selection_1, alt.value(100), alt.value(60)),
            tooltip=[
                alt.Tooltip(c, format=".2g") if per_rep_titers[c].dtype == float else c
                for c in per_rep_titers.columns
                if c not in {"group", "serum", "titer_as"}
            ],
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
    Write the individual per-replicate titers to a file, this is before any QC has been applied:
    """)
    return


@app.cell
def _(mo, per_rep_titers, per_rep_titers_csv):
    mo.output.append(
        mo.md(
            f"Writing per-replicate titers (without QC filtering) to `{per_rep_titers_csv}`"
        )
    )
    per_rep_titers.to_csv(per_rep_titers_csv, index=False, float_format="%.4g")
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
def _(mo, numpy, pd, per_rep_titers, qc_thresholds, viruses):
    mo.output.append(mo.md(f"Using the following `qc_thresholds={qc_thresholds!r}`"))

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
    per_rep_chart = (
        alt.Chart(per_rep_titers_w_fc)
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
                (
                    alt.Tooltip(c, format=".3g")
                    if per_rep_titers_w_fc[c].dtype == float
                    else c
                )
                for c in per_rep_titers_w_fc
            ],
        )
        .mark_point(
            size=35, filled=True, fillOpacity=0.5, strokeOpacity=1, stroke="black"
        )
    )
    median_chart = (
        alt.Chart(median_titers_noqc)
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
                (
                    alt.Tooltip(c, format=".3g")
                    if median_titers_noqc[c].dtype == float
                    else c
                )
                for c in median_titers_noqc
            ],
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

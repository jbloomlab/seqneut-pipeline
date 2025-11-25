# /// script
# [tool.marimo.runtime]
# auto_instantiate = false
# ///

import marimo

__generated_with = "0.17.6"
app = marimo.App(width="full")


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

    import marimo as mo

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
        # if running in edit mode, set `context_pickle_path` to valid pickle
        context_pickle_path = None
        # context_pickle_path = pathlib.Path("test_example/results/aggregated_titers/aggregate_titers_context.pickle")

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
    return context, mo, pickle


@app.cell
def _(context, mo):
    # Extract variables from context - raises KeyError if required keys missing
    input_pickles = context["input"]["pickles"]
    input_titers = context["input"]["titers"]
    output_pickles = context["output"]["pickles"]
    output_titers = context["output"]["titers"]
    titers_chart_html = context["output"]["titers_chart"]
    groups_sera = context["params"]["groups_sera"]
    groups = context["params"]["groups"]
    viral_strain_plot_order = context["params"]["viral_strain_plot_order"]

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
    return (
        groups,
        groups_sera,
        input_pickles,
        input_titers,
        output_pickles,
        output_titers,
        titers_chart_html,
        viral_strain_plot_order,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    # Aggregate titers across all sera
    Aggregate the titers across all sera.
    """
    )
    return


@app.cell
def _():
    # Setup and read data

    import altair as alt

    import neutcurve

    import pandas as pd

    _ = alt.data_transformers.disable_max_rows()
    return alt, neutcurve, pd


@app.cell
def _(groups, groups_sera):
    # Get the groups ordered by number of sera in each:

    ordered_groups = [
        g
        for (_, g) in sorted(
            [(sum(g == group for (g, _) in groups_sera), group) for group in groups],
            reverse=True,
        )
    ]

    ordered_groups_sera = [
        f"{group} {serum}"
        for group in ordered_groups
        for (g, serum) in sorted(groups_sera)
        if g == group
    ]
    return ordered_groups, ordered_groups_sera


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## Write merged titers for each group
    Get the merged titers and merged `CurveFits` object:
    """
    )
    return


@app.cell
def _(
    groups,
    groups_sera,
    input_pickles,
    input_titers,
    mo,
    neutcurve,
    output_pickles,
    output_titers,
    pd,
    pickle,
):
    assert len(input_titers) == len(input_pickles) == len(groups_sera)
    assert len(groups) == len(output_titers) == len(output_pickles)

    titers = pd.concat(
        [pd.read_csv(f_titer) for f_titer in input_titers], ignore_index=True
    )
    assert len(titers) == len(titers.groupby(["group", "serum", "virus"]))
    for group, f_out_titer in zip(groups, output_titers):
        mo.output.append(
            mo.md(f"Writing aggregated titers for `{group}` to `{f_out_titer}`")
        )
        titers.query("group == @group").to_csv(
            f_out_titer, index=False, float_format="%.4g"
        )

    for group, f_out_pickle in zip(groups, output_pickles):
        fits_list = []
        for (g, serum), pickle_f in zip(groups_sera, input_pickles):
            if g == group:
                with open(pickle_f, "rb") as f_in_pickle:
                    fits_list.append(pickle.load(f_in_pickle))
        curvefits = neutcurve.CurveFits.combineCurveFits(fits_list)
        mo.output.append(
            mo.md(f"Pickling aggregated `CurveFits` for `{group}` to `{f_out_pickle}`")
        )
        with open(f_out_pickle, "wb") as f_out_pickle_file:
            pickle.dump(curvefits, f_out_pickle_file)
    return curvefits, titers


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## Plot all the titers
    """
    )
    return


@app.cell
def _(
    alt,
    curvefits,
    mo,
    ordered_groups,
    ordered_groups_sera,
    titers,
    titers_chart_html,
    viral_strain_plot_order,
):
    if viral_strain_plot_order is None:
        viral_strain_plot_order_used = sorted(curvefits.allviruses)
    elif set(curvefits.allviruses) - set(viral_strain_plot_order):
        raise ValueError(
            "not all viruses in `viral_strain_plot_order`:\n"
            + f"{set(curvefits.allviruses) - set(viral_strain_plot_order)=}"
        )
    else:
        viral_strain_plot_order_used = viral_strain_plot_order

    viruses = [v for v in viral_strain_plot_order_used if v in curvefits.allviruses]
    assert set(viruses) == set(curvefits.allviruses)

    serum_selection = alt.selection_point(
        fields=["group_serum"],
        bind="legend",
        toggle="true",
        empty=False,
        on="mouseover",
        clear="mouseout",
    )

    group_selection = alt.selection_point(
        fields=["group"],
        value=None,
        bind=alt.binding_select(
            options=[None] + ordered_groups,
            labels=["all"] + ordered_groups,
            name="group",
        ),
    )

    base_chart = (
        alt.Chart(titers)
        .transform_calculate(group_serum=alt.datum["group"] + " " + alt.datum["serum"])
        .add_params(serum_selection, group_selection)
        .transform_filter(group_selection)
        .encode(
            alt.X(
                "titer",
                title="neutralization titer",
                scale=alt.Scale(nice=False, padding=5, type="log"),
                axis=alt.Axis(labelOverlap=True),
            ),
        )
        .properties(height=alt.Step(10))
    )

    median_titers_chart = (
        base_chart.encode(
            alt.Y("virus", title=None, sort=viruses, axis=alt.Axis(labelLimit=500)),
            color=alt.value("black"),
            tooltip=["virus", alt.Tooltip("titer", format=".3g")],
        )
        .transform_aggregate(titer="median(titer)", groupby=["virus"])
        .mark_line(point={"filled": True, "size": 50})
        .properties(width=125, title="median titer")
    )

    per_titers_chart = (
        base_chart.encode(
            alt.Y(
                "virus",
                title=None,
                sort=viruses,
                axis=alt.Axis(labels=False, ticks=False),
            ),
            alt.Detail("group_serum:N"),
            alt.Color(
                "group_serum:N",
                scale=alt.Scale(range=["black", "black"]),
                legend=alt.Legend(
                    orient="right",
                    columns=1 + int(1.6 * len(ordered_groups_sera)) // len(viruses),
                    symbolLimit=0,
                    symbolFillColor="black",
                    title="serum (click to select)",
                ),
                sort=ordered_groups_sera,
            ),
            opacity=alt.condition(serum_selection, alt.value(1), alt.value(0.15)),
            strokeWidth=alt.condition(serum_selection, alt.value(3), alt.value(1)),
            tooltip=[
                alt.Tooltip(c, format=".3g") if titers[c].dtype == float else c
                for c in titers.columns
            ],
        )
        .mark_line(point={"filled": True, "size": 20})
        .properties(width=225, title="per-serum titers")
    )

    titers_chart = (
        alt.hconcat(median_titers_chart, per_titers_chart, spacing=5)
        .configure_axis(grid=False)
        .configure_legend(padding=10, labelOffset=2, columnPadding=8, labelLimit=400)
        .properties(
            title=alt.TitleParams(
                "Interactive chart of serum neutralization titers",
                subtitle=[
                    "Median titers are at left; per-serum titers are at right.",
                    "Mouseover points for details.",
                    "Click points on serum legend at right to select sera to bold on per-serum chart.",
                    "Use dropdown at bottom to select serum groups to show.",
                ],
                fontSize=15,
                dx=100,
                dy=-5,
            ),
            autosize=alt.AutoSizeParams(resize=True),
        )
    )

    mo.output.append(mo.md(f"Saving chart to `{titers_chart_html}`"))
    titers_chart.save(titers_chart_html)
    mo.output.append(titers_chart)
    return


if __name__ == "__main__":
    app.run()

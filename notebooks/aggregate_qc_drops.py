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
        # set `context_pickle_path` to valid pickle if running via marimo edit
        context_pickle_path = None
        # context_pickle_path = pathlib.Path("test_example/results/qc_drops/aggregate_qc_drops_context.pickle")

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
    return context, mo


@app.cell
def _(context, mo):
    # Extract variables from context - raises KeyError if required keys missing
    input_plate_qc_drops = context["input"]["plate_qc_drops"]
    input_groups_sera_qc_drops = context["input"]["groups_sera_qc_drops"]
    output_plate_qc_drops = context["output"]["plate_qc_drops"]
    output_barcode_qc_drops = context["output"]["barcode_qc_drops"]
    output_groups_sera_qc_drops = context["output"]["groups_sera_qc_drops"]
    plates = context["params"]["plates"]
    groups_sera = context["params"]["groups_sera"]

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
        groups_sera,
        input_groups_sera_qc_drops,
        input_plate_qc_drops,
        output_barcode_qc_drops,
        output_groups_sera_qc_drops,
        output_plate_qc_drops,
        plates,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    # Aggregate and analyze the drops from QC-ing the plates and sera
    """
    )
    return


@app.cell
def _():
    import altair as alt

    import pandas as pd

    from ruamel.yaml import YAML

    yaml = YAML(typ="rt")

    _ = alt.data_transformers.disable_max_rows()
    return alt, pd, yaml


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## Analyze plate QC drops
    Read QC drops for individual plates into a merged dictionary, write it to YAML, and also convert to a DataFrame.
    If you really want to look into the details of what is being dropped, you will want to look at that merged YAML file.
    """
    )
    return


@app.cell
def _(input_plate_qc_drops, mo, output_plate_qc_drops, pd, plates, yaml):
    # read dictionary of QC drops
    assert len(plates) == len(input_plate_qc_drops)
    plate_qc_drops = {}
    for _plate, _qc_drops_yaml in zip(plates, input_plate_qc_drops):
        with open(_qc_drops_yaml) as _f:
            plate_qc_drops[_plate] = yaml.load(_f)
    assert len(plate_qc_drops) == len(input_plate_qc_drops)

    mo.output.append(mo.md(f"Writing merged plate drops to {output_plate_qc_drops}"))
    with open(output_plate_qc_drops, "w") as _f:
        yaml.dump(plate_qc_drops, stream=_f)

    # convert dictionary of QC drops into list of tuples
    plate_qc_drop_tups = [
        (plate_key, droptype_key, drop_key, reason)
        for plate_key, plate_val in plate_qc_drops.items()
        for droptype_key, droptype_val in plate_val.items()
        for drop_key, reason in droptype_val.items()
    ]

    # create data frame of QC drops
    plate_qc_drops_df = pd.DataFrame(
        plate_qc_drop_tups, columns=["plate", "drop type", "drop", "reason"]
    )
    return plate_qc_drops, plate_qc_drops_df


@app.cell
def _(pd, plate_qc_drops_df):
    plate_qc_drop_counts = plate_qc_drops_df.groupby(
        ["plate", "drop type", "reason"], as_index=False
    ).aggregate(n_drops=pd.NamedAgg("drop", "nunique"))
    assert plate_qc_drop_counts["n_drops"].sum() == len(plate_qc_drops_df)
    return (plate_qc_drop_counts,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    Now plot the number of drops for each plate.
    You should be worried (maybe re-do or discard) any plates with a very large number of drops:
    """
    )
    return


@app.cell
def _(alt, plate_qc_drop_counts, plates):
    plate_selection = alt.selection_point(fields=["plate"], on="mouseover", empty=False)

    plate_qc_drop_counts_chart = (
        alt.Chart(plate_qc_drop_counts)
        .add_params(plate_selection)
        .encode(
            alt.X(
                "n_drops",
                title="number of drops",
            ),
            alt.Y(
                "plate",
                sort=plates,
                title=None,
                axis=alt.Axis(labelFontStyle="bold", labelFontSize=11),
            ),
            alt.Column(
                "drop type",
                title=None,
                spacing=5,
                header=alt.Header(
                    labelFontSize=12, labelFontStyle="bold", labelPadding=1
                ),
            ),
            alt.Color(
                "reason",
                legend=alt.Legend(
                    orient="top", columns=1, labelLimit=230, title=None, padding=1
                ),
            ),
            strokeWidth=alt.condition(plate_selection, alt.value(3), alt.value(0.5)),
            tooltip=plate_qc_drop_counts.columns.tolist(),
        )
        .mark_bar(height={"band": 0.8}, stroke="black")
        .properties(
            width=230,
            height=alt.Step(16),
            title=alt.TitleParams(
                "Number of QC drops when processing plates", anchor="middle", dy=-2
            ),
        )
        .configure_axis(grid=False)
        .resolve_scale(color="independent", x="independent")
    )

    plate_qc_drop_counts_chart
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## Analyze barcode QC drops
    If a barcode is dropped especially often across plates, that could indicate something problematic with that barcode such that it should be removed altogether from the library analysis.

    First, write a YAML with this information:
    """
    )
    return


@app.cell
def _(mo, output_barcode_qc_drops, plate_qc_drops, yaml):
    barcode_qc_drops = {}

    for _plate, plate_d in plate_qc_drops.items():
        for key, val in [
            tup for tup in plate_d.items() if tup[0].startswith("barcode") and tup[1]
        ]:
            for desc, reason in val.items():
                desc_entries = desc.split(None, maxsplit=1)
                barcode = desc_entries[0]
                if len(desc_entries) == 2:
                    description = f"{_plate} {desc_entries[1]}"
                elif len(desc_entries) == 1:
                    description = _plate
                else:
                    raise RuntimeError("should not get here")
                if barcode not in barcode_qc_drops:
                    barcode_qc_drops[barcode] = {}
                if key not in barcode_qc_drops[barcode]:
                    barcode_qc_drops[barcode][key] = {}
                barcode_qc_drops[barcode][key][description] = reason

    # sort keys
    def sort_nested(d):
        if isinstance(d, dict):
            return dict(sorted((key, sort_nested(val)) for key, val in d.items()))
        else:
            return d

    barcode_qc_drops = sort_nested(barcode_qc_drops)

    mo.output.append(
        mo.md((f"Writing merged barcode drops to {output_barcode_qc_drops}"))
    )
    with open(output_barcode_qc_drops, "w") as _f:
        yaml.dump(barcode_qc_drops, stream=_f)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    Now make a plot showing how often each barcode is dropped for each reason:
    """
    )
    return


@app.cell
def _(pd, plate_qc_drops_df):
    barcode_drops = (
        plate_qc_drops_df.query("`drop type`.str.startswith('barcode')")
        .assign(barcode=lambda x: x["drop"].str.split().str[0])
        .groupby(["drop type", "barcode"], as_index=False)
        .aggregate(
            plates_where_dropped=pd.NamedAgg("plate", "nunique"),
            total_drops=pd.NamedAgg("plate", "count"),
        )
    )
    return (barcode_drops,)


@app.cell
def _(alt, barcode_drops):
    barcode_selection = alt.selection_point(
        fields=["barcode"], on="mouseover", empty=False
    )

    barcode_drops_chart = (
        alt.Chart(barcode_drops)
        .add_params(barcode_selection)
        .encode(
            alt.X(
                "total_drops",
                title="times barcode dropped",
            ),
            alt.Y(
                "barcode",
                sort=alt.SortField("total_drops", order="descending"),
                axis=alt.Axis(labelFontSize=9),
            ),
            alt.Column(
                "drop type",
                title=None,
                spacing=8,
                header=alt.Header(
                    labelFontSize=12, labelFontStyle="bold", labelPadding=1
                ),
            ),
            strokeWidth=alt.condition(barcode_selection, alt.value(3), alt.value(0.5)),
            tooltip=barcode_drops.columns.tolist(),
        )
        .mark_bar(height={"band": 0.8}, stroke="black")
        .properties(
            width=200,
            height=alt.Step(10),
            title=alt.TitleParams(
                "Number of QC drops when processing plates", anchor="middle", dy=-2
            ),
        )
        .configure_axis(grid=False)
        .resolve_scale(color="independent", x="independent", y="independent")
    )

    barcode_drops_chart
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## Analyze the groups/sera QC
    Analyze the QC performed on the groups/sera, which involves completely dropping titers for certain virus-sera pairs.

    Read the QC for different groups/sera into a merged dictionary, write it to YAML, and also convert to a DataFrame.
    If you really want to look into the details of what is being dropped, you will want to look at that merged YAML file.
    """
    )
    return


@app.cell
def _(
    groups_sera,
    input_groups_sera_qc_drops,
    mo,
    output_groups_sera_qc_drops,
    pd,
    yaml,
):
    # read dictionary of QC drops
    assert len(groups_sera) == len(input_groups_sera_qc_drops)
    groups_sera_qc_drops = {}
    for (group, serum), _qc_drops_yaml in zip(groups_sera, input_groups_sera_qc_drops):
        if group not in groups_sera_qc_drops:
            groups_sera_qc_drops[group] = {}
        with open(_qc_drops_yaml) as _f:
            groups_sera_qc_drops[group][serum] = yaml.load(_f)

    mo.output.append(
        mo.md(f"Writing merged groups/sera drops to {output_groups_sera_qc_drops}")
    )
    with open(output_groups_sera_qc_drops, "w") as _f:
        yaml.dump(groups_sera_qc_drops, stream=_f)

    # convert dictionary of QC drops into list of tuples
    groups_sera_qc_drop_tups = [
        (group_key, serum_key, virus, reason)
        for group_key, group_val in groups_sera_qc_drops.items()
        for serum_key, serum_val in group_val.items()
        for virus, reason in serum_val.items()
    ]

    # create data frame of QC drops
    groups_sera_qc_drops_df = pd.DataFrame(
        groups_sera_qc_drop_tups, columns=["group", "serum", "virus", "reason"]
    )
    return (groups_sera_qc_drops_df,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    Plot the number of viruses dropped for each group/serum.
    If a group/serum has many missed viruses, then you will lack a lot of titers and so it may be worth reviewing the cause of the drops.
    """
    )
    return


@app.cell
def _(alt, groups_sera_qc_drops_df, pd):
    groups_sera_n_drops = groups_sera_qc_drops_df.groupby(
        ["group", "serum", "reason"], as_index=False
    ).aggregate(n_viruses=pd.NamedAgg("virus", "nunique"))
    assert groups_sera_n_drops["n_viruses"].sum() == len(groups_sera_qc_drops_df)

    groups_sera_n_drops_chart = (
        alt.Chart(groups_sera_n_drops)
        .encode(
            alt.X("n_viruses", title="number of viruses dropped"),
            alt.Y("serum"),
            alt.Row("group"),
            alt.Color(
                "reason", title="reason dropped", legend=alt.Legend(labelLimit=350)
            ),
            tooltip=groups_sera_n_drops.columns.tolist(),
        )
        .mark_bar(height={"band": 0.8})
        .properties(
            width=250,
            height=alt.Step(13),
            title="Number of viruses dropped at serum QC for each serum",
        )
        .configure_axis(grid=False)
        .resolve_scale(y="independent", x="independent")
    )

    groups_sera_n_drops_chart
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    Plot the number of sera for which each virus is dropped during serum QC.
    If a virus is dropped for many sera, that may indicate some issue with that virus in assays:
    """
    )
    return


@app.cell
def _(alt, groups_sera_qc_drops_df, pd):
    virus_n_drops = groups_sera_qc_drops_df.groupby(
        ["group", "virus", "reason"], as_index=False
    ).aggregate(n_sera=pd.NamedAgg("serum", "nunique"))
    assert virus_n_drops["n_sera"].sum() == len(groups_sera_qc_drops_df)

    virus_n_drops_chart = (
        alt.Chart(virus_n_drops)
        .encode(
            alt.X("n_sera", title="number of sera for which virus is dropped"),
            alt.Y("virus", sort=alt.SortField("n_sera", order="descending")),
            alt.Row("group"),
            alt.Color(
                "reason", title="reason dropped", legend=alt.Legend(labelLimit=350)
            ),
            tooltip=virus_n_drops.columns.tolist(),
        )
        .mark_bar(height={"band": 0.8})
        .properties(
            width=250,
            height=alt.Step(13),
            title="Number of sera for which each virus is dropped at serum QC",
        )
        .configure_axis(grid=False)
        .resolve_scale(y="independent", x="independent")
    )

    virus_n_drops_chart
    return


if __name__ == "__main__":
    app.run()

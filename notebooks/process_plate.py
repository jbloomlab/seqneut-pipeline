# /// script
# [tool.marimo.runtime]
# auto_instantiate = false
# ///

import marimo

__generated_with = "0.17.6"
app = marimo.App(width="full")


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Process plate counts to get fraction infectivities and fit curves
    This notebook analyzes a plate of sequencing-based neutralization assays.

    The plots are interactive, so you can mouseover points for details, use the mouse-scroll to zoom and pan, and use interactive dropdowns at the bottom of the plots.
    """)
    return


@app.cell
def _():
    # Import Python modules

    import io
    import pickle
    import sys
    import warnings

    import altair as alt

    import matplotlib.style

    import neutcurve
    from neutcurve.colorschemes import CBPALETTE, CBMARKERS

    import numpy

    import pandas as pd

    import ruamel.yaml as yaml

    import marimo as mo

    _ = alt.data_transformers.disable_max_rows()

    # avoid clutter w RuntimeWarning during curve fitting
    warnings.filterwarnings("ignore", category=RuntimeWarning)

    # faster plotting of neut curves
    matplotlib.style.use("fast")

    # get marimo PDF from matplotlib figure
    def mo_pdf_from_fig(fig):
        buf = io.BytesIO()
        with matplotlib.rc_context({
            "pdf.use14corefonts": True,
            "pdf.compression": 7,
            "path.simplify": True,
            "path.simplify_threshold": 0.2,
        }):
            fig.savefig(buf, format="pdf", metadata={})
        buf.seek(0)
        return mo.pdf(src=buf, width="100%")
    return (
        CBMARKERS,
        CBPALETTE,
        alt,
        io,
        mo,
        mo_pdf_from_fig,
        neutcurve,
        numpy,
        pd,
        pickle,
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
        # If running in edit mode, make sure `context_pickle_path` below specifies a valid pickle
        context_pickle_path = None
        # context_pickle_path = pathlib.Path("test_example/results/plates/plate11/process_plate11_context.pickle")

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
def _(context, io, mo, pd, yaml):
    # Extract variables from context - raises KeyError if required keys missing
    count_csvs = context["input"]["count_csvs"]
    fate_csvs = context["input"]["fate_csvs"]
    qc_drops_yaml = context["output"]["qc_drops"]
    frac_infectivity_csv = context["output"]["frac_infectivity_csv"]
    fits_csv = context["output"]["fits_csv"]
    fits_pickle = context["output"]["fits_pickle"]
    viral_barcodes = context["params"]["viral_barcodes"]
    neut_standard_barcodes = context["params"]["neut_standard_barcodes"]
    samples = context["params"]["samples"]
    plate = context["wildcards"]["plate"]
    plate_params = context["params"]["plate_params"]

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

    # Process plate_params only if we have real data
    if plate_params:
        # get thresholds turning lists into tuples as needed
        manual_drops = {
            filter_type: [tuple(w) if isinstance(w, list) else w for w in filter_drops]
            for (filter_type, filter_drops) in plate_params["manual_drops"].items()
        }
        group = plate_params["group"]
        qc_thresholds = plate_params["qc_thresholds"]
        curvefit_params = plate_params["curvefit_params"]
        curvefit_qc = plate_params["curvefit_qc"]
        if "barcode_serum_replicates_ignore_curvefit_qc" in curvefit_qc:
            curvefit_qc["barcode_serum_replicates_ignore_curvefit_qc"] = [
                tuple(w)
                for w in curvefit_qc["barcode_serum_replicates_ignore_curvefit_qc"]
            ]

        mo.output.append(mo.md(f"Processing `{plate}`"))

        samples_df = pd.DataFrame(plate_params["samples"])
        mo.output.append(mo.md(f"Plate has {len(samples)} samples (wells)"))
        assert all(
            (len(samples_df) == samples_df[c].nunique())
            for c in ["well", "sample", "sample_noplate"]
        )
        assert len(samples_df) == len(
            samples_df.groupby(["serum_replicate", "dilution_factor"])
        )
        assert len(samples) == len(count_csvs) == len(fate_csvs) == len(samples_df)

        for d, key, title in [
            (manual_drops, "manual_drops", "Data manually specified to drop:"),
            (qc_thresholds, "qc_thresholds", "QC thresholds applied to data:"),
            (curvefit_params, "curvefit_params", "Curve-fitting parameters:"),
            (curvefit_qc, "curvefit_qc", "Curve-fitting QC:"),
        ]:
            mo.output.append(mo.md(f"{title}"))
            yaml_buffer_params = io.StringIO()
            yaml.YAML(typ="rt").dump({key: d}, stream=yaml_buffer_params)
            mo.output.append(mo.md(f"```yaml\n{yaml_buffer_params.getvalue()}```"))
    else:
        # Stub context - initialize with empty values
        manual_drops = {}
        group = "unknown"
        qc_thresholds = {}
        curvefit_params = {}
        curvefit_qc = {}
        samples_df = pd.DataFrame()
    return (
        count_csvs,
        curvefit_params,
        curvefit_qc,
        fate_csvs,
        fits_csv,
        fits_pickle,
        frac_infectivity_csv,
        group,
        manual_drops,
        neut_standard_barcodes,
        plate,
        qc_drops_yaml,
        qc_thresholds,
        samples,
        samples_df,
        viral_barcodes,
    )


@app.cell
def _(manual_drops):
    # Set up dictionary to keep track of wells, barcodes, well-barcodes, and serum-replicates that are dropped:
    qc_drops = {
        "wells": {},
        "barcodes": {},
        "barcode_wells": {},
        "barcode_serum_replicates": {},
        "serum_replicates": {},
    }

    assert set(manual_drops).issubset(
        qc_drops
    ), f"{manual_drops.keys()=}, {qc_drops.keys()}"
    return (qc_drops,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Statistics on barcode-parsing for each sample
    Make interactive chart of the "fates" of the sequencing reads parsed for each sample on the plate.

    If most sequencing reads are not "valid barcodes", this could potentially indicate some problem in the sequencing or barcode set you are parsing.

    Potential fates are:
     - *valid barcode*: barcode that matches a known virus or neutralization standard, we hope most reads are this.
     - *invalid barcode*: a barcode with proper flanking sequences, but does not match a known virus or neutralization standard. If you  have a lot of reads of this type, it is probably a good idea to look at the invalid barcode CSVs (in the `./results/barcode_invalid/` subdirectory created by the pipeline) to see what these invalid barcodes are.
     - *unparseable barcode*: could not parse a barcode from this read as there was not a sequence of the correct length with the appropriate flanking sequence.
     - *invalid outer flank*: if using an outer upstream or downstream region (`upstream2` or `downstream2` for the [illuminabarcodeparser](https://jbloomlab.github.io/dms_variants/dms_variants.illuminabarcodeparser.html)), reads that are otherwise valid except for this outer flank. Typically you would be using `upstream2` if you have a plate index embedded in your primer, and reads with this classification correspond to a different index than the one for this plate.
     - *low quality barcode*: low-quality or `N` nucleotides in barcode, could indicate problem with sequencing.
     - *failed chastity filter*: reads that failed the Illumina chastity filter, if these are reported in the FASTQ (they may not be).

    Also, if the number of reads per sample is very uneven, that could indicate that you did not do a good job of balancing the different samples in the Illumina sequencing.
    """)
    return


@app.cell
def _(alt, fate_csvs, pd, plate, samples, samples_df):
    fates = (
        pd.concat([pd.read_csv(f).assign(sample=s) for f, s in zip(fate_csvs, samples)])
        .merge(samples_df, validate="many_to_one", on="sample")
        .assign(
            fate_counts=lambda x: x.groupby("fate")["count"].transform("sum"),
            sample_well=lambda x: x["sample_noplate"] + " (" + x["well"] + ")",
        )
        .query("fate_counts > 0")[  # only keep fates with at least one count
            [
                "fate",
                "count",
                "well",
                "serum_replicate",
                "sample_well",
                "dilution_factor",
            ]
        ]
    )

    assert len(fates) == len(fates.drop_duplicates())

    serum_replicates = sorted(fates["serum_replicate"].unique())
    sample_wells = list(
        fates.sort_values(["serum_replicate", "dilution_factor"])["sample_well"]
    )

    serum_selection = alt.selection_point(
        fields=["serum_replicate"],
        bind=alt.binding_select(
            options=[None] + serum_replicates,
            labels=["all"] + serum_replicates,
            name="serum",
        ),
    )

    fates_chart = (
        alt.Chart(fates)
        .add_params(serum_selection)
        .transform_filter(serum_selection)
        .encode(
            alt.X("count", scale=alt.Scale(nice=False, padding=3)),
            alt.Y(
                "sample_well",
                title=None,
                sort=sample_wells,
            ),
            alt.Color("fate", sort=sorted(fates["fate"].unique(), reverse=True)),
            alt.Order("fate", sort="descending"),
            tooltip=fates.columns.tolist(),
        )
        .mark_bar(height={"band": 0.85})
        .properties(
            height=alt.Step(10),
            width=200,
            title=f"Barcode parsing for {plate}",
        )
        .configure_axis(grid=False)
    )

    fates_chart
    return sample_wells, serum_replicates, serum_selection


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Read barcode counts and apply manually specified drops
    Read the counts per barcode, then apply any manually specified drops.
    """)
    return


@app.cell
def _(
    count_csvs,
    neut_standard_barcodes,
    pd,
    sample_wells,
    samples,
    samples_df,
    serum_replicates,
    viral_barcodes,
):
    # get barcode counts
    counts = (
        pd.concat(
            [pd.read_csv(c).assign(sample=s) for c, s in zip(count_csvs, samples)]
        )
        .merge(samples_df, validate="many_to_one", on="sample")
        .drop(columns=["replicate", "plate", "fastq"])
        .assign(sample_well=lambda x: x["sample_noplate"] + " (" + x["well"] + ")")
    )

    # classify barcodes as viral or neut standard
    barcode_class = pd.concat(
        [
            pd.DataFrame(viral_barcodes).assign(neut_standard=False),
            pd.DataFrame(neut_standard_barcodes).assign(
                neut_standard=True, strain=pd.NA
            ),
        ],
        ignore_index=True,
    )

    # merge counts and classification of barcodes
    assert set(counts["barcode"]) == set(barcode_class["barcode"])
    counts = counts.merge(barcode_class, on="barcode", validate="many_to_one")
    assert set(sample_wells) == set(counts["sample_well"])
    assert set(serum_replicates) == set(counts["serum_replicate"])
    return (counts,)


@app.cell
def _(counts, manual_drops, mo, qc_drops):
    counts_qc_1 = counts.copy()
    for filter_type, filter_drops in manual_drops.items():
        mo.output.append(
            mo.md(
                f"Dropping {len(filter_drops)} {filter_type} specified in manual_drops"
            )
        )
        assert filter_type in qc_drops
        qc_drops[filter_type].update(
            {w: "manual_drop" for w in filter_drops if not isinstance(w, list)}
        )
        if filter_type == "barcode_wells":
            counts_qc_1 = counts_qc_1[
                ~counts.assign(
                    barcode_well=lambda x: x.apply(
                        lambda r: (r["barcode"], r["well"]), axis=1
                    )
                )["barcode_well"].isin(qc_drops[filter_type])
            ]
        elif filter_type == "barcode_serum_replicates":
            counts_qc_1 = counts_qc_1[
                ~counts_qc_1.assign(
                    barcode_serum_replicate=lambda x: x.apply(
                        lambda r: (r["barcode"], r["serum_replicate"]), axis=1
                    )
                )["barcode_serum_replicate"].isin(qc_drops[filter_type])
            ]
        elif filter_type == "wells":
            counts_qc_1 = counts_qc_1[~counts_qc_1["well"].isin(qc_drops[filter_type])]
        elif filter_type == "barcodes":
            counts_qc_1 = counts_qc_1[
                ~counts_qc_1["barcode"].isin(qc_drops[filter_type])
            ]
        elif filter_type == "serum_replicates":
            counts_qc_1 = counts_qc_1[
                ~counts_qc_1["serum_replicate"].isin(qc_drops[filter_type])
            ]
        elif filter_type == "barcode_serum_replicates":
            counts_qc_1 = counts_qc_1[
                ~counts_qc_1["barcode_serum_replicate"].isin(qc_drops[filter_type])
            ]
        else:
            assert filter_type in set(counts_qc_1.columns)
            counts_qc_1 = counts_qc_1[
                ~counts_qc_1[filter_type].isin(qc_drops[filter_type])
            ]

    if not manual_drops:
        mo.output.append(mo.md("No drops specified in manual_drops"))
    return (counts_qc_1,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Average counts per barcode in each well
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Plot average counts per barcode.
    If a sample has inadequate barcode counts, it may not have good enough statistics for accurate analysis, and a QC-threshold is applied:
    """)
    return


@app.cell
def _(
    alt,
    counts_qc_1,
    mo,
    pd,
    plate,
    qc_drops,
    qc_thresholds,
    sample_wells,
    serum_selection,
):
    # Compute average barcode counts per well
    avg_barcode_counts = (
        counts_qc_1.groupby(
            ["well", "serum_replicate", "sample_well"],
            dropna=False,
            as_index=False,
        )
        .aggregate(avg_count=pd.NamedAgg("count", "mean"))
        .assign(
            fails_qc=lambda x: (
                x["avg_count"] < qc_thresholds["avg_barcode_counts_per_well"]
            )
        )
    )

    # Create chart
    avg_barcode_counts_chart = (
        alt.Chart(avg_barcode_counts)
        .add_params(serum_selection)
        .transform_filter(serum_selection)
        .encode(
            alt.X(
                "avg_count",
                title="average barcode counts per well",
                scale=alt.Scale(nice=False, padding=3),
            ),
            alt.Y("sample_well", sort=sample_wells),
            alt.Color(
                "fails_qc",
                title=f"fails qc_thresholds['avg_barcode_counts_per_well']={qc_thresholds['avg_barcode_counts_per_well']!r}",
                legend=alt.Legend(titleLimit=500),
            ),
            tooltip=[
                (
                    alt.Tooltip(c, format=".3g")
                    if avg_barcode_counts[c].dtype == float
                    else c
                )
                for c in avg_barcode_counts.columns
            ],
        )
        .mark_bar(height={"band": 0.85})
        .properties(
            height=alt.Step(10),
            width=250,
            title=f"Average barcode counts per well for {plate}",
        )
        .configure_axis(grid=False)
    )

    # Display the chart
    mo.output.append(avg_barcode_counts_chart)

    # Drop wells failing QC
    avg_barcode_counts_per_well_drops = list(
        avg_barcode_counts.query("fails_qc")["well"]
    )
    mo.output.append(
        mo.md(
            f"Dropping {len(avg_barcode_counts_per_well_drops)} wells for failing "
            f"`qc_thresholds['avg_barcode_counts_per_well']={qc_thresholds['avg_barcode_counts_per_well']!r}`: "
            f"{avg_barcode_counts_per_well_drops}"
        )
    )
    qc_drops["wells"].update(
        {w: "avg_barcode_counts_per_well" for w in avg_barcode_counts_per_well_drops}
    )
    counts_qc_2 = counts_qc_1[~counts_qc_1["well"].isin(qc_drops["wells"])]
    return (counts_qc_2,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Fraction of counts from neutralization standard
    Determine the fraction of counts from the neutralization standard in each sample, and make sure this fraction passess the QC threshold.
    """)
    return


@app.cell
def _(
    alt,
    counts_qc_2,
    pd,
    plate,
    qc_thresholds,
    sample_wells,
    serum_selection,
):
    # Compute neutralization standard fractions per well
    neut_standard_fracs = (
        counts_qc_2.assign(
            neut_standard_count=lambda x: x["count"] * x["neut_standard"].astype(int)
        )
        .groupby(
            ["well", "serum_replicate", "sample_well"],
            dropna=False,
            as_index=False,
        )
        .aggregate(
            total_count=pd.NamedAgg("count", "sum"),
            neut_standard_count=pd.NamedAgg("neut_standard_count", "sum"),
        )
        .assign(
            neut_standard_frac=lambda x: x["neut_standard_count"] / x["total_count"],
            fails_qc=lambda x: (
                x["neut_standard_frac"]
                < qc_thresholds["min_neut_standard_frac_per_well"]
            ),
        )
    )

    # Create chart
    neut_standard_fracs_chart = (
        alt.Chart(neut_standard_fracs)
        .add_params(serum_selection)
        .transform_filter(serum_selection)
        .encode(
            alt.X(
                "neut_standard_frac",
                title="frac counts from neutralization standard per well",
                scale=alt.Scale(nice=False, padding=3),
            ),
            alt.Y("sample_well", sort=sample_wells),
            alt.Color(
                "fails_qc",
                title=f"fails qc_thresholds['min_neut_standard_frac_per_well']={qc_thresholds['min_neut_standard_frac_per_well']!r}",
                legend=alt.Legend(titleLimit=500),
            ),
            tooltip=[
                (
                    alt.Tooltip(c, format=".3g")
                    if neut_standard_fracs[c].dtype == float
                    else c
                )
                for c in neut_standard_fracs.columns
            ],
        )
        .mark_bar(height={"band": 0.85})
        .properties(
            height=alt.Step(10),
            width=250,
            title=f"Neutralization-standard fracs per well for {plate}",
        )
        .configure_axis(grid=False)
        .configure_legend(titleLimit=1000)
    )

    # Display the chart
    neut_standard_fracs_chart
    return (neut_standard_fracs,)


@app.cell
def _(counts_qc_2, mo, neut_standard_fracs, qc_drops, qc_thresholds):
    # drop wells failing QC
    min_neut_standard_frac_per_well_drops = list(
        neut_standard_fracs.query("fails_qc")["well"]
    )
    mo.output.append(
        mo.md(
            f"Dropping {len(min_neut_standard_frac_per_well_drops)} wells for failing `qc_thresholds['min_neut_standard_frac_per_well']={qc_thresholds['min_neut_standard_frac_per_well']!r}`: {min_neut_standard_frac_per_well_drops}"
        )
    )
    qc_drops["wells"].update(
        {
            w: "min_neut_standard_frac_per_well"
            for w in min_neut_standard_frac_per_well_drops
        }
    )
    counts_qc_3 = counts_qc_2[~counts_qc_2["well"].isin(qc_drops["wells"])]
    return (counts_qc_3,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Consistency and minimum fractions for barcodes
    We examine the fraction of counts attributable to each barcode. We do this splitting the data two ways:

     1. Looking at all viral (but not neut-standard) barcodes only for the no-serum samples (wells).

     2. Looking at just the neut-standard barcodes for all samples (wells).

    The reason is that if the experiment is set up perfectly, these fractions should be the same across all samples for each barcode.
    (We do not expect viral barcodes to have consistent fractions across no-serum samples as they will be neutralized differently depending on strain).

    We plot these fractions in interactive plots (you can mouseover points and zoom) so you can identify barcodes that fail the expected consistency QC thresholds.

    We also make sure the barcodes meet specified QC minimum thresholds for all samples, and flag any that do not.
    """)
    return


@app.cell
def _(
    alt,
    counts_qc_3,
    mo,
    numpy,
    pd,
    plate,
    qc_drops,
    qc_thresholds,
    sample_wells,
):
    # Create barcode selection parameter for mouseover highlighting
    barcode_selection = alt.selection_point(
        fields=["barcode"],
        on="mouseover",
        empty=False,
    )

    # Analyze barcode evenness for two groups: neut-standard barcodes and viral barcodes
    counts_qc_4 = counts_qc_3.copy()

    for is_neut_standard, df in counts_qc_3.groupby("neut_standard"):
        # Determine which QC to apply based on barcode type
        if is_neut_standard:
            mo.output.append(
                mo.md(
                    f"{'=' * 89}\n\nAnalyzing neut-standard barcodes from all samples (wells)"
                )
            )
            qc_name = "per_neut_standard_barcode_filters"
        else:
            mo.output.append(
                mo.md(
                    f"{'=' * 89}\n\nAnalyzing all barcodes from no-serum samples (wells)"
                )
            )
            qc_name = "no_serum_per_viral_barcode_filters"
            df = df.query("serum == 'none'")

        # Compute barcode fractions and fold changes
        df = df.assign(
            sample_counts=lambda x: x.groupby("sample")["count"].transform("sum"),
            count_frac=lambda x: x["count"] / x["sample_counts"],
            median_count_frac=lambda x: x.groupby("barcode")["count_frac"].transform(
                "median"
            ),
            fold_change_from_median=lambda x: numpy.where(
                x["count_frac"] > x["median_count_frac"],
                x["count_frac"] / x["median_count_frac"],
                x["median_count_frac"] / x["count_frac"],
            ),
        )[
            ["barcode", "count", "sample_well", "count_frac", "fold_change_from_median"]
            + ([] if is_neut_standard else ["strain"])
        ]

        # Apply QC thresholds
        qc = qc_thresholds[qc_name]
        mo.output.append(mo.md(f"Apply QC `{qc_name}`: `{qc}`"))

        fails_qc = (
            df.assign(
                fails_qc=lambda x: ~(
                    (x["count_frac"] >= qc["min_frac"])
                    & (x["fold_change_from_median"] <= qc["max_fold_change"])
                )
            )
            .groupby("barcode", as_index=False)
            .aggregate(n_wells_fail_qc=pd.NamedAgg("fails_qc", "sum"))
            .assign(fails_qc=lambda x: x["n_wells_fail_qc"] >= qc["max_wells"])[
                ["barcode", "fails_qc"]
            ]
        )

        df = df.merge(fails_qc, on="barcode", validate="many_to_one")

        # Create evenness chart
        evenness_chart = (
            alt.Chart(df)
            .add_params(barcode_selection)
            .encode(
                alt.X(
                    "count_frac",
                    title=(
                        "barcode's fraction of neut standard counts"
                        if is_neut_standard
                        else "barcode's fraction of non-neut standard counts"
                    ),
                    scale=alt.Scale(nice=False, padding=5),
                ),
                alt.Y("sample_well", sort=sample_wells),
                alt.Fill(
                    "fails_qc",
                    title=f"fails {qc_name}",
                    legend=alt.Legend(titleLimit=500),
                ),
                strokeWidth=alt.condition(
                    barcode_selection, alt.value(2), alt.value(0)
                ),
                size=alt.condition(barcode_selection, alt.value(60), alt.value(35)),
                tooltip=[
                    alt.Tooltip(c, format=".2g") if df[c].dtype == float else c
                    for c in df.columns
                ],
            )
            .mark_circle(fillOpacity=0.45, stroke="black", strokeOpacity=1)
            .properties(
                height=alt.Step(10),
                width=300,
                title=alt.TitleParams(
                    (
                        f"{plate} all samples, neut-standard barcodes"
                        if is_neut_standard
                        else f"{plate} no-serum samples, all barcodes"
                    ),
                    subtitle="x-axis is zoomable (use mouse scroll/pan)",
                ),
            )
            .configure_axis(grid=False)
            .configure_legend(titleLimit=1000)
            .interactive()
        )

        # Display the chart
        mo.output.append(evenness_chart)

        # Drop barcodes failing QC
        barcode_drops = list(fails_qc.query("fails_qc")["barcode"])
        mo.output.append(
            mo.md(
                f"Dropping {len(barcode_drops)} barcodes for failing `qc={qc!r}`: {barcode_drops}"
            )
        )
        qc_drops["barcodes"].update(
            {bc: "min_neut_standard_frac_per_well" for bc in barcode_drops}
        )
        counts_qc_4 = counts_qc_4[~counts_qc_4["barcode"].isin(qc_drops["barcodes"])]
    return barcode_selection, counts_qc_4


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Compute fraction infectivity

    The fraction infectivity for viral barcode $v_b$ in sample $s$ is computed as:

    $$F_{v_b,s} = \frac{c_{v_b,s} / \left(\sum_{n_b} c_{n_b,s}\right)}{{\rm median}_{s_0}\left[ c_{v_b,s_0} / \left(\sum_{n_b} c_{n_b,s_0}\right)\right]}$$

    where
     - $c_{v_b,s}$ is the counts of viral barcode $v_b$ in sample $s$.
     - $\sum_{n_b} c_{n_b,s}$ is the sum of the counts for all neutralization standard barcodes $n_b$ for sample $s$.
     - $c_{v_b,s_0}$ is the counts of viral barcode $v_b$ in no-serum sample $s_0$.
     - $\sum_{n_b} c_{n_b,s_0}$ is the sum of the counts for all neutralization standard barcodes $n_b$ for no-serum sample $s_0$.
     - ${\rm median}_{s_0}\left[ c_{v_b,s_0} / \left(\sum_{n_b} c_{n_b,s_0}\right)\right]$ is the median taken across all no-serum samples of the counts of viral barcode $v_b$ versus the total counts for all neutralization standard barcodes.

    First, compute the total neutralization-standard counts for each sample (well).
    Plot these, and drop any wells that do not meet the QC threshold.
    """)
    return


@app.cell
def _(
    alt,
    counts_qc_4,
    pd,
    plate,
    qc_thresholds,
    sample_wells,
    serum_selection,
):
    # Compute neutralization standard counts per well
    neut_standard_counts = (
        counts_qc_4.query("neut_standard")
        .groupby(
            ["well", "serum_replicate", "sample_well", "dilution_factor"],
            dropna=False,
            as_index=False,
        )
        .aggregate(neut_standard_count=pd.NamedAgg("count", "sum"))
        .assign(
            fails_qc=lambda x: (
                x["neut_standard_count"]
                < qc_thresholds["min_neut_standard_count_per_well"]
            )
        )
    )

    # Create chart
    neut_standard_counts_chart = (
        alt.Chart(neut_standard_counts)
        .add_params(serum_selection)
        .transform_filter(serum_selection)
        .encode(
            alt.X(
                "neut_standard_count",
                title="counts from neutralization standard",
                scale=alt.Scale(nice=False, padding=3),
            ),
            alt.Y("sample_well", sort=sample_wells),
            alt.Color(
                "fails_qc",
                title=f"fails qc_thresholds['min_neut_standard_count_per_well']={qc_thresholds['min_neut_standard_count_per_well']!r}",
                legend=alt.Legend(titleLimit=500),
            ),
            tooltip=[
                (
                    alt.Tooltip(c, format=".3g")
                    if neut_standard_counts[c].dtype == float
                    else c
                )
                for c in neut_standard_counts.columns
            ],
        )
        .mark_bar(height={"band": 0.85})
        .properties(
            height=alt.Step(10),
            width=250,
            title=f"Neutralization-standard counts for {plate}",
        )
        .configure_axis(grid=False)
        .configure_legend(titleLimit=1000)
    )

    # Display the chart
    neut_standard_counts_chart
    return (neut_standard_counts,)


@app.cell
def _(counts_qc_4, mo, neut_standard_counts, qc_drops, qc_thresholds):
    # drop wells failing QC
    min_neut_standard_count_per_well_drops = list(
        neut_standard_counts.query("fails_qc")["well"]
    )
    mo.output.append(
        mo.md(
            f"Dropping {len(min_neut_standard_count_per_well_drops)} wells for failing `qc_thresholds['min_neut_standard_count_per_well']={qc_thresholds['min_neut_standard_count_per_well']!r}`: {min_neut_standard_count_per_well_drops}"
        )
    )
    qc_drops["wells"].update(
        {
            w: "min_neut_standard_count_per_well"
            for w in min_neut_standard_count_per_well_drops
        }
    )
    neut_standard_counts_1 = neut_standard_counts[
        ~neut_standard_counts["well"].isin(qc_drops["wells"])
    ]
    counts_qc_5 = counts_qc_4[~counts_qc_4["well"].isin(qc_drops["wells"])]
    return counts_qc_5, neut_standard_counts_1


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Compute and plot the no-serum sample viral barcode counts and check if they pass the QC filters.
    """)
    return


@app.cell
def _(
    alt,
    barcode_selection,
    counts_qc_5,
    neut_standard_counts_1,
    plate,
    qc_thresholds,
    sample_wells,
):
    # Compute no-serum viral barcode counts with QC
    no_serum_counts = (
        counts_qc_5.query("serum == 'none'")
        .query("not neut_standard")
        .merge(neut_standard_counts_1, validate="many_to_one")[
            ["barcode", "strain", "well", "sample_well", "count", "neut_standard_count"]
        ]
        .assign(
            fails_qc=lambda x: (
                x["count"] <= qc_thresholds["min_no_serum_count_per_viral_barcode_well"]
            )
        )
    )

    # Create strain selection dropdown
    strains = sorted(no_serum_counts["strain"].unique())
    strain_selection_dropdown = alt.selection_point(
        fields=["strain"],
        bind=alt.binding_select(
            options=[None] + strains,
            labels=["all"] + strains,
            name="virus strain",
        ),
    )

    # Prepare data for plotting
    no_serum_counts_plot_df = no_serum_counts.drop(
        columns=["well", "neut_standard_count"]
    )

    # Create chart
    no_serum_counts_chart = (
        alt.Chart(no_serum_counts_plot_df)
        .add_params(barcode_selection, strain_selection_dropdown)
        .transform_filter(strain_selection_dropdown)
        .encode(
            alt.X(
                "count",
                title="viral barcode count",
                scale=alt.Scale(nice=False, padding=5),
            ),
            alt.Y("sample_well", sort=sample_wells),
            alt.Fill(
                "fails_qc",
                title=f"fails qc_thresholds['min_no_serum_count_per_viral_barcode_well']={qc_thresholds['min_no_serum_count_per_viral_barcode_well']!r}",
                legend=alt.Legend(titleLimit=500),
            ),
            strokeWidth=alt.condition(barcode_selection, alt.value(2), alt.value(0)),
            size=alt.condition(barcode_selection, alt.value(60), alt.value(35)),
            tooltip=no_serum_counts_plot_df.columns.tolist(),
        )
        .mark_circle(fillOpacity=0.6, stroke="black", strokeOpacity=1)
        .properties(
            height=alt.Step(10),
            width=400,
            title=f"{plate} viral barcode counts in no-serum samples",
        )
        .configure_axis(grid=False)
        .configure_legend(titleLimit=1000)
        .interactive()
    )

    # Display the chart
    no_serum_counts_chart
    return no_serum_counts, strain_selection_dropdown


@app.cell
def _(counts_qc_5, mo, no_serum_counts, qc_drops, qc_thresholds):
    # drop barcode / wells failing QC
    min_no_serum_count_per_viral_barcode_well_drops = list(
        no_serum_counts.query("fails_qc")[["barcode", "well"]].itertuples(
            index=False, name=None
        )
    )
    mo.output.append(
        mo.md(
            f"Dropping {len(min_no_serum_count_per_viral_barcode_well_drops)} barcode-wells for failing `qc_thresholds['min_no_serum_count_per_viral_barcode_well']={qc_thresholds['min_no_serum_count_per_viral_barcode_well']!r}`: {min_no_serum_count_per_viral_barcode_well_drops}"
        )
    )
    qc_drops["barcode_wells"].update(
        {
            w: "min_no_serum_count_per_viral_barcode_well"
            for w in min_no_serum_count_per_viral_barcode_well_drops
        }
    )
    no_serum_counts_1 = no_serum_counts[
        ~no_serum_counts.assign(
            barcode_well=lambda x: x.apply(lambda r: (r["barcode"], r["well"]), axis=1)
        )["barcode_well"].isin(qc_drops["barcode_wells"])
    ]
    counts_qc_6 = counts_qc_5[
        ~counts_qc_5.assign(
            barcode_well=lambda x: x.apply(lambda r: (r["barcode"], r["well"]), axis=1)
        )["barcode_well"].isin(qc_drops["barcode_wells"])
    ]
    return counts_qc_6, no_serum_counts_1


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Compute and plot the median ratio of viral barcode count to neut standard counts across no-serum samples.
    If library composition is equal, all of these values should be similar:
    """)
    return


@app.cell
def _(alt, no_serum_counts_1, pd, plate):
    # Compute median ratio of viral barcode to neut standard counts
    median_no_serum_ratio = (
        no_serum_counts_1.assign(ratio=lambda x: x["count"] / x["neut_standard_count"])
        .groupby(["barcode", "strain"], as_index=False)
        .aggregate(median_no_serum_ratio=pd.NamedAgg("ratio", "median"))
    )

    # Create strain selection for mouseover
    strain_selection = alt.selection_point(
        fields=["strain"],
        on="mouseover",
        empty=False,
    )

    # Create chart
    median_no_serum_ratio_chart = (
        alt.Chart(median_no_serum_ratio)
        .add_params(strain_selection)
        .encode(
            alt.X(
                "median_no_serum_ratio",
                title="median ratio of counts",
                scale=alt.Scale(nice=False, padding=5),
            ),
            alt.Y(
                "barcode",
                sort=alt.SortField("median_no_serum_ratio", order="descending"),
                axis=alt.Axis(labelFontSize=5),
            ),
            color=alt.condition(
                strain_selection,
                alt.value("orange"),
                alt.value("gray"),
            ),
            tooltip=[
                (
                    alt.Tooltip(c, format=".3g")
                    if median_no_serum_ratio[c].dtype == float
                    else c
                )
                for c in median_no_serum_ratio.columns
            ],
        )
        .mark_bar(height={"band": 0.85})
        .properties(
            height=alt.Step(5),
            width=250,
            title=f"{plate} no-serum median ratio viral barcode to neut-standard barcode",
        )
        .configure_axis(grid=False)
        .configure_legend(titleLimit=1000)
    )

    # Display the chart
    median_no_serum_ratio_chart
    return (median_no_serum_ratio,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Compute and plot the actual fraction infectivities.
    We compute both the raw fraction infectivities and the ones with the ceiling applied:
    """)
    return


@app.cell
def _(
    counts_qc_6,
    curvefit_params,
    median_no_serum_ratio,
    neut_standard_counts_1,
):
    frac_infectivity = (
        counts_qc_6.query("not neut_standard")
        .query("serum != 'none'")
        .merge(median_no_serum_ratio, validate="many_to_one")
        .merge(neut_standard_counts_1, validate="many_to_one")
        .assign(
            frac_infectivity_raw=lambda x: (
                x["count"] / x["neut_standard_count"] / x["median_no_serum_ratio"]
            ),
            frac_infectivity_ceiling=lambda x: (
                x["frac_infectivity_raw"].clip(
                    upper=curvefit_params["frac_infectivity_ceiling"]
                )
            ),
            concentration=lambda x: 1 / x["dilution_factor"],
            plate_barcode=lambda x: x["plate_replicate"] + "-" + x["barcode"],
        )[
            [
                "barcode",
                "plate_barcode",
                "well",
                "strain",
                "serum",
                "serum_replicate",
                "dilution_factor",
                "concentration",
                "frac_infectivity_raw",
                "frac_infectivity_ceiling",
            ]
        ]
    )
    assert len(
        frac_infectivity.groupby(["serum", "plate_barcode", "dilution_factor"])
    ) == len(frac_infectivity)
    assert frac_infectivity["dilution_factor"].notnull().all()
    assert frac_infectivity["frac_infectivity_raw"].notnull().all()
    assert frac_infectivity["frac_infectivity_ceiling"].notnull().all()
    return (frac_infectivity,)


@app.cell
def _(curvefit_params, frac_infectivity, qc_thresholds):
    frac_infectivity_cols = {
        "frac_infectivity_raw": "raw fraction infectivity",
        "frac_infectivity_ceiling": f"fraction infectivity with ceiling at {curvefit_params['frac_infectivity_ceiling']}",
    }

    frac_infectivity_chart_df = frac_infectivity.assign(
        fails_qc=lambda x: (
            x["frac_infectivity_raw"]
            > qc_thresholds["max_frac_infectivity_per_viral_barcode_well"]
        ),
    )[
        [
            "barcode",
            "strain",
            "well",
            "serum_replicate",
            "dilution_factor",
            "fails_qc",
            *list(frac_infectivity_cols),
        ]
    ].rename(
        columns=frac_infectivity_cols
    )

    # some manipulations to shrink data frame plotted by altair below by putting
    # them in smaller data frames that are used via transform_lookup
    barcode_lookup_df = frac_infectivity[["barcode", "strain"]].drop_duplicates()
    assert len(barcode_lookup_df) == barcode_lookup_df["barcode"].nunique()
    well_lookup_df = frac_infectivity[
        ["well", "serum_replicate", "dilution_factor"]
    ].drop_duplicates()
    assert len(well_lookup_df) == well_lookup_df["well"].nunique()

    frac_infectivity_chart_df = frac_infectivity_chart_df.drop(
        columns=["strain", "serum_replicate", "dilution_factor"]
    )
    return (
        barcode_lookup_df,
        frac_infectivity_chart_df,
        frac_infectivity_cols,
        well_lookup_df,
    )


@app.cell
def _(
    alt,
    barcode_lookup_df,
    barcode_selection,
    frac_infectivity_chart_df,
    frac_infectivity_cols,
    mo,
    plate,
    qc_thresholds,
    strain_selection_dropdown,
    well_lookup_df,
):
    frac_infectivity_chart = (
        alt.Chart(frac_infectivity_chart_df)
        .transform_lookup(
            lookup="barcode",
            from_=alt.LookupData(barcode_lookup_df, key="barcode", fields=["strain"]),
        )
        .transform_lookup(
            lookup="well",
            from_=alt.LookupData(
                well_lookup_df,
                key="well",
                fields=["serum_replicate", "dilution_factor"],
            ),
        )
        .transform_fold(
            frac_infectivity_cols.values(), ["ceiling_applied", "frac_infectivity"]
        )
        .add_params(strain_selection_dropdown, barcode_selection)
        .transform_filter(strain_selection_dropdown)
        .encode(
            alt.X(
                "dilution_factor:Q",
                title="dilution factor",
                scale=alt.Scale(nice=False, padding=5, type="log"),
            ),
            alt.Y(
                "frac_infectivity:Q",
                title="fraction infectivity",
                scale=alt.Scale(nice=False, padding=5),
            ),
            alt.Column(
                "ceiling_applied:N",
                sort="descending",
                title=None,
                header=alt.Header(
                    labelFontSize=13, labelFontStyle="bold", labelPadding=2
                ),
            ),
            alt.Row(
                "serum_replicate:N",
                title=None,
                spacing=3,
                header=alt.Header(labelFontSize=13, labelFontStyle="bold"),
            ),
            alt.Detail("barcode"),
            alt.Shape(
                "fails_qc",
                title=f"fails {qc_thresholds['max_frac_infectivity_per_viral_barcode_well']=}",
                legend=alt.Legend(titleLimit=500, orient="bottom"),
            ),
            color=alt.condition(
                barcode_selection, alt.value("black"), alt.value("MediumBlue")
            ),
            strokeWidth=alt.condition(barcode_selection, alt.value(3), alt.value(1)),
            opacity=alt.condition(barcode_selection, alt.value(1), alt.value(0.25)),
            tooltip=[
                (
                    alt.Tooltip(c, format=".3g")
                    if frac_infectivity_chart_df[c].dtype == float
                    else c
                )
                for c in frac_infectivity_chart_df.columns
            ]
            + [
                alt.Tooltip("strain:N"),
                alt.Tooltip("serum_replicate:N"),
                alt.Tooltip("dilution_factor:Q"),
            ],
        )
        .mark_line(point=True)
        .properties(
            height=150,
            width=250,
            title=f"Fraction infectivities for {plate}",
        )
        .interactive(bind_x=False)
        .configure_axis(grid=False)
        .configure_legend(titleLimit=1000)
        .configure_point(size=50)
        .resolve_scale(x="independent", y="independent")
    )

    # Display the chart
    mo.output.append(frac_infectivity_chart)
    return


@app.cell
def _(
    frac_infectivity,
    frac_infectivity_chart_df,
    mo,
    qc_drops,
    qc_thresholds,
):
    # drop barcode / wells failing QC
    max_frac_infectivity_per_viral_barcode_well_drops = list(
        frac_infectivity_chart_df.query("fails_qc")[["barcode", "well"]]
        .drop_duplicates()
        .itertuples(index=False, name=None)
    )
    mo.output.append(
        mo.md(
            f"Dropping {len(max_frac_infectivity_per_viral_barcode_well_drops)} barcode-wells for failing `qc_thresholds['max_frac_infectivity_per_viral_barcode_well']={qc_thresholds['max_frac_infectivity_per_viral_barcode_well']!r}`: {max_frac_infectivity_per_viral_barcode_well_drops}"
        )
    )
    qc_drops["barcode_wells"].update(
        {
            w: "max_frac_infectivity_per_viral_barcode_well"
            for w in max_frac_infectivity_per_viral_barcode_well_drops
        }
    )
    frac_infectivity_1 = frac_infectivity[
        ~frac_infectivity.assign(
            barcode_well=lambda x: x.apply(lambda r: (r["barcode"], r["well"]), axis=1)
        )["barcode_well"].isin(qc_drops["barcode_wells"])
    ]
    return (frac_infectivity_1,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Check how many dilutions we have per barcode / serum-replicate:
    """)
    return


@app.cell
def _(
    alt,
    barcode_selection,
    frac_infectivity_1,
    mo,
    pd,
    qc_drops,
    qc_thresholds,
):
    # Count number of dilutions per barcode/serum-replicate
    n_dilutions = (
        frac_infectivity_1.groupby(
            ["serum_replicate", "strain", "barcode"],
            as_index=False,
        )
        .aggregate(**{"number of dilutions": pd.NamedAgg("dilution_factor", "nunique")})
        .assign(
            fails_qc=lambda x: (
                x["number of dilutions"]
                < qc_thresholds["min_dilutions_per_barcode_serum_replicate"]
            )
        )
    )

    # Create chart
    n_dilutions_chart = (
        alt.Chart(n_dilutions)
        .add_params(barcode_selection)
        .encode(
            alt.X(
                "number of dilutions",
                scale=alt.Scale(nice=False, padding=4),
            ),
            alt.Y("strain", title=None),
            alt.Column(
                "serum_replicate",
                title=None,
                header=alt.Header(
                    labelFontSize=12,
                    labelFontStyle="bold",
                    labelPadding=0,
                ),
            ),
            alt.Fill(
                "fails_qc",
                title=f"fails qc_thresholds['min_dilutions_per_barcode_serum_replicate']={qc_thresholds['min_dilutions_per_barcode_serum_replicate']!r}",
                legend=alt.Legend(titleLimit=500, orient="bottom"),
            ),
            strokeWidth=alt.condition(barcode_selection, alt.value(2), alt.value(0)),
            size=alt.condition(barcode_selection, alt.value(55), alt.value(35)),
            tooltip=[
                alt.Tooltip(c, format=".3g") if n_dilutions[c].dtype == float else c
                for c in n_dilutions.columns
            ],
        )
        .mark_circle(stroke="black", strokeOpacity=1, fillOpacity=0.45)
        .properties(
            height=alt.Step(10),
            width=120,
            title=alt.TitleParams(
                "number of dilutions for each barcode for each serum-replicate",
                dy=-2,
            ),
        )
    )

    # Display the chart
    mo.output.append(n_dilutions_chart)

    # Drop barcode/serum-replicates failing QC
    min_dilutions_per_barcode_serum_replicate_drops = list(
        n_dilutions.query("fails_qc")[["barcode", "serum_replicate"]].itertuples(
            index=False, name=None
        )
    )
    mo.output.append(
        mo.md(
            f"Dropping {len(min_dilutions_per_barcode_serum_replicate_drops)} barcode/serum-replicates for failing `qc_thresholds['min_dilutions_per_barcode_serum_replicate']={qc_thresholds['min_dilutions_per_barcode_serum_replicate']!r}`: {min_dilutions_per_barcode_serum_replicate_drops}"
        )
    )
    qc_drops["barcode_serum_replicates"].update(
        {
            w: "min_dilutions_per_barcode_serum_replicate"
            for w in min_dilutions_per_barcode_serum_replicate_drops
        }
    )
    frac_infectivity_2 = frac_infectivity_1[
        ~frac_infectivity_1.assign(
            barcode_serum_replicate=lambda x: x.apply(
                lambda r: (r["barcode"], r["serum_replicate"]), axis=1
            )
        )["barcode_serum_replicate"].isin(qc_drops["barcode_serum_replicates"])
    ]
    return (frac_infectivity_2,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Fit neutralization curves without applying QC to curves
    First fit curves to all serum replicates, then we will apply QC on the curve fits.
    Note that the fitting is done to the fraction infectivities **with** the ceiling.
    """)
    return


@app.cell
def _(curvefit_params, frac_infectivity_2, neutcurve):
    fits_noqc = neutcurve.CurveFits(
        frac_infectivity_2.rename(
            columns={
                "frac_infectivity_ceiling": "fraction infectivity",
                "concentration": "serum concentration",
            }
        ),
        conc_col="serum concentration",
        fracinf_col="fraction infectivity",
        virus_col="strain",
        serum_col="serum_replicate",
        replicate_col="barcode",
        fixtop=curvefit_params["fixtop"],
        fixbottom=curvefit_params["fixbottom"],
        fixslope=curvefit_params["fixslope"],
    )
    return (fits_noqc,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Determine which fits fail the curve fitting QC, and plot them.
    Note the plot indicates as failing QC any barcode / serum-replicate that fails, even if we are also specified to ignore the QC for that one (so it will not be removed later):
    """)
    return


@app.cell
def _(curvefit_qc, fits_noqc, frac_infectivity_2, pd):
    goodness_of_fit = curvefit_qc["goodness_of_fit"]
    fit_params_noqc = (
        frac_infectivity_2.groupby(["serum_replicate", "barcode"], as_index=False)
        .aggregate(max_frac_infectivity=pd.NamedAgg("frac_infectivity_ceiling", "max"))
        .merge(
            fits_noqc.fitParams(average_only=False, no_average=True)[
                ["serum", "virus", "replicate", "r2", "rmsd"]
            ].rename(columns={"serum": "serum_replicate", "replicate": "barcode"}),
            validate="one_to_one",
        )
        .assign(
            fails_max_frac_infectivity_at_least=lambda x: (
                x["max_frac_infectivity"] < curvefit_qc["max_frac_infectivity_at_least"]
            ),
            fails_goodness_of_fit=lambda x: (
                (x["r2"] < goodness_of_fit["min_R2"])
                & (x["rmsd"] > goodness_of_fit["max_RMSD"])
            ),
            fails_qc=lambda x: (
                x["fails_max_frac_infectivity_at_least"] | x["fails_goodness_of_fit"]
            ),
            ignore_qc=lambda x: x.apply(
                lambda r: (
                    r["serum_replicate"]
                    in curvefit_qc["serum_replicates_ignore_curvefit_qc"]
                    or (r["barcode"], r["serum_replicate"])
                    in curvefit_qc["barcode_serum_replicates_ignore_curvefit_qc"]
                ),
                axis=1,
            ),
        )
    )
    return (fit_params_noqc,)


@app.cell
def _(alt, barcode_selection, curvefit_qc, fit_params_noqc, mo):
    mo.output.append(
        mo.md(
            f"Plotting barcode / serum-replicates that fail `curvefit_qc={curvefit_qc!r}`"
        )
    )

    fit_params_noqc_base_chart = alt.Chart(fit_params_noqc).add_params(
        barcode_selection
    )
    for prop, col in [
        ("max frac infectivity", "max_frac_infectivity"),
        ("curve fit R2", "r2"),
        ("curve fit RMSD", "rmsd"),
    ]:
        chart = (
            fit_params_noqc_base_chart.encode(
                alt.X(col, title=prop, scale=alt.Scale(nice=False, padding=4)),
                alt.Y("virus", title=None),
                alt.Fill("fails_qc"),
                alt.Column(
                    "serum_replicate",
                    title=None,
                    header=alt.Header(
                        labelFontSize=12, labelFontStyle="bold", labelPadding=0
                    ),
                ),
                strokeWidth=alt.condition(
                    barcode_selection, alt.value(2), alt.value(0)
                ),
                size=alt.condition(barcode_selection, alt.value(55), alt.value(35)),
                tooltip=[
                    (
                        alt.Tooltip(c, format=".3g")
                        if fit_params_noqc[c].dtype == float
                        else c
                    )
                    for c in fit_params_noqc.columns
                ],
            )
            .mark_circle(stroke="black", strokeOpacity=1, fillOpacity=0.55)
            .properties(
                height=alt.Step(10),
                width=90,
                title=alt.TitleParams(
                    f"{prop} for each barcode serum-replicate", dy=-2
                ),
            )
        )

        # Display the chart
        mo.output.append(chart)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Now plot curves for all virus vs serum-replicates that have a barcode that fails any of the QC.
    In these plots, the suffix on the barcode name in the color key indicates if it passed or failed QC:
    """)
    return


@app.cell
def _(
    CBMARKERS,
    CBPALETTE,
    curvefit_qc,
    fit_params_noqc,
    fits_noqc,
    mo,
    mo_pdf_from_fig,
    pd,
):
    barcode_serum_replicates_fail_qc = fit_params_noqc.query("fails_qc").reset_index(
        drop=True
    )
    mo.output.append(
        mo.md(
            f"Here are barcode / serum-replicates that fail `curvefit_qc={curvefit_qc!r}`"
        )
    )
    mo.output.append(barcode_serum_replicates_fail_qc)

    if len(barcode_serum_replicates_fail_qc):
        mo.output.append(
            mo.md(
                "Curves for virus vs serum-replicates with at least one failed barcode.\n\nColor key labels indicate if barcodes failed or passed QC."
            )
        )
        plots = {}
        ncol = 6
        for iplot, (serum, virus, failed_barcodes) in enumerate(
            barcode_serum_replicates_fail_qc.groupby(
                ["serum_replicate", "virus"], as_index=False
            )
            .aggregate(barcodes=pd.NamedAgg("barcode", list))
            .itertuples(index=False)
        ):
            passed_barcodes = [
                bc
                for bc in fits_noqc.replicates[serum, virus]
                if bc not in failed_barcodes and bc != "average"
            ]
            curvelist = []
            assert len(CBMARKERS) >= len(failed_barcodes + passed_barcodes)
            assert len(CBPALETTE) >= len(failed_barcodes + passed_barcodes)
            for replicate, marker, color in zip(
                failed_barcodes + passed_barcodes, CBMARKERS, CBPALETTE
            ):
                curvelist.append(
                    {
                        "serum": serum,
                        "virus": virus,
                        "replicate": replicate,
                        "label": replicate
                        + ("-fail" if replicate in failed_barcodes else "-pass"),
                        "color": color,
                        "marker": marker,
                    }
                )
            plots[iplot // ncol, iplot % ncol] = (f"{serum} vs {virus}", curvelist)
        _fig_fail_qc, _ = fits_noqc.plotGrid(
            plots,
            attempt_shared_legend=False,
            legendfontsize=8,
            titlesize=9,
            ticksize=10,
            draw_in_bounds=True,
        )
        mo.output.append(mo_pdf_from_fig(_fig_fail_qc))
    else:
        mo.output.append(mo.md("No serum-replicates fail QC."))
    return


@app.cell
def _(curvefit_qc, fit_params_noqc, frac_infectivity_2, mo, qc_drops):
    # drop barcode / serum-replicates failing QC
    frac_infectivity_3 = frac_infectivity_2.copy()
    fit_params_noqc_1 = fit_params_noqc.copy()
    for qc_filter in ["max_frac_infectivity_at_least", "goodness_of_fit"]:
        fits_qc_drops = list(
            fit_params_noqc_1.query(f"fails_{qc_filter} and (not ignore_qc)")[
                ["barcode", "serum_replicate"]
            ].itertuples(index=False, name=None)
        )
        mo.output.append(
            mo.md(
                f"Dropping {len(fits_qc_drops)} barcode/serum-replicates for failing `{qc_filter}={curvefit_qc[qc_filter]}`: {fits_qc_drops}"
            )
        )
        qc_drops["barcode_serum_replicates"].update(
            {w: qc_filter for w in fits_qc_drops}
        )
        frac_infectivity_3 = frac_infectivity_3[
            ~frac_infectivity_3.assign(
                barcode_serum_replicate=lambda x: x.apply(
                    lambda r: (r["barcode"], r["serum_replicate"]), axis=1
                )
            )["barcode_serum_replicate"].isin(qc_drops["barcode_serum_replicates"])
        ]
        fit_params_noqc_1 = fit_params_noqc_1[
            ~fit_params_noqc_1.assign(
                barcode_serum_replicate=lambda x: x.apply(
                    lambda r: (r["barcode"], r["serum_replicate"]), axis=1
                )
            )["barcode_serum_replicate"].isin(qc_drops["barcode_serum_replicates"])
        ]
    return (frac_infectivity_3,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Fit neutralization curves after applying QC
    No we re-fit and plot curves after applying all the QC.
    """)
    return


@app.cell
def _(curvefit_params, fits_noqc, frac_infectivity_3, group, mo, neutcurve):
    fits_qc = neutcurve.CurveFits(
        frac_infectivity_3.rename(
            columns={
                "frac_infectivity_ceiling": "fraction infectivity",
                "concentration": "serum concentration",
            }
        ),
        conc_col="serum concentration",
        fracinf_col="fraction infectivity",
        virus_col="strain",
        serum_col="serum",
        replicate_col="plate_barcode",
        fixtop=curvefit_params["fixtop"],
        fixbottom=curvefit_params["fixbottom"],
        fixslope=curvefit_params["fixslope"],
    )
    fit_params_qc = fits_qc.fitParams(average_only=False, no_average=True)
    assert len(fit_params_qc) <= len(
        fits_noqc.fitParams(average_only=False, no_average=True)
    )
    mo.output.append(mo.md(f"Assigning fits for this plate to `{group}`"))
    fit_params_qc.insert(0, "group", group)
    return fit_params_qc, fits_qc


@app.cell
def _(fits_qc, mo, mo_pdf_from_fig):
    if fits_qc.sera:
        _fig_passed_qc, _ = fits_qc.plotReplicates(
            attempt_shared_legend=False,
            legendfontsize=8,
            titlesize=9,
            ticksize=10,
            ncol=6,
            draw_in_bounds=True,
        )
        mo.output.append(mo_pdf_from_fig(_fig_passed_qc))
    else:
        mo.output.append(mo.md("No sera passed QC."))
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Save results to files
    """)
    return


@app.cell
def _(
    fit_params_qc,
    fits_csv,
    fits_pickle,
    fits_qc,
    frac_infectivity_3,
    frac_infectivity_csv,
    io,
    mo,
    pickle,
    qc_drops,
    qc_drops_yaml,
    yaml,
):
    mo.output.append(
        mo.md(f"Writing fraction infectivities to `{frac_infectivity_csv}`")
    )
    (
        frac_infectivity_3[
            [
                "serum",
                "strain",
                "plate_barcode",
                "dilution_factor",
                "frac_infectivity_raw",
                "frac_infectivity_ceiling",
            ]
        ]
        .sort_values(["serum", "plate_barcode", "dilution_factor"])
        .to_csv(frac_infectivity_csv, index=False, float_format="%.4g")
    )
    mo.output.append(mo.md(f"Writing fit parameters to `{fits_csv}`"))
    fit_params_qc.drop(columns=["nreplicates", "ic50_str"]).to_csv(
        fits_csv, index=False, float_format="%.4g"
    )
    mo.output.append(
        mo.md(f"Pickling neutcurve.CurveFits object for these data to `{fits_pickle}`")
    )
    with open(fits_pickle, "wb") as f_pickle:
        pickle.dump(fits_qc, f_pickle)
    mo.output.append(mo.md(f"Writing QC drops to `{qc_drops_yaml}`"))

    def tup_to_str(x):
        return " ".join(x) if isinstance(x, tuple) else x

    qc_drops_for_yaml = {
        key: {tup_to_str(key2): val2 for key2, val2 in val.items()}
        for key, val in qc_drops.items()
    }
    with open(qc_drops_yaml, "w") as f_yaml:
        yaml.YAML(typ="rt").dump(qc_drops_for_yaml, f_yaml)
    mo.output.append(mo.md("Here are the QC drops:"))
    yaml_buffer_qc_drops = io.StringIO()
    yaml.YAML(typ="rt").dump(qc_drops_for_yaml, stream=yaml_buffer_qc_drops)
    mo.output.append(mo.md(f"```yaml\n{yaml_buffer_qc_drops.getvalue()}```"))
    return


if __name__ == "__main__":
    app.run()

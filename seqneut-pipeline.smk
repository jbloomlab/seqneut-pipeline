"""Snakemake file for the ``seqneut-pipeline`.

Designed to be included in another ``Snakefile`` that specifies the config.

"""

import pandas as pd


snakemake.utils.min_version("9.0")


include: "funcs.smk"  # include functions


# --- Process configuration ------------------------------------------------------------

pipeline_subdir = config["seqneut-pipeline"]

viral_libraries = {
    lib: pd.read_csv(f) for (lib, f) in config["viral_libraries"].items()
}

if ("viral_strain_plot_order" not in config) or (
    config["viral_strain_plot_order"] is None
):
    viral_strain_plot_order = None
else:
    viral_strain_plot_order = pd.read_csv(config["viral_strain_plot_order"])[
        "strain"
    ].tolist()
    assert len(viral_strain_plot_order) == len(set(viral_strain_plot_order))

neut_standard_sets = {
    s: pd.read_csv(f) for (s, f) in config["neut_standard_sets"].items()
}

plates = {
    str(plate): process_plate(str(plate), plate_params)
    for (plate, plate_params) in config["plates"].items()
}

groups = sorted(set(plate_params["group"] for plate_params in plates.values()))
groups_cannot_contain = ["|", "_"]  # wildcard problems if group contains these
if any(s in group for s in groups_cannot_contain for group in groups):
    raise ValueError(f"found {groups_cannot_contain=} character in {groups=}")


wildcard_constraints:
    group="|".join(groups),


if not set(config["sera_override_defaults"]).issubset(groups):
    raise ValueError(f"{config['sera_override_defaults']=} keyed by invalid groups")


if plates == {}:
    samples = {}
else:
    samples = pd.concat(
        [plate_d["samples"] for plate_d in plates.values()],
        ignore_index=True,
    )
    assert samples["sample"].nunique() == len(samples)
    samples = samples.set_index("sample").to_dict(orient="index")


if "miscellaneous_plates" in config:
    miscellaneous_plates = process_miscellaneous_plates(config["miscellaneous_plates"])
else:
    miscellaneous_plates = {}

if "curve_display_method" in config:
    curve_display_method = config["curve_display_method"]
else:
    curve_display_method = "png8"

# define `add_htmls_to_docs` if not already defined.
try:
    add_htmls_to_docs
except NameError:  # if not defined
    add_htmls_to_docs = {}


# --- Snakemake rules -------------------------------------------------------------------

if plates:

    rule count_barcodes:
        """Count barcodes for a sample."""
        input:
            fastq=lambda wc: samples[wc.sample]["fastq"],
        output:
            counts="results/barcode_counts/{sample}.csv",
            invalid="results/barcode_invalid/{sample}.csv",
            fates="results/barcode_fates/{sample}.csv",
        params:
            viral_barcodes=lambda wc: sorted(
                viral_libraries[plates[samples[wc.sample]["plate"]]["viral_library"]][
                    "barcode"
                ]
            ),
            neut_standard_barcodes=lambda wc: sorted(
                neut_standard_sets[
                    plates[samples[wc.sample]["plate"]]["neut_standard_set"]
                ]["barcode"]
            ),
            illumina_barcode_parser_params=lambda wc: plates[
                samples[wc.sample]["plate"]
            ]["illumina_barcode_parser_params"],
        conda:
            "envs/count_barcodes.yml"
        log:
            "results/logs/count_barcodes_{sample}.txt",
        script:
            "scripts/count_barcodes.py"

    rule process_plate:
        """Process a plate to QC and convert counts to fraction infectivity."""
        input:
            marimo_nb=os.path.join(pipeline_subdir, "notebooks/process_plate.py"),
            count_csvs=lambda wc: expand(
                rules.count_barcodes.output.counts,
                sample=plates[wc.plate]["samples"]["sample"],
            ),
            fate_csvs=lambda wc: expand(
                rules.count_barcodes.output.fates,
                sample=plates[wc.plate]["samples"]["sample"],
            ),
        output:
            marimo_html="results/plates/{plate}/process_{plate}.html",
            context_pickle="results/plates/{plate}/process_{plate}_context.pickle",
            qc_drops="results/plates/{plate}/qc_drops.yml",
            frac_infectivity_csv="results/plates/{plate}/frac_infectivity.csv",
            fits_csv="results/plates/{plate}/curvefits.csv",
            fits_pickle="results/plates/{plate}/curvefits.pickle",
        log:
            "results/logs/process_{plate}.txt",
        params:
            # pass DataFrames/Series as dict/list for snakemake params rerun triggers
            viral_barcodes=lambda wc: (
                viral_libraries[plates[wc.plate]["viral_library"]][
                    ["barcode", "strain"]
                ]
                .sort_values("barcode")
                .reset_index(drop=True)
                .to_dict()
            ),
            neut_standard_barcodes=lambda wc: (
                neut_standard_sets[plates[wc.plate]["neut_standard_set"]][["barcode"]]
                .sort_values("barcode")
                .reset_index(drop=True)
                .to_dict()
            ),
            samples=lambda wc: plates[wc.plate]["samples"]["sample"].tolist(),
            plate_params=lambda wc: {
                param: (val if param != "samples" else val.to_dict())
                for (param, val) in plates[wc.plate].items()
            },
            curve_display_method=curve_display_method,
        conda:
            "environment.yml"
        script:
            "scripts/run_marimo_w_context_pickle.py"

    checkpoint groups_sera_by_plate:
        """Get list of all groups/sera and plates they are on."""
        input:
            csvs=expand(rules.process_plate.output.fits_csv, plate=plates),
        output:
            csv="results/sera/groups_sera_by_plate.csv",
        params:
            plates=list(plates),
        log:
            "results/logs/groups_sera_by_plate.txt",
        conda:
            "environment.yml"
        script:
            "scripts/groups_sera_by_plate.py"

    rule group_serum_titers:
        """Aggregate and analyze titers for a group / serum."""
        input:
            marimo_nb=os.path.join(pipeline_subdir, "notebooks/group_serum_titers.py"),
            pickles=lambda wc: [
                rules.process_plate.output.fits_pickle.format(plate=plate)
                for plate in groups_sera_plates()[(wc.group, wc.serum)]
            ],
        output:
            marimo_html="results/sera/{group}_{serum}/{group}_{serum}_titers.html",
            context_pickle="results/sera/{group}_{serum}/{group}_{serum}_titers_context.pickle",
            per_rep_titers="results/sera/{group}_{serum}/titers_per_replicate.csv",
            titers="results/sera/{group}_{serum}/titers.csv",
            curves_pdf="results/sera/{group}_{serum}/curves.pdf",
            pickle="results/sera/{group}_{serum}/curvefits.pickle",
            qc_drops="results/sera/{group}_{serum}/qc_drops.yml",
        params:
            viral_strain_plot_order=viral_strain_plot_order,
            serum_titer_as=lambda wc: (
                config["sera_override_defaults"][wc.group][wc.serum]["titer_as"]
                if (
                    (wc.group in config["sera_override_defaults"])
                    and (wc.serum in config["sera_override_defaults"][wc.group])
                    and (
                        "titer_as"
                        in config["sera_override_defaults"][wc.group][wc.serum]
                    )
                )
                else config["default_serum_titer_as"]
            ),
            qc_thresholds=lambda wc: (
                config["sera_override_defaults"][wc.group][wc.serum]["qc_thresholds"]
                if (
                    (wc.group in config["sera_override_defaults"])
                    and (wc.serum in config["sera_override_defaults"][wc.group])
                    and (
                        "qc_thresholds"
                        in config["sera_override_defaults"][wc.group][wc.serum]
                    )
                )
                else config["default_serum_qc_thresholds"]
            ),
            curve_display_method=curve_display_method,
        log:
            "results/logs/group_serum_titers_{group}_{serum}.txt",
        conda:
            "environment.yml"
        script:
            "scripts/run_marimo_w_context_pickle.py"

    rule aggregate_titers:
        """Aggregate all serum titers."""
        input:
            marimo_nb=os.path.join(pipeline_subdir, "notebooks/aggregate_titers.py"),
            pickles=lambda wc: [
                rules.group_serum_titers.output.pickle.format(group=group, serum=serum)
                for (group, serum) in groups_sera_plates()
            ],
            titers=lambda wc: [
                rules.group_serum_titers.output.titers.format(group=group, serum=serum)
                for (group, serum) in groups_sera_plates()
            ],
        output:
            marimo_html="results/aggregated_titers/aggregate_titers.html",
            context_pickle="results/aggregated_titers/aggregate_titers_context.pickle",
            pickles=[
                f"results/aggregated_titers/curvefits_{group}.pickle"
                for group in groups
            ],
            titers=[f"results/aggregated_titers/titers_{group}.csv" for group in groups],
            titers_chart="results/aggregated_titers/titers.html",
        params:
            viral_strain_plot_order=viral_strain_plot_order,
            groups_sera=lambda wc: list(groups_sera_plates()),
            groups=groups,
        conda:
            "environment.yml"
        log:
            "results/logs/aggregate_titers.txt",
        script:
            "scripts/run_marimo_w_context_pickle.py"

    rule aggregate_qc_drops:
        """Aggregate all QC drops."""
        input:
            marimo_nb=os.path.join(pipeline_subdir, "notebooks/aggregate_qc_drops.py"),
            plate_qc_drops=expand(rules.process_plate.output.qc_drops, plate=plates),
            groups_sera_qc_drops=lambda wc: [
                rules.group_serum_titers.output.qc_drops.format(
                    group=group, serum=serum
                )
                for (group, serum) in groups_sera_plates()
            ],
        output:
            marimo_html="results/qc_drops/aggregate_qc_drops.html",
            context_pickle="results/qc_drops/aggregate_qc_drops_context.pickle",
            plate_qc_drops="results/qc_drops/plate_qc_drops.yml",
            barcode_qc_drops="results/qc_drops/barcode_qc_drops.yml",
            groups_sera_qc_drops="results/qc_drops/groups_sera_qc_drops.yml",
        params:
            plates=list(plates),
            groups_sera=lambda wc: list(groups_sera_plates()),
        conda:
            "environment.yml"
        log:
            "results/logs/aggregate_qc_drops.txt",
        script:
            "scripts/run_marimo_w_context_pickle.py"

    rule build_docs:
        """Build the HTML documentation."""
        input:
            lambda wc: [f for d in add_htmls_to_docs.values() for f in d.values()],
            titers_chart=rules.aggregate_titers.output.titers_chart,
            serum_titers_htmls=lambda wc: [
                f"results/sera/{group}_{serum}/{group}_{serum}_titers.html"
                for (group, serum) in groups_sera_plates()
            ],
            process_plates_htmls=expand(
                "results/plates/{plate}/process_{plate}.html",
                plate=plates,
            ),
            qc_drops_html="results/qc_drops/aggregate_qc_drops.html",
        output:
            docs=directory(config["docs"]),
        params:
            description=config["description"],
            groups_sera=lambda wc: list(groups_sera_plates()),
            plates={plate: plates[plate]["group"] for plate in plates},
            add_htmls_to_docs=lambda wc: {
                key: {key2: str(val2) for (key2, val2) in val.items()}
                for (key, val) in add_htmls_to_docs.items()
            },
        conda:
            "environment.yml"
        log:
            "results/logs/build_docs.txt",
        script:
            "scripts/build_docs.py"


rule miscellaneous_plate_count_barcodes:
    """Count barcodes for a well in a miscellaneous plate."""
    input:
        fastq=lambda wc: miscellaneous_plates[wc.misc_plate]["wells"][wc.well],
    output:
        counts="results/miscellaneous_plates/{misc_plate}/{well}_counts.csv",
        invalid="results/miscellaneous_plates/{misc_plate}/{well}_invalid.csv",
        fates="results/miscellaneous_plates/{misc_plate}/{well}_fates.csv",
    params:
        viral_barcodes=lambda wc: sorted(
            viral_libraries[miscellaneous_plates[wc.misc_plate]["viral_library"]][
                "barcode"
            ]
        ),
        neut_standard_barcodes=lambda wc: sorted(
            neut_standard_sets[
                miscellaneous_plates[wc.misc_plate]["neut_standard_set"]
            ]["barcode"]
        ),
        illumina_barcode_parser_params=lambda wc: miscellaneous_plates[wc.misc_plate][
            "illumina_barcode_parser_params"
        ],
    conda:
        "envs/count_barcodes.yml"
    log:
        "results/logs/miscellaneous_plate_count_barcodes_{misc_plate}_{well}.txt",
    script:
        "scripts/count_barcodes.py"


if plates:
    seqneut_pipeline_outputs = [
        rules.aggregate_titers.output.titers,
        rules.aggregate_titers.output.pickles,
        rules.aggregate_qc_drops.output.plate_qc_drops,
        rules.aggregate_qc_drops.output.barcode_qc_drops,
        rules.aggregate_qc_drops.output.groups_sera_qc_drops,
        rules.build_docs.output.docs,
        *[
            f"results/miscellaneous_plates/{plate}/{well}_{suffix}"
            for plate in miscellaneous_plates
            for well in miscellaneous_plates[plate]["wells"]
            for suffix in ["counts.csv", "fates.csv"]
        ],
    ]

else:
    seqneut_pipeline_outputs = [
        *[
            f"results/miscellaneous_plates/{plate}/{well}_{suffix}"
            for plate in miscellaneous_plates
            for well in miscellaneous_plates[plate]["wells"]
            for suffix in ["counts.csv", "fates.csv"]
        ],
    ]

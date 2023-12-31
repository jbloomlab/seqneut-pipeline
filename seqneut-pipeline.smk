"""Snakemake file for the ``seqneut-pipeline`.

Designed to be included in another ``Snakefile`` that specifies the config.

"""


import pandas as pd


snakemake.utils.min_version("8.0")


include: "funcs.smk"  # include functions


# --- Process configuration ------------------------------------------------------------

pipeline_subdir = config["seqneut-pipeline"]

viral_libraries = config["viral_libraries"]

if ("viral_strain_plot_order" not in config) or (
    config["viral_strain_plot_order"] is None
):
    viral_strain_plot_order = None
else:
    viral_strain_plot_order = pd.read_csv(config["viral_strain_plot_order"])[
        "strain"
    ].tolist()
    assert len(viral_strain_plot_order) == len(set(viral_strain_plot_order))

neut_standard_sets = config["neut_standard_sets"]

plates = {
    str(plate): process_plate(str(plate), plate_params)
    for (plate, plate_params) in config["plates"].items()
}

samples = pd.concat(
    [plate_d["samples"] for plate_d in plates.values()],
    ignore_index=True,
)
assert samples["sample"].nunique() == samples["fastq"].nunique() == len(samples)
samples = samples.set_index("sample").to_dict(orient="index")


# --- Snakemake rules -------------------------------------------------------------------


rule count_barcodes:
    """Count barcodes for a sample."""
    input:
        fastq=lambda wc: samples[wc.sample]["fastq"],
        viral_library=lambda wc: (
            viral_libraries[plates[samples[wc.sample]["plate"]]["viral_library"]]
        ),
        neut_standard_set=lambda wc: (
            neut_standard_sets[
                plates[samples[wc.sample]["plate"]]["neut_standard_set"]
            ]
        ),
    output:
        counts="results/barcode_counts/{sample}.csv",
        invalid="results/barcode_invalid/{sample}.csv",
        fates="results/barcode_fates/{sample}.csv",
    params:
        illumina_barcode_parser_params=config["illumina_barcode_parser_params"],
    conda:
        "envs/count_barcodes.yml"
    log:
        "results/logs/count_barcodes_{sample}.txt",
    script:
        "scripts/count_barcodes.py"


rule process_plate:
    """Process a plate to QC and convert counts to fraction infectivity."""
    input:
        count_csvs=lambda wc: expand(
            rules.count_barcodes.output.counts,
            sample=plates[wc.plate]["samples"]["sample"],
        ),
        fate_csvs=lambda wc: expand(
            rules.count_barcodes.output.fates,
            sample=plates[wc.plate]["samples"]["sample"],
        ),
        viral_library_csv=lambda wc: (
            viral_libraries[plates[wc.plate]["viral_library"]]
        ),
        neut_standard_set_csv=lambda wc: (
            neut_standard_sets[plates[wc.plate]["neut_standard_set"]]
        ),
    output:
        qc_drops="results/plates/{plate}/qc_drops.yml",
        frac_infectivity_csv="results/plates/{plate}/frac_infectivity.csv",
        fits_csv="results/plates/{plate}/curvefits.csv",
        fits_pickle="results/plates/{plate}/curvefits.pickle",
    log:
        notebook="results/plates/{plate}/process_{plate}.ipynb",
    params:
        # pass DataFrames/Series as dict/list for snakemake params rerun triggers
        samples=lambda wc: plates[wc.plate]["samples"]["sample"].tolist(),
        plate_params=lambda wc: {
            param: (val if param != "samples" else val.to_dict())
            for (param, val) in plates[wc.plate].items()
        },
    conda:
        "environment.yml"
    notebook:
        "notebooks/process_plate.py.ipynb"


checkpoint sera_by_plate:
    """Get list of all sera and plates they are on."""
    input:
        csvs=expand(rules.process_plate.output.fits_csv, plate=plates),
    output:
        csv="results/sera/sera_by_plate.csv",
    params:
        plates=list(plates),
    log:
        "results/logs/sera_by_plate.txt",
    conda:
        "environment.yml"
    script:
        "scripts/sera_by_plate.py"


rule serum_titers:
    """Aggregate and analyze titers for a serum."""
    input:
        pickles=lambda wc: [
            rules.process_plate.output.fits_pickle.format(plate=plate)
            for plate in sera_plates()[wc.serum]
        ],
    output:
        per_rep_titers="results/sera/{serum}/titers_per_replicate.csv",
        titers="results/sera/{serum}/titers.csv",
        curves_pdf="results/sera/{serum}/curves.pdf",
        pickle="results/sera/{serum}/curvefits.pickle",
        qc_drops="results/sera/{serum}/qc_drops.yml",
    params:
        viral_strain_plot_order=viral_strain_plot_order,
        serum_titer_as=lambda wc: (
            config["sera_override_defaults"][wc.serum]["titer_as"]
            if (
                (wc.serum in config["sera_override_defaults"])
                and ("titer_as" in config["sera_override_defaults"][wc.serum])
            )
            else config["default_serum_titer_as"]
        ),
        qc_thresholds=lambda wc: (
            config["sera_override_defaults"][wc.serum]["qc_thresholds"]
            if (
                (wc.serum in config["sera_override_defaults"])
                and ("qc_thresholds" in config["sera_override_defaults"][wc.serum])
            )
            else config["default_serum_qc_thresholds"]
        ),
    log:
        notebook="results/sera/{serum}/{serum}_titers.ipynb",
    conda:
        "environment.yml"
    notebook:
        "notebooks/serum_titers.py.ipynb"


rule aggregate_titers:
    """Aggregate all serum titers."""
    input:
        pickles=lambda wc: expand(rules.serum_titers.output.pickle, serum=sera_plates()),
        titers=lambda wc: expand(rules.serum_titers.output.titers, serum=sera_plates()),
    output:
        pickle="results/aggregated_titers/curvefits.pickle",
        titers="results/aggregated_titers/titers.csv",
        titers_chart="results/aggregated_titers/titers.html",
    params:
        viral_strain_plot_order=viral_strain_plot_order,
    conda:
        "environment.yml"
    log:
        notebook="results/aggregated_titers/aggregate_titers.ipynb",
    notebook:
        "notebooks/aggregate_titers.py.ipynb"


rule aggregate_qc_drops:
    """Aggregate all QC drops."""
    input:
        plate_qc_drops=expand(rules.process_plate.output.qc_drops, plate=plates),
        sera_qc_drops=lambda wc: expand(
            rules.serum_titers.output.qc_drops,
            serum=sera_plates(),
        ),
    output:
        plate_qc_drops="results/qc_drops/plate_qc_drops.yml",
        sera_qc_drops="results/qc_drops/sera_qc_drops.yml",
    params:
        plates=list(plates),
        sera=lambda wc: list(sera_plates()),
    conda:
        "environment.yml"
    log:
        notebook="results/qc_drops/aggregate_qc_drops.ipynb",
    notebook:
        "notebooks/aggregate_qc_drops.py.ipynb"


rule notebook_to_html:
    """Convert Jupyter notebook to HTML"""
    input:
        notebook="{notebook}.ipynb",
    output:
        html="{notebook}.html",
    log:
        "results/logs/notebook_to_html_{notebook}.txt",
    conda:
        "environment.yml"
    shell:
        "jupyter nbconvert --to html {input.notebook} &> {log}"


rule build_docs:
    """Build the HTML documentation."""
    input:
        titers_chart=rules.aggregate_titers.output.titers_chart,
        serum_titers_htmls=lambda wc: expand(
            "results/sera/{serum}/{serum}_titers.html",
            serum=sera_plates(),
        ),
        process_plates_htmls=expand(
            "results/plates/{plate}/process_{plate}.html",
            plate=plates,
        ),
        qc_drops_html="results/qc_drops/aggregate_qc_drops.html",
    output:
        docs=directory(config["docs"]),
    params:
        description=config["description"],
        sera=lambda wc: list(sera_plates()),
        plates=list(plates),
    conda:
        "environment.yml"
    log:
        "results/logs/build_docs.txt",
    script:
        "scripts/build_docs.py"


seqneut_pipeline_outputs = [
    rules.aggregate_titers.output.titers,
    rules.aggregate_titers.output.pickle,
    rules.aggregate_qc_drops.output.plate_qc_drops,
    rules.aggregate_qc_drops.output.sera_qc_drops,
    rules.build_docs.output.docs,
]

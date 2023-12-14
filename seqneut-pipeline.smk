"""Snakemake file for the ``seqneut-pipeline`.

Designed to be included in another ``Snakefile`` that specifies the config.

"""


import pandas as pd


snakemake.utils.min_version("7.32")


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


rule process_counts:
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
        qc_failures="results/plates/{plate}/process_counts_qc_failures.txt",
        frac_infectivity_csv="results/plates/{plate}/frac_infectivity.csv",
    log:
        notebook="results/plates/{plate}/process_counts_{plate}.ipynb",
    params:
        samples=lambda wc: plates[wc.plate]["samples"]["sample"],
        plate_params=lambda wc: {
            param: val
            for (param, val) in plates[wc.plate].items()
            if param not in {"curvefit_params"}
        },
    conda:
        "environment.yml"
    notebook:
        "notebooks/process_counts.py.ipynb"


rule qc_process_counts:
    """Check QC results on `process_counts` rule."""
    input:
        qc_failures=expand(rules.process_counts.output.qc_failures, plate=plates),
        process_counts_htmls=expand(
            "results/plates/{plate}/process_counts_{plate}.html",
            plate=plates,
        ),
    output:
        qc_summary="results/plates/qc_process_counts_summary.txt",
    conda:
        "environment.yml"
    params:
        plates=list(plates),
    log:
        "results/logs/qc_process_counts.txt",
    script:
        "scripts/qc_process_counts.py"


rule curvefits:
    """Fit neutralization curves for a plate."""
    input:
        qc_failures=rules.process_counts.output.qc_failures,
        frac_infectivity_csv=rules.process_counts.output.frac_infectivity_csv,
    output:
        csv="results/plates/{plate}/curvefits.csv",
        pdf="results/plates/{plate}/curvefits.pdf",
        pickle="results/plates/{plate}/curvefits.pickle",
    log:
        notebook="results/plates/{plate}/curvefits_{plate}.ipynb",
    params:
        curvefit_params=lambda wc: plates[wc.plate]["curvefit_params"],
    conda:
        "environment.yml"
    notebook:
        "notebooks/curvefits.py.ipynb"


checkpoint sera_by_plate:
    """Get list of all sera and plates they are on."""
    input:
        csvs=expand(rules.curvefits.output.csv, plate=plates),
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
        plate_fits=lambda wc: [
            rules.curvefits.output.csv.format(plate=plate)
            for plate in sera_plates()[wc.serum]
        ],
        pickles=lambda wc: [
            rules.curvefits.output.pickle.format(plate=plate)
            for plate in sera_plates()[wc.serum]
        ],
    output:
        per_rep_titers="results/sera/{serum}/titers_per_replicate.csv",
        median_titers="results/sera/{serum}/titers_median.csv",
        curves_pdf="results/sera/{serum}/curves.pdf",
        qc_failures="results/sera/{serum}/qc_failures.txt",
        pickle="results/sera/{serum}/curvefits.pickle",
    params:
        viral_strain_plot_order=viral_strain_plot_order,
        qc_thresholds=config["serum_titers_qc_thresholds"],
        qc_exclusions=lambda wc: (
            config["serum_titers_qc_exclusions"][wc.serum]
            if wc.serum in config["serum_titers_qc_exclusions"]
            else {}
        ),
    log:
        notebook="results/sera/{serum}/serum_titers_{serum}.ipynb",
    conda:
        "environment.yml"
    notebook:
        "notebooks/serum_titers.py.ipynb"


rule qc_serum_titers:
    """Check QC serum titeres from `serum_titers` rule."""
    input:
        qc_failures=lambda wc: expand(
            rules.serum_titers.output.qc_failures,
            serum=sera_plates(),
        ),
        serum_titers_htmls=lambda wc: expand(
            "results/sera/{serum}/serum_titers_{serum}.html",
            serum=sera_plates(),
        ),
    output:
        qc_summary="results/sera/qc_serum_titers_summary.txt",
    conda:
        "environment.yml"
    params:
        sera=lambda wc: list(sera_plates()),
    log:
        "results/logs/qc_serum_titers.txt",
    script:
        "scripts/qc_serum_titers.py"


rule aggregate_titers:
    """Aggregate all serum titers."""
    input:
        qc_serum_titer_failures=rules.qc_serum_titers.output.qc_summary,
        pickles=lambda wc: expand(rules.serum_titers.output.pickle, serum=sera_plates()),
        titers=lambda wc: expand(
            rules.serum_titers.output.median_titers, serum=sera_plates(),
        ),
    output:
        pickle="results/aggregated_titers/curvefits.pickle",
        titers="results/aggregated_titers/titers.csv",
        titers_chart="results/aggregated_titers/titers.html",
    params:
        viral_strain_plot_order=viral_strain_plot_order,
    conda:
        "environment.ymml"
    log:
        notebook="results/aggregated_titers/aggregate_titers.ipynb",
    notebook:
        "notebooks/aggregate_titers.py.ipynb"


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


seqneut_pipeline_outputs = [
    rules.qc_process_counts.output.qc_summary,
    rules.qc_serum_titers.output.qc_summary,
    rules.aggregate_titers.output.titers,
    rules.aggregate_titers.output.pickle,
]

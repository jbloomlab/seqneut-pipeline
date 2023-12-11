"""Snakemake file for the ``seqneut-pipeline`.

Designed to be included in another ``Snakefile`` that specifies the config.

"""


import copy
import os

import pandas as pd 


snakemake.utils.min_version("7.32")


include: "funcs.smk"  # include functions


# --- Process configuration ------------------------------------------------------------

pipeline_subdir = config["seqneut-pipeline"]

viral_libraries = config["viral_libraries"]

viral_strain_plot_order = get_viral_strain_plot_order(viral_libraries, config)

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
            neut_standard_sets[plates[samples[wc.sample]["plate"]]["neut_standard_set"]]
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
        "results/logs/count_barcodes_{sample}.txt"   
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
    log:
        notebook="results/plates/{plate}/curvefits_{plate}.ipynb",
    params:
        curvefit_params=lambda wc: plates[wc.plate]["curvefit_params"],
    conda:
        "environment.yml"
    notebook:
        "notebooks/curvefits.py.ipynb"


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
    expand(rules.count_barcodes.output.counts, sample=samples),
    expand(rules.process_counts.output.frac_infectivity_csv, plate=plates),
    expand(rules.curvefits.output.csv, plate=plates),
    rules.qc_process_counts.output.qc_summary,
]

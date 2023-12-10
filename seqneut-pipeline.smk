"""Snakemake file for the ``seqneut-pipeline`.

Designed to be included in another ``Snakefile`` that specifies the config.

"""


import copy
import os

import pandas as pd 


snakemake.utils.min_version("7.32")


# --- Process configuration ------------------------------------------------------------

pipeline_subdir = config["seqneut-pipeline"]

viral_libraries = config["viral_libraries"]

# process viral_strain_plot_order
viral_strain_plot_order = {viral_library: None for viral_library in viral_libraries}
if "viral_strain_plot_order" in config:
    viral_strain_plot_order = viral_strain_plot_order | config["viral_strain_plot_order"]
    if set(viral_strain_plot_order) != set(viral_libraries):
        raise ValueError(
            f"{viral_strain_plot_order.keys()=} != {viral_libraries.keys()=}"
        )
for _viral_library, _csv in viral_strain_plot_order.items():
    _viral_library_strains = sorted(
        set(pd.read_csv(viral_libraries[_viral_library])["strain"])
    )
    if _csv:
        _viral_order = pd.read_csv(_csv)["strain"].tolist()
        if len(_viral_order) != len(set(_viral_order)):
            raise ValueError(f"duplicate strains in viral_strain_order CSV {_csv}")
        if set(_viral_order) != set(_viral_library_strains):
            raise ValueError(
                f"viral_strain_order does not have correct strains for {_viral_library}"
            )
        viral_strain_plot_order[_viral_library] = _viral_order
    else:
        viral_strain_plot_order[_viral_library] = _viral_library_strains


neut_standard_sets = config["neut_standard_sets"]

def process_plate(plate, plate_params):
    """Process a plot from the configuration."""

    # Process plate parameters
    req_plate_params = {
        "date",
        "viral_library",
        "neut_standard_set",
        "samples_csv",
        "process_counts_qc_thresholds",
        "barcodes_to_drop",
        "wells_to_drop",
        "curvefit_params",
    }
    if not req_plate_params.issubset(plate_params):
        raise ValueError(f"{plate=} {plate_params=} lacks {req_plate_params=}")
    if plate_params["viral_library"] not in viral_libraries:
        raise ValueError(
            f"{plate=} {plate_params['viral_library']=} not in {viral_libraries=}"
        )
    if plate_params["neut_standard_set"] not in neut_standard_sets:
        raise ValueError(
            f"{plate=} {plate_params['neut_standard_set']=} not in {neut_standard_sets=}"
        )
    plate_d = copy.deepcopy(plate_params)
    plate_d["date"] = str(plate_d["date"])
    if not re.fullmatch("\d{4}\-\d{2}\-\d{2}", str(plate_d["date"])):
        raise ValueError(f"{plate=} {plate_d['date']=} not in YYYY-MM-DD format")

    # Process samples_csv to create the sample data frame
    req_sample_cols = ["well", "serum", "dilution_factor", "replicate", "fastq"]
    samples_df = pd.read_csv(plate_params["samples_csv"], comment="#")
    if not set(req_sample_cols).issubset(samples_df.columns):
        raise ValueError(f"{plate=} {samples_df.columns=} lacks {req_sample_cols=}")

    if samples_df["serum"].isnull().any():
        raise ValueError(f"{plate=} 'samples_csv' has null values in 'serum' column")

    # try to turn columns of ints and NAs into Int64 to avoid ints appearing as flaots
    for col in ["replicate", "dilution_factor"]:
        try:
            samples_df[col] = samples_df[col].astype("Int64")
        except TypeError:
            pass     

    # make serum_replicate that defines serum and replicate if needed
    samples_df = samples_df[req_sample_cols].assign(
        one_serum_replicate=lambda x: (
            x.groupby("serum")["replicate"].transform("nunique", dropna=False) == 1
        ),
        serum_replicate=lambda x: x.apply(
            lambda row: (
                str(row["serum"]) + (
                    "" if row["one_serum_replicate"] == 1 else f"-{row['replicate']}"
                )
            ),
            axis=1,
        ),
        sample_noplate=lambda x: x.apply(
            lambda row: (
                row["serum_replicate"] + (
                    ""
                    if pd.isnull(row["dilution_factor"])
                    else f"_{row['dilution_factor']}"
                )
            ),
            axis=1,
        ),
        sample=lambda x: plate + "_" + x["sample_noplate"],
        plate=plate,
        plate_replicate=lambda x: x.apply(
            lambda row: (
                plate + ("" if row["one_serum_replicate"] else f"{-row['replicate']}")
            ),
            axis=1,
        ),
    ).drop(columns="one_serum_replicate")

    assert len(samples_df) == samples_df["sample"].nunique(), plate

    # make sure serum_replicate and dilution_factor are unique
    dup_rows = (
        samples_df.assign(
            duplicates=lambda x: (
                x.groupby(["serum_replicate", "dilution_factor"], dropna=False)
                ["sample"]
                .transform("count")
            ),
        )
        .query("duplicates > 1")
        .drop(columns="duplicates")
    )
    if len(dup_rows):
        raise ValueError(f"{plate=} has duplicated serum / replicates:\n{dup_rows}")

    # make sure dilution_factor is valid
    if not ((samples_df["dilution_factor"] >= 1) | (samples_df["serum"] == "none")).all():
        raise ValueError(f"{plate=} has dilution factors not >= 1 for non-none serum")

    # make sure there is at least one "none" sample
    if "none" not in set(samples_df["serum"]):
        raise ValueError(f"{plate=} has no samples with serum set to 'none'")

    # make sure fastqs are unique
    dup_fastqs = (
        samples_df
        .assign(duplicates=lambda x: x.groupby("fastq")["fastq"].transform("count"))
        .query("duplicates > 1")
        .drop(columns="duplicates")
    )
    if len(dup_fastqs):
        raise ValueError(f"{plate=} has duplicate FASTQs:\n{dup_fastqs}")

    plate_d["samples"] = samples_df

    return plate_d


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

"""Snakemake file for the ``seqneut-pipeline`.

Designed to be included in another ``Snakefile`` that specifies the config.

"""


import copy
import os

import pandas as pd 


# --- Process configuration ------------------------------------------------------------

pipeline_subdir = config["seqneut-pipeline"]

viral_libraries = config["viral_libraries"]

neut_standard_sets = config["neut_standard_sets"]

def process_plate(plate, plate_params):
    """Process a plot from the configuration."""

    # Process plate parameters
    req_plate_params = {"date", "viral_library", "neut_standard_set", "samples_csv"}
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
    req_sample_cols = {"serum", "dilution_factor", "replicate", "fastq"}
    samples_df = pd.read_csv(plate_params["samples_csv"])
    if not req_sample_cols.issubset(samples_df.columns):
        raise alueError(f"{plate=} {samples_df.columns=} lacks {req_sample_cols=}")

    if samples_df["serum"].isnull().any():
        raise ValueError(f"{plate=} 'samples_csv' has null values in 'serum' column")

    # try to turn columns of ints and NAs into Int64 to avoid ints appearing as flaots
    for col in ["replicate", "dilution_factor"]:
        try:
            samples_df[col] = samples_df[col].astype("Int64")
        except TypeError:
            pass     

    # make serum_replicate that defines serum and replicate if needed
    samples_df = samples_df.assign(
        n_serum_replicates=lambda x: (
            x.groupby("serum")["replicate"].transform("nunique", dropna=False)
        ),
        serum_replicate=lambda x: x.apply(
            lambda row: (
                str(row["serum"]) + (
                    "" if row["n_serum_replicates"] == 1 else f"-{row['replicate']}"
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
    ).drop(columns="n_serum_replicates")

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
        "environment.yml"
    log:
        "results/logs/count_barcodes_{sample}.txt"   
    script:
        "scripts/count_barcodes.py"


rule process_plate:
    """Process a plate."""
    input:
        count_csvs=lambda wc: expand(
            rules.count_barcodes.output.counts,
            sample=plates[wc.plate]["samples"]["sample"],
        ),
        fate_csvs=lambda wc: expand(
            rules.count_barcodes.output.fates,
            sample=plates[wc.plate]["samples"]["sample"],
        ),
        invalid_count_csvs=lambda wc: expand(
            rules.count_barcodes.output.invalid,
            sample=plates[wc.plate]["samples"]["sample"],
        ),
        viral_library_csv=lambda wc: (
            viral_libraries[plates[wc.plate]["viral_library"]]
        ),
        neut_standard_set_csv=lambda wc: (
            neut_standard_sets[plates[wc.plate]["neut_standard_set"]]
        ),
    output:
        frac_infectivity_csv="results/plates/{plate}/frac_infectivity.csv",
    log:
        notebook="results/plates/{plate}/process_{plate}.ipynb",
    params:
        samples=lambda wc: plates[wc.plate]["samples"]["sample"],
        plate_params=lambda wc: plates[wc.plate],
    conda:
        "environment.yml"
    notebook:
        "notebooks/process_plate.py.ipynb"


rule calculate_neutralization_potency:
    """Process fraction infectivity files to calculate NT50s and generate plots."""
    input:
#        fractioninfectivity=expand(rules.calculate_fraction_infectivity.output.fraction_infectivity, plate=plates),
        ipynb=os.path.join(pipeline_subdir, "notebooks/calculate-neutralization-potency.ipynb"),
    output:
        ipynb="results/notebooks/calculate-neutralization-potency.ipynb",
        html="results/notebooks/calculate-neutralization-potency.html",
        median_ic50s="results/selections/nt50_measurements_by_strain.csv",
    conda:
        "envs/calculate-neutralization-potency.yml"
    log:
        "results/logs/calculate_neutralization_potency.txt",
    shell:
        """
        papermill {input.ipynb} {output.ipynb} \
            -p snakemake True \
            -p median_ic50s {output.median_ic50s} \
            &> {log}

        jupyter nbconvert --to html {output.ipynb}
        """


#### ----------------------- Imports ----------------------- ####

import pandas as pd 
from os.path import join

#### ------------------------ Input ------------------------ ####

if "pipeline_subdir" in config:
    pipeline_subdir = config["pipeline_subdir"]
else:
    pipeline_subdir = "seqneut_pipeline"

# Read in the barcode runs 
barcode_runs = pd.read_csv(config["barcode_runs"])

# Check for the necessary columns
required_columns = {
    "date",
    "plate",
    "library",
    "standard_set",
    "fastq",
}
required_columns.update(config["id_columns"])
missing_columns = required_columns - set(barcode_runs.columns)
if missing_columns:
    raise ValueError(
        f"The following columns must exist in the `barcode_runs` dataframe: {list(missing_columns)}"
    )

# Make the 'sample' column using the id columns
if "sample" in barcode_runs.columns:
    raise ValueError(
        f"The `barcode_runs` dataframe already has a column called 'sample', please rename the column or remove it."
    )
barcode_runs["sample"] = barcode_runs[config["id_columns"]].apply(
    lambda x: "-".join(x.astype(str)), axis=1
)

# Check that the sample names are unique
if len(barcode_runs["sample"].unique()) != len(barcode_runs):
    raise ValueError(f"The sample names derived from the provided id columns are not unique.")


# Check that the neutralization standard file has the 'barcode' and 'standard_set' columns
standards_required_columns = {"barcode", "standard_set"}
standards_missing_columns = standards_required_columns - set(pd.read_csv(config["neut_standards"]).columns)
if standards_missing_columns:
    raise ValueError(
        f"The following columns must exist in the `neut_standards` dataframe: {list(standards_missing_columns)}"
    )

# Check that the strain to barcode file has the 'strain', 'barcode' and 'library' columns
library_required_columns = {"strain", "barcode", "library"}
library_missing_columns = library_required_columns - set(pd.read_csv(config["strain_to_barcode"]).columns)
if library_missing_columns:
    raise ValueError(
        f"The following columns must exist in the `strain_to_barcode` dataframe: {list(library_missing_columns)}"
    )


# Get a list of the sample names
samples = barcode_runs["sample"].unique().tolist()
# Get a list of plate names
plates = barcode_runs["plate"].unique().tolist()


#### ------------------------ Rules ------------------------ ####

rule clean:
    shell:
        """
        rm -rf logs/
        rm -rf tmp/
        rm -f slurm*.out
        """

rule count_barcodes:
    """Count barcodes for each sample."""
    input:
        fastq=lambda wildcards: barcode_runs.set_index("sample").at[wildcards.sample, "fastq"],
        variants=config["strain_to_barcode"],
        standards=config["neut_standards"],
    output:
        counts="results/barcode_counts/{sample}/{sample}_counts.csv",
        invalid="results/barcode_counts/{sample}/{sample}_invalid.csv",
        fates="results/barcode_counts/{sample}/{sample}_fates.csv",
    params:
        library=lambda wildcards: barcode_runs.set_index("sample").at[wildcards.sample, "library"],
        standard_set=lambda wildcards: barcode_runs.set_index("sample").at[wildcards.sample, "standard_set"],
        parser_params=config["illumina_barcode_parser_params"]
    conda:
        "envs/count_barcodes.yml"
    log:
        "results/logs/count_barcodes_{sample}.txt"   
    script:
        "scripts/count_barcodes.py"


rule analyze_barcode_counts:
    """Process barcode counts and perform basic quality control."""
    input:
        counts=expand(rules.count_barcodes.output.counts, sample=samples),
        invalid=expand(rules.count_barcodes.output.invalid, sample=samples),
        fates=expand(rules.count_barcodes.output.fates, sample=samples),
        ipynb=os.path.join(pipeline_subdir, "notebooks/analyze-barcode-counts.ipynb"),
    output:
        ipynb="results/notebooks/analyze-barcode-counts.ipynb",
        html="results/notebooks/analyze-barcode-counts.html",
        joined_counts="results/barcode_counts/barcode_counts.csv",
        counts_by_plate=expand("results/barcode_counts/{plate}_barcode_counts.csv", plate=plates),
    conda:
        "envs/count_barcodes.yml"
    log:
        "results/logs/analyze_barcode_counts.txt",
    shell:
        """
        papermill {input.ipynb} {output.ipynb} \
            -p joined_counts {output.joined_counts} \
            -p snakemake True \
            &> {log}

        jupyter nbconvert --to html {output.ipynb}
        """

rule calculate_fraction_infectivity:
    """Process counts files by plate and calculate fraction infectivity."""
    input:
        variants=config["strain_to_barcode"],
        standards=config["neut_standards"],
        counts="results/barcode_counts/{plate}_barcode_counts.csv",
    output:
        fraction_infectivity="results/fraction_infectivity/{plate}_fractioninfectivity.csv",
    log:
        "results/logs/calculate_fraction_infectivity_{plate}.txt",
    script:
        "scripts/calculate_fraction_infectivity-perplate.py"


rule calculate_neutralization_potency:
    """Process fraction infectivity files to calculate NT50s and generate plots."""
    input:
        fractioninfectivity=expand(rules.calculate_fraction_infectivity.output.fraction_infectivity, plate=plates),
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


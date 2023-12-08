# `seqneut-pipeline` for analyzing sequencing-based neutralization assays

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/charliermarsh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![Snakemake](https://img.shields.io/badge/snakemake-≥7.32-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)

---

This is a modular analysis pipeline for analyzing high-throughput sequencing-based neutralization assays of the type developed in the [Bloom lab](https://research.fredhutch.org/bloom/en.html).
See **[add link to Loes et al when available]** for a description of these assays.
That paper is also the scientific citation for this analysis pipeline.

## Using this pipeline
This pipeline is designed to be included as a modular portion of a larger [snakemake](https://snakemake.readthedocs.io/) analysis.
This pipeline processes FASTQ files to get the counts of each barcode, analyzes those to determine the fraction infectivity at each serum concentration, and then fits and plots neutralization curves.

To use the pipeline, first create another repository specific to your project.
The include this pipeline as a [git submodule](https://git-scm.com/book/en/v2/Git-Tools-Submodules) in your repository, where it will be present in a subdirectory named `seqneut-pipeline`.
So the overall structure will look like this:

```
<your_project_repo>
├── seqneut-pipeline [added as git submodule]
│   ├── seqneut-pipeline.smk [snakemake rules for pipeline]
│   ├── environment.yml [conda environment for pipeline]
│   └── ... [other contents of seqneut-pipeline]
├── README.md [README for your project]
├── config.yml [YAML configuration for project]
├── Snakefile [top-level snakemake file that includes `seqneut-pipeline.smk`]
├── data [subdirectory with input data for your project]
├── results [subdirectory with results created by pipeline]
├── docs [HTML summary of results created by pipeline]
└── <other files / subdirectories that are part of project>
```

So after you have created your project repo, add this pipeline as a git submodule with:

      git submodule add https://github.com/jbloomlab/seqneut-pipeline

This creates a file called `gitmodules` and the `seqneut-pipeline` submodule, which can then be committed to the repo.
If at some point you want to updated the version of the pipeline, simply `cd` into the `seqneut-pipeline` subdirectory and pull or checkout the version you want.

To use the pipeline, you then need to add a few things to your main-level repo.
The first is a top-level `Snakefile` that includes `seqneut-pipeline`, reads in your configuration, file, and contains as some of its targets the outputs of `seqneut-pipeline`.
So minimally that top-level `Snakefile` should contain the following lines (it can also contain additional stuff if you also have it running project-specific analyses on the output of the pipeline):
```
configfile: "config.yml"

include: "seqneut-pipeline/seqneut-pipeline.smk"

rule all:
    input:
        seqneut_pipeline_outputs
```

In addition, you need to create the configuration file `config.yml` and ensure it includes the appropriate configuration for `seqneut-pipeline` as described below.

Finally, you need to create a `conda` environment that minimally includes the packages needed to run the pipeline, which are a recent version of [snakemake](https://snakemake.readthedocs.io/) and [pandas](https://pandas.pydata.org/).
You can either create your own environment containing these, or simply build and use the one specified in [environment.yml](environment.yml) file of `seqneut-pipeline`, which is named `seqneut-pipeline`. So if you are using that environment, you can simply run the pipeline with:
```
conda activate seqneut-pipeline
snakemake -j <n_jobs> --use-conda --keep-going
```

The use of `--keep-going` is recommended for the QC steps below as it will create the notebooks helpful for manually doing the QC.

## Configuring the pipeline
The configuration for the pipeline is in a file called `config.yml`.
An example configuration file is in [./test_example/config.yml](test_example/config.yml).
Here we describe the required keys in this YAML file (you can also add additional information specific to your repo, but we advise adding comments to put that in a separate part of the YAML from the `seqneut-pipeline` configuration).
For background on YAML format, including on the anchor (`&`) and merge (`<<: *`) operators that can be helpful to simplify the YAML file, see [here](https://spacelift.io/blog/yaml) and [here](https://ktomk.github.io/writing/yaml-anchor-alias-and-merge-key.html).

The top-level keys in the YAML are:

### seqneut-pipeline
Location of the `seqneut-pipeline` relative to the top-level repo.
This will almost always be a subdirectory of the same name, so this key will be as shown below unless you have a good reason to do otherwise:

        seqneut-pipeline: seqneut-pipeline

### viral_libraries
A dictionary (mapping) of viral library names to CSV files holding these libraries.
So in general this key will look like:
```
viral_libraries:
  pdmH1N1_lib2023_loes: data/viral_libraries/pdmH1N1_lib2023_loes.csv
  <potentially more viral libraries specified as name: CSV pairs>
```
Note how the recommended way to organize the viral libraries (as indicated above) is to put them in a `./data/viral_libraries/` subdirectory.

The CSV files themselves will have columns specifying the viral barcode and the strain it corresponds to, such as:
```
barcode,strain
ACGGAATCCCCTGAGA,A/Washington/23/2020
GCATGGATCCTTTACT,A/Togo/845/2020
<additional lines>
```

### neut_standard_sets
A dictionary (mapping) of neutralization-standard set names to CSV files holding the barcodes for the neutralization standard set.
So in general, this key will look like:
```
neut_standard_sets:
  loes2023: data/neut_standard_sets/loes2023_neut_standards.csv
  <potentially more neut standard sets specified as name: CSV pairs>
```
Note how the recommended way to organize the viral libraries (as indicated above) is to put them in a `./data/viral_libraries/` subdirectory.

The CSV files need just a single column specifying the neutralization standard barcode, such as:
```
barcode
CTTTAAATTATAGTCT
CATACAGAGTTTGTTG
<additional lines>
```

### illumina_barcode_parser_params
A dictionary (mapping) specifying how to parse the Illumina FASTQ files to barcode counts.
This mapping should just specify key-word arguments that can be passed to [dms_variants.illuminabarcodeparser.IlluminaBarcodeParser](https://jbloomlab.github.io/dms_variants/dms_variants.illuminabarcodeparser.html#dms_variants.illuminabarcodeparser.IlluminaBarcodeParser); note that the barcode-length (`bclen`) parameter should not be specified as it is inferred from the length of the barcodes.

So in general, this key will look like this:
```
illumina_barcode_parser_params:
  upstream: CTCCCTACAATGTCGGATTTGTATTTAATAG
  downstream: ''
  minq: 20
  upstream_mismatch: 4
  bc_orientation: R2
```

#### plates
This dictionary (mapping) contains the heart of the configuration, and may be quite large.
Essentially, it specifies what samples are contained in each plate, how those samples should be processed, QC thresholds, and any specific barcodes or samples that should be excluded for failing QC.

The basic structure is that `plates` specifies additional mappings of plate names to configuration for that plate.
Specifically, this key should look like this:
```
plates:

  plate1:
    date: 2023-08-01
    viral_library: pdmH1N1_lib2023_loes
    neut_standard_set: loes2023
    samples_csv: data/plates/plate1_samples.csv
    process_counts_qc_thresholds:
      <<: *default_process_counts_qc_thresholds
    barcodes_to_drop:
      - AAAACAGTATAGAAGA  # viral barcode w no counts in no-serum samples
      - AATCTCCTCACGCAGC  # viral barcode w no counts in no-serum samples
      - ATGCAATATTAAGGAA  # viral barcode w no counts in no-serum samples
      - GGTCCATCTCAGATCG  # neut standard barcode w inconsistently low counts in Y106d0_4860
    wells_to_drop: []

  <additional_plates>
```
The above example shows the configuration of a plate called `plate1`, and there may be many additional plates.
The elements under each plate-mapping are in turn as follows:

#### date
The `date` key specifies the date on which the plate was processed in `YYYY-MM-DD` format.

#### viral_library
The `viral_library` key gives the name of a key in `viral_libraries` that specifies the barcodes / strains for the viral library used for this plate.

#### neut_standard_set
The `neut_standard_set` key gives the name of a key in `neut_standard_sets` that specifies the neutralization-standard set barcodes used for this plate.

#### samples_csv
The `samples_csv` key gives the name of a CSV file specifying the samples for that plate.
The recommended way to organize these sample CSVs is to put them in `./data/plates/` subdirectory.
The CSV file must have the following columns:
 - *well*: well in plate in which sample was run, typically names like "A1", "B1", etc.
 - *serum*: name of the serum in this well, or "none" if it is a no-serum sample.
 - *dilution_factor*: dilution factor of the serum (should be a number > 1), leave blank for the no-serum samples (*serum* of "none")
 - *replicate*: the replicate of this serum, which you only need to specify if there are multiple different samples with the same *serum* and *dilution_factor* in the plate.
 - *fastq*: path to the FASTQ file, can be gzipped
 - other columns (e.g., *notes*, etc) are allowed but are ignored by the pipeline

Here are a few lines of an example CSV file:
```
well,serum,dilution_factor,replicate,fastq
A1,none,,1,/fh/fast/bloom_j/SR/ngs/illumina/aloes/230801_VH00319_391_AACWKHTM5/Unaligned/Project_aloes/Plate1_Noserum1_S1_R1_001.fastq.gz
A2,Y106d182,20.0,1,/fh/fast/bloom_j/SR/ngs/illumina/aloes/230801_VH00319_391_AACWKHTM5/Unaligned/Project_aloes/Y106_d182_conc1_S9_R1_001.fastq.gz
<additional lines>
```

#### process_counts_qc_thresholds
This key defines a mapping of the quality-control thresholds for processing the sequencing counts to get fraction infectivities.
These thresholds are used for the QC in the `process_count` rule (see section below on quality-control for more details).

Since it is a bit complex, you may want to use the [YAML anchor / merge](https://ktomk.github.io/writing/yaml-anchor-alias-and-merge-key.html) syntax to define a default that you then merge for specific plates.
The default can be defined like this:
```
default_process_counts_qc_thresholds: &default_process_counts_qc_thresholds
  avg_barcode_counts: 500
  min_neut_standard_frac: 0.01
  max_neut_standard_frac_no_serum: 0.1
  barcode_frac_consistency: 3
  min_viral_barcode_frac: 0.001
  min_neut_standard_barcode_frac: 0.01
  min_neut_standard_count: 1000
  min_no_serum_viral_barcode_count: 10
  min_dilutions_per_serum_replicate: 4
  max_frac_infectivity: 5
```
and then for specific plates you can merge it in and overwrite.
For instance, below would merge the above defaults but then overwrite the `min_viral_barcode_frac` to a different value:
```
plates:

  plate1:
    <other keys>
    process_counts_qc_thresholds:
      <<: *default_process_counts_qc_thresholds  # merge in defaults
      min_viral_barcode_frac: 0.0002  # overwrite default for this key for this plate

  <other plates>
```

The keys for `counts_qc_thresholds` define the following thresholds:

##### avg_barcode_counts
Require at least this many average counts per barcode (averaged across barcodes) for each sample.
Designed to ensure enough sequencing for barcodes for each sample.

##### min_neut_standard_frac
Require at least this fraction of counts to come from the neutralization standard for each sample.
Designed to ensure all samples have adequate use of neutralization standard.

##### max_neut_standard_frac_no_serum
Require that no more than this fraction of counts comes from the neutralization standard in each no-serum sample.
Designed to check against excessive neutralization standard use.

##### barcode_consistency_frac
Require that the frequency (fraction of counts) for each barcode has a fold-change from the median across samples that does not exceed this amount in the samples where it should be consistent.
Specifically, in the no-serum samples all viral barcodes should be at consistent relative levels across samples, and this is checked using this threshold.
Likewise, in all samples the neutralization standard barcodes should be at consistent levels, and this is checked using this threshold.

##### min_viral_barcode_frac
Require that each viral barcode comprises at least this fraction of total counts in the no-serum sample.
Designed to make sure no viral barcodes are too rare.

##### min_neut_standard_barcode_frac
Require that each neut-standard barcode is at least this fraction of all counts in all samples.
Designed to make sure no neut-standard barcode is too rare.

##### min_neut_standard_counts
Require that each sample has at least this number of neut-standard barcode counts (summed across all neut-standard barcodes).
Designed to make sure all samples have enough of these barcodes to accurately estimate the fraction infectivity.

##### min_no_serum_viral_barcode_count
Require each viral barcode has at least this number of counts in each no-serum sample.
Designed to make sure each viral barcode adequately represented to get good estimates.

##### min_dilutions_per_serum_replicate
Require at least this many different dilutions for each serum-replicate on a plate.
Designed to make sure we have enough dilutions to actually fit a curve.

##### max_frac_infectivity
Require all fraction infectivity values to be <= this number.
They should really all be <= one, but sometimes they will be larger due to noise but it's a bad sign if they get too large.

#### barcodes_to_drop
If there are barcodes that are failing some of the QC above, you can specify them here and they will be excluded.
As you add barcodes to drop for a plate, you should add a comment in the YAML on why the barcode is being dropped.
Note that if a barcode is entirely missing from the viral library for many plates, you may want to update the `viral_barcodes` entry instead to just remove it.

#### wells_to_drop
If there are samples that are failing some of the QC above, you can specify their well identifiers in a list here and they will be dropped.
As you add wells (samples) to drop for a plate, you should add a comment in the YAML on why the sample is being dropped.

## Output of the pipeline
The results of running the pipeline are put in the `./results/` subdirectory of your main repo.
We recommend using the `.gitignore` file in [./test_example/.gitignore] in your main repo to only track key results in your GitHub repo.
The set of full created outputs are as follows (note only some will be tracked depending on your `.gitignore`):

  - `./results/barcode_counts/`: files giving the barcode counts for each sample. You should track this in the repo.
  - `./results/barcode_fates/`: files giving the statistics (fates) of reads in the barcode counting for each sample. You do not need to track this in the repo as the results are plotted.
  - `./results/barcode_invalid/`: files giving counts of invalid barcodes for each sample. You do not need to track this in the repo, but it could be helpful to look at these identities in counts if QC shows you are getting many invalid barcodes.
  - `./results/plates/{plate}/frac_infectivity.csv`: fraction infectivity for viral barcodes for a plate. You should track this in the repo.
  - `./results/plates/{plate}/process_counts_{plate}.ipynb`: Jupyter notebook processing counts for a plate. You do not need to track this, look at the HTMl version of notebook instead.
  - `./results/plates/{plate}/process_counts_{plate}.html`: HTML of Jupyter notebook processing counts for a plate. You do not need to track this as it will be rendered in docs when pipeline runs successfully.
  - `./results/plates/{plate}/process_counts_qc_failures.txt`: List of QC failures when processing counts for plate. You do not need to track this as the summary for all plates is tracked instead.
  - `./results/plates/qc_process_counts_summary.txt`: summary of QC for processing counts for all plates. You should track this in the repo.
  - `./logs/`: logs from `snakemake` rules, you may want to look at these if there are rule failures.

## Running pipeline to identify QC failures and fixing them
If you run the pipeline via `snakemake` with the `--keep-going` flag as recommended above, the pipeline will run as far as possible.
However, if there are any QC failures that will keep it from running to completion.

For the processing of counts to fraction infectivity, the file `./results/plates/qc_process_counts_summary.txt` will summarize the QC failures and tell you which HTML notebooks to look at for details.
You then need to address these QC failures by doing one of the following:

 - Removing the offending barcodes by adding them to `barcodes_to_drop` for that plate. (If the barcode is missing in all plates, you might remove from `viral_barcodes` or `neut_standard_sets`.)

 - Remove the offending samples by adding them to `wells_to_drop` for that plate.

 - Adjusting the `process_counts_qc_thresholds` for that plate to be more lenient.

It is expected that you may have to perform several iterations of running and fixing QC failures.

# `seqneut-pipeline` for analyzing sequencing-based neutralization assays

[![Release](https://img.shields.io/github/v/release/jbloomlab/seqneut-pipeline?logo=github)](https://github.com/jbloomlab/seqneut-pipeline/releases)
[![Build Status](https://github.com/jbloomlab/seqneut-pipeline/actions/workflows/test.yaml/badge.svg)](https://github.com/jbloomlab/seqneut-pipeline/actions/workflows/test.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/charliermarsh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)

---

This is a modular analysis pipeline for analyzing high-throughput sequencing-based neutralization assays of the type developed in the [Bloom lab](https://jbloomlab.org).
See [Loes et al (2024)](https://doi.org/10.1128/jvi.00689-24) for a description of these assays.


Please cite [Loes et al (2024)](https://doi.org/10.1128/jvi.00689-24) if you use this pipeline for your scientific study.

See [here](https://github.com/jbloomlab/seqneut-pipeline/graphs/contributors) for a list of the contributors to this pipeline.

## Overview

This pipeline goes from the FASTQ files that represent the counts of each barcoded viral variant to the computed neutralization titers for each sera.
The titers are computed by fitting Hill-curve style neutralization curves using the [neutcurve](https://jbloomlab.github.io/neutcurve/) package; see the documentation for the details of these curves.
The titers represent the reciprocal serum dilutions at which half the viral infectivity is neutralized.
The pipeline provides options to compute these titers as either:
 - the reciprocal of the inhibitory concentration 50% (IC50), namely as the neutralization titer 50% (NT50)
 - the reciprocal of the midpoint of the neutralization curve

When the curves are fit with a top plateau of 1 and a bottom plateau of zero, these two different ways of calculating the titers are identical, and represent the serum dilution factor at which the serum neutralizes half of the viral infectivity.
Typically there will be multiple replicates, and the final reported titer is the median among replicates.

The pipeline also performs extensive quality control at different steps using configurable options described below.

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
If at some point you want to update the version of the pipeline, simply `cd` into the `seqneut-pipeline` subdirectory and pull or checkout the version you want.

To use the pipeline, you then need to add a few things to your main-level repo.
The first is a top-level `Snakefile` that includes your configuration, `seqneut-pipeline`, and outputs of `seqneut-pipeline` as targets of the `all` rule.
So minimally that top-level `Snakefile` should contain the following lines (it can also contain additional stuff if you also have it running project-specific analyses on the output of the pipeline):
```
import os

configfile: "config.yml"

include: os.path.join(config["seqneut-pipeline"], "seqneut-pipeline.smk")

rule all:
    input:
        seqneut_pipeline_outputs
```

In this `Snakefile`, the `seqneut_pipeline_outputs` specify files created by the pipeline.
Several of of these may be of special interest for you to use in additional rules you define in `Snakefile`:
  - `./results/aggregated_titers/titers_{group}.csv`: CSV with the final (median across replicates) titers for each serum-virus pair in a group after applying quality control filters.
  - `./results/aggregated_titers/curvefits_{group}.pickle`: a pickled [neutcurve.CurveFits](https://jbloomlab.github.io/neutcurve/neutcurve.curvefits.html#neutcurve.curvefits.CurveFits) object with all of the curve fits for all serum-virus-replicates in a group after applying the QC filters. You can use the methods of this object to make plots of neutralization curves for specific sera / viruses / replicates.

In addition, you need to create the configuration file `config.yml` and ensure it includes the appropriate configuration for `seqneut-pipeline` as described below.

To track the correct files in the created results, we suggest you copy the [./test_example/.gitignore](test_example/.gitignore) file to be the `.gitignore` for your main repo.
This will track key results files, but not an excessive number of non-essential files.

Finally, you need to create a `conda` environment that minimally includes the packages needed to run the pipeline, which are a recent version of [snakemake](https://snakemake.readthedocs.io/) and [pandas](https://pandas.pydata.org/).
You can either create your own environment containing these, or simply build and use the one specified in [environment.yml](environment.yml) file of `seqneut-pipeline`, which is named `seqneut-pipeline`. So if you are using that environment, you can simply run the pipeline with:
```
conda activate seqneut-pipeline
snakemake -j <n_jobs> --software-deployment-method conda
```

Note also that a few rules have rule-specific `conda` environments in [./envs/](envs).

## Configuring the pipeline
The configuration for the pipeline is in a file called `config.yml`.
An example configuration file is in [./test_example/config.yml](test_example/config.yml) (although some of the QC thresholds are set more leniently to make the test example work for small data as described in the comments in that YAML).

Here we describe the required keys in this YAML file (you can also add additional information specific to your repo, but we advise adding comments to put that in a separate part of the YAML from the `seqneut-pipeline` configuration).
For background on YAML format, including on the anchor (`&`) and merge (`<<: *`) operators that can be helpful to simplify the YAML file, see [here](https://spacelift.io/blog/yaml) and [here](https://ktomk.github.io/writing/yaml-anchor-alias-and-merge-key.html).

The top-level keys in the YAML are:

### seqneut-pipeline
Location of the `seqneut-pipeline` relative to the top-level repo.
This will almost always be a subdirectory of the same name, so this key will be as shown below unless you have a good reason to do otherwise:

        seqneut-pipeline: seqneut-pipeline

### docs
Location where we create the `./docs/` subdirectory with HTMLs for rendering on GitHub pages.
This will almost always be `docs`, so this key will be as shown below unless you have a good reason to do otherwise:

        docs: docs

### description
Description of pipeline, used in the HTML docs rendering.
Should include title (with markdown `#` heading, authors and/or citation, and link to GitHub repo.
For instance:
```
description: |
  # <title>
  <short description>

  <authors and/or link to citation>

  See <GitHub repo link> for code and numerical data.
```

### viral_libraries
A dictionary (mapping) of viral library names to CSV files holding these libraries.
So in general this key will look like:
```
viral_libraries:
  pdmH1N1_lib2023_loes: data/viral_libraries/pdmH1N1_lib2023_loes.csv
  <potentially more viral libraries specified as name: CSV pairs>
```
The recommended way to organize the viral libraries (as indicated above) is to put them in a `./data/viral_libraries/` subdirectory.

The CSV files themselves will have columns specifying the viral barcode and the strain it corresponds to, such as:
```
barcode,strain
ACGGAATCCCCTGAGA,A/Washington/23/2020
GCATGGATCCTTTACT,A/Togo/845/2020
<additional lines>
```

### viral_strain_plot_order
A a CSV with a column named "strain" that lists the strains in the order they should be plotted.
If not specified or set to "null", then plotting is just alphabetical.
Must include all strains being used if specified.
So should look like this:
```
viral_strain_plot_order: data/viral_strain_plot_order.csv
```

The CSV file itself will just have a column named "strain" specifying the order, such as:
```
strain
A/California/07/2009
A/Michigan/45/2015
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
The recommended way to organize the neutralization-standard sets (as indicated above) is to put them in a `./data/neut_standard_sets/` subdirectory.

The CSV files need just a single column specifying the neutralization standard barcode, such as:
```
barcode
CTTTAAATTATAGTCT
CATACAGAGTTTGTTG
<additional lines>
```

### illumina_barcode_parser_params
A dictionary (mapping) specifying how to parse the Illumina FASTQ files to barcode counts.
This is a global dictionary that is applied to all plates, but can be augmented or overriden on a per-plate basis by specifying plate specific `illumina_barcode_parser_params` as described in the plate configuration below.
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

### plates
This dictionary (mapping) contains the heart of the configuration, and may be quite large.
Essentially, it specifies what samples are contained in each plate, how those samples should be processed, QC thresholds, and any specific barcodes or samples that should be dropped.
In addition, each plate is assigned to a *group*, which might be "serum" or "pilot" (if you are mixing analyses of your sera with pilot experiments), or could be additional groups if you have two distinct sets of sera.

The basic structure is that `plates` maps plate names to configurations for the plates.
Specifically, it should look like this:
```
plates:

  plate1:
    group: serum
    date: 2023-08-01
    viral_library: pdmH1N1_lib2023_loes
    neut_standard_set: loes2023
    samples_csv: data/plates/plate1_samples.csv
    manual_drops: {}
    qc_thresholds:
      <<: *default_process_plate_qc_thresholds
    curvefit_params:
      <<: *default_process_plate_curvefit_params
    curvefit_qc:
      <<: *default_process_plate_curvefit_qc
    illumina_barcode_parser_params:  # optional argument
        upstream2: GCTACA

  <additional_plates>
```
The above example shows the configuration of a plate called `plate1`, and there may be many additional plates.
The elements under each plate-mapping are in turn as follows:

#### group
The group that this plate is assigned to (cannot contain any underscores). Typically this might be "serum" or "pilot" or however you are categorizing the runs.

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

#### manual_drops
As you analyze plates, you may find specific barcodes, wells, etc that you want to drop even if they don't fail the QC if they appear problematic to you for some reason.
If so, specify them using this key (if you don't want to manually drop any data for the plate, then just set this key to an empy dictionary `{}`).

The manual drops can have the following keys:

 - `wells`: list of wells to drop

 - `barcodes`: list of barcodes to drop from all wells

 - `barcode_wells`: list of `[barcode, well]` lists to drop specific barcodes in specific wells

 - `barcode_serum_replicates`: list of `[barcode, serum_replicate]` to drop specific barcodes for specific serum-replicates

 - `serum_replicates`: list of serum-replicates to drop

So for instance, you could have this specification if you wanted to drop barcode `AGTCCTATCCTCAAAT` for all wells of serum-replicate `M099d0`
```
manual_drops:
  barcode_serum_replicates:
    - [AGTCCTATCCTCAAAT, M099d0]
```

#### qc_thresholds
This key defines a mapping of the quality-control thresholds for processing the sequencing counts to get fraction infectivities.
These thresholds are used for the QC in the `process_count` rule (see section below on quality-control for more details).

Since it is a bit complex, you typically will want to use the [YAML anchor / merge](https://ktomk.github.io/writing/yaml-anchor-alias-and-merge-key.html) syntax to define a default that you then merge for specific plates.
The default can be defined like this:
```
default_process_plate_qc_thresholds: &default_process_plate_qc_thresholds
  avg_barcode_counts_per_well: 500
  min_neut_standard_frac_per_well: 0.005
  no_serum_per_viral_barcode_filters:
    min_frac: 0.0005
    max_fold_change: 3
    max_wells: 2
  per_neut_standard_barcode_filters:
    min_frac: 0.005
    max_fold_change: 3
    max_wells: 2
  min_neut_standard_count_per_well: 1000
  min_no_serum_count_per_viral_barcode_well: 30
  max_frac_infectivity_per_viral_barcode_well: 5
  min_dilutions_per_barcode_serum_replicate: 6
```
and then for specific plates you can merge this default it in and overwrite any specific keys if needed.
For instance, below would merge the above defaults but then overwrite the `min_viral_barcode_frac` to a different value:
```
plates:

  plate1:
    <other keys>
    qc_thresholds:
      <<: *default_process_plate_qc_thresholds  # merge in defaults
      avg_barcode_counts_per_well: 1000  # overwrite default for this key for this plate

  <other plates>
```

The QC-thresholds defined here are applied in order to drop data (wells, barcodes, etc) when processing the plates.
Specifically:

 - `avg_barcode_counts_per_well`: drop any well that does not have at least this many average counts per barcode.

 - `min_neut_standard_frac_per_well`: drop any well where the neutralization standard is not at least this fraction of the counts in the well.

 - `no_serum_per_viral_barcode_filters`: has subkeys `min_frac`, `max_fold_change`, and `max_wells`. The QC analyzes the fraction of all viral-barcode counts in the no-serum samples (wells) that are attributable to each viral barcode, and checks that this fraction is at least `min_frac` and is not more than `max_fold_change` different from the median fraction for this viral barcode across no-serum samples. If a given viral barcode fails either of these filters in at least `max_wells` wells, it is dropped entirely from that plate.

 - `per_neut_standard_barcode_filters`: has subkeys `min_frac`, `max_fold_change`, and `max_wells`. The QC analyzes the fraction of all neutralization-standard barcode counts in all samples (wells) that are attributable to each neutralization-standard barcode, and checks that this fraction is at least `min_frac` and is not more than `max_fold_change` different from the median fraction for this  barcode across all wells samples. If a given neutralization-standard barcode fails either of these filters in at least `max_wells` wells, it is dropped entirely from that plate.

 - `min_neut_standard_count_per_well`: drop any well where the total counts for neutralization standard barcodesis not at least this large.

 - `min_no_serum_count_per_viral_barcode_well`: drop any viral-barcode / no-serum well combination where the viral barcode does not have at least this many counts.

 - `max_frac_infectivity_per_viral_barcode_well`: drop any viral-barcode / well combinations where the viral barcode in that well has a computed fraction infectivity exceeding this.

 - `min_dilutions_per_barcode_serum_replicate`: drop any viral-barcode / serum-replicate combinations where the serum-replicate does not have at least this many dilutions for the viral barcode.


#### curvefit_params
This key defines some parameters specifying how the neutralization curves are fit, which is done using the Hill curves defined in the [neutcurve](https://jbloomlab.github.io/neutcurve/) package.

You typically want to use the [YAML anchor/merge](https://ktomk.github.io/writing/yaml-anchor-alias-and-merge-key.html) syntax to define a default that you then merge for specific plates.
The default can be defined like this:
```
default_process_plate_curvefit_params: &default_process_plate_curvefit_params
  frac_infectivity_ceiling: 1
  fixtop: [0.75, 1.0]
  fixbottom: 0
  fixslope: [0.8, 10]
```

The specific meaning of these curve-fitting parameters are as follows:

 - `frac_infectivity_ceiling`: ceiling to apply to all fraction infectivities before fitting curves. You may want to this to one to put a ceiling on all values >1. In principle, no values should be >1 in the absence of experimental noise.

 - `fixtop`: how to set the top plateau of the neutralization curve. You can set it to:
   - A list of length two specifying a reasonable range, such as `[0.75, 1.0]`, in which case the top is optimized within that range. **This is typically the recommended setting.**
   - A fixed value to fix to a specific number, typically 1. This is recommended if you want to force all curves to have a top plateau of one.
   - The value `false` if you want it to be a totally free parameter. This is not recommended as you can sometimes get spurious fits of a very large value when the data don't capture fully neutralization.

 - `fixbottom`: how to set the bottom plateau of the neutralization curve to this value. Like `fixtop`, it can be a length-two list, a fixed value, or `false`. Typically you should set it to 0 unless you have a good reason otherwise.

 - `fixslope`: how to set the slope of the neutralization curve. Like `fixtop`, it can be a length-two list, a fixed value, or `false`. If you don't know the "slope" of the neutralization curve, setting to `false` is a reasonable choice. However, in many cases it is preferable to set to a range that encompasses "reasonable" slopes. Note that what is "reasonable" will depend on the units of the concentration, but when they are serum dilutions a "reasonable" range is often `[0.8, 10]`.

#### curvefit_qc
This key defines some parameters on quality-control performed after the curve-fitting; viral-barcode / serum-replicate combinations that fail this QC are dropped.

You typically want to use the [YAML anchor/merge](https://ktomk.github.io/writing/yaml-anchor-alias-and-merge-key.html) syntax to define a default that you then merge for specific plates.
The default can be defined like this:
```
default_process_plate_curvefit_qc:  &default_process_plate_curvefit_qc
  max_frac_infectivity_at_least: 0
  goodness_of_fit:
    min_R2: 0.75
    max_RMSD: 0.05
  serum_replicates_ignore_curvefit_qc: []
  barcode_serum_replicates_ignore_curvefit_qc: []
```

The specific meanings of these QC parameters are:

 - `max_frac_infectivity_at_least`: drop any viral-barcode / serum-replicate combination that does not have a maximum frac infectivity across all concentrations of at least this value. Typically if you want to allow curves where the sera neutralize at all tested concentrations then you should set a value of 0. But you should set a value >0.5 if you want to require all sera to have a midpoint within the dilution range.

 - `goodness_of_fit`: drop any viral-barcode / serum-replicate combination where the curve fit does not have reasonable goodness of fit. A curve is dropped if it fails **both** of `min_R2` and `max_RMSD` (passing one is enough). The reason for using both is that when the data has more variation, we can tolerate a higher RMSD if the R2 is still good. There are two keys specified under `goodness_of_fit`:
  - `min_R2`: does curve fit have a [coefficient of determination](https://en.wikipedia.org/wiki/Coefficient_of_determination) at least this large (a coefficient of determination of 1 is a perfect fit). Used to drop very poor fitting curves. Reasonable values might be in the 0.6 to 0.8 range, although you should also just look at the curves being dropped to see if they look good.
  - `max_RMSD`: does curve fit have a root-mean square deviation (square root of mean residuals) no larger than this? Reasonable values might be in the 0.05 to 0.1 range.

 - `serum_replicates_ignore_curvefit_qc`: list of any serum replicates for which we ignore the curve-fitting QC for all viral barcodes.

 - `barcode_serum_replicates_ignore_curvefit_qc`: list (as `[barcode, serum_replicate]`) of viral-barcodes / serum-replicates where we ignore the curve-fitting QC.

#### illumina_barcode_parser_params
This key defines parameters for the [illuminabarcodeparser](https://jbloomlab.github.io/dms_variants/dms_variants.illuminabarcodeparser.html) that override anything set in the global `illumina_barcode_parser_params` above.
It is optional, and if not defined just the global params are used.
If this is defined, it is used to update the global params (adding new params and overriding any shared ones).
The main anticipated use case is if you add plate-specific indices in the round 1 PCR and want to specify those indices here using `upstream2` and `upstream2_mismatch`.

### default_serum_titer_as
Specifies how we compute the final titers in `serum_titers`.
Can be either `midpoint` or `nt50` depending on whether you want to report the value where the fraction infectivity gets to 50%, or the midpoint of the curve, so should be either
```
default_serum_titer_as: midpoint
```
or
```
default_serum_titer_as: nt50
```
The difference only becomes relevant if some your curves have plateaus substantially different than zero and one.

If you want to handle specific sera different, see `sera_override_defaults`.

### default_serum_qc_thresholds
Default QC we apply to each serum-virus pair when reporting out the final titers (medians across replicaes) in `serum_titers`.
Any serum-virus pair that fails this QC does not have a titer reported unless it is specified in `sera_override_defaults`.

Should look like this:
```
default_serum_qc_thresholds: &default_serum_qc_thresholds
  min_replicates: 2
  max_fold_change_from_median: 3
  viruses_ignore_qc: []
```
where:

 - `min_replicates`: drop any virus-serum titer that is not supported by at least this many replicates.
 - `max_fold_change_from_median`: drop any virus-serum titer where any replicate differs by more than this from the median across replicates.
 - `viruses_ignore_qc`: list of viruses for which you want to ignore the above QC. Specifying a virus here will ignore the QC for **all** sera, if you want to make a virus-serum specific exclusion then instead specify this in `sera_override_defaults`.

### sera_override_defaults
Override `default_serum_titer_as` or `default_serum_qc_thresholds` for specific sera in each group (recall groups are assigned per-plate).
For instance, this could look like:
```
sera_override_defaults:
  serum:
    M099d30:
      qc_thresholds:
        <<: *default_serum_qc_thresholds
        viruses_ignore_qc:
          - A/Belgium/H0017/2022
    Y044d30:
      qc_thresholds:
        <<: *default_serum_qc_thresholds
        max_fold_change_from_median: 4
      titer_as: nt50
```

The above means that in the group called *serum*, for serum `M099d30` we override the `default_serum_qc_thresholds` to exclude virus ` A/Belgium/H0017/2022`, and for serum `Y044d30` we override the defaults to allow a greater fold-change from median for individual replicates, and compute the titer as `nt50`.
Anything not listed here gets handled by the defaults in `default_serum_titer_as` and `default_serum_qc_thresholds`.

### miscellaneous_plates
This is an optional key that can be used specify plates that you just want to count barcodes for, and then analyze those counts outside the main pipeline.
This might be useful for library pooling or QC, for instance---or if you want to look at some failed plates that you don't actually want to fit curves for.

If you do not want to specify any miscellaneous plates either leave this key out or set it to an empty dictionary (`{}`).

The key should look like this:

```
miscellaneous_plates:

  <plate_name_1>:
    date: <date>
    viral_library: <viral library>
    neut_standard_set: <standard set>
    samples_csv: <filename>
    illumina_barcode_parser_params:  # optional key
        <parser params to override global>

  <plate_name_2>:
    ...
```

The plate name is just the name assigned to the plate.
The `date`, `viral_library`, `neut_standard_set`, and `illumina_barcode_parser_params` keys have the same meaning as for the plates specified under `plates`.

The `samples_csv` should specify the samples to analyze in a CSV that has columns named "well" and "fastq", and optionally other columns as well.

The output is that for each plate, the following files are created:

 - `results/miscellaneous_plates/<plate_name>/<well>_counts.csv`: counts of each viral barcode in that well of that plate.
 - `results/miscellaneous_plates/<plate_name>/<well>_invalid.csv`: counts of each invalid barcode in that well of that plate.
 - `results/miscellaneous_plates/<plate_name>/<well>_fates.csv`: summarizing number of reads that are valid and various types of invalid for each well of that plate.


## Results of running the pipeline
The results of running the pipeline are put in the `./results/` subdirectory of your main repo.
We recommend using the `.gitignore` file in [./test_example/.gitignore] in your main repo to only track key results in your GitHub repo.
The key results if the pipeline runs to completion are in `./results/aggregated_titers/titers_{group}.csv` for each group of sera.
The set of full created outputs are as follows (note only some will be tracked depending on your `.gitignore`):

  - Outputs related to barcode counting:
    - `./results/barcode_counts/`: files giving the barcode counts for each sample. You should track this in the repo.
    - `./results/barcode_fates/`: files giving the statistics (fates) of reads in the barcode counting for each sample. You do not need to track this in the repo as the results are plotted.
    - `./results/barcode_invalid/`: files giving counts of invalid barcodes for each sample. You do not need to track this in the repo, but it could be helpful to look at these identities in counts if QC shows you are getting many invalid barcodes.

  - Outputs related to processing each plate:
    - `./results/plates/{plate}/frac_infectivity.csv`: fraction infectivity for viral barcodes for a plate. You should track this in the repo.
    - `./results/plates/{plate}/process_{plate}.ipynb`: Jupyter notebook processing counts for a plate. You do not need to track this as an HTML version will be rendered in `./docs/` when pipeline runs successfully.
    - `./results/plates/{plate}/process_{plate}.html`: HTML of Jupyter notebook processing counts for a plate. You do not need to track this as it will be rendered in `./docs/` when pipeline runs successfully.
    - `./results/plates/{plate}/qc_drops.yml`: details on data (barcodes, wells, etc) dropped for failing QC when processing this plate.
    - `./results/plates/{plate}/curvefits.csv`: the neutralization curve fits to each serum on each plate. You should track this in repo.
    - `./results/plates/{plate}/curvefits.pickle`: pickle file with the `neutcurve.CurveFits` object for the plate. You do not need to track this in the repo as both the plots and numerical data are rendered elsewhere.

  - Output related to per-serum titers (aggregated across replicates potentially run on different plates); note that serum are organized per-group as specified in the plates:
    - `./results/sera/groups_sera_by_plate.csv` summarizes which plate(s) each group/serum was run on.
    - `./results/sera/{group}_{serum}/titers.csv`: titer for each virus against the group/serum, reported as the median across replicates, and only keeping those that pass QC. You should track this file in the repo.
    - `./results/sera/{group}_{serum}/titers_per_replicate.csv`: titers for each replicate of each virus against the group/serum. You should track this file in the repo.
    - `./results/sera/{group}_{serum}/curves.pdf`: PDF rendering of the neutralization curves for the group/serum. You do not need to track this in the repo as a HTML version of a notebook containing the plots is tracked in `./docs/`.
    - `./results/sera/{group}_{serum}/curvefits.pickle`: pickle file with the `neutcurve.CurveFits` object for this group/serum, after applying QC filters. You do not need to track this in the repo as both the plots and numerical data are rendered elsewhere.
    - `./results/sera/{group}_{serum}/{group}_{serum}_titers.ipynb`: Jupyter notebook that aggregates titers for a group/serum across all plates. You do not need to track this in the repo as a HTML version of the notebook is tracked in `./docs/`.
    - `./results/sera/{group}_{serum}/{group}_{serum}_titers.html`: HTML rendering of the Jupyter notebook that aggregates titers for a group/serum across all plates. You do not need to track this in the repo as it will be rendered in `./docs/` when the pipeline runs successfully.
    - `./results/sera/{group}_{serum}/qc_drops.yml`: virus-group/serum titers dropped due to QC when processing this serum's titers.

  - Results related to aggregated titers across all sera in a group after applying all quality control:
    - `./results/aggregated_titers/titers_{group}.csv`: titers for all sera / virus in a group (median of replicates). You should track this file as it has the final processed results.
    - `./results/aggregated_titers/curvefits_{group}.pickle`: pickle file with the `neutcurve.CurveFits` object holding all final curves for a group. You do not need to track this in the repo, but if you have further code that makes specific plots you may want to use this.
    - `./results/aggregated_titers/titers.html`: interactive plot of titers for all sera. You do not need to track this in the repo as it is rendered in `./docs/` when the pipeline runs successfully.
    - `./results/aggregated_titers/aggregate_titers.ipynb`: Jupyter notebook that aggregates all the titers. You do not need to track this in the repo.

  - Results summarizing data dropped due to QC:
    - `./results/plate_qc_drops.yml`: YAML file summarizing all data (barcodes, wells, etc) dropped during the plate-processing QC. You should track this in repo.
    - `./results/groups_sera_qc_drops.yml`: YAML file summarizing all group/serum-virus titers dropped during the serum titers QC. You should track this in repo.
    - `./results/aggregate_qc_drops.ipynb`: Jupypter notebook summarizing the QC drops. You do not need to track as an HTML version is rendered in `./docs/`
    - `./results/aggregate_qc_drops.html`: HTML version Jupypter notebook summarizing the QC drops. You do not need to track as it is rendered in `./docs/`

## Examining the output and setting appropriate QC values in the configuration
When you run the pipeline, the QC values in the configuration will be automatically applied, and HTML notebooks summarizing the processing of each plate and sera are rendered in `./docs`, alongside a summary of all QC across all plates / sera.
YAML summaries of the QC are also created.

While the QC is designed to hopefully make reasonable default choices, you should **always** carefully look through these notebooks after adding new data, and potentially adjust the QC in the configuration and re-run.

## Rendering HTML plots and notebooks in docs
If the pipeline runs to completion, it will create HTML documentation with plots of the overall titers, per-serum titer analyses, per-plate analyses and overall QC summary in a docs subdirectory, which will typically named be `./docs/` (if you use suggested key in configuration YAML).
This HTML documentation can be rendered via [GitHub Pages](https://docs.github.com/en/pages/getting-started-with-github-pages/configuring-a-publishing-source-for-your-github-pages-site) from the `./docs/` directory.

Looking at this documentation is a good way to QC the data and understand the results.

The documentation for the test example for this pipeline is at [https://jbloomlab.github.io/seqneut-pipeline/](https://jbloomlab.github.io/seqneut-pipeline/).

If you want to add additional HTML files to the docs, specify a dict in the top-level `Snakefile` with the name `add_htmls_to_docs` like this:
```
add_htmls_to_docs = {
    "Additional files": {
        "Example HTML file":  "results/extra_htmls/example_html.html",
        <other keys specifying file names and their paths>
    },
    <other nested dicts with a heading and then name: file key-value pairs>
}
```

## Test example and testing via GitHub Actions
The [./test_example](test_example) subdirectory contains a small test example that illustrates use of the pipeline.

The code is tested by running this example, as well as formatted with [black](https://github.com/psf/black) and [snakefmt](https://github.com/snakemake/snakefmt) and linted with [ruff](https://github.com/astral-sh/ruff) and [snakemake --lint](https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html) via the GitHub Action specified in [.github/workflows/test.yaml](.github/workflows/test.yaml).

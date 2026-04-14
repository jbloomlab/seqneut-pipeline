# CHANGELOG

### version 6.3.0
- Update package versions (helps fix `snakemake` error when using `nextstrain-prot-titers-tree` by matching versions):
  + marimo: 0.17.6 -> 0.23
  + markdown: 3.9 -> 3.10
  + pillow: 12.0 -> 12.2
  + ruamel.yaml: 0.18.15 -> 0.19
  + snakemake: 9.13.2 -> 9.19

#### version 6.2.1
- Set target version for `black` to Python 3.13 in `pyproject.toml`.
- Fix divide by zero error if no barcode counts ([this pull request](https://github.com/jbloomlab/seqneut-pipeline/pull/85/changes) addressing [this issue](https://github.com/jbloomlab/seqneut-pipeline/issues/83)).

### version 6.2.0
Allow one level of nesting in `add_htmls_to_docs`.

### version 6.1.0
Add more information to the invalid barcodes CSV produced by *count_barcodes*, specifying the closest valid barcode, the Hamming distance to the closest valid barcode, and the counts of the closest valid barcode. Also for valid and invalid barcodes report fraction of all barcode reads (valid and invalid) that map to that barcode. Addresses [this issue](https://github.com/jbloomlab/seqneut-pipeline/issues/81).

## version 6.0.0
Change how docs summarizing the output are tracked on GitHub Pages.
This change is designed to solve the issue where the main repository is bloated by committing of the docs, which both makes the repo very large and makes it hard to see code changes vi git diffs.
The following changes were made:
  - Docs are no longer written to `./docs` that is tracked via git; instead they are written to `./results/docs` which is not tracked via git.
  - Anytime you want to commit the current docs, you run the script [publish_docs_gh-pages.sh](publish_docs_gh-pages.sh) which commits the current contents of `./results/docs` to a *gh-pages* branch with no history; you then set GitHub Pages to serve from the *gh-pages* branch and the `/root` directory on GitHub.

So to migrate an existing repository, do the following:
  1. Remove any existing `./docs/` directory (`rm -rf docs`).
  2. Remove any lines in you `.gitignore` that involves `docs` (such as `docs/*` or `!docs/*.html`).
  3. Remove the *docs* key from your `config.yml`, it is no longer needed.
  4. (If you are migrating from version 4 or earlier) Set the *curve_display_method* in `config.yml` to a valid value (suggestion is *svg*); create this option if it does not already exist.
  5. Anytime you want to commit the current docs, run `./seqneut-pipeline/publish_docs_gh-pages.sh`.
  6. On your GitHub Repository, go to *Settings* then *Pages* on the left toolbar and set GitHub Pages to serve from the *gh-pages* branch and the `/root` directory.

## version 5.0.0
Migrated to use `marimo` notebooks rather than Jupyter notebooks. Main changes:
  - all notebooks are now `marimo`, not Jupyter; these notebooks are now run by a driver script from a context pickle.
  - in the docs, the notebooks render much better, including hiding code by default (you can click on an option in upper left to show code).
  - the `conda` environment is updated to newer package versions and to add `marimo` and drop `jupyterlab`.
  - the *curve_display_method* in `config.yaml` is now required, and has possible values of *svg*, *pdf*, *png8*, and *inline* (see repo README for details on each).
  - project now has a single consolidated `pyproject.toml` that includes version number.

#### version 4.0.2
Provide a better error message when there are duplicated samples.

#### version 4.0.1
Fix bug in how docs show plots added via `add_htmls_to_docs`.

## version 4.0.0
Primarily includes updates to make pipeline run more efficiently and produce more manageable output on large pipelines. Specifically:
 - Update to [newer version of `neutcurve`](https://github.com/jbloomlab/neutcurve/pull/63) that plots large panels of neutralization curves faster.
 - Some changes to `process_plate` and `group_serum_titers` notebooks to accelerate plotting and improve plotted neutralization curves.
 - In `process_plate`, only plot virus-serum pairs that fail QC in that part, rather than all viruses or sera with any pair failing QC (see [here](https://github.com/jbloomlab/seqneut-pipeline/issues/65)).
 - Fix bug in `agregate_titers` when *viral_strain_plot_order* is *None*.
 - Add `curve_display_method` option to display neutralization curves in notebooks at smaller size; this involves adding the `notebook_funcs.py` file as part of the pipeline internally. By default the method is now `png8` which is lower resolution but leads to much smaller notebooks in the results and docs; use `inline` for the older bigger but higher resolution plots.
 - Shrink size of `altair` plots in `process_plate` by making them use data more efficiently.
 - Add a file summarizing QC drops at barcode level (`results/qc_drops/barcode_qc_drops.yml`) (addresses [this issue](https://github.com/jbloomlab/seqneut-pipeline/issues/69))
 - Make the aggregated titers chart produced by pipeline in docs size in a way that is reasonable when many sera and viruses: rather than a different plot for each serum, there is just a median chart and an all-sera chart with individual sera selectable on legend.

### version 3.4.0
- Make the pipeline run OK even if FASTQs are not available if the counts have already been computed (see [here](https://github.com/jbloomlab/seqneut-pipeline/pull/62)). This is designed to make it more portable. Specifically:
  - default `.gitignore` and README now suggest tracking the barcode fates files as well as the counts files, as they are needed by downstream pipeline steps.
  - the invalid barcodes are no longer considered an output of the pipeline as they are not tracked.

### version 3.3.0
- In internal functioning of `seqneut-pipeline.smk`, make *viral_libraries* and *neut_standard_sets* Python variables passed as params rather than CSVs passed as input file (see [here](https://github.com/jbloomlab/seqneut-pipeline/pull/60)). The advantage of this is altering the CSVs holding the viral libraries and neutralization standards now only triggers a pipeline re-run if the relevant content has changed rather than always. Specifically:
  - *viral_libraries* is now a dict of data frames rather than a dict of CSVs.
  - *neut_standard_sets* is now a dict of dataframes rather than a dict of CSVs.

### version 3.2.0
- Update software versions in `conda` environments (see [here](https://github.com/jbloomlab/seqneut-pipeline/pull/58)).

#### version 3.1.4
- Minor bug fix to make it possible to run `miscellaneous_plates` without running any `plates`

#### version 3.1.3
- Remove `defaults` channel from `conda` envs to avoid issue with Anaconda licensing.
- Update `snakemake` version.

#### version 3.1.2
- Minor bug fix to resolve issue with specifying serum_replicates and barcode_serum_replicates to manual drops, similar issue resolved in 3.1.1. Addresses [this issue](https://github.com/jbloomlab/seqneut-pipeline/issues/49).

#### version 3.1.1
- Minor bug fix to make it possible to add `wells` or `barcodes` to the `manual_drops` specified for each plate. Addresses [this issue](https://github.com/jbloomlab/seqneut-pipeline/issues/45).

### version 3.1.0
- Configured to enable plate-level indices to be embedded in the round-1 PCR primers (see [this issue](https://github.com/jbloomlab/seqneut-pipeline/issues/40)). Essentially, this amounts to allowing a per-plate flanking sequence to be specified for each plate, and only FASTQ reads with that flanking sequence are read for that plate. Typically this index would be specified as `upstream2` in the [illuminabarcodeparser](https://jbloomlab.github.io/dms_variants/dms_variants.illuminabarcodeparser.html). To enable this change, altered the configuration from the previous setup of just having a single global `illumina_barcode_parser_params` applied to all plates. Now such a global parser is still specified that has default values that you want to apply to all plates. But in addition, in the per-plate configuration you can specify `illumina_barcode_parser_params` that are added to (and override) anything in the global parser params, and can contain plate specific `upstream2` and other relevant setting (eg, `upstream2_mismatch`). The test example was modified to use this option for plate2 and plate11.

- Update software versions:
  - `dms_variants` to 1.6.0
  - `neutcurve` to 2.1.0
  - `altair` to 5.3
  - `python` to 3.12

- Draw neutralization curves using `draw_in_bounds=True` with `neutcurve` to avoid lines extrapolating beyond data. Addresses [this issue](https://github.com/jbloomlab/neutcurve/issues/59).

## version 3.0.0
- In `curvefit_params` in the YAML configuration, now `fixslope` should be specified in addition `fixtop` and `fixbottom`. In addition, all three of these can be set to constraint ranges rather than just totally free or to fixed values. Alongside this change, the slope of curve fits are now reported in key output files. Addresses [this issue](https://github.com/jbloomlab/neutcurve/issues/53) and [this issue](https://github.com/jbloomlab/seqneut-pipeline/issues/32).
  - This is a **backward-incompatible change** in the configuration YAML, now you must specify `fixslope` under `curvefit_params`.

- In `process_plate_curvefit_qc` in the YAML configuration, there is a new key called `goodness_of_fit` and now both `min_R2` (the minimum coefficient of determination) and `max_RMSD` (the maximum mean square deviation) for each curve fit are specified as keys under that. The curves are then filtered to retain only those that meet *either* of these criteria (so must fail both to be dropped). Addresses [this issue](https://github.com/jbloomlab/seqneut-pipeline/issues/33) and [this issue](https://github.com/jbloomlab/neutcurve/issues/55#issuecomment-2016975219). Alongside this change, the `rmsd` is now reported in key output files. Also, in the tabulation of failures, `fails_min_R2` now becomes `fails_goodness_of_fit`.
  - This is a **backward-incompatible change** in the configuration YAML. Previously `min_R2` was a standalone key under `process_plate_curvefit_qc`; now `goodness_of_fit` is the required key and `min_R2` and `max_RMSD` are required keys under it.

- Handle titers that are outside the range of the dilutions series by reporting them as upper or lower bounds rather than as interpolated, and marking them appropriately on plots. This change helps with low potency or high potency sera, where there may be no neutralization or high neutralization at all tested concentrations. Addresses [this issue](https://github.com/jbloomlab/seqneut-pipeline/issues/30).

- Each plate is now assigned to a *group*, which makes it possible to have separate groups (for instance, "serum" and "pilot" if you have serum samples of interest and pilot experiments, although it can be everything). **This is a backward-incompatible change** that requires you to update the configuration YAML and changes the names of some output files (so you will need to update your `.gitignore` to be similar to the new one in the `test_example`). Specifically:
  - For each plate under `plates` in the configuration YAML, you now specify a `group` as one of the keys (eg, serum, pilot, etc)
  - For `sera_override_defaults` in the configuration YAML, the keys for individual sera are now nested under keys for their groups.
  - The sera are processed by group, so "group" is now a column in the output CSVs and the serum results files are now in subdirectories named `./results/sera/{group}_{serum}` rather than `./results/sera/{serum}` as before.
  - The aggregated titers are now in per-group CSVs with names like `./results/aggregated_titers/titers_{group}.csv` rather than in the single `./results/aggregated_titers/titers.csv` from before.
  - The final aggregated output plot allows you to select by group.
  - The docs are organized by group in the per-plate and per-sera plots.

- Added another plate (of H3N2 rather than H1N1) to the `test_example` to test some of the changes introduced in this version.

- Update `seqneut-pipeline` conda environment in `environment.yml`. Update `neutcurve` 2.0.1, also update other packages (`pandas`, `snakemake`, `markdown`, `papermill`) to latest versions.

### version 2.2.0
- Add the `add_htmls_to_docs` option, which can be specified in `Snakefile` to add additional HTML documentation to pipeline.
- Update to `snakemake` 8.5.4.

### version 2.1.0
- Add an option to specify `miscellaneous_plates` which are plates that just have their barcodes counted (addresses [this issue](https://github.com/jbloomlab/seqneut-pipeline/issues/26)).

#### version 2.0.1
- Update to `dms_variants` 1.5.0 (addresses [this issue](https://github.com/jbloomlab/seqneut-pipeline/issues/24)).

## version 2.0.0
Full re-write that changes how configuration is specified to automatically do the QC, and uses a newer version of `neutcurve` that fits better. Completely backward-incompatible with version 1.*.

### version 1.1.0
- Update software environment to include `neutcurve` 1.0 and `snakemake` 8.0.
- Fix linting of notebooks with `ruff`

## version 1.0.0
The initial version of this pipeline was created from prior code from Andrea Loes and Will Hannon.
The pipeline was re-factored into this initial version with the following design goals:

  - Separate the pipeline from the configuration and code for specific analyses. This is important because we want to be able to have an portable pipeline that can be used for multiple different studies, and can be independently tested and versioned. This will be especially important of we eventually imagine this assay being used widely in large-scale studies where standardization is important.

  - As part of the modularization of the pipeline from the specific analyses, the configuration and code is designed now solely to focus on calculating neutralization titers for sera. Additional configuration such as grouping sera by individuals for different timepoints etc will not be a universal feature of such studies, so is moved to the upstream project-specific code that runs the pipeline.

  - Rename from names like "NGS neuts" to "seqneut", because in paper we are choosing to describe as sequencing-based neutralization assays.

  - Perform analyses and provide configuration on a per-plate basis, and (except for last output steps) do not aggregate configuration or analyses across plates. This is because we envision typical studies as sequentially running and analyzing plates, each of which will be QC-ed separately.

  - Move all configuration, sample exclusion, etc in a YAML configuration file. Samples will not be dropped by the code or have thresholds applied that are not specifically delineated in this configuration. This makes it easier to see what is being done in the analysis since a user can look at just the configuration YAML file without examining the code. **The principle here is that all QC and/or sample or barcode exclusions must be transparently specified here, and is not done silently in the code.** This ensures any problems with a plate will be flagged at time of analysis and an intentional decision has to be made whether to accept or reject the data based on the QC thresholds.

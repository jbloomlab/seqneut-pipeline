# CHANGELOG

## version 8.0.0
- Update software versions, most notably `python` to 3.14:
  + `python`: 3.13 -> 3.14
  + `snakemake`: 9.19 -> 9.23
  + `pandas`: pin loosened from 2.3.3 to 2.3 for consistency with the other pins (still resolves to 2.3.3). Note that `pandas` cannot yet be updated to 3.0 because the `bioconda` `snakemake` package pins `pandas <3`.
  + Drop `mamba` from `environment.yml`, as `snakemake` has deprecated the `mamba` conda frontend and now relies on `conda` (>24.7.1) for conda-environment deployment.
  + In `pyproject.toml`, set `requires-python` to `>=3.13` (was `>=3.9`) and add `py314` to the `black` target versions.
  + Update the `snakemake` version badge in the README to `>=9`.

- Update the GitHub Actions used for testing in `.github/workflows/test.yaml` to their latest versions:
  + `actions/checkout`: v5 -> v7
  + `conda-incubator/setup-miniconda`: v3 -> v4, and rename the now-deprecated `auto-activate-base` input to `auto-activate`
  + `actions/upload-artifact`: v4 -> v7
  + Also drop the `rm -rf docs` from the test-example step, as it is a leftover from version 6.0.0 when the docs moved from `./docs` to `./results/docs`.

- Convert the `date` of each plate in the config to a string (in the new `stringify_plate_dates` in `funcs.smk`). YAML parses unquoted dates such as `2023-08-01` into `datetime.date` objects, and the updated `snakemake` hashes the config with `json.dumps` at startup, which fails on those objects with `TypeError: Object of type date is not JSON serializable`. Previously the dates were only converted to strings on a copy of the config, and not at all for `miscellaneous_plates`.

- Adopt the much broader default rule selection of `ruff` 0.16, and fix or explicitly ignore everything it newly reports in the notebooks and scripts. In `pyproject.toml`, `B018` and `PLR1711` are ignored for `notebooks/*.py`, as they fire only on the file format that `marimo` writes rather than on anything actually wrong: `PLR1711` flags the `return` that `marimo` generates at the end of every cell, and `B018` flags a cell's trailing bare expression, which is what that cell renders.

- Add plate-to-plate correlations of titers to the per-serum notebook. For a serum measured on more than one plate, there is now a scatter plot for each pair of plates showing the titer of each viral strain on one plate versus the other, labeled with the Pearson correlation of the log titers and arranged in a grid of up to four columns when there are several pairs. The titer for a strain on a plate is the median over that plate's barcodes. All strains are shown colored by whether the per-serum QC on replicate-to-replicate variation retained or dropped them, and a `show strains` dropdown below the plot restricts it to just the retained or just the dropped strains; the Pearson correlation is reported both over all strains and over just the retained ones, so the effect of the QC is visible without a separate plot.

- Add a per-plate view of the plate-to-plate correlations, organized by plate rather than by serum. A new `plate_to_plate_corr` rule creates a notebook and a CSV of Pearson correlations for each plate, comparing its titers with those on each other plate of its group, and adds the notebooks to the docs under a subheading per group. The correlations annotating the plots are computed by `altair` transforms from the points currently shown, so they follow the dropdowns that select on the per-serum QC status and on serum.

- Add a `plate` column to `./results/plates/{plate}/curvefits.csv` and `./results/sera/{group}_{serum}/titers_per_replicate.csv`. The per-serum analysis needs to know which plate each replicate came from, and recording it upstream avoids parsing it out of the replicate name, which is ambiguous when one plate name is a prefix of another (as it now is in the test example).

- Add a `dropped_by_qc` column to `./results/sera/{group}_{serum}/titers_per_replicate.csv` indicating whether the virus was dropped by the per-serum QC on replicate-to-replicate variation. That file remains unfiltered by this QC, so the column is what makes it possible to tell which of its titers made it into the final results. Note that the file is still filtered by the QC applied when each plate was processed, and that this QC drops a virus rather than an individual replicate, so the column has the same value for all replicates of a virus. Writing the file now happens after the per-serum QC is determined rather than before it.

- Add a `plate2_exact_dup` plate to the test example, an exact duplicate of `plate2` that gives its sera a third plate so that the multi-panel correlation grid is tested. Being an exact duplicate also makes two of the three panels serve as checks: its titers must correlate perfectly with `plate2`, and its correlation with `plate11` must match that of `plate2` with `plate11`. This changes the number of replicates per titer for that group, and so also the expected titers in `expected_titers_for_test.csv`.

- Simplify [./test_example/.gitignore](test_example/.gitignore), which the README recommends copying into your own repo, using `results/**` plus `!results/**/` to ignore all of `results` while still allowing the key results to be re-included. This removes the need to re-include each subdirectory before ignoring its contents. Which files are tracked is unchanged (verified over all files in the test example's `results`), so you do not need to update your own copy.

- Add an optional `collapse_strain_barcodes` configuration key. When it is set to true, the counts of all of a strain's barcodes are summed in the processing of each plate before any of the QC or curve fitting, so that a plate gives one neutralization curve per strain rather than one per barcode.

  The rationale against collapsing is that having multiple barcodes per strain provides internal replicates.
  The rationale for collapsing barcodes is that the noise in the measurements is largely determined by the finite number of virions of each barcoded strain that go in each well, so collapsing barcodes can increase the total input virions for each strain and so reduce experimental noise.
  Which rationale wins out can be experiment-dependent, depending for instance on the evenness of the barcode distribution and the number of strains and barcodes in the library, so we do not have a universal recommendation on whether collapsing is or is not a good idea.

  This option defaults to false, which analyzes each barcode separately as before the option was introduced. Note that `min_replicates` in the serum QC thresholds then counts only the plates a serum was measured on, so it typically has to be set to one; see the README for this and the other implications of collapsing. The test example is also run a second time with the barcodes collapsed, using the new `config_collapse_strain_barcodes.yml` overrides.

### version 7.1.0
Update software versions (reason is to update `neutcurve` to fix bug causing scrunched panels in curve plots):
 - `altair`: 5.5 -> 6.2
 - `pillow`: 12.2 -> 12.3
 - `neutcurve`: 2.3.0 -> 2.3.1

## version 7.0.0
Add support for titrating a `concentration` (e.g. an antibody at some µg/ml) as an alternative to a serum `dilution_factor`.
A plate can now set `dilution_factor_or_concentration: concentration` (default is `dilution_factor`, so existing configs are unchanged) together with a `concentration_units` value (e.g. `ug/ml`).
Such a plate's `samples_csv` must have a `concentration` column instead of `dilution_factor`.
For a concentration, the value is used directly as the curve-fitting concentration (rather than taking the reciprocal), and the reported titer is the fitted IC50 (or midpoint) in the given units, with `serum_titer_as` accepting `ic50` or `midpoint`.
All plates in a group must share the same `dilution_factor_or_concentration` and `concentration_units`.

This is a major version bump because of two output changes that affect **all** projects (not just those using concentrations):
  - The aggregated interactive titers plot is now written as **one plot per group** at `results/aggregated_titers/titers_{group}.html`, replacing the previous single combined `results/aggregated_titers/titers.html`. The GitHub Pages docs now link to one titers plot per group. (Besides being cleaner, this avoids plotting incomparable titers, such as reciprocal dilutions vs. concentrations, on a shared axis.)
  - Every titers CSV now has a `titer_units` column (`reciprocal_dilution` in dilution-factor mode, or the plate's `concentration_units` in concentration mode).

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

# CHANGELOG

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

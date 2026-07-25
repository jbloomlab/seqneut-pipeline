# CLAUDE.md - Context for seqneut-pipeline

The [README](README.md) is the canonical reference for what this pipeline is, how to
configure it, and what it produces: every config key, every QC threshold, the
curve-fitting parameters, the full list of result files, how to include the pipeline as
a submodule in a project repo, and the format/lint commands. Read it first, and do not
duplicate it here.

This file covers only what the README does not: how the code is organized internally,
and the conventions and invariants that cannot be inferred from the config
documentation.

## Rule map

`seqneut-pipeline.smk` defines the rules below, and project repos include it from their
top-level `Snakefile`. Most of the actual analysis lives in the `marimo` notebooks
rather than in the rules.

| rule | implemented by | does |
|---|---|---|
| `count_barcodes` | `scripts/count_barcodes.py` | parses FASTQs into viral and neut-standard barcode counts using `dms_variants.illuminabarcodeparser` |
| `process_plate` | `notebooks/process_plate.py` | per plate: applies `qc_thresholds`, computes fraction infectivity, fits curves with `neutcurve` |
| `groups_sera_by_plate` | `scripts/groups_sera_by_plate.py` | maps each group/serum to the plate(s) it ran on. A checkpoint because its input is the plates' `curvefits.csv`, so which group/serum jobs exist is only known once the plates have been processed and their QC applied |
| `group_serum_titers` | `notebooks/group_serum_titers.py` | per group/serum: medians titers across replicates from all of its plates, applies the serum-level QC |
| `aggregate_titers` | `notebooks/aggregate_titers.py` | combines all sera in a group and builds the interactive titer plots |
| `aggregate_qc_drops` | `notebooks/aggregate_qc_drops.py` | summarizes every QC drop across all plates and sera |
| `build_docs` | `scripts/build_docs.py` | collects the notebook HTML into `results/docs` |
| `miscellaneous_plate_count_barcodes` | `scripts/count_barcodes.py` | counting only, for `miscellaneous_plates` |

`funcs.smk` holds the config-processing helpers that run while the rules are being
defined (`process_plate`, `process_miscellaneous_plates`, `stringify_plate_dates`, and
the sample-table validation). It relies on `pd` and `re` being imported by
`seqneut-pipeline.smk`.

## Marimo notebooks

The analysis notebooks are `marimo`, not Jupyter (since v5.0.0).
`scripts/run_marimo_w_context_pickle.py` pickles the `snakemake` context (input, output,
params, wildcards, threads, resources) and runs `marimo export`, so a notebook reads its
inputs from that pickle rather than from a global `snakemake` object. The same notebooks
can be opened with `marimo edit` for interactive work.

Two properties of the format explain why `pyproject.toml` ignores `B018` and `PLR1711`
for `notebooks/*.py`:

 - `marimo` writes every cell as a function ending in a `return` of the variables that
   cell exports (see `to_functiondef` in `marimo/_ast/codegen.py`), which trips
   `PLR1711`. Deleting those returns is futile, as `marimo` regenerates them on save.
 - A cell's trailing bare expression is what that cell renders, which trips `B018`.
   Deleting them would remove plots from the docs.

## Naming conventions

Built in `funcs.smk` (~lines 155-195) and relied on throughout the notebooks:

 - `serum_replicate`: `{serum}`, or `{serum}-{replicate}` when a serum has more than one
   replicate on the plate.
 - `sample`: `{plate}_{serum_replicate}_{dilution_factor}` (or `_{concentration}`), with
   the trailing part omitted for no-serum samples.
 - `plate_replicate`: `{plate}`, or `{plate}-{replicate}` for multi-replicate sera.
 - `plate_barcode`: `{plate_replicate}-{barcode}`, which lets one barcode be tracked as a
   distinct replicate across plates. This is what is passed to `neutcurve` as
   `replicate_col`.

## Fraction infectivity

Computed in `notebooks/process_plate.py` as:

    frac_infectivity_raw = count / neut_standard_count / median_no_serum_ratio

That is, each viral barcode's count is normalized by its well's total neut-standard
counts, then divided by the median of that same ratio across the plate's no-serum wells.
The `frac_infectivity_ceiling` from `curvefit_params` is then applied as a separate
column; the curve fits use the ceilinged values, while the raw values are also saved.

## Invariants enforced in code

These raise errors, and are not all stated in the README:

 - Groups cannot contain `_` or `|`, as both break wildcards (`seqneut-pipeline.smk`).
 - Every plate needs at least one `serum: none` sample (`funcs.smk`). No-serum wells are
   designated by `serum: none` with a blank dilution factor or concentration, never by a
   value of 0.
 - All barcodes on a plate, viral library and neut standard set together, must be the
   same length, since `bclen` is inferred from the first one
   (`scripts/count_barcodes.py`).
 - Plate `date` values must be strings by the time `snakemake` starts executing, because
   it hashes the config with `json.dumps` and YAML parses an unquoted `2023-08-01` into a
   `datetime.date`, which is not JSON serializable. `stringify_plate_dates` in
   `funcs.smk` normalizes this for both `plates` and `miscellaneous_plates`; do not
   remove that call.

## Working in this repo

 - Verify changes with the format and lint commands given in the README, plus
   `cd test_example && snakemake -n`, which catches config and DAG errors cheaply.
 - Ask before running the full `test_example` pipeline. It deletes and regenerates
   ~300 git-tracked files under `test_example/results/`, and a rerun yields small
   floating-point differences that should not be committed.
 - Describe changes in `CHANGELOG.md` at a high level, including the rationale for
   anything surprising, but leave per-file and per-rule detail to the commit message.
   Bump the version in `pyproject.toml`, which is the single source of truth for it.

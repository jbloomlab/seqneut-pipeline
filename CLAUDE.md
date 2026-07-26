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

Two conventions of `seqneut_pipeline_outputs`, the list of final targets at the bottom of
`seqneut-pipeline.smk`:

 - Notebook HTML is not listed. Every notebook HTML is an input to `build_docs`, whose
   `docs` directory is listed, so requesting the HTML as well is redundant.
 - The non-HTML results of a rule are listed before `build_docs.output.docs`, so that the
   numerical outputs read as a group ahead of the documentation build.

Prefer letting a notebook handle a degenerate case over precomputing which jobs are
worth running. A rule that runs per group runs for every group, and its notebook reports
that there is nothing to show when, say, the group has one plate; deriving the subset of
groups in the `.smk` instead would add config-processing code for no gain. The same
applies to adding a config key to turn an analysis off: do not add one unless asked.

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

## Keeping Altair plots small

`altair` embeds the whole data frame as inline JSON in the chart spec, so the size of the
notebook HTML scales as rows x columns, with every string value repeated on every row.
The plots here are large enough that this matters, so keep the plotted data frame narrow
and let `vega-lite` reconstruct the rest client-side:

 - Any column that is functionally determined by a key already in the data frame belongs
   in a small lookup table, not in the plotted frame. Drop it and re-join it with
   `.transform_lookup(lookup=key, from_=alt.LookupData(lookup_df, key=key, fields=[...]))`.
   Assert that the key is unique in the lookup table so the join is 1:1 and lossless.
 - Ship columns wide and fold them with `.transform_fold()` rather than melting in
   `pandas`, which would multiply the row count of the embedded data.
 - A column that exists only to filter rows to a subplot of an `alt.concat` is not needed:
   slice the frame in `pandas` and give each subplot its own. That costs no extra rows,
   because `altair` consolidates every frame into the chart's top-level `datasets` and the
   rows are partitioned among the subplots either way, and it drops both the column and the
   `transform_filter`.
 - Round floats that came from arithmetic before plotting. A median over an even number of
   values gives `210.35000000000002`, which costs 18 characters on every row for precision
   far below a pixel.
 - A column that is the same on every row (`group` and `serum` in a per-serum notebook, or
   `titer_units`, which follows from `titer_as`) belongs in a `transform_calculate` that
   returns a literal, so the value appears once in the chart rather than once per row. Note
   that the quoting is `transform_calculate(titer_units=f"'{units}'")`: the argument is a
   `vega` expression, so a string constant needs quotes inside the Python string. As with a
   lookup, the tooltip must then name the type explicitly (`alt.Tooltip("titer_units:N")`).
   `notebooks/group_serum_titers.py` does all three of these things, via `narrow_for_altair`.

`notebooks/process_plate.py` (~lines 1237-1285) is the worked example: `strain` is
determined by `barcode`, and `serum_replicate` and the dilution factor / concentration
are determined by `well`, so all three are dropped in favor of a barcode lookup and a
well lookup. Note that tooltips do not follow automatically; a looked-up field has no
`pandas` dtype for `altair` to infer from, so it must be listed explicitly with its type
suffix (`alt.Tooltip("strain:N")`).

## Interaction in concatenated Altair charts

A param belongs to whichever subplot declares it, so in an `alt.concat` of panels it must
be declared on the concatenated chart with `.add_params(...)` to be in scope for all of
them. Two consequences, both worked through in the plate-to-plate correlation plots of
`notebooks/group_serum_titers.py` and `notebooks/plate_to_plate_corr.py`:

 - `bind="legend"` cannot drive a selection across panels, as the legend is drawn by one
   panel. Use an input binding (`alt.binding_select`) instead, which the concatenated
   chart can declare; note it renders below the plot, not above it.
 - A mouseover selection needs `empty=False`, or every point matches when nothing is
   hovered. `vega-lite` puts this on the predicate rather than the selection, so it
   appears in the emitted spec as `{"param": ..., "empty": false}` in each condition and
   not in the param definition.

Since interaction cannot be exercised headlessly, verify it by parsing the emitted spec
out of the exported HTML, checking that each param appears once with a `views` list
covering every panel.

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
 - A new or changed notebook can be run against the committed `test_example/results`
   without running the pipeline: build the context dict that
   `scripts/run_marimo_w_context_pickle.py` would build (`workdir` set to
   `test_example`, `input` / `output` / `params` / `wildcards` filled in, outputs pointed
   somewhere scratch), pickle it, and run
   `marimo export html <notebook> -o <html> -- --context-pickle <pickle>` with
   `test_example` as the working directory. Synthesizing input CSVs this way also reaches
   cases the test example does not contain. `scripts/build_docs.py` can likewise be
   exercised by supplying a stub object as the global `snakemake`. Prefer this over
   copying cells into a separate script, which tests a copy rather than the notebook.
 - Ask before running the full `test_example` pipeline. It deletes and regenerates
   ~300 git-tracked files under `test_example/results/`, and a rerun yields small
   floating-point differences that should not be committed.
 - Describe changes in `CHANGELOG.md` at a high level, including the rationale for
   anything surprising, but leave per-file and per-rule detail to the commit message.
   Bump the version in `pyproject.toml`, which is the single source of truth for it.
 - A version accumulates changes over several pull requests before it is released, so
   its `CHANGELOG.md` heading carries an `(unreleased)` marker while it is in progress.
   Add new entries under that heading rather than starting a new version, and do not
   record the release state anywhere else: not in commit messages, which would go stale
   if the version were renumbered, and not in `pyproject.toml`, whose version is already
   the one being worked toward. Releasing means removing the `(unreleased)` marker and
   tagging the version; the tag is what actually marks a version as released.

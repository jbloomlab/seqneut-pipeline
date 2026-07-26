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
| `plate_to_plate_corr` | `notebooks/plate_to_plate_corr.py` | per plate: correlates its titers with those on the other plates of its group, one panel per other plate |
| `aggregate_titers` | `notebooks/aggregate_titers.py` | combines all sera in a group and builds the interactive titer plots |
| `aggregate_qc_drops` | `notebooks/aggregate_qc_drops.py` | summarizes every QC drop across all plates and sera |
| `build_docs` | `scripts/build_docs.py` | collects the notebook HTML into `results/docs` |
| `miscellaneous_plate_count_barcodes` | `scripts/count_barcodes.py` | counting only, for `miscellaneous_plates` |

`funcs.smk` holds the config-processing helpers that run while the rules are being
defined (`process_plate`, `process_miscellaneous_plates`, `stringify_plate_dates`,
`get_plate_comparators`, `get_collapse_strain_barcodes`, and the sample-table
validation). It relies on `pd` and `re` being imported by `seqneut-pipeline.smk`.

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

`plate_to_plate_corr` is the exception, and shows what justifies one: it runs per plate,
and a project can have many plates that share no serum with any other, so running it for
all of them would put a notebook per plate in the docs with nothing in it.
`get_plate_comparators` therefore derives the plates worth running from the sample tables
in the config, which works only because it is a superset — QC can only remove sera, so a
plate excluded there could never have had anything to show. The notebook still handles the
degenerate case, since QC can empty out a plate that was worth running. Precompute a job
list only when it changes what appears in the docs like this, not merely to skip work.

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
 - A number annotating a plot that must agree with what a selection is showing (such as a
   correlation coefficient) has to be computed by transforms in the layer that draws it,
   with the same `transform_filter` calls ahead of the aggregate, so `vega` recomputes it
   whenever the selection changes. It cannot go in the title: a title is a fixed string
   that can refer to a param but not to a value aggregated from the data. So it is drawn
   as a `mark_text` layer positioned with `x=alt.value(...)`, `y=alt.value(...)`, one
   layer per line of text, since a text mark draws one string and the aggregate of a layer
   cannot be shared with another layer. `notebooks/plate_to_plate_corr.py` builds a shared
   `_base` with the filters and then chains both the points and those layers onto it,
   which is also how to keep the layers from disagreeing (`altair` copies a chart when a
   method is chained onto it, so a shared base is safe).

Since interaction cannot be exercised headlessly, verify it by parsing the emitted spec
out of the exported HTML and checking that each param appears exactly once, at the top
level of the concatenated spec rather than repeated per panel, with a non-empty `views`
list that grows with the number of panels. To get the spec, undo the JavaScript `\uXXXX`
escapes in the HTML, then HTML-unescape the `data-data` attribute of the
`marimo-mime-renderer` whose `data-mime` names `vegalite` *twice*, then `json.loads` it
twice (it is a JSON string holding the JSON spec). Do not expect the `views` entries to
match the panels one for one: `altair` names only some layers, and the point layers of a
layered panel are typically left unnamed, so compare against a chart already known to
work rather than against the panel count.

A value computed by transforms can be checked further than that, because a `vega`
expression is a subset of JavaScript and `node` is available. Pull the transform list and
the dataset for a layer out of the spec, apply the selection by hand (a `{"param": ...}`
filter needs `vega`'s selection machinery, so drop those and filter the rows yourself),
run the remaining `calculate` / `joinaggregate` / `aggregate` transforms in `node` with the
expressions taken verbatim from the spec, and compare the result with the same quantity
computed in `pandas`. Shim only what the expressions use (`format` for a `.Nf` spec is
`toFixed`; `sqrt`, `pow`, `log`, and `LN10` come from `Math`). This catches a bad
expression, which no amount of reading the spec will, and it exercises selection states
that the initial rendering does not.

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
 - the barcode of a strain collapsed by `collapse_strain_barcodes`: the strain name,
   rather than a shared literal such as `collapsed`. `barcode` is a key that identifies a
   strain and not just a label, since the per-barcode QC drops are keyed on
   `(barcode, serum_replicate)` and `(barcode, well)` and the fraction-infectivity chart
   looks up `strain` from `barcode`, so one shared name would make dropping a single
   curve drop every strain.

## Fraction infectivity

Computed in `notebooks/process_plate.py` as:

    frac_infectivity_raw = count / neut_standard_count / median_no_serum_ratio

That is, each viral barcode's count is normalized by its well's total neut-standard
counts, then divided by the median of that same ratio across the plate's no-serum wells.
Under `collapse_strain_barcodes`, `count` is instead the sum over the strain's barcodes,
which are collapsed right after the `manual_drops` are applied so that everything from
the QC onward is per strain.
The `frac_infectivity_ceiling` from `curvefit_params` is then applied as a separate
column; the curve fits use the ceilinged values, while the raw values are also saved.

## Invariants enforced in code

These raise errors, and are not all stated in the README:

 - Groups cannot contain `_` or `|`, as both break wildcards (`seqneut-pipeline.smk`).
   A plate can contain `_`, so a path with both wildcards (as in `plate_to_plate_corr`)
   splits correctly only because the group is constrained to an alternation of the known
   groups. That rule adds a `plate` constraint too. Note that `snakefmt` 2.0.3 wants a
   rule's `wildcard_constraints` block after `log:`; put it anywhere else and `snakefmt`
   silently reorders the rule into something it can no longer parse, moving the docstring
   down with it.
 - Every plate needs at least one `serum: none` sample (`funcs.smk`). No-serum wells are
   designated by `serum: none` with a blank dilution factor or concentration, never by a
   value of 0.
 - All barcodes on a plate, viral library and neut standard set together, must be the
   same length, since `bclen` is inferred from the first one
   (`scripts/count_barcodes.py`).
 - `collapse_strain_barcodes` cannot be combined with `manual_drops` of `barcode_wells` or
   `barcode_serum_replicates` (`funcs.smk`), which remove a barcode from only some wells
   and so would leave a strain's summed counts over a different set of barcodes in those
   wells than in the rest of the plate.
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
   floating-point differences that should not be committed. The run with the barcodes
   collapsed (`--configfile config_collapse_strain_barcodes.yml`, checked against its own
   expected-titers CSV) writes over the same results, so run it in a copy of the
   subdirectory named `_test_example_collapsed`, which `.gitignore` already excludes via
   `_*`.
 - `snakemake --config key=true` yields the *string* `"true"`, not a bool, since
   `snakemake` only evaluates a value capitalized as in Python. A config value that flips
   behavior therefore has to accept both, as `get_collapse_strain_barcodes` does.
 - Describe changes in `CHANGELOG.md` at a high level, including the rationale for
   anything surprising, but leave per-file and per-rule detail to the commit message.
   Bump the version in `pyproject.toml`, which is the single source of truth for it.
 - Keep `CHANGELOG.md` and README entries succinct, and mostly about what a change adds
   or what a result shows rather than how it is implemented. State implementation only
   where it is needed to understand or use the thing, and briefly: for the README, that
   means what a file contains and whether to track it, not how a plot is laid out, and
   for `CHANGELOG.md`, one entry describing what a version ends up with rather than the
   steps taken to get there. Steps that never shipped are not changes and do not belong
   in `CHANGELOG.md` at all, so rewrite an `(unreleased)` entry in place rather than
   adding a second entry describing how the first one was revised.
 - A version accumulates changes over several pull requests before it is released, so
   its `CHANGELOG.md` heading carries an `(unreleased)` marker while it is in progress.
   Add new entries under that heading rather than starting a new version, and do not
   record the release state anywhere else: not in commit messages, which would go stale
   if the version were renumbered, and not in `pyproject.toml`, whose version is already
   the one being worked toward. Releasing means removing the `(unreleased)` marker and
   tagging the version; the tag is what actually marks a version as released.

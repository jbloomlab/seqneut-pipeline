# Instructions for Claude Code

## Bloom lab coding standards

@bloomlab-coding-standards/CLAUDE.md

The [standards](https://github.com/jbloomlab/bloomlab-coding-standards) are pinned at a
commit of that submodule; update periodically with
`git submodule update --remote bloomlab-coding-standards`.

## Additional guidelines for this project

- This pipeline is included as a git submodule by project repositories, which supply the
  `config.yml`, `Snakefile`, `data/`, and `results/` that the lab layout places at a
  repository root; here they exist only under `test_example/`. The rules live in
  `seqneut-pipeline.smk` and `funcs.smk` rather than in `rules/`; their current length is
  fine, but raise splitting them if either passes 1000 lines.
- Maintain backward compatibility, except where breaking it is well justified by a gain in
  clarity, simplicity, or performance; discuss such a change before making it.
- [README](README.md) is the canonical reference for what this pipeline is, how to
  configure it, and what it produces — every config key, QC threshold, and results file.
  Read it and do not duplicate it here.
- A configuration key that can be empty is read as `config.get(key) or {}` rather than with
  a `get` default: emptying such a key from a `--configfile` override requires setting it
  to null, since `snakemake.utils.update_config` merges recursively and so leaves the
  original value untouched given an empty dict.

## Working in this repo

- Verify changes with the format and lint commands in the README, plus
  `cd test_example && snakemake -n`, which catches config and DAG errors cheaply.
- A changed analysis script can be re-run on its own against the committed
  `test_example/results`, without running the whole pipeline: copy `test_example` to a
  `_`-prefixed directory and run `snakemake --force <one of the rule's outputs>
  --rerun-triggers mtime`, which re-runs just that rule from results that are already
  there. Synthesizing input CSVs this way also reaches cases the test example does not
  contain.
- `scripts/seqneut_report.py` and `scripts/seqneut_funcs.py` are imported by the analysis
  scripts and so are tested directly, with `pytest`.
- Ask before running the full `test_example` pipeline: it rewrites ~300 git-tracked files
  under `test_example/results/`, and a rerun yields small floating-point differences that
  should not be committed. The collapsed-barcode run
  (`--configfile config_collapse_strain_barcodes.yml`) overwrites the same results, so run
  it in a copy of the subdirectory named `_test_example_collapsed`, which `.gitignore`
  excludes via `_*`.

## Keeping Altair plots small

The plots here are large enough that the size of the inline JSON matters, so beyond the lab
rules on narrowing a chart's data frame:

- Assert that the key of a `transform_lookup` table is unique, so the join is 1:1 and
  lossless. A looked-up field has no `pandas` dtype for `altair` to infer, so its tooltip
  must name the type explicitly (`alt.Tooltip("strain:N")`).
- A column that exists only to filter rows to a subplot of an `alt.concat` is not needed:
  slice the frame in `pandas` and give each subplot its own. This costs no extra rows,
  since `altair` consolidates every frame into the chart's top-level `datasets` and the
  rows are partitioned among the subplots either way, and it drops both the column and the
  `transform_filter`.
- Round floats that came from arithmetic before plotting. A median over an even number of
  values gives `210.35000000000002`, 18 characters on every row for precision far below a
  pixel.
- A column with the same value on every row (`group` and `serum` in a per-serum notebook,
  or `titer_units`, which follows from `titer_as`) belongs in a `transform_calculate`
  returning a literal, so the value appears once rather than once per row. The quoting is
  `transform_calculate(titer_units=f"'{units}'")`: the argument is a `vega` expression, so
  a string constant needs quotes inside the Python string. As with a lookup, the tooltip
  must then name the type explicitly.
- A param belongs to whichever subplot declares it, so in an `alt.concat` of panels it must
  be declared on the concatenated chart with `.add_params(...)` to be in scope for all of
  them.
- `bind="legend"` cannot drive a selection across panels, as the legend is drawn by one
  panel. Use an input binding (`alt.binding_select`) instead, which the concatenated chart
  can declare; note it renders below the plot, not above it.
- A mouseover selection needs `empty=False`, or every point matches when nothing is
  hovered. `vega-lite` puts this on the predicate rather than the selection, so it appears
  in the emitted spec as `{"param": ..., "empty": false}` in each condition and not in the
  param definition.
- A number annotating a plot that must agree with what a selection shows (a correlation
  coefficient, say) has to be computed by transforms in the layer that draws it, with the
  same `transform_filter` calls ahead of the aggregate, so `vega` recomputes it whenever
  the selection changes. It cannot go in the title, which is a fixed string that can refer
  to a param but not to an aggregate of the data. Draw it instead as a `mark_text` layer
  positioned with `x=alt.value(...)` and `y=alt.value(...)`, one layer per line of text,
  since a text mark draws one string and layers cannot share an aggregate. Chain both the
  points and the text layers onto a shared base holding the filters, which keeps them from
  disagreeing (`altair` copies a chart when a method is chained onto it, so a shared base
  is safe).

Interaction cannot be exercised headlessly, so verify it against the spec parsed out of the
exported HTML: each param should appear exactly once, at the top level of the concatenated
spec rather than repeated per panel, with a non-empty `views` list. Do not expect the
`views` entries to match the panels one for one, as `altair` names only some layers;
compare against a chart already known to work. A value computed by transforms can be
checked further, since a `vega` expression is a subset of JavaScript and `node` is
available: pull a layer's transforms and dataset out of the spec, filter the rows by hand
in place of the `{"param": ...}` predicates, run the remaining expressions verbatim in
`node`, and compare with the same quantity computed in `pandas`. This catches a bad
expression, which no amount of reading the spec will, and exercises selection states that
the initial rendering does not.

## New versions and CHANGELOG

- Describe changes in `CHANGELOG.md` at a high level, including the rationale for anything
  surprising, but leave per-file and per-rule detail to the commit message. Say what a
  version ends up with rather than the steps taken to get there: rewrite an `(unreleased)`
  entry in place rather than adding a second describing how the first was revised. Bump the
  version in `pyproject.toml`, the single source of truth for it.
- A version accumulates changes over several pull requests, so its `CHANGELOG.md` heading
  carries an `(unreleased)` marker while in progress; add new entries under that heading
  rather than starting a new version. Do not record the release state anywhere else.
  Releasing means removing the marker and tagging the version, and the tag is what actually
  marks a version as released.

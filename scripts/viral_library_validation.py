"""Validation of a viral library CSV against the checks configured for it.

The checks are declared per column in `viral_library_validations`, so that only `strain`
and `barcode` are named here: they are the sole columns `seqneut-pipeline` itself reads,
and everything else a library carries is described entirely by the configuration.

Nothing here reads or writes a file beyond the library CSV, and nothing prints, so that
every check can be exercised directly by the tests.

"""

import collections
import re

import Bio.Data.CodonTable
import Bio.Seq
import pandas as pd

# the only columns `seqneut-pipeline` itself reads from a viral library
REQUIRED_COLS = {"strain", "barcode"}

# the keys a `columns` entry may have, the subset a `groups` entry may have, and the
# values `unique` may take
SPEC_KEYS = {
    "nulls",
    "unique",
    "unique_within",
    "values",
    "non_null_when",
    "matches",
    "length",
    "translates_to",
    "numeric",
    "by",
    "groups",
}
GROUP_KEYS = {"matches", "length"}
UNIQUE_VALUES = {"barcode", "strain"}

# how many failing values a failed check names before it just gives the count
MAX_OFFENDERS_SHOWN = 10

# the sections of the report, in the order it presents them
SECTIONS = ["columns", "numbers", "patterns", "uniqueness", "sequences"]

# a check of a library, in the section of the report it belongs to. `offenders` empty
# means the check passed; `id` names the check stably, as `description` interpolates
# configured values and so is not something to match on.
Check = collections.namedtuple(
    "Check", ["id", "section", "description", "examined", "offenders"]
)


class SpecError(ValueError):
    """The `viral_library_validations` entry for a library is itself malformed."""


class LibraryError(ValueError):
    """A viral library failed the checks configured for it.

    `failures` holds the failed `Check`s and `report` the whole report, so that a failing
    run still says everything the library got through.

    """

    def __init__(self, message, failures, report):
        super().__init__(message)
        self.failures = failures
        self.report = report


def _require(condition, message):
    """Raise `SpecError` describing a malformed spec unless `condition` holds."""
    if not condition:
        raise SpecError(f"`viral_library_validations` is invalid: {message}")


def _validate_patterns(spec, where):
    """Validate the `matches` and `length` of a column or group spec, described by `where`."""
    if "matches" in spec:
        _require(
            isinstance(spec["matches"], list) and spec["matches"],
            f"`{where}: matches` must be a non-empty list of regexes",
        )
        for pattern in spec["matches"]:
            _require(
                isinstance(pattern, str),
                f"`{where}: matches` holds {pattern!r}, which is not a regex",
            )
            try:
                re.compile(pattern)
            except re.error as e:
                _require(
                    False, f"`{where}: matches` holds invalid regex {pattern!r}: {e}"
                )
    if "length" in spec:
        length = spec["length"]
        _require(
            isinstance(length, dict) and set(length) == {"min", "max"},
            f"`{where}: length` must have exactly `min` and `max`, but has "
            f"{sorted(length) if isinstance(length, dict) else length!r}",
        )
        _require(
            all(isinstance(length[k], int) for k in ("min", "max"))
            and length["min"] <= length["max"],
            f"`{where}: length` must be integers with `min` no greater than `max`, "
            f"but is {length}",
        )


def _require_non_null_grouping(col_specs, grouping_col, where):
    """Require `grouping_col` to be non-null, as it partitions the rows of a check.

    A row whose grouping column is empty belongs to no group, and so would be quietly
    dropped from the check the grouping keys rather than failing it.

    """
    _require(
        not col_specs[grouping_col]["nulls"],
        f"`{where}` groups on {grouping_col!r}, which must be declared `nulls: false`: a "
        "row that does not set it belongs to no group and so would go unchecked",
    )


def validate_spec(spec):
    """Raise `SpecError` describing `spec` unless it is valid checks for one library.

    Called before the library is read, so that a configuration error is never reported as
    though it were a problem with the data.

    """
    _require(
        isinstance(spec, dict) and set(spec) == {"columns"},
        "each library's entry must have exactly the key `columns`, but has "
        f"{sorted(spec) if isinstance(spec, dict) else spec!r}",
    )
    col_specs = spec["columns"]
    _require(
        isinstance(col_specs, dict) and col_specs,
        "`columns` must be a non-empty mapping of column name to its checks",
    )
    _require(
        REQUIRED_COLS <= set(col_specs),
        f"`columns` must list {sorted(REQUIRED_COLS)}, the columns `seqneut-pipeline` "
        f"itself reads, but is missing {sorted(REQUIRED_COLS - set(col_specs))}",
    )

    for col, col_spec in col_specs.items():
        where = f"columns: {col}"
        _require(
            isinstance(col_spec, dict),
            f"`{where}` must be a mapping of check name to its value",
        )
        _require(
            not set(col_spec) - SPEC_KEYS,
            f"`{where}` has unknown keys {sorted(set(col_spec) - SPEC_KEYS)}, "
            f"expected some of {sorted(SPEC_KEYS)}",
        )
        _require(
            isinstance(col_spec.get("nulls"), bool),
            f"`{where}` must declare `nulls` as true or false",
        )

        if "unique" in col_spec:
            _require(
                col_spec["unique"] in UNIQUE_VALUES,
                f"`{where}` has `unique: {col_spec['unique']}`, expected one of "
                f"{sorted(UNIQUE_VALUES)}",
            )
            # a strain names itself, so the check would hold no matter what the data said
            _require(
                not (col == "strain" and col_spec["unique"] == "strain"),
                "`columns: strain` cannot have `unique: strain`, which asks whether a "
                "strain names one strain and so is true of any library",
            )
        if "unique_within" in col_spec:
            _require(
                col_spec.get("unique") == "strain",
                f"`{where}` has `unique_within` without `unique: strain`, which is the "
                "only check it scopes",
            )
            _require(
                col_spec["unique_within"] in col_specs
                and col_spec["unique_within"] != col,
                f"`{where}` has `unique_within: {col_spec['unique_within']}`, which is "
                "not another column listed in `columns`",
            )
            _require_non_null_grouping(col_specs, col_spec["unique_within"], where)

        if "values" in col_spec:
            _require(
                isinstance(col_spec["values"], list) and col_spec["values"],
                f"`{where}: values` must be a non-empty list",
            )
            # the CSV is read as text, so an unquoted YAML boolean or number would be
            # compared as a Python object and match nothing in the library
            _require(
                all(isinstance(v, str) for v in col_spec["values"]),
                f"`{where}: values` must all be strings, but holds "
                f"{[v for v in col_spec['values'] if not isinstance(v, str)]}; quote "
                "them, as the library is read as text",
            )

        if "non_null_when" in col_spec:
            when = col_spec["non_null_when"]
            _require(
                isinstance(when, dict) and len(when) == 1,
                f"`{where}: non_null_when` must name exactly one column, but names "
                f"{sorted(when) if isinstance(when, dict) else when!r}",
            )
            ((when_col, when_val),) = when.items()
            _require(
                when_col in col_specs,
                f"`{where}: non_null_when` names {when_col!r}, which `columns` does not "
                "list",
            )
            when_values = col_specs[when_col].get("values")
            _require(
                when_values is None or when_val in when_values,
                f"`{where}: non_null_when` expects {when_col} == {when_val!r}, which is "
                f"not among the {when_values} that column may hold",
            )

        if "translates_to" in col_spec:
            _require(
                col_spec["translates_to"] in col_specs
                and col_spec["translates_to"] != col,
                f"`{where}: translates_to` names {col_spec['translates_to']!r}, which "
                "is not another column listed in `columns`",
            )
        if "numeric" in col_spec:
            _require(
                isinstance(col_spec["numeric"], bool),
                f"`{where}: numeric` must be true or false",
            )
        _validate_patterns(col_spec, where)

        if "by" in col_spec or "groups" in col_spec:
            _require(
                "by" in col_spec and "groups" in col_spec,
                f"`{where}` must have `by` and `groups` together, as `by` names the "
                "column whose values key `groups`",
            )
            by_col = col_spec["by"]
            _require(
                by_col in col_specs and by_col != col,
                f"`{where}: by` names {by_col!r}, which is not another column listed in "
                "`columns`",
            )
            by_values = col_specs[by_col].get("values")
            _require(
                by_values is not None,
                f"`{where}: by` names {by_col!r}, which must declare the `values` that "
                "key `groups`",
            )
            _require_non_null_grouping(col_specs, by_col, where)
            _require(
                set(col_spec["groups"]) == set(by_values),
                f"`{where}: groups` is keyed by {sorted(col_spec['groups'])}, but "
                f"`{by_col}` holds {sorted(by_values)}; every value needs an entry and "
                "an entry for anything else would never be applied",
            )
            for group, group_spec in col_spec["groups"].items():
                group_where = f"{where}: groups: {group}"
                _require(
                    isinstance(group_spec, dict) and group_spec,
                    f"`{group_where}` must be a non-empty mapping",
                )
                _require(
                    not set(group_spec) - GROUP_KEYS,
                    f"`{group_where}` has unknown keys "
                    f"{sorted(set(group_spec) - GROUP_KEYS)}, expected some of "
                    f"{sorted(GROUP_KEYS)}",
                )
                _validate_patterns(group_spec, group_where)


def read_library(csv):
    """Read a viral library CSV, every column as text so no identifier becomes a number."""
    return pd.read_csv(csv, dtype=str)


def describe_column(col, spec):
    """Human-readable summary of the checks `spec` configures on column `col`."""
    col_spec = spec["columns"].get(col)
    if col_spec is None:
        return "NOT VALIDATED"
    parts = ["nulls allowed" if col_spec["nulls"] else "non-null"]
    if col_spec.get("unique") == "barcode":
        parts.append("distinct per row")
    elif col_spec.get("unique") == "strain":
        within = col_spec.get("unique_within")
        parts.append(
            f"names one strain within its {within}" if within else "names one strain"
        )
    if "values" in col_spec:
        parts.append(f"values {col_spec['values']}")
    if "non_null_when" in col_spec:
        ((when_col, when_val),) = col_spec["non_null_when"].items()
        parts.append(f"set when {when_col} is {when_val}")
    for pattern in col_spec.get("matches", []):
        parts.append(f"matches {pattern}")
    if "length" in col_spec:
        parts.append(f"length {col_spec['length']['min']}-{col_spec['length']['max']}")
    if "by" in col_spec:
        parts.append(f"further checks per {col_spec['by']}")
    if "translates_to" in col_spec:
        parts.append(f"translates to {col_spec['translates_to']}")
    if col_spec.get("numeric"):
        parts.append("a number")
    return "; ".join(parts)


def _table_lines(header, rows, last_free=True):
    """Lines of a table, the first column left-aligned and the rest right-aligned.

    The last column is left as written where `last_free`, as it holds free text too wide
    to pad every other row out to.

    """
    n = len(header) - 1 if last_free else len(header)
    widths = [max(len(row[i]) for row in [header, *rows]) for i in range(n)]
    return [
        "  "
        + row[0].ljust(widths[0])
        + "".join(f"  {v.rjust(w)}" for v, w in zip(row[1:n], widths[1:]))
        + (f"  {row[-1]}" if last_free else "")
        for row in [header, *rows]
    ]


def summary_lines(df, spec, csv):
    """Lines describing the library, its categorical columns, and every column in it."""
    col_specs = spec["columns"]
    lines = [
        f"=== summary of {csv} ===",
        f"{len(df)} barcodes for {df['strain'].nunique()} strains",
    ]

    # a breakdown per value of every column declared categorical, so that the summary
    # names no particular column of its own
    strains = df.drop_duplicates("strain")
    for col, col_spec in col_specs.items():
        if "values" not in col_spec:
            continue
        lines.append("")
        lines.extend(
            _table_lines(
                (col, "strains", "barcodes"),
                [
                    (
                        f"  {value}",
                        str(int((strains[col] == value).sum())),
                        str(int((df[col] == value).sum())),
                    )
                    for value in col_spec["values"]
                ],
                last_free=False,
            )
        )

    # one row per column of the CSV, in the order the CSV has them, so that a column the
    # configuration does not mention still shows up rather than passing unnoticed
    lines.append("")
    lines.extend(
        _table_lines(
            ("column", "null", "non-null", "distinct", "max per strain", "checks"),
            [
                (
                    col,
                    str(int(df[col].isnull().sum())),
                    str(int(df[col].notnull().sum())),
                    str(int(df[col].nunique())),  # counts distinct non-null values
                    str(int(df.groupby("strain")[col].nunique(dropna=False).max())),
                    describe_column(col, spec),
                )
                for col in df.columns
            ],
        )
    )

    unlisted = [col for col in df.columns if col not in col_specs]
    lines.append(
        "\ncolumns not listed in `viral_library_validations: columns`: "
        + (str(unlisted) if unlisted else "none")
    )
    return lines


def _offenders(df, failed, col):
    """The strains and values of the rows of `df` where the boolean `failed` holds."""
    rows = df.loc[failed]
    return [f"{strain}: {value!r}" for strain, value in zip(rows["strain"], rows[col])]


def schema_checks(df, spec):
    """Checks that the library has the columns the spec lists."""
    col_specs = spec["columns"]
    yield Check(
        id="schema",
        section="schema",
        description="all columns listed in `columns` are present",
        examined=f"{len(col_specs)} columns",
        offenders=[col for col in col_specs if col not in df.columns],
    )


def column_checks(df, spec):
    """Checks of `nulls`, `values`, `non_null_when`, and `numeric`."""
    col_specs = spec["columns"]
    non_null_cols = [col for col, s in col_specs.items() if not s["nulls"]]
    yield Check(
        id="nulls",
        section="columns",
        description="columns declared non-null have no null values",
        examined=f"{len(non_null_cols)} columns",
        offenders=[
            f"{col}: {int(df[col].isnull().sum())} null"
            for col in non_null_cols
            if df[col].isnull().any()
        ],
    )
    for col, col_spec in col_specs.items():
        set_rows = df[df[col].notnull()]
        if "values" in col_spec:
            yield Check(
                id=f"values:{col}",
                section="columns",
                description=f"`{col}` holds only {col_spec['values']}",
                examined=f"{len(set_rows)} non-null rows",
                offenders=sorted(
                    set(set_rows.loc[~set_rows[col].isin(col_spec["values"]), col])
                ),
            )
        if "non_null_when" in col_spec:
            ((when_col, when_val),) = col_spec["non_null_when"].items()
            yield Check(
                id=f"non_null_when:{col}",
                section="columns",
                description=(
                    f"`{col}` is set for exactly the rows with {when_col} "
                    f"== {when_val}"
                ),
                examined=f"{len(df)} rows",
                offenders=sorted(
                    set(
                        df.loc[
                            df[col].notnull() != (df[when_col] == when_val), "strain"
                        ]
                    )
                ),
            )
        if col_spec.get("numeric"):
            yield Check(
                id=f"numeric:{col}",
                section="numbers",
                description=f"`{col}` is a number",
                examined=f"{len(set_rows)} non-null rows",
                offenders=_offenders(
                    set_rows, pd.to_numeric(set_rows[col], errors="coerce").isna(), col
                ),
            )


def _pattern_checks_over(df, col, patterns, length, id_suffix, described):
    """Checks of `patterns` and `length` over the rows of `df` that set `col`."""
    set_rows = df[df[col].notnull()]
    examined = f"{len(set_rows)} non-null rows"
    for pattern in patterns:
        yield Check(
            id=f"matches:{col}:{id_suffix}{pattern}",
            section="patterns",
            description=f"`{col}`{described} matches {pattern}",
            examined=examined,
            offenders=_offenders(
                set_rows,
                ~set_rows[col].map(re.compile(pattern).search).astype(bool),
                col,
            ),
        )
    if length is not None:
        lengths = set_rows[col].str.len()
        yield Check(
            id=f"length:{col}:{id_suffix}".rstrip(":"),
            section="patterns",
            description=(
                f"`{col}`{described} is {length['min']} to {length['max']} characters"
                if length["min"] != length["max"]
                else f"`{col}`{described} is {length['min']} characters"
            ),
            examined=examined,
            offenders=_offenders(
                set_rows, (lengths < length["min"]) | (lengths > length["max"]), col
            ),
        )


def pattern_checks(df, spec):
    """Checks of `matches` and `length`, including those keyed by `by` and `groups`."""
    for col, col_spec in spec["columns"].items():
        yield from _pattern_checks_over(
            df, col, col_spec.get("matches", []), col_spec.get("length"), "", ""
        )
        for group, group_spec in col_spec.get("groups", {}).items():
            yield from _pattern_checks_over(
                df[df[col_spec["by"]] == group],
                col,
                group_spec.get("matches", []),
                group_spec.get("length"),
                f"{group}:",
                f" where {col_spec['by']} is {group}",
            )


def uniqueness_checks(df, spec):
    """Checks of `unique`, `unique_within`, and that a strain's barcodes agree."""
    col_specs = spec["columns"]
    for col, col_spec in col_specs.items():
        # nulls are excluded, so a column that is optional per barcode can still be
        # required to identify the barcodes that do set it
        set_rows = df[df[col].notnull()]
        if col_spec.get("unique") == "barcode":
            yield Check(
                id=f"unique_barcode:{col}",
                section="uniqueness",
                description=f"`{col}` is distinct in every row that sets it",
                examined=f"{len(set_rows)} non-null rows",
                offenders=sorted(set(set_rows.loc[set_rows[col].duplicated(), col])),
            )
        elif col_spec.get("unique") == "strain":
            within = col_spec.get("unique_within")
            keys = ([within] if within else []) + [col]
            per_value = set_rows.drop_duplicates([*keys, "strain"]).groupby(keys)[
                "strain"
            ]
            yield Check(
                id=f"unique_strain:{col}",
                section="uniqueness",
                description=(
                    f"each `{col}` value names one strain within its {within}"
                    if within
                    else f"each `{col}` value names one strain"
                ),
                examined=f"{len(per_value)} distinct values",
                offenders=per_value.apply(list)[per_value.nunique() > 1].tolist(),
            )

    # every column but a `unique: barcode` one describes the strain, so it must hold the
    # same value on each of that strain's barcodes
    # `strain` groups this check, so it stays even where it is declared per barcode --
    # a library of one barcode per strain may legitimately say so
    per_barcode_cols = [
        col
        for col, s in col_specs.items()
        if s.get("unique") == "barcode" and col != "strain"
    ]
    per_strain = df[list(col_specs)].drop(columns=per_barcode_cols).drop_duplicates()
    yield Check(
        id="constant_per_strain",
        section="uniqueness",
        description=(
            f"each strain has one set of values for all listed columns but "
            f"{per_barcode_cols}"
        ),
        examined=f"{df['strain'].nunique()} strains",
        offenders=sorted(
            set(per_strain.loc[per_strain["strain"].duplicated(keep=False), "strain"])
        ),
    )


def _translation(nt):
    """The protein `nt` translates to, or `None` where it is not a coding sequence.

    What the sequence column may actually hold is the `matches` regex's business, so this
    tries to translate anything of a whole number of codons and reports what will not
    translate as a failure of `translates_to` rather than as a crash. Duplicating the
    format check here would instead make one malformed sequence fail both.

    """
    if len(nt) % 3:  # a partial codon, which `Bio` only warns about
        return None
    try:
        return str(Bio.Seq.Seq(nt).translate())
    except Bio.Data.CodonTable.TranslationError:
        return None


def sequence_checks(df, spec):
    """Checks of `translates_to`."""
    for col, col_spec in spec["columns"].items():
        if "translates_to" not in col_spec:
            continue
        protein_col = col_spec["translates_to"]
        rows = df[df[col].notnull() & df[protein_col].notnull()]
        yield Check(
            id=f"translates_to:{col}",
            section="sequences",
            description=f"`{col}` translates to `{protein_col}`",
            examined=f"{len(rows)} rows setting both",
            offenders=sorted(
                {
                    strain
                    for strain, nt, protein in zip(
                        rows["strain"], rows[col], rows[protein_col]
                    )
                    if _translation(nt) != protein
                }
            ),
        )


def library_checks(df, spec):
    """Every check `spec` configures on library `df`, as `Check`s in report order."""
    yield from column_checks(df, spec)
    yield from pattern_checks(df, spec)
    yield from uniqueness_checks(df, spec)
    yield from sequence_checks(df, spec)


def numeric_lines(df, spec):
    """Lines reporting the range of each column declared `numeric`."""
    lines = []
    for col, col_spec in spec["columns"].items():
        if not col_spec.get("numeric"):
            continue
        # one value per strain, so a strain carrying several barcodes is not weighted by
        # them; two decimals is the precision a decimal-year date is recorded to
        values = pd.to_numeric(
            df.drop_duplicates("strain")[col], errors="coerce"
        ).dropna()
        lines.append(
            f"      `{col}` over {len(values)} strains:  min {values.min():.2f}"
            f"  median {values.median():.2f}  max {values.max():.2f}"
            if len(values)
            else f"      `{col}`: no numeric values over {len(values)} strains"
        )
    return lines


def _check_lines(checks, width):
    """Lines reporting `checks`, each naming what it examined and any failures."""
    lines = []
    for c in checks:
        status = "FAIL" if c.offenders else " OK "
        lines.append(f"  {status}  {c.description.ljust(width)}  ({c.examined})")
        if c.offenders:
            shown = c.offenders[:MAX_OFFENDERS_SHOWN]
            more = len(c.offenders) - len(shown)
            lines.append(
                f"        {len(c.offenders)} failures: {shown}"
                + (f" and {more} more" if more else "")
            )
    return lines


def report(csv, spec):
    """Validate the viral library `csv` against `spec`.

    Returns the report as text. Raises `SpecError` if `spec` is malformed, and
    `LibraryError` carrying the whole report if any check fails.

    """
    validate_spec(spec)
    df = read_library(csv)

    schema = list(schema_checks(df, spec))
    # nothing else can be checked against columns that are not there, so a failed schema
    # ends the report rather than raising a `KeyError` from every check after it
    if any(c.offenders for c in schema):
        text = "\n".join(
            ["=== schema ===", *_check_lines(schema, len(schema[0].description))]
        )
        raise LibraryError(f"{csv} is missing listed columns", schema, text)

    checks = list(library_checks(df, spec))
    width = max(len(c.description) for c in [*schema, *checks])
    lines = ["=== schema ===", *_check_lines(schema, width), ""]
    lines.extend(summary_lines(df, spec, csv))
    # a fixed order, so the report does not depend on the order the checks are yielded in
    for section in SECTIONS:
        in_section = [c for c in checks if c.section == section]
        # the range of a numeric column is reported alongside the check that it is one
        extra = numeric_lines(df, spec) if section == "numbers" else []
        if in_section or extra:
            lines.extend(["", f"=== {section} ===", *_check_lines(in_section, width)])
            lines.extend(extra)
    text = "\n".join(lines)

    failures = [c for c in checks if c.offenders]
    if failures:
        raise LibraryError(
            f"{len(failures)} check{'s' if len(failures) > 1 else ''} failed for "
            f"{csv}: {[c.description for c in failures]}",
            failures,
            text,
        )
    return text

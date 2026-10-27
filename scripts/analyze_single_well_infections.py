"""Analyze the plates that infect each well with a single virus of the library.

Writes a self-contained HTML report, a CSV of what every well held, and a CSV of how each
well's reads split across the barcodes it was expected to hold. The question these plates
answer is whether each well holds only what it should, and whether virus from one well ends
up in another; they also give a loose estimate of how well each virus grows on its own.

What a well should hold is read from the plates' samples CSVs. Where they have a `barcode`
column each well was infected with one barcoded variant and should hold that barcode alone
(with any other barcode carried in the same tube, where the viral library has a `stock_id`
column); otherwise their `strain` column names a strain and the well should hold every
barcode of it. That is the well's *unit*, and its *expected barcodes* are what belongs in
it, along with the neutralization standard spiked in. Every parsed read is put into one of
five categories that partition it:

    expected            one of the well's expected barcodes
    neut standard       a barcode of the neutralization standard
    other library       a library barcode that is not expected here, so it does not belong
    unmatched, near     matched neither whitelist, but is within `MAX_MISREAD_HAMMING`
                        of an expected barcode or of the neutralization standard, so it
                        is the well's own material misread
    unmatched, far      matched neither whitelist and is further than that, so it cannot be
                        attributed to what belongs in the well

Before that, the reads of any `discarded_barcodes`, and those within `MAX_MISREAD_HAMMING`
of one, are set aside: these are barcodes physically in every well that no analysis should
count. They are reported but enter no category or fraction.

When each well holds one barcode, another barcode of the well's own strain is *not*
expected: the well was supposed to hold one variant, so a sibling barcode is material that
does not belong, exactly like any other.

The titer is the expected reads against the neutralization standard in the well, `expected /
neut standard`. Only those two are in that ratio, so a contaminated well's titer is not
inflated by whatever else was in the well; that is reported separately instead.

Because `well` is the physical position on the plate, contamination can be traced rather
than only counted: a read that belongs to another well is attributed to that well and to how
far away it is, counting a step onto a diagonal neighbor as one. Material from a neighboring
well is spillover; material from far away, from a well on another plate of the group, or of
a unit that no well here was infected with, arrived some other way. This works because each
unit has at most one home well per plate, which is checked below.

Two things to hold on to when reading the report:

  - A well whose own virus barely grew shows a large foreign *fraction* from a small
    absolute spill, so every fraction is reported beside the counts behind it.
  - When a well holds a whole strain, its reads are summed over the strain's barcodes, of
    which strains have differing numbers, so a strain with more barcodes contributes more
    reads at the same titer. That count is an upper bound on the inflation rather than a
    correction for it, since a barcode the virus does not carry contributes nothing; the
    last section reports which expected barcodes produced no reads at all.

One group of plates is analyzed together and reported in one HTML. Their barcodes are
counted against one whitelist, the group's `viral_library`. Where `final_viral_library` is
set, the units whose barcodes were all dropped from it are marked wherever they appear; it
is never used in a calculation. Everything specific to a group comes from the
`analyze_single_well_infections` section of the configuration via `snakemake.params`, as
does the top-level `discarded_barcodes`.

"""

import sys
from pathlib import Path

import altair as alt
import numpy as np
import pandas as pd
from miscellaneous_plate_charts import (
    CARRIED_FORWARD_LABELS,
    DROPPED_COLOR,
    LABEL_LIMIT,
    POINT_STROKE_WIDTH,
    POINT_STROKE_WIDTH_HOVERED,
    STROKE_WIDTH,
    STROKE_WIDTH_HOVERED,
    WellPanels,
    carried_forward_stroke,
    lookup_tooltip,
    stroke_width,
    well_hover,
)
from miscellaneous_plate_funcs import (
    ANOTHER_PLATE,
    FATE_EXPLANATIONS,
    MAX_MISREAD_HAMMING,
    MAX_NAMED_DISTANCE,
    NOT_ASSAYED,
    PlateReport,
    all_origins,
    check_discarded_barcodes,
    check_fates_add_up,
    distance_origin,
    min_hamming,
    read_per_well_csvs,
    reads_per_well,
    split_discarded,
    trace_origin,
    well_position,
)
from seqneut_funcs import narrow_for_altair

# `noqa: SIM115` as this log file must stay open for the life of the script
sys.stderr = sys.stdout = open(snakemake.log[0], "w")  # noqa: SIM115

_ = alt.data_transformers.disable_max_rows()

# Other units and unmatched barcodes listed for each well the two reporting thresholds pick
# out. A display limit, affecting only how much of a well is broken out and never what is
# computed, so it lives here rather than in the configuration.
N_TOP_BARCODES = 3

# the non-neut-standard parsed reads the categories above are taken as fractions of, spelled
# with hyphens as it is a column of the output CSV as well as of the frame here
NON_NEUT_STANDARD_READS = "non-neut-standard_parsed_reads"

# What an empty well is called on the charts' axis. A well whose unit is blank was left
# uninfected, so it holds only the neutralization standard spiked into it and every library
# barcode in it is material that does not belong.
EMPTY_WELL = "empty"


# --- read the configuration -------------------------------------------------------------

single_well = snakemake.wildcards.single_well
single_well_config = snakemake.params.single_well_config
plates = single_well_config["miscellaneous_plates"]
addtl_viral_cols = single_well_config["addtl_viral_cols_to_show"]
has_final_library = single_well_config["final_viral_library"] is not None

if not isinstance(single_well_config["description"], str):
    raise TypeError(f"'description' of {single_well} must be text")
if not (
    isinstance(addtl_viral_cols, list)
    and all(isinstance(col, str) for col in addtl_viral_cols)
):
    raise TypeError(
        f"'addtl_viral_cols_to_show' of {single_well} must be a list of column names, "
        f"not {addtl_viral_cols!r}"
    )

# The thresholds a well is flagged and warned about by. None of them changes what is
# computed: `max_other_library_frac` and `max_unmatched_frac` pick out the wells whose
# foreign and unmatched material are broken out, and `min_parsed_reads_per_well` and
# `min_neut_standard_reads_for_titer` the wells too shallow for their fractions, or with too
# little of the standard for their titer, to rest on much.
_thresholds = {}
for _key in [
    "max_other_library_frac",
    "max_unmatched_frac",
    "min_parsed_reads_per_well",
    "min_neut_standard_reads_for_titer",
]:
    _value = single_well_config[_key]
    if isinstance(_value, bool) or not isinstance(_value, (int, float)) or _value < 0:
        raise ValueError(
            f"{_key!r} of {single_well} must be a non-negative number, not {_value!r}"
        )
    _thresholds[_key] = _value
max_other_library_frac = _thresholds["max_other_library_frac"]
max_unmatched_frac = _thresholds["max_unmatched_frac"]
min_parsed_reads_per_well = _thresholds["min_parsed_reads_per_well"]
min_neut_standard_reads_for_titer = _thresholds["min_neut_standard_reads_for_titer"]

# The four categories that partition the virus in a well: the column each is computed in,
# mapped to the label it carries in the report, in the order they stack.
#
# The neutralization standard is not among them, and every fraction in this report is a
# fraction of the *non-neut-standard* parsed reads, named `NON_NEUT_STANDARD_READS` below.
# The standard is a fixed spike rather than part of what a well grew, so including it would
# make each of these fractions depend on how much virus the well happened to hold.
COMPOSITION_CATEGORIES = {
    "expected_reads": "expected",
    "other_library_reads": "other library",
    "unmatched_near_expected_reads": f"≤{MAX_MISREAD_HAMMING} nt from expected/standard",
    "unmatched_far_expected_reads": f">{MAX_MISREAD_HAMMING} nt from expected/standard",
}


# the axis the wells run along in every chart, and so how the panels of each are laid out
panels = WellPanels(single_well_config["strain_axis"])

print(f"Analyzing single virus per well infections {single_well} over plates {plates}")

report = PlateReport(title=f"Single virus per well infections: {single_well}")


report.md(f"""
    This document analyzes the single virus per well infection experiment -
    {single_well}. The purpose of this analysis is a quality control check on the
    individual viral stocks to validate that each stock contains the expected barcode(s)
    and to assess the level of barcode contamination. Additionally, it provides estimates
    of the titer of each viral stock. These titers should be interpreted carefully as they
    assume each infection was performed in the linear range (where there is a linear
    relationship between the amount of input virus and resulting spike-in normalized
    counts), which may not always be true for every virus in these experiments. The plots
    below are interactive: mouse over for details.

    ## Overview
    """)
report.md(single_well_config["description"])  # added separately as it is multi-line

# --- read the samples and the libraries ------------------------------------------------

samples_by_plate = {
    plate: pd.read_csv(csv_file, dtype={"well": str, "strain": str, "barcode": str})
    for (plate, csv_file) in zip(plates, snakemake.input.samples_csvs, strict=True)
}

# What each well was infected with, read off the samples CSVs: a `barcode` column means one
# barcoded variant to a well, and otherwise the `strain` column names a whole strain. Every
# plate of a group has to agree, as the report puts all of their wells on one axis.
_holds = {
    plate: "barcode" if "barcode" in df.columns else "strain"
    for (plate, df) in samples_by_plate.items()
}
if len(set(_holds.values())) != 1:
    raise ValueError(
        f"the plates of {single_well} do not agree on whether a well holds a strain or a "
        f"barcode: {_holds}. A samples CSV with a 'barcode' column names a barcode per "
        "well, and one without it a strain per well"
    )
well_holds = next(iter(_holds.values()))

# A dilution factor is optional, but if any plate gives it then every plate must, as the
# titers of all of them are ranked against each other below.
_with_dilution = {
    plate for (plate, df) in samples_by_plate.items() if "dilution_factor" in df
}
if _with_dilution and _with_dilution != set(plates):
    raise ValueError(
        f"only some of the plates of {single_well} give a 'dilution_factor': "
        f"{sorted(_with_dilution)}"
    )
has_dilution = bool(_with_dilution)

for plate, df in samples_by_plate.items():
    missing = {"well", well_holds} - set(df.columns)
    if missing:
        raise ValueError(
            f"the samples CSV of {plate} lacks the columns {sorted(missing)}"
        )
    # A well is uninfected only where nothing names what infected it, so a strain beside a
    # blank barcode is a barcode left out rather than an empty well.
    if well_holds == "barcode" and "strain" in df.columns:
        _no_barcode = df[df["barcode"].isna() & df["strain"].notna()]
        if len(_no_barcode):
            raise ValueError(
                f"the samples CSV of {plate} names a strain but no barcode for the wells "
                f"{_no_barcode['well'].tolist()}, so what infected them is not known"
            )

samples = pd.concat(
    [
        df[["well", well_holds, *(["dilution_factor"] if has_dilution else [])]].assign(
            plate=plate
        )
        for (plate, df) in samples_by_plate.items()
    ],
    ignore_index=True,
).rename(columns={well_holds: "expected"})

# A well that names no unit was left uninfected. A blank is the record here rather than a
# value that has to be interpreted.
samples["is_empty"] = samples["expected"].isna()
samples[["row", "column"]] = [well_position(well) for well in samples["well"]]

if samples.duplicated(["plate", "well"]).any():
    raise ValueError("a plate has two samples in the same well")
if samples[~samples["is_empty"]].duplicated(["plate", "expected"]).any():
    raise ValueError(
        f"a plate has the same {well_holds} in more than one well, so a read of that "
        f"{well_holds} cannot be traced back to the well it came from"
    )

# Where each well falls in the order the plate was filled, down each column and then across.
# The charts order their rows by what they are about, so this is for the output CSV.
samples = samples.sort_values(["plate", "column", "row"]).reset_index(drop=True)
samples["position"] = samples.groupby("plate").cumcount()

if has_dilution:
    # over the infected wells only, as an empty well has no titer to compare
    dilutions = sorted(samples.loc[~samples["is_empty"], "dilution_factor"].unique())
    # Every titer below is a ratio to the same fixed spike, so two wells are comparable only
    # if they were at the same dilution, and they are all ranked against each other.
    if len(dilutions) != 1:
        raise ValueError(
            f"the wells of {single_well} are at more than one dilution factor "
            f"{dilutions}, so their titers are not comparable as computed"
        )

# identifiers read as strings, so that one made of digits is not taken for a number
viral_library = pd.read_csv(
    snakemake.input.viral_library,
    dtype={"barcode": str, "strain": str, "stock_id": str},
)
# A stock is one tube, which can carry more than one barcode. Where the library does not
# record stocks, each barcode is its own.
if "stock_id" not in viral_library.columns:
    viral_library["stock_id"] = viral_library["barcode"]
viral_library["stock_id"] = viral_library["stock_id"].fillna(viral_library["barcode"])
# A tube holds virus of one strain, and a well infected from it is labeled by that strain.
_strains_of_stock = viral_library.groupby("stock_id")["strain"].unique()
_mixed_stocks = _strains_of_stock[_strains_of_stock.map(len) > 1]
if len(_mixed_stocks):
    raise ValueError(
        f"these stocks of {snakemake.input.viral_library} carry barcodes of more than one "
        "strain, but a stock is one tube of one strain: "
        + ", ".join(
            f"{stock} ({', '.join(sorted(strains))})"
            for (stock, strains) in _mixed_stocks.items()
        )
    )
neut_standard_barcodes = pd.read_csv(
    snakemake.input.neut_standard_set, dtype={"barcode": str}
)["barcode"].tolist()
barcode_class = pd.concat(
    [
        viral_library[["barcode", "strain"]].assign(neut_standard=False),
        pd.DataFrame({"barcode": neut_standard_barcodes}).assign(
            neut_standard=True, strain=pd.NA
        ),
    ],
    ignore_index=True,
)

# The units of the library, which are what a well can be infected with, and the barcodes
# each brings to its well. A well holds a stock and a stock brings every barcode it carries,
# so this resolves through the stock: a passenger barcode sharing a tube belongs to the same
# well rather than to none.
barcodes_of_stock = viral_library.groupby("stock_id")["barcode"].apply(set).to_dict()
unit_stocks = viral_library.groupby(well_holds)["stock_id"].apply(set).to_dict()
unit_barcodes = {
    unit: sorted({b for s in stocks for b in barcodes_of_stock[s]})
    for (unit, stocks) in unit_stocks.items()
}
unit_strain = dict(zip(viral_library[well_holds], viral_library["strain"], strict=True))

missing_viral_cols = sorted(set(addtl_viral_cols) - set(viral_library.columns))
if missing_viral_cols:
    raise ValueError(
        f"'addtl_viral_cols_to_show' of {single_well} names {missing_viral_cols}, which "
        f"are not columns of {snakemake.input.viral_library}"
    )

discarded_barcodes = check_discarded_barcodes(
    snakemake.params.discarded_barcodes,
    barcode_class["barcode"].tolist(),
)

assayed_units = set(samples.loc[~samples["is_empty"], "expected"])
unknown_units = sorted(assayed_units - set(unit_barcodes))
if unknown_units:
    raise ValueError(
        f"the samples CSVs of {single_well} infect wells with {well_holds}s that are not "
        f"in {snakemake.input.viral_library}, so they have no barcodes to count: "
        f"{unknown_units}"
    )

# Where each well holds a barcode, a samples CSV may also name the strain beside it. That
# column is not what the analysis reads, so it is checked against the library rather than
# being left to disagree with it unnoticed.
if well_holds == "barcode":
    _named_strains = pd.concat(
        [
            df.loc[
                df["barcode"].notna() & df["strain"].notna(),
                ["well", "barcode", "strain"],
            ].assign(plate=plate)
            for (plate, df) in samples_by_plate.items()
            if "strain" in df.columns
        ]
        or [pd.DataFrame(columns=["plate", "well", "barcode", "strain"])],
        ignore_index=True,
    )
    _disagree = _named_strains[
        _named_strains["strain"] != _named_strains["barcode"].map(unit_strain)
    ]
    if len(_disagree):
        raise ValueError(
            f"the samples CSVs of {single_well} name a strain for these barcodes that is "
            f"not the strain {snakemake.input.viral_library} gives them: "
            + ", ".join(
                f"{r.plate} {r.well} {r.barcode} ({r.strain!r}, library "
                f"{unit_strain.get(r.barcode)!r})"
                for r in _disagree.itertuples()
            )
        )

# The unit each library barcode belongs to, which is the well it should have come from. A
# barcode expected in an assayed unit's well belongs to that unit; any other belongs to its
# own, which no well here was infected with.
barcode_unit = dict(
    zip(viral_library["barcode"], viral_library[well_holds], strict=True)
)
_claimed = {}
for unit in sorted(assayed_units):
    for barcode in unit_barcodes[unit]:
        if barcode in _claimed:
            raise ValueError(
                f"the barcode {barcode} is expected in the wells of both {_claimed[barcode]} "
                f"and {unit}, sharing a stock, so a read of it cannot be traced to one well"
            )
        _claimed[barcode] = unit
barcode_unit.update(_claimed)
n_expected_barcodes = pd.Series({u: len(bs) for (u, bs) in unit_barcodes.items()})

# How a well is labeled on the charts' shared axis. When wells hold strains the label is the
# strain; when they hold barcodes the strain names the well too, a strain with several
# stocks being told apart by the stock's position among its strain's in library order.
if well_holds == "strain":
    unit_label = {strain: strain for strain in viral_library["strain"]}
else:
    _stocks = (
        viral_library[["strain", "stock_id"]].drop_duplicates().reset_index(drop=True)
    )
    _by_strain = _stocks.groupby("strain")
    _stock_label = {
        stock_id: strain if n_stocks == 1 else f"{strain} [bc {number}]"
        for (stock_id, strain, n_stocks, number) in zip(
            _stocks["stock_id"],
            _stocks["strain"],
            _by_strain["stock_id"].transform("nunique"),
            _by_strain.cumcount() + 1,
            strict=True,
        )
    }
    unit_label = {
        barcode: _stock_label[stock_id]
        for (barcode, stock_id) in zip(
            viral_library["barcode"], viral_library["stock_id"], strict=True
        )
    }

# Extra per-strain columns of the viral library carried into the tooltips, the tables and
# the output CSV. They describe the strain rather than the barcode, so a strain carrying
# several barcodes must give the same value on each of them.
strain_annotations = viral_library[["strain", *addtl_viral_cols]].drop_duplicates()
if strain_annotations["strain"].nunique() != len(strain_annotations):
    raise ValueError(
        f"the 'addtl_viral_cols_to_show' {addtl_viral_cols} are not constant across the "
        f"barcodes of every strain in {snakemake.input.viral_library}"
    )

# How each unit fared in the library carried forward: kept with all of its expected
# barcodes, kept with some, or dropped. By barcode, which is what the two libraries share.
unit_status = pd.DataFrame(
    {
        "unit": list(unit_barcodes),
        "assayed": [u in assayed_units for u in unit_barcodes],
        "n_expected_barcodes": [len(bs) for bs in unit_barcodes.values()],
    }
)
if has_final_library:
    (final_viral_library_csv,) = snakemake.input.final_viral_library
    _final_barcodes = set(
        pd.read_csv(final_viral_library_csv, dtype={"barcode": str})["barcode"]
    )
    unit_status = unit_status.assign(
        n_expected_barcodes_final=[
            len(set(bs) & _final_barcodes) for bs in unit_barcodes.values()
        ],
        final_library=lambda x: np.where(
            x["n_expected_barcodes_final"] == 0,
            "dropped from the library",
            np.where(
                x["n_expected_barcodes_final"] < x["n_expected_barcodes"],
                "kept, some barcodes dropped",
                "kept",
            ),
        ),
        carried_forward=lambda x: (x["n_expected_barcodes_final"] > 0).map(
            CARRIED_FORWARD_LABELS
        ),
    )
final_cols = (
    ["n_expected_barcodes_final", "final_library", "carried_forward"]
    if has_final_library
    else []
)

# --- read the counts ---------------------------------------------------------------------

plate_wells = list(zip(samples["plate"], samples["well"], strict=True))
counts = read_per_well_csvs(
    snakemake.input.counts, "_counts", ["barcode", "count"], plate_wells
)[["plate", "well", "barcode", "count"]]
if set(counts["barcode"]) != set(barcode_class["barcode"]):
    raise ValueError(
        "barcodes in the counts do not match those in the viral library and "
        "neutralization standard set"
    )
counts = counts.merge(barcode_class, on="barcode", validate="many_to_one").merge(
    samples[["plate", "well", "expected"]], on=["plate", "well"], validate="many_to_one"
)
# the unit a library barcode belongs to, null for a neutralization-standard barcode, which
# belongs to no well; a read belongs where this equals the well's own `expected`
counts["barcode_unit"] = counts["barcode"].map(barcode_unit)

fates = read_per_well_csvs(
    snakemake.input.fates, "_fates", ["fate", "count"], plate_wells
)
unexplained_fates = sorted(
    set(fates[fates["count"] > 0]["fate"]) - set(FATE_EXPLANATIONS)
)
if unexplained_fates:
    raise ValueError(
        f"the fates {unexplained_fates} have no explanation in `FATE_EXPLANATIONS`"
    )

invalid, discarded = split_discarded(
    read_per_well_csvs(
        snakemake.input.invalid,
        "_invalid",
        [
            "barcode",
            "count",
            "closest_valid_barcode",
            "closest_valid_barcode_hamming_distance",
        ],
        plate_wells,
    ),
    discarded_barcodes,
)

# --- the wells assayed ---------------------------------------------------------------------


def warned_wells(rows, reads_column=None):
    """The wells of `well_reads` rows a warning applies to, as the table it shows."""
    columns = {"plate": "plate", "well": "well", "unit_label": well_holds}
    if reads_column:
        columns[reads_column] = reads_column.replace("_", " ")
    return rows[list(columns)].rename(columns=columns)


_n_units = len(unit_barcodes)
_n_empty = int(samples["is_empty"].sum())
# A unit that no well was infected with, other than one whose barcodes are all in the tube
# of a unit that was, the passenger barcode of a shared tube being assayed with it.
_not_assayed = unit_status[
    ~unit_status["assayed"]
    & ~unit_status["unit"].map(lambda u: set(unit_barcodes[u]) <= set(_claimed))
]
report.md("""
    ## Experimental summary

    The tiles below summarize key experimental parameters.
    """)
report.stat_tiles(
    {
        "plates": len(plates),
        "infected wells": len(samples) - _n_empty,
        "blank wells": _n_empty,
        **({"dilution factor": dilutions[0]} if has_dilution else {}),
        f"{well_holds}s assayed": f"{int(unit_status['assayed'].sum())} / {_n_units}",
    },
    caption=f"Viral library: {Path(snakemake.input.viral_library).name}",
)
report.reserve_warnings()


def library_units(rows, columns=None):
    """The units of `unit_status` rows a library warning applies to, as the table it shows.

    Where wells hold barcodes each is shown beside its label, which the barcodes of one tube
    share. `columns` maps any further columns to show to their headers.

    """
    named = (
        {"label": "strain", "unit": "barcode"}
        if well_holds == "barcode"
        else {"label": well_holds}
    )
    named.update(columns or {})
    return (
        rows.assign(label=rows["unit"].map(unit_label))
        .sort_values(["label", "unit"])[list(named)]
        .rename(columns=named)
    )


# Warnings about the library rather than the wells, listed after those about the wells.
library_warnings = []
if len(_not_assayed):
    library_warnings.append(
        (
            "of the viral library not assayed here",
            library_units(
                _not_assayed,
                {"final_library": "final library"} if has_final_library else None,
            ),
        )
    )
if has_final_library:
    _assayed = unit_status.query("assayed")
    _dropped = _assayed.query("n_expected_barcodes_final == 0")
    _partial = _assayed.query("final_library == 'kept, some barcodes dropped'")
    if len(_dropped):
        library_warnings.append(
            (
                f"assayed here not carried through to {Path(final_viral_library_csv).name}",
                library_units(_dropped),
            )
        )
    if len(_partial):
        library_warnings.append(
            (
                (
                    f"carried through to {Path(final_viral_library_csv).name} with fewer "
                    "barcodes than expected here"
                ),
                library_units(
                    _partial,
                    {
                        "n_expected_barcodes": "barcodes expected here",
                        "n_expected_barcodes_final": "barcodes carried through",
                    },
                ),
            )
        )

# --- tally the reads of each well --------------------------------------------------------

# How far each unmatched barcode is from what belongs in the well it turned up in: the
# barcodes of the unit infected into it, and those of the neutralization standard spiked into
# it. Within `MAX_MISREAD_HAMMING` of either, the read is a misread of something that
# belongs there; further away it cannot be put down to either. The distance to what the well
# expects is computed a well at a time, so each barcode is measured against its own well's.
invalid = invalid.merge(
    samples[["plate", "well", "expected", "is_empty"]],
    on=["plate", "well"],
    validate="many_to_one",
)
invalid["expected_barcode_distance"] = pd.Series(np.nan, index=invalid.index)
for _expected, _group in invalid.groupby("expected", sort=False):
    invalid.loc[_group.index, "expected_barcode_distance"] = min_hamming(
        _group["barcode"], unit_barcodes[_expected]
    )
# an empty well expects no barcodes of the library, so it is left out of the loop above by
# its null `expected`, and the standard is the only thing its barcodes are measured against
if invalid.loc[~invalid["is_empty"], "expected_barcode_distance"].isnull().any():
    raise ValueError(
        "an unmatched barcode was not measured against the barcodes expected in its well"
    )
invalid["neut_standard_distance"] = min_hamming(
    invalid["barcode"], neut_standard_barcodes
)
invalid["expected_distance"] = invalid[
    ["expected_barcode_distance", "neut_standard_distance"]
].min(axis=1)

well_reads = (
    samples.assign(
        unit=well_holds,
        # an empty well expects none, rather than an unknown number
        n_expected_barcodes=lambda x: x["expected"]
        .map(n_expected_barcodes)
        .fillna(0)
        .astype(int),
        strain=lambda x: x["expected"].map(unit_strain),
        unit_label=lambda x: x["expected"]
        .map(unit_label)
        .where(~x["is_empty"], f"({EMPTY_WELL}) " + x["plate"] + ":" + x["well"]),
    )
    .merge(strain_annotations, on="strain", how="left", validate="many_to_one")
    .merge(
        unit_status[["unit", *final_cols]].rename(columns={"unit": "expected"}),
        on="expected",
        how="left",
        validate="many_to_one",
    )
)
# Every chart puts one well on each step of the shared axis, so two wells sharing a label
# would be drawn on one step. The same unit can be in a well of more than one plate, so
# those labels say which plate.
_repeated = well_reads["unit_label"].duplicated(keep=False)
well_reads.loc[_repeated, "unit_label"] = (
    well_reads["unit_label"] + " (" + well_reads["plate"] + ")"
)
if well_reads["unit_label"].nunique() != len(well_reads):
    raise ValueError("two wells of the group have the same label for the charts' axis")

_well_index = pd.MultiIndex.from_frame(well_reads[["plate", "well"]])


def reads_where(frame, mask):
    """Sum `frame`'s counts per well where `mask` holds, as a column of `well_reads`."""
    return reads_per_well(frame[mask], _well_index)


_is_expected = counts["barcode_unit"] == counts["expected"]
_misread = invalid["expected_distance"] <= MAX_MISREAD_HAMMING
well_reads["neut_standard_reads"] = reads_where(counts, counts["neut_standard"])
well_reads["expected_reads"] = reads_where(counts, _is_expected)
well_reads["other_library_reads"] = reads_where(
    counts, ~counts["neut_standard"] & ~_is_expected
)
well_reads["unmatched_near_expected_reads"] = reads_where(invalid, _misread)
well_reads["unmatched_far_expected_reads"] = reads_where(invalid, ~_misread)
# Thrown out rather than put in a category, but kept as a count so that the reads of each
# well still add up to what the pipeline parsed.
well_reads["discarded_reads"] = reads_where(
    discarded, pd.Series(True, index=discarded.index)
)
well_reads[NON_NEUT_STANDARD_READS] = well_reads[list(COMPOSITION_CATEGORIES)].sum(
    axis=1
)
if well_reads.loc[well_reads["is_empty"], "expected_reads"].any():
    raise ValueError("an empty well was credited with reads of its own material")
well_reads["parsed_reads"] = (
    well_reads[NON_NEUT_STANDARD_READS] + well_reads["neut_standard_reads"]
)

check_fates_add_up(
    fates,
    {
        fate: well_reads[columns].sum(axis=1).set_axis(_well_index)
        for (fate, columns) in [
            (
                "valid barcode",
                ["neut_standard_reads", "expected_reads", "other_library_reads"],
            ),
            (
                "invalid barcode",
                [
                    "unmatched_near_expected_reads",
                    "unmatched_far_expected_reads",
                    "discarded_reads",
                ],
            ),
        ]
    },
)

well_reads["too_shallow"] = well_reads["parsed_reads"] < min_parsed_reads_per_well
if well_reads["too_shallow"].any():
    report.warning(
        f"with fewer than {min_parsed_reads_per_well} parsed reads",
        warned_wells(well_reads.query("too_shallow"), "parsed_reads"),
    )

# every fraction is of the non-neut-standard parsed reads, so the categories sum to one
for column in COMPOSITION_CATEGORIES:
    well_reads[f"frac_{column}"] = well_reads[column] / well_reads[
        NON_NEUT_STANDARD_READS
    ].where(well_reads[NON_NEUT_STANDARD_READS] > 0)
_category_sums = well_reads[[f"frac_{c}" for c in COMPOSITION_CATEGORIES]].sum(axis=1)
if not np.allclose(_category_sums[well_reads[NON_NEUT_STANDARD_READS] > 0], 1):
    raise ValueError("the categories of a well do not account for all of its reads")

# What does not belong in the well at all, which is what the wells are ordered by: the other
# library material plus the unmatched reads too far from anything the well should hold.
well_reads["frac_foreign"] = well_reads[
    ["frac_other_library_reads", "frac_unmatched_far_expected_reads"]
].sum(axis=1, min_count=1)
well_reads["frac_unmatched"] = well_reads[
    ["frac_unmatched_near_expected_reads", "frac_unmatched_far_expected_reads"]
].sum(axis=1, min_count=1)
well_reads["flagged_other_library"] = (
    well_reads["frac_other_library_reads"] > max_other_library_frac
)
well_reads["flagged_unmatched"] = well_reads["frac_unmatched"] > max_unmatched_frac

# --- trace where material that belongs elsewhere came from --------------------------------

# The wells each unit was grown in, at most one per plate, which is what lets a read that
# belongs elsewhere be traced back to where it was.
_unit_homes = {}
for _row in samples[~samples["is_empty"]].itertuples():
    _unit_homes.setdefault(_row.expected, []).append((_row.plate, _row.well))

_other_unit_reads = counts[
    ~counts["neut_standard"] & ~_is_expected & (counts["count"] > 0)
].copy()
_origins = [
    trace_origin(row.plate, row.well, _unit_homes.get(row.barcode_unit, []))
    for row in _other_unit_reads.itertuples()
]
_other_unit_reads[["origin", "source_well", "wells_away", "came_from"]] = pd.DataFrame(
    _origins,
    index=_other_unit_reads.index,
    columns=["origin", "source_well", "wells_away", "came_from"],
)

# The reads of each foreign unit in each well, with where that unit was grown and its share
# of the well. Grouped by the unit rather than by where it came from, since `wells_away` is
# null for a unit not grown on this plate and pandas drops rows with a null grouping key: the
# origin columns are constant within a (well, unit) as a unit has one home well per plate.
other_unit_totals = (
    _other_unit_reads.groupby(["plate", "well", "barcode_unit"], as_index=False)
    .aggregate(
        reads=pd.NamedAgg("count", "sum"),
        origin=pd.NamedAgg("origin", "first"),
        source_well=pd.NamedAgg("source_well", "first"),
        wells_away=pd.NamedAgg("wells_away", "first"),
        came_from=pd.NamedAgg("came_from", "first"),
    )
    .rename(columns={"barcode_unit": "other_unit"})
    .merge(
        well_reads[["plate", "well", NON_NEUT_STANDARD_READS]],
        on=["plate", "well"],
        validate="many_to_one",
    )
    .assign(
        other_unit_label=lambda x: x["other_unit"].map(unit_label),
        fraction_of_reads=lambda x: x["reads"] / x[NON_NEUT_STANDARD_READS],
    )
    .sort_values(["reads", "other_unit"], ascending=[False, True])
)


def top_other_units(wells, n_top):
    """Get the `n_top` foreign units contributing most of each of `wells`' reads."""
    return (
        other_unit_totals.merge(
            wells[["plate", "well"]], on=["plate", "well"], validate="many_to_one"
        )
        .groupby(["plate", "well"])
        .head(n_top)
    )


# nearest first, with the origins that are not a distance last, so the stacking and the
# legend both read from "came from next door" towards "came from somewhere else"
origin_order = [o for o in all_origins() if o in set(other_unit_totals["origin"])]

# --- the reads-per-well chart -------------------------------------------------------------

# The infected wells, ordered by how much of each is not its own material, largest first.
# Computed here and handed to every panel as an explicit sort, so the panels need not carry
# the quantity they are sorted by.
infected = well_reads[~well_reads["is_empty"]]
well_order = infected.sort_values(
    ["frac_foreign", "unit_label"], ascending=[False, True]
)["unit_label"].tolist()

# What describes a well rather than a bar, one row per well, looked up by each panel rather
# than repeated on each of its rows. The key is unique, so the join is one to one.
well_info = well_reads[
    [
        "unit_label",
        "plate",
        "well",
        "expected",
        "strain",
        *addtl_viral_cols,
        *final_cols,
    ]
]
assert well_info["unit_label"].is_unique
well_lookup = alt.LookupData(
    data=well_info,
    key="unit_label",
    fields=[c for c in well_info.columns if c != "unit_label"],
)


well_tooltips = [
    lookup_tooltip(well_info, c)
    for c in ["plate", "well", "expected", "strain", *addtl_viral_cols]
]
if has_final_library:
    well_tooltips.append(lookup_tooltip(well_info, "final_library"))

_hover = well_hover()
_stroke = carried_forward_stroke(has_final_library, panels.legend_orient)


report.md(f"""
    ## Barcode composition of each well

    The multi-panel plot below can be used to examine the barcode composition of each well
    and help flag any QC issues that may need further attention. Each bar is one infected
    well, labeled by the {well_holds} it was infected with. Hover over a bar to outline
    that well in all four panels.
    """)

# The invalid barcode reads were checked above to be the unmatched plus the discarded ones,
# so they can be shown split in two.
fates_shown = pd.concat(
    [
        fates.loc[
            fates["fate"] != "invalid barcode", ["plate", "well", "fate", "count"]
        ],
        *(
            well_reads[["plate", "well"]].assign(
                fate=fate, count=well_reads[columns].sum(axis=1)
            )
            for (fate, columns) in [
                (
                    "invalid barcode",
                    ["unmatched_near_expected_reads", "unmatched_far_expected_reads"],
                ),
                ("discarded barcode", ["discarded_reads"]),
            ]
        ),
    ],
    ignore_index=True,
)
_infected_plate_wells = infected[["plate", "well", "unit_label"]]
fates_plotted = narrow_for_altair(
    fates_shown[fates_shown.groupby("fate")["count"].transform("sum") > 0].merge(
        _infected_plate_wells, on=["plate", "well"], validate="many_to_one"
    )[["unit_label", "fate", "count"]]
)
# the legend lists the fates in this order, and the notes explain them in it too
fate_order = sorted(fates_plotted["fate"].unique(), reverse=True)
composition_long = infected.melt(
    id_vars=["unit_label", NON_NEUT_STANDARD_READS],
    value_vars=list(COMPOSITION_CATEGORIES),
    var_name="category",
    value_name="reads",
).assign(
    fraction=lambda x: x["reads"] / x[NON_NEUT_STANDARD_READS],
    # the bars stack in the order the categories are declared in, which the label alone
    # cannot convey to Vega, so the position in that order rides along to sort on
    category_order=lambda x: x["category"].map(
        {c: i for (i, c) in enumerate(COMPOSITION_CATEGORIES)}
    ),
    category=lambda x: x["category"].map(COMPOSITION_CATEGORIES),
)
composition_long = narrow_for_altair(
    composition_long.drop(columns=NON_NEUT_STANDARD_READS)
)
origins_per_well = (
    other_unit_totals.groupby(["plate", "well", "origin"], as_index=False)
    .aggregate(reads=pd.NamedAgg("reads", "sum"))
    .merge(
        infected[["plate", "well", "unit_label", NON_NEUT_STANDARD_READS]],
        on=["plate", "well"],
        validate="many_to_one",
    )
    .assign(
        fraction_of_reads=lambda x: x["reads"] / x[NON_NEUT_STANDARD_READS],
        origin_order=lambda x: x["origin"].map(
            {origin: i for (i, origin) in enumerate(origin_order)}
        ),
    )[["unit_label", "origin", "origin_order", "reads", "fraction_of_reads"]]
)
origins_per_well = narrow_for_altair(origins_per_well)
# The reads of the "valid barcode" fate are the neutralization standard and every library
# barcode; a well with none of them has no fraction.
viral_fraction = narrow_for_altair(
    infected[["unit_label", "expected_reads", "neut_standard_reads"]].assign(
        valid_reads=lambda x: (
            x["expected_reads"]
            + infected["other_library_reads"]
            + x["neut_standard_reads"]
        ),
        viral_fraction=lambda x: (
            x["expected_reads"] / x["valid_reads"].where(x["valid_reads"] > 0)
        ),
    )
)

# a bar is a band of each well's step thick, along the wells' axis
_bar = panels.bar_band(0.85)
_panels = [
    alt.Chart(fates_plotted)
    .transform_lookup(lookup="unit_label", from_=well_lookup)
    .mark_bar(**_bar)
    .encode(
        panels.value_channel(
            "count:Q", title="reads", scale=alt.Scale(nice=False, padding=3)
        ),
        panels.well_axis(well_order, labeled=panels.labels_wells(0, 4)),
        alt.Fill(
            "fate:N",
            sort=fate_order,
            scale=alt.Scale(scheme="tableau10"),
            legend=alt.Legend(orient=panels.legend_orient, columns=1),
        ),
        alt.Order("fate:N", sort="descending"),
        stroke=_stroke,
        strokeWidth=stroke_width(_hover, STROKE_WIDTH_HOVERED, STROKE_WIDTH),
        tooltip=[*well_tooltips, "fate:N", "count:Q"],
    )
    .properties(**panels.panel_size(170)),
    alt.Chart(viral_fraction)
    .transform_lookup(lookup="unit_label", from_=well_lookup)
    .mark_bar(**_bar)
    .encode(
        panels.value_channel(
            "viral_fraction:Q",
            title="viral fraction of valid reads",
            scale=alt.Scale(domain=[0, 1], nice=False),
            axis=alt.Axis(tickCount=4),
        ),
        panels.well_axis(well_order, labeled=panels.labels_wells(1, 4)),
        stroke=_stroke,
        strokeWidth=stroke_width(_hover, STROKE_WIDTH_HOVERED, STROKE_WIDTH),
        tooltip=[
            *well_tooltips,
            alt.Tooltip("viral_fraction:Q", format=".3g"),
            "expected_reads:Q",
            "neut_standard_reads:Q",
            "valid_reads:Q",
        ],
    )
    .properties(**panels.panel_size(100)),
    alt.Chart(composition_long)
    .transform_lookup(lookup="unit_label", from_=well_lookup)
    .mark_bar(**_bar)
    .encode(
        panels.value_channel(
            "fraction:Q",
            title="fraction of the non-neut-standard reads",
            stack=True,
            # to one, so that a bar's length is read against the whole of a well
            scale=alt.Scale(domain=[0, 1], nice=False),
            axis=alt.Axis(tickCount=4),
        ),
        panels.well_axis(well_order, labeled=panels.labels_wells(2, 4)),
        alt.Fill(
            "category:N",
            title="barcode category",
            sort=list(COMPOSITION_CATEGORIES.values()),
            scale=alt.Scale(scheme="dark2"),
            legend=alt.Legend(
                labelLimit=LABEL_LIMIT, orient=panels.legend_orient, columns=1
            ),
        ),
        alt.Order("category_order:Q", sort="ascending"),
        stroke=_stroke,
        strokeWidth=stroke_width(_hover, STROKE_WIDTH_HOVERED, STROKE_WIDTH),
        tooltip=[
            *well_tooltips,
            "category:N",
            "reads:Q",
            alt.Tooltip("fraction:Q", format=".3g"),
        ],
    )
    .properties(**panels.panel_size(200)),
    alt.Chart(origins_per_well)
    .transform_lookup(lookup="unit_label", from_=well_lookup)
    .mark_bar(**_bar)
    .encode(
        panels.value_channel(
            "fraction_of_reads:Q",
            title="fraction from elsewhere in the library",
            stack=True,
            scale=alt.Scale(nice=False, padding=3),
            axis=alt.Axis(tickCount=3),
        ),
        panels.well_axis(well_order, labeled=panels.labels_wells(3, 4)),
        alt.Fill(
            "origin:N",
            title="came from",
            sort=origin_order,
            scale=alt.Scale(scheme="set1"),
            legend=alt.Legend(
                labelLimit=LABEL_LIMIT, orient=panels.legend_orient, columns=1
            ),
        ),
        alt.Order("origin_order:Q", sort="ascending"),
        stroke=_stroke,
        strokeWidth=stroke_width(_hover, STROKE_WIDTH_HOVERED, STROKE_WIDTH),
        tooltip=[
            *well_tooltips,
            "origin:N",
            "reads:Q",
            alt.Tooltip("fraction_of_reads:Q", format=".3g"),
        ],
    )
    .properties(**panels.panel_size(130)),
]
report.chart(
    panels.configure(
        panels.concat(*_panels, spacing=20)
        .add_params(_hover)
        .resolve_scale(**{panels.strain_axis: "shared"}, fill="independent")
    )
)
with report.notes("How to read this plot"):
    # the discarded barcodes are split out of the invalid ones, so they are not among them
    _fate_explanations = FATE_EXPLANATIONS | (
        {
            "invalid barcode": (
                "the read contains a barcode sequence not present in the viral library, "
                "the neutralization standard set, or the discarded barcode set."
            ),
            "discarded barcode": (
                "the read contains a barcode present in the discarded barcode set or "
                f"within {MAX_MISREAD_HAMMING} nt of one."
            ),
        }
        if discarded_barcodes
        else {}
    )
    report.md("""
        The **first panel** classifies every read sequenced for the well into
        the following categories:
        """)
    report.md(
        "\n".join(f" - **{fate}**: {_fate_explanations[fate]}" for fate in fate_order)
    )
    report.md(f"""
        The **second panel** takes every valid barcode read (from the "valid barcode"
        category in the first panel) and plots the fraction of these that map to the
        expected {well_holds}. This is a rough proxy for the titer of the expected virus and
        can help with further QC interpretation, especially for low titer strains.
        """)
    _category_labels = list(COMPOSITION_CATEGORIES.values())
    report.md(f"""
        The **third panel** classifies every barcode found in the well into
        the following categories:

         - **{_category_labels[0]}**: {
            "the barcode maps correctly to the infecting strain"
            if well_holds == "strain"
            else "the barcode is the infecting barcode, or another barcode carried in the "
            "same tube"}.
         - **{_category_labels[1]}**: the barcode maps to another {well_holds} in the
           library.
         - **{_category_labels[2]}**: the barcode is within {MAX_MISREAD_HAMMING} nt
           of an expected barcode or a neutralization standard barcode. These are likely
           sequencing errors as opposed to contamination.
         - **{_category_labels[3]}**: the barcode is further than
           {MAX_MISREAD_HAMMING} nt from every expected and neutralization standard
           barcode.

        Note that this excludes neutralization standard barcodes{
            " and discarded barcodes" if discarded_barcodes else ""}, and values are plotted
        as the fraction of non-neutralization-standard reads.

        The **fourth panel** takes every barcode in the "{_category_labels[1]}"
        category of the third panel and classifies it into the following categories:
        """)
    _belongs = (
        "the barcode belongs to a strain expected in a well"
        if well_holds == "strain"
        else "the barcode is expected in a well"
    )
    _origin_explanations = {
        distance_origin(1): (
            f"{_belongs} one well away in any direction on the same plate, including "
            "diagonally."
        ),
        **{
            distance_origin(d): f"{_belongs} on the same plate, {d} wells away."
            for d in range(2, MAX_NAMED_DISTANCE + 1)
        },
        distance_origin(MAX_NAMED_DISTANCE + 1): (
            f"{_belongs} on the same plate but over {MAX_NAMED_DISTANCE} wells away."
        ),
        ANOTHER_PLATE: f"{_belongs} on a different plate.",
        NOT_ASSAYED: (
            "the barcode belongs to a strain in the library but not assayed on the plates "
            "measured here."
            if well_holds == "strain"
            else "the barcode is in the library but not assayed on the plates measured "
            "here."
        ),
    }
    assert list(_origin_explanations) == all_origins()
    report.md(
        "\n".join(
            f" - **{origin}**: {_origin_explanations[origin]}"
            for origin in origin_order
            # only what the legend shows, which leaves out what is in empty wells alone
            if origin in set(origins_per_well["origin"])
        )
    )
    if samples["is_empty"].any():
        report.md(
            "Empty wells are not plotted but appear in the tables below if over a "
            f"threshold of {100 * max_other_library_frac:g}% (library contamination) or "
            f"{100 * max_unmatched_frac:g}% (unknown barcodes)."
        )


def dropped_rows(frame):
    """Mask of the rows of `frame` whose unit was dropped, or None with no final library."""
    if not has_final_library:
        return None
    return frame["carried_forward"] == CARRIED_FORWARD_LABELS[False]


# --- discarded barcodes --------------------------------------------------------------------

if discarded_barcodes:
    report.md("""
        ## Discarded barcodes

        The following barcodes were intentionally marked for discard:
        """)
    report.table(
        pd.DataFrame(
            {
                "barcode": list(discarded_barcodes),
                "reason": list(discarded_barcodes.values()),
            }
        )
    )
    _with_discards = well_reads[well_reads["discarded_reads"] > 0]
    if not len(_with_discards):
        report.md(
            "No wells contained these barcodes or barcodes within "
            f"{MAX_MISREAD_HAMMING} nt of them."
        )
    else:
        report.md(f"""
            Any reads matching these barcodes or within {MAX_MISREAD_HAMMING} nt
            were discarded, for a total of {int(well_reads["discarded_reads"].sum())} reads
            across the {len(well_reads)} wells analyzed. A summary table of this is below:
            """)
        report.table(
            _with_discards.assign(
                fraction=lambda x: x["discarded_reads"]
                / (x["parsed_reads"] + x["discarded_reads"]).where(
                    x["parsed_reads"] + x["discarded_reads"] > 0
                )
            )[["plate", "well", "unit_label", "discarded_reads", "fraction"]].rename(
                columns={
                    "unit_label": well_holds,
                    "discarded_reads": "discarded reads",
                    "fraction": "fraction of all parsed reads",
                }
            ),
            marked=dropped_rows(_with_discards),
            marked_color=DROPPED_COLOR,
        )


# --- wells with library contamination -----------------------------------------------------

flagged_other_library = well_reads.query("flagged_other_library").sort_values(
    "frac_other_library_reads", ascending=False
)
report.md(f"""
    ## Wells with library contamination

    Below is a table showing wells where more than {100 * max_other_library_frac:g}% of
    the non-neutralization-standard reads are unexpected library barcodes. Read counts are
    included and should be considered, especially in cases where the expected {well_holds}
    failed to grow.
    """)
if not len(flagged_other_library):
    report.md("No well is over the threshold.")
else:
    _table = flagged_other_library[
        [
            "plate",
            "well",
            "unit_label",
            *addtl_viral_cols,
            *(["carried_forward"] if has_final_library else []),
            NON_NEUT_STANDARD_READS,
            "expected_reads",
            "other_library_reads",
            "frac_other_library_reads",
        ]
    ].merge(
        top_other_units(flagged_other_library, N_TOP_BARCODES).assign(
            # a well on another plate is named with its plate, being no distance away
            source=lambda x: x["came_from"].where(
                x["origin"] == ANOTHER_PLATE, x["source_well"]
            )
        )[
            [
                "plate",
                "well",
                "other_unit_label",
                "fraction_of_reads",
                "source",
                "origin",
            ]
        ],
        on=["plate", "well"],
        validate="one_to_many",
    )
    _top_other = {
        "other_unit_label": well_holds,
        "fraction_of_reads": "fraction of non-neut-standard reads",
        "source": "well",
        "origin": "distance",
    }
    report.table(
        _table.drop(columns=["carried_forward"], errors="ignore").rename(
            columns={
                "unit_label": well_holds,
                NON_NEUT_STANDARD_READS: "non-neut-standard reads",
                "expected_reads": "expected reads",
                "other_library_reads": "other library reads",
                "frac_other_library_reads": "fraction unexpected library barcodes",
            }
        ),
        marked=dropped_rows(_table),
        marked_color=DROPPED_COLOR,
        group=(f"top other {well_holds}s", _top_other),
    )

# --- wells with unknown barcodes ----------------------------------------------------------

flagged_unmatched = well_reads.query("flagged_unmatched").sort_values(
    "frac_unmatched", ascending=False
)
report.md(f"""
    ## Wells with unknown barcodes

    Below is a table showing every well where more than {100 * max_unmatched_frac:g}% of
    the non-neutralization-standard reads contain barcodes not specified in the library{
        ", the discarded barcode set," if discarded_barcodes else ""} or the
    neutralization standard set.
    """)
if not len(flagged_unmatched):
    report.md("No well is over the threshold.")
else:
    _top_unmatched = (
        invalid.merge(
            flagged_unmatched[["plate", "well"]],
            on=["plate", "well"],
            validate="many_to_one",
        )
        .sort_values(["count", "barcode"], ascending=[False, True])
        .groupby(["plate", "well"])
        .head(N_TOP_BARCODES)
        .merge(
            barcode_class.rename(columns={"barcode": "closest_valid_barcode"}),
            on="closest_valid_barcode",
            how="left",
            validate="many_to_one",
        )
        .merge(
            well_reads[
                [
                    "plate",
                    "well",
                    "unit_label",
                    NON_NEUT_STANDARD_READS,
                    "frac_unmatched",
                    *final_cols,
                ]
            ],
            on=["plate", "well"],
            validate="many_to_one",
        )
        .assign(
            fraction_of_reads=lambda x: x["count"] / x[NON_NEUT_STANDARD_READS],
            # Close enough to a known barcode and it is that barcode misread, so the thing
            # it came off is named; further away the closest match is a coincidence rather
            # than an origin, so only the distance is given.
            relation=lambda x: (
                "sequencing error off "
                + x["strain"].where(~x["neut_standard"].astype(bool), "neut standard")
            ).where(
                x["closest_valid_barcode_hamming_distance"] <= MAX_MISREAD_HAMMING,
                x["closest_valid_barcode_hamming_distance"].astype(str)
                + " nt from any known barcode",
            ),
        )
        .sort_values(
            ["frac_unmatched", "plate", "well", "count", "barcode"],
            ascending=[False, True, True, False, True],
        )
    )
    report.table(
        _top_unmatched[
            [
                "plate",
                "well",
                "unit_label",
                NON_NEUT_STANDARD_READS,
                "frac_unmatched",
                "barcode",
                "fraction_of_reads",
                "relation",
            ]
        ].rename(
            columns={
                "unit_label": f"{well_holds} in well",
                NON_NEUT_STANDARD_READS: "non-neut-standard reads",
                "frac_unmatched": "fraction of well unknown",
                "fraction_of_reads": "fraction of non-neut-standard reads",
                "relation": "relation to closest barcode",
            }
        ),
        marked=dropped_rows(_top_unmatched),
        marked_color=DROPPED_COLOR,
        group=(
            "top unknown barcodes",
            [
                "barcode",
                "fraction of non-neut-standard reads",
                "relation to closest barcode",
            ],
        ),
    )

# --- estimated titer of what each well held -----------------------------------------------

report.md(f"""
    ## Estimated titer of each {well_holds}

    These single virus per well infection experiments allow for an estimated titer of each
    {well_holds}. This is simply computed as the ratio of reads containing the expected
    barcode to the neutralization standard reads in a well{
        ", summed over all of the strain's barcodes (the number of barcodes is shown "
        "when hovering over points)"
        if well_holds == "strain"
        else ""}. The titer is relative, in arbitrary units rather than infectious units.
    These titers should be
    interpreted carefully as they assume each infection was performed in the linear range
    (where there is a linear relationship between the amount of input virus and resulting
    spike-in normalized counts), which may not always be true for every virus in these
    experiments. Additionally, the total counts for both the expected {well_holds} and the
    neutralization standard are displayed when hovering over points. These can be used to
    help assess the precision of these estimates.
    """)

well_reads["titer_imprecise"] = ~well_reads["is_empty"] & (
    well_reads["neut_standard_reads"] < min_neut_standard_reads_for_titer
)
if well_reads["titer_imprecise"].any():
    report.warning(
        f"with fewer than {min_neut_standard_reads_for_titer} neutralization "
        "standard reads",
        warned_wells(well_reads.query("titer_imprecise"), "neut_standard_reads"),
    )

# A well with no standard reads has no titer, there being nothing to divide by, and an empty
# well has none to measure rather than one of zero.
well_reads["titer"] = (
    well_reads["expected_reads"]
    / well_reads["neut_standard_reads"].where(well_reads["neut_standard_reads"] > 0)
).where(~well_reads["is_empty"])
_no_titer = well_reads[~well_reads["is_empty"] & well_reads["titer"].isna()]
if len(_no_titer):
    report.warning("with no neutralization standard reads", warned_wells(_no_titer))

titers = well_reads.loc[
    well_reads["titer"].notna(),
    [
        "unit_label",
        "titer",
        "expected_reads",
        "neut_standard_reads",
        "parsed_reads",
        "n_expected_barcodes",
        "too_shallow",
        "titer_imprecise",
        "flagged_other_library",
        "flagged_unmatched",
    ],
]
_zero_titer = well_reads["titer"] == 0
report.section_warnings(
    [
        {
            "text": "with zero titer",
            "details": warned_wells(well_reads[_zero_titer], "expected_reads"),
            "note": (
                "These wells have no reads of their own, and a logarithmic axis cannot "
                f"place zero, so they are absent from the {panels.positions[-1]} panel "
                "of the plot below."
            ),
        }
    ]
    if _zero_titer.any()
    else []
)

if not len(titers):
    report.md("No well has a titer to report.")
else:
    _titer_order = titers.sort_values(["titer", "unit_label"], ascending=[False, True])[
        "unit_label"
    ].tolist()
    _titer_base = alt.Chart(narrow_for_altair(titers)).transform_lookup(
        lookup="unit_label", from_=well_lookup
    )
    # One frame for both panels, so that each row is in the chart once. A logarithmic scale
    # cannot place zero, and a zero in its data breaks its domain, so that panel drops them.
    titer_panel_bases = {
        "linear": _titer_base,
        "log": _titer_base.transform_filter("datum.titer > 0"),
    }
    report.chart(
        panels.configure(
            panels.concat(
                *(
                    titer_panel_bases[scale_type]
                    .mark_point(filled=True, fill="gray", size=25)
                    .encode(
                        panels.value_channel(
                            "titer:Q",
                            title=f"estimated titer, {scale_type} scale",
                            scale=alt.Scale(type=scale_type),
                        ),
                        panels.well_axis(
                            _titer_order, labeled=panels.labels_wells(i_panel, 2)
                        ),
                        stroke=_stroke,
                        strokeWidth=stroke_width(
                            _hover, POINT_STROKE_WIDTH_HOVERED, POINT_STROKE_WIDTH
                        ),
                        tooltip=[
                            *well_tooltips,
                            alt.Tooltip("titer:Q", format=".3g"),
                            "expected_reads:Q",
                            "neut_standard_reads:Q",
                            "parsed_reads:Q",
                            "n_expected_barcodes:Q",
                            "too_shallow:N",
                            "titer_imprecise:N",
                            "flagged_other_library:N",
                            "flagged_unmatched:N",
                        ],
                    )
                    .properties(**panels.panel_size(220))
                    for (i_panel, scale_type) in enumerate(["linear", "log"])
                ),
                spacing=20,
            )
            .add_params(_hover)
            .resolve_scale(**{panels.strain_axis: "shared"})
        )
    )

# --- expected barcode recovery -----------------------------------------------------------

report.md("""
    ## Expected barcode recovery

    This section examines the expected barcode recovery in each well. This can be helpful
    to examine the representation of each barcode in the viral stock when there are
    multiple barcodes mapping to an individual strain.
    """)

# Built by reindexing onto the barcodes each well expects rather than by filtering the
# counts, so that a barcode with no row in a well reads as zero rather than vanishing.
expected_barcode_reads = (
    pd.DataFrame(
        [
            {
                "plate": row.plate,
                "well": row.well,
                "expected": row.expected,
                "barcode": barcode,
            }
            for row in samples[~samples["is_empty"]].itertuples()
            for barcode in unit_barcodes[row.expected]
        ],
        columns=["plate", "well", "expected", "barcode"],
    )
    .merge(
        counts[["plate", "well", "barcode", "count"]],
        on=["plate", "well", "barcode"],
        how="left",
        validate="one_to_one",
    )
    .rename(columns={"count": "reads"})
    .fillna({"reads": 0})
    .astype({"reads": int})
    .merge(
        well_reads[["plate", "well", "unit_label", "expected_reads", *final_cols]],
        on=["plate", "well"],
        validate="many_to_one",
    )
    .assign(
        unit=well_holds,
        frac_of_expected=lambda x: x["reads"]
        / x["expected_reads"].where(x["expected_reads"] > 0),
    )
    .sort_values(["plate", "well", "barcode"])
    .reset_index(drop=True)
)

if len(expected_barcode_reads):
    report.chart(
        panels.configure(
            alt.Chart(
                narrow_for_altair(
                    expected_barcode_reads.assign(
                        # colored by rank in its well, as which barcode a segment is means
                        # nothing from one well to the next and is in the tooltip
                        rank=lambda x: x.sort_values(
                            ["reads", "barcode"], ascending=[False, True]
                        )
                        .groupby(["plate", "well"])
                        .cumcount()
                        + 1
                    )[["unit_label", "barcode", "rank", "reads", "frac_of_expected"]]
                )
            )
            .transform_lookup(lookup="unit_label", from_=well_lookup)
            .mark_bar(**_bar)
            .encode(
                panels.value_channel(
                    "frac_of_expected:Q", title="fraction of the well's expected reads"
                ),
                panels.well_axis(well_order, labeled=True),
                alt.Fill(
                    "rank:O",
                    title="barcode rank in well",
                    scale=alt.Scale(scheme="tableau10"),
                    legend=alt.Legend(orient=panels.legend_orient),
                ),
                alt.Order("rank:O", sort="ascending"),
                tooltip=[
                    *well_tooltips,
                    "barcode:N",
                    "reads:Q",
                    alt.Tooltip("frac_of_expected:Q", format=".3g"),
                ],
            )
            .properties(**panels.panel_size(300))
        )
    )

_absent = expected_barcode_reads[expected_barcode_reads["reads"] == 0]
if len(_absent):
    report.warning(
        "with no reads in expected well",
        _absent[["barcode", "plate", "well"]],
        shows="barcode",
    )
    report.md(
        "Additionally, the following table shows barcodes with no reads in the "
        "expected well:"
    )
    report.table(
        _absent[["plate", "well", "unit_label", "barcode"]].rename(
            columns={"unit_label": well_holds}
        ),
        marked=dropped_rows(_absent),
        marked_color=DROPPED_COLOR,
    )
else:
    report.md("Every barcode a well should hold produced at least one read.")

# --- write the outputs -----------------------------------------------------------------

for _text, _details in library_warnings:
    report.warning(_text, _details, shows=well_holds)
if report.warnings:
    print(f"Finished with {len(report.warnings)} warnings")

_worst_other_unit = top_other_units(well_reads, 1).rename(
    columns={
        "other_unit": "top_other_unit",
        "reads": "top_other_unit_reads",
        "source_well": "top_other_unit_well",
        "wells_away": "top_other_unit_wells_away",
        "origin": "top_other_unit_origin",
    }
)[
    [
        "plate",
        "well",
        "top_other_unit",
        "top_other_unit_reads",
        "top_other_unit_well",
        "top_other_unit_wells_away",
        "top_other_unit_origin",
    ]
]

well_composition = (
    well_reads.merge(
        _worst_other_unit, on=["plate", "well"], how="left", validate="one_to_one"
    )
    .sort_values(["plate", "position"])[
        [
            "plate",
            "well",
            "unit",
            "expected",
            "strain",
            *addtl_viral_cols,
            *(["dilution_factor"] if has_dilution else []),
            "n_expected_barcodes",
            *(c for c in final_cols if c != "carried_forward"),
            "parsed_reads",
            "discarded_reads",
            "neut_standard_reads",
            NON_NEUT_STANDARD_READS,
            *COMPOSITION_CATEGORIES,
            *(f"frac_{column}" for column in COMPOSITION_CATEGORIES),
            "frac_foreign",
            "frac_unmatched",
            "titer",
            "too_shallow",
            "titer_imprecise",
            "flagged_other_library",
            "flagged_unmatched",
            "top_other_unit",
            "top_other_unit_reads",
            "top_other_unit_well",
            "top_other_unit_wells_away",
            "top_other_unit_origin",
        ]
    ]
    .reset_index(drop=True)
)
well_composition.to_csv(
    snakemake.output.well_composition, index=False, float_format="%.4g"
)
expected_barcode_reads[
    ["plate", "well", "unit", "expected", "barcode", "reads", "frac_of_expected"]
].to_csv(snakemake.output.expected_barcode_reads, index=False, float_format="%.4g")
report.write(snakemake.output.html)

print(f"Wrote the report to {snakemake.output.html}")

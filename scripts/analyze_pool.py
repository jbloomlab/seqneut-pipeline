"""Measure the composition of a library pool sequenced over a serial dilution.

The pool was used to infect cells over a serial dilution and the barcodes sequenced. This
measures how much of the pool each unit of the library contributes, over the dilutions
where the counts still respond to being diluted, and how much virus the undiluted pool
holds against the neutralization standard.

What to do about an unbalanced pool is left to the project: nothing here computes a volume,
a ratio to add, or which units to leave out of a re-pool.

"""

import math
import sys
from pathlib import Path

import altair as alt
import pandas as pd
from miscellaneous_plate_charts import (
    LABEL_LIMIT,
    fate_fill,
    lookup_tooltip,
    order_fates,
)
from miscellaneous_plate_funcs import (
    FATE_EXPLANATIONS,
    MAX_MISREAD_HAMMING,
    PlateReport,
    check_discarded_barcodes,
    check_fates_add_up,
    plate_and_well,
    read_per_well_csvs,
    reads_per_well,
    split_discarded,
)
from pool_funcs import (
    IDEAL_SLOPE,
    analysis_wells,
    barcode_fractions,
    diverging_stocks,
    fit_dilution_windows,
    member_barcodes,
    pool_samples,
    pool_units,
    unit_shares,
)
from seqneut_funcs import narrow_for_altair, round_sig

# `noqa: SIM115` as this log file must stay open for the life of the script
sys.stderr = sys.stdout = open(snakemake.log[0], "w")  # noqa: SIM115

# factor by which the odds axis is padded beyond the wells plotted on it
ODDS_AXIS_PAD = 1.5

# height in pixels of the charts that plot a value against the dilution series
SERIES_CHART_HEIGHT = 300

# width in pixels of the dilution-series charts and of each well-QC panel
SERIES_CHART_WIDTH = 400
QC_PANEL_WIDTH = 250

# pixels per unit along the axis of the composition chart
UNIT_STEP = 10

# Whether a well passes a QC threshold, fixed per value so that passing reads the same in
# every panel and whether or not that panel holds a well of the other kind. From `tableau10`,
# as the other charts' colors are.
QC_PASS = "#4c78a8"
QC_FAIL = "#e45756"

# Pixels allowed for a line of a legend's title, enough for the longest configuration key
# named in one. The default truncates well short of it, and a title is drawn past the width
# of its panel rather than wrapped, so the limit is what decides whether the key reads whole.
LEGEND_TITLE_LIMIT = 400


# --- read the configuration ---------------------------------------------------------------

pool = snakemake.wildcards.pool
pool_config = snakemake.params.pool_config
date = pool_config["date"]
plates = pool_config["miscellaneous_plates"]
pool_by = pool_config["analyze_by"]
equal_volume = pool_config["pooled_by_equal_volume"]
max_stock_share_ratio = pool_config["max_stock_share_ratio"]

print(f"Analyzing pool {pool} from {date}, by {pool_by}")

report = PlateReport(f"Composition of library pool {pool}")
report.md(f"""
    This document analyzes the pooling experiment - {pool}. In these experiments a pool of
    individual {pool_by}s is made and then infected on cells over a serial dilution. This
    was performed here on the miscellaneous
    {"plates" if len(plates) > 1 else "plate"}
    {", ".join(f"`{p}`" for p in plates)} from {date}. There are several goals with this
    analysis:

     - to verify that the pool contains all intended viral {pool_by}s
     - to assess how evenly the included viral {pool_by}s are balanced in this pool
     - to determine the dilutions over which the pool is in the linear range and recommend
       a dilution factor for future experiments that use this pool

    This analysis outputs CSV files with the information needed for any subsequent
    re-pooling calculations. The plots in the report below are interactive: mouse over for
    details.

    ## Overview
    """)
report.md(pool_config["description"])  # added separately as it is multi-line


# --- read the samples and the libraries ---------------------------------------------------

# The plates are replicates of one dilution series, so their wells are read as one series
# and labeled by the plate they came from. One plate can also hold several pools, each down
# its own series, which `pool_wells` picks this pool's wells out of.
samples = pd.concat(
    [
        pool_samples(
            pd.read_csv(csv, dtype={"well": str}).drop(columns="fastq"),
            pool_config["pool_wells"],
        ).assign(plate=plate)
        for (plate, csv) in zip(plates, snakemake.input.samples_csvs, strict=True)
    ],
    ignore_index=True,
)
samples_csvs = ", ".join(snakemake.input.samples_csvs)
if "dilution_factor" not in samples.columns:
    raise ValueError(
        f"{samples_csvs} lacks a 'dilution_factor' column, which says how much the pool "
        "was diluted in each well and so is what the linear range is found over"
    )
if samples["dilution_factor"].isna().any():
    raise ValueError(
        f"{samples_csvs} gives no 'dilution_factor' for the wells "
        f"{sorted(samples.loc[samples['dilution_factor'].isna(), 'well'])}"
    )
if "replicate" not in samples.columns:
    raise ValueError(
        f"{samples_csvs} lacks a 'replicate' column, which says which replicate of the "
        "dilution series each well belongs to"
    )
if samples["replicate"].isna().any():
    raise ValueError(
        f"{samples_csvs} gives no 'replicate' for the wells "
        f"{sorted(samples.loc[samples['replicate'].isna(), 'well'])}"
    )
plate_wells = list(samples[["plate", "well"]].itertuples(index=False, name=None))


def well_label(plate, well):
    """How a well is named in the report and the outputs.

    The plates of a pool are replicates of one dilution series and so name their wells
    alike, which makes a well's own name ambiguous as soon as there is more than one plate.
    A pool measured from a single plate is the ordinary case, so there the plate is left off
    rather than repeated down every axis.

    """
    return f"{plate}:{well}" if len(plates) > 1 else well


def for_this_pool(csv_files, suffix):
    """The per-well CSVs of this pool's wells, out of every well of its plates.

    The rule passes whole plates, of which only this pool's wells are wanted. They are
    filtered by path before being read rather than dropped after, as `read_per_well_csvs`
    requires the files it reads to be exactly the wells it is told to expect. The plate is
    matched as well as the well, as two plates of a pool name their wells alike.

    """
    wanted = set(plate_wells)
    return [f for f in csv_files if plate_and_well(f, suffix) in wanted]


viral_library = pd.read_csv(snakemake.input.viral_library)
neut_standard_barcodes = pd.read_csv(snakemake.input.neut_standard_set)[
    "barcode"
].tolist()

addtl_viral_cols = list(pool_config["addtl_viral_cols_to_show"])
missing_viral_cols = [c for c in addtl_viral_cols if c not in viral_library.columns]
if missing_viral_cols:
    raise ValueError(
        f"'addtl_viral_cols_to_show' names {missing_viral_cols}, which are not columns of "
        f"{snakemake.input.viral_library}"
    )
strain_annotations = viral_library[["strain", *addtl_viral_cols]].drop_duplicates()
if strain_annotations["strain"].nunique() != len(strain_annotations):
    raise ValueError(
        f"the 'addtl_viral_cols_to_show' {addtl_viral_cols} are not constant across the "
        f"barcodes of every strain in {snakemake.input.viral_library}"
    )

members = member_barcodes(viral_library, pool_config["pool_members"])
member_library = viral_library[viral_library["barcode"].isin(members)]
n_units = pool_units(member_library, pool_by)
even_share = 1 / n_units
report.md("""
    ## Experimental summary

    The tiles below summarize key experimental parameters.
    """)
# A replicate of the series is a plate's replicate: the plates of a pool are themselves
# replicates of one series, and each numbers its own replicates from one, so the pair is what
# distinguishes them rather than the column on its own.
replicates = sorted(
    set(samples[["plate", "replicate"]].itertuples(index=False, name=None))
)
report.stat_tiles(
    {
        "replicates": len(replicates),
        "dilutions": len(set(samples["dilution_factor"])),
        "strains pooled": f"{member_library['strain'].nunique()} / "
        f"{viral_library['strain'].nunique()}",
        "barcodes pooled": f"{len(members)} / {len(viral_library)}",
        # `pool_by` is the strains or the barcodes themselves unless the pool is by stock,
        # in which case the units are counted on their own tile
        **({f"{pool_by}s pooled": n_units} if pool_by == "stock" else {}),
    },
    captions=[
        f"Viral library: {Path(snakemake.input.viral_library).name}",
        f"Pooled by: {pool_by}",
    ],
)
# A dilution measured in fewer replicates than another is not an error, as a well can be
# left out or fail, but it is worth seeing: the series is then uneven across the range.
_replicates_per_dilution = samples.groupby("dilution_factor").apply(
    lambda g: g[["plate", "replicate"]].drop_duplicates().shape[0], include_groups=False
)
_ragged = _replicates_per_dilution[_replicates_per_dilution != len(replicates)]
if len(_ragged):
    report.warning(
        f"not measured in all {len(replicates)} replicates",
        _ragged.rename("replicates").reset_index(),
        shows="dilution",
        note=(
            "Certain dilution values were measured in fewer total replicates. Depending "
            "on where these fall, this could result in less confident estimation of the "
            "linear range or library composition."
        ),
    )
# What the units are is worth stating, as which one the pool is balanced by decides what
# an even share means; it is hidden, being background rather than a result.
with report.notes("What is the difference between strains, barcodes, and stocks?"):
    report.md("""
        - **barcode**: the base unit of sequence-based counting.
        - **strain**: a viral strain in the library. Each strain carries one or more
          barcodes.
        - **stock**: a physical tube of virus from the library. Each stock contains a
          single strain, but may carry several of that strain's barcodes if they were
          rescued together. Useful when you have a mix of individual barcode rescues and
          combined barcode rescues for some strains.
        """)

# Reserved here and filled at the end, once every warning has been collected: a warning is
# found as the analysis runs but is wanted where it will be read first.
report.reserve_warnings()


# --- read the counts ----------------------------------------------------------------------

counts = read_per_well_csvs(
    for_this_pool(snakemake.input.counts, "_counts"),
    "_counts",
    ["barcode", "count"],
    plate_wells,
)[["plate", "well", "barcode", "count"]]
whitelisted = set(viral_library["barcode"]) | set(neut_standard_barcodes)
if set(counts["barcode"]) != whitelisted:
    raise ValueError(
        "barcodes in the counts do not match those in the viral library and "
        "neutralization standard set"
    )

fates = read_per_well_csvs(
    for_this_pool(snakemake.input.fates, "_fates"),
    "_fates",
    ["fate", "count"],
    plate_wells,
)
unexplained_fates = sorted(
    set(fates[fates["count"] > 0]["fate"]) - set(FATE_EXPLANATIONS)
)
if unexplained_fates:
    raise ValueError(
        f"the fates {unexplained_fates} have no explanation in `FATE_EXPLANATIONS`"
    )

discarded_barcodes = check_discarded_barcodes(
    snakemake.params.discarded_barcodes, sorted(whitelisted)
)
invalid, discarded = split_discarded(
    read_per_well_csvs(
        for_this_pool(snakemake.input.invalid, "_invalid"),
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

# The counts and the invalid CSVs have to add back up to what the fates say became of the
# reads, or a tally below has dropped or double-counted one. `reads_per_well` returns a
# tally indexed from zero in the order of `well_index`, so each is labeled by plate and
# well for the comparison.
well_index = pd.MultiIndex.from_tuples(plate_wells, names=["plate", "well"])
check_fates_add_up(
    fates,
    {
        "valid barcode": reads_per_well(counts, well_index).set_axis(well_index),
        "invalid barcode": (
            reads_per_well(invalid, well_index) + reads_per_well(discarded, well_index)
        ).set_axis(well_index),
    },
)


# --- fates of the sequencing reads --------------------------------------------------------

report.md("""
    ## Read level QC

    Every sequencing read gets assigned to one of the "fate" categories below. This can be
    used to identify if something went wrong in the sequencing or parsing of the barcodes.
    """)

well_order = [
    well_label(p, w)
    for (p, w) in samples.sort_values(["dilution_factor", "plate", "well"])[
        ["plate", "well"]
    ].itertuples(index=False, name=None)
]
drawn_fates = fates[fates.groupby("fate")["count"].transform("sum") > 0].assign(
    well=lambda x: [
        well_label(p, w) for (p, w) in zip(x["plate"], x["well"], strict=True)
    ]
)
fate_order, _stack_order = order_fates(drawn_fates["fate"])
report.chart(
    alt.Chart(
        narrow_for_altair(drawn_fates, drop=["plate"]).assign(fate_order=_stack_order)
    )
    .encode(
        alt.X("count", title="reads", scale=alt.Scale(nice=False, padding=3)),
        alt.Y("well:N", title=None, sort=well_order),
        fate_fill(fate_order, "right"),
        alt.Order("fate_order:Q"),
        tooltip=["well:N", "fate:N", "count:Q"],
    )
    .mark_bar(height={"band": 0.85})
    .properties(height=alt.Step(UNIT_STEP), width=QC_PANEL_WIDTH)
    .configure_axis(grid=False)
    .configure_axisX(labelOverlap="greedy")
)
with report.notes("What each fate means"):
    report.md(
        "\n".join(f" - **{fate}**: {FATE_EXPLANATIONS[fate]}" for fate in fate_order)
    )


# --- which wells measure the pool ---------------------------------------------------------

report.md("""
    ## Well level QC
    """)

# before the section's text, which a reader otherwise has to get past to learn that
# something needs their attention
_well_qc_warnings_slot = report.reserve_section_warnings()

report.md(f"""
    Each well needs to have enough barcode counts to accurately assess the composition of
    the pool. A well fails if it averages fewer than
    {pool_config["min_avg_barcode_count_per_well"]} counts per barcode.

    There also needs to be sufficient counts of the neutralization standard barcodes for
    normalization and a well fails if fewer than
    {round_sig(100 * pool_config["min_neut_standard_frac_per_well"])}% of the total barcode
    counts are from neutralization standard barcodes. Additionally, as the dilution series
    progresses, the fraction of barcode counts from the neutralization standard should
    generally increase as the viral library becomes relatively more dilute. Discrepancies
    here could indicate pipetting errors when performing the serial dilution.
    """)

is_neut = counts["barcode"].isin(neut_standard_barcodes)
is_member = counts["barcode"].isin(members)
well_qc = (
    counts.groupby(["plate", "well"], as_index=False)
    .aggregate(
        avg_count=pd.NamedAgg("count", "mean"),
        total_count=pd.NamedAgg("count", "sum"),
    )
    .merge(
        counts[is_neut]
        .groupby(["plate", "well"], as_index=False)["count"]
        .sum()
        .rename(columns={"count": "neut_standard_count"}),
        on=["plate", "well"],
        how="left",
        validate="one_to_one",
    )
    .merge(
        counts[is_member]
        .groupby(["plate", "well"], as_index=False)["count"]
        .sum()
        .rename(columns={"count": "member_count"}),
        on=["plate", "well"],
        how="left",
        validate="one_to_one",
    )
    .fillna({"neut_standard_count": 0, "member_count": 0})
    .merge(
        samples[["plate", "well", "dilution_factor", "replicate"]],
        on=["plate", "well"],
        validate="one_to_one",
    )
    .assign(
        # What a well is called on an axis and in a table. `plate` and `well` are kept as
        # they are, being what a well is selected by.
        well_name=lambda x: [
            well_label(p, w) for (p, w) in zip(x["plate"], x["well"], strict=True)
        ],
        viral_count=lambda x: x["total_count"] - x["neut_standard_count"],
        # a well with no counts at all gets a null fraction rather than a zero division
        neut_standard_frac=lambda x: (
            x["neut_standard_count"] / x["total_count"].where(x["total_count"] > 0)
        ),
        fails_count_qc=lambda x: (
            x["avg_count"] < pool_config["min_avg_barcode_count_per_well"]
        ),
        # a null fraction cannot be compared, so treat it as failing
        fails_neut_standard_qc=lambda x: ~(
            x["neut_standard_frac"] >= pool_config["min_neut_standard_frac_per_well"]
        ),
        fails_qc=lambda x: x["fails_count_qc"] | x["fails_neut_standard_qc"],
        # The odds are what a dilution step should move: diluting leaves the standard alone
        # while dividing the viral counts, so the standard's counts over the viral counts
        # rise along a line of slope IDEAL_SLOPE on log-log axes for as long as the pool
        # responds. The fraction itself just flattens towards one.
        neut_standard_odds=lambda x: (
            x["neut_standard_count"] / x["viral_count"].where(x["viral_count"] > 0)
        ),
        # The counts also have to be positive for the odds and their log to be defined at
        # all, which is arithmetic rather than a further threshold.
        fittable=lambda x: ~x["fails_qc"]
        & (x["viral_count"] > 0)
        & (x["neut_standard_count"] > 0),
    )
)

empty_wells = sorted(well_qc.loc[well_qc["total_count"] <= 0, "well_name"])
# Collected as the section is built and written into the slot left at the top of it, these
# wells being what a reader needs to know about before reading the thresholds against them.
well_qc_warnings = []
if empty_wells:
    well_qc_warnings.append(
        {
            "text": "with zero barcode counts",
            "details": well_qc.loc[
                well_qc["total_count"] <= 0, ["well_name", "dilution_factor"]
            ].rename(columns={"well_name": "well"}),
        }
    )

qc_panels = []
# Each panel is colored by its own threshold, so each legend names the configuration key
# that sets it and the value it is set to: the two legends are drawn side by side and a
# shared title would read as one. A title is given as its lines, a single string being
# truncated at the width of a panel rather than wrapped.
for i, (qc_col, x_col, x_title, qc_title) in enumerate(
    [
        (
            "fails_count_qc",
            "avg_count",
            "average counts per barcode",
            [
                "fails min_avg_barcode_count_per_well",
                f"({pool_config['min_avg_barcode_count_per_well']} counts per barcode)",
            ],
        ),
        (
            "fails_neut_standard_qc",
            "neut_standard_frac",
            "fraction of counts from the standard",
            [
                "fails min_neut_standard_frac_per_well",
                (
                    f"({round_sig(100 * pool_config['min_neut_standard_frac_per_well'])}"
                    "% of counts from the standard)"
                ),
            ],
        ),
    ]
):
    # `dict.fromkeys` to drop a repeat while keeping the order: the panel's own `x_col` is
    # also one of the columns every panel shows, and a frame of duplicated columns is not
    # one `altair` can read
    panel_columns = list(
        dict.fromkeys(
            [
                "well_name",
                "dilution_factor",
                x_col,
                qc_col,
                "avg_count",
                "neut_standard_count",
                "member_count",
            ]
        )
    )
    qc_panels.append(
        alt.Chart(narrow_for_altair(well_qc[panel_columns]))
        .encode(
            alt.X(
                x_col,
                title=x_title,
                scale=alt.Scale(nice=False, padding=3),
                # few enough ticks that their labels cannot collide at this panel width
                axis=alt.Axis(tickCount=4),
            ),
            alt.Y(
                "well_name:N",
                title=None,
                sort=well_order,
                axis=alt.Axis(labelLimit=LABEL_LIMIT) if i == 0 else None,
            ),
            # Pinned to the value rather than taken from the scheme's order, the panels
            # coloring by their own threshold: a panel whose wells all pass or all fail
            # holds one value, which would otherwise take the first color of the scheme and
            # so mean the opposite of the same color in the panel beside it.
            alt.Fill(
                qc_col,
                title=qc_title,
                scale=alt.Scale(domain=[False, True], range=[QC_PASS, QC_FAIL]),
                # wide enough for the configuration key, which is longer than the default
                # limit and would otherwise be drawn with its end cut off
                legend=alt.Legend(titleLimit=LEGEND_TITLE_LIMIT),
            ),
            tooltip=[
                "well_name:N",
                "dilution_factor:Q",
                "avg_count:Q",
                "neut_standard_count:Q",
                "member_count:Q",
            ],
        )
        .mark_bar(height={"band": 0.85})
        .properties(height=alt.Step(UNIT_STEP), width=QC_PANEL_WIDTH)
    )
report.chart(
    # spaced apart so the last tick label of one panel cannot collide with the first of
    # the next, which `labelOverlap` cannot prevent as it works within one axis
    alt.hconcat(*qc_panels, spacing=45)
    # Each panel colors by its own threshold, so each keeps its own legend. The two fields
    # are both boolean, and `vega-lite` merges the legends of scales it can share, which
    # would draw one legend titled after whichever panel came first.
    .resolve_scale(y="shared", fill="independent")
    .resolve_legend(fill="independent")
    .configure_axis(grid=False)
    .configure_axisX(labelOverlap="greedy")
)
report.fill_section_warnings(_well_qc_warnings_slot, well_qc_warnings)


# --- identifying the linear range ---------------------------------------------------------

report.md("""
    ## Identifying the linear range
    """)

_linear_range_warnings_slot = report.reserve_section_warnings()

report.md(f"""
    The linear range for a pool is determined by examining how the share of neutralization
    standard counts changes across the dilution series. At high concentrations of the viral
    pool (low dilution factor), it is possible to saturate the cellular capacity to produce
    barcoded nucleic acid. In this regime, adding less virus doesn't necessarily result in a
    corresponding decrease in viral barcode counts and so the fraction of neutralization
    standard barcodes remains relatively consistent. As the pool is diluted further, there
    should be a point where a linear relationship (on a log-scale) is observed between the
    share of neutralization standard counts and the dilution factor.

    To determine this point, lines are fit to a series of
    {pool_config["linearity_window_dilutions"]} consecutive dilutions with a sliding window
    of 1 dilution. In ideal conditions this slope would be {IDEAL_SLOPE:+.0f}, meaning a
    2-fold dilution of the library pool would result in a 2-fold increase in the amount of
    neutralization standard. The first (most concentrated) dilution window with a slope of
    {IDEAL_SLOPE:+.0f} +/- {pool_config["linearity_slope_tol"]} is selected and the middle
    dilution from this window is the recommended dilution factor. This is quantified in the
    following table:
    """)

fittable_wells = well_qc[well_qc["fittable"]]
if not len(fittable_wells):
    raise ValueError(
        f"no well of pool {pool} passes QC with countable reads, so there is nothing to "
        "measure the pool from"
    )

windows = fit_dilution_windows(
    fittable_wells,
    window_dilutions=pool_config["linearity_window_dilutions"],
    slope_tol=pool_config["linearity_slope_tol"],
)
chosen = analysis_wells(
    windows, fittable_wells, given=pool_config["linear_range_wells"]
)
analysis_plate_wells = chosen.wells
analysis_well_labels = [well_label(p, w) for (p, w) in chosen.wells]
# No window of the series responds to dilution as it should, so there is no dilution at
# which what the pool holds can be read. The report says so and stops short of a
# composition rather than reading one from a dilution the series does not support.
no_linear_range = not chosen.wells
if no_linear_range:
    print(
        "No window of dilutions is within tolerance of the ideal slope, so no composition "
        "is read"
    )
else:
    print(
        f"Reading the composition from wells {analysis_well_labels} at "
        f"{chosen.dilution}-fold"
    )

if windows:
    report.table(
        pd.DataFrame(
            [
                {
                    "dilutions": f"{w['window'][0]}-{w['window'][-1]}",
                    "slope": round_sig(w["slope"]),
                    "standard error": round_sig(w["slope_se"]),
                    "off ideal by": round_sig(w["deviation"]),
                    "within tolerance": w["within_tol"],
                }
                for w in windows
            ]
        ),
        # the window the dilution read from is the middle of, which a reader looks for first
        bold=pd.Series([w["window"] == chosen.titer_window for w in windows]),
    )

# Limits from the wells rather than from everything drawn, as the ideal line runs the
# length of the series and reaches odds the wells never take.
odds_domain = [
    round_sig(fittable_wells["neut_standard_odds"].min() / ODDS_AXIS_PAD),
    round_sig(fittable_wells["neut_standard_odds"].max() * ODDS_AXIS_PAD),
]
dilution_axis = alt.X(
    "dilution_factor",
    scale=alt.Scale(type="log"),
    title="pool reciprocal dilution factor",
)
odds_axis = alt.Y(
    "neut_standard_odds",
    scale=alt.Scale(type="log", domain=odds_domain, nice=False),
    title="neutralization standard odds, frac / (1 - frac)",
)


def ideal_height_over(window):
    """The height a line of the ideal slope takes to pass through `window`'s wells.

    The intercept of `IDEAL_SLOPE` through the mean of the wells of `window` on log-log
    axes, which is what the odds of the undiluted pool would be if the pool responded to
    dilution ideally over that window.

    """
    in_window = fittable_wells[fittable_wells["dilution_factor"].isin(window)]
    return (
        in_window["neut_standard_odds"].map(math.log10)
        - IDEAL_SLOPE * in_window["dilution_factor"].map(math.log10)
    ).mean()


# The titer: how much virus the undiluted pool holds against the standard. The standard is
# the same spike-in in every well, so the viral-to-standard ratio measures the pool in units
# that do not depend on sequencing depth, and multiplying by the dilution recovers the
# undiluted pool. Only the line's height is fit, with the slope fixed at the ideal, over the
# window that gives the titer.
if chosen.titer_window is None:
    ideal_intercept = None
    titer_vs_neut_standard = None
else:
    ideal_intercept = ideal_height_over(chosen.titer_window)
    # The intercept is the odds of the undiluted pool, and the titer is how much virus that
    # is per unit of standard, which is the reciprocal.
    titer_vs_neut_standard = 10**-ideal_intercept


def ideal_line(dilutions, intercept):
    """Where a line of the ideal slope and height `intercept` sits at each of `dilutions`.

    The slope is fixed at `IDEAL_SLOPE` rather than fit, so this is what the odds would do
    if every dilution step divided the viral counts by that step; only the height is fit.
    Every dilution asked for is returned, the line being cut to the axis by the mark's own
    clipping: dropping the points that fall outside the odds the wells set would shorten
    the line to its endpoints inside them instead of running it to the edge of the plot.

    """
    line = pd.DataFrame({"dilution_factor": sorted(set(dilutions))}).assign(
        neut_standard_odds=lambda x: x["dilution_factor"].map(
            lambda d: 10 ** (intercept + IDEAL_SLOPE * math.log10(d))
        )
    )
    return narrow_for_altair(line)


series_layers = []

if ideal_intercept is not None:
    # A checkbox rather than an entry in the legend the windows are selected from: the
    # ideal line belongs to whichever window is selected rather than being one more window
    # to pick, so it is shown or hidden on its own. Off by default, so the plot opens on
    # the wells and the fits alone.
    ideal_param_name = "ideal_slope_line"
    ideal_shown = alt.param(
        name=ideal_param_name,
        value=False,
        bind=alt.binding_checkbox(name="ideal slope line "),
    )
    # Dashed, where the fits to the windows are solid: it is the slope the assay should
    # have rather than anything measured from the wells, and a reader comparing a fit
    # against it should not have to tell them apart by color alone. The dash is longer than
    # the chosen-dilution rule's, so the two do not read as the same annotation.
    ideal_dash = [6, 3]
    # One line per window rather than one path through every window's points. A line mark
    # draws a single connected path through the rows behind it, and only a grouping channel
    # splits it; the fits are split by the color they carry, which these have none of, being
    # black. Without this, selecting two windows joins their lines end to end, running the
    # path back across the series to start the next.
    ideal_group = alt.Detail("window:N")

# the ideal lines, built per window once the windows are fit and the selection exists
ideal_layers = []

series_layers += [
    # every well is drawn, but the limits come from the wells that are fit, so a failing
    # well far off the scale is clipped rather than compressing all the others
    alt.Chart(
        narrow_for_altair(
            well_qc[
                [
                    "well_name",
                    "dilution_factor",
                    "neut_standard_odds",
                    "neut_standard_frac",
                    "viral_count",
                    "neut_standard_count",
                    "fails_qc",
                ]
            ]
        )
    )
    .mark_circle(size=60, clip=True)
    .encode(
        dilution_axis,
        odds_axis,
        # Pinned rather than taken from the scheme's order, as the fitted lines now color
        # by their own scale: a plot whose wells all pass holds one value, which would
        # otherwise take the first color of the scheme and so mean the opposite of the
        # same color in the well-QC panels above.
        alt.Color(
            "fails_qc",
            title="fails a threshold",
            scale=alt.Scale(domain=[False, True], range=[QC_PASS, QC_FAIL]),
        ),
        tooltip=[
            "well_name:N",
            "dilution_factor:Q",
            "neut_standard_frac:Q",
            "neut_standard_odds:Q",
            "viral_count:Q",
            "neut_standard_count:Q",
        ],
    ),
    # the dilution the composition is read from, where there is one
    alt.Chart(
        pd.DataFrame({"dilution_factor": [] if no_linear_range else [chosen.dilution]})
    )
    .mark_rule(color="black", strokeDash=[2, 2])
    .encode(dilution_axis),
]


# Each window's own fitted line, chosen by clicking the legend. The selection opens on the
# window the composition is read from, so the view on load is the fit that decided the
# dilution; clicking another isolates it, and shift-clicking adds to the selection, which is
# what lets a few windows be compared without the rest of a long series in the way. The
# windows not selected stay on the plot faintly rather than vanishing, so that the shape of
# the whole series is still legible: a pool that breaks cleanly into responding to dilution
# has its later fits along one line, where one that only creeps toward it has them fanning.
def window_label(w):
    """How a window of dilutions is named in the plot.

    The window the composition is read from is marked, that being what a reader looks for
    first and what the colors otherwise leave them to work out. It is only a recommendation
    where the analysis made it: given `linear_range_wells`, the window is the one that
    happens to center the wells that were set.

    """
    dilutions = f"{w['window'][0]}-{w['window'][-1]}"
    if w["window"] != chosen.titer_window:
        return dilutions
    return f"{dilutions} ({'analyzed' if chosen.wells_given else 'recommended'})"


chosen_label = next(
    (window_label(w) for w in windows if w["window"] == chosen.titer_window), None
)
if windows and chosen_label is not None:
    all_dilutions = sorted(set(fittable_wells["dilution_factor"]))

    def fitted(w, dilutions):
        """Where window `w`'s fitted line sits at each of `dilutions`."""
        return pd.DataFrame(
            {
                "dilution_factor": dilutions,
                "neut_standard_odds": [
                    10 ** (w["intercept"] + w["slope"] * math.log10(d))
                    for d in dilutions
                ],
            }
        ).assign(window=window_label(w), slope=round_sig(w["slope"]))

    def fits_frame(dilutions_of):
        """Every window's fitted line, over the dilutions `dilutions_of` gives each.

        Clipped to the odds the wells set: a saturated window's line extended across the
        series runs far off the scale, and would otherwise compress every well into a band.

        """
        frame = pd.concat(
            [fitted(w, dilutions_of(w)) for w in windows], ignore_index=True
        )
        return frame[frame["neut_standard_odds"].between(*odds_domain)]

    window_fits = fits_frame(lambda w: w["window"])
    # The same line carried across the whole series, drawn faintly: it says where the
    # window's fit would put the pool at every other dilution, which is what shows a
    # saturated window's slope failing to describe the dilute end. Drawn for each window
    # selected, so that comparing two compares where each would put the whole series.
    window_extensions = fits_frame(lambda w: all_dilutions)

    # The ideal line of each window, at the height that carries the ideal slope through
    # that window's wells, so that selecting a window draws the reference its own fit is to
    # be read against rather than one pinned to the window the titer came from.
    def ideal_frame(dilutions_of):
        """Every window's ideal line, over the dilutions `dilutions_of` gives each."""
        return pd.concat(
            [
                ideal_line(dilutions_of(w), ideal_height_over(w["window"])).assign(
                    window=window_label(w)
                )
                for w in windows
            ],
            ignore_index=True,
        )

    ideal_fits = ideal_frame(lambda w: w["window"])
    ideal_extensions = ideal_frame(lambda w: all_dilutions)

    # the domain names every window, so that a window keeps its color however many are
    # selected and the legend lists them all to click
    window_fill = alt.Color(
        "window:N",
        title="window of dilutions",
        scale=alt.Scale(domain=[window_label(w) for w in windows]),
        legend=alt.Legend(
            # a line rather than the circle a symbol legend draws by default, the lines
            # being what it names and what it selects
            symbolType="stroke",
            symbolStrokeWidth=3,
        ),
    )
    window_selection = alt.selection_point(
        fields=["window"],
        # Clicking the legend, which is both the key to the colors and the control: there
        # is no multi-select input binding, and a window at a time is not enough on a long
        # series. Shift-click adds to the selection.
        bind="legend",
        # the chosen window, so the plot opens on the fit that placed the linear range
        value=[{"window": chosen_label}],
        name="window_fit",
    )
    # Which windows to draw. A legend binding always clears the selection on a click that
    # lands on no mark, which cannot be switched off, and an empty selection would
    # otherwise draw every window; this puts the plot back to the window the composition
    # is read from instead, so that clicking beside a legend entry returns to where the
    # plot started rather than to something never asked for.
    selected_windows = (
        "length(data('window_fit_store')) "
        "? indexof(window_fit.window, datum.window) >= 0 "
        # a `vega` expression, so the string constant needs quotes inside the Python
        # string; the label is dilutions and a word, so it holds none of its own
        f": datum.window === '{chosen_label}'"
    )
    # Only while a single window is drawn. The line is the reference that one fit is read
    # against, and one per window at once is several dashed lines at several heights
    # competing with the fits they are meant to qualify rather than a reference. An empty
    # store falls back to the window the composition is read from, which is a single window
    # too, so the line is there on load.
    one_window_selected = "length(data('window_fit_store')) <= 1"

    if ideal_intercept is not None:
        ideal_layers = [
            # Behind the fits and the wells, being the reference they are read against.
            # Three filters apply: the checkbox shows it at all, the count keeps it to a
            # single window, and the selection picks whose it is.
            alt.Chart(narrow_for_altair(ideal_extensions))
            .mark_line(
                color="black",
                opacity=0.25,
                strokeWidth=2,
                strokeDash=ideal_dash,
                clip=True,
            )
            .encode(dilution_axis, odds_axis, ideal_group, tooltip=["window:N"])
            .add_params(ideal_shown)
            # the param's own name, as a `vega` expression: passing the param itself makes
            # a selection predicate, which tests a row against a selection's values rather
            # than reading the boolean the checkbox sets
            .transform_filter(ideal_param_name)
            .transform_filter(one_window_selected)
            .transform_filter(selected_windows),
            # Full strength over the window's own dilutions, where its wells are claimed to
            # follow the ideal, and faint across the rest of the series, where they are
            # not. Both are dashed, so the opacity is what separates the two.
            alt.Chart(narrow_for_altair(ideal_fits))
            .mark_line(color="black", strokeWidth=2, strokeDash=ideal_dash, clip=True)
            .encode(dilution_axis, odds_axis, ideal_group, tooltip=["window:N"])
            .transform_filter(ideal_param_name)
            .transform_filter(one_window_selected)
            .transform_filter(selected_windows),
        ]

    series_layers += [
        # the extension first, so the fit over the window's own dilutions draws on top
        alt.Chart(narrow_for_altair(window_extensions))
        .mark_line(strokeWidth=2, opacity=0.3)
        .encode(
            dilution_axis,
            odds_axis,
            # no legend of its own, the layer below carrying the one that selects
            window_fill.copy().legend(None),
            tooltip=["window:N", "slope:Q"],
        )
        .transform_filter(selected_windows),
        alt.Chart(narrow_for_altair(window_fits))
        .mark_line(strokeWidth=2)
        .encode(
            dilution_axis,
            odds_axis,
            window_fill,
            tooltip=["window:N", "slope:Q"],
        )
        # Declared on this layer rather than on the layered chart: a legend binding
        # attaches to the legend of the scale its own view colors by, and the chart's
        # first view colors by the well QC. Declared there, the binding reaches the
        # wrong legend and silently does nothing.
        .add_params(window_selection)
        # Only the selected windows, rather than drawing the rest faintly: a line still
        # on the plot reads as one being shown, and the point of selecting one is to see
        # it against the wells by itself.
        .transform_filter(selected_windows),
    ]

# The dashed vertical marks the dilution read from, which the wells were set to on the
# explicit path and the fit recommends on the automatic one.
if no_linear_range:
    report.md("""
        The plot below displays the neutralization standard odds (fraction of standard / 1
        - fraction of standard) versus the reciprocal dilution factor on a log-log scale,
        with the linear fit to each dilution window. No window is within tolerance of the
        ideal slope, so no dilution is marked and the pool's composition is not read; the
        warning above this section gives the slope of every window.
        """)
else:
    _marked_dilution = (
        "the dilutions set in `linear_range_wells`"
        if chosen.wells_given
        else "the recommended dilution"
    )
    report.md(f"""
        The plot below displays the neutralization standard odds (fraction of standard / 1
        - fraction of standard) versus the reciprocal dilution factor on a log-log scale.
        The linear fit to the dilution windows are drawn with {_marked_dilution} marked by
        a dashed vertical line. The downstream compositional analysis of the pool is based
        on the wells at this chosen dilution factor.
        """)

if not chosen.wells_given and not no_linear_range:
    report.md(f"**Recommended dilution factor: {chosen.dilution}-fold**")

series_chart = (
    # the ideal lines first, being the reference the wells and the fits are read against
    alt.layer(*ideal_layers, *series_layers)
    .properties(height=SERIES_CHART_HEIGHT, width=SERIES_CHART_WIDTH)
    # The wells color by whether they pass a threshold and the fitted lines by which window
    # they belong to, which are unrelated. Both fields are nominal, and `vega-lite` merges
    # the legends of scales it can share, which would list a dilution window among the
    # pass and fail values under both titles at once.
    .resolve_scale(color="independent")
    .resolve_legend(color="independent")
)
report.chart(series_chart.configure_axisX(labelOverlap="greedy"))

with report.notes("How can I interact with this plot?"):
    report.md("""
        By default, the linear fit for the window containing the recommended dilution
        factor is displayed. The line is solid across the window, continuing fainter
        through the rest of the dilution range. Multiple windows can be selected by shift
        clicking the windows in the legend. The ideal slope (+1) line can be overlayed if a
        single window is selected for comparison purposes.
        """)

# Collected as the section is built and written into the slot left at the top of it, each
# being about the fit this section tabulates and so read before rather than after it.
linear_range_warnings = []

# A window spans consecutive dilutions of those that reach the fit, so a dilution whose every
# well is dropped closes the window over the gap: the window then covers a wider range of the
# series than its length says, and the composition is read from a middle that is not the
# middle of the dilutions the plate holds. The ragged-replicate warning above counts the
# replicates of a dilution and so does not catch one that is missing altogether.
_dropped_dilutions = sorted(
    set(samples["dilution_factor"]) - set(fittable_wells["dilution_factor"])
)
if _dropped_dilutions:
    linear_range_warnings.append(
        {
            "text": "left out of the fit, so a window of dilutions may span the gap",
            "details": pd.DataFrame(
                [
                    {
                        "dilution": dilution,
                        "wells": ", ".join(
                            sorted(
                                well_label(p, w)
                                for (p, w) in samples.loc[
                                    samples["dilution_factor"] == dilution,
                                    ["plate", "well"],
                                ].itertuples(index=False, name=None)
                            )
                        ),
                    }
                    for dilution in _dropped_dilutions
                ]
            ),
            "shows": "dilution",
            "note": (
                "Every well at these dilutions either failed a QC threshold or had no "
                "countable reads, so they are not fit. A window is a run of consecutive "
                "dilutions among those that are fit, so a window reported here can cover "
                "a wider span of the series than its number of dilutions suggests. Check "
                "the windows above against the dilutions the plate actually holds."
            ),
        }
    )

if chosen.wells_given and chosen.titer_window is None:
    linear_range_warnings.append(
        {
            "text": "reported without a titer",
            "details": pd.DataFrame(
                [
                    {
                        "wells": ", ".join(analysis_well_labels),
                        "dilution": chosen.dilution,
                    }
                ]
            ),
            "shows": "pool",
            "note": (
                f"No window of {pool_config['linearity_window_dilutions']} dilutions "
                f"that is within {pool_config['linearity_slope_tol']} of the ideal "
                f"slope has {chosen.dilution}-fold as its middle dilution, so there is "
                "no fit over these wells to take a titer from and "
                "`titer_vs_neut_standard` is empty. The composition below is still "
                "measured from them."
            ),
        }
    )

# a well kept by the fit may still have failed QC where the wells were given
failing_analysis = well_qc[
    pd.MultiIndex.from_frame(well_qc[["plate", "well"]]).isin(set(chosen.wells))
    & well_qc["fails_qc"]
]
if len(failing_analysis):
    linear_range_warnings.append(
        {
            "text": "used for library composition failed a QC threshold",
            "details": failing_analysis.assign(
                # the threshold beside the value it is compared with, so that the table says
                # which of the two the well failed and by how much
                min_avg_count=pool_config["min_avg_barcode_count_per_well"],
                min_neut_standard_frac=pool_config["min_neut_standard_frac_per_well"],
            )[
                [
                    "well_name",
                    "dilution_factor",
                    "avg_count",
                    "min_avg_count",
                    "neut_standard_frac",
                    "min_neut_standard_frac",
                ]
            ].rename(
                columns={"well_name": "well"}
            ),
        }
    )

if no_linear_range:
    _closest = min(windows, key=lambda w: w["deviation"], default=None)
    linear_range_warnings.append(
        {
            "text": "{verb} outside tolerance of the ideal slope, so the pool is "
            "not measured",
            "details": pd.DataFrame(
                [
                    {
                        "dilutions": f"{w['window'][0]}-{w['window'][-1]}",
                        "slope": round_sig(w["slope"]),
                        "standard error": round_sig(w["slope_se"]),
                        "off ideal by": round_sig(w["deviation"]),
                        "closest": w is _closest,
                    }
                    for w in windows
                ]
            ),
            "shows": "dilution window",
            "verb": "falls",
            "note": (
                "No window of dilutions has a slope within "
                f"{pool_config['linearity_slope_tol']} of the ideal "
                f"{IDEAL_SLOPE:+.0f}, so the pool never clearly responds to being diluted "
                "and there is no dilution at which its composition can be read. This "
                "report therefore stops at the wells and their QC: it gives no "
                "recommended dilution factor, no composition and no titer. Check the "
                "pool and the dilution series, and set `linear_range_wells` to read it "
                "from particular wells anyway."
            ),
        }
    )

report.fill_section_warnings(_linear_range_warnings_slot, linear_range_warnings)


# --- a pool with no linear range stops here -----------------------------------------------

# Everything below reads the wells at the chosen dilution, of which there are none. The
# outputs are still written, being declared by the rule: the summary records the pool and
# what was fit, with the analysis fields empty, and the per-barcode composition is written
# with its header and no rows.
if no_linear_range:
    pd.DataFrame(
        columns=[
            "barcode",
            "mean_fraction",
            "strain",
            *(["stock"] if "stock" in member_library.columns else []),
        ]
    ).to_csv(snakemake.output.barcode_composition, index=False)
    pd.DataFrame(
        [
            {
                "pool": pool,
                "date": date,
                "miscellaneous_plates": " ".join(plates),
                "neut_standard_set": pool_config["neut_standard_set"],
                "analyze_by": pool_by,
                "pooled_by_equal_volume": equal_volume,
                f"n_{pool_by}s": n_units,
                "n_strains": member_library["strain"].nunique(),
                "n_member_barcodes": len(members),
                "analysis_dilution": None,
                "analysis_wells": None,
                "linear_range_dilutions": None,
                "linear_range_strict": False,
                "analysis_wells_given": False,
                "titer_vs_neut_standard": None,
            }
        ]
    ).to_csv(snakemake.output.pool_summary, index=False, float_format="%.4g")
    report.write(snakemake.output.html)
    sys.exit(0)


# --- invalid and discarded barcodes -------------------------------------------------------

report.md("""
    ## Pool level QC
    """)

_pool_qc_warnings_slot = report.reserve_section_warnings()

report.md(f"""
    The following table classifies the barcodes of all reads that came from the chosen
    dilution factor wells. Reads with barcodes within {MAX_MISREAD_HAMMING} nucleotide of
    an expected barcode is likely sequencing error where as those further away may
    indicate contamination.
    """)


def in_analysis_wells(df):
    """Rows of `df` from the wells the composition is read from."""
    return df[
        pd.MultiIndex.from_frame(df[["plate", "well"]]).isin(set(analysis_plate_wells))
    ]


analysis_invalid = in_analysis_wells(invalid)
analysis_discarded = in_analysis_wells(discarded)
analysis_counts = in_analysis_wells(counts)
valid_reads = int(analysis_counts["count"].sum())
near = analysis_invalid["closest_valid_barcode_hamming_distance"] <= MAX_MISREAD_HAMMING
# every read of those wells, so that the rows below partition them and the fractions sum to
# one: a discarded barcode's reads are split out of the unmatched ones rather than counted
total_reads = (
    valid_reads + analysis_invalid["count"].sum() + analysis_discarded["count"].sum()
)
report.table(
    pd.DataFrame(
        [
            {
                "reads": reads,
                "count": count,
                "fraction": count / total_reads if total_reads else None,
            }
            for (reads, count) in [
                ("known barcode or neut standard", valid_reads),
                (
                    f"<= {MAX_MISREAD_HAMMING} nt of known barcode",
                    int(analysis_invalid.loc[near, "count"].sum()),
                ),
                (
                    f"> {MAX_MISREAD_HAMMING} nt of known barcode",
                    int(analysis_invalid.loc[~near, "count"].sum()),
                ),
            ]
            # only where such a barcode is configured, there being no such category of read
            # otherwise rather than none of them
            + (
                [("discarded barcode", int(analysis_discarded["count"].sum()))]
                if discarded_barcodes
                else []
            )
        ]
    )
)

if discarded_barcodes:
    report.md("""
        The following barcodes were intentionally marked for discard:
        """)
    report.table(
        pd.DataFrame(
            [
                {"barcode": barcode, "reason discarded": reason}
                for (barcode, reason) in sorted(discarded_barcodes.items())
            ]
        )
    )


# --- reads of strains the pool was not made from ------------------------------------------

non_members = (
    analysis_counts[
        ~analysis_counts["barcode"].isin(members)
        & ~analysis_counts["barcode"].isin(neut_standard_barcodes)
    ]
    .merge(viral_library[["barcode", "strain"]], on="barcode", validate="many_to_one")
    .groupby("strain", as_index=False)["count"]
    .sum()
    .query("count > 0")
    .sort_values("count", ascending=False)
)
if len(non_members):
    member_reads = int(
        analysis_counts.loc[analysis_counts["barcode"].isin(members), "count"].sum()
    )
    total = member_reads + int(non_members["count"].sum())
    # slotted into this section rather than the report's block: contamination is what the
    # section classifies the reads for, so the warning is read with that classification
    report.fill_section_warnings(
        _pool_qc_warnings_slot,
        [
            {
                "text": "detected",
                "details": non_members.assign(
                    fraction_of_viral_reads=lambda x: (x["count"] / total).map(
                        round_sig
                    )
                ).rename(columns={"count": "reads"}),
                "shows": "unexpected strain",
                "note": (
                    "There are some strains accounting for "
                    f"{round_sig(100 * non_members['count'].sum() / total)}% of all viral reads "
                    "that were not supposed to be included in this pool. Flagging them here and "
                    "leaving out of composition analysis. They are not excluded from the titer "
                    "fit, so `titer_vs_neut_standard` is biased slightly high by the "
                    "fraction above."
                ),
            }
        ],
    )


# --- each unit's share of the pool --------------------------------------------------------

fractions = barcode_fractions(counts, members, analysis_plate_wells)
shares = unit_shares(fractions, member_library, pool_by, equal_volume=equal_volume)
share_cols = ["measured_share"] + (
    ["divided_share"] if "divided_share" in shares.columns else []
)

unit_order = shares["unit"].tolist()
plotted_shares = shares.merge(
    strain_annotations, on="strain", how="left", validate="many_to_one"
)

# A bar reads its magnitude off its length, which is only proportional to the share on a
# linear axis, so the axis is linear even though a weak unit's bar is then short.
unit_axis = alt.Y(
    "unit:N", title=None, sort=unit_order, scale=alt.Scale(domain=unit_order)
)
# hovering a bar outlines that unit in both panels
unit_hover = alt.selection_point(
    name="unit_hover",
    fields=["unit"],
    on="pointerover",
    clear="pointerout",
    empty=False,
)
hover_stroke = {
    "stroke": alt.value("black"),
    "strokeWidth": alt.condition(unit_hover, alt.value(2), alt.value(0)),
}

share_panel = alt.layer(
    alt.Chart(narrow_for_altair(plotted_shares))
    .encode(
        alt.X(
            "measured_share",
            title=f"share of the pool per {pool_by}",
            scale=alt.Scale(nice=False, padding=3),
        ),
        unit_axis.copy().axis(alt.Axis(labelLimit=LABEL_LIMIT)),
        tooltip=[
            "unit:N",
            "strain:N",
            *[f"{c}:Q" for c in share_cols],
            *[lookup_tooltip(strain_annotations, c) for c in addtl_viral_cols],
        ],
        **hover_stroke,
    )
    .mark_bar(height={"band": 0.85}, fill="gray"),
    # one row, so the rule is drawn once rather than once per unit
    alt.Chart(pd.DataFrame({"even_share": [round_sig(even_share)]}))
    .mark_rule(color="red", strokeDash=[2, 2], strokeWidth=2)
    .encode(alt.X("even_share")),
)

# How a unit's share splits among its barcodes, which says whether one barcode is poorly
# represented in the rescued virus. Only where a unit has more than one barcode: under a
# `pool_by` of `barcode` every unit is one barcode and the panel is a column of full bars.
within_unit = fractions.merge(
    member_library[["barcode", "strain", *(["stock"] if pool_by == "stock" else [])]],
    on="barcode",
    validate="one_to_one",
).assign(unit=lambda x: x[pool_by])
within_unit = within_unit.assign(
    unit_share=lambda x: x.groupby("unit")["mean_fraction"].transform("sum"),
    barcode_share=lambda x: (
        x["mean_fraction"] / x["unit_share"].where(x["unit_share"] > 0)
    ).map(round_sig),
    barcode_number=lambda x: (x.groupby("unit").cumcount() + 1).astype(str),
)
several_barcodes = (
    pool_by != "barcode" and within_unit.groupby("unit")["barcode"].size().max() > 1
)

report.md("""
    ## Pool composition
    """)

_composition_warnings_slot = report.reserve_section_warnings()

report.md(f"""
    A perfectly balanced pool would have equal representation from each {pool_by},
    represented by the red dashed line on the barplot below. The degree of variation from
    this equal share should be assessed for each {pool_by} to determine whether and how to
    perform any subsequent re-pooling attempts.
    """)
if several_barcodes:
    report.md(f"""
        The relative representation of each barcode for any given {pool_by} can also be
        examined for cases where there are multiple barcodes.
        """)

# How far from even the pool actually is, which the plot shows but does not quantify: the
# spread decides whether a re-pool is worth doing and by how much. Over the units with reads,
# a unit with none having no finite fold change to the rest; they are counted beside it, and
# the warning below the plot names them.
_present = plotted_shares.loc[plotted_shares["measured_share"] > 0, "measured_share"]
_n_absent = int((plotted_shares["measured_share"] <= 0).sum())
if len(_present) > 1:
    _within_2fold = int(
        ((_present / even_share < 2) & (even_share / _present < 2)).sum()
    )
    _absent_caption = (
        f"Computed over the {len(_present)} {pool_by}s with reads; {_n_absent} with none "
        "are left out, having no finite fold change."
    )
    report.stat_tiles(
        {
            f"highest : lowest {pool_by}": (
                f"{_present.max() / _present.min():.1f}-fold"
            ),
            "highest vs even share": f"{_present.max() / even_share:.1f}x over",
            "lowest vs even share": f"{even_share / _present.min():.1f}x under",
            f"{pool_by}s within 2-fold of even": f"{_within_2fold} / {len(_present)}",
        },
        captions=[_absent_caption] if _n_absent else [],
    )

panels = [share_panel]
if several_barcodes:
    panels.append(
        alt.Chart(narrow_for_altair(within_unit))
        .encode(
            alt.X("barcode_share", title=f"share of the {pool_by} per barcode"),
            # the labels come from the share panel this axis is shared with
            unit_axis.copy().axis(None),
            alt.Fill("barcode_number:N", title="barcode"),
            tooltip=["unit:N", "barcode:N", "barcode_share:Q"],
            **hover_stroke,
        )
        .mark_bar(height={"band": 0.85})
    )

report.chart(
    alt.hconcat(
        *(
            p.properties(height=alt.Step(UNIT_STEP), width=QC_PANEL_WIDTH)
            for p in panels
        )
    )
    .add_params(unit_hover)
    .resolve_scale(y="shared")
    .configure_axis(grid=False)
    .configure_axisX(labelOverlap="greedy")
)

# A unit with no reads draws a bar of no length, so it is named here as well: on the plot it
# is indistinguishable from one too weak to see, and it is the case that matters most.
absent = plotted_shares[plotted_shares["measured_share"] <= 0]
if len(absent):
    report.md(f"""
        The following {pool_by}s were not detected at all in the chosen dilution factor
        wells.
        """)
    report.table(
        absent[
            ["unit", *(["strain"] if pool_by != "strain" else []), *addtl_viral_cols]
        ].rename(columns={"unit": pool_by})
    )

if pool_by == "stock" and equal_volume:
    diverging = diverging_stocks(shares, member_library, max_stock_share_ratio)
    if len(diverging):
        # slotted into this subsection rather than the report's block: it is about the
        # shares the subsection plots, and is read against them
        report.fill_section_warnings(
            _composition_warnings_slot,
            [
                {
                    "text": (
                        "whose stocks differ in measured share by more than "
                        f"{max_stock_share_ratio}-fold"
                    ),
                    "details": diverging.assign(
                        min_share=lambda x: x["min_share"].map(round_sig),
                        max_share=lambda x: x["max_share"].map(round_sig),
                        ratio=lambda x: x["ratio"].map(
                            lambda v: round_sig(v) if pd.notna(v) else None
                        ),
                    ),
                    "shows": "strain",
                    "note": (
                        "At least one viral strain is present in multiple stocks that "
                        "were pooled here. This warning is triggered when the difference "
                        "between stocks of the same strain in their measured share of "
                        f"the pool exceeds {max_stock_share_ratio}-fold. Some possible "
                        "explanations include: mistakes when pooling, inaccuracies in "
                        "titering, differences in viral growth. A blank ratio means the "
                        "weakest stock has no reads at all."
                    ),
                }
            ],
        )


# --- write the outputs --------------------------------------------------------------------

composition = (
    fractions.merge(
        member_library[
            [
                "barcode",
                "strain",
                *(["stock"] if "stock" in member_library.columns else []),
            ]
        ],
        on="barcode",
        validate="one_to_one",
    )
    .merge(strain_annotations, on="strain", how="left", validate="many_to_one")
    .sort_values("mean_fraction", ascending=False)
    .reset_index(drop=True)
)
if pool_by == "stock":
    composition = composition.merge(
        shares.rename(
            columns={"unit": "stock", **{c: f"stock_{c}" for c in share_cols}}
        )[["stock", *[f"stock_{c}" for c in share_cols]]],
        on="stock",
        how="left",
        validate="many_to_one",
    )
composition.to_csv(
    snakemake.output.barcode_composition, index=False, float_format="%.4g"
)

pd.DataFrame(
    [
        {
            "pool": pool,
            "date": date,
            "miscellaneous_plates": " ".join(plates),
            "neut_standard_set": pool_config["neut_standard_set"],
            "analyze_by": pool_by,
            "pooled_by_equal_volume": equal_volume,
            f"n_{pool_by}s": n_units,
            "n_strains": member_library["strain"].nunique(),
            "n_member_barcodes": len(members),
            "analysis_dilution": chosen.dilution,
            "analysis_wells": " ".join(analysis_well_labels),
            "linear_range_dilutions": (
                " ".join(str(d) for d in chosen.titer_window)
                if chosen.titer_window
                else None
            ),
            "linear_range_strict": chosen.strict,
            "analysis_wells_given": chosen.wells_given,
            "titer_vs_neut_standard": titer_vs_neut_standard,
        }
    ]
).to_csv(snakemake.output.pool_summary, index=False, float_format="%.4g")

report.write(snakemake.output.html)

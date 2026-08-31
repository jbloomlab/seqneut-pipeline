"""Process a plate to QC and convert counts to fraction infectivity, then fit curves."""

import pickle
import sys
import warnings

import altair as alt
import matplotlib.style
import neutcurve
import numpy
import pandas as pd
from neutcurve.colorschemes import CBMARKERS, CBPALETTE
from seqneut_funcs import yaml_rt
from seqneut_report import Report

# `noqa: SIM115` as this log file must stay open for the life of the script
sys.stderr = sys.stdout = open(snakemake.log[0], "w")  # noqa: SIM115

_ = alt.data_transformers.disable_max_rows()

# avoid clutter w RuntimeWarning during curve fitting
warnings.filterwarnings("ignore", category=RuntimeWarning)

# faster plotting of neut curves
matplotlib.style.use("fast")

count_csvs = snakemake.input.count_csvs
fate_csvs = snakemake.input.fate_csvs
qc_drops_yaml = snakemake.output.qc_drops
frac_infectivity_csv = snakemake.output.frac_infectivity_csv
fits_csv = snakemake.output.fits_csv
fits_pickle = snakemake.output.fits_pickle
viral_barcodes = snakemake.params.viral_barcodes
neut_standard_barcodes = snakemake.params.neut_standard_barcodes
samples = snakemake.params.samples
plate = snakemake.wildcards.plate
plate_params = snakemake.params.plate_params
curve_display_method = snakemake.params.curve_display_method
collapse_strain_barcodes = snakemake.params.collapse_strain_barcodes

report = Report(
    title="Process plate counts to get fraction infectivities and fit curves"
)
report.md("""
    This report analyzes a plate of sequencing-based neutralization assays.

    The plots are interactive, so you can mouseover points for details, use the
    mouse-scroll to zoom and pan, and use interactive dropdowns at the bottom of the
    plots.
    """)

# get thresholds turning lists into tuples as needed
manual_drops = {
    filter_type: [tuple(w) if isinstance(w, list) else w for w in filter_drops]
    for (filter_type, filter_drops) in plate_params["manual_drops"].items()
}
group = plate_params["group"]
qc_thresholds = plate_params["qc_thresholds"]
curvefit_params = plate_params["curvefit_params"]
curvefit_qc = plate_params["curvefit_qc"]
if "barcode_serum_replicates_ignore_curvefit_qc" in curvefit_qc:
    curvefit_qc["barcode_serum_replicates_ignore_curvefit_qc"] = [
        tuple(w) for w in curvefit_qc["barcode_serum_replicates_ignore_curvefit_qc"]
    ]

# whether this plate titrates a 'dilution_factor' or a 'concentration'; the
# value is also the name of the corresponding column in samples_df
dilution_factor_or_concentration = plate_params["dilution_factor_or_concentration"]
assert dilution_factor_or_concentration in {"dilution_factor", "concentration"}
concentration_units = plate_params["concentration_units"]
x_label = (
    "dilution factor"
    if dilution_factor_or_concentration == "dilution_factor"
    else f"concentration ({concentration_units})"
)

report.md(f"Processing `{plate}`")

samples_df = pd.DataFrame(plate_params["samples"])
report.md(f"Plate has {len(samples)} samples (wells)")
assert all(
    (len(samples_df) == samples_df[c].nunique())
    for c in ["well", "sample", "sample_noplate"]
)
assert len(samples_df) == len(
    samples_df.groupby(["serum_replicate", dilution_factor_or_concentration])
)
assert len(samples) == len(count_csvs) == len(fate_csvs) == len(samples_df)

for d, key, title in [
    (manual_drops, "manual_drops", "Data manually specified to drop:"),
    (qc_thresholds, "qc_thresholds", "QC thresholds applied to data:"),
    (curvefit_params, "curvefit_params", "Curve-fitting parameters:"),
    (curvefit_qc, "curvefit_qc", "Curve-fitting QC:"),
]:
    report.md(f"{title}")
    report.yaml({key: d})

# Set up dictionary to keep track of wells, barcodes, well-barcodes, and
# serum-replicates that are dropped:
qc_drops = {
    "wells": {},
    "barcodes": {},
    "barcode_wells": {},
    "barcode_serum_replicates": {},
    "serum_replicates": {},
}

assert set(manual_drops).issubset(
    qc_drops
), f"{manual_drops.keys()=}, {qc_drops.keys()}"

report.md("""
    ## Statistics on barcode-parsing for each sample
    Make interactive chart of the "fates" of the sequencing reads parsed for each sample
    on the plate.

    If most sequencing reads are not "valid barcodes", this could potentially indicate
    some problem in the sequencing or barcode set you are parsing.

    Potential fates are:

     - *valid barcode*: barcode that matches a known virus or neutralization standard, we
       hope most reads are this.
     - *invalid barcode*: a barcode with proper flanking sequences, but does not match a
       known virus or neutralization standard. If you  have a lot of reads of this type,
       it is probably a good idea to look at the invalid barcode CSVs (in the
       `./results/barcode_invalid/` subdirectory created by the pipeline) to see what
       these invalid barcodes are.
     - *unparseable barcode*: could not parse a barcode from this read as there was not a
       sequence of the correct length with the appropriate flanking sequence.
     - *invalid outer flank*: if using an outer upstream or downstream region
       (`upstream2` or `downstream2` for the
       [illuminabarcodeparser](https://jbloomlab.github.io/dms_variants/dms_variants.illuminabarcodeparser.html)),
       reads that are otherwise valid except for this outer flank. Typically you would be
       using `upstream2` if you have a plate index embedded in your primer, and reads with
       this classification correspond to a different index than the one for this plate.
     - *low quality barcode*: low-quality or `N` nucleotides in barcode, could indicate
       problem with sequencing.
     - *failed chastity filter*: reads that failed the Illumina chastity filter, if these
       are reported in the FASTQ (they may not be).

    Also, if the number of reads per sample is very uneven, that could indicate that you
    did not do a good job of balancing the different samples in the Illumina sequencing.
    """)

fates = (
    pd.concat([pd.read_csv(f).assign(sample=s) for f, s in zip(fate_csvs, samples)])
    .merge(samples_df, validate="many_to_one", on="sample")
    .assign(
        fate_counts=lambda x: x.groupby("fate")["count"].transform("sum"),
        sample_well=lambda x: x["sample_noplate"] + " (" + x["well"] + ")",
    )
    .query("fate_counts > 0")[  # only keep fates with at least one count
        [
            "fate",
            "count",
            "well",
            "serum_replicate",
            "sample_well",
            dilution_factor_or_concentration,
        ]
    ]
)

assert len(fates) == len(fates.drop_duplicates())

serum_replicates = sorted(fates["serum_replicate"].unique())
sample_wells = list(
    fates.sort_values(["serum_replicate", dilution_factor_or_concentration])[
        "sample_well"
    ]
)

serum_selection = alt.selection_point(
    fields=["serum_replicate"],
    bind=alt.binding_select(
        options=[None] + serum_replicates,
        labels=["all"] + serum_replicates,
        name="serum",
    ),
)

fates_chart = (
    alt.Chart(fates)
    .add_params(serum_selection)
    .transform_filter(serum_selection)
    .encode(
        alt.X("count", scale=alt.Scale(nice=False, padding=3)),
        alt.Y(
            "sample_well",
            title=None,
            sort=sample_wells,
        ),
        alt.Color("fate", sort=sorted(fates["fate"].unique(), reverse=True)),
        alt.Order("fate", sort="descending"),
        tooltip=fates.columns.tolist(),
    )
    .mark_bar(height={"band": 0.85})
    .properties(
        height=alt.Step(10),
        width=200,
        title=f"Barcode parsing for {plate}",
    )
    .configure_axis(grid=False)
)

report.chart(fates_chart)

report.md("""
    ## Read barcode counts and apply manually specified drops
    Read the counts per barcode, then apply any manually specified drops.
    """)

# get barcode counts
counts = (
    pd.concat([pd.read_csv(c).assign(sample=s) for c, s in zip(count_csvs, samples)])
    .merge(samples_df, validate="many_to_one", on="sample")
    .drop(columns=["replicate", "plate", "fastq"])
    .assign(sample_well=lambda x: x["sample_noplate"] + " (" + x["well"] + ")")
)

# classify barcodes as viral or neut standard
barcode_class = pd.concat(
    [
        pd.DataFrame(viral_barcodes).assign(neut_standard=False),
        pd.DataFrame(neut_standard_barcodes).assign(neut_standard=True, strain=pd.NA),
    ],
    ignore_index=True,
)

# merge counts and classification of barcodes
assert set(counts["barcode"]) == set(barcode_class["barcode"])
counts = counts.merge(barcode_class, on="barcode", validate="many_to_one")
assert set(sample_wells) == set(counts["sample_well"])
assert set(serum_replicates) == set(counts["serum_replicate"])

counts_qc_1 = counts.copy()
for filter_type, filter_drops in manual_drops.items():
    report.md(f"Dropping {len(filter_drops)} {filter_type} specified in manual_drops")
    assert filter_type in qc_drops
    qc_drops[filter_type].update({w: "manual_drop" for w in filter_drops})
    if filter_type == "barcode_wells":
        counts_qc_1 = counts_qc_1[
            ~counts_qc_1.assign(
                barcode_well=lambda x: x.apply(
                    lambda r: (r["barcode"], r["well"]), axis=1
                )
            )["barcode_well"].isin(qc_drops[filter_type])
        ]
    elif filter_type == "barcode_serum_replicates":
        counts_qc_1 = counts_qc_1[
            ~counts_qc_1.assign(
                barcode_serum_replicate=lambda x: x.apply(
                    lambda r: (r["barcode"], r["serum_replicate"]), axis=1
                )
            )["barcode_serum_replicate"].isin(qc_drops[filter_type])
        ]
    elif filter_type == "wells":
        counts_qc_1 = counts_qc_1[~counts_qc_1["well"].isin(qc_drops[filter_type])]
    elif filter_type == "barcodes":
        counts_qc_1 = counts_qc_1[~counts_qc_1["barcode"].isin(qc_drops[filter_type])]
    elif filter_type == "serum_replicates":
        counts_qc_1 = counts_qc_1[
            ~counts_qc_1["serum_replicate"].isin(qc_drops[filter_type])
        ]
    else:
        assert filter_type in set(counts_qc_1.columns)
        counts_qc_1 = counts_qc_1[~counts_qc_1[filter_type].isin(qc_drops[filter_type])]

if not manual_drops:
    report.md("No drops specified in manual_drops")

report.md("""
    ## Collapse the barcodes for each strain
    If `collapse_strain_barcodes` is set in the configuration, sum the counts of all of a
    strain's barcodes, so that everything below is per strain rather than per barcode.
    """)

if collapse_strain_barcodes:
    # Sum the counts of a strain's viral barcodes in each sample, naming the
    # collapsed barcode for the strain so that it still identifies the strain
    # uniquely in the QC and curve fits below. The neut-standard barcodes have no
    # strain and their counts are only ever used summed over each well, so they are
    # kept as they are.
    viral = counts_qc_1.query("not neut_standard")
    # the counts columns are summed over a strain's barcodes, while all of the other
    # columns describe the sample and so are just carried through
    count_cols = ["count", "fraction_all_valid_and_invalid_counts"]
    sample_info = viral[
        [c for c in counts_qc_1.columns if c not in {"barcode", "strain", *count_cols}]
    ].drop_duplicates()
    assert len(sample_info) == sample_info["sample"].nunique()

    counts_collapsed = pd.concat(
        [
            viral.groupby(["sample", "strain"], as_index=False)[count_cols]
            .sum(min_count=1)
            .assign(barcode=lambda x: x["strain"])
            .merge(sample_info, on="sample", validate="many_to_one"),
            counts_qc_1.query("neut_standard"),
        ],
        ignore_index=True,
    )[counts_qc_1.columns]

    barcodes_per_strain = (
        viral.groupby("strain", as_index=False)
        .aggregate(n_barcodes=pd.NamedAgg("barcode", "nunique"))
        .sort_values(["n_barcodes", "strain"], ascending=[False, True])
        .set_index("strain")  # so the display has no meaningless integer index
    )
    report.md(
        f"Collapsed the {viral['barcode'].nunique()} viral barcodes into "
        f"{len(barcodes_per_strain)} strains, with "
        f"{barcodes_per_strain['n_barcodes'].min()} to "
        f"{barcodes_per_strain['n_barcodes'].max()} barcodes per strain.\n\n"
        "Everything below is therefore per strain even where it is labeled as "
        "being per barcode, since each strain now has a single barcode named "
        "for the strain. In particular, any QC threshold on the counts or the "
        "fraction of counts for a barcode now applies to the summed counts of "
        "all of a strain's barcodes."
    )
    report.table(barcodes_per_strain)
else:
    counts_collapsed = counts_qc_1
    report.md(
        "Not collapsing the barcodes for each strain, as "
        "`collapse_strain_barcodes` is not set in the configuration."
    )

report.md("## Average counts per barcode in each well")
report.md("""
    Plot average counts per barcode.
    If a sample has inadequate barcode counts, it may not have good enough statistics for
    accurate analysis, and a QC-threshold is applied:
    """)

# Compute average barcode counts per well
avg_barcode_counts = (
    counts_collapsed.groupby(
        ["well", "serum_replicate", "sample_well"],
        dropna=False,
        as_index=False,
    )
    .aggregate(avg_count=pd.NamedAgg("count", "mean"))
    .assign(
        fails_qc=lambda x: (
            x["avg_count"] < qc_thresholds["avg_barcode_counts_per_well"]
        )
    )
)

avg_barcode_counts_chart = (
    alt.Chart(avg_barcode_counts)
    .add_params(serum_selection)
    .transform_filter(serum_selection)
    .encode(
        alt.X(
            "avg_count",
            title="average barcode counts per well",
            scale=alt.Scale(nice=False, padding=3),
        ),
        alt.Y("sample_well", sort=sample_wells),
        alt.Color(
            "fails_qc",
            title=(
                "fails qc_thresholds['avg_barcode_counts_per_well']="
                f"{qc_thresholds['avg_barcode_counts_per_well']!r}"
            ),
            legend=alt.Legend(titleLimit=500),
        ),
        tooltip=[
            (
                alt.Tooltip(c, format=".3g")
                if avg_barcode_counts[c].dtype == float
                else c
            )
            for c in avg_barcode_counts.columns
        ],
    )
    .mark_bar(height={"band": 0.85})
    .properties(
        height=alt.Step(10),
        width=250,
        title=f"Average barcode counts per well for {plate}",
    )
    .configure_axis(grid=False)
)

report.chart(avg_barcode_counts_chart)

# Drop wells failing QC
avg_barcode_counts_per_well_drops = list(avg_barcode_counts.query("fails_qc")["well"])
report.md(
    f"Dropping {len(avg_barcode_counts_per_well_drops)} wells for failing "
    f"`qc_thresholds['avg_barcode_counts_per_well']="
    f"{qc_thresholds['avg_barcode_counts_per_well']!r}`: "
    f"{avg_barcode_counts_per_well_drops}"
)
qc_drops["wells"].update(
    {w: "avg_barcode_counts_per_well" for w in avg_barcode_counts_per_well_drops}
)
counts_qc_2 = counts_collapsed[~counts_collapsed["well"].isin(qc_drops["wells"])]

report.md("""
    ## Fraction of counts from neutralization standard
    Determine the fraction of counts from the neutralization standard in each sample, and
    make sure this fraction passess the QC threshold.
    """)

# Compute neutralization standard fractions per well
neut_standard_fracs = (
    counts_qc_2.assign(
        neut_standard_count=lambda x: x["count"] * x["neut_standard"].astype(int)
    )
    .groupby(
        ["well", "serum_replicate", "sample_well"],
        dropna=False,
        as_index=False,
    )
    .aggregate(
        total_count=pd.NamedAgg("count", "sum"),
        neut_standard_count=pd.NamedAgg("neut_standard_count", "sum"),
    )
    .assign(
        neut_standard_frac=lambda x: x["neut_standard_count"] / x["total_count"],
        fails_qc=lambda x: (
            x["neut_standard_frac"] < qc_thresholds["min_neut_standard_frac_per_well"]
        ),
    )
)

neut_standard_fracs_chart = (
    alt.Chart(neut_standard_fracs)
    .add_params(serum_selection)
    .transform_filter(serum_selection)
    .encode(
        alt.X(
            "neut_standard_frac",
            title="frac counts from neutralization standard per well",
            scale=alt.Scale(nice=False, padding=3),
        ),
        alt.Y("sample_well", sort=sample_wells),
        alt.Color(
            "fails_qc",
            title=(
                "fails qc_thresholds['min_neut_standard_frac_per_well']="
                f"{qc_thresholds['min_neut_standard_frac_per_well']!r}"
            ),
            legend=alt.Legend(titleLimit=500),
        ),
        tooltip=[
            (
                alt.Tooltip(c, format=".3g")
                if neut_standard_fracs[c].dtype == float
                else c
            )
            for c in neut_standard_fracs.columns
        ],
    )
    .mark_bar(height={"band": 0.85})
    .properties(
        height=alt.Step(10),
        width=250,
        title=f"Neutralization-standard fracs per well for {plate}",
    )
    .configure_axis(grid=False)
    .configure_legend(titleLimit=1000)
)

report.chart(neut_standard_fracs_chart)

# drop wells failing QC
min_neut_standard_frac_per_well_drops = list(
    neut_standard_fracs.query("fails_qc")["well"]
)
report.md(
    f"Dropping {len(min_neut_standard_frac_per_well_drops)} wells for failing "
    f"`qc_thresholds['min_neut_standard_frac_per_well']="
    f"{qc_thresholds['min_neut_standard_frac_per_well']!r}`: "
    f"{min_neut_standard_frac_per_well_drops}"
)
qc_drops["wells"].update(
    {
        w: "min_neut_standard_frac_per_well"
        for w in min_neut_standard_frac_per_well_drops
    }
)
counts_qc_3 = counts_qc_2[~counts_qc_2["well"].isin(qc_drops["wells"])]

report.md("""
    ## Consistency and minimum fractions for barcodes
    We examine the fraction of counts attributable to each barcode. We do this splitting
    the data two ways:

     1. Looking at all viral (but not neut-standard) barcodes only for the no-serum
        samples (wells).

     2. Looking at just the neut-standard barcodes for all samples (wells).

    The reason is that if the experiment is set up perfectly, these fractions should be
    the same across all samples for each barcode.
    (We do not expect viral barcodes to have consistent fractions across no-serum samples
    as they will be neutralized differently depending on strain).

    We plot these fractions in interactive plots (you can mouseover points and zoom) so
    you can identify barcodes that fail the expected consistency QC thresholds.

    We also make sure the barcodes meet specified QC minimum thresholds for all samples,
    and flag any that do not.
    """)

# Create barcode selection parameter for mouseover highlighting
barcode_selection = alt.selection_point(
    fields=["barcode"],
    on="mouseover",
    empty=False,
)

# Analyze barcode evenness for two groups: neut-standard barcodes and viral barcodes
counts_qc_4 = counts_qc_3.copy()

for is_neut_standard, df in counts_qc_3.groupby("neut_standard"):
    # Determine which QC to apply based on barcode type
    if is_neut_standard:
        report.md(
            f"{'=' * 89}\n\nAnalyzing neut-standard barcodes from all samples (wells)"
        )
        qc_name = "per_neut_standard_barcode_filters"
    else:
        report.md(f"{'=' * 89}\n\nAnalyzing all barcodes from no-serum samples (wells)")
        qc_name = "no_serum_per_viral_barcode_filters"
        df = df.query("serum == 'none'")

    # Compute barcode fractions and fold changes
    df = df.assign(
        sample_counts=lambda x: x.groupby("sample")["count"].transform("sum"),
        count_frac=lambda x: x["count"] / x["sample_counts"],
        median_count_frac=lambda x: x.groupby("barcode")["count_frac"].transform(
            "median"
        ),
        fold_change_from_median=lambda x: numpy.where(
            x["count_frac"] > x["median_count_frac"],
            x["count_frac"] / x["median_count_frac"],
            x["median_count_frac"] / x["count_frac"],
        ),
    )[
        ["barcode", "count", "sample_well", "count_frac", "fold_change_from_median"]
        + ([] if is_neut_standard else ["strain"])
    ]

    # Apply QC thresholds
    qc = qc_thresholds[qc_name]
    report.md(f"Apply QC `{qc_name}`: `{qc}`")

    # the lambdas below take `qc` as a keyword argument so that they bind the
    # current loop iteration's `qc` rather than closing over the loop variable
    fails_qc = (
        df.assign(
            fails_qc=lambda x, qc=qc: ~(
                (x["count_frac"] >= qc["min_frac"])
                & (x["fold_change_from_median"] <= qc["max_fold_change"])
            )
        )
        .groupby("barcode", as_index=False)
        .aggregate(n_wells_fail_qc=pd.NamedAgg("fails_qc", "sum"))
        .assign(fails_qc=lambda x, qc=qc: x["n_wells_fail_qc"] >= qc["max_wells"])[
            ["barcode", "fails_qc"]
        ]
    )

    df = df.merge(fails_qc, on="barcode", validate="many_to_one")

    evenness_chart = (
        alt.Chart(df)
        .add_params(barcode_selection)
        .encode(
            alt.X(
                "count_frac",
                title=(
                    "barcode's fraction of neut standard counts"
                    if is_neut_standard
                    else "barcode's fraction of non-neut standard counts"
                ),
                scale=alt.Scale(nice=False, padding=5),
            ),
            alt.Y("sample_well", sort=sample_wells),
            alt.Fill(
                "fails_qc",
                title=f"fails {qc_name}",
                legend=alt.Legend(titleLimit=500),
            ),
            strokeWidth=alt.condition(barcode_selection, alt.value(2), alt.value(0)),
            size=alt.condition(barcode_selection, alt.value(60), alt.value(35)),
            tooltip=[
                alt.Tooltip(c, format=".2g") if df[c].dtype == float else c
                for c in df.columns
            ],
        )
        .mark_circle(fillOpacity=0.45, stroke="black", strokeOpacity=1)
        .properties(
            height=alt.Step(10),
            width=300,
            title=alt.TitleParams(
                (
                    f"{plate} all samples, neut-standard barcodes"
                    if is_neut_standard
                    else f"{plate} no-serum samples, all barcodes"
                ),
                subtitle="x-axis is zoomable (use mouse scroll/pan)",
            ),
        )
        .configure_axis(grid=False)
        .configure_legend(titleLimit=1000)
        .interactive()
    )

    report.chart(evenness_chart)

    # Drop barcodes failing QC
    barcode_drops = list(fails_qc.query("fails_qc")["barcode"])
    report.md(
        f"Dropping {len(barcode_drops)} barcodes for failing `qc={qc!r}`: {barcode_drops}"
    )
    qc_drops["barcodes"].update({bc: qc_name for bc in barcode_drops})
    counts_qc_4 = counts_qc_4[~counts_qc_4["barcode"].isin(qc_drops["barcodes"])]

report.md(r"""
    ## Compute fraction infectivity

    The fraction infectivity for viral barcode $v_b$ in sample $s$ is computed as:

    $$F_{v_b,s} = \frac{c_{v_b,s} / \left(\sum_{n_b} c_{n_b,s}\right)}{{\rm median}\_{s_0}\left[ c_{v_b,s_0} / \left(\sum_{n_b} c_{n_b,s_0}\right)\right]}$$

    where

     - $c_{v_b,s}$ is the counts of viral barcode $v_b$ in sample $s$.
     - $\sum_{n_b} c_{n_b,s}$ is the sum of the counts for all neutralization standard
       barcodes $n_b$ for sample $s$.
     - $c_{v_b,s_0}$ is the counts of viral barcode $v_b$ in no-serum sample $s_0$.
     - $\sum_{n_b} c_{n_b,s_0}$ is the sum of the counts for all neutralization standard
       barcodes $n_b$ for no-serum sample $s_0$.
     - ${\rm median}\_{s_0}\left[ c_{v_b,s_0} / \left(\sum_{n_b} c_{n_b,s_0}\right)\right]$
       is the median taken across all no-serum samples of the counts of viral barcode
       $v_b$ versus the total counts for all neutralization standard barcodes.

    First, compute the total neutralization-standard counts for each sample (well).
    Plot these, and drop any wells that do not meet the QC threshold.
    """)

# Compute neutralization standard counts per well
neut_standard_counts = (
    counts_qc_4.query("neut_standard")
    .groupby(
        [
            "well",
            "serum_replicate",
            "sample_well",
            dilution_factor_or_concentration,
        ],
        dropna=False,
        as_index=False,
    )
    .aggregate(neut_standard_count=pd.NamedAgg("count", "sum"))
    .assign(
        fails_qc=lambda x: (
            x["neut_standard_count"] < qc_thresholds["min_neut_standard_count_per_well"]
        )
    )
)

neut_standard_counts_chart = (
    alt.Chart(neut_standard_counts)
    .add_params(serum_selection)
    .transform_filter(serum_selection)
    .encode(
        alt.X(
            "neut_standard_count",
            title="counts from neutralization standard",
            scale=alt.Scale(nice=False, padding=3),
        ),
        alt.Y("sample_well", sort=sample_wells),
        alt.Color(
            "fails_qc",
            title=(
                "fails qc_thresholds['min_neut_standard_count_per_well']="
                f"{qc_thresholds['min_neut_standard_count_per_well']!r}"
            ),
            legend=alt.Legend(titleLimit=500),
        ),
        tooltip=[
            (
                alt.Tooltip(c, format=".3g")
                if neut_standard_counts[c].dtype == float
                else c
            )
            for c in neut_standard_counts.columns
        ],
    )
    .mark_bar(height={"band": 0.85})
    .properties(
        height=alt.Step(10),
        width=250,
        title=f"Neutralization-standard counts for {plate}",
    )
    .configure_axis(grid=False)
    .configure_legend(titleLimit=1000)
)

report.chart(neut_standard_counts_chart)

# drop wells failing QC
min_neut_standard_count_per_well_drops = list(
    neut_standard_counts.query("fails_qc")["well"]
)
report.md(
    f"Dropping {len(min_neut_standard_count_per_well_drops)} wells for failing "
    f"`qc_thresholds['min_neut_standard_count_per_well']="
    f"{qc_thresholds['min_neut_standard_count_per_well']!r}`: "
    f"{min_neut_standard_count_per_well_drops}"
)
qc_drops["wells"].update(
    {
        w: "min_neut_standard_count_per_well"
        for w in min_neut_standard_count_per_well_drops
    }
)
neut_standard_counts_1 = neut_standard_counts[
    ~neut_standard_counts["well"].isin(qc_drops["wells"])
]
counts_qc_5 = counts_qc_4[~counts_qc_4["well"].isin(qc_drops["wells"])]

report.md("""
    Compute and plot the no-serum sample viral barcode counts and check if they pass the
    QC filters.
    """)

# Compute no-serum viral barcode counts with QC
no_serum_counts = (
    counts_qc_5.query("serum == 'none'")
    .query("not neut_standard")
    .merge(neut_standard_counts_1, validate="many_to_one")[
        ["barcode", "strain", "well", "sample_well", "count", "neut_standard_count"]
    ]
    .assign(
        fails_qc=lambda x: (
            x["count"] <= qc_thresholds["min_no_serum_count_per_viral_barcode_well"]
        )
    )
)

# Create strain selection dropdown
strains = sorted(no_serum_counts["strain"].unique())
strain_selection_dropdown = alt.selection_point(
    fields=["strain"],
    bind=alt.binding_select(
        options=[None] + strains,
        labels=["all"] + strains,
        name="virus strain",
    ),
)

# Prepare data for plotting
no_serum_counts_plot_df = no_serum_counts.drop(columns=["well", "neut_standard_count"])

no_serum_counts_chart = (
    alt.Chart(no_serum_counts_plot_df)
    .add_params(barcode_selection, strain_selection_dropdown)
    .transform_filter(strain_selection_dropdown)
    .encode(
        alt.X(
            "count",
            title="viral barcode count",
            scale=alt.Scale(nice=False, padding=5),
        ),
        alt.Y("sample_well", sort=sample_wells),
        alt.Fill(
            "fails_qc",
            title=(
                "fails qc_thresholds['min_no_serum_count_per_viral_barcode_well']="
                f"{qc_thresholds['min_no_serum_count_per_viral_barcode_well']!r}"
            ),
            legend=alt.Legend(titleLimit=500),
        ),
        strokeWidth=alt.condition(barcode_selection, alt.value(2), alt.value(0)),
        size=alt.condition(barcode_selection, alt.value(60), alt.value(35)),
        tooltip=no_serum_counts_plot_df.columns.tolist(),
    )
    .mark_circle(fillOpacity=0.6, stroke="black", strokeOpacity=1)
    .properties(
        height=alt.Step(10),
        width=400,
        title=f"{plate} viral barcode counts in no-serum samples",
    )
    .configure_axis(grid=False)
    .configure_legend(titleLimit=1000)
    .interactive()
)

report.chart(no_serum_counts_chart)

# drop barcode / wells failing QC
min_no_serum_count_per_viral_barcode_well_drops = list(
    no_serum_counts.query("fails_qc")[["barcode", "well"]].itertuples(
        index=False, name=None
    )
)
report.md(
    f"Dropping {len(min_no_serum_count_per_viral_barcode_well_drops)} barcode-wells for "
    f"failing `qc_thresholds['min_no_serum_count_per_viral_barcode_well']="
    f"{qc_thresholds['min_no_serum_count_per_viral_barcode_well']!r}`: "
    f"{min_no_serum_count_per_viral_barcode_well_drops}"
)
qc_drops["barcode_wells"].update(
    {
        w: "min_no_serum_count_per_viral_barcode_well"
        for w in min_no_serum_count_per_viral_barcode_well_drops
    }
)
no_serum_counts_1 = no_serum_counts[
    ~no_serum_counts.assign(
        barcode_well=lambda x: x.apply(lambda r: (r["barcode"], r["well"]), axis=1)
    )["barcode_well"].isin(qc_drops["barcode_wells"])
]
counts_qc_6 = counts_qc_5[
    ~counts_qc_5.assign(
        barcode_well=lambda x: x.apply(lambda r: (r["barcode"], r["well"]), axis=1)
    )["barcode_well"].isin(qc_drops["barcode_wells"])
]

report.md("""
    Compute and plot the median ratio of viral barcode count to neut standard counts
    across no-serum samples.
    If library composition is equal, all of these values should be similar:
    """)

# Compute median ratio of viral barcode to neut standard counts
median_no_serum_ratio = (
    no_serum_counts_1.assign(ratio=lambda x: x["count"] / x["neut_standard_count"])
    .groupby(["barcode", "strain"], as_index=False)
    .aggregate(median_no_serum_ratio=pd.NamedAgg("ratio", "median"))
)

# Create strain selection for mouseover
strain_selection = alt.selection_point(
    fields=["strain"],
    on="mouseover",
    empty=False,
)

median_no_serum_ratio_chart = (
    alt.Chart(median_no_serum_ratio)
    .add_params(strain_selection)
    .encode(
        alt.X(
            "median_no_serum_ratio",
            title="median ratio of counts",
            scale=alt.Scale(nice=False, padding=5),
        ),
        alt.Y(
            "barcode",
            sort=alt.SortField("median_no_serum_ratio", order="descending"),
            axis=alt.Axis(labelFontSize=5),
        ),
        color=alt.condition(
            strain_selection,
            alt.value("orange"),
            alt.value("gray"),
        ),
        tooltip=[
            (
                alt.Tooltip(c, format=".3g")
                if median_no_serum_ratio[c].dtype == float
                else c
            )
            for c in median_no_serum_ratio.columns
        ],
    )
    .mark_bar(height={"band": 0.85})
    .properties(
        height=alt.Step(5),
        width=250,
        title=f"{plate} no-serum median ratio viral barcode to neut-standard barcode",
    )
    .configure_axis(grid=False)
    .configure_legend(titleLimit=1000)
)

report.chart(median_no_serum_ratio_chart)

report.md("""
    Compute and plot the actual fraction infectivities.
    We compute both the raw fraction infectivities and the ones with the ceiling applied:
    """)

frac_infectivity = (
    counts_qc_6.query("not neut_standard")
    .query("serum != 'none'")
    .merge(median_no_serum_ratio, validate="many_to_one")
    .merge(neut_standard_counts_1, validate="many_to_one")
    .assign(
        frac_infectivity_raw=lambda x: (
            x["count"] / x["neut_standard_count"] / x["median_no_serum_ratio"]
        ),
        frac_infectivity_ceiling=lambda x: (
            x["frac_infectivity_raw"].clip(
                upper=curvefit_params["frac_infectivity_ceiling"]
            )
        ),
        # concentration used for curve fitting: reciprocal of the dilution factor,
        # or the specified concentration used directly
        concentration=lambda x: (
            1 / x["dilution_factor"]
            if dilution_factor_or_concentration == "dilution_factor"
            else x["concentration"]
        ),
        plate_barcode=lambda x: x["plate_replicate"] + "-" + x["barcode"],
    )[
        [
            "barcode",
            "plate_barcode",
            "well",
            "strain",
            "serum",
            "serum_replicate",
        ]
        + (
            ["dilution_factor", "concentration"]
            if dilution_factor_or_concentration == "dilution_factor"
            else ["concentration"]
        )
        + [
            "frac_infectivity_raw",
            "frac_infectivity_ceiling",
        ]
    ]
)
assert len(
    frac_infectivity.groupby(
        ["serum", "plate_barcode", dilution_factor_or_concentration]
    )
) == len(frac_infectivity)
assert frac_infectivity[dilution_factor_or_concentration].notnull().all()
assert frac_infectivity["frac_infectivity_raw"].notnull().all()
assert frac_infectivity["frac_infectivity_ceiling"].notnull().all()

frac_infectivity_cols = {
    "frac_infectivity_raw": "raw fraction infectivity",
    "frac_infectivity_ceiling": (
        f"fraction infectivity with ceiling at "
        f"{curvefit_params['frac_infectivity_ceiling']}"
    ),
}

frac_infectivity_chart_df = frac_infectivity.assign(
    fails_qc=lambda x: (
        x["frac_infectivity_raw"]
        > qc_thresholds["max_frac_infectivity_per_viral_barcode_well"]
    ),
)[
    [
        "barcode",
        "strain",
        "well",
        "serum_replicate",
        dilution_factor_or_concentration,
        "fails_qc",
        *list(frac_infectivity_cols),
    ]
].rename(
    columns=frac_infectivity_cols
)

# some manipulations to shrink data frame plotted by altair below by putting
# them in smaller data frames that are used via transform_lookup
barcode_lookup_df = frac_infectivity[["barcode", "strain"]].drop_duplicates()
assert len(barcode_lookup_df) == barcode_lookup_df["barcode"].nunique()
well_lookup_df = frac_infectivity[
    ["well", "serum_replicate", dilution_factor_or_concentration]
].drop_duplicates()
assert len(well_lookup_df) == well_lookup_df["well"].nunique()

frac_infectivity_chart_df = frac_infectivity_chart_df.drop(
    columns=["strain", "serum_replicate", dilution_factor_or_concentration]
)

frac_infectivity_chart = (
    alt.Chart(frac_infectivity_chart_df)
    .transform_lookup(
        lookup="barcode",
        from_=alt.LookupData(barcode_lookup_df, key="barcode", fields=["strain"]),
    )
    .transform_lookup(
        lookup="well",
        from_=alt.LookupData(
            well_lookup_df,
            key="well",
            fields=["serum_replicate", dilution_factor_or_concentration],
        ),
    )
    .transform_fold(
        frac_infectivity_cols.values(), ["ceiling_applied", "frac_infectivity"]
    )
    .add_params(strain_selection_dropdown, barcode_selection)
    .transform_filter(strain_selection_dropdown)
    .encode(
        alt.X(
            f"{dilution_factor_or_concentration}:Q",
            title=x_label,
            scale=alt.Scale(nice=False, padding=5, type="log"),
        ),
        alt.Y(
            "frac_infectivity:Q",
            title="fraction infectivity",
            scale=alt.Scale(nice=False, padding=5),
        ),
        alt.Column(
            "ceiling_applied:N",
            sort="descending",
            title=None,
            header=alt.Header(labelFontSize=13, labelFontStyle="bold", labelPadding=2),
        ),
        alt.Row(
            "serum_replicate:N",
            title=None,
            spacing=3,
            header=alt.Header(labelFontSize=13, labelFontStyle="bold"),
        ),
        alt.Detail("barcode"),
        alt.Shape(
            "fails_qc",
            title=(
                "fails "
                f"{qc_thresholds['max_frac_infectivity_per_viral_barcode_well']=}"
            ),
            legend=alt.Legend(titleLimit=500, orient="bottom"),
        ),
        color=alt.condition(
            barcode_selection, alt.value("black"), alt.value("MediumBlue")
        ),
        strokeWidth=alt.condition(barcode_selection, alt.value(3), alt.value(1)),
        opacity=alt.condition(barcode_selection, alt.value(1), alt.value(0.25)),
        tooltip=[
            (
                alt.Tooltip(c, format=".3g")
                if frac_infectivity_chart_df[c].dtype == float
                else c
            )
            for c in frac_infectivity_chart_df.columns
        ]
        + [
            alt.Tooltip("strain:N"),
            alt.Tooltip("serum_replicate:N"),
            alt.Tooltip(f"{dilution_factor_or_concentration}:Q", title=x_label),
        ],
    )
    .mark_line(point=True)
    .properties(
        height=150,
        width=250,
        title=f"Fraction infectivities for {plate}",
    )
    .interactive(bind_x=False)
    .configure_axis(grid=False)
    .configure_legend(titleLimit=1000)
    .configure_point(size=50)
    .resolve_scale(x="independent", y="independent")
)

report.chart(frac_infectivity_chart)

# drop barcode / wells failing QC
max_frac_infectivity_per_viral_barcode_well_drops = list(
    frac_infectivity_chart_df.query("fails_qc")[["barcode", "well"]]
    .drop_duplicates()
    .itertuples(index=False, name=None)
)
report.md(
    f"Dropping {len(max_frac_infectivity_per_viral_barcode_well_drops)} barcode-wells "
    f"for failing `qc_thresholds['max_frac_infectivity_per_viral_barcode_well']="
    f"{qc_thresholds['max_frac_infectivity_per_viral_barcode_well']!r}`: "
    f"{max_frac_infectivity_per_viral_barcode_well_drops}"
)
qc_drops["barcode_wells"].update(
    {
        w: "max_frac_infectivity_per_viral_barcode_well"
        for w in max_frac_infectivity_per_viral_barcode_well_drops
    }
)
frac_infectivity_1 = frac_infectivity[
    ~frac_infectivity.assign(
        barcode_well=lambda x: x.apply(lambda r: (r["barcode"], r["well"]), axis=1)
    )["barcode_well"].isin(qc_drops["barcode_wells"])
]

report.md("Check how many dilutions we have per barcode / serum-replicate:")

# Count number of dilutions (or concentrations) per barcode/serum-replicate
n_dilutions = (
    frac_infectivity_1.groupby(
        ["serum_replicate", "strain", "barcode"],
        as_index=False,
    )
    .aggregate(
        **{
            "number of dilutions": pd.NamedAgg(
                dilution_factor_or_concentration, "nunique"
            )
        }
    )
    .assign(
        fails_qc=lambda x: (
            x["number of dilutions"]
            < qc_thresholds["min_dilutions_per_barcode_serum_replicate"]
        )
    )
)

n_dilutions_chart = (
    alt.Chart(n_dilutions)
    .add_params(barcode_selection)
    .encode(
        alt.X(
            "number of dilutions",
            scale=alt.Scale(nice=False, padding=4),
        ),
        alt.Y("strain", title=None),
        alt.Column(
            "serum_replicate",
            title=None,
            header=alt.Header(
                labelFontSize=12,
                labelFontStyle="bold",
                labelPadding=0,
            ),
        ),
        alt.Fill(
            "fails_qc",
            title=(
                "fails qc_thresholds['min_dilutions_per_barcode_serum_replicate']="
                f"{qc_thresholds['min_dilutions_per_barcode_serum_replicate']!r}"
            ),
            legend=alt.Legend(titleLimit=500, orient="bottom"),
        ),
        strokeWidth=alt.condition(barcode_selection, alt.value(2), alt.value(0)),
        size=alt.condition(barcode_selection, alt.value(55), alt.value(35)),
        tooltip=[
            alt.Tooltip(c, format=".3g") if n_dilutions[c].dtype == float else c
            for c in n_dilutions.columns
        ],
    )
    .mark_circle(stroke="black", strokeOpacity=1, fillOpacity=0.45)
    .properties(
        height=alt.Step(10),
        width=120,
        title=alt.TitleParams(
            "number of dilutions for each barcode for each serum-replicate",
            dy=-2,
        ),
    )
)

report.chart(n_dilutions_chart)

# Drop barcode/serum-replicates failing QC
min_dilutions_per_barcode_serum_replicate_drops = list(
    n_dilutions.query("fails_qc")[["barcode", "serum_replicate"]].itertuples(
        index=False, name=None
    )
)
report.md(
    f"Dropping {len(min_dilutions_per_barcode_serum_replicate_drops)} "
    f"barcode/serum-replicates for failing "
    f"`qc_thresholds['min_dilutions_per_barcode_serum_replicate']="
    f"{qc_thresholds['min_dilutions_per_barcode_serum_replicate']!r}`: "
    f"{min_dilutions_per_barcode_serum_replicate_drops}"
)
qc_drops["barcode_serum_replicates"].update(
    {
        w: "min_dilutions_per_barcode_serum_replicate"
        for w in min_dilutions_per_barcode_serum_replicate_drops
    }
)
frac_infectivity_2 = frac_infectivity_1[
    ~frac_infectivity_1.assign(
        barcode_serum_replicate=lambda x: x.apply(
            lambda r: (r["barcode"], r["serum_replicate"]), axis=1
        )
    )["barcode_serum_replicate"].isin(qc_drops["barcode_serum_replicates"])
]

report.md("""
    ## Fit neutralization curves without applying QC to curves
    First fit curves to all serum replicates, then we will apply QC on the curve fits.
    Note that the fitting is done to the fraction infectivities **with** the ceiling.
    """)

fits_noqc = neutcurve.CurveFits(
    frac_infectivity_2.rename(
        columns={
            "frac_infectivity_ceiling": "fraction infectivity",
        }
    ),
    conc_col="concentration",
    fracinf_col="fraction infectivity",
    virus_col="strain",
    serum_col="serum_replicate",
    replicate_col="barcode",
    fixtop=curvefit_params["fixtop"],
    fixbottom=curvefit_params["fixbottom"],
    fixslope=curvefit_params["fixslope"],
)

report.md("""
    Determine which fits fail the curve fitting QC, and plot them.
    Note the plot indicates as failing QC any barcode / serum-replicate that fails, even
    if we are also specified to ignore the QC for that one (so it will not be removed
    later):
    """)

goodness_of_fit = curvefit_qc["goodness_of_fit"]
fit_params_noqc = (
    frac_infectivity_2.groupby(["serum_replicate", "barcode"], as_index=False)
    .aggregate(max_frac_infectivity=pd.NamedAgg("frac_infectivity_ceiling", "max"))
    .merge(
        fits_noqc.fitParams(average_only=False, no_average=True)[
            ["serum", "virus", "replicate", "r2", "rmsd"]
        ].rename(columns={"serum": "serum_replicate", "replicate": "barcode"}),
        validate="one_to_one",
    )
    .assign(
        fails_max_frac_infectivity_at_least=lambda x: (
            x["max_frac_infectivity"] < curvefit_qc["max_frac_infectivity_at_least"]
        ),
        fails_goodness_of_fit=lambda x: (
            (x["r2"] < goodness_of_fit["min_R2"])
            & (x["rmsd"] > goodness_of_fit["max_RMSD"])
        ),
        fails_qc=lambda x: (
            x["fails_max_frac_infectivity_at_least"] | x["fails_goodness_of_fit"]
        ),
        ignore_qc=lambda x: x.apply(
            lambda r: (
                r["serum_replicate"]
                in curvefit_qc["serum_replicates_ignore_curvefit_qc"]
                or (r["barcode"], r["serum_replicate"])
                in curvefit_qc["barcode_serum_replicates_ignore_curvefit_qc"]
            ),
            axis=1,
        ),
    )
)

report.md(
    f"Plotting barcode / serum-replicates that fail `curvefit_qc={curvefit_qc!r}`"
)

fit_params_noqc_base_chart = alt.Chart(fit_params_noqc).add_params(barcode_selection)
for prop, col in [
    ("max frac infectivity", "max_frac_infectivity"),
    ("curve fit R2", "r2"),
    ("curve fit RMSD", "rmsd"),
]:
    chart = (
        fit_params_noqc_base_chart.encode(
            alt.X(col, title=prop, scale=alt.Scale(nice=False, padding=4)),
            alt.Y("virus", title=None),
            alt.Fill("fails_qc"),
            alt.Column(
                "serum_replicate",
                title=None,
                header=alt.Header(
                    labelFontSize=12, labelFontStyle="bold", labelPadding=0
                ),
            ),
            strokeWidth=alt.condition(barcode_selection, alt.value(2), alt.value(0)),
            size=alt.condition(barcode_selection, alt.value(55), alt.value(35)),
            tooltip=[
                (
                    alt.Tooltip(c, format=".3g")
                    if fit_params_noqc[c].dtype == float
                    else c
                )
                for c in fit_params_noqc.columns
            ],
        )
        .mark_circle(stroke="black", strokeOpacity=1, fillOpacity=0.55)
        .properties(
            height=alt.Step(10),
            width=90,
            title=alt.TitleParams(f"{prop} for each barcode serum-replicate", dy=-2),
        )
    )

    report.chart(chart)

report.md("""
    Now plot curves for all virus vs serum-replicates that have a barcode that fails any
    of the QC.
    In these plots, the suffix on the barcode name in the color key indicates if it
    passed or failed QC:
    """)

barcode_serum_replicates_fail_qc = fit_params_noqc.query("fails_qc").reset_index(
    drop=True
)

if len(barcode_serum_replicates_fail_qc):
    report.md(
        f"Here are barcode / serum-replicates that fail `curvefit_qc={curvefit_qc!r}`"
    )
    report.table(barcode_serum_replicates_fail_qc, index=False)
    report.md(
        "Curves for virus vs serum-replicates with at least one failed barcode.\n\n"
        "Color key labels indicate if barcodes failed or passed QC."
    )
    plots = {}
    ncol = 6
    for iplot, (serum, virus, failed_barcodes) in enumerate(
        barcode_serum_replicates_fail_qc.groupby(
            ["serum_replicate", "virus"], as_index=False
        )
        .aggregate(barcodes=pd.NamedAgg("barcode", list))
        .itertuples(index=False)
    ):
        passed_barcodes = [
            bc
            for bc in fits_noqc.replicates[serum, virus]
            if bc not in failed_barcodes and bc != "average"
        ]
        curvelist = []
        assert len(CBMARKERS) >= len(failed_barcodes + passed_barcodes)
        assert len(CBPALETTE) >= len(failed_barcodes + passed_barcodes)
        for replicate, marker, color in zip(
            failed_barcodes + passed_barcodes, CBMARKERS, CBPALETTE
        ):
            curvelist.append(
                {
                    "serum": serum,
                    "virus": virus,
                    "replicate": replicate,
                    "label": replicate
                    + ("-fail" if replicate in failed_barcodes else "-pass"),
                    "color": color,
                    "marker": marker,
                }
            )
        plots[iplot // ncol, iplot % ncol] = (f"{serum} vs {virus}", curvelist)
    fig_fail_qc, _ = fits_noqc.plotGrid(
        plots,
        attempt_shared_legend=False,
        legendfontsize=8,
        titlesize=9,
        ticksize=10,
        draw_in_bounds=True,
    )
    report.figure(fig_fail_qc, curve_display_method)
else:
    report.md("No serum-replicates fail QC.")

# drop barcode / serum-replicates failing QC
frac_infectivity_3 = frac_infectivity_2.copy()
fit_params_noqc_1 = fit_params_noqc.copy()
for qc_filter in ["max_frac_infectivity_at_least", "goodness_of_fit"]:
    fits_qc_drops = list(
        fit_params_noqc_1.query(f"fails_{qc_filter} and (not ignore_qc)")[
            ["barcode", "serum_replicate"]
        ].itertuples(index=False, name=None)
    )
    report.md(
        f"Dropping {len(fits_qc_drops)} barcode/serum-replicates for failing "
        f"`{qc_filter}={curvefit_qc[qc_filter]}`: {fits_qc_drops}"
    )
    qc_drops["barcode_serum_replicates"].update({w: qc_filter for w in fits_qc_drops})
    frac_infectivity_3 = frac_infectivity_3[
        ~frac_infectivity_3.assign(
            barcode_serum_replicate=lambda x: x.apply(
                lambda r: (r["barcode"], r["serum_replicate"]), axis=1
            )
        )["barcode_serum_replicate"].isin(qc_drops["barcode_serum_replicates"])
    ]
    fit_params_noqc_1 = fit_params_noqc_1[
        ~fit_params_noqc_1.assign(
            barcode_serum_replicate=lambda x: x.apply(
                lambda r: (r["barcode"], r["serum_replicate"]), axis=1
            )
        )["barcode_serum_replicate"].isin(qc_drops["barcode_serum_replicates"])
    ]

report.md("""
    ## Fit neutralization curves after applying QC
    No we re-fit and plot curves after applying all the QC.
    """)

fits_qc = neutcurve.CurveFits(
    frac_infectivity_3.rename(
        columns={
            "frac_infectivity_ceiling": "fraction infectivity",
        }
    ),
    conc_col="concentration",
    fracinf_col="fraction infectivity",
    virus_col="strain",
    serum_col="serum",
    replicate_col="plate_barcode",
    fixtop=curvefit_params["fixtop"],
    fixbottom=curvefit_params["fixbottom"],
    fixslope=curvefit_params["fixslope"],
)
fit_params_qc = fits_qc.fitParams(average_only=False, no_average=True)
assert len(fit_params_qc) <= len(
    fits_noqc.fitParams(average_only=False, no_average=True)
)
report.md(f"Assigning fits for this plate to `{group}`")
fit_params_qc.insert(0, "group", group)
# record the plate so downstream analyses do not have to parse it out of the
# replicate name (which is a `plate_barcode`), where one plate name being a prefix
# of another would make the parse ambiguous
fit_params_qc.insert(1, "plate", plate)

if fits_qc.sera:
    fig_passed_qc, _ = fits_qc.plotReplicates(
        attempt_shared_legend=False,
        legendfontsize=8,
        titlesize=9,
        ticksize=10,
        ncol=6,
        draw_in_bounds=True,
    )
    report.figure(fig_passed_qc, curve_display_method)
else:
    report.md("No sera passed QC.")

report.md("## Save results to files")

report.md(f"Writing fraction infectivities to `{frac_infectivity_csv}`", log=True)
(
    frac_infectivity_3[
        [
            "serum",
            "strain",
            "plate_barcode",
            dilution_factor_or_concentration,
            "frac_infectivity_raw",
            "frac_infectivity_ceiling",
        ]
    ]
    .sort_values(["serum", "plate_barcode", dilution_factor_or_concentration])
    .to_csv(frac_infectivity_csv, index=False, float_format="%.4g")
)
report.md(f"Writing fit parameters to `{fits_csv}`", log=True)
fit_params_qc.drop(columns=["nreplicates", "ic50_str"]).to_csv(
    fits_csv, index=False, float_format="%.4g"
)
report.md(
    f"Pickling neutcurve.CurveFits object for these data to `{fits_pickle}`", log=True
)
with open(fits_pickle, "wb") as f_pickle:
    pickle.dump(fits_qc, f_pickle)
report.md(f"Writing QC drops to `{qc_drops_yaml}`", log=True)


def tup_to_str(x):
    """Join a tuple key into the space-separated string written to the YAML."""
    return " ".join(x) if isinstance(x, tuple) else x


qc_drops_for_yaml = {
    key: {tup_to_str(key2): val2 for key2, val2 in val.items()}
    for key, val in qc_drops.items()
}
with open(qc_drops_yaml, "w") as f_yaml:
    yaml_rt().dump(qc_drops_for_yaml, f_yaml)
report.md("Here are the QC drops:")
report.yaml(qc_drops_for_yaml)

report.write(snakemake.output.html)
print(f"Wrote the report to {snakemake.output.html}")

"""Analyze titers for a serum assigned to a group, aggregating replicates across plates."""

import itertools
import pickle
import sys

import altair as alt
import matplotlib
import matplotlib.pyplot as plt
import neutcurve
import numpy
import pandas as pd
from ruamel import yaml
from seqneut_funcs import (
    get_median_bound,
    narrow_for_altair,
    padded_log_domain,
    pearson_r_log10,
)
from seqneut_report import Report

# `noqa: SIM115` as this log file must stay open for the life of the script
sys.stderr = sys.stdout = open(snakemake.log[0], "w")  # noqa: SIM115

_ = alt.data_transformers.disable_max_rows()

# faster plotting of neut curves
matplotlib.style.use("fast")

pickle_fits = snakemake.input.pickles
fits_csvs = snakemake.input.fits_csvs
per_rep_titers_csv = snakemake.output.per_rep_titers
titers_csv = snakemake.output.titers
curves_pdf = snakemake.output.curves_pdf
output_pickle = snakemake.output.pickle
qc_drops_file = snakemake.output.qc_drops
viral_strain_plot_order = snakemake.params.viral_strain_plot_order
serum_titer_as = snakemake.params.serum_titer_as
qc_thresholds = snakemake.params.qc_thresholds
curve_display_method = snakemake.params.curve_display_method
dilution_factor_or_concentration = snakemake.params.dilution_factor_or_concentration
concentration_units = snakemake.params.concentration_units
serum = snakemake.wildcards.serum
group = snakemake.wildcards.group

report = Report(title="Titers for a serum in a group")

report.md("""
    Analyze titers for a serum assigned to a group, aggregating replicates which may be
    across multiple plates.
    """)
report.md(f"Processing `{group=}`, `{serum=}`")

report.md("""
    ## Get all titers for this plate
    Combine all the pickled `neutcurve.CurveFits` from plates for this serum into a
    single `neutcurve.CurveFits`:
    """)

report.md(
    f"Combining the curve fits for `group={group!r}`, `serum={serum!r}` from "
    f"`pickle_fits={pickle_fits!r}`"
)
fits_to_combine = []
for fname in pickle_fits:
    with open(fname, "rb") as f_combine:
        fits_to_combine.append(pickle.load(f_combine))
fits_noqc = neutcurve.CurveFits.combineCurveFits(fits_to_combine, sera=[serum])

report.md("""
    Get the plate that each replicate was measured on.
    Each plate records this in its `curvefits.csv`, so it does not have to be parsed out
    of the replicate name (a `plate_barcode`), which would be ambiguous when one plate
    name is a prefix of another:
    """)

plate_reps = pd.concat(
    [
        pd.read_csv(f, dtype={"plate": str, "replicate": str})[["plate", "replicate"]]
        for f in fits_csvs
    ],
    ignore_index=True,
).drop_duplicates()
assert (
    len(plate_reps) == plate_reps["replicate"].nunique()
), f"a replicate is assigned to more than one plate:\n{plate_reps}"
plate_of_replicate = plate_reps.set_index("replicate")["plate"].to_dict()
# plates in the order they are listed for this serum, used to order the plate pairs
serum_plates = list(dict.fromkeys(plate_reps["plate"]))
report.md(f"`{serum}` was measured on {len(serum_plates)} plate(s): {serum_plates}")

report.md("Indicate how we are calculating the titer:")

report.md(f"Calculating with `serum_titer_as={serum_titer_as!r}`")
if dilution_factor_or_concentration == "dilution_factor":
    assert serum_titer_as in {"nt50", "midpoint"}, serum_titer_as
else:
    assert serum_titer_as in {"ic50", "midpoint"}, serum_titer_as

report.md("""
    Get all the per-replicate fit params with the titers.
    We also convert the IC50 to NT50, and take inverse of midpoint to get it on same
    scale as NT50s:
    """)

fit_params = fits_noqc.fitParams(average_only=False, no_average=True).assign(
    group=group,
    plate=lambda x: x["replicate"].map(plate_of_replicate),
)
assert fit_params["plate"].notnull().all(), (
    "could not assign a plate to every replicate:\n"
    f"{fit_params[fit_params['plate'].isnull()]}"
)
if dilution_factor_or_concentration == "dilution_factor":
    # titer is a reciprocal serum dilution: convert IC50 to NT50 and take the
    # inverse of the midpoint to put it on the same scale. Taking the reciprocal
    # flips the direction of a bound, so map lower <-> upper.
    per_rep_titers = fit_params.assign(
        nt50=lambda x: 1 / x["ic50"],
        midpoint=lambda x: 1 / x["midpoint_bound"],
        titer=lambda x: (x["midpoint"] if serum_titer_as == "midpoint" else x["nt50"]),
        titer_bound=lambda x: (
            x["midpoint_bound_type"]
            if serum_titer_as == "midpoint"
            else x["ic50_bound"]
        ).map({"lower": "upper", "upper": "lower", "interpolated": "interpolated"}),
        titer_as=serum_titer_as,
        titer_units="reciprocal_dilution",
    )[
        [
            "group",
            "serum",
            "virus",
            "plate",
            "replicate",
            "titer",
            "titer_bound",
            "titer_as",
            "titer_units",
            "nt50",
            "midpoint",
            "top",
            "bottom",
            "slope",
        ]
    ]
else:
    # titer is the fitted concentration directly (e.g. IC50 in `concentration_units`);
    # no reciprocal is taken, so bounds are reported as-is (no lower <-> upper swap).
    per_rep_titers = fit_params.assign(
        midpoint=lambda x: x["midpoint_bound"],
        titer=lambda x: (x["midpoint"] if serum_titer_as == "midpoint" else x["ic50"]),
        titer_bound=lambda x: (
            x["midpoint_bound_type"]
            if serum_titer_as == "midpoint"
            else x["ic50_bound"]
        ),
        titer_as=serum_titer_as,
        titer_units=concentration_units,
    )[
        [
            "group",
            "serum",
            "virus",
            "plate",
            "replicate",
            "titer",
            "titer_bound",
            "titer_as",
            "titer_units",
            "ic50",
            "midpoint",
            "top",
            "bottom",
            "slope",
        ]
    ]
assert per_rep_titers.notnull().all().all()

if len(
    invalid_titer_as := per_rep_titers.query(
        "(titer_as in ['nt50', 'ic50']) and top <= 0.5"
    )
):
    raise ValueError(
        "There are titers computed as nt50/ic50 when curve top <= 0.5:\n"
        f"{invalid_titer_as}"
    )
assert len(per_rep_titers) == per_rep_titers["replicate"].nunique()

# get viruses in the order to plot them
viruses = sorted(per_rep_titers["virus"].unique())
if viral_strain_plot_order is not None:
    if not set(viruses).issubset(viral_strain_plot_order):
        raise ValueError(
            "`viral_strain_plot_order` lacks some viruses with titers:\n"
            + str(set(viruses) - set(viral_strain_plot_order))
        )
    viruses = [v for v in viral_strain_plot_order if v in viruses]
report.md(f"`{serum}` has titers for a total of {len(viruses)} viruses")

report.md("""
    ## Correlate NT50s with midpoints of curves
    Plot the correlation of the NT50s with the midpoint (this is an interactive plot,
    mouse over points for details).
    This plot can help you determine if you made the correct choice of `serum_titer_as`
    when choosing to use the midpoint or NT50 for the titer.
    For titers where they are well correlated it should not matter which you chose.
    But if there are titers far from the correlation line, you should look at those
    measurements and curves to make sure you made the correct choice of calculating the
    titer as the NT50 versus midpoint:
    """)

virus_selection_1 = alt.selection_point(fields=["virus"], on="mouseover", empty=False)
# the per-replicate potency column is 'nt50' (dilution) or 'ic50' (concentration)
potency_col = (
    "nt50" if dilution_factor_or_concentration == "dilution_factor" else "ic50"
)
# `group`, `serum`, and `titer_as` are dropped because they are not in the tooltips
# below, and `titer_units` because it is the same on every row and so is put back as a
# constant by the `transform_calculate`
titer_units = per_rep_titers["titer_units"].unique()
assert len(titer_units) == 1, titer_units
chart_df = narrow_for_altair(
    per_rep_titers, drop=["group", "serum", "titer_as", "titer_units"]
)
midpoint_vs_nt50_chart = (
    alt.Chart(chart_df)
    .add_params(virus_selection_1)
    .transform_calculate(titer_units=f"'{titer_units[0]}'")
    .encode(
        alt.X(potency_col, scale=alt.Scale(type="log", nice=False, padding=8)),
        alt.Y("midpoint", scale=alt.Scale(type="log", nice=False, padding=8)),
        alt.Color("titer_bound"),
        strokeWidth=alt.condition(virus_selection_1, alt.value(3), alt.value(0)),
        size=alt.condition(virus_selection_1, alt.value(100), alt.value(60)),
        tooltip=[
            alt.Tooltip(c, format=".2g") if chart_df[c].dtype == float else c
            for c in chart_df.columns
        ]
        # a calculated field has no `pandas` dtype, so its type is given explicitly
        + [alt.Tooltip("titer_units:N")],
    )
    .mark_circle(stroke="black", fillOpacity=0.45, color="black")
    .properties(
        width=350,
        height=350,
        title=f"{potency_col.upper()} versus midpoint for {group} {serum}",
    )
    .configure_axis(grid=False)
)
report.chart(midpoint_vs_nt50_chart)

report.md("""
    ## Plot median titers and determine if they pass QC
    Get the median titers for each virus across replicates, then add these median titers
    to the per-replicate titers and calculate the fold-change in titer between each
    replicate and its median.
    Finally, for each virus indicate whether it passes the QC:
    """)

report.md(f"Using the following `qc_thresholds={qc_thresholds!r}`")

median_titers_noqc = (
    per_rep_titers.sort_values("titer")  # for getting median bound
    .groupby(["group", "serum", "virus", "titer_as", "titer_units"], as_index=False)
    .aggregate(
        titer=pd.NamedAgg("titer", "median"),
        n_replicates=pd.NamedAgg("replicate", "count"),
        titer_sem=pd.NamedAgg("titer", "sem"),
        titer_bound=pd.NamedAgg("titer_bound", get_median_bound),
    )
)

per_rep_titers_w_fc = (
    per_rep_titers.merge(
        median_titers_noqc[["group", "serum", "virus", "titer"]].rename(
            columns={"titer": "median_titer"}
        ),
        validate="many_to_one",
        on=["group", "serum", "virus"],
    )
    .assign(
        fc_from_median=lambda x: numpy.where(
            x["titer"] > x["median_titer"],
            x["titer"] / x["median_titer"],
            x["median_titer"] / x["titer"],
        ),
    )
    .drop(columns=["group", "serum", "titer_as", "median_titer"])
)

median_titers_noqc = median_titers_noqc.merge(
    per_rep_titers_w_fc.groupby("virus", as_index=False).aggregate(
        max_fc_from_median=pd.NamedAgg("fc_from_median", "max")
    ),
    on="virus",
    validate="one_to_one",
).assign(
    fails_min_reps=lambda x: x["n_replicates"] < qc_thresholds["min_replicates"],
    fails_max_fc=lambda x: (
        x["max_fc_from_median"] >= qc_thresholds["max_fold_change_from_median"]
    ),
    fails_qc=lambda x: x["fails_min_reps"] | x["fails_max_fc"],
    fails_qc_reason=lambda x: (
        x.apply(
            lambda r: ", ".join(
                (["min_replicates"] if r["fails_min_reps"] else [])
                + (["max_fold_change_from_median"] if r["fails_max_fc"] else [])
            ),
            axis=1,
        )
    ),
)

# get viruses failing QC in order to plot
viruses_failing_qc = (
    median_titers_noqc.query("fails_qc").set_index("virus")["fails_qc_reason"].to_dict()
)
viruses_failing_qc = {
    v: viruses_failing_qc[v] for v in viruses if v in viruses_failing_qc
}

median_titers_noqc = median_titers_noqc.drop(
    columns=["fails_min_reps", "fails_max_fc", "fails_qc_reason"]
)

per_rep_titers_w_fc = per_rep_titers_w_fc.merge(
    median_titers_noqc[["virus", "fails_qc"]],
    on="virus",
    validate="many_to_one",
)

report.md("""
    Now plot the per-replicate and median titers, indicating any viruses that failed QC.
    Note that potentially some of these titers may still be retained if the viruses in
    question are specified in `viruses_ignore_qc` of `qc_thresholds`.
    """)

virus_selection_2 = alt.selection_point(fields=["virus"], on="mouseover", empty=False)
titer_axis_title = (
    "titer"
    if dilution_factor_or_concentration == "dilution_factor"
    else f"titer ({concentration_units})"
)
# the columns dropped here are the same on every row, and are put back as constants by
# the `transform_calculate` of each layer so that the tooltips still show them
constants = {
    "group": group,
    "serum": serum,
    "titer_as": median_titers_noqc["titer_as"].unique()[0],
    "titer_units": median_titers_noqc["titer_units"].unique()[0],
}
calculate = {col: f"'{val}'" for col, val in constants.items()}
per_rep_df = narrow_for_altair(per_rep_titers_w_fc, drop=["titer_units"])
median_df = narrow_for_altair(median_titers_noqc, drop=list(constants))
per_rep_chart = (
    alt.Chart(per_rep_df)
    .transform_calculate(titer_units=calculate["titer_units"])
    .encode(
        alt.X(
            "titer",
            title=titer_axis_title,
            scale=alt.Scale(nice=False, padding=5, type="log"),
        ),
        alt.Y("virus", sort=viruses),
        alt.Fill(
            "fails_qc",
            title=(
                f"fails qc_thresholds['min_replicates']={qc_thresholds['min_replicates']!r}, "
                f"qc_thresholds['max_fold_change_from_median']="
                f"{qc_thresholds['max_fold_change_from_median']!r}"
            ),
            legend=alt.Legend(titleLimit=500),
        ),
        alt.Shape("titer_bound"),
        strokeWidth=alt.condition(virus_selection_2, alt.value(2), alt.value(0)),
        tooltip=[
            (alt.Tooltip(c, format=".3g") if per_rep_df[c].dtype == float else c)
            for c in per_rep_df
        ]
        # a calculated field has no `pandas` dtype, so its type is given explicitly
        + [alt.Tooltip("titer_units:N")],
    )
    .mark_point(size=35, filled=True, fillOpacity=0.5, strokeOpacity=1, stroke="black")
)
median_chart = (
    alt.Chart(median_df)
    .transform_calculate(**calculate)
    .encode(
        alt.X(
            "titer",
            title=titer_axis_title,
            scale=alt.Scale(nice=False, padding=5, type="log"),
        ),
        alt.Y("virus", sort=viruses),
        alt.Fill("fails_qc"),
        alt.Shape("titer_bound"),
        strokeWidth=alt.condition(virus_selection_2, alt.value(2), alt.value(0.5)),
        tooltip=[
            (alt.Tooltip(c, format=".3g") if median_df[c].dtype == float else c)
            for c in median_df
        ]
        + [alt.Tooltip(f"{c}:N") for c in calculate],
    )
    .mark_point(size=75, filled=True, fillOpacity=0.9, strokeOpacity=1, stroke="black")
)
titer_chart = (
    (per_rep_chart + median_chart)
    .add_params(virus_selection_2)
    .properties(
        height=alt.Step(11),
        width=250,
        title=(
            f"{group} {serum} median (large points) and per-replicate (small points) titers"
        ),
    )
    .configure_axis(grid=False)
)
report.chart(titer_chart)

report.md("""
    ## Plot individual curves for any viruses failing QC
    Plot individual curves for viruses failing QC.
    Note that potentially some of these titers may still be retained if the viruses in
    question are specified in `viruses_ignore_qc` of `qc_thresholds`.
    """)

report.md(
    f"Neutralization curves for the {len(viruses_failing_qc)} viruses failing QC:"
)
if len(viruses_failing_qc):
    fig_failing, _ = fits_noqc.plotReplicates(
        viruses=viruses_failing_qc,
        attempt_shared_legend=False,
        legendfontsize=8,
        titlesize=12,
        ncol=4,
        heightscale=1.2,
        widthscale=1.2,
        subplot_titles="{virus}",
        draw_in_bounds=True,
    )
    _ = fig_failing.suptitle(
        f"neutralization curves for viruses failing QC for {group} {serum}",
        y=1,
        fontsize=18,
        fontweight="bold",
    )
    fig_failing.tight_layout()
    report.figure(fig_failing, curve_display_method)
    plt.close(fig_failing)
else:
    report.md("No curves fail QC")

report.md("""
    ## Get the viruses to drop for QC failures
    Drop any viruses that fail QC and are not specified in `viruses_ignore_qc` of
    `qc_thresholds`.
    """)

viruses_to_drop = {
    v: reason
    for v, reason in viruses_failing_qc.items()
    if v not in qc_thresholds["viruses_ignore_qc"]
}
report.md(f"Dropping {len(viruses_to_drop)} viruses for failing QC:")
report.yaml(viruses_to_drop)
if nkept := (len(viruses_failing_qc) - len(viruses_to_drop)):
    kept_viruses = {
        v: reason
        for v, reason in viruses_failing_qc.items()
        if v in qc_thresholds["viruses_ignore_qc"]
    }
    report.md(
        f"Retaining {nkept} viruses that fail QC because they are in "
        f"`viruses_ignore_qc`: {kept_viruses}"
    )
report.md(f"Writing QC drops to `{qc_drops_file}`", log=True)
with open(qc_drops_file, "w") as f_qc_drops:
    yaml.YAML(typ="rt").dump(viruses_to_drop, f_qc_drops)

report.md("""
    Write the individual per-replicate titers to a file.
    These are all of the replicates that passed the QC applied when their plate was
    processed, and so are not filtered by the per-serum QC just applied above.
    Instead, the `dropped_by_qc` column indicates whether a replicate's virus was dropped
    by that per-serum QC, making this file a superset of the titers written below.
    Note that the per-serum QC drops a virus rather than an individual replicate, so this
    column has the same value for all replicates of a virus:
    """)

report.md(f"Writing per-replicate titers to `{per_rep_titers_csv}`", log=True)
per_rep_titers.assign(
    dropped_by_qc=lambda x: x["virus"].isin(set(viruses_to_drop))
).to_csv(per_rep_titers_csv, index=False, float_format="%.4g")

report.md("""
    ## Plate-to-plate correlation of titers
    If this serum was measured on more than one plate, compare the titers between each
    pair of plates to assess plate-to-plate reproducibility.
    The titer for a strain on a plate is the median over all of that plate's barcodes,
    all of which have already passed the per-plate QC applied when the plate was
    processed.

    All strains are shown, colored by whether they are retained in the final titers or
    were dropped just above by the per-serum QC on replicate-to-replicate variation.
    **Use the `show strains` dropdown at the bottom of the plots to show only the
    retained or only the dropped strains**, and set it back to `all` to show all of them.
    Note that a strain listed in `viruses_ignore_qc` of `qc_thresholds` is colored as
    retained even if it fails the QC, because it is kept in the final titers; mouse over
    a point for the details.
    The Pearson correlation is reported both over all strains and over just the retained
    strains, and the dashed line is `y = x`.
    """)

# the titer for a strain on a plate is the median over that plate's barcodes
plate_titers = (
    per_rep_titers.sort_values("titer")  # for getting median bound
    .groupby(["virus", "plate"], as_index=False)
    .aggregate(
        titer=pd.NamedAgg("titer", "median"),
        n_barcodes=pd.NamedAgg("replicate", "count"),
        titer_bound=pd.NamedAgg("titer_bound", get_median_bound),
    )
)

dropped = set(viruses_to_drop)
pair_viruses = set()
plate_pair_info = []
for plate_x, plate_y in itertools.combinations(serum_plates, 2):
    pair_df = (
        plate_titers[plate_titers["plate"] == plate_x]
        .drop(columns="plate")
        .merge(
            plate_titers[plate_titers["plate"] == plate_y].drop(columns="plate"),
            on="virus",
            suffixes=("_x", "_y"),
            validate="one_to_one",
        )
        .assign(
            censored=lambda x: (x["titer_bound_x"] != "interpolated")
            | (x["titer_bound_y"] != "interpolated"),
        )
    )
    if not len(pair_df):
        continue
    kept_df = pair_df[~pair_df["virus"].isin(dropped)]
    plate_pair_info.append(
        {
            "pair": f"{plate_x} vs {plate_y}",
            "plate_x": plate_x,
            "plate_y": plate_y,
            "n": len(pair_df),
            "r": pearson_r_log10(pair_df),
            "n_kept": len(kept_df),
            "r_kept": pearson_r_log10(kept_df),
            # both axes of the subplot get this domain, so that the y = x line drawn
            # across it is a true diagonal
            "domain": padded_log_domain(pair_df["titer_x"], pair_df["titer_y"]),
            # keep the plotted data frame narrow: each subplot gets just its own pair
            # of plates, so no column names the pair. That costs no extra rows, since
            # `altair` consolidates every frame into the chart's `datasets` and the
            # rows are partitioned among the subplots either way. The plates are named
            # in the subplot's axis titles and tooltips rather than in a column, and
            # the QC status of a strain is determined by `virus` so it goes in the
            # `virus_qc_lookup` below.
            "titers": narrow_for_altair(
                pair_df[
                    [
                        "virus",
                        "titer_x",
                        "titer_y",
                        "n_barcodes_x",
                        "n_barcodes_y",
                        "censored",
                    ]
                ].reset_index(drop=True)
            ),
        }
    )
    pair_viruses.update(pair_df["virus"])

if plate_pair_info:
    # end points of each subplot's y = x line, which are the ends of its axes rather
    # than of its data, so that the line spans the whole subplot
    for info in plate_pair_info:
        info["diagonal"] = pd.DataFrame({"titer": info["domain"]})

    def qc_reason(v):
        """Explain the QC status of a strain."""
        if v in viruses_to_drop:
            return f"fails {viruses_to_drop[v]}"
        elif v in viruses_failing_qc:
            return f"fails {viruses_failing_qc[v]} but retained via `viruses_ignore_qc`"
        else:
            return "passes QC"

    virus_qc_lookup = pd.DataFrame(
        {"virus": [v for v in viruses if v in pair_viruses]}
    ).assign(
        qc_status=lambda x: numpy.where(
            x["virus"].isin(dropped),
            "dropped by per-serum QC",
            "retained in final titers",
        ),
        qc_reason=lambda x: x["virus"].map(qc_reason),
    )
    assert len(virus_qc_lookup) == virus_qc_lookup["virus"].nunique()

    report.md(
        f"Comparing titers between {len(plate_pair_info)} pair(s) of plates for `{serum}`"
    )
else:
    virus_qc_lookup = None
    if len(serum_plates) < 2:
        report.md(
            f"`{serum}` was measured on only one plate (`{serum_plates[0]}`), so there "
            "is no plate-to-plate correlation to show."
        )
    else:
        report.md(
            "No pair of plates has any strain measured on both plates, so there is no "
            "plate-to-plate correlation to show."
        )

if plate_pair_info:
    titer_label = (
        "titer"
        if dilution_factor_or_concentration == "dilution_factor"
        else f"titer ({concentration_units})"
    )

    # A dropdown is used to select which strains to show rather than binding the
    # selection to the QC-status legend. Binding to the legend does not work here: the
    # legend of a concatenated chart belongs to whichever subplot draws it, so the
    # binding never reaches the other subplots. An input binding does work, because the
    # param is declared on the concatenated chart itself (see `add_params` below) and
    # so is in scope for every subplot. Selecting "all" shows every strain.
    qc_statuses = ["retained in final titers", "dropped by per-serum QC"]
    qc_selection = alt.selection_point(
        fields=["qc_status"],
        bind=alt.binding_select(
            options=[None, *qc_statuses],
            labels=["all", *qc_statuses],
            name="show strains",
        ),
        name="qc_status_selection",
    )

    # mousing over a strain highlights it on every subplot, not just the one being
    # moused over, because this param is also declared on the concatenated chart
    virus_hover = alt.selection_point(
        fields=["virus"], on="mouseover", empty=False, name="virus_hover"
    )

    def scale_for(info):
        """Log scale shared by both axes of a subplot."""
        return alt.Scale(type="log", nice=False, domain=info["domain"])

    def fmt_r(r):
        """Format a Pearson R that is `None` when there were too few strains."""
        return "not computed (< 3 strains)" if r is None else f"{r:.2f}"

    charts = []
    for i, info in enumerate(plate_pair_info):
        base = alt.Chart(info["titers"]).transform_lookup(
            lookup="virus",
            from_=alt.LookupData(
                virus_qc_lookup,
                key="virus",
                fields=["qc_status", "qc_reason"],
            ),
        )
        x_title = f"{titer_label} on {info['plate_x']}"
        y_title = f"{titer_label} on {info['plate_y']}"
        points = (
            base.transform_filter(qc_selection)
            .encode(
                alt.X("titer_x:Q", title=x_title, scale=scale_for(info)),
                alt.Y("titer_y:Q", title=y_title, scale=scale_for(info)),
                alt.Color(
                    "qc_status:N",
                    scale=alt.Scale(
                        domain=qc_statuses,
                        range=["MediumBlue", "OrangeRed"],
                    ),
                    legend=(
                        alt.Legend(
                            title=[
                                "strain QC status",
                                "(use the 'show strains'",
                                "dropdown below the plots",
                                "to show only one of these)",
                            ],
                            titleLimit=400,
                            symbolLimit=0,
                        )
                        if i == 0
                        else None
                    ),
                ),
                alt.Shape(
                    "censored:N",
                    title=["is either median titer", "censored at a bound?"],
                    legend=alt.Legend(titleLimit=400) if i == 0 else None,
                ),
                tooltip=[
                    alt.Tooltip("virus:N", title="strain"),
                    alt.Tooltip(
                        "titer_x:Q",
                        title=f"{titer_label} on {info['plate_x']}",
                        format=".3g",
                    ),
                    alt.Tooltip(
                        "titer_y:Q",
                        title=f"{titer_label} on {info['plate_y']}",
                        format=".3g",
                    ),
                    alt.Tooltip(
                        "n_barcodes_x:Q", title=f"barcodes on {info['plate_x']}"
                    ),
                    alt.Tooltip(
                        "n_barcodes_y:Q", title=f"barcodes on {info['plate_y']}"
                    ),
                    alt.Tooltip("qc_status:N", title="QC status"),
                    alt.Tooltip("qc_reason:N", title="QC detail"),
                ],
                # the moused-over strain gets a larger point with a thick red stroke
                size=alt.condition(virus_hover, alt.value(180), alt.value(60)),
                stroke=alt.condition(virus_hover, alt.value("red"), alt.value("black")),
                strokeWidth=alt.condition(virus_hover, alt.value(3), alt.value(0.5)),
                fillOpacity=alt.condition(virus_hover, alt.value(1), alt.value(0.5)),
            )
            .mark_point(filled=True, strokeOpacity=1)
        )
        # y = x, drawn from the ends of the axes rather than from the titers so that
        # the line spans the whole subplot. This layer repeats the axis titles of the
        # points layer rather than leaving them unset, because the layers share their
        # axes and a title of `None` here blanks the shared title.
        diagonal = (
            alt.Chart(info["diagonal"])
            .encode(
                alt.X("titer:Q", title=x_title, scale=scale_for(info)),
                alt.Y("titer:Q", title=y_title, scale=scale_for(info)),
            )
            .mark_line(color="gray", strokeDash=[4, 4], strokeWidth=1)
        )
        charts.append(
            (diagonal + points).properties(
                width=200,
                height=200,
                title=alt.Title(
                    info["pair"],
                    subtitle=[
                        f"R = {fmt_r(info['r'])} (all {info['n']} strains)",
                        (
                            f"R = {fmt_r(info['r_kept'])} "
                            f"({info['n_kept']} retained strains)"
                        ),
                    ],
                    fontSize=14,
                    subtitleFontSize=10,
                ),
            )
        )
    plate_correlation_chart = (
        alt.concat(*charts, columns=min(4, len(charts)))
        # both params are declared here rather than on the subplots so that they are
        # in scope for all of them
        .add_params(qc_selection, virus_hover).properties(
            title=alt.Title(
                f"plate-to-plate correlation for {group} {serum}",
                subtitle=(
                    "each point is a viral strain; titers are medians over each "
                    "plate's barcodes"
                ),
                fontSize=15,
                anchor="middle",
            )
        )
        # `titleFontWeight` is set because axis titles are bold by default, and here
        # only the overall and subplot titles should be bold. `labelOverlap` drops
        # tick labels rather than letting them overlap, and is set here rather than on
        # each axis so that the chart spec carries it once instead of once per axis.
        .configure_axis(
            grid=False,
            labelOverlap=True,
            titleFontSize=13,
            titleFontWeight="normal",
        )
        # `limit=0` means no limit, so that neither the titles nor the subtitles of the
        # plot or its subplots are truncated. It is set in the chart's config rather
        # than on each title so that the spec carries it once.
        .configure_title(limit=0)
    )
    report.chart(plate_correlation_chart)

report.md("""
    ## Get and plot the neutralization curves for all retained viruses
    First, get the `CurveFits` for just those retained viruses (dropping ones that fail
    QC), and plot:
    """)

fits_qc = neutcurve.CurveFits.combineCurveFits(
    [fits_noqc],
    viruses=[v for v in viruses if v not in viruses_to_drop],
)
assert len(viruses) == len(fits_qc.viruses[serum]) + len(viruses_to_drop)
fig_retained, _ = fits_qc.plotReplicates(
    attempt_shared_legend=False,
    legendfontsize=8,
    titlesize=12,
    ncol=4,
    heightscale=1.2,
    widthscale=1.2,
    subplot_titles="{virus}",
    viruses=[v for v in viruses if v not in viruses_to_drop],
    draw_in_bounds=True,
)
_ = fig_retained.suptitle(
    f"neutralization curves for retained viruses for {group} {serum}",
    y=1,
    fontsize=18,
    fontweight="bold",
)
fig_retained.tight_layout()
report.figure(fig_retained, curve_display_method)
report.md(f"Saving plot of curves to `{curves_pdf}`", log=True)
# `CreationDate` is omitted so that the same curves always give the same PDF; matplotlib
# otherwise stamps the current time into it, as `neutcurve.fig_utils.fig_html` also avoids
fig_retained.savefig(curves_pdf, metadata={"CreationDate": None})
plt.close(fig_retained)

report.md("Save the `CurveFits` to a pickle file:")

with open(output_pickle, "wb") as f_out_pickle:
    pickle.dump(fits_qc, f_out_pickle)

report.md(f"Writing curve fits to {output_pickle}", log=True)

report.md("Write the titers (excluding QC dropped viruses) to a CSV:")

report.md(f"Writing titers to `{titers_csv}`", log=True)
(
    median_titers_noqc.query("virus not in @viruses_to_drop")[
        [
            "group",
            "serum",
            "virus",
            "titer",
            "titer_bound",
            "titer_sem",
            "n_replicates",
            "titer_as",
            "titer_units",
        ]
    ].to_csv(titers_csv, index=False, float_format="%.4g")
)

report.write(snakemake.output.html)
print(f"Wrote the report to {snakemake.output.html}")

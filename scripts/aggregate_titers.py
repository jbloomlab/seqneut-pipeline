"""Aggregate titers across all sera."""

import pickle
import sys

import altair as alt
import neutcurve
import pandas as pd
from seqneut_report import Report

# `noqa: SIM115` as this log file must stay open for the life of the script
sys.stderr = sys.stdout = open(snakemake.log[0], "w")  # noqa: SIM115

_ = alt.data_transformers.disable_max_rows()

input_pickles = snakemake.input.pickles
input_titers = snakemake.input.titers
output_pickles = snakemake.output.pickles
output_titers = snakemake.output.titers
titers_charts = snakemake.output.titers_charts
groups_sera = snakemake.params.groups_sera
groups = snakemake.params.groups
groups_dilution_factor_or_concentration = (
    snakemake.params.groups_dilution_factor_or_concentration
)
groups_concentration_units = snakemake.params.groups_concentration_units
viral_strain_plot_order = snakemake.params.viral_strain_plot_order

report = Report(title="Aggregate titers across all sera")

report.md("""
    ## Write merged titers for each group
    Get the merged titers and merged `CurveFits` object:
    """)

assert len(input_titers) == len(input_pickles) == len(groups_sera)
assert len(groups) == len(output_titers) == len(output_pickles)

titers = pd.concat(
    [pd.read_csv(f_titer) for f_titer in input_titers], ignore_index=True
)
assert len(titers) == len(titers.groupby(["group", "serum", "virus"]))
for group, f_out_titer in zip(groups, output_titers):
    report.md(f"Writing aggregated titers for `{group}` to `{f_out_titer}`")
    print(f"Writing aggregated titers for {group} to {f_out_titer}")
    titers.query("group == @group").to_csv(
        f_out_titer, index=False, float_format="%.4g"
    )

for group, f_out_pickle in zip(groups, output_pickles):
    fits_list = []
    for (g, serum), pickle_f in zip(groups_sera, input_pickles):
        if g == group:
            with open(pickle_f, "rb") as f_in_pickle:
                fits_list.append(pickle.load(f_in_pickle))
    curvefits = neutcurve.CurveFits.combineCurveFits(fits_list)
    report.md(f"Pickling aggregated `CurveFits` for `{group}` to `{f_out_pickle}`")
    print(f"Pickling aggregated CurveFits for {group} to {f_out_pickle}")
    with open(f_out_pickle, "wb") as f_out_pickle_file:
        pickle.dump(curvefits, f_out_pickle_file)

report.md("## Plot all the titers")

if viral_strain_plot_order is None:
    viral_strain_plot_order_used = sorted(curvefits.allviruses)
elif set(curvefits.allviruses) - set(viral_strain_plot_order):
    raise ValueError(
        "not all viruses in `viral_strain_plot_order`:\n"
        + f"{set(curvefits.allviruses) - set(viral_strain_plot_order)=}"
    )
else:
    viral_strain_plot_order_used = viral_strain_plot_order

# one chart per group; each group has a single titer unit and so a single
# correct x-axis label (avoids mixing reciprocal dilutions and concentrations)
for group, titers_chart_html in zip(groups, titers_charts):
    group_titers = titers[titers["group"] == group]
    viruses = [
        v for v in viral_strain_plot_order_used if v in set(group_titers["virus"])
    ]
    group_sera = sorted(group_titers["serum"].unique())

    if groups_dilution_factor_or_concentration[group] == "dilution_factor":
        x_title = "neutralization titer"
    else:
        x_title = f"titer ({groups_concentration_units[group]})"

    serum_selection = alt.selection_point(
        fields=["serum"],
        bind="legend",
        toggle="true",
        empty=False,
    )

    base_chart = (
        alt.Chart(group_titers)
        .add_params(serum_selection)
        .encode(
            alt.X(
                "titer",
                title=x_title,
                scale=alt.Scale(nice=False, padding=5, type="log"),
                axis=alt.Axis(labelOverlap=True),
            ),
        )
        .properties(height=alt.Step(10))
    )

    median_titers_chart = (
        base_chart.encode(
            alt.Y("virus", title=None, sort=viruses, axis=alt.Axis(labelLimit=500)),
            color=alt.value("black"),
            tooltip=["virus", alt.Tooltip("titer", format=".3g")],
        )
        .transform_aggregate(titer="median(titer)", groupby=["virus"])
        .mark_line(point={"filled": True, "size": 50})
        .properties(width=125, title="median titer")
    )

    per_titers_chart = (
        base_chart.encode(
            alt.Y(
                "virus",
                title=None,
                sort=viruses,
                axis=alt.Axis(labels=False, ticks=False),
            ),
            alt.Detail("serum:N"),
            alt.Color(
                "serum:N",
                scale=alt.Scale(range=["black", "black"]),
                legend=alt.Legend(
                    orient="right",
                    columns=1 + int(1.6 * len(group_sera)) // max(len(viruses), 1),
                    symbolLimit=0,
                    symbolFillColor="black",
                    title="serum (click to select)",
                ),
                sort=group_sera,
            ),
            opacity=alt.condition(serum_selection, alt.value(1), alt.value(0.15)),
            strokeWidth=alt.condition(serum_selection, alt.value(3), alt.value(1)),
            tooltip=[
                (alt.Tooltip(c, format=".3g") if group_titers[c].dtype == float else c)
                for c in group_titers.columns
            ],
        )
        .mark_line(point={"filled": True, "size": 20})
        .properties(width=225, title="per-serum titers")
    )

    titers_chart = (
        alt.hconcat(median_titers_chart, per_titers_chart, spacing=5)
        .configure_axis(grid=False)
        .configure_legend(padding=10, labelOffset=2, columnPadding=8, labelLimit=400)
        .properties(
            title=alt.TitleParams(
                f"Interactive chart of neutralization titers for {group}",
                subtitle=[
                    "Median titers are at left; per-serum titers are at right.",
                    "Mouseover points for details.",
                    (
                        "Click points on serum legend at right to select sera "
                        "to bold on per-serum chart."
                    ),
                ],
                fontSize=15,
                dx=100,
                dy=-5,
            ),
            autosize=alt.AutoSizeParams(resize=True),
        )
    )

    report.md(f"Saving chart for `{group}` to `{titers_chart_html}`")
    print(f"Saving chart for {group} to {titers_chart_html}")
    titers_chart.save(titers_chart_html)
    report.chart(titers_chart)

report.write(snakemake.output.html)
print(f"Wrote the report to {snakemake.output.html}")

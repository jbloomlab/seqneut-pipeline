"""Snakemake file for the ``seqneut-pipeline`.

Designed to be included in another ``Snakefile`` that specifies the config.

"""

import pandas as pd
import re

snakemake.utils.min_version("9.0")


include: "funcs.smk"  # include functions


# --- Process configuration ------------------------------------------------------------

pipeline_subdir = config["seqneut-pipeline"]

viral_libraries = {
    lib: pd.read_csv(f) for (lib, f) in config["viral_libraries"].items()
}

if ("viral_strain_plot_order" not in config) or (
    config["viral_strain_plot_order"] is None
):
    viral_strain_plot_order = None
else:
    viral_strain_plot_order = pd.read_csv(config["viral_strain_plot_order"])[
        "strain"
    ].tolist()
    assert len(viral_strain_plot_order) == len(set(viral_strain_plot_order))

neut_standard_sets = {
    s: pd.read_csv(f) for (s, f) in config["neut_standard_sets"].items()
}

stringify_plate_dates(config)  # so `snakemake` can JSON serialize the config

# whether to sum the counts of all barcodes for a strain before any curve fitting; read
# before the plates are processed, as their validation depends on it
collapse_strain_barcodes = get_collapse_strain_barcodes(config)

for lib, lib_df in viral_libraries.items():
    check_no_whitespace(lib_df["barcode"], f"barcodes in viral library {lib!r}")
    if collapse_strain_barcodes:
        # collapsing names a strain's barcode for the strain, so the strain name is
        # then subject to the same constraint as a barcode
        check_no_whitespace(lib_df["strain"], f"strains in viral library {lib!r}")
for neut_standard_set, neut_standard_df in neut_standard_sets.items():
    check_no_whitespace(
        neut_standard_df["barcode"],
        f"barcodes in neut standard set {neut_standard_set!r}",
    )

# `plates` may be absent, null, or empty for a project that only has
# `miscellaneous_plates`, and so just counts barcodes rather than fitting any curves
plates = {
    str(plate): process_plate(str(plate), plate_params)
    for (plate, plate_params) in (config.get("plates") or {}).items()
}
check_no_whitespace(plates, "plate names")

groups = sorted(set(plate_params["group"] for plate_params in plates.values()))
groups_cannot_contain = ["|", "_"]  # wildcard problems if group contains these
if any(s in group for s in groups_cannot_contain for group in groups):
    raise ValueError(f"found {groups_cannot_contain=} character in {groups=}")

# all plates in a group must share 'dilution_factor_or_concentration' and
# 'concentration_units'; build per-group maps of these values (fail fast otherwise)
group_dilution_factor_or_concentration = {}
group_concentration_units = {}
for plate_d in plates.values():
    plate_group = plate_d["group"]
    plate_dfc = plate_d["dilution_factor_or_concentration"]
    plate_units = plate_d["concentration_units"]
    if plate_group not in group_dilution_factor_or_concentration:
        group_dilution_factor_or_concentration[plate_group] = plate_dfc
        group_concentration_units[plate_group] = plate_units
    elif (group_dilution_factor_or_concentration[plate_group] != plate_dfc) or (
        group_concentration_units[plate_group] != plate_units
    ):
        raise ValueError(
            f"plates in group {plate_group!r} do not all share the same "
            "'dilution_factor_or_concentration' and 'concentration_units'"
        )


# the plates each plate's titers can be correlated with, and the subset of plates that
# have at least one such comparator and so get a plate-to-plate correlation report
# (keyed to their group). Both are determined by the configuration rather than by the
# QC, so the set of correlation reports is known without the `groups_sera_by_plate`
# checkpoint even though the inputs to each of them are not.
plate_comparators = get_plate_comparators(plates)
corr_plates = {
    plate: plates[plate]["group"]
    for (plate, comparators) in plate_comparators.items()
    if comparators
}


wildcard_constraints:
    group="|".join(groups),


# like `plates`, may be absent or null, both of which mean no sera override the defaults
sera_override_defaults = config.get("sera_override_defaults") or {}

if not set(sera_override_defaults).issubset(groups):
    raise ValueError(f"{sera_override_defaults=} keyed by invalid groups")


if plates == {}:
    samples = {}
else:
    samples = pd.concat(
        [plate_d["samples"] for plate_d in plates.values()],
        ignore_index=True,
    )
    assert samples["sample"].nunique() == len(samples)
    samples = samples.set_index("sample").to_dict(orient="index")


miscellaneous_plates = process_miscellaneous_plates(
    config.get("miscellaneous_plates") or {}
)

# define `add_htmls_to_docs` if not already defined.
try:
    add_htmls_to_docs
except NameError:  # if not defined
    add_htmls_to_docs = {}


# --- Snakemake rules -------------------------------------------------------------------

if plates:

    rule count_barcodes:
        """Count barcodes for a sample."""
        input:
            fastq=lambda wc: samples[wc.sample]["fastq"],
        output:
            counts="results/barcode_counts/{sample}.csv",
            invalid="results/barcode_invalid/{sample}.csv",
            fates="results/barcode_fates/{sample}.csv",
        log:
            "results/logs/count_barcodes_{sample}.txt",
        conda:
            "environment.yml"
        params:
            viral_barcodes=lambda wc: sorted(
                viral_libraries[plates[samples[wc.sample]["plate"]]["viral_library"]][
                    "barcode"
                ]
            ),
            neut_standard_barcodes=lambda wc: sorted(
                neut_standard_sets[
                    plates[samples[wc.sample]["plate"]]["neut_standard_set"]
                ]["barcode"]
            ),
            illumina_barcode_parser_params=lambda wc: plates[
                samples[wc.sample]["plate"]
            ]["illumina_barcode_parser_params"],
        script:
            "scripts/count_barcodes.py"

    rule process_plate:
        """Process a plate to QC and convert counts to fraction infectivity."""
        input:
            report_module=os.path.join(pipeline_subdir, "scripts/seqneut_report.py"),
            count_csvs=lambda wc: expand(
                rules.count_barcodes.output.counts,
                sample=plates[wc.plate]["samples"]["sample"],
            ),
            fate_csvs=lambda wc: expand(
                rules.count_barcodes.output.fates,
                sample=plates[wc.plate]["samples"]["sample"],
            ),
        output:
            html="results/plates/{plate}/process_{plate}.html",
            qc_drops="results/plates/{plate}/qc_drops.yml",
            frac_infectivity_csv="results/plates/{plate}/frac_infectivity.csv",
            fits_csv="results/plates/{plate}/curvefits.csv",
            fits_pickle="results/plates/{plate}/curvefits.pickle",
        log:
            "results/logs/process_{plate}.txt",
        conda:
            "environment.yml"
        params:
            # pass DataFrames/Series as dict/list for snakemake params rerun triggers
            viral_barcodes=lambda wc: (
                viral_libraries[plates[wc.plate]["viral_library"]][
                    ["barcode", "strain"]
                ]
                .sort_values("barcode")
                .reset_index(drop=True)
                .to_dict()
            ),
            neut_standard_barcodes=lambda wc: (
                neut_standard_sets[plates[wc.plate]["neut_standard_set"]][["barcode"]]
                .sort_values("barcode")
                .reset_index(drop=True)
                .to_dict()
            ),
            samples=lambda wc: plates[wc.plate]["samples"]["sample"].tolist(),
            plate_params=lambda wc: {
                param: (val if param != "samples" else val.to_dict())
                for (param, val) in plates[wc.plate].items()
            },
            curve_display_method=config["curve_display_method"],
            collapse_strain_barcodes=collapse_strain_barcodes,
        script:
            "scripts/process_plate.py"

    checkpoint groups_sera_by_plate:
        """Get list of all groups/sera and plates they are on."""
        input:
            csvs=expand(rules.process_plate.output.fits_csv, plate=plates),
        output:
            csv="results/sera/groups_sera_by_plate.csv",
        log:
            "results/logs/groups_sera_by_plate.txt",
        conda:
            "environment.yml"
        params:
            plates=list(plates),
        script:
            "scripts/groups_sera_by_plate.py"

    rule group_serum_titers:
        """Aggregate and analyze titers for a group / serum."""
        input:
            report_module=os.path.join(pipeline_subdir, "scripts/seqneut_report.py"),
            funcs_module=os.path.join(pipeline_subdir, "scripts/seqneut_funcs.py"),
            pickles=lambda wc: [
                rules.process_plate.output.fits_pickle.format(plate=plate)
                for plate in groups_sera_plates()[(wc.group, wc.serum)]
            ],
            fits_csvs=lambda wc: [
                rules.process_plate.output.fits_csv.format(plate=plate)
                for plate in groups_sera_plates()[(wc.group, wc.serum)]
            ],
        output:
            html="results/sera/{group}_{serum}/{group}_{serum}_titers.html",
            per_rep_titers="results/sera/{group}_{serum}/titers_per_replicate.csv",
            titers="results/sera/{group}_{serum}/titers.csv",
            curves_pdf="results/sera/{group}_{serum}/curves.pdf",
            pickle="results/sera/{group}_{serum}/curvefits.pickle",
            qc_drops="results/sera/{group}_{serum}/qc_drops.yml",
        log:
            "results/logs/group_serum_titers_{group}_{serum}.txt",
        conda:
            "environment.yml"
        params:
            viral_strain_plot_order=viral_strain_plot_order,
            curve_display_method=config["curve_display_method"],
            dilution_factor_or_concentration=lambda wc: (
                group_dilution_factor_or_concentration[wc.group]
            ),
            concentration_units=lambda wc: group_concentration_units[wc.group],
            serum_titer_as=lambda wc: (
                sera_override_defaults[wc.group][wc.serum]["titer_as"]
                if (
                    (wc.group in sera_override_defaults)
                    and (wc.serum in sera_override_defaults[wc.group])
                    and ("titer_as" in sera_override_defaults[wc.group][wc.serum])
                )
                else config["default_serum_titer_as"]
            ),
            qc_thresholds=lambda wc: (
                sera_override_defaults[wc.group][wc.serum]["qc_thresholds"]
                if (
                    (wc.group in sera_override_defaults)
                    and (wc.serum in sera_override_defaults[wc.group])
                    and ("qc_thresholds" in sera_override_defaults[wc.group][wc.serum])
                )
                else config["default_serum_qc_thresholds"]
            ),
        script:
            "scripts/group_serum_titers.py"

    rule plate_to_plate_corr:
        """Correlate the titers on a plate with those on the other plates of its group.

        Only created for the plates in `corr_plates`, which are those sharing a serum
        with another plate of their group. Each pair of plates is therefore compared
        twice, once in each plate's report, which keeps a report from growing with
        the square of the number of plates in the group.

        The script handles the case where the QC has dropped every titer that this
        plate shared with its comparators, by reporting that there is nothing to
        correlate.

        """
        input:
            report_module=os.path.join(pipeline_subdir, "scripts/seqneut_report.py"),
            funcs_module=os.path.join(pipeline_subdir, "scripts/seqneut_funcs.py"),
            # The sera measured on this plate. Each of these files holds that serum's
            # titers from every plate it was measured on, which is what the correlations
            # need. All of the plate's sera are read rather than just those also measured
            # elsewhere, so that the report can list the ones that are not compared.
            per_rep_titers=lambda wc: [
                rules.group_serum_titers.output.per_rep_titers.format(
                    group=wc.group, serum=serum
                )
                for ((group, serum), serum_plates) in groups_sera_plates().items()
                if (group == wc.group) and (wc.plate in serum_plates)
            ],
        output:
            html="results/plate_to_plate_corrs/plate_to_plate_corr_{group}_{plate}.html",
            corrs_csv="results/plate_to_plate_corrs/plate_to_plate_corr_{group}_{plate}.csv",
        log:
            "results/logs/plate_to_plate_corr_{group}_{plate}.txt",
        wildcard_constraints:
            # a group cannot contain '_' but a plate can, so without this the plate part
            # of the path would be split at its first '_'
            plate="|".join(re.escape(plate) for plate in plates),
        conda:
            "environment.yml"
        params:
            plate_date=lambda wc: plates[wc.plate]["date"],
            # keyed by plate in config order, which orders the panels of the plots
            comparator_dates=lambda wc: {
                other: plates[other]["date"] for other in plate_comparators[wc.plate]
            },
            dilution_factor_or_concentration=lambda wc: (
                group_dilution_factor_or_concentration[wc.group]
            ),
            concentration_units=lambda wc: group_concentration_units[wc.group],
        script:
            "scripts/plate_to_plate_corr.py"

    rule aggregate_titers:
        """Aggregate all serum titers."""
        input:
            report_module=os.path.join(pipeline_subdir, "scripts/seqneut_report.py"),
            pickles=lambda wc: [
                rules.group_serum_titers.output.pickle.format(group=group, serum=serum)
                for (group, serum) in groups_sera_plates()
            ],
            titers=lambda wc: [
                rules.group_serum_titers.output.titers.format(group=group, serum=serum)
                for (group, serum) in groups_sera_plates()
            ],
        output:
            html="results/aggregated_titers/aggregate_titers.html",
            pickles=[
                f"results/aggregated_titers/curvefits_{group}.pickle"
                for group in groups
            ],
            titers=[f"results/aggregated_titers/titers_{group}.csv" for group in groups],
            titers_charts=[
                f"results/aggregated_titers/titers_{group}.html" for group in groups
            ],
        log:
            "results/logs/aggregate_titers.txt",
        conda:
            "environment.yml"
        params:
            viral_strain_plot_order=viral_strain_plot_order,
            groups_sera=lambda wc: list(groups_sera_plates()),
            groups=groups,
            groups_dilution_factor_or_concentration=group_dilution_factor_or_concentration,
            groups_concentration_units=group_concentration_units,
        script:
            "scripts/aggregate_titers.py"

    rule aggregate_qc_drops:
        """Aggregate all QC drops."""
        input:
            # `snakemake` re-runs a `script:` rule when the script itself changes, but
            # does not follow its imports, so the modules it imports are inputs
            report_module=os.path.join(pipeline_subdir, "scripts/seqneut_report.py"),
            plate_qc_drops=expand(rules.process_plate.output.qc_drops, plate=plates),
            groups_sera_qc_drops=lambda wc: [
                rules.group_serum_titers.output.qc_drops.format(
                    group=group, serum=serum
                )
                for (group, serum) in groups_sera_plates()
            ],
        output:
            html="results/qc_drops/aggregate_qc_drops.html",
            plate_qc_drops="results/qc_drops/plate_qc_drops.yml",
            barcode_qc_drops="results/qc_drops/barcode_qc_drops.yml",
            groups_sera_qc_drops="results/qc_drops/groups_sera_qc_drops.yml",
        log:
            "results/logs/aggregate_qc_drops.txt",
        conda:
            "environment.yml"
        params:
            plates=list(plates),
            groups_sera=lambda wc: list(groups_sera_plates()),
        script:
            "scripts/aggregate_qc_drops.py"


rule miscellaneous_plate_count_barcodes:
    """Count barcodes for a well in a miscellaneous plate."""
    input:
        fastq=lambda wc: miscellaneous_plates[wc.misc_plate]["wells"][wc.well],
    output:
        counts="results/miscellaneous_plates/{misc_plate}/{well}_counts.csv",
        invalid="results/miscellaneous_plates/{misc_plate}/{well}_invalid.csv",
        fates="results/miscellaneous_plates/{misc_plate}/{well}_fates.csv",
    log:
        "results/logs/miscellaneous_plate_count_barcodes_{misc_plate}_{well}.txt",
    conda:
        "environment.yml"
    params:
        viral_barcodes=lambda wc: sorted(
            viral_libraries[miscellaneous_plates[wc.misc_plate]["viral_library"]][
                "barcode"
            ]
        ),
        neut_standard_barcodes=lambda wc: sorted(
            neut_standard_sets[
                miscellaneous_plates[wc.misc_plate]["neut_standard_set"]
            ]["barcode"]
        ),
        illumina_barcode_parser_params=lambda wc: miscellaneous_plates[wc.misc_plate][
            "illumina_barcode_parser_params"
        ],
    script:
        "scripts/count_barcodes.py"


rule build_docs:
    """Build the HTML documentation.

    Defined outside the `if plates` block above so that a project with only
    `miscellaneous_plates` still gets docs of the HTMLs it adds with
    `add_htmls_to_docs`. Everything derived from the plates is then empty, and
    `scripts/build_docs.py` omits the sections that describe them.

    """
    input:
        lambda wc: [
            f
            for d in add_htmls_to_docs.values()
            for v in d.values()
            for f in (v.values() if isinstance(v, dict) else [v])
        ],
        titers_charts=(rules.aggregate_titers.output.titers_charts if plates else []),
        serum_titers_htmls=lambda wc: (
            [
                f"results/sera/{group}_{serum}/{group}_{serum}_titers.html"
                for (group, serum) in groups_sera_plates()
            ]
            if plates
            else []
        ),
        process_plates_htmls=expand(
            "results/plates/{plate}/process_{plate}.html",
            plate=plates,
        ),
        plate_corr_htmls=[
            rules.plate_to_plate_corr.output.html.format(group=group, plate=plate)
            for (plate, group) in corr_plates.items()
        ],
        qc_drops_html=(rules.aggregate_qc_drops.output.html if plates else []),
    output:
        docs=directory("results/docs"),
    log:
        "results/logs/build_docs.txt",
    conda:
        "environment.yml"
    params:
        description=config["description"],
        groups_sera=lambda wc: (list(groups_sera_plates()) if plates else []),
        plates={plate: plates[plate]["group"] for plate in plates},
        # keyed by plate in the same order as `plate_corr_htmls`
        corr_plates=corr_plates,
        add_htmls_to_docs=lambda wc: {
            key: {
                key2: (
                    {k3: str(v3) for k3, v3 in val2.items()}
                    if isinstance(val2, dict)
                    else str(val2)
                )
                for (key2, val2) in val.items()
            }
            for (key, val) in add_htmls_to_docs.items()
        },
    script:
        "scripts/build_docs.py"


seqneut_pipeline_outputs = [
    *(
        [
            rules.aggregate_titers.output.titers,
            rules.aggregate_titers.output.pickles,
            rules.aggregate_qc_drops.output.plate_qc_drops,
            rules.aggregate_qc_drops.output.barcode_qc_drops,
            rules.aggregate_qc_drops.output.groups_sera_qc_drops,
            *[
                rules.plate_to_plate_corr.output.corrs_csv.format(
                    group=group, plate=plate
                )
                for (plate, group) in corr_plates.items()
            ],
        ]
        if plates
        else []
    ),
    *[
        f"results/miscellaneous_plates/{plate}/{well}_{suffix}"
        for plate in miscellaneous_plates
        for well in miscellaneous_plates[plate]["wells"]
        for suffix in ["counts.csv", "fates.csv"]
    ],
    # the report HTMLs are not listed here, as they are inputs to `build_docs`
    rules.build_docs.output.docs,
]

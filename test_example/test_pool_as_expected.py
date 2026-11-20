"""Test the pool analyses measure what was planted in their synthetic plates.

The plates are written by `fastqs/pool/make_pool_fastqs.py`. Every read of them was built to
parse, so each barcode's share of the pool, the dilution the composition is read from, and
the linear range found must match the plan exactly, to the precision the CSVs are written at.

"""

import pandas as pd

# the analyses of `config.yml`, each keyed by the plate whose plan it should reproduce
ANALYSES = {
    "pool_whole_library": "pool_plate_1",
    "pool_subgroup_by_stock": "pool_plate_2",
    "pool_no_linear_range": "pool_plate_3",
}

# the shares are written to four significant digits, so they agree to about that
SHARE_TOLERANCE = 1e-3

expected_composition = pd.read_csv("expected_pool_composition.csv")
expected_summary = pd.read_csv("expected_pool_summary.csv")

for analysis, plate in ANALYSES.items():
    actual = pd.read_csv(f"results/pool_composition/{analysis}_barcode_composition.csv")
    planted = expected_composition[expected_composition["pool"] == plate]
    assert len(planted), f"nothing planted for {plate}"

    # every member barcode is reported, including any that got no reads
    assert set(actual["barcode"]) == set(planted["barcode"]), (
        f"{analysis} reports barcodes {sorted(set(actual['barcode']) ^ set(planted['barcode']))} "
        "that were not planted, or the other way round"
    )

    merged = actual.merge(planted, on="barcode", suffixes=("_found", "_planted"))
    for column in ["mean_fraction", "stock_measured_share"]:
        if f"{column}_planted" not in merged.columns:
            continue
        diff = (merged[f"{column}_found"] - merged[f"{column}_planted"]).abs()
        assert (
            diff < SHARE_TOLERANCE
        ).all(), f"{analysis} {column} differs from what was planted by up to {diff.max():.3g}:\n" + merged.loc[
            diff >= SHARE_TOLERANCE,
            ["barcode", f"{column}_found", f"{column}_planted"],
        ].to_string(
            index=False
        )

    # the divided share is reported only where the pool was made with equal volumes
    planted_divided = "stock_divided_share" in planted.columns and (
        planted["stock_divided_share"].notna().any()
    )
    assert planted_divided == ("stock_divided_share" in actual.columns), (
        f"{analysis} reports a divided per-stock share where none was planted, or the "
        "other way round"
    )

    summary = pd.read_csv(f"results/pool_composition/{analysis}_pool_summary.csv")
    assert len(summary) == 1, f"{analysis} summary is not one row"
    summary = summary.iloc[0]
    planted_summary = expected_summary[expected_summary["pool"] == plate].iloc[0]
    for column in [
        "analysis_dilution",
        "linear_range_dilutions",
        "linear_range_strict",
        "analysis_wells_given",
        "n_strains",
    ]:
        found, want = summary[column], planted_summary[column]
        if pd.isna(want):
            assert pd.isna(
                found
            ), f"{analysis} {column} is {found!r}, but nothing should be reported"
        else:
            assert str(found) == str(
                want
            ), f"{analysis} {column} is {found!r} but {want!r} was planted"
    # a pool with no linear range has no fit to take a titer from
    if pd.isna(planted_summary["linear_range_dilutions"]):
        assert pd.isna(summary["titer_vs_neut_standard"]), (
            f"{analysis} reports a titer of {summary['titer_vs_neut_standard']} although "
            "no run of dilutions qualified"
        )
    else:
        assert summary["titer_vs_neut_standard"] > 0, (
            f"{analysis} reports no titer although the run "
            f"{planted_summary['linear_range_dilutions']} qualified"
        )

print(f"pool analyses match what was planted in all {len(ANALYSES)} of them")

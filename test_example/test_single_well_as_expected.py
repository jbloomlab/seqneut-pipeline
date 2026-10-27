"""Test the single-well analyses find what was planted in their synthetic plates.

The plates are written by `fastqs/single_well/make_single_well_fastqs.py`. Every read of
them was built to parse, so the categories of each well, the flags they raise, and where
its contaminant came from must match the plan exactly. The titer is compared to the
precision the composition CSV is written at.

"""

import glob

import numpy
import pandas as pd

EXACT_COLUMNS = [
    "n_expected_barcodes",
    "expected_reads",
    "neut_standard_reads",
    "other_library_reads",
    "unmatched_near_expected_reads",
    "unmatched_far_expected_reads",
    "discarded_reads",
    "flagged_other_library",
    "flagged_unmatched",
    "too_shallow",
    "titer_imprecise",
    "top_other_unit_origin",
    # absent from the CSV, and so null, for an analysis with no final library
    "n_expected_barcodes_final",
    "final_library",
]

expected = pd.read_csv("expected_single_well_composition.csv")
actual = pd.concat(
    [
        pd.read_csv(f)
        for f in glob.glob("results/single_well_infections/*_well_composition.csv")
    ],
    ignore_index=True,
)

merged = expected.merge(
    actual,
    on=["plate", "well"],
    how="outer",
    suffixes=("_expected", "_actual"),
    validate="one_to_one",
    indicator=True,
)
missing = merged[merged["_merge"] != "both"]
assert not len(missing), f"wells in only one of the two:\n{missing[['plate', 'well']]}"

for column in EXACT_COLUMNS:
    e, a = merged[f"{column}_expected"], merged[f"{column}_actual"]
    differ = ~((e == a) | (e.isna() & a.isna()))
    assert not differ.any(), (
        f"{column} differs:\n"
        + merged.loc[
            differ, ["plate", "well", f"{column}_expected", f"{column}_actual"]
        ].to_string()
    )

# the expected barcodes given no reads, which the analysis must list and warn about
barcode_reads = pd.concat(
    [
        pd.read_csv(f)
        for f in glob.glob(
            "results/single_well_infections/*_expected_barcode_reads.csv"
        )
    ],
    ignore_index=True,
)
absent_found = set(
    map(tuple, barcode_reads.query("reads == 0")[["plate", "well", "barcode"]].values)
)
absent_planted = {
    (row.plate, row.well, barcode)
    for row in expected.dropna(subset=["absent_barcodes"]).itertuples()
    for barcode in row.absent_barcodes.split(";")
}
assert absent_found == absent_planted, f"{absent_found=} but {absent_planted=}"

# The two wells of one strain on the plate recording stocks hold different tubes of it, so
# the charts must tell them apart by tube.
with open(
    "results/single_well_infections/single_well_stocks_analyze_single_well_infections.html"
) as f:
    stocks_html = f.read()
for label in ["[bc 1]", "[bc 2]"]:
    assert label in stocks_html, f"{label} is not a label of the stocks report"

numpy.testing.assert_allclose(
    merged["titer_actual"], merged["titer_expected"], rtol=1e-3, equal_nan=True
)

print(f"All {len(merged)} wells of the single-well plates are as expected.")

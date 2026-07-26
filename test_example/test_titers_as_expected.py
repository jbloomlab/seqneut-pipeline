"""Test script to test titers as expected.

The expected titers are read from the CSV specified as the first command-line argument,
which defaults to ``expected_titers_for_test.csv``. Pass
``expected_titers_for_test_collapse_strain_barcodes.csv`` for a run that collapses the
barcodes for each strain.

"""

import glob
import sys

import numpy
import pandas as pd

expected_titers_csv = (
    sys.argv[1] if len(sys.argv) > 1 else "expected_titers_for_test.csv"
)
print(f"Reading expected titers from {expected_titers_csv}")

expected_titers = (
    pd.read_csv(expected_titers_csv)
    .assign(log10_titer=lambda x: numpy.log10(x["titer"]))[
        ["group", "serum", "virus", "log10_titer", "n_replicates", "titer_bound"]
    ]
    .sort_values(["group", "serum", "virus"])
    .reset_index(drop=True)
)

actual_titers = (
    pd.concat(
        [pd.read_csv(f) for f in glob.glob("results/aggregated_titers/titers_*.csv")]
    )
    .assign(log10_titer=lambda x: numpy.log10(x["titer"]))[
        ["group", "serum", "virus", "log10_titer", "n_replicates", "titer_bound"]
    ]
    .sort_values(["group", "serum", "virus"])
    .reset_index(drop=True)
)

rtol = 0.02
atol = 0.1

print(f"Comparing {expected_titers} and {actual_titers} with {atol=} and {rtol=}")

pd.testing.assert_frame_equal(
    expected_titers,
    actual_titers,
    check_exact=False,
    atol=atol,
    rtol=rtol,
)

print("Titers are sufficiently similar.")

# The `abconc` group holds the same data as the `pilot` group but titrated as a
# `concentration` set to 1 / dilution_factor rather than as a `dilution_factor`.
# The fits are therefore identical, and the concentration-mode titer (the IC50
# reported directly, with no reciprocal) must equal the reciprocal of the
# corresponding dilution-mode `pilot` titer. This is checked at the per-replicate
# (per-barcode) level, where the relationship is exact -- it does not hold for
# medians because the median of reciprocals is not the reciprocal of the median.


def _read_per_rep(group):
    dfs = [
        pd.read_csv(f)
        for f in glob.glob(f"results/sera/{group}_*/titers_per_replicate.csv")
    ]
    assert dfs, f"no per-replicate titers found for group {group!r}"
    return pd.concat(dfs, ignore_index=True).assign(
        barcode=lambda x: x["replicate"].str.rsplit("-", n=1).str[-1]
    )


pilot_rep = _read_per_rep("pilot")
abconc_rep = _read_per_rep("abconc")

# units are reported per the mode
assert (pilot_rep["titer_units"] == "reciprocal_dilution").all()
assert (abconc_rep["titer_units"] == "ug/ml").all()

reciprocal_check = pilot_rep.merge(
    abconc_rep,
    on=["serum", "virus", "barcode"],
    suffixes=("_pilot", "_abconc"),
    validate="one_to_one",
)
assert len(reciprocal_check) == len(pilot_rep) == len(abconc_rep), (
    f"barcode mismatch: {len(pilot_rep)=}, {len(abconc_rep)=}, "
    f"{len(reciprocal_check)=}"
)
numpy.testing.assert_allclose(
    reciprocal_check["titer_pilot"] * reciprocal_check["titer_abconc"],
    1.0,
    rtol=0.02,
)

print("Concentration-mode titers are reciprocals of dilution-mode titers.")

"""Test script to test titers as expected."""

import glob

import numpy

import pandas as pd


expected_titers = (
    pd.read_csv("expected_titers_for_test.csv")
    .assign(log10_titer=lambda x: numpy.log10(x["titer"]))[
        ["group", "serum", "virus", "log10_titer", "n_replicates", "titer_bound"]
    ]
    .sort_values(["serum", "virus"])
    .reset_index(drop=True)
)

actual_titers = (
    pd.concat(
        [pd.read_csv(f) for f in glob.glob("results/aggregated_titers/titers_*.csv")]
    )
    .assign(log10_titer=lambda x: numpy.log10(x["titer"]))[
        ["group", "serum", "virus", "log10_titer", "n_replicates", "titer_bound"]
    ]
    .sort_values(["serum", "virus"])
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

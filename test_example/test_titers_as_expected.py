"""Test script to test titers as expected."""

import numpy

import pandas as pd


expected_titers = pd.read_csv("expected_titers_for_test.csv").assign(
    log10_titer=lambda x: numpy.log10(x["titer"])
)[["serum", "virus", "log10_titer", "n_replicates", "titer_bound"]]

actual_titers = pd.read_csv("results/aggregated_titers/titers.csv").assign(
    log10_titer=lambda x: numpy.log10(x["titer"])
)[["serum", "virus", "log10_titer", "n_replicates", "titer_bound"]]

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

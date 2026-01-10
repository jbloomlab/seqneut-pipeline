"""Count barcodes from FASTQ barcodes."""

import sys

import dms_variants.illuminabarcodeparser

import pandas as pd


def hamming_distance(s1, s2):
    """Compute Hamming distance between two equal-length strings."""
    return sum(c1 != c2 for c1, c2 in zip(s1, s2, strict=True))


def find_closest_valid_barcode(invalid_bc, valid_barcodes_df):
    """Find closest valid barcode by Hamming distance with tie-breaking."""
    closest = (
        valid_barcodes_df.assign(
            hamming_distance=lambda df: df["barcode"].apply(
                lambda vb: hamming_distance(invalid_bc, vb)
            )
        )
        .sort_values(
            ["hamming_distance", "count", "barcode"], ascending=[True, False, True]
        )
        .iloc[0]
    )
    return closest["barcode"], closest["hamming_distance"], closest["count"]


sys.stderr = sys.stdout = log = open(snakemake.log[0], "w")

viral_barcodes = snakemake.params.viral_barcodes
neut_standard_barcodes = snakemake.params.neut_standard_barcodes

valid_barcodes = list({bc.upper() for bc in viral_barcodes + neut_standard_barcodes})
if len(valid_barcodes) != (len(viral_barcodes) + len(neut_standard_barcodes)):
    raise ValueError("barcodes not all unique")
if not valid_barcodes:
    raise ValueError("no barcodes specified")
bclen = len(valid_barcodes[0])
if not all(bclen == len(bc) for bc in valid_barcodes):
    raise ValueError("barcodes not all the same length")

parser = dms_variants.illuminabarcodeparser.IlluminaBarcodeParser(
    bclen=bclen,
    **snakemake.params.illumina_barcode_parser_params,
)

counts, fates = parser.parse(snakemake.input.fastq, outer_flank_fates=True)

total_counts = counts["count"].sum()

# create valid barcodes DataFrame with counts (needed for closest barcode lookup)
valid_barcodes_df = (
    pd.concat([counts, pd.DataFrame({"barcode": valid_barcodes, "count": 0})])
    .query("barcode in @valid_barcodes")
    .groupby("barcode", as_index=False)
    .aggregate({"count": "sum"})
)

# write valid barcode counts including a 0 count for missing ones
(
    valid_barcodes_df.sort_values(["count", "barcode"], ascending=[False, True])
    .assign(fraction_all_valid_and_invalid_counts=lambda x: x["count"] / total_counts)
    .to_csv(snakemake.output.counts, index=False, float_format="%.3g")
)

# process invalid barcodes and find closest valid barcode for each
invalid_barcodes_df = (
    counts.query("barcode not in @valid_barcodes")
    .assign(fraction_all_valid_and_invalid_counts=lambda x: x["count"] / total_counts)
    .copy()
)

# compute closest valid barcodes
closest_valid_barcode_cols = [
    "closest_valid_barcode",
    "closest_valid_bacode_hamming_distance",
    "closest_valid_barcode_count",
]
if len(invalid_barcodes_df) > 0:
    closest_info = invalid_barcodes_df["barcode"].apply(
        lambda bc: find_closest_valid_barcode(bc, valid_barcodes_df)
    )
    invalid_barcodes_df[closest_valid_barcode_cols] = pd.DataFrame(
        closest_info.tolist(), index=invalid_barcodes_df.index
    )
else:
    for col in closest_valid_barcode_cols:
        invalid_barcodes_df[col] = None

# write invalid barcodes with closest valid barcode information
(
    invalid_barcodes_df.sort_values(
        ["count", "barcode"], ascending=[False, True]
    ).to_csv(snakemake.output.invalid, index=False, float_format="%.3g")
)

# write fates, note we have to add valid and invalid fates
(
    pd.concat(
        [
            fates.query("fate not in ['valid barcode', 'invalid barcode']"),
            pd.DataFrame(
                {
                    "fate": ["valid barcode", "invalid barcode"],
                    "count": [
                        counts.query("barcode in @valid_barcodes")["count"].sum(),
                        counts.query("barcode not in @valid_barcodes")["count"].sum(),
                    ],
                }
            ),
        ]
    )
    .sort_values(["count", "fate"], ascending=False)
    .to_csv(snakemake.output.fates, index=False)
)

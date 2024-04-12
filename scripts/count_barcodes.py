"""Count barcodes from FASTQ barcodes."""

import sys

import dms_variants.illuminabarcodeparser

import pandas as pd


sys.stderr = sys.stdout = log = open(snakemake.log[0], "w")

viral_barcodes = pd.read_csv(snakemake.input.viral_library)["barcode"].tolist()
neut_standard_barcodes = pd.read_csv(snakemake.input.neut_standard_set)[
    "barcode"
].tolist()

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

# write valide barcode counts including a 0 count for missing ones
(
    pd.concat([counts, pd.DataFrame({"barcode": valid_barcodes, "count": 0})])
    .query("barcode in @valid_barcodes")
    .groupby("barcode", as_index=False)
    .aggregate({"count": "sum"})
    .sort_values(["count", "barcode"], ascending=[False, True])
    .to_csv(snakemake.output.counts, index=False)
)

# write invalid barcodes
(
    counts.query("barcode not in @valid_barcodes")
    .sort_values(["count", "barcode"], ascending=[False, True])
    .to_csv(snakemake.output.invalid, index=False)
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

"""Count variants from Illumina barcodes."""

import os
import sys
import pandas as pd
import dms_variants.illuminabarcodeparser

# Write the log file to the specified location
sys.stderr = sys.stdout = log = open(snakemake.log[0], "w")

## ====== Inputs ====== ##
fastq = snakemake.input.fastq
variants = snakemake.input.variants
standards = snakemake.input.standards
sample = snakemake.wildcards.sample

## ====== Params ====== ##
library = snakemake.params.library
standard_set = snakemake.params.standard_set
parser_params = snakemake.params.parser_params

## ====== Output ====== ##
counts_csv = snakemake.output.counts
invalid_csv = snakemake.output.invalid
fates_csv = snakemake.output.fates

# Get the library barcodes
variants_df = pd.read_csv(variants)
# Get the neturalization standard barcodes
standards_df = pd.read_csv(standards)

# Check that the library for the sample matches the variants
if library not in variants_df["library"].unique():
    raise ValueError(f"{library=} not in {variants['library'].unique()=}")
# Check that the standard set for the sample matches the standards
if standard_set not in standards_df["standard_set"].unique():
    raise ValueError(f"{standard_set=} not in {standards['standard_set'].unique()=}")

# Get the valid barcodes for the library
library_barcodes = variants_df.query("library == @library")["barcode"].tolist()
# Get the valid barcodes for the standard set
standard_barcodes = standards_df.query("standard_set == @standard_set")[
    "barcode"
].tolist()

# Check that there are barcodes for the library and standard set
if len(library_barcodes) < 1:
    raise ValueError(f"no barcodes for {library=}")
if len(standard_barcodes) < 1:
    raise ValueError(f"no barcodes for {standard_set=}")
# Check that the barcodes are unique in both the library and standard set
if len(library_barcodes) != len(set(library_barcodes)):
    raise ValueError("There are non-unique barcodes in library")
if len(standard_barcodes) != len(set(standard_barcodes)):
    raise ValueError("There are non-unique barcodes in the neutralization standard set")

print(
    f"There are {len(library_barcodes)} valid barcodes for {library=} and there are {len(standard_barcodes)} valid barcodes for {standard_set=}"
)

# Join the library and standard barcodes
valid_barcodes = library_barcodes + standard_barcodes
if len(valid_barcodes) != len(set(valid_barcodes)):
    raise ValueError(
        "There is overlap between the library and neutralization standard barcodes"
    )
# Get the barcode length for the library
bclen = len(valid_barcodes[0])
# Check that all barcodes are the same length
if not all(bclen == len(bc) for bc in valid_barcodes):
    raise ValueError("Not all barcodes are of the same length")

valid_barcodes = set(valid_barcodes)

# Parse barcodes from the FASTQ files
print(f"Parsing barcodes from {fastq}")

# Check that the FASTQ files exist
if not os.path.isfile(fastq):
    raise ValueError(f"Cannot find FASTQ {fastq}")

# Initialize the parser
parser = dms_variants.illuminabarcodeparser.IlluminaBarcodeParser(
    bclen=bclen,
    **parser_params,
)
counts, fates = parser.parse(fastq)
counts["valid"] = counts["barcode"].isin(valid_barcodes)
print("Counts of valid and invalid barcodes:")
print(
    counts.groupby("valid").aggregate(
        n_barcodes=pd.NamedAgg("barcode", "nunique"),
        n_counts=pd.NamedAgg("count", "sum"),
    )
)

# Write the barcodes to their corresponding CSV
print(f"Writing valid barcode counts to {counts_csv}")
counts_valid = counts.query("valid").drop(columns="valid")
missing_valid_barcodes = sorted(list(valid_barcodes - set(counts_valid["barcode"])))
counts_valid = counts_valid = pd.concat(
    [
        counts_valid[["barcode", "count"]],
        pd.DataFrame({"barcode": missing_valid_barcodes, "count": 0}),
    ]
)
counts_valid.to_csv(counts_csv, index=False)

print(f"Writing invalid barcode counts to {invalid_csv}")
counts_invalid = (
    counts.query("not valid")
    .drop(columns="valid")
    .assign(library=library, sample=sample)
)
counts_invalid[["barcode", "count"]].to_csv(invalid_csv, index=False)

print(f"Writing barcode fates to {fates_csv}")
fates = pd.concat(
    [
        fates.query("fate not in ['valid barcode', 'invalid barcode']"),
        pd.DataFrame(
            {
                "fate": ["valid barcode", "invalid barcode"],
                "count": [
                    counts_valid["count"].sum(),
                    counts_invalid["count"].sum(),
                ],
            }
        ),
    ]
).assign(library=library, sample=sample)
fates[["fate", "count"]].to_csv(fates_csv, index=False)


print("\nFinished!")

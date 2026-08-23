"""Validate a viral library CSV against the checks configured for it in `config.yml`."""

import sys

import viral_library_validation

# `noqa: SIM115` as this log file must stay open for the life of the script
sys.stderr = sys.stdout = open(snakemake.log[0], "w")  # noqa: SIM115

try:
    report = viral_library_validation.report(
        snakemake.input.csv, snakemake.params.validations
    )
except viral_library_validation.LibraryError as e:
    # the report of a failed run is written too, so that it says what the library did get
    # through rather than only naming the first thing wrong with it
    print(e.report)
    raise

print(report)

with open(snakemake.output.validation, "w") as f:
    f.write(f"{report}\n")

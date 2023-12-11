"""Process QC for processing counts for all plates."""


import sys


sys.stderr = sys.stdout = log = open(snakemake.log[0], "w")

plates = snakemake.params.plates
qc_failures = snakemake.input.qc_failures
process_counts_htmls = snakemake.input.process_counts_htmls

assert len(plates) == len(qc_failures) == len(process_counts_htmls)

qc_summary = []

for plate, qc_failure, html in zip(
    plates,
    qc_failures,
    process_counts_htmls,
):
    with open(qc_failure) as f:
        failures = ["\t" + line.strip() for line in f if line.strip()]
    if failures:
        qc_summary.append(f"{plate} process_counts QC failures, see {html} for details")
        qc_summary += failures
    else:
        qc_summary.append(f"{plate} process_counts passed all QC")

with open(snakemake.output.qc_summary, "w") as f:
    f.write("\n".join(qc_summary))

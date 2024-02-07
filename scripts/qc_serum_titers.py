"""Process QC for serum titers."""

import sys


sys.stderr = sys.stdout = log = open(snakemake.log[0], "w")

sera = snakemake.params.sera
qc_failures = snakemake.input.qc_failures
serum_titers_htmls = snakemake.input.serum_titers_htmls

assert len(sera) == len(qc_failures) == len(serum_titers_htmls)

qc_summary = []

for serum, qc_failure, html in zip(
    sera,
    qc_failures,
    serum_titers_htmls,
):
    with open(qc_failure) as f:
        failures = ["\t" + line.strip() for line in f if line.strip()]
    if failures:
        qc_summary.append(f"{serum} serum_titers QC failures, see {html} for details")
        qc_summary += failures
    else:
        qc_summary.append(f"{serum} serum passed all QC")

with open(snakemake.output.qc_summary, "w") as f:
    f.write("\n".join(qc_summary))

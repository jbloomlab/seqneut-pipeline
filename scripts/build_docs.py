"""Implements ``snakemake`` rule to translate gene sequence."""


import os
import shutil
import sys

import markdown
import markdown.extensions.toc


sys.stderr = sys.stdout = log = open(snakemake.log[0], "w")

copied_files = {
    f: os.path.join(snakemake.output.docs, os.path.basename(f)) for f in snakemake.input
}
assert len(copied_files) == len(set(copied_files.values())) == len(snakemake.input)
os.makedirs(snakemake.output.docs, exist_ok=True)
for f in snakemake.input:
    shutil.copy(f, snakemake.output.docs)
assert all(os.path.isfile(f) for f in copied_files.values())

md_text = [
    snakemake.params.description,
    "",
    # table of contents: https://python-markdown.github.io/extensions/toc/
    "[TOC]",
    "",
    "## Plot of titers for all sera",
    f"[Interactive chart of titers]({copied_files[snakemake.input.titers_chart]})",
    "",
    "## Analyses of per-serum neutralization titers",
]

assert len(snakemake.params.sera) == len(snakemake.input.serum_titers_htmls)
for serum, f in zip(snakemake.params.sera, snakemake.input.serum_titers_htmls):
    md_text.append(f"  - [{serum}]({copied_files[f]})")

md_text += ["", "## Analyses of per-plate counts and curve fits"]
assert len(snakemake.params.plates) == len(snakemake.input.process_counts_htmls)
for plate, f in zip(snakemake.params.plates, snakemake.input.process_counts_htmls):
    md_text.append(f"  - [{plate}]({copied_files[f]})")

md_text = "\n".join(md_text)

print(f"Rendering the following markdown text:\n\n{md_text}\n\n")

html = markdown.markdown(
    md_text,
    extensions=[
        markdown.extensions.toc.TocExtension(
            title="Contents",
            toc_depth="2-2",
        ),
    ],
)

index_html = os.path.join(snakemake.output.docs, "index.html")
with open(index_html, "w") as f:
    f.write(html)

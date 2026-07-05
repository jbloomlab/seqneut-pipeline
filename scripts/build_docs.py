"""Implements ``snakemake`` rule to translate gene sequence."""

import os
import time
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

# get sera in each group, groups ordered by how many sera are in each
groups = {group: [] for group in set(snakemake.params.plates.values())}
for group, serum in snakemake.params.groups_sera:
    groups[group].append(serum)
group_order = [
    group
    for _, group in sorted(
        [(len(sera), group) for (group, sera) in groups.items()], reverse=True
    )
]
assert len(group_order) == len(set(group_order))

# map each per-group titers chart (named titers_{group}.html) to its group
titers_charts_by_group = {}
for f in snakemake.input.titers_charts:
    base = os.path.basename(f)
    assert base.startswith("titers_") and base.endswith(".html"), base
    titers_charts_by_group[base[len("titers_") : -len(".html")]] = f

md_text = [
    snakemake.params.description,
    "",
    f"Documentation of results rendered as of {time.asctime()}",
    "",
    # table of contents: https://python-markdown.github.io/extensions/toc/
    "[TOC]",
    "",
    "## Titers for all sera",
]
for group in group_order:
    md_text.append(
        f"- [Interactive chart of titers (`{group}`)]"
        f"({os.path.basename(copied_files[titers_charts_by_group[group]])})"
    )
md_text += [
    "",
    "## Per-serum neutralization titers",
]

assert len(snakemake.params.groups_sera) == len(snakemake.input.serum_titers_htmls)
for group in group_order:
    md_text.append(f"\n### {group}")
    for serum in groups[group]:
        f = snakemake.input.serum_titers_htmls[
            snakemake.params.groups_sera.index((group, serum))
        ]
        md_text.append(f"- [{serum}]({os.path.basename(copied_files[f])})")

md_text += ["", "## Per-plate counts and curve fits"]
assert len(snakemake.params.plates) == len(snakemake.input.process_plates_htmls)
for group in group_order:
    md_text.append(f"\n### {group}")
    for plate in [p for (p, g) in snakemake.params.plates.items() if g == group]:
        f = snakemake.input.process_plates_htmls[
            list(snakemake.params.plates).index(plate)
        ]
        md_text.append(f"- [{plate}]({os.path.basename(copied_files[f])})")

md_text += [
    "",
    "## Summary of data dropped during quality control",
    f"[Notebook summarizing QC drops]({os.path.basename(copied_files[snakemake.input.qc_drops_html])})",
]

for heading, heading_d in snakemake.params.add_htmls_to_docs.items():
    md_text += ["", f"## {heading}"]
    for name, value in heading_d.items():
        if isinstance(value, dict):
            md_text += ["", f"### {name}"]
            for sub_name, sub_fname in value.items():
                md_text.append(
                    f"- [{sub_name}]({os.path.basename(copied_files[sub_fname])})"
                )
        else:
            md_text.append(f"- [{name}]({os.path.basename(copied_files[value])})")

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

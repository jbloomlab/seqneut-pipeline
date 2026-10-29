"""Functions shared by the analyses of miscellaneous plates.

The miscellaneous plates are counted by ``miscellaneous_plate_count_barcodes``, which
writes a ``counts``, ``invalid``, and ``fates`` CSV for each well. The analyses built on
them read those per-well CSVs, set aside reads of barcodes that no analysis should count,
measure how far unmatched barcodes are from known ones, and place wells on the plate.

"""

import contextlib
import html
import re
from pathlib import Path

import numpy
import pandas
from seqneut_report import Report

# The plate layout that a well name is read as, so that `B6` becomes a row and a column and
# two wells can be told how far apart they are. A 96-well plate; a well outside it raises
# rather than being placed somewhere arbitrary.
PLATE_ROWS = "ABCDEFGH"
PLATE_COLUMNS = 12

# Material that belongs in another well of the same plate is attributed to how far away
# that well is, and distances further than this are reported as one group rather than each
# on its own: what matters is whether it came from a neighboring well, not whether it came
# from five wells away or from nine.
MAX_NAMED_DISTANCE = 2

# How many nucleotides an unmatched barcode can differ from a barcode and still be a misread
# of it. Fixed rather than configured, as it is part of what a misread means rather than an
# experimental choice; it does change which category a read falls in.
MAX_MISREAD_HAMMING = 1

# where material came from when it is not a distance on the same plate
ANOTHER_PLATE = "another plate here"
NOT_ASSAYED = "not assayed here"

# What each fate of a read means. Only the fates that have counts are explained in a
# report, and a fate the pipeline writes that is missing here raises: a fate nobody can
# interpret is worse than one nobody sees. "discarded barcode" is not one the pipeline
# writes, but split out of "invalid barcode" for the report.
FATE_EXPLANATIONS = {
    "valid barcode": (
        "the read contains a valid barcode present in either the viral library or the "
        "neutralization standard set."
    ),
    "invalid barcode": (
        "the read contains a barcode sequence not present in either the viral library or "
        "the neutralization standard set."
    ),
    "discarded barcode": (
        "the read contains a barcode present in the discarded barcode set."
    ),
    "low quality barcode": (
        "the read contains poor-quality nucleotide sequencing within the barcode region."
    ),
    "unparseable barcode": (
        "the read contains the sequence flanking the barcode, but the barcode could not "
        "be read from its expected position."
    ),
    "invalid outer flank": (
        "the read contains a different outer flanking sequence than specified. This outer "
        "barcode is used to separate samples on a shared sequencing run, so reads "
        "belonging to other samples on a multiplexed run are counted here and are expected "
        "in this sequencing regime."
    ),
    "read too short": "the read is too short to contain the barcode and its flanks.",
    "failed chastity filter": (
        "the read did not pass the Illumina chastity filter for signal purity."
    ),
}


def well_position(well):
    """Get the (row index, column) of `well`, a well name such as ``B6``."""
    row, column = well[:1], well[1:]
    if row not in PLATE_ROWS or not column.isdigit():
        raise ValueError(
            f"well {well!r} is not a row letter followed by a column number, which a "
            f"{len(PLATE_ROWS)} x {PLATE_COLUMNS} plate needs it to be"
        )
    if not 1 <= int(column) <= PLATE_COLUMNS:
        raise ValueError(f"well {well!r} is outside a {PLATE_COLUMNS}-column plate")
    return PLATE_ROWS.index(row), int(column)


def wells_apart(well_1, well_2):
    """Get how many wells apart `well_1` and `well_2` are on the plate.

    A step onto any of the eight surrounding wells counts as one, diagonals included, so
    two wells one apart are the wells that touch. That is the right measure for liquid
    carried between wells, which has no reason to prefer a row to a column or either to a
    diagonal.

    """
    (row_1, column_1), (row_2, column_2) = well_position(well_1), well_position(well_2)
    return max(abs(row_1 - row_2), abs(column_1 - column_2))


def distance_origin(distance):
    """The origin of material from a well `distance` wells away on the same plate."""
    if distance <= MAX_NAMED_DISTANCE:
        return f"{distance} well{'s' if distance > 1 else ''} away"
    return f"over {MAX_NAMED_DISTANCE} wells away"


def all_origins():
    """Every origin `trace_origin` can give, nearest first."""
    return [
        *(distance_origin(d) for d in range(1, MAX_NAMED_DISTANCE + 2)),
        ANOTHER_PLATE,
        NOT_ASSAYED,
    ]


def trace_origin(plate, well, homes):
    """Where material found in `well` of `plate` came from.

    Args:
        `plate`, `well` (str)
            Where the material was found.
        `homes` (list)
            The (plate, well) pairs the material belongs in; empty for material that no
            well was infected with.

    Returns:
        `(origin, source_well, wells_away, came_from)`: the origin, one of `all_origins`;
        the well the material belongs in, or null if it belongs in none; how many wells
        away that is on the same plate, or null if it is not on this plate; and a short
        description of the source for tables. Of several homes on this plate the source
        is the nearest, as the likeliest to have spilled, and the first given of those
        equally near; of homes only on other plates it is the first given.

    """
    same_plate = [w for (p, w) in homes if p == plate]
    if same_plate:
        source = min(same_plate, key=lambda w: wells_apart(well, w))
        distance = wells_apart(well, source)
        origin = distance_origin(distance)
        return origin, source, distance, f"{source}, {origin}"
    if homes:
        home_plate, home_well = homes[0]
        return ANOTHER_PLATE, home_well, None, f"{home_well} of {home_plate}"
    return NOT_ASSAYED, None, None, NOT_ASSAYED


def check_discarded_barcodes(discarded_barcodes, counted_barcodes):
    """Check barcodes to discard against the barcodes that are counted.

    Each discarded barcode is given with the reason it is discarded, which the reports show.
    A discarded barcode must be in neither the viral library nor the neutralization standard
    set, whose reads are counted before any could be discarded; it can be any distance from
    them otherwise, as only its exact reads are discarded.

    Args:
        `discarded_barcodes` (dict)
            The barcodes to discard, each to the reason it is discarded, as configured.
        `counted_barcodes` (list)
            The viral library and neutralization standard barcodes.

    Returns:
        The discarded barcodes upper-cased, as the counted barcodes are compared, each to
        its reason.

    """
    if not isinstance(discarded_barcodes, dict):
        raise TypeError(
            "'discarded_barcodes' must map each barcode to the reason it is discarded, "
            f"not be {discarded_barcodes!r}"
        )
    counted = [bc.upper() for bc in counted_barcodes]
    checked = {}
    for barcode, reason in discarded_barcodes.items():
        if not isinstance(barcode, str) or len(barcode) != len(counted[0]):
            raise ValueError(
                f"discarded barcode {barcode!r} is not a barcode of the length of those it "
                f"is counted against, {len(counted[0])}"
            )
        barcode = barcode.upper()
        if set(barcode) - set("ACGT"):
            raise ValueError(
                f"discarded barcode {barcode} is not made of A, C, G, and T"
            )
        if not isinstance(reason, str) or not reason.strip():
            raise ValueError(
                f"discarded barcode {barcode} needs the reason it is discarded, not "
                f"{reason!r}"
            )
        if barcode in checked:
            raise ValueError(f"discarded barcode {barcode} is given more than once")
        if barcode in counted:
            raise ValueError(
                f"discarded barcode {barcode} is in the viral library or neutralization "
                "standard set, so its reads are counted rather than discarded"
            )
        checked[barcode] = reason.strip()
    return checked


def plate_and_well(csv_file, suffix):
    """Get the (plate, well) that a per-well `suffix` CSV is for, from its path.

    The counting rule writes ``results/miscellaneous_plates/<plate>/<well><suffix>.csv``.

    """
    path = Path(csv_file)
    return path.parent.name, path.stem.removesuffix(suffix)


def read_per_well_csvs(csv_files, suffix, required_columns, plate_wells):
    """Read per-well CSVs of several plates into one frame labeled by plate and well.

    Args:
        `csv_files` (list)
            The CSVs, one per well, as written by the counting rule.
        `suffix` (str)
            What follows the well in each file's name, such as ``_counts``.
        `required_columns` (list)
            Columns each file must have. Checked on each file as it is read rather than on
            the concatenated frame, where a column missing from one file would survive as
            nulls and then read as a measurement of zero.
        `plate_wells` (iterable)
            The (plate, well) pairs that must be read, exactly. Checked from the file
            *paths* rather than from the rows read, since a CSV of this kind can
            legitimately be empty: a well with no unmatched barcodes has nothing to list.

    Returns:
        `pandas.DataFrame` of the rows of every file, with columns ``plate`` and ``well``.

    """
    read_for = [plate_and_well(csv_file, suffix) for csv_file in csv_files]
    if len(set(read_for)) != len(read_for):
        raise ValueError(f"two of the {suffix!r} CSVs are for the same plate and well")
    expected = set(map(tuple, plate_wells))
    if set(read_for) != expected:
        raise ValueError(
            f"the {suffix!r} CSVs are for wells the samples CSVs do not list, or the "
            f"other way round: {sorted(set(read_for) ^ expected)}"
        )
    frames = []
    for csv_file, (plate, well) in zip(csv_files, read_for, strict=True):
        frame = pandas.read_csv(csv_file)
        missing = set(required_columns) - set(frame.columns)
        if missing:
            raise ValueError(f"{csv_file} lacks the columns {sorted(missing)}")
        frames.append(frame.assign(plate=plate, well=well))
    # the empty files are left out where there are others, as `pandas` would otherwise
    # take their columns' dtypes into account when typing the concatenated frame
    return pandas.concat([f for f in frames if len(f)] or frames, ignore_index=True)


def reads_per_well(frame, well_index):
    """Sum the counts of `frame` in each well, a well without rows summing to zero.

    Args:
        `frame` (pandas.DataFrame)
            Rows with ``plate``, ``well``, and ``count`` columns.
        `well_index` (pandas.MultiIndex)
            The ``(plate, well)`` of every well, in the order to return them.

    Returns:
        pandas.Series of the reads of each well of `well_index`, indexed from zero.

    Reads are counts, so the result is an integer however `frame` is typed: summing a
    frame of no rows gives a column of Python objects, which carries through arithmetic
    into anything derived from it and is written to a CSV at full precision.

    """
    return (
        frame.groupby(["plate", "well"])["count"]
        .sum()
        .reindex(well_index, fill_value=0)
        .astype(int)
        .reset_index(drop=True)
    )


def as_barcode_array(barcodes):
    """Get `barcodes` as a 2D array of one nucleotide per cell, for comparing them."""
    lengths = {len(barcode) for barcode in barcodes}
    if len(lengths) != 1:
        raise ValueError(f"the barcodes are not all one length: {sorted(lengths)}")
    return numpy.frombuffer("".join(barcodes).encode(), dtype="S1").reshape(
        len(barcodes), -1
    )


def min_hamming(barcodes, targets):
    """Get the smallest Hamming distance from each of `barcodes` to any of `targets`."""
    barcodes, targets = list(barcodes), list(targets)
    if not barcodes:
        return numpy.array([], dtype=int)
    if not targets:
        raise ValueError("no target barcodes to measure a distance to")
    if len(barcodes[0]) != len(targets[0]):
        raise ValueError("the barcodes and the targets are of different lengths")
    return (
        (
            as_barcode_array(barcodes)[:, None, :]
            != as_barcode_array(targets)[None, :, :]
        )
        .sum(axis=2)
        .min(axis=1)
    )


def split_discarded(invalid, discarded_barcodes):
    """Split unmatched barcodes into those kept and those of discarded barcodes.

    A discarded barcode is one physically in every well that no analysis should count, such
    as a spike-in barcode left out of the neutralization standard set because it sequences
    badly. It is in neither whitelist, so its reads are among the unmatched ones. Only an
    exact match is set aside, as only an exact match of a counted barcode is counted; a
    misread of a discarded barcode stays among the unmatched reads.

    Args:
        `invalid` (pandas.DataFrame)
            Unmatched barcodes, with a ``barcode`` column.
        `discarded_barcodes` (list or dict)
            The barcodes to discard, possibly empty: a dict, as from
            `check_discarded_barcodes`, is taken for its keys.

    Returns:
        `(kept, discarded)`, the rows of `invalid` split in two.

    """
    if not discarded_barcodes or not len(invalid):
        return invalid, invalid.iloc[:0]
    is_discarded = invalid["barcode"].isin(list(discarded_barcodes))
    return invalid[~is_discarded], invalid[is_discarded]


def check_fates_add_up(fates, tallies):
    """Check that reads tallied from the counts add back up to the per-well fates.

    The tallies are built from the ``counts`` and ``invalid`` CSVs, so they have to agree
    with what the ``fates`` CSV says became of the reads; a tally that double-counted or
    dropped a read would otherwise go unnoticed.

    Args:
        `fates` (pandas.DataFrame)
            Columns ``plate``, ``well``, ``fate``, ``count``.
        `tallies` (dict)
            Maps a fate to a `pandas.Series` of reads indexed by (plate, well), which must
            equal that fate's count in each well, a fate absent from a well counting as 0.

    """
    for fate, tally in tallies.items():
        from_fates = (
            fates[fates["fate"] == fate]
            .groupby(["plate", "well"])["count"]
            .sum()
            .reindex(tally.index, fill_value=0)
        )
        disagree = from_fates[from_fates != tally]
        if len(disagree):
            raise ValueError(
                f"the reads tallied do not match the {fate!r} fate in the wells "
                f"{list(disagree.index)}: {tally[disagree.index].tolist()} tallied, "
                f"{disagree.tolist()} in the fates"
            )


def table_html(df, index=False, marked=None, marked_color="#cc00cc", group=None):
    """HTML for `df` in a scrolling box, with floats to four significant digits.

    The arguments are those of `PlateReport.table`.

    """
    text = df.map(
        lambda v: (
            ""
            if pandas.isna(v)
            else html.escape(f"{v:.4g}" if isinstance(v, float) else str(v))
        )
    )
    if marked is not None and marked.any():
        text.loc[marked.to_numpy()] = text.loc[marked.to_numpy()].map(
            lambda v: f'<span style="color: {marked_color}">{v}</span>'
        )
    if group is None:
        table = text.rename(columns=_header_html).to_html(
            index=index, border=0, escape=False
        )
    else:
        if index:
            raise ValueError("a table with a group of columns cannot show its index")
        table = _grouped_table_html(text, *group)
    return '<div class="scroll-box" style="width: fit-content">' + table + "</div>"


def _header_html(name):
    """A column header that wraps rather than widening a column of shorter entries.

    It wraps between words only, its hyphens being made non-breaking, so that a term such
    as "non-neut-standard" stays on one line.

    """
    name = html.escape(str(name)).replace("-", "&#8209;")
    return f'<div class="head">{name}</div>'


def _grouped_table_html(text, label, columns):
    """HTML for `text` with `columns` under a header `label`, the rest spanning rows.

    `columns` is a list of columns, or a dict mapping each to its sub-header, which lets a
    sub-column share a name with a column outside the group.

    """
    if not isinstance(columns, dict):
        columns = {c: c for c in columns}
    keys = [c for c in text.columns if c not in columns]
    if not keys or set(columns) - set(text.columns):
        raise ValueError(
            f"cannot group {columns} out of the columns {list(text.columns)}"
        )
    # Stuck as a whole rather than cell by cell, as its second row sits below a first whose
    # height depends on how its headers wrap.
    head = (
        "<tr>"
        + "".join(f'<th rowspan="2">{_header_html(c)}</th>' for c in keys)
        + f'<th colspan="{len(columns)}" style="text-align: center">'
        + f"{html.escape(label)}</th></tr><tr>"
        + "".join(f"<th>{_header_html(c)}</th>" for c in columns.values())
        + "</tr>"
    )
    body = []
    for _, rows in text.groupby(keys, sort=False):
        for i, row in enumerate(rows.itertuples(index=False)):
            row = dict(zip(text.columns, row))
            spanned = "".join(
                f'<td rowspan="{len(rows)}" style="vertical-align: top">{row[c]}</td>'
                for c in keys
            )
            body.append(
                "<tr>"
                + (spanned if i == 0 else "")
                + "".join(f"<td>{row[c]}</td>" for c in columns)
                + "</tr>"
            )
    return (
        f'<table class="dataframe"><thead class="grouped">{head}</thead>'
        f'<tbody>{"".join(body)}</tbody></table>'
    )


# Each warning is a line with a button that opens its table, closed until clicked.
_WARNINGS_CSS = """
details.warning {
  background: #fff8e1; border-left: 4px solid #f0ad4e; margin: 0.4em 0;
  padding: 0.4em 0.8em;
}
details.warning summary { cursor: pointer; list-style: none; }
details.warning .scroll-box { background: white; }
details.warning summary::-webkit-details-marker { display: none; }
details.warning .toggle {
  border: 1px solid #0366d6; border-radius: 1em; color: #0366d6; font-size: 0.85em;
  margin-left: 0.5em; padding: 0.05em 0.6em; white-space: nowrap;
}
details.warning summary:hover .toggle { background: #0366d6; color: white; }
details.warning .toggle::before { content: "\\25B8  "; }
details.warning[open] .toggle::before { content: "\\25BE  "; }
details.warning .hide, details.warning[open] .show { display: none; }
details.warning[open] .hide { display: inline; }
"""

# Lines between a table's columns; a row's last cell is always its last column, as the cells
# a grouped table spans over rows come first. A header's lines are drawn as shadows, as a
# sticky cell's borders are not drawn in a collapsed table. Headers sit at the bottom of
# their row, and wrap past a width cap, so a column is as wide as the wider of its entries
# and its wrapped header. A grouped table's header sticks as a whole, its rows' heights
# depending on how they wrap, with a line under the header spanning its group.
#
# A table or chart wider than the text runs on past its right edge up to the window's,
# scrolling only beyond that, the body's padding being kept on the right. Beside the sidebar
# the body's left edge is the sidebar's width in from the window's; on a narrow window,
# where the sidebar is above the body and the body is centered, it is half of what the
# window has beyond the text's width, `(100vw - 100%) / 2`.
_TABLE_CSS = """
.scroll-box td { border-right: 1px solid #ddd; }
.scroll-box td:last-child { border-right: none; }
.scroll-box th {
  border-right: none; box-shadow: inset -1px 0 #ddd, inset 0 -1px #ddd; vertical-align: bottom;
}
.scroll-box th:last-child { box-shadow: inset 0 -1px #ddd; }
.scroll-box th .head {
  margin-left: auto; max-width: 7em; min-width: min-content; white-space: normal;
}
.scroll-box thead.grouped { position: sticky; top: 0; z-index: 1; }
.scroll-box thead.grouped th { position: static; }
.scroll-box, .chart-block {
  max-width: calc(100vw - var(--sidebar-width) - 3em); width: fit-content;
}
@media (max-width: 60em) {
  .scroll-box, .chart-block { max-width: calc((100vw + 100%) / 2 - 1.5em); }
}
"""

# A sidebar of the report's sections, fixed at the left of the window and with the section
# in view picked out, which sits above the body instead on a narrow window. It links back
# to the docs' index, which `build_docs` puts beside every report it copies there.
_SIDEBAR_CSS = """
:root { --sidebar-width: 15rem; }
body { margin: 0 auto 0 var(--sidebar-width); }
nav.sidebar {
  background: #fafafa; border-right: 1px solid #ddd; bottom: 0; box-sizing: border-box;
  font-size: 0.9rem; left: 0; overflow-y: auto; padding: 1.5rem 1rem; position: fixed;
  top: 0; width: var(--sidebar-width);
}
nav.sidebar a.home {
  border: 1px solid #ccc; border-radius: 4px; display: inline-block; margin-bottom: 1.2em;
  padding: 0.25em 0.7em;
}
nav.sidebar a.home:hover { background: #eee; text-decoration: none; }
nav.sidebar .contents {
  border-bottom: 1px solid #ddd; color: #666; font-size: 0.85em; font-weight: 600;
  letter-spacing: 0.06em; margin: 0 0 0.6em; padding-bottom: 0.2em;
  text-transform: uppercase;
}
nav.sidebar ul { list-style: none; margin: 0.2em 0; padding-left: 0; }
nav.sidebar a { color: #245; display: block; padding: 0.16em 0; text-decoration: none; }
nav.sidebar a:hover { text-decoration: underline; }
nav.sidebar a.active { color: #111; font-weight: 600; }
@media (max-width: 60em) {
  :root { --sidebar-width: 0px; }
  body { margin: 0 auto; }
  nav.sidebar {
    border: 1px solid #ddd; border-radius: 4px; max-height: 40vh; position: static;
    width: auto;
  }
}
"""

# Picks out the sidebar's entry for the section being read: the last whose heading is above
# a fifth of the way down the window, or the first before any is, or the last once the page
# is scrolled to its foot, which a short last section's heading may never get that high.
# Only the links to sections are followed, not the one back to the index.
_SIDEBAR_JS = """
const links = [...document.querySelectorAll("nav.sidebar a[href^='#']")];
const headings = links.map((a) => document.getElementById(a.hash.slice(1)));
function pickSection() {
  const atFoot = window.innerHeight + window.scrollY
    >= document.documentElement.scrollHeight - 2;
  let current = 0;
  headings.forEach((h, i) => {
    if (h.getBoundingClientRect().top <= window.innerHeight / 5) current = i;
  });
  if (atFoot) current = headings.length - 1;
  links.forEach((a, i) => a.classList.toggle("active", i === current));
}
window.addEventListener("scroll", pickSection, { passive: true });
window.addEventListener("resize", pickSection);
window.addEventListener("load", pickSection);
pickSection();
"""

# The headings the sidebar lists: the report's sections, but not the warnings or anything
# else headed within one.
_SIDEBAR_HEADING_RE = re.compile(r"<h2>(.*?)</h2>", re.DOTALL)


def _with_sidebar(blocks):
    """`blocks` with an id on each section's heading, after a sidebar linking to them."""
    entries = []

    def add_id(match):
        inner = match.group(1)
        text = html.unescape(re.sub(r"<[^>]+>", "", inner))
        base = re.sub(r"[^a-z0-9]+", "-", text.lower()).strip("-") or "section"
        # a heading used twice gets the next free id
        taken = {anchor for (anchor, _) in entries}
        anchor, n = base, 1
        while anchor in taken:
            n += 1
            anchor = f"{base}-{n}"
        entries.append((anchor, inner))
        return f'<h2 id="{anchor}">{inner}</h2>'

    blocks = [_SIDEBAR_HEADING_RE.sub(add_id, block) for block in blocks]
    toc = "".join(
        f'<li><a href="#{anchor}">{inner}</a></li>' for (anchor, inner) in entries
    )
    sidebar = (
        f'<style>{_SIDEBAR_CSS}</style><nav class="sidebar">'
        '<a class="home" href="index.html">&larr; All results</a>'
        f'<p class="contents">Contents</p><ul>{toc}</ul></nav>'
    )
    return [
        sidebar,
        *blocks,
        f"<script>{_SIDEBAR_JS}</script>",
    ]


# A notes block is a link-like line that opens the text under it, closed until clicked.
_NOTES_CSS = """
details.notes { border-left: 4px solid #ddd; margin: 0.6em 0; padding: 0.2em 0.8em; }
details.notes summary { color: #0366d6; cursor: pointer; list-style: none; }
details.notes summary::-webkit-details-marker { display: none; }
details.notes summary::before { content: "\\25B8  "; }
details.notes[open] summary::before { content: "\\25BE  "; }
"""


def _logged_warning(text, details, shows="well", note=None):
    """Print a warning to the log, and return it as `_warnings_html` takes it."""
    print(
        f"WARNING: {text}\n"
        + (f"{note}\n" if note else "")
        + details.to_string(index=False)
    )
    return (text, details, shows, note)


def _warnings_html(warnings):
    """HTML for `warnings` under a heading, each a line that opens its table."""
    return f"<style>{_WARNINGS_CSS}</style><h3>Warnings</h3>" + "".join(
        f'<details class="warning"><summary><b>{len(details)}</b> {noun} '
        f'{html.escape(text)} <span class="toggle"><span class="show">show'
        f'</span><span class="hide">hide</span> {noun}</span></summary>'
        + (f"<p>{html.escape(note)}</p>" if note else "")
        + f"{table_html(details)}</details>"
        for (text, details, shows, note) in warnings
        for noun in [shows if len(details) == 1 else shows + "s"]
    )


class PlateReport(Report):
    """A `Report` with its warnings gathered in one slot, formatted tables, and a sidebar.

    Warnings are found as the analysis runs, but are wanted where they will be read first, so
    each is logged as it is added and they are all written in the slot `reserve_warnings`
    leaves for them.

    """

    def __init__(self, title):
        super().__init__(title)
        self._blocks.append(f"<style>{_TABLE_CSS}</style>")
        self._warnings = []
        self._section_warnings = []
        self._warnings_index = None
        self._notes_styled = False

    def reserve_warnings(self):
        """Leave the slot the warnings are written into, at the current position."""
        if self._warnings_index is not None:
            raise ValueError("the warnings slot is already reserved")
        self._warnings_index = len(self._blocks)
        self._blocks.append("")

    def warning(self, text, details, shows="well", note=None):
        """Record a warning, printing it to the log as well.

        Args:
            `text` (str)
                What the warning is about, following the count of rows of `details`
                and `shows`: "with no reads" to read as "2 wells with no reads".
            `details` (pandas.DataFrame)
                What it applies to, shown in a dropdown beside `text`.
            `shows` (str)
                What a row of `details` is.
            `note` (str or None)
                A line shown in the dropdown above `details`.

        """
        self._warnings.append(_logged_warning(text, details, shows, note))

    def section_warnings(self, warnings):
        """Append `warnings` here, under their own heading, rather than in the slot.

        For warnings about one section, which are read with it. Each is a dict of the
        arguments of `warning`, and nothing is appended if there are none.

        """
        warnings = [_logged_warning(**w) for w in warnings]
        if warnings:
            self._section_warnings += warnings
            self._blocks.append(_warnings_html(warnings))

    @property
    def warnings(self):
        """The warnings recorded so far, those of the slot and then those of sections."""
        return self._warnings + self._section_warnings

    def table(self, df, index=False, marked=None, marked_color="#cc00cc", group=None):
        """Append `df` as a table, with floats to four significant digits.

        Args:
            `df` (pandas.DataFrame)
                The frame to draw.
            `index` (bool)
                Write the index.
            `marked` (pandas.Series or None)
                Boolean mask of rows whose text is written in `marked_color`.
            `marked_color` (str)
                The color of marked rows.
            `group` (tuple or None)
                `(label, columns)`: draw `columns` as sub-columns under the header
                `label`, where `columns` is a list, or a dict mapping each column to its
                sub-header; and draw the rows that agree on all the other columns together,
                where they first appear, with those columns as one cell spanning them.

        """
        self._blocks.append(table_html(df, index, marked, marked_color, group))

    @contextlib.contextmanager
    def notes(self, summary):
        """Hide whatever is appended within the block behind `summary`, closed by default."""
        style = "" if self._notes_styled else f"<style>{_NOTES_CSS}</style>"
        self._notes_styled = True
        self._blocks.append(
            f'{style}<details class="notes"><summary>{html.escape(summary)}</summary>'
        )
        yield
        self._blocks.append("</details>")

    def stat_tiles(self, tiles, caption):
        """Append a row of tiles, each a value over its label, with `caption` below.

        Args:
            `tiles` (dict)
                The value of each tile, keyed by its label.
            `caption` (str)
                A line of text under the tiles.

        """
        self._blocks.append(
            '<div style="display: flex; flex-wrap: wrap; gap: 0.6em; margin: 1em 0 0.4em">'
            + "".join(
                '<div style="border: 1px solid #ddd; border-radius: 6px; min-width: 6em; '
                'padding: 0.5em 1em; text-align: center">'
                f'<div style="font-size: 1.6em; font-weight: bold">'
                f"{html.escape(str(value))}</div>"
                f'<div style="color: #555; font-size: 0.85em">{html.escape(label)}</div>'
                "</div>"
                for (label, value) in tiles.items()
            )
            + f'</div><p style="color: #555">{html.escape(caption)}</p>'
        )

    def write(self, path):
        """Write the report to `path`, with the warnings in their slot and a sidebar."""
        if self._warnings:
            if self._warnings_index is None:
                raise ValueError("there are warnings but no slot was reserved for them")
            self._blocks[self._warnings_index] = _warnings_html(self._warnings)
        blocks = self._blocks
        self._blocks = _with_sidebar(blocks)
        try:
            super().write(path)
        finally:
            self._blocks = blocks

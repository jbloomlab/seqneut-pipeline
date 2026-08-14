"""Build the HTML report that a pipeline script writes alongside its data outputs.

A report is a sequence of blocks appended in the order the script computes them, then
written to one self-contained HTML file. The blocks are markdown prose, ``altair``
charts, ``pandas`` data frames, and ``matplotlib`` figures.

``altair`` builds the whole chart embedding itself, so nothing here writes a script tag
or names a ``vega`` version: ``Chart.to_html`` with ``template="universal"`` returns a
fragment that loads ``vega``, ``vega-lite``, and ``vega-embed`` on demand, guarded by a
global recording which versions are already present, so the libraries load once however
many charts a page holds.

"""

import io
import textwrap

import markdown
from neutcurve.fig_utils import fig_html
from ruamel.yaml import YAML

# Blocks that can be taller than the window are put in a box of this height that scrolls
# rather than pushing the rest of the page down. `neutcurve.fig_utils.fig_html` uses the
# same height for the figures it wraps.
SCROLL_BOX_HEIGHT = "80vh"

# `KaTeX` renders the LaTeX in the markdown, which `markdown` passes through untouched.
# It is restricted to the prose blocks so that it cannot rewrite anything a chart drew.
KATEX_VERSION = "0.16.11"

_PAGE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@{katex}/dist/katex.min.css">
<style>
body {{
  font-family: system-ui, -apple-system, "Segoe UI", Helvetica, Arial, sans-serif;
  line-height: 1.5;
  margin: 0 auto;
  max-width: 80em;
  padding: 1.5em;
}}
h1, h2, h3 {{ line-height: 1.25; }}
code, pre {{ font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; }}
code {{ background: #f4f4f4; border-radius: 3px; padding: 0.1em 0.3em; }}
pre {{ background: #f4f4f4; border-radius: 4px; overflow-x: auto; padding: 0.6em 0.8em; }}
pre code {{ background: none; padding: 0; }}
/* a chart is drawn at its natural size, so a wide one scrolls rather than being clipped */
.chart-block {{ margin: 1em 0; overflow-x: auto; }}
.scroll-box {{
  border: 1px solid #ddd;
  border-radius: 4px;
  margin: 1em 0;
  max-height: {scroll_height};
  overflow: auto;
}}
.scroll-box table {{ border-collapse: collapse; font-size: 0.85em; }}
.scroll-box th, .scroll-box td {{
  border-bottom: 1px solid #eee;
  padding: 0.25em 0.6em;
  text-align: right;
  white-space: nowrap;
}}
.scroll-box th {{ background: #fafafa; position: sticky; top: 0; }}
</style>
</head>
<body>
<h1>{title}</h1>
{body}
<script defer src="https://cdn.jsdelivr.net/npm/katex@{katex}/dist/katex.min.js"></script>
<script defer src="https://cdn.jsdelivr.net/npm/katex@{katex}/dist/contrib/auto-render.min.js"
  onload="document.querySelectorAll('.md-block').forEach(
    (el) => renderMathInElement(el, {{delimiters: [
      {{left: '$$', right: '$$', display: true}},
      {{left: '$', right: '$', display: false}}
    ]}})
  )"></script>
</body>
</html>
"""


class Report:
    """An HTML report built up a block at a time.

    Args:
        `title` (str)
            Title of the page, also written as its top-level heading.

    """

    def __init__(self, title):
        self.title = title
        self._blocks = []
        self._n_charts = 0  # gives each chart a unique `altair` `output_div`

    def md(self, text):
        """Append a block of markdown prose, which may contain LaTeX.

        The text is dedented first, so that it can be written as an indented triple-
        quoted string without markdown reading the indentation as a code block.

        """
        html = markdown.markdown(
            textwrap.dedent(text),
            extensions=["fenced_code", "tables", "sane_lists"],
        )
        self._blocks.append(f'<div class="md-block">{html}</div>')

    def yaml(self, obj):
        """Append `obj` as a fenced YAML block."""
        buffer = io.StringIO()
        YAML(typ="rt").dump(obj, buffer)
        self.md(f"```yaml\n{buffer.getvalue()}```")

    def chart(self, chart):
        """Append an `altair` chart."""
        self._n_charts += 1
        fragment = chart.to_html(
            fullhtml=False,
            output_div=f"chart_{self._n_charts}",
            template="universal",
        )
        self._blocks.append(f'<div class="chart-block">{fragment}</div>')

    def table(self, df, index=True):
        """Append a `pandas` data frame as a table in a scrollable box.

        Args:
            `df` (pandas.DataFrame)
                The frame to draw.
            `index` (bool)
                Write the index. Pass `False` where it is a meaningless range.

        """
        self._blocks.append(
            f'<div class="scroll-box">{df.to_html(index=index, border=0)}</div>'
        )

    def figure(self, fig, display_method):
        """Append a `matplotlib` figure, rendered as by `neutcurve.fig_utils.fig_html`."""
        self._blocks.append(fig_html(fig, display_method))

    def write(self, path):
        """Write the report to `path`."""
        with open(path, "w") as f:
            f.write(
                _PAGE.format(
                    title=self.title,
                    katex=KATEX_VERSION,
                    scroll_height=SCROLL_BOX_HEIGHT,
                    body="\n".join(self._blocks),
                )
            )

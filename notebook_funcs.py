"""Functions used inside some Jupyter notebooks."""

import base64
import io

import IPython.display

import PIL


def display_curve_fig(fig, curve_display_method):
    """How to display a ``matplotlib`` curve figure.

    Displays a ``matplotlib`` figure in one of several possible ways in Jupyter notebook.

    Parameters
    -----------
    fig : matplotlib.Figure
        The figure to display
    curve_display_method : {"inline", "png8", "no_display"}
        If "inline" just displays the figure inline. If "png8" then converts to
        ``*.png8`` and displays; will be lower resolution but smaller. If
        "no_display" then no curve shown.

    """
    if curve_display_method == "inline":
        display(fig)

    elif curve_display_method == "png8":
        # show smaller image using PNG-8
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=80)
        im = PIL.Image.open(io.BytesIO(buf.getvalue())).quantize(
            colors=48, dither=PIL.Image.NONE
        )
        out = io.BytesIO()
        im.save(out, format="PNG", optimize=True, pnginfo=PIL.PngImagePlugin.PngInfo())
        data = base64.b64encode(out.getvalue()).decode("ascii")
        display(
            IPython.display.HTML(
                f'<img src="data:image/png;base64,{data}" alt="figure">'
            )
        )

    elif curve_display_method != "no_display":
        raise ValueError(f"invalid {curve_display_method=}")

"""Shared matplotlib style for DeepUnitMatch paper figures.

Makes saved .svg files Inkscape-friendly:
  - svg.fonttype = "none": text is written as real <text> elements (editable
    in Inkscape) instead of being converted to glyph outlines (paths).
  - Arial as the font, so what Inkscape shows matches the rest of the figure.
  - Figure sizes in mm (mm_figsize), so panels paste in at their final size.

Usage:
    import paper_style
    paper_style.apply()
    fig, ax = plt.subplots(figsize=paper_style.mm_figsize(45, 30))
    ...
    paper_style.save_svg(fig, "panel.svg")
"""

import os

import matplotlib as mpl

FONT_SIZE = 7  # pt; Nature Methods asks for 5-7 pt

# Colours used across Figure 2
COLOURS = {
    "half1": "#000000",        # first half / reference unit
    "half2": "#808080",        # second half of the same unit
    "diff_within": "#EE6C21",  # different unit within session
    "match_across": "#209120", # match across sessions
    "DUM": "#E41A1C",
    "UM": "#1F5FBF",
}

MM_PER_INCH = 25.4


def mm_figsize(width_mm, height_mm):
    return (width_mm / MM_PER_INCH, height_mm / MM_PER_INCH)


def apply(font_size=FONT_SIZE):
    """Set rcParams for editable-text svg output. font_size=None leaves sizes alone."""
    params = {
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "svg.hashsalt": "deepunitmatch",  # deterministic ids between re-saves
    }
    if font_size is not None:
        params.update({
            "font.size": font_size,
            "axes.labelsize": font_size,
            "axes.titlesize": font_size,
            "xtick.labelsize": font_size,
            "ytick.labelsize": font_size,
            "legend.fontsize": font_size,
            "axes.linewidth": 0.6,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.minor.width": 0.4,
            "ytick.minor.width": 0.4,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "xtick.minor.size": 1.5,
            "ytick.minor.size": 1.5,
            "lines.linewidth": 1.0,
        })
    mpl.rcParams.update(params)


def offset_spines(ax, points=4):
    ax.spines[["right", "top"]].set_visible(False)
    ax.spines["left"].set_position(("outward", points))
    ax.spines["bottom"].set_position(("outward", points))


def save_svg(fig, path, png=False):
    """Save fig as editable svg (and optionally a png preview next to it)."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with mpl.rc_context({"svg.fonttype": "none"}):
        fig.savefig(path, format="svg", bbox_inches="tight", transparent=True,
                    metadata={"Date": None})
    if png:
        fig.savefig(os.path.splitext(path)[0] + ".png", dpi=300, bbox_inches="tight")

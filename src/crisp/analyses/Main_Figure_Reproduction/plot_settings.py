"""Matplotlib settings for the supplied quantitative figures."""

import matplotlib as mpl


def apply_figure_style(font="Arial", sizes=(8, 7, 6)):
    """Apply explicit font, line and export settings before creating axes."""
    base, annotation, tick = sizes
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": [font, "DejaVu Sans", "sans-serif"],
        "font.size": base,
        "axes.titlesize": base,
        "axes.titlelocation": "left",
        "axes.labelsize": base,
        "legend.fontsize": annotation,
        "xtick.labelsize": tick,
        "ytick.labelsize": tick,
        "axes.linewidth": 0.6,
        "patch.linewidth": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "legend.frameon": False,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.dpi": 300,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })

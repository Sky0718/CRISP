"""Shared vector and preview export for the CRISP scientific figures."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager


def apply_figure_style():
    """Select an installed sans-serif family and editable vector text."""
    families = {font.name for font in font_manager.fontManager.ttflist}
    family = next((name for name in ("Arial", "Helvetica", "DejaVu Sans")
                   if name in families), "DejaVu Sans")
    plt.rcParams.update({"font.family": family, "font.sans-serif": [family],
                         "pdf.fonttype": 42, "ps.fonttype": 42,
                         "svg.fonttype": "none", "savefig.dpi": 300,
                         "savefig.bbox": None, "legend.frameon": False})
    return family


def save_figure(figure, output, stem, title, min_font=None):
    """Save one fixed-size PDF, editable SVG and 300-dpi PNG."""
    from matplotlib.text import Text
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    targets = [output / (stem + suffix) for suffix in (".pdf", ".svg", ".png")]
    if any(path.exists() for path in targets):
        raise FileExistsError("Choose an output directory without existing figure files.")
    family = plt.rcParams["font.sans-serif"][0]
    for item in figure.findobj(match=Text):
        item.set_fontfamily(family)
        if min_font is not None:
            item.set_fontsize(max(min_font, item.get_fontsize()))
    figure.canvas.draw()
    figure.savefig(targets[0], dpi=100, metadata={"Title": title, "Author": "Anonymized",
                                                "Creator": "Matplotlib"})
    figure.savefig(targets[1], dpi=100, metadata={"Title": title, "Creator": "Matplotlib",
                                                "Date": None})
    figure.savefig(targets[2], dpi=300)
    plt.close(figure)
    return targets

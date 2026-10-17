"""Render cancer-role and supplementary figures from their numerical source data."""
from pathlib import Path
import argparse
from importlib import import_module
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIGURES = {
    "2": ("Figure_2.json", "_figure_2", "Figure_2"),
    "S1": ("Supplementary_Figure_1.json", "_figure_s1", "Supplementary_Figure_1"),
    "S2": ("Supplementary_Figure_2.json", "_figure_s2", "Supplementary_Figure_2"),
    "S3": ("Supplementary_Figure_3.json", "_figure_s3", "Supplementary_Figure_3"),
    "S4": ("Supplementary_Figure_4.json", "_figure_s4", "Supplementary_Figure_4"),
}


def run(data: Path, output: Path, figures=None):
    """Read the supplied figure records and render PDF, SVG and PNG files.

    ``data`` is the extracted Supplementary_Data root or its Figure_Source_Data
    directory. ``figures`` optionally selects labels from 2, S1, S2, S3 and S4.
    Existing figure files are never overwritten. No input data are changed.
    """
    data, output = Path(data), Path(output)
    sources = data / "Figure_Source_Data" if (data / "Figure_Source_Data").is_dir() else data
    labels = list(FIGURES) if figures is None else list(figures)
    if len(labels) != len(set(labels)) or not labels or any(label not in FIGURES for label in labels):
        raise ValueError("Select unique figure labels from 2, S1, S2, S3 and S4.")
    records = {}
    for label in labels:
        filename, module, stem = FIGURES[label]
        records[label] = json.loads((sources / filename).read_text(encoding="utf-8-sig"))
        for suffix in (".pdf", ".svg", ".png"):
            if (output / (stem + suffix)).exists():
                raise FileExistsError("Choose an output directory without existing figure files.")
    output.mkdir(parents=True, exist_ok=True)
    paths = []
    for label in labels:
        module = import_module("." + FIGURES[label][1], __package__)
        with plt.rc_context():
            paths.extend(module.draw(records[label], output))
    return paths


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--figures", nargs="+", choices=tuple(FIGURES))
    args = parser.parse_args()
    run(args.data, args.output, args.figures)


if __name__ == "__main__":
    main()

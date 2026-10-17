




import argparse
import csv
import hashlib
import json
from collections import Counter
from fractions import Fraction as F
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.text import Text

PARSER = argparse.ArgumentParser(description="Render the fixed-record donor-composition figure.")
PARSER.add_argument("--source", type=Path, required=True)
PARSER.add_argument("--output", type=Path, required=True)
PARSER.add_argument("--work", type=Path, required=True)
PARSER.add_argument("--font-dir", type=Path)
PARSER.add_argument("--font-family", default="Arial")
ARGS = PARSER.parse_args()
WORK = ARGS.work.resolve()
DATA = ARGS.source.resolve()
OUTPUT = ARGS.output.resolve()
assert WORK.is_dir() and DATA.is_dir()
assert not OUTPUT.exists()
assert not any((WORK / name).exists() for name in ("Fig_Composition_source.json", "Fig_Composition_layout_check.json", "Fig_Composition.png"))
KEYS = [("PRODIGY", "BOTH_CORE"), ("DawnRank", "BOTH_CORE"),
        ("PRODIGY", "ONCOGENE"), ("DawnRank", "ONCOGENE")]
LABELS = ["Shared core\nPRODIGY", "Shared core\nDawnRank",
          "Oncogene only\nPRODIGY", "Oncogene only\nDawnRank"]
YS = [0, 1, 2.45, 3.45]
METHOD_COLOURS = {"PRODIGY": "#326B8C", "DawnRank": "#AB6554"}
SUPPORT_COLOURS = {"both": "#788B98", "entrant": "#528777",
                   "exit": "#B78054", "neither": "#E1E5E7"}
COMPONENTS = ["within_common", "reweight_common", "entrant", "exit"]
COMPONENT_COLOURS = ["#586B79", "#9A8DA9", "#528777", "#B78054"]
CHANGE_COLOURS = {"increase": "#528777", "tie": "#D7DDE1", "decrease": "#AB6554"}


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def prepare():
    summary_path = DATA / "ratio_decomposition_summary.csv"
    donors_path = DATA / "ratio_decomposition_donors.csv"
    summaries = read_csv(summary_path)
    donors = read_csv(donors_path)
    output = []
    for key in KEYS:
        row = next(r for r in summaries if (r["policy"], r["category"]) == key
                   and r["measure"] == "binary_fraction")
        records = [r for r in donors if (r["policy"], r["category"]) == key
                   and r["measure"] == "binary_fraction"]
        assert len(records) == 83 and len({r["donor_id"] for r in records}) == 83
        support = Counter(r["support"] for r in records)
        assert support == Counter({"both": int(row["both"]), "entrant": int(row["entrants"]),
                                   "exit": int(row["exits"]), "neither": int(row["neither"])})
        total0 = sum((F(r["R0"]) for r in records), F(0))
        total1 = sum((F(r["R1"]) for r in records), F(0))
        p0 = sum((F(r["numerator0"]) for r in records), F(0)) / total0
        p1 = sum((F(r["numerator1"]) for r in records), F(0)) / total1
        assert p0 == F(row["ratio0"]) and p1 == F(row["ratio1"])
        c = (p0 + p1) / 2
        components = {k: F(0) for k in COMPONENTS}
        changes = Counter({"increase": 0, "tie": 0, "decrease": 0})
        for r in records:
            r0, r1 = F(r["R0"]), F(r["R1"])
            a0, a1 = r0 / total0, r1 / total1
            assert bool(r0) == bool(r["donor_ratio0"])
            assert bool(r1) == bool(r["donor_ratio1"])
            if r0 and r1:
                f0, f1 = F(r["donor_ratio0"]), F(r["donor_ratio1"])
                assert f0 == F(r["numerator0"]) / r0 and f1 == F(r["numerator1"]) / r1
                components["within_common"] += (a0 + a1) / 2 * (f1 - f0)
                components["reweight_common"] += ((f0 + f1) / 2 - c) * (a1 - a0)
                changes["increase" if f1 > f0 else "decrease" if f1 < f0 else "tie"] += 1
            elif r1:
                components["entrant"] += a1 * (F(r["donor_ratio1"]) - c)
            elif r0:
                components["exit"] -= a0 * (F(r["donor_ratio0"]) - c)
        assert all(components[k] == F(row[k]) for k in COMPONENTS)
        assert sum(components.values(), F(0)) == p1 - p0 == F(row["difference"])
        assert sum(changes.values()) == support["both"]
        if key == ("PRODIGY", "BOTH_CORE"):
            common = [r for r in records if r["support"] == "both"]
            assert len(common) == 14
            assert all(F(r["donor_ratio0"]) == F(r["donor_ratio1"]) == 1 for r in common)
        output.append({"policy": key[0], "category": key[1], "support": dict(support),
                       "common_donor_change_counts": dict(changes),
                       "ratio0": str(p0), "ratio1": str(p1), "difference": str(p1 - p0),
                       "components": {k: str(v) for k, v in components.items()}})
    result = {"status": "PASS", "models": 85, "donors": 83,
              "checks": ["All 332 focal donor rows retained", "Exact pooled ratios reconstructed",
                         "Exact midpoint decomposition independently reconstructed",
                         "All common-support changes counted from exact ratios",
                         "All 14 PRODIGY shared-core common-support pairs equal one"],
              "sources": [{"path": p.name,
                           "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
                          for p in (summary_path, donors_path)], "rows": output}
    return result


def panel(ax, letter, title):
    ax.text(0, 1.08, letter, transform=ax.transAxes, fontweight="bold", fontsize=11,
            ha="left", va="bottom")
    ax.annotate(title, xy=(0, 1.08), xycoords="axes fraction", xytext=(15, 0),
                textcoords="offset points", fontsize=10.5, ha="left", va="bottom")


def four_rows(ax):
    ax.set_yticks(YS, LABELS)
    ax.set_ylim(4.05, -.65)
    ax.tick_params(axis="y", length=0, pad=4)
    for label, key in zip(ax.get_yticklabels(), KEYS):
        label.set_color(METHOD_COLOURS[key[0]])
    ax.spines["left"].set_visible(False)


def base_style(ax):
    ax.set_axisbelow(True)
    ax.grid(axis="x", color="#E5E9EC", linewidth=.5)
    ax.tick_params(length=3, width=.6, pad=3)


def render(data):
    if ARGS.font_dir:
        for font in sorted(ARGS.font_dir.glob("*.ttf")):
            font_manager.fontManager.addfont(str(font))
    font_manager.findfont(ARGS.font_family, fallback_to_default=False)
    plt.rcParams.update({"font.family": ARGS.font_family, "font.size": 9.4, "axes.labelsize": 9.4,
                         "xtick.labelsize": 9.0, "ytick.labelsize": 9.2, "legend.fontsize": 8.8,
                         "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": .65,
                         "axes.unicode_minus": False, "pdf.fonttype": 42, "svg.fonttype": "none",
                         "legend.frameon": False, "savefig.facecolor": "white"})
    fig = plt.figure(figsize=(220 / 25.4, 200 / 25.4), facecolor="white")
    a = fig.add_axes((.125, .655, .335, .260))
    b = fig.add_axes((.665, .655, .305, .260))
    c = fig.add_axes((.125, .200, .335, .290))
    d = fig.add_axes((.665, .200, .305, .290))
    rows = data["rows"]

                                                                                 
    panel(a, "a", "PRODIGY shared core")
    p = [100 * float(F(rows[0][k])) for k in ("ratio0", "ratio1")]
    a.plot([0, 1], p, color=METHOD_COLOURS["PRODIGY"], lw=1.8)
    a.scatter([0, 1], p, s=32, facecolors=["white", METHOD_COLOURS["PRODIGY"]],
              edgecolors=METHOD_COLOURS["PRODIGY"], lw=1.3, zorder=4)
    a.plot([0, 1], [100, 100], color="#555E64", lw=1.3)
    a.scatter([0, 1], [100, 100], s=30, facecolors=["white", "#555E64"],
              edgecolors="#555E64", lw=1.2, zorder=4)
    for x, value in enumerate(p):
        a.text(x, value - .65, f"{value:.2f}", ha="center", va="top",
               color=METHOD_COLOURS["PRODIGY"], fontsize=9)
    a.text(.5, 100.75, "14/14 common-support donors", ha="center", va="bottom", fontsize=9)
    a.text(.5, 99.25, "Each 100% in both states", ha="center", va="top", fontsize=8.8, color="#555E64")
    a.text(.50, 87.8, "Pooled fraction", color=METHOD_COLOURS["PRODIGY"], ha="center", fontsize=9)
    a.text(.50, 86.4, f"+{100 * float(F(rows[0]['difference'])):.3f} percentage points",
           color=METHOD_COLOURS["PRODIGY"], ha="center", fontsize=8.8)
    a.set_xlim(-.18, 1.18)
    a.set_ylim(85, 103.2)
    a.set_xticks([0, 1], ["Old", "New"])
    a.set_yticks([85, 90, 95, 100])
    a.set_ylabel("Positive fraction (%)")
    a.grid(axis="y", color="#E5E9EC", lw=.5)
    a.set_axisbelow(True)

                                                                               
    panel(b, "b", "Category support")
    base_style(b)
    for y, row in zip(YS, rows):
        left = 0
        for key in ("both", "entrant", "exit", "neither"):
            count = row["support"].get(key, 0)
            b.barh(y, count, left=left, height=.52, color=SUPPORT_COLOURS[key], linewidth=0)
            if count:
                b.text(left + count / 2, y, str(count), ha="center", va="center", fontsize=8.8,
                       color="white" if key in ("both", "entrant") else "#222222")
            left += count
        assert left == 83
    four_rows(b)
    b.set_xlim(0, 83)
    b.set_xticks([0, 20, 40, 60, 83])
    b.set_xlabel("Donors (all 83 retained)")
    b.legend(handles=[Patch(facecolor=SUPPORT_COLOURS[k], label=lab) for k, lab in
                      [("both", "Both states"), ("entrant", "New only"),
                       ("exit", "Old only"), ("neither", "Neither")]],
             loc="upper center", bbox_to_anchor=(.5, -.29), ncol=2, borderaxespad=0,
             handlelength=1.1, handletextpad=.45, columnspacing=1.15, labelspacing=.55)

                                                                  
    panel(c, "c", "Common-support donor changes")
    base_style(c)
    for y, row in zip(YS, rows):
        left = 0
        for key in ("decrease", "tie", "increase"):
            count = row["common_donor_change_counts"][key]
            c.barh(y, count, left=left, height=.52, color=CHANGE_COLOURS[key], linewidth=0)
            if count:
                c.text(left + count / 2, y, str(count), ha="center", va="center", fontsize=8.8,
                       color="white" if key in ("decrease", "increase") else "#222222")
            left += count
        assert left == row["support"]["both"]
    four_rows(c)
    c.set_xlim(0, 60)
    c.set_xticks([0, 20, 40, 60])
    c.set_xlabel("Common-support donors")
    c.legend(handles=[Patch(facecolor=CHANGE_COLOURS[k], label=lab) for k, lab in
                      [("decrease", "Decrease"), ("tie", "Unchanged"), ("increase", "Increase")]],
             loc="upper center", bbox_to_anchor=(.5, -.30), ncol=3, borderaxespad=0,
             handlelength=1.0, handletextpad=.4, columnspacing=.8)

                                                                                
    panel(d, "d", "Sources of pooled change")
    base_style(d)
    for y, row in zip(YS, rows):
        pos, neg = 0., 0.
        for key, colour in zip(COMPONENTS, COMPONENT_COLOURS):
            value = 100 * float(F(row["components"][key]))
            start = pos if value >= 0 else neg
            d.barh(y, value, left=start, height=.52, color=colour, linewidth=0)
            if value >= 0:
                pos += value
            else:
                neg += value
        net = 100 * float(F(row["difference"]))
        d.scatter(net, y, marker="D", s=27, facecolor="white", edgecolor="#20272C", lw=1.0, zorder=5)
        d.text(13.7, y, f"{net:+.2f}", ha="right", va="center", fontsize=8.8, color="#222222")
    four_rows(d)
    d.set_xlim(-14.0, 14.0)
    d.set_xticks([-10, -5, 0, 5, 10])
    d.axvline(0, color="#86919A", lw=.7, ls=(0, (3, 3)), zorder=0)
    d.text(13.7, -.72, "Net", ha="right", va="bottom", fontsize=8.8)
    d.set_xlabel("Contribution (percentage points)")
    d.legend(handles=[Patch(facecolor=col, label=lab) for col, lab in
                      zip(COMPONENT_COLOURS, ["Within donor", "Reweighting", "New only", "Old only"])] +
                     [Line2D([], [], marker="D", color="#20272C", mfc="white", linestyle="", label="Net change", markersize=4.5)],
             loc="upper center", bbox_to_anchor=(.5, -.29), ncol=2, borderaxespad=0,
             handlelength=1.1, handletextpad=.45, columnspacing=1.0, labelspacing=.55)

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    bounds = fig.bbox
    outside = []
    visible_text = 0
    for item in fig.findobj(match=Text):
        if not item.get_visible() or not item.get_text():
            continue
        visible_text += 1
        box = item.get_window_extent(renderer)
        if box.x0 < bounds.x0 - 1 or box.y0 < bounds.y0 - 1 or box.x1 > bounds.x1 + 1 or box.y1 > bounds.y1 + 1:
            outside.append(item.get_text())
    assert not outside, outside
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, metadata={"Title": "Donor composition and pooled category fractions",
                                 "Author": "Anonymized", "Creator": "Matplotlib",
                                 "Subject": "Fixed-record descriptive donor-composition accounting"})
    fig.savefig(WORK / "Fig_Composition.png", dpi=180)
    (WORK / "Fig_Composition_layout_check.json").write_text(json.dumps({
        "status": "PASS", "size_mm": [220, 200], "visible_text_items_checked": visible_text,
        "outside_canvas": outside, "figure_sha256": hashlib.sha256(OUTPUT.read_bytes()).hexdigest(),
        "runtime": {"matplotlib": matplotlib.__version__, "font": ARGS.font_family + "; embedded TrueType"},
        "caption": "Main_Donor_Composition_Figure.tex",
        "scope": "Visualisation only; exact analysis records retained"
    }, indent=2) + "\n", encoding="utf-8")
    plt.close(fig)


if __name__ == "__main__":
    data = prepare()
    (WORK / "Fig_Composition_source.json").write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    render(data)
    print(json.dumps({"status": "PASS", "figure": str(OUTPUT), "source_checks": data["checks"]}, indent=2))

"""Cancer-role allocations, functional fractions and oncogene identity flows."""
from fractions import Fraction
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from ._figure_common import apply_figure_style, save_figure

COLOURS = {"PRODIGY": "#326B8C", "DawnRank": "#AB6554"}
PALE = {"PRODIGY": "#C9D8E1", "DawnRank": "#E7D2CB"}
CATEGORIES = ("BOTH_CORE", "ONCOGENE")
KEYS = [(cat, method) for cat in CATEGORIES for method in COLOURS]
LABELS = ["Shared core\nPRODIGY", "Shared core\nDawnRank", "Oncogene only\nPRODIGY", "Oncogene only\nDawnRank"]
POSITIONS = (0, 1, 2.4, 3.4)

def number(value):
    return float(Fraction(str(value)))

def panel(ax, letter, title, colour="#222222"):
    ax.text(0, 1.10, letter, transform=ax.transAxes, fontweight="bold", va="bottom", fontsize=11)
    ax.annotate(title, (0, 1.10), xycoords="axes fraction", xytext=(15, 0), textcoords="offset points",
                va="bottom", fontsize=10.5, color=colour)


def style(ax):
    ax.set_axisbelow(True)
    ax.grid(axis="x", color="#E7EAEC", linewidth=.55)
    ax.tick_params(length=3, width=.6, pad=3)


def four_rows(ax):
    ax.set_yticks(POSITIONS, LABELS)
    ax.set_ylim(4.05, -.65)
    for label, (_, method) in zip(ax.get_yticklabels(), KEYS):
        label.set_color(COLOURS[method])


def draw(data, output):
    family = apply_figure_style()
    plt.rcParams.update({"font.family": family, "font.size": 9.4, "axes.labelsize": 9.4,
                         "xtick.labelsize": 9, "ytick.labelsize": 9.1, "legend.fontsize": 9,
                         "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": .65,
                         "axes.unicode_minus": False, "pdf.fonttype": 42, "svg.fonttype": "none",
                         "legend.frameon": False, "savefig.facecolor": "white"})
    fig = plt.figure(figsize=(220 / 25.4, 236 / 25.4), facecolor="white")
    boxes = [(.155, .701, .315, .221), (.650, .701, .315, .221),
             (.155, .370, .315, .205), (.650, .370, .315, .205),
             (.155, .060, .315, .210), (.650, .060, .315, .210)]
    axes = [fig.add_axes(box) for box in boxes]
    for ax in axes:
        style(ax)
    levels = {(x["category"], x["method"], x["state"]): x for x in data["levels"]}
    a, b, c, d, e, f = axes
    for y, (category, method) in zip(POSITIONS, KEYS):
        colour = COLOURS[method]
        for state, offset in (("0", -.17), ("1", .17)):
            row = levels[(category, method, state)]
            positive, zero, total = map(number, (row["H"], row["Z"], row["R"]))
            a.barh(y + offset, positive, height=.25, color=colour, linewidth=0)
            a.barh(y + offset, zero, left=positive, height=.25, color=PALE[method], linewidth=0)
            a.text(total + .035, y + offset, f"{total:.2f}", va="center", fontsize=8.8)
        old, new = [100 * number(levels[(category, method, state)]["pooled_positive_fraction"]) for state in ("0", "1")]
        b.plot([old, new], [y, y], color=colour, linewidth=1.5)
        b.scatter(old, y, s=33, facecolors="white", edgecolors=colour, linewidths=1.2, zorder=3)
        b.scatter(new, y, s=25, facecolors=colour, edgecolors=colour, zorder=4)
        b.text(old, y - .24, f"{old:.1f}", color=colour, ha="center", va="bottom", fontsize=8.8)
        b.text(new, y + .22, f"{new:.1f}", color=colour, ha="center", va="top", fontsize=8.8)
    for ax in (a, b, c, d):
        four_rows(ax)
    a.set_xlim(0, 1.87)
    a.set_xticks([0, .5, 1, 1.5])
    a.set_xlabel("Allocated slots per 10 offered")
    a.text(0, 1.01, "Upper bar: old    Lower bar: new", transform=a.transAxes, fontsize=8.8)
    a.legend(handles=[Patch(facecolor="#4F626C", label="Positive"), Patch(facecolor="#CCD4D8", label="Measured zero")],
             loc="upper left", bbox_to_anchor=(-.01, -.22), ncol=2, borderaxespad=0, handlelength=1.1,
             columnspacing=1.0, handletextpad=.45)
    b.set_xlim(38, 102)
    b.set_xticks([40, 60, 80, 100])
    b.set_xlabel("Pooled positive fraction (%)")
    b.legend(handles=[Line2D([], [], color="#596167", marker="o", markerfacecolor="white", linestyle="", label="Old"),
                      Line2D([], [], color="#596167", marker="o", linestyle="", label="New")],
             loc="upper left", bbox_to_anchor=(.16, -.22), ncol=2, borderaxespad=0, handlelength=1.1, columnspacing=1.2)
    panel(a, "a", "Category allocation")
    panel(b, "b", "Within-category positive fraction")
    for ax, endpoint, scale in ((c, "delta_R", 1), (d, "delta_p", 100)):
        for y, (category, method) in zip(POSITIONS, KEYS):
            row = next(x for x in data["envelopes"] if x["category"] == category and x["method"] == method and x["endpoint"] == endpoint)
            low, high, estimate = scale * row["lower"], scale * row["upper"], scale * number(row["value"])
            ax.plot([low, high], [y, y], color=COLOURS[method], linewidth=1.8)
            ax.plot([low, high], [y, y], "|", color=COLOURS[method], markersize=6, markeredgewidth=1.1)
            ax.scatter(estimate, y, s=26, color=COLOURS[method], zorder=3)
        ax.axvline(0, color="#69747B", linestyle="--", linewidth=.8)
        ax.text(0, 1.01, "Donor-composition sensitivity envelopes", transform=ax.transAxes, fontsize=8.8)
    c.set_xlim(-1.4, .8)
    c.set_xticks([-1, -.5, 0, .5])
    c.set_xlabel("Allocation change per 10 offered slots")
    d.set_xlim(-27, 27)
    d.set_xticks([-20, -10, 0, 10, 20])
    d.set_xlabel("Positive-fraction change (percentage points)")
    panel(c, "c", "Allocation change")
    panel(d, "d", "Fraction change")
    flow_rows = {(x["method"], x["label_group"]): x for x in data["oncogene_identity_flows"]}
    stages = ("old", "retained", "loss", "gain", "new")
    stage_labels = ("Old", "Retained", "Removed", "Added", "New")
    for ax, method, letter in ((e, "PRODIGY", "e"), (f, "DawnRank", "f")):
        positive, zero = flow_rows[(method, "positive")], flow_rows[(method, "zero")]
        for y, stage in enumerate(stages):
            pos, zer = number(positive[stage]), number(zero[stage])
            ax.barh(y, pos, height=.56, color=COLOURS[method], linewidth=0)
            ax.barh(y, zer, left=pos, height=.56, color=PALE[method], linewidth=0)
            ax.text(pos + zer + .035, y, f"{pos + zer:.2f}", va="center", fontsize=8.8)
        ax.set_yticks(range(5), stage_labels)
        ax.set_ylim(4.6, -.6)
        ax.set_xlim(0, 1.87)
        ax.set_xticks([0, .5, 1, 1.5])
        ax.set_xlabel("Oncogene-only allocation per 10 slots")
        panel(ax, letter, method + " identity flows", COLOURS[method])
    for ax in axes:
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)
    fig.canvas.draw()
    return save_figure(fig, output, "Figure_2", "Cancer-role allocation and identity flows")

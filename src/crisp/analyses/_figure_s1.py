"""Numerical stopping, strict cutoffs and functional allocation bounds."""
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter
import numpy as np
from ._figure_common import apply_figure_style, save_figure

COLORS = {'first': '#225F85', 'second': '#AE6243', 'negative': '#326A8B', 'positive': '#A75E43', 'zero': '#B7BEC1', 'unknown': '#E4CFA1', 'menu': '#CCA858'}

def panel(ax, letter, title):
    ax.text(-0.17, 1.055, '(' + letter + ')', transform=ax.transAxes, fontsize=11.5, fontweight='bold', ha='left', va='bottom')
    ax.set_title(title, loc='left', pad=13)
    ax.tick_params(direction='out', length=3)

def note_axes(ax, direction='y'):
    ax.grid(axis=direction, color='#DFE3E5', linewidth=0.55, zorder=0)
    ax.set_axisbelow(True)

def draw(source, output):
    family = apply_figure_style()
    matplotlib.rcParams.update({'font.family': family, 'font.size': 10, 'axes.labelsize': 10.5, 'axes.titlesize': 11.5, 'xtick.labelsize': 10, 'ytick.labelsize': 10, 'legend.fontsize': 10, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.edgecolor': '#555555', 'axes.labelcolor': '#222222', 'text.color': '#222222', 'xtick.color': '#333333', 'ytick.color': '#333333', 'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none', 'axes.linewidth': 0.7, 'xtick.major.width': 0.7, 'ytick.major.width': 0.7, 'savefig.facecolor': 'white', 'mathtext.default': 'regular', 'mathtext.fontset': 'custom', 'mathtext.rm': family, 'mathtext.it': family + ':italic', 'mathtext.bf': family + ':bold'})
    figure, axes = plt.subplots(3, 2, figsize = (220 / 25.4, 292 / 25.4))
    figure.subplots_adjust(left = 0.125, right = 0.975, bottom = 0.12, top = 0.95, wspace = 0.46, hspace = 0.65)

    ax = axes[0, 0]
    panel(ax, 'a', 'Iterations to stopping')
    for key, offset in [('first', -0.12), ('second', 0.12)]:
        selected = [row for row in source['iterations'] if row['transition'] == key]
        x = np.array([-np.log10(float(row['epsilon'])) for row in selected]) + offset
        midpoint = np.array([(row['minimum'] + row['maximum']) / 2 for row in selected])
        errors = np.array([(row['maximum'] - row['minimum']) / 2 for row in selected])
        ax.errorbar(x, midpoint, yerr = errors, fmt = 'o' if key == 'first' else 's', ms = 5, capsize = 4, color = COLORS[key], label = 'First transition' if key == 'first' else 'Second transition', lw = 1.4, zorder = 3)
    ax.set_xticks([4, 6, 8, 12], ['1e-4', '1e-6', '1e-8', '1e-12'])
    ax.set_xlim(3.4, 12.6)
    ax.set_ylim(0, 530)
    ax.set_xlabel('Stopping tolerance')
    ax.set_ylabel('Iteration range across 340 cells')
    note_axes(ax)
    ax.legend(frameon = False, loc = 'upper left', handlelength = 1, borderaxespad = 0.2)

    ax = axes[0, 1]
    panel(ax, 'b', 'Changed ten-slot allocations')
    cells = ['G0_H0', 'G0_H1', 'G1_H0', 'G1_H1']
    for key, offset in [('first', -0.18), ('second', 0.18)]:
        lookup = {row['cell']: row['count'] for row in source['changed'] if row['transition'] == key}
        values = [lookup[cell] for cell in cells]
        x = np.arange(4) + offset
        ax.bar(x, values, width = 0.33, color = COLORS[key], zorder = 3)
        for position, value in zip(x, values):
            ax.text(position, value + 0.7, str(value), ha = 'center', fontsize = 10)
    ax.set_xticks(range(4), ['G0H0', 'G0H1', 'G1H0', 'G1H1'])
    ax.set_ylim(0, 24)
    ax.set_yticks([0, 5, 10, 15, 20])
    ax.set_xlabel('Graph and evidence cell')
    ax.set_ylabel('Changed models out of 85')
    note_axes(ax)

    ax = axes[1, 0]
    panel(ax, 'c', 'Separation at strict cutoffs')
    for key, marker in [('first', 'o'), ('second', 's')]:
        selected = [row for row in source['ratios'] if row['transition'] == key]
        values = [row['minimum_gap_over_twice_error'] for row in selected]
        ax.plot(range(4), values, marker = marker, color = COLORS[key], lw = 1.2, ms = 5, zorder = 3)
    ax.axhline(1, color = '#6B6B6B', ls = '--', lw = 0.9)
    ax.set_yscale('log')
    ax.set_ylim(0.65, 1.7 * max(row['minimum_gap_over_twice_error'] for row in source['ratios']))
    ax.set_xticks(range(4), ['1', '5', '10', '20'])
    ax.set_xlabel('Offered candidate budget')
    ax.set_ylabel('Minimum gap / (2 × error bound)')
    ax.yaxis.set_major_locator(LogLocator(base = 10, numticks = 5))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, position: '1e' + str(int(round(np.log10(value))))))
    ax.yaxis.set_minor_formatter(NullFormatter())
    note_axes(ax)

    ax = axes[1, 1]
    panel(ax, 'd', 'Joint functional bounds')
    labels = []
    for index, row in enumerate(source['intervals']):
        position = 5 - index
        color = COLORS[row['transition']]
        ax.plot([row['lower'], row['upper']], [position, position], color = color, lw = 2.2, solid_capstyle = 'butt')
        ax.plot([row['lower'], row['upper']], [position, position], '|', color = color, ms = 8, markeredgewidth = 1.1)
        ax.scatter([row['known']], [position], color = color, s = 22, zorder = 3)
        prefix = 'First' if row['transition'] == 'first' else 'Second'
        exponent = int(-np.log10(float(row['epsilon'])))
        labels.append(prefix + ' 1e-' + str(exponent))
    ax.axvline(0, color = '#777777', lw = 0.8, ls = '--')
    ax.set_yticks(range(5, -1, -1), labels)
    ax.set_ylim(-0.6, 5.6)
    ax.set_xlim(-1.2, 0.32)
    ax.set_xticks([-1, -0.5, 0])
    ax.set_xlabel('Positive allocation change / 10 slots')
    note_axes(ax, 'x')

    ax = axes[2, 0]
    panel(ax, 'e', 'Donor direction classifications')
    labels = []
    left = np.zeros(4)
    ordered = source['directions']
    for category, label in [('negative', 'Negative'), ('positive', 'Positive'), ('zero', 'Zero'), ('unknown', 'Unresolved')]:
        values = np.array([row[category] for row in ordered])
        ax.barh(range(4), values, left = left, height = 0.57, color = COLORS[category], label = label, zorder = 3)
        for index, value in enumerate(values):
            ax.text(left[index] + value / 2, index, str(value), color = 'white' if category in ['negative', 'positive'] else '#222222', ha = 'center', va = 'center', fontsize = 10)
        left += values
    for row in ordered:
        labels.append(('First' if row['transition'] == 'first' else 'Second') + (' original' if row['epsilon'] == '1e-04' else ' strict'))
    ax.set_yticks(range(4), labels)
    ax.invert_yaxis()
    ax.set_xlim(0, 83)
    ax.set_xticks([0, 20, 40, 60, 83])
    ax.set_xlabel('Donors')
    ax.legend(frameon = False, loc = 'upper left', bbox_to_anchor = (-0.02, -0.23), ncol = 2, columnspacing = 1.3, handlelength = 1.1)
    note_axes(ax, 'x')

    ax = axes[2, 1]
    panel(ax, 'f', 'Strict membership accounting')
    labels = []
    for index, row in enumerate(source['certification']):
        position = 7 - index
        color = COLORS[row['transition']]
        ax.barh(position, row['score_gap_certified'], height = 0.6, color = color, zorder = 3)
        if (row['complete_menu']):
            ax.barh(position, row['complete_menu'], left = row['score_gap_certified'], height = 0.6, color = COLORS['menu'], zorder = 3)
            ax.text(343, position, '(+4)', fontsize = 10, va = 'center')
        ax.text(168, position, str(row['score_gap_certified']), color = 'white', fontsize = 10, ha = 'center', va = 'center')
        labels.append(('First' if row['transition'] == 'first' else 'Second') + ' k=' + str(row['budget']))
    ax.set_yticks(range(7, -1, -1), labels)
    ax.set_xlim(0, 372)
    ax.set_xticks([0, 100, 200, 340])
    ax.set_xlabel('Model–cells with resolved membership')
    ax.set_ylim(-0.7, 7.7)
    note_axes(ax, 'x')
    return save_figure(figure, output, "Supplementary_Figure_1", "Numerical stopping and strict candidate cutoffs")

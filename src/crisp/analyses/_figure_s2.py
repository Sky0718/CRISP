"""Crossed input states and shared-label bounds."""
from ._figure_common import apply_figure_style, save_figure
from fractions import Fraction
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.text import Text
import numpy as np

COLOURS = {'PRODIGY':'#326B8C','DawnRank':'#AB6554','positive':'#528777','zero':'#AAB4BE','unknown':'#D9CBA7','vacant':'#F0F1F2'}
CELLS = ['G0H0','G1H0','G0H1','G1H1']
ORDER = ['graph_at_old_evidence','graph_at_new_evidence','evidence_at_old_graph','evidence_at_new_graph','interaction','joint']
LABELS = ['Graph, old evidence','Graph, new evidence','Evidence, old graph','Evidence, new graph','Graph × evidence','Joint update']

def value(number):
    return float(Fraction(str(number)))

def panel(axis, letter, title):
    axis.text(0,1.065,'(' + letter + ')',transform=axis.transAxes,fontweight='bold',va='bottom',fontsize=11)
    axis.annotate(title,(0,1.065),xycoords='axes fraction',xytext=(23,0),textcoords='offset points',va='bottom',fontsize=11)

def draw(data, output):
    family = apply_figure_style()
    plt.rcParams.update({'font.family':family,'font.size':10.5,'axes.labelsize':10.5,'xtick.labelsize':10,'ytick.labelsize':10,'legend.fontsize':10,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.65,'pdf.fonttype':42,'svg.fonttype':'none','axes.unicode_minus':False,'legend.frameon':False,'figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white'})
    figure = plt.figure(figsize=(220/25.4,250/25.4))
    diagram = figure.add_axes([.065,.765,.88,.18])
    diagram.set_axis_off()
    diagram.set_xlim(0,1)
    diagram.set_ylim(0,1)
    panel(diagram,'a','Colorectal cancer models: crossed input states')
    diagram.text(.5,-.17,'85 models, 83 donors; 10 candidate slots per model',ha='center',va='top',fontsize=11)
    for x in (.23,.615):
        diagram.plot([x,x],[.04,.96],color='black',linewidth=.8)
    for y in (.35,.65):
        diagram.plot([0,1],[y,y],color='black',linewidth=.8)
    for x,label in [(.4225,'Old evidence'),(.8075,'New evidence')]:
        diagram.text(x,.805,label,ha='center',va='center',fontsize=12)
    for y,label in [(.50,'Old graph'),(.195,'New graph')]:
        diagram.text(.115,y,label,ha='center',va='center',fontsize=12)
    for x,y,label in [(.4225,.50,'G0H0'),(.8075,.50,'G0H1'),(.4225,.195,'G1H0'),(.8075,.195,'G1H1')]:
        diagram.text(x,y,label,ha='center',va='center',fontsize=12)
    for index,method in enumerate(('PRODIGY','DawnRank')):
        axis = figure.add_axes([.29+index*.385,.43,.27,.225])
        rows = {row['contrast']:row for row in data['contrasts'] if row['method'] == method}
        for position,name in enumerate(ORDER):
            row = rows[name]
            lower,upper,known = [value(row[key]) for key in ('lower','upper','known')]
            axis.plot([lower,upper],[position,position],color=COLOURS[method],linewidth=2)
            axis.scatter([known],[position],s=20,color=COLOURS[method],zorder=3)
        axis.set_yticks(range(6),LABELS if index == 0 else [])
        axis.set_ylim(5.6,-.6)
        axis.set_xlim(-1.2,.20)
        axis.set_xticks([-1.2,-.6,0])
        axis.set_xlabel('Depletion-positive allocation change\n(per 10 candidate slots)')
        axis.axvline(0,color='#89969F',linestyle=(0,(3,3)),linewidth=.7)
        axis.set_axisbelow(True)
        axis.grid(axis='x',color='#E2E7EA',linewidth=.4)
        panel(axis,chr(98+index),method)
    identity = figure.add_axes([.29,.10,.54,.215])
    labels = []
    for index,method in enumerate(('PRODIGY','DawnRank')):
        row = next(row for row in data['contrasts'] if (row['method'] == method and row['contrast'] == 'joint'))
        for offset,shared in enumerate((False,True)):
            position = index * 2.5 + offset
            lower,upper = [value(row[key]) for key in (('lower','upper') if shared else ('independent_lower','independent_upper'))]
            colour = COLOURS[method] if shared else '#87949D'
            identity.plot([lower,upper],[position,position],color=colour,linewidth=2.3 if shared else 1.7)
            identity.scatter([value(row['known'])],[position],s=24,facecolors=colour if shared else 'white',edgecolors=colour,zorder=3)
            identified = (upper < 0)
            identity.text(1.025,position,'Decrease' if identified else 'Unresolved',transform=identity.get_yaxis_transform(),va='center',color=colour,fontsize=10)
            labels.append(method + (', shared labels' if shared else ', separate labels'))
    identity.set_yticks([0,1,2.5,3.5],labels)
    identity.set_ylim(4.15,-.65)
    identity.set_xlim(-1.4,.2)
    identity.set_xticks([-1.2,-.6,0])
    identity.set_xlabel('Joint depletion-positive allocation change\n(per 10 candidate slots)')
    identity.axvline(0,color='#89969F',linestyle=(0,(3,3)),linewidth=.7)
    identity.grid(axis='x',color='#E2E7EA',linewidth=.4)
    identity.set_axisbelow(True)
    panel(identity,'d','Matching gene identities resolves update direction')
    return save_figure(figure, output, "Supplementary_Figure_2", "Crossed input states and shared-label bounds", min_font=10)

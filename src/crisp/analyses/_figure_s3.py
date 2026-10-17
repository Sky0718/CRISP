"""Donor directions, candidate overlap and seed variation."""
from ._figure_common import apply_figure_style, save_figure
from fractions import Fraction
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.text import Text
import numpy as np

COLOURS = {'PRODIGY':'#326B8C','DawnRank':'#AB6554','other':'#82929B','opposed':'#A45B46'}
SIGNS = ['negative','positive','zero','unresolved']
CELLS = ['G0H0','G1H0','G0H1','G1H1']

def number(value):
    return float(Fraction(str(value)))

def panel(axis, letter, title):
    axis.text(0,1.08,'(' + letter + ')',transform=axis.transAxes,fontweight='bold',va='bottom',fontsize=11)
    axis.annotate(title,(0,1.08),xycoords='axes fraction',xytext=(23,0),textcoords='offset points',va='bottom',fontsize=11)

def draw(data, output):
    family = apply_figure_style()
    plt.rcParams.update({'font.family':family,'font.size':10.5,'axes.labelsize':10.5,'xtick.labelsize':10,'ytick.labelsize':10,'legend.fontsize':10,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.65,'pdf.fonttype':42,'svg.fonttype':'none','axes.unicode_minus':False,'legend.frameon':False,'figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white'})
    figure = plt.figure(figsize=(220/25.4,285/25.4))
    a = figure.add_axes([.14,.755,.30,.175])
    matrix = np.array([[data['paired_directions']['matrix'].get(first+'_'+second,0) for second in SIGNS] for first in SIGNS])
    a.imshow(matrix,cmap=LinearSegmentedColormap.from_list('density',['#F2F5F5','#326B8C']),vmin=0,vmax=matrix.max(),aspect='auto')
    for row in range(4):
        for column in range(4):
            a.text(column,row,str(matrix[row,column]),ha='center',va='center',color='white' if matrix[row,column] > matrix.max()*.65 else '#24353D')
    for row,column in [(0,1),(1,0)]:
        a.add_patch(Rectangle((column-.49,row-.49),.98,.98,fill=False,edgecolor=COLOURS['opposed'],linewidth=1.6))
    a.set_xticks(range(4),['Negative','Positive','Zero','Unresolved'],rotation=25,ha='right')
    a.set_yticks(range(4),['Negative','Positive','Zero','Unresolved'])
    a.tick_params(length=0)
    a.set_xlabel('DawnRank')
    a.set_ylabel('PRODIGY')
    panel(a,'a','Depletion-capture direction by donor')
    b = figure.add_axes([.65,.755,.30,.175])
    histogram = data['seed_identity_persistence']['frequency_histogram']
    counts = [histogram[str(seed_count)] for seed_count in range(6)]
    colours = ['#AAB4BE']*5 + [COLOURS['opposed']]
    b.bar(range(6),counts,color=colours,width=.7)
    for position,count in enumerate(counts):
        b.text(position,count+1.8,str(count),ha='center',va='bottom',fontsize=10)
    b.set_xticks(range(6))
    b.set_yticks([0,20,40,60,80])
    b.set_xlim(-.65,5.65)
    b.set_ylim(0,83)
    b.set_xlabel('PRODIGY seeds with opposite direction')
    b.set_ylabel('Donors')
    b.set_axisbelow(True)
    b.grid(axis='y',color='#E2E7EA',linewidth=.4)
    panel(b,'b','Opposite directions across seeds')
    c = figure.add_axes([.14,.435,.30,.17])
    for position,direction,colour in [(0,'removal',COLOURS['DawnRank']),(1,'addition',COLOURS['PRODIGY'])]:
        rows = sorted((row for row in data['donor_overlap'] if row['direction'] == direction),key=lambda row:row['donor_id'])
        values = np.array([100*number(row['ratio']) for row in rows if row['ratio'] != ''])
        offsets = np.array([((index*37)%83)/82-.5 for index in range(len(values))])*.36
        c.scatter(np.full(len(values),position)+offsets,values,s=9,color=colour,alpha=.58,linewidths=0)
        c.boxplot([values],positions=[position],widths=.48,showfliers=False,patch_artist=False,medianprops={'color':'#24353D','linewidth':1.4},boxprops={'color':'#24353D','linewidth':.8},whiskerprops={'color':'#24353D','linewidth':.8},capprops={'color':'#24353D','linewidth':.8})
    c.set_xticks([0,1],['Removal','Addition'])
    c.set_xlim(-.6,1.6)
    c.set_ylim(-3,103)
    c.set_yticks([0,25,50,75,100])
    c.set_ylabel('Same-gene overlap within donor (%)')
    c.set_axisbelow(True)
    c.grid(axis='y',color='#E2E7EA',linewidth=.4)
    panel(c,'c','Genes removed or added by both methods')
    d = figure.add_axes([.73,.435,.22,.17])
    seed = {row['cell']:row for row in data['seed_variability']}
    values = [number(seed[cell]['donor_weighted_mean_pair_l1']) for cell in CELLS]+[number(data['turnover_cancellation']['mean_seed_l1'])]
    for position,value in enumerate(values):
        colour = COLOURS['PRODIGY'] if position == 4 else '#AAB8C0'
        d.barh(position,value,height=.6,color=colour)
        d.text(value+.17,position,f'{value:.2f}',va='center',fontsize=10)
    d.set_yticks(range(5),['Within '+cell for cell in CELLS]+['Paired update'])
    d.set_ylim(4.6,-.6)
    d.set_xlim(0,14.3)
    d.set_xticks([0,6,12])
    d.set_xlabel('Mean allocation distance (L1)')
    panel(d,'d','Update and seed variation')
    e = figure.add_axes([.14,.105,.30,.175])
    for method,offset,colour in [('PRODIGY',-.16,COLOURS['PRODIGY']),('DawnRank',.16,COLOURS['DawnRank'])]:
        rows = {row['contrast']:row for row in data['evidence_allocation'] if row['method'] == method}
        values = [rows[name]['allocation_changed_models'] for name in ('evidence_at_old_graph','evidence_at_new_graph')]
        bars = e.bar(np.arange(2)+offset,values,width=.30,color=colour,label=method)
        for bar,value in zip(bars,values):
            e.text(bar.get_x()+bar.get_width()/2,value+2,str(value),ha='center',va='bottom',fontsize=10)
    e.set_xticks([0,1],['Old graph','New graph'])
    e.set_ylabel('Models with changed allocation (of 85)')
    e.set_ylim(0,85)
    e.set_yticks([0,20,40,60,80])
    e.legend(loc='upper left',bbox_to_anchor=(-.07,-.26),ncol=2,columnspacing=.75,handlelength=1.0)
    e.set_axisbelow(True)
    e.grid(axis='y',color='#E2E7EA',linewidth=.4)
    panel(e,'e','Effect of evidence construction')
    f = figure.add_axes([.70,.105,.25,.175])
    rows = data['seed_joint']+[data['mean_PRO_joint']]
    for position,row in enumerate(rows):
        colour = '#24353D' if position == 5 else COLOURS['PRODIGY']
        f.plot([number(row['lower']),number(row['upper'])],[position,position],color=colour,linewidth=2 if position == 5 else 1.5)
        f.scatter([number(row['known'])],[position],s=20 if position == 5 else 13,color=colour,zorder=3)
    f.axvline(0,color='#ABB7BE',linewidth=.7,linestyle=(0,(3,3)))
    f.set_yticks(range(6),[str(row['seed']) for row in data['seed_joint']]+['Five-seed mean'])
    f.set_ylim(5.7,-.7)
    f.set_xlim(-.65,.08)
    f.set_xticks([-.6,-.3,0])
    f.set_xlabel('Depletion-positive allocation change\n(per 10 candidate slots)')
    f.set_axisbelow(True)
    f.grid(axis='x',color='#E2E7EA',linewidth=.4)
    panel(f,'f','PRODIGY update in each seed')
    return save_figure(figure, output, "Supplementary_Figure_3", "Donor directions, candidate overlap and seed variation", min_font=10)

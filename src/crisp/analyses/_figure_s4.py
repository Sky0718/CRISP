"""Oesophageal and construction comparisons."""
from ._figure_common import apply_figure_style, save_figure
from fractions import Fraction
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.text import Text
from matplotlib.ticker import MaxNLocator
import numpy as np

COLOURS={'DawnRank':'#AB6554','Degree':'#586D7B','positive':'#326B8C','zero':'#C1CBD1','unknown':'#DBCDA8','first':'#326B8C','second':'#AB6554','bridge':'#628778'}

def number(value):
    return float(Fraction(str(value)))

def panel(axis,letter,title):
    axis.text(0,1.08,'('+letter+')',transform=axis.transAxes,fontweight='bold',va='bottom',fontsize=11)
    axis.annotate(title,(0,1.08),xycoords='axes fraction',xytext=(23,0),textcoords='offset points',va='bottom',fontsize=11)

def zero(axis):
    axis.axvline(0,color='#ABB7BE',linewidth=.7,linestyle=(0,(3,3)))
    axis.set_axisbelow(True)
    axis.grid(axis='x',color='#E2E7EA',linewidth=.4)

def draw(data, output):
    family = apply_figure_style()
    plt.rcParams.update({'font.family':family,'font.size':10.5,'axes.labelsize':10.5,'xtick.labelsize':10,'ytick.labelsize':10,'legend.fontsize':10,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.65,'pdf.fonttype':42,'svg.fonttype':'none','axes.unicode_minus':False,'legend.frameon':False,'figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white'})
    full_data = data
    data=full_data['oesophageal_transport']
    summary=full_data['oesophageal_identity_analysis']
    figure=plt.figure(figsize=(220/25.4,285/25.4))
    a=figure.add_axes([.19,.755,.29,.175])
    levels=data['population_levels']
    positions=np.arange(6)
    left=np.zeros(6)
    for field in ['positive','zero','unknown']:
        values=np.array([number(row[field]) for row in levels])
        a.barh(positions,values,left=left,height=.64,color=COLOURS[field],label={'positive':'Depletion-positive','zero':'Measured zero','unknown':'Unmeasured'}[field],edgecolor='white',linewidth=.4)
        left+=values
    if not np.allclose(left,10,rtol=0,atol=1e-12):
        raise ValueError('Complete ten-slot source ledger required')
    a.set_yticks(positions,['G0H0','G1H0','G0H1','G1H1','Degree old','Degree new'])
    a.set_ylim(5.6,-.6)
    a.set_xlim(0,10)
    a.set_xticks([0,5,10])
    a.set_xlabel('Allocation per 10 candidate slots')
    a.legend(loc='upper left',bbox_to_anchor=(-.08,-.25),ncol=2,columnspacing=.8,handlelength=1,handletextpad=.4)
    panel(a,'a','Oesophageal cancer: complete allocation')
    b=figure.add_axes([.67,.755,.285,.175])
    comparisons={row['contrast']:row for row in data['population_contrasts']}
    for position,contrast in enumerate(['Delta','Theta']):
        row=comparisons[contrast]
        outer=data['composition']['outer_envelopes'][contrast]
        b.hlines(position,outer[0],outer[1],color='#7A8992',linewidth=1.1,label='Composition envelope' if position==0 else None)
        b.vlines(outer,position-.07,position+.07,color='#7A8992',linewidth=1.1)
        b.plot([number(row['lower']),number(row['upper'])],[position,position],color=COLOURS['positive'],linewidth=3,label='Shared-label bounds' if position==0 else None)
        b.scatter([number(row['observed'])],[position],s=24,color='#173B50',edgecolors='white',linewidths=.45,zorder=4)
    zero(b)
    b.set_yticks([0,1],['Joint update','Degree-relative'])
    b.set_ylim(1.65,-.65)
    b.set_xlim(-1.5,.25)
    b.set_xticks([-1.5,-1,-.5,0])
    b.set_xlabel('Depletion-positive allocation change\n(per 10 candidate slots)')
    comparison_handles=[plt.Line2D([0],[0],color='#7A8992',linewidth=1.1,label='Donor-composition envelope'),plt.Line2D([0],[0],color=COLOURS['positive'],linewidth=3,label='Missing-label bounds')]
    b.legend(handles=comparison_handles,loc='upper left',bbox_to_anchor=(-.10,-.25),handlelength=1.5,handletextpad=.5)
    panel(b,'b','Joint and relative updates')
    c=figure.add_axes([.19,.42,.125,.20])
    proportion=figure.add_axes([.365,.42,.115,.20],sharey=c)
    for index,row in enumerate(summary['quantity_composition']):
        position=[0,.48,1.55,2.03][index]
        colour=COLOURS[row['method']]
        for axis,key,scale in [(c,'R',1),(proportion,'p',100)]:
            old,new=[scale*number(row[key+state]) for state in ['0','1']]
            axis.plot([old,new],[position,position],color=colour,linewidth=1.3)
            axis.scatter(old,position,s=36,facecolors='white',edgecolors=colour,linewidth=1.15,zorder=3)
            axis.scatter(new,position,s=16,color=colour,zorder=4)
    c.set_yticks([.24,1.79],['Shared core-fitness','Oncogene only'])
    c.set_ylim(2.65,-.6)
    c.set_xlim(0,1.7)
    c.set_xticks([0,.8,1.6])
    c.set_xlabel('Allocation\n(per 10 slots)')
    proportion.set_xlim(0,105)
    proportion.set_xticks([0,50,100])
    proportion.tick_params(axis='y',left=False,labelleft=False)
    proportion.spines['left'].set_visible(False)
    proportion.set_xlabel('Depletion-positive\nfraction (%)')
    handles=[plt.Line2D([0],[0],color=COLOURS['DawnRank'],linewidth=1.5,label='DawnRank'),plt.Line2D([0],[0],color='#4D565B',marker='o',markerfacecolor='white',linewidth=0,markersize=4.5,label='Old'),plt.Line2D([0],[0],color=COLOURS['Degree'],linewidth=1.5,label='Degree'),plt.Line2D([0],[0],color='#4D565B',marker='o',markerfacecolor='#4D565B',linewidth=0,markersize=4.5,label='New')]
    c.legend(handles=handles,loc='upper left',bbox_to_anchor=(-.18,-.23),ncol=2,columnspacing=1,handlelength=1.3,handletextpad=.4)
    panel(c,'c','Allocation and depletion-positive fraction')
    d=figure.add_axes([.67,.42,.285,.20])
    categories=[('core','BOTH_CORE','Shared core'),('core','ORGANOID_ONLY','Organoid only'),('core','CELL_LINE_ONLY','Cell-line only'),('core','NOT_LISTED_IN_EITHER','Neither core'),('role','ONCOGENE','Oncogene only'),('role','TSG','TSG only'),('role','DUAL_ONCOGENE_TSG','Dual role'),('role','OTHER_OR_UNSPECIFIED','Other role'),('role','UNLISTED','Unlisted')]
    positions=[0,1,2,3,4.7,5.7,6.7,7.7,8.7]
    lookup={(row['partition'],row['category']):row for row in summary['category_relative_contrasts']}
    limits=[]
    for position,(partition,category,label) in zip(positions,categories):
        row=lookup[(partition,category)]
        lower,upper,known=[number(row[key]) for key in ['L','U','K']]
        colour=COLOURS['positive'] if partition=='core' else COLOURS['DawnRank']
        d.plot([lower,upper],[position,position],color=colour,linewidth=2)
        d.scatter(known,position,s=19,color=colour,zorder=3)
        limits.extend([lower,upper])
    zero(d)
    d.axhline(3.85,color='#DEE3E6',linewidth=.5)
    d.set_yticks(positions,[row[2] for row in categories])
    d.set_ylim(9.3,-.65)
    span=max(limits)-min(limits)
    d.set_xlim(min(min(limits)-.09*span,-.025),max(max(limits)+.09*span,.025))
    d.xaxis.set_major_locator(MaxNLocator(nbins=4))
    d.set_xlabel('Update minus degree update')
    panel(d,'d','Category changes relative to degree')
    e=figure.add_axes([.19,.085,.29,.205])
    matrix=np.asarray(summary['identity_gain']['matrix'],dtype=int)
    palette=LinearSegmentedColormap.from_list('count',['#FFFFFF','#B8CAD5','#326B8C'])
    e.imshow(matrix,cmap=palette,vmin=0,vmax=int(matrix.max()),aspect='auto')
    for i in range(4):
        for j in range(4):
            e.text(j,i,str(matrix[i,j]),ha='center',va='center',fontsize=11,fontweight='bold' if i!=j and matrix[i,j]>0 else 'normal',color='white' if matrix[i,j]>matrix.max()*.7 else '#233842')
    e.set_xticks(range(4),['Negative','Positive','Zero','Unresolved'],rotation=28,ha='right')
    e.set_yticks(range(4),['Negative','Positive','Zero','Unresolved'])
    e.set_xlabel('Shared-label direction')
    e.set_ylabel('Separate-arm direction')
    e.set_xticks(np.arange(-.5,4,1),minor=True)
    e.set_yticks(np.arange(-.5,4,1),minor=True)
    e.grid(which='minor',color='#D5DFE4',linewidth=.5)
    e.tick_params(which='minor',bottom=False,left=False)
    for spine in e.spines.values():
        spine.set_visible(False)
    panel(e,'e','Direction assessment for all 58 donors')
    f=figure.add_axes([.67,.085,.285,.205])
    for index,(row,key) in enumerate(zip(full_data['contrasts'],['first','second','bridge'])):
        bounds=row['population']
        f.plot([number(bounds['lower']),number(bounds['upper'])],[index,index],color=COLOURS[key],linewidth=2)
        f.scatter([number(bounds['known'])],[index],s=23,color=COLOURS[key],zorder=3)
    zero(f)
    f.set_yticks(range(3),['11.0 to 11.5','11.5 to 12.0','11.5 bridge'])
    f.set_ylim(2.55,-.55)
    f.set_xlim(-1.2,.35)
    f.set_xticks([-1,-.5,0])
    f.set_xlabel('Depletion-positive allocation change\n(per 10 candidate slots)')
    panel(f,'f','Colorectal comparison contexts')
    return save_figure(figure, output, "Supplementary_Figure_4", "Oesophageal and construction comparisons", min_font=10)

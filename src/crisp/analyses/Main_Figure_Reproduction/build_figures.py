"""Render main figures from their numerical source records.

Run: python build_figures.py [--output OUTPUT] [--inputs INPUTS]
Optional --manuscript copies the three final PDFs into a manuscript directory.
Print-scale canvas with editable SVG/PDF output.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
import argparse, csv, json, hashlib, sys, shutil
ROOT = Path(__file__).resolve().parent
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle, FancyArrowPatch
from matplotlib.lines import Line2D
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap
from plot_settings import apply_figure_style

ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,default=ROOT/'output');ap.add_argument('--inputs',type=Path,default=ROOT/'inputs');ap.add_argument('--manuscript',type=Path);args=ap.parse_args()
OUT=args.output;OUT.mkdir(parents=True,exist_ok=True)
apply_figure_style(font='Arial',sizes=(8,7,6))
plt.rcParams.update({'axes.labelpad':3,'axes.titlepad':6,'svg.fonttype':'none','pdf.fonttype':42})
BLUE='#397DA5';ORANGE='#BD6B42';GREY='#ABB3BB';INK='#27333E';PALE='#E6EBEF';ZERO='#C8CED4'
METHODS=('PRODIGY','DawnRank','degree'); LABEL={'PRODIGY':'PRODIGY','DawnRank':'DawnRank','degree':'Degree'}
read=lambda name:list(csv.DictReader((args.inputs/name).open(encoding='utf-8-sig',newline='')))
LEVELS=read('union_levels.csv');ENDPOINTS=read('union_endpoints.csv');BOUNDS=read('fraction_bounds.csv');CAND=read('candidate_identities.csv');ATTR=read('attribution.csv');UNITS=read('decomposition_units.csv')
CRC=json.loads((args.inputs/'crc_categories.json').read_text());ESCA=json.loads((args.inputs/'esca_categories.json').read_text())
RECORDS=[]
def val(x):return float(F(str(x)))
def unique(rows,**key):
    found=[r for r in rows if all(r[k]==v for k,v in key.items())];assert len(found)==1,(key,len(found));return found[0]
def record(fig,panel,**kw):RECORDS.append(dict(figure=fig,panel=panel,**kw))
def axes(fig,rect,title,letter):
    ax=fig.add_axes(rect);ax.set_title(title,loc='left');ax.text(-.033/rect[2],1+.018/rect[3],letter,transform=ax.transAxes,weight='bold',size=9,ha='left',va='bottom');return ax
def quiet(ax):ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0)
def export(fig,name,panels):
    folder=OUT/name;folder.mkdir(exist_ok=True)
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    texts=[t for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text().strip()]
    boxes=[(t,t.get_window_extent(renderer)) for t in texts]
    outside=[t.get_text() for t,b in boxes if b.x0<-.5 or b.y0<-.5 or b.x1>fig.bbox.x1+.5 or b.y1>fig.bbox.y1+.5]
    overlap=[]
    for i,(t,b) in enumerate(boxes):
        for u,c in boxes[i+1:]:
            if b.overlaps(c):overlap.append([t.get_text(),u.get_text()])
    for ext in ('pdf','svg','png'):fig.savefig(folder/f'{name}.{ext}',dpi=300,facecolor='white')
    geom=[]
    for letter,aa in panels.items():
        aa=aa if isinstance(aa,list) else [aa]
        boxes=[a.get_tightbbox(renderer) for a in aa];bb=matplotlib.transforms.Bbox.union(boxes).expanded(1.02,1.05).transformed(fig.dpi_scale_trans.inverted())
        for ext in ('pdf','png'):fig.savefig(folder/f'panel_{letter}.{ext}',dpi=300,bbox_inches=bb)
        geom.append(dict(panel=letter,axes=[list(a.get_position().bounds) for a in aa]))
    (folder/'geometry_review.json').write_text(json.dumps(dict(dimensions_mm=list(fig.get_size_inches()*25.4),panels=geom,outside=outside,text_bbox_overlaps=overlap),indent=2))
    plt.close(fig)
    return folder/f'{name}.pdf'

def overview():
    fig=plt.figure(figsize=(183/25.4,144/25.4))
    anchor=fig.add_axes([.16,.855,.63,.09])
    anchor.text(-.175,1.16,'a',weight='bold',size=9,transform=anchor.transAxes)
    donor_rows=[r for r in UNITS if r['context']=='CRC' and r['partition']=='donor' and r['policy'] in METHODS[:2]]
    donor_ids=sorted({r['unit'] for r in donor_rows},key=int);assert len(donor_ids)==83
    cmap=LinearSegmentedColormap.from_list('donor_fraction_change',[ORANGE,'#F7F7F7',BLUE]);norm=TwoSlopeNorm(vmin=-100,vcenter=0,vmax=100)
    for mi,m in enumerate(METHODS[:2]):
        current=[r for r in donor_rows if r['policy']==m];assert len(current)==83
        assert sum(r['support']=='common' for r in current)==(78 if m=='PRODIGY' else 82)
        for j,donor in enumerate(donor_ids):
            r=unique(current,unit=donor);common=r['support']=='common'
            delta=100*(F(r['f1'])-F(r['f0'])) if common else None
            anchor.add_patch(Rectangle((j+.5,mi-.5),1,1,facecolor=cmap(norm(float(delta))) if common else 'white',edgecolor='white' if common else '#7B838B',linewidth=.15 if common else .35,hatch=None if common else 'xx'))
            record('Figure1','a',context='CRC',method=m,donor_id=donor,donor_index=j+1,support=r['support'],endpoint='donor_fraction_delta_pp',f0_exact=r['f0'] if F(r['R0'])>0 else '',f1_exact=r['f1'] if F(r['R1'])>0 else '',R0_exact=r['R0'],R1_exact=r['R1'],exact=str(delta) if common else '',value=float(delta) if common else '',undefined=not common)
    anchor.set_xlim(.5,83.5);anchor.set_ylim(1.5,-.5);anchor.set_yticks([0,1],METHODS[:2],fontsize=7)
    anchor.set_xticks([1,20,40,60,83]);anchor.tick_params(axis='both',length=0,pad=3);anchor.set_xlabel('Donor',fontsize=7,labelpad=1)
    for spine in anchor.spines.values():spine.set_visible(False)
    color_axis=fig.add_axes([.845,.899,.11,.018],label='<colorbar>')
    cb=fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm,cmap=cmap),cax=color_axis,orientation='horizontal',ticks=[-100,0,100]);cb.outline.set_linewidth(.35)
    color_axis.set_title('Change in fraction (pp)',fontsize=7,pad=4);color_axis.tick_params(labelsize=6,length=2,pad=2)
    fig.legend(handles=[Patch(facecolor='white',edgecolor='#7B838B',hatch='xx',linewidth=.35,label='Undefined')],loc='lower left',bbox_to_anchor=(.832,.837),fontsize=6,handlelength=1.3,handletextpad=.5)
    b=axes(fig,[.085,.465,.38,.285],'Cancer-listed selection expands','b')
    xs=[0,1,3,4,6,7]
    for mi,m in enumerate(METHODS):
        for si in range(2):
            r=unique(LEVELS,context='CRC',policy=m,state=str(si));bottom=0
            for field,color in [('H',BLUE),('Z',ZERO),('U','white')]:
                v=val(r[field]);b.bar(xs[2*mi+si],v,bottom=bottom,width=.72,color=color,edgecolor=GREY if field=='U' else 'white',lw=.45);bottom+=v
                record('Figure1','b',context='CRC',method=m,state=si,endpoint=field,exact=r[field],value=v)
            assert abs(bottom-val(r['R']))<1e-10
        b.text(3*mi+.5,-.21,LABEL[m],ha='center',va='top',transform=b.get_xaxis_transform(),size=7)
    b.set_xticks(xs,['11.0','11.5']*3);b.tick_params(axis='x',length=0);b.set_xlim(-.7,7.8);b.set_ylim(0,5.3);b.set_yticks([0,1,2,3,4,5]);b.set_ylabel('Listed slots per ten offered')
    b.legend(handles=[Patch(facecolor=BLUE,label='Positive'),Patch(facecolor=ZERO,label='Measured zero'),Patch(facecolor='white',edgecolor=GREY,label='Unknown')],loc='upper left',bbox_to_anchor=(-.08,-.27),ncol=3,fontsize=6,handlelength=1,columnspacing=.9)
    c=axes(fig,[.635,.465,.33,.285],'Positive fraction declines','c')
    for j,m in enumerate(METHODS):
        e=unique(ENDPOINTS,context='CRC',contrast=m,endpoint='delta_H_over_R');r=unique(BOUNDS,context='CRC',policy=m);y=2-j;col=BLUE if m!='degree' else '#7B838B'
        lo,hi=100*val(e['envelope_lower']),100*val(e['envelope_upper']);c.plot([lo,hi],[y,y],color=col,lw=.9);c.plot([lo,lo],[y-.06,y+.06],color=col,lw=.7);c.plot([hi,hi],[y-.06,y+.06],color=col,lw=.7)
        c.plot([100*val(r['sharp_lower']),100*val(r['sharp_upper'])],[y,y],color=col,lw=4,solid_capstyle='butt');c.scatter([100*val(e['estimate'])],[y],color=col,edgecolor='white',linewidth=.55,s=24,zorder=4)
        record('Figure1','c',context='CRC',method=m,endpoint='delta_fraction_pp',estimate=100*val(e['estimate']),donor_lower=lo,donor_upper=hi,sharp_lower=100*val(r['sharp_lower']),sharp_upper=100*val(r['sharp_upper']))
    c.axvline(0,color=GREY,lw=.65,ls='--');c.set_xlim(-18,1);c.set_xticks([-15,-10,-5,0]);c.set_ylim(-.45,2.75);c.set_yticks([2,1,0],[LABEL[m] for m in METHODS]);quiet(c);c.set_xlabel('Change in positive fraction\n(percentage points)')
    c.legend(handles=[Line2D([0],[0],lw=.9,color=INK,label='Donor envelope'),Line2D([0],[0],lw=4,color=INK,label='Label bounds')],loc='upper left',bbox_to_anchor=(-.22,-.29),ncol=2,handlelength=1.3,fontsize=6,columnspacing=1)
    ds=[]
    for mi,m in enumerate(METHODS[:2]):
        d=fig.add_axes([.19+mi*.45,.070,.30,.205]);ds.append(d);d.set_title(m,loc='left',fontsize=7)
        if mi==0:fig.text(.05,.310,'d',size=9,weight='bold');fig.text(.09,.310,'Listed positive gains accompany larger unlisted losses',size=8)
        for j,(cat,partition,key,color) in enumerate([('Cancer-listed','role',None,BLUE),('Unlisted','role','UNLISTED',GREY),('All genes','ALL','ALL',INK)]):
            found=[r for r in ATTR if r['method']==('DawnRank_strict' if m=='DawnRank' else m) and r['partition']==partition and (r['category']!='UNLISTED' if key is None else r['category']==key)]
            vv=[sum((F(r[f]) for r in found),F(0)) for f in ('known','lower','upper')];obs,lo,hi=map(float,vv);y=2-j
            d.plot([lo,hi],[y,y],color=color,lw=3,solid_capstyle='butt');d.scatter([obs],[y],s=24,color=color,edgecolor='white',lw=.5,zorder=3)
            record('Figure1','d',context='CRC',method=m,category=cat,endpoint='delta_positive_slots',exact=str(vv[0]),value=obs,sharp_lower=lo,sharp_upper=hi)
        d.axvline(0,color=GREY,lw=.6,ls='--');d.set_xlim(-1.17,.25);d.set_ylim(-.45,2.55);d.set_xticks([-1,-.5,0]);d.set_yticks([2,1,0],['Cancer-listed','Unlisted','All genes']);quiet(d);d.set_xlabel('Change in positive slots per ten offered',fontsize=7)
    return export(fig,'Fig_Overview',dict(a=[anchor,color_axis],b=b,c=c,d=ds))

def candidates():
    pos=[r for r in CAND if r['label_group']=='positive' and F(r['net'])!=0];ids=sorted({(r['gene'],r['model_id']) for r in pos});genes=sorted({g for g,s in ids});assert len(ids)==34 and len(genes)==12
    vals={(m,g,k):sum((F(r[k])*F(r['donor_weight']) for r in pos if r['method']==m and r['gene']==g),F(0)) for m in METHODS[:2] for g in genes for k in ('gain','loss')}
    assert sum(vals['PRODIGY',g,'gain'] for g in genes)==F(159,830)
    assert sum(vals['PRODIGY',g,'loss'] for g in genes)==F(53,415)
    assert sum(vals['DawnRank',g,'gain'] for g in genes)==F(8,83)
    assert sum(vals['DawnRank',g,'loss'] for g in genes)==F(1,83)
    fig=plt.figure(figsize=(183/25.4,184/25.4))
    method_colors={'PRODIGY':'#665D8C','DawnRank':'#35424B'}
    a=axes(fig,[.10,.39,.28,.50],'Gross gains and losses','a')
    for mi,m in enumerate(METHODS[:2]):
        for j,g in enumerate(genes):
            y=11-j+(.16 if mi==0 else -.16);gain=float(vals[m,g,'gain']);loss=float(vals[m,g,'loss']);color=method_colors[m]
            a.barh(y,gain,height=.25,color=color,alpha=.9);a.barh(y,-loss,height=.25,color=color,alpha=.35,edgecolor=color,lw=.5)
            if gain==0 and loss==0:a.scatter([0],[y],s=6,color=color,zorder=4)
            record('Figure3','a',method=m,gene=g,endpoint='gross_positive_allocation',gain_exact=str(vals[m,g,'gain']),loss_exact=str(vals[m,g,'loss']),gain=gain,loss=loss)
    a.axvline(0,color=GREY,lw=.6);a.set_xlim(-.085,.065);a.set_xticks([-.08,-.04,0,.04]);a.set_ylim(-.6,12);a.set_yticks(range(11,-1,-1),genes,fontstyle='italic');quiet(a);a.set_xlabel('Positive slots per ten offered')
    a.text(.0,.995,'Loss',transform=a.transAxes,size=7,ha='left',va='top');a.text(1,.995,'Gain',transform=a.transAxes,size=7,ha='right',va='top')
    fig.legend(handles=[Patch(facecolor=method_colors['PRODIGY'],label='PRODIGY'),Patch(facecolor=method_colors['DawnRank'],label='DawnRank')],loc='upper left',bbox_to_anchor=(.095,.965),fontsize=7,ncol=2,handlelength=1.0,columnspacing=1.2)
    # All 34 changed identities, not a selected best-case subset.
    b=fig.add_axes([.50,.105,.47,.785]);b.set_xlim(0,1);b.set_ylim(-.5,34.8);b.set_axis_off();fig.text(.46,.919,'b',size=9,weight='bold');fig.text(.50,.919,'Each changed model–gene experiment',size=8)
    head=[(.0,'Gene / model'),(.39,'P'),(.475,'D'),(.57,'Protein change'),(.96,'LFC')]
    for x,t in head:b.text(x,34.3,t,size=7,ha='right' if t=='LFC' else 'left',va='center')
    cmap=LinearSegmentedColormap.from_list('loss_gain',['#BA6745','#F7F7F7','#337EA5']);norm=TwoSlopeNorm(vmin=-1,vcenter=0,vmax=1)
    lookup={(r['method'],r['gene'],r['model_id']):r for r in CAND}
    for i,(g,s) in enumerate(ids):
        y=33-i;source=next(r for r in pos if r['gene']==g and r['model_id']==s)
        if i%2==0:b.add_patch(Rectangle((-.01,y-.47),1.02,.94,color='#F3F5F6',lw=0))
        b.text(0,y,g,size=6,ha='left',va='center',fontstyle='italic');b.text(.195,y,s[4:],size=6,ha='left',va='center')
        for mi,m in enumerate(METHODS[:2]):
            r=lookup.get((m,g,s));v=F(r['net']) if r else F(0);x=.39+mi*.085
            b.add_patch(Rectangle((x-.008,y-.40),.070,.80,facecolor=cmap(norm(float(v))),edgecolor='white',lw=.5))
            txt=f'{float(v):g}' if r else '–';b.text(x+.027,y,txt,ha='center',va='center',size=6,color='white' if abs(float(v))>.7 else INK)
            record('Figure3','b',method=m,gene=g,model=s,endpoint='model_allocation_change',exact=str(v),value=float(v),both_states_absent=(r is None),protein_changes=source['source_protein_changes'],LFC=source['supplied_LFC'])
        protein=source['source_protein_changes'].replace('p.','').replace('|',' / ')
        b.text(.57,y,protein,size=6,ha='left',va='center');b.text(.96,y,f"{float(source['supplied_LFC']):.2f}",size=6,ha='right',va='center')
    b.text(.0,-1.35,'P: PRODIGY   D: DawnRank   Model: SIDM suffix',size=6,ha='left');b.text(.0,-2.35,'Cells: allocation change (−1 loss, +1 gain)',size=6,ha='left');b.text(.0,-3.35,'–: absent in both lists; 0: unchanged allocation',size=6,ha='left')
    c=axes(fig,[.10,.10,.28,.20],'Opposite KRAS screen contexts','c')
    for i,(s,marker) in enumerate([('SIDM01266','o'),('SIDM02231','s')]):
        for mi,m in enumerate(METHODS[:2]):
            r=unique(CAND,model_id=s,gene='KRAS',method=m);y=3-i*2-mi*.55;old,new=val(r['old']),val(r['new']);color=method_colors[m]
            c.plot([old,new],[y,y],color=color,lw=1.1);c.scatter([old],[y],s=21,facecolor='white',edgecolor=color,lw=.8,zorder=4);c.scatter([new],[y],s=21,color=color,marker='o',zorder=5)
            record('Figure3','c',method=m,gene='KRAS',model=s,old=r['old'],new=r['new'],LFC=r['supplied_LFC'],label_group=r['label_group'])
        r=unique(CAND,model_id=s,gene='KRAS',method='PRODIGY');c.text(-.08,3-i*2+.50,f"{s} / {r['source_protein_changes'].replace('p.','')}\nLFC {float(r['supplied_LFC']):+.2f}; {r['label_group']}",size=6,ha='left',va='bottom')
    c.set_xlim(-.08,1.10);c.set_ylim(-.10,4.1);c.set_xticks([0,.5,1]);c.set_yticks([]);quiet(c);c.set_xlabel('Allocation in a ten-slot list',fontsize=7)
    c.legend(handles=[Line2D([],[],marker='o',color=INK,markerfacecolor='white',ls='',label='Old'),Line2D([],[],marker='o',color=INK,ls='',label='Updated')],loc='upper left',bbox_to_anchor=(-.02,-.23),ncol=2,fontsize=6,handletextpad=.3)
    return export(fig,'Fig_Candidates',dict(a=a,b=b,c=c))

def context():
    fig=plt.figure(figsize=(183/25.4,145/25.4))
    controls=read('graph_controls.csv');controls=[r for r in controls if r['mode'] in ('H0','H1')];assert len(controls)==18
    a=axes(fig,[.105,.20,.43,.64],'More listed genes, lower positive fraction','a')
    groupcols={'original':INK,'low':BLUE,'high':ORANGE,'threshold':'#81719C'}
    for r in controls:
        group=r['pair'].split('_')[0];color=groupcols[group];x=val(r['delta_R']);y=100*val(r['observed_delta']);lo=100*val(r['sharp_lower']);hi=100*val(r['sharp_upper'])
        a.plot([x,x],[lo,hi],color=color,lw=.75,alpha=.65,zorder=1);a.scatter([x],[y],color=color,marker='o' if r['mode']=='H0' else 's',s=23,edgecolor='white',linewidth=.4,zorder=3)
        record('Figure5','a',context='CRC',method='DawnRank',graph_pair=r['pair'],fixed_inputs=r['mode'],endpoint='union_graph_contrast',delta_listed_slots=x,observed_fraction_delta_pp=y,sharp_lower_pp=lo,sharp_upper_pp=hi)
    a.axhline(0,color=GREY,lw=.75,ls='--');a.set_xlim(.85,1.39);a.set_xticks([.9,1.0,1.1,1.2,1.3]);a.set_ylim(-8.4,.45);a.set_yticks([-8,-6,-4,-2,0]);a.set_xlabel('Change in listed slots per ten offered');a.set_ylabel('Change in positive fraction (percentage points)')
    a.legend(handles=[Line2D([],[],color=groupcols[g],marker='o',ls='',label=label) for g,label in [('original','Original threshold'),('low','65,783 edges'),('high','119,423 edges'),('threshold','Alternative thresholds')]],loc='upper left',fontsize=6,handletextpad=.3)
    fig.legend(handles=[Line2D([],[],color=INK,marker='o',ls='',label='Old inputs fixed'),Line2D([],[],color=INK,marker='s',ls='',label='Updated inputs fixed'),Line2D([],[],color=INK,lw=.8,label='Shared-label bounds')],loc='lower left',bbox_to_anchor=(.075,.06),ncol=1,fontsize=7,handletextpad=.4)
    top=[]
    for j,(endpoint,title,xlab,scale,lims,ticks) in enumerate([
        ('delta_R','Listed selection expands','Change in listed slots',1,(-.1,1.85),[0,.5,1,1.5]),
        ('delta_H_over_R','Observed-positive fraction','Change in fraction (pp)',100,(-12,6),[-10,-5,0,5]),
        ('delta_H','Positive gains remain modest','Change in positive slots',1,(-.19,.36),[-.1,0,.1,.2,.3])]):
        ax=axes(fig,[.73,.67-j*.27,.24,.17],title,chr(98+j));top.append(ax)
        for i,(co,color) in enumerate([('CRC',BLUE),('ESCA',ORANGE)]):
            r=unique(ENDPOINTS,context=co,contrast='DawnRank',endpoint=endpoint);e=scale*val(r['estimate']);lo=scale*val(r['envelope_lower']);hi=scale*val(r['envelope_upper']);y=1-i
            ax.plot([lo,hi],[y,y],lw=1.2,color=color);ax.scatter([e],[y],s=24,color=color,edgecolor='white',lw=.5,zorder=4)
            record('Figure5',chr(98+j),context=co,method='DawnRank',endpoint=endpoint,estimate=e,donor_lower=lo,donor_upper=hi,scale=scale)
        ax.set_xlim(*lims);ax.set_xticks(ticks);ax.set_ylim(-.55,1.55);ax.set_yticks([1,0],['Colorectal','Oesophageal']);quiet(ax);ax.axvline(0,color=GREY,lw=.6,ls='--');ax.set_xlabel(xlab,fontsize=7)
    return export(fig,'Fig_4',dict(a=a,b=top[0],c=top[1],d=top[2]))

paths=[overview(),candidates(),context()]
keys=sorted(set().union(*(r.keys() for r in RECORDS)))
with (OUT/'source_data.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(RECORDS)
manifest={'backend':'Python/matplotlib','records':len(RECORDS),'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(args.inputs.iterdir()) if p.is_file()},'row_selection':{'Figure1':'Panel a: all 83 CRC donors, in ascending numeric donor-ID order, for both primary methods (166 cells). Plot 100*(f1-f0) only on common listed-allocation support (78 PRODIGY, 82 DawnRank donors); six non-common cells are hatched as undefined, not assigned zero. Panels b-d: CRC all three union policies; two primary methods for whole-list partition.','Figure3':'Every positive oncogene-only identity with nonzero old/new change for either primary method: 34 identities, 12 genes. Two KRAS source cases also displayed.','Figure5':'All18 fixed-input DawnRank graph contrasts; DawnRank same union endpoints in both contexts. The9degree controls retained in the supporting data but not plotted in mainFigure5.'},'intervals':'Donor-composition sensitivity envelopes retain original endpoint families. Fixed-cohort shared-label completion bounds are distinct intervals, not confidence intervals.'}
def table_summary(path):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        rows = sum(1 for _ in csv.DictReader(stream))
    return {'sha256': hashlib.sha256(path.read_bytes()).hexdigest(), 'rows_input': rows}
manifest['source_tables']={p.name:table_summary(p) for p in args.inputs.glob('*.csv')}
(OUT/'source_manifest.json').write_text(json.dumps(manifest,indent=2))
if args.manuscript:
    for p in paths:
        target=args.manuscript/'Figure'/p.stem/(p.stem+'.pdf');target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
print(json.dumps({'figures':[str(p) for p in paths],'records':len(RECORDS)},indent=2))

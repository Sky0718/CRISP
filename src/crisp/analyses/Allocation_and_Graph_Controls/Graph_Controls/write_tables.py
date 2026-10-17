import argparse
import csv
from fractions import Fraction
from pathlib import Path

parser = argparse.ArgumentParser(description='Generate graph-control LaTeX tables from numerical records.')
parser.add_argument('--source', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
ROOT, OUTPUT = args.source.resolve(), args.output.resolve()
if OUTPUT.exists() or OUTPUT.is_relative_to(ROOT):
    raise ValueError('Choose a new output directory outside the input directory.')
OUTPUT.mkdir(parents=True)
with (ROOT/'population_contrasts.csv').open(encoding='utf-8',newline='') as handle:
    rows = list(csv.DictReader(handle))
lookup = {(row['contrast'],row['category']):row for row in rows}


def val(row,key,scale=1):
    return float(Fraction(row[key]))*scale


def number(value,digits=3):
    return f'{value:+.{digits}f}'.replace('+0.000','0.000')


def interval(low,high,digits=3):
    return '$['+number(low,digits)+', '+number(high,digits)+']$'


def span(family,category,key,scale=1,digits=3):
    data = [val(row,key,scale) for row in rows if row['pair'].startswith(family) and row['mode']!='DEGREE' and row['category']==category]
    return interval(min(data),max(data),digits)


def table(caption,label,headers,body,notes,columns):
    return '\n'.join([r'\begin{table}[H]',r'\begin{minipage}{\textwidth}',r'\caption{'+caption+r'}\label{'+label+'}',r'\begin{singlespace}\fontsize{9}{11}\selectfont',r'\setlength{\tabcolsep}{2.5pt}\renewcommand{\arraystretch}{1.10}',r'\begin{tabular*}{\textwidth}{@{}'+columns+r'@{}}',r'\toprule',' & '.join(headers)+r' \\ \midrule',*body,r'\bottomrule',r'\end{tabular*}',r'\par\vspace{3pt}\noindent\textit{Note.} '+notes,r'\end{singlespace}',r'\end{minipage}',r'\end{table}',''])


families=[('original_T700',r'Original $>700$'),('low_count','65,783 edges each'),('high_count','119,423 edges each'),('threshold_T650',r'Threshold $>650$'),('threshold_T750',r'Threshold $>750$')]
mainbody=[]
for family,name in families:
    values=[span(family,category,key,scale,digits) for category,key,scale,digits in [('BOTH_CORE','delta_R',1,3),('ONCOGENE','delta_R',1,3),('BOTH_CORE','delta_fraction_lower',100,2),('ONCOGENE','delta_fraction_lower',100,2)]]
    mainbody.append(' & '.join([name]+values)+r' \\')
main=table('DawnRank category changes under finite graph-construction controls.','tab:graph_controls_main',[r'\textbf{Construction}',r'\textbf{$\Delta R_{\rm core}$}',r'\textbf{$\Delta R_{\rm onco}$}',r'\textbf{$\Delta p_{\rm core}$ (pp)}',r'\textbf{$\Delta p_{\rm onco}$ (pp)}'],mainbody,r'New minus old, within each fixed consumed-input state. Entries are descriptive minima and maxima across $H_0/H_1$ and, for equal-count comparisons, three fixed boundary-tie choices; they are not confidence intervals. $R$ is donor-weighted allocation per ten offered slots. $p$ is the pooled within-category measured-positive fraction; pp denotes percentage points. Every selected identity in these two categories was measured. Core is the shared core-fitness category; onco is oncogene-only. All 85 models and 83 donors are retained. Equal edge counts leave edge identities, degrees and isolated-node counts unequal.',r'l@{\extracolsep{\fill}}rrrr')
(OUTPUT/'Main_Graph_Controls_Table.tex').write_text(main,encoding='utf-8')

metadata=list(csv.DictReader((ROOT/'graph_metadata.csv').open(encoding='utf-8',newline='')))
graphbody=[]
for row in metadata:
    name=row['graph'].replace('_',r'\_')
    graphbody.append(' & '.join([r'\texttt{'+name+'}',f"{int(row['edges']):,}",f"{int(row['isolates']):,}",row['maximum_degree']])+r' \\')
graph=table('Complete graph constructions retained in the finite control experiment.','tab:graph_control_constructions',[r'\textbf{Graph}',r'\textbf{Edges}',r'\textbf{Isolates}',r'\textbf{Maximum degree}'],graphbody,r'$G_0/G_1$ identify the original/new physical sources. T650, T700 and T750 retain integer scores strictly above the indicated threshold. M65783 and M119423 retain exactly that many highest-score edges; S identifies the fixed boundary-tie seed. All graphs retain 7,399 genes. G0\_T700 and G1\_T700 reuse the accepted original rankings. No graph was removed based on its results.',r'l@{\extracolsep{\fill}}rrr')


def control_table(mode):
    selected=[]
    for row in rows:
        if row['category']=='ALL' and ((row['mode']=='DEGREE') if mode=='degree' else (row['mode']!='DEGREE')):
            selected.append(row)
    body=[]
    for row in selected:
        name=row['contrast']
        core=lookup[(name,'BOTH_CORE')]
        onco=lookup[(name,'ONCOGENE')]
        label=row['pair'].replace('original_T700','Original700').replace('low_count_S','Low/').replace('high_count_S','High/').replace('threshold_T','T')
        if mode!='degree':
            label+=r'/$'+row['mode'][0]+'_'+row['mode'][1]+'$'
        fields=[label,interval(val(row,'lower'),val(row,'upper')), '$'+number(val(core,'delta_R'))+'$', '$'+number(val(onco,'delta_R'))+'$', '$'+number(val(core,'delta_fraction_lower',100),2)+'$', '$'+number(val(onco,'delta_fraction_lower',100),2)+'$']
        assert core['delta_fraction_lower']==core['delta_fraction_upper'] and onco['delta_fraction_lower']==onco['delta_fraction_upper']
        body.append(' & '.join(fields)+r' \\')
    method='Same-query degree' if mode=='degree' else 'DawnRank'
    return table(method+' results for every prespecified graph comparison.','tab:graph_controls_'+mode,[r'\textbf{Pair/input}',r'\textbf{Positive-mass bounds}',r'\textbf{$\Delta R_{\rm core}$}',r'\textbf{$\Delta R_{\rm onco}$}',r'\textbf{$\Delta p_{\rm core}$}',r'\textbf{$\Delta p_{\rm onco}$}'],body,r'All contrasts are new minus old at fixed input. Low pairs the old source restricted to 65,783 edges with unchanged G1\_T700; High pairs unchanged G0\_T700 with the new source expanded to 119,423 edges. Numbers following Low/High are tie seeds. Positive-mass intervals are finite-cohort sharp shared-label bounds, not confidence intervals. $\Delta p$ is in percentage points and is point-identified here because both focal categories have fully measured selected labels. All rows use equal-donor weights, complete original mutation queries and ten slots. Degree has no input-state dependence.',r'l@{\extracolsep{\fill}}rrrrr')

(OUTPUT/'Supplementary_Graph_Control_Tables.tex').write_text(
    graph+control_table('dawn')+control_table('degree'), encoding='utf-8')
print('Wrote graph-control tables from numerical endpoint records.')

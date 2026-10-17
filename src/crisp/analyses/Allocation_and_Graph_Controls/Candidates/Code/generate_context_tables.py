"""Format complete candidate tables from validated CSVs, without case selection."""
import argparse
import csv
from fractions import Fraction as F
from pathlib import Path


def rows(path):
    with path.open(encoding='utf-8',newline='') as f:
        return list(csv.DictReader(f))


def esc(x):
    return str(x).replace('_',r'\_').replace('&',r'\&').replace('%',r'\%')


def wt(old,new):
    return '$'+format(float(F(old)),'.2g')+r'\to'+format(float(F(new)),'.2g')+'$'


def table(body,caption,label,cols,heads):
    return '\n'.join([r'\begin{table}[H]',r'\centering',r'\begin{singlespace}',
        r'\small',r'\setlength{\tabcolsep}{4pt}',r'\begin{tabular}{'+cols+'}',r'\toprule',
        ' & '.join(heads)+r' \\',r'\midrule',*body,r'\bottomrule',r'\end{tabular}',
        r'\caption{'+caption+'}',r'\label{'+label+'}',r'\end{singlespace}',r'\end{table}',''])


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--source-data',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args()
    rr=rows(a.source_data/'three_node_full_cohort_context.csv')
    changed=[r for r in rr if any(F(r[m+'_old'])!=F(r[m+'_new']) for m in ['PRODIGY','DawnRank'])]
    assert len(changed)==12
    body=[]
    for r in sorted(changed,key=lambda r:(r['gene'],r['model_id'])):
        body.append(' & '.join([esc(r['gene']),esc(r['model_id']),esc(r['source_protein_changes'].replace('|','; ')),
            r['binary_source_value'],format(float(r['supplied_LFC']),'.3f'),
            wt(r['PRODIGY_old'],r['PRODIGY_new']),wt(r['DawnRank_old'],r['DawnRank_new'])])+r' \\')
    text=table(body,r'\textbf{Every changed allocation in the source-anchored three-node diagnostic.} '+
        r'The table includes all 12 changed model--gene combinations among 255 combinations, including the zero-labelled KRAS gain. '+
        r'Call is the original BAGEL2 depletion call; LFC is the supplied continuous log-fold change. P and D show old-to-new PRODIGY and DawnRank allocations, respectively; PRODIGY weights average the original five seeds. '+
        r'Protein changes are copied from the source variant annotation without assigning an activating or inactivating function. Full assay-sample and donor identifiers, variant allele fractions and genomic alleles are retained in the data.',
        'tab:node_changed_context','lllrrrr',['Gene','Model','Source allele','Call','LFC','P','D'])
    (a.output/'Supplementary_Node_Changes_Table.tex').write_text(text,encoding='utf-8')
    rr=rows(a.source_data/'changed_positive_34_context.csv')
    assert len(rr)==34 and len({(r['model_id'],r['gene']) for r in rr})==34
    body=[]
    for r in sorted(rr,key=lambda r:(r['gene'],r['model_id'])):
        body.append(' & '.join([esc(r['gene']),esc(r['model_id']),esc(r['source_protein_changes'].replace('|','; ')),
            format(float(r['supplied_LFC']),'.3f'),
            wt(r['PRODIGY_old_weight'],r['PRODIGY_new_weight']),wt(r['DawnRank_old_weight'],r['DawnRank_new_weight'])])+r' \\')
    text=table(body,r'\textbf{Complete allele context for the 34 changed-positive identities.} '+
        r'Every row is query-eligible and has an original positive BAGEL2 call. A change under either method admits the identity; the other method may be unchanged. '+
        r'P and D show old-to-new PRODIGY and DawnRank allocations, respectively. Protein annotations, including unresolved p.? entries, are retained verbatim; their functional direction remains unassigned. '+
        r'LFC values are descriptive source measurements, rounded here to three decimals and retained at full supplied precision in the data. The two CTNNB1 annotations in SIDM01268 are two source records for one model--gene identity.',
        'tab:changed_positive_alleles','lllrrr',['Gene','Model','Source allele','LFC','P','D'])
    (a.output/'Supplementary_Changed_Positive_Alleles_Table.tex').write_text(text,encoding='utf-8')
    print('Formatted all 12 changed node identities and all 34 changed-positive identities.')


if __name__=='__main__':
    main()

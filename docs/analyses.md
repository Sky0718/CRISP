# Analyses from supplied data

The commands in this guide calculate allocation, functional-screen and
candidate summaries from the accompanying `Supplementary_Data` records. Extract
all 16 supplementary ZIP parts together, keeping their shared directory tree.
Run the examples from the code repository root and choose a new output directory
outside both the code and input-data directories.

```console
python crisp.py reproduce --data ../Supplementary_Data --output ../CRISP_results
python crisp.py reproduce --data ../Supplementary_Data --output ../CRISP_results_with_continuous --union-readouts
python crisp.py analysis union --data ../Supplementary_Data --output ../CRISP_union
```

The first command derives cancer-listed binary summaries, exact shared-label
bounds and donor/role decompositions at budget ten, and checks figure-source hashes.
`--union-readouts` adds continuous depletion summaries, paired donor resampling
and exact measured-label sensitivity. Use the `union` component for the additional
budget and graph-control comparisons. The installed `crisp` command accepts the
same arguments as `python crisp.py`.

## Components

Use `python crisp.py analysis COMPONENT --data DATA --output OUTPUT`.

| Component | Calculation and required supplementary component |
|---|---|
| `union` | Completed-fraction bounds, donor/role decomposition, budgets 1/5/10/20 and matched graph controls; `Union_Sensitivity_Analysis` |
| `union-ledger` | Complete CRC and oesophageal union allocation and donor-resampling tables; `Reader_Data/Data_10/Union_Analysis` |
| `continuous` | Continuous depletion, observed-label flips and four-role summaries; `Union_Validation` and the allocation inputs in `Reader_Data/Data_10/Union_Analysis` |
| `categories` | Four-cell category levels/contrasts and oncogene identity flows; `Colorectal_Functional_Analysis`, `functional_content` and allocation inputs |
| `composition` | Symmetric within-donor, reweighting, entry and exit terms plus focal continuous envelopes; `functional_content`, `Category_Sensitivity` and `Colorectal_Functional_Analysis` |
| `candidates` | Candidate identities, source variants, depletion values and three-node full-cohort context; `Allocation_and_Graph_Controls/Candidates` |
| `structure` | Category connectivity, query opportunities and allocation flows across degree-change strata; `Allocation_and_Graph_Controls/Graph_Controls` |
| `opportunity` | Opportunity and selected-flow denominator comparison; `Variant_and_Network_Context/Structure` |
| `alleles` | Exact allele/CIViC matching from cached public source records; `Variant_and_Network_Context/Alleles` |
| `graph-endpoints` | Exact reconstruction of graph-control allocation and endpoint tables; `Allocation_and_Graph_Controls/Graph_Controls` |
| `graph-numerical` | Sparse fixed-point, convergence, ordering and cutoff-bound checks; native matrices and reference vectors in `Allocation_and_Graph_Controls/Graph_Controls` |
| `figures` | Main and supplementary figure panels from their numerical source records; `Main_Figure_Reproduction`, `Figure_Source_Data` and donor-composition records |

Each component writes results and its numerical comparisons below the selected
output directory. Hash checks compare inputs with their supplied manifests;
reference-table comparisons use the tolerances specified by each analysis.

The `alleles` component copies the cached source records into its output directory
and summarises them without network access. The separate `prepare` command in
`src/crisp/analyses/Variant_and_Network_Context/Alleles/annotate_coverage.py`
reads original provider VCFs using their source-path mappings. Its `query` command
uses the resulting HGVS queries to retrieve public allele records. The supplied
cache is sufficient for the `alleles` command shown here.

## Dependencies

Use Python 3.10 or later and install `requirements.txt` for the numerical
components. Continuous, composition, category and union-ledger calculations use
NumPy; `graph-numerical` also uses SciPy. Install `requirements_figures.txt` for
figure rendering. Arial is preferred, with an installed sans-serif fallback when
it is unavailable. Native R calculations have a separate environment described
in [Native calculations](native.md).

## Native graph-control calculations

`src/crisp/analyses/Allocation_and_Graph_Controls/Graph_Controls/prepare_graphs.py`
constructs the prespecified graphs from separately obtained
STRING releases. Its required options are `--genes`, `--old-links`, `--old-info`,
`--new-links`, `--new-info`, `--edge-reference` and `--output`. Inputs include the
ordered common gene domain and complete source network records.

The R entry point
`src/crisp/analyses/Allocation_and_Graph_Controls/Graph_Controls/Reproduce_Native.R`
accepts four positional arguments:

1. Graph work directory containing `GRAPH_MANIFEST.json`, `graph_metadata.csv`
   and `graphs/`.
2. External graph-control data directory containing `reference_allocations/numerical/`.
3. Prepared input RDS file at
   `Supplementary_Data/Allocation_and_Graph_Controls/Graph_Controls/runtime/prepared_input.rds`,
   containing the gene domain, models, input states, mutation queries, graph
   operators and seed aliases.
4. `src/crisp/analyses/Allocation_and_Graph_Controls/Graph_Controls`, which contains
   `runtime/native/functions.R` and `runtime/native/original_adapter.R`.

Use a separate graph work directory for the newly generated matrices and score
vectors. This entry point requires R 4.5.2 and Matrix 1.7.4. The parameters retain all 85
models, 83 donors, 7,399 genes, both processed-input states and the five specified
seed aliases. The numerical checks use the fixed tolerances and complete query
orders defined in the source. Raw provider datasets, prepared inputs and derived
matrices are external data, not GitHub source files.

## Source layout

`src/crisp/analyses/runner.py` dispatches the components. Component directories
contain calculation code, independent numerical checks and scientific protocols.
`src/crisp/statistics/` supplies the shared exact-allocation and resampling
primitives. Source extraction, table generation and graph-construction scripts
also retain their individual command-line arguments for analyses using original
provider records.

# Native calculations and allocation analysis

CRISP connects network-based cancer-gene rankings to model-matched screen profiles. This guide covers molecular preparation, native PRODIGY and DawnRank calculations, exact allocation summaries, numerical precision, and cohort-level statistics. The [main README](../README.md) describes the supplied-record analysis interface.

## Command and path conventions

From the repository root:

```bash
python crisp.py native COMMAND CONFIG.json
```

The installed `crisp` command accepts the same arguments. Configuration paths are relative to the JSON file containing them. Executable names such as `Rscript` resolve through `PATH`; explicit executable and library paths are also configuration-relative. Select a new output directory or file for each calculation.

For the accompanying data, use the configurations in `CRISP_Base_Data/config`.
Their command index is `CRISP_Base_Data/CONFIGURATIONS.json`. The examples below
assume that the code checkout is named `CRISP` and is beside `CRISP_Base_Data`
and `Supplementary_Data`, as shown in [Data and dependencies](data.md).

The files in this repository's [config](../config) directory specify software
interfaces and shared rules. Their input and output paths are templates for a
corresponding dataset. Before invoking a native command, check all input paths
and choose a fresh `output` destination. Copying a configuration to another
directory also changes the base of its relative paths. Preserve model and donor
identifiers, source eligibility rules, ordered gene axes and numerical parameters
when adapting a study comparison.

## Environments

Molecular preparation uses a separate numerical environment with NumPy `1.26.4` and pandas `2.2.3`. The `prepare-molecular` command checks these versions. Activate that environment before invoking this command. Python post-processing uses Python `3.10+`, NumPy and SciPy; its dependencies are listed in [requirements.txt](../requirements.txt). Molecular preparation does not need to share its environment with figure generation.

Native R calculations use R `4.5.2`, Matrix `1.7.4`, `jsonlite`, `data.table`, `plyr` and the dependencies of the selected method. PRODIGY additionally uses DESeq2 `1.50.2`, igraph `2.3.3`, PCSF `0.99.1`, MASS `7.3.65` and mixtools `2.0.0.1`. Oesophageal count preparation uses edgeR `4.8.2`. Supply `rscript`, the required library paths and the package `versions` map in native job configurations. Precision configurations use `library_paths`; oesophageal preparation uses `r_libraries`.

The supplied R entry points set single-thread numerical-library settings. Keep the same numerical environment across the graph and evidence states of a comparison.

## External method sources

Obtain these distributions separately, retaining their licences and notices:

- [TARGET-SL](https://github.com/RhysGillman/TARGET-SL/tree/ec057f299b3b7017e743660b3a96bcaab9e7fb4f), commit `ec057f299b3b7017e743660b3a96bcaab9e7fb4f`, GPL-3.0. PRODIGY uses `data/PRODIGY_pathways.rds` and eight files from `scripts/PRODIGY/`, in this order: `PRODIGY.R`, `run_single_PCSF.R`, `get_pathway_network.R`, `get_enrichment_pvalue.R`, `get_enriched_pathways.R`, `get_diff_expressed_genes.R`, `get_DEGs.R`, `analyze_PRODIGY_results.R`.
- DawnRank input normalisation uses `scripts/DawnRank/DawnNormalize.R` from that TARGET-SL distribution. The [DawnRank source](https://github.com/MartinFXP/DawnRank/tree/e34fefd2e0d52ff6fa8b3384cbf6d24eac4766fd) identifies the upstream method and its GPL licence.
- [PCSF](https://github.com/IOR-Bioinformatics/PCSF/tree/4e5f2707fdd096b760e5bb4ca60938a3c9836d78), commit `4e5f2707fdd096b760e5bb4ca60938a3c9836d78`, supplies the `0.99.1` package, solver and associated licence.

Source-data access and redistribution conditions apply separately from the CRISP code licence. The oesophageal source specification provides the public molecular-data URLs and source identities.

## Input structures

The colorectal population contains 85 models from 83 donors. The oesophageal population contains 59 models from 58 donors. Model-level contributions are averaged within donor before equal-donor averaging. PRODIGY uses all five master seeds: `104729`, `130363`, `155921`, `196613`, `228017`. Deterministic seed aliases describe the same allocation.

CSV matrices have unique ordered gene rows and sample columns, with the gene identifier in the first column. Mutation matrices are binary. Native retained-network TSVs have columns `gene1`, `gene2`, `combined_score`; scores are integers satisfying `700 < score <= 1000`, and undirected non-self edges are unique.

Allocation JSON contains a `models` roster with `model_id`, `sample_id`, `donor_id` and `query`, together with allocation records carrying `model_id`, `cell`, `budget` and `weights`. The four states are `G0H0`, `G1H0`, `G0H1`, `G1H1`; the corresponding underscored forms are accepted by the relevant readers. Weights are gene-to-rational mappings or unique gene/value pairs. PRODIGY records also identify the seed and complete emitted menu. Cutoff ties share allocation fractionally. Underfilled menus retain vacancies, and degree ranking includes eligible zero-degree genes.

Functional labels use `calls` records with `sample_id`, `gene`, `state` and `value`. `MEASURED` requires integer `0` or `1`. Unknown states retain null values. Calls shared across allocations must have the same biological identity and value. Diagnostics additionally require `source_binding.source_id`.

The colorectal category file supplies `sets.organoid`, `sets.cell_line` and `roles[gene].cosmic_tokens`. The study sets contain 751 organoid-core and 1,121 cell-line-core symbols. Oesophageal scope records provide `categories.core` and `categories.role` for the complete query union. Core-fitness and cancer-role partitions are evaluated separately.

## Molecular preparation and native rankings

```bash
python crisp.py native prepare-molecular molecular.json
python crisp.py native prepare common_domain.json
python crisp.py native compute native_job.json
python crisp.py native export-native export.json
```

Create these job-specific configurations with the full inputs and roster:

| Command | Required configuration fields |
| --- | --- |
| `prepare-molecular` | `contract`, `source_root`, `observed_genes`, `mutation`, `mutation_sample_map`, `network`, `expected_models`, `output` |
| `prepare` | Two ordered paths each for `counts` and `mutations`; `models`, `expected_models`, `expected_genes`, `output`; optional inherited `evidence` records |
| `compute` | `method`, `mode`, `models`, `expected_models`, `context`, `seeds`, `rscript`, `versions`, `output`, and method-specific fields below |
| `export-native` | `operation` and `output`, with the operation-specific inputs below |

The source-column and transformation contract for molecular preparation is implemented in [molecular.py](../src/crisp/core/molecular.py). Common-domain preparation retains the ordered shared genes and the binary64 values of inherited evidence.

For **PRODIGY**, provide `original_sources`, `counts`, `mutations`, `pathways`, `sample_origins`, `alpha`, `cells` and `tasks`. Sample origins use the literal values `normal` and `tumor`. Each task identifies `cell`, `sample`, `master_seed`, `patient_seed`, `network_id` and `network`. Tasks cover every model, cell and seed exactly once. The primary `alpha` is `0.05`. Four-cell jobs use the order `G0_H0`, `G1_H0`, `G0_H1`, `G1_H1`; each task supplies inherited `explicit_evidence`. Its available object contains `schema="cdf_explicit_evidence_v1"`, ordered genes and little-endian `binary64_le` values. Keep the corresponding integer patient seeds.

For **DawnRank**, native-mode inputs are `tumour`, `normal`, `mutation`, `network`, `network_id`, `normalisation_names`, `normalise_source` and `expected_normals`. The name CSV has `role`, `canonical` and `original` columns. Four-cell inputs instead provide `common_genes`, `mutation`, `edges_G0`, `edges_G1`, `inherited_H0`, `inherited_H1`, `native_genes_H0`, `native_genes_H1` and `expected_genes`. Set `mu=3` and explicit `epsilon` and `max_iterations`.

Native export operations are:

- `dawn_cells`: `directory` identifies a complete four-cell DawnRank output.
- `prodigy_carriers`: provide `directory`, `scope`, `cell_ids` and `source_id`. Each cell-identity row contains `model_id`, `graph`, `evidence` and `signature_id`.
- `native_ranks`: provide a `jobs` list, with `directory`, `pool`, `release` and `source_id` for each native output.

The `analyse` command accepts operation configurations for `prodigy_freeze`, `prodigy_evaluate`, `dawn_freeze`, `dawn_evaluate`, `factorial_decision_freeze`, `factorial_decision_evaluate`, `cross_method`, `bridge_freeze` and `bridge_evaluate`. Their argument schemas follow the functions registered in [analysis/run.py](../src/crisp/analysis/run.py). Allocation definitions precede functional-label lookup.

## Numerical precision and component comparisons

```bash
python crisp.py native precision ../CRISP_Base_Data/config/precision_crc_first.json
python crisp.py native verify-precision ../CRISP_Base_Data/config/precision_crc_first_check.json
python crisp.py native precision ../CRISP_Base_Data/config/precision_components.json
python crisp.py native mechanism ../CRISP_Base_Data/config/base_mechanism_freeze.json
python crisp.py native mechanism ../CRISP_Base_Data/config/base_mechanism_evaluate.json
python crisp.py native mechanism ../CRISP_Base_Data/config/base_mechanism_score.json
```

The first precision configuration reads the prepared RDS from
`Supplementary_Data/Allocation_and_Graph_Controls/Graph_Controls/runtime/prepared_input.rds`.
It contains ordered `genes`, model sample names, a named binary `mutation` matrix,
`graphs$G0` and `graphs$G1` with `transition` and `damping`, and `evidence$H0` and
`evidence$H1` with values and evaluability. Graph `degree` or `adjacency` supports
degree-allocation export. The roster supplies model, sample and donor identifiers
in the same order. The component precision configuration consumes the preceding
precision output, including its prepared R objects.

The strict setting is `epsilon=1e-12` with a 10,000-iteration ceiling. The `1e-4` comparison uses a ceiling of 100; intermediate `1e-6` and `1e-8` settings use 10,000. Stopping thresholds and cutoff-membership bounds are separate quantities. Score arrays use little-endian binary64, column-major order, with `genes.txt` and `models.txt` axes.

Component states `Pa_Db_Hh` cross propagation, damping/restart and evidence. `reuse_directory` identifies matching strict diagonal states. The mechanism calculations retain all 22 specified contrasts, shared model–gene outcomes, unknown-label bounds and donor weights. Score-reference calculations preserve signed reference scores and exact cutoff ties.

The following filenames refer to the software templates in `CRISP/config`.
The corresponding data-specific files and their dependencies are listed in
`CRISP_Base_Data/CONFIGURATIONS.json`.

| Configuration | Command | Purpose |
| --- | --- | --- |
| `precision_first.json` | `precision` | First colorectal transition |
| `precision_second.json` | `precision` | Second colorectal transition |
| `precision_components.json` | `precision` | Propagation–damping–evidence cube |
| `precision_verify.json` | `verify-precision` | Fixed-point and cutoff-bound calculation |
| `mechanism_freeze.json` | `mechanism` | Allocation and shared-label requests |
| `mechanism_evaluate.json` | `mechanism` | Complete component contrasts |
| `mechanism_score.json` | `mechanism` | Score-reference decomposition |

## Oesophageal calculations

The [source specification](../config/oesophageal_sources.json) contains the 59-model roster, 403 GTEx gastro-oesophageal-junction references, 118 VCFs, source URLs and source file names. The [mutation rule](../config/oesophageal_mutation_rule.json) defines qualifying calls and exclusions. These two files are inputs to the configuration chain, rather than standalone commands.

```bash
python crisp.py native oesophageal ../CRISP_Base_Data/config/oesophageal_inputs.json
python crisp.py native oesophageal ../CRISP_Base_Data/config/oesophageal_prepare.json
python crisp.py native precision ../CRISP_Base_Data/config/oesophageal_compute.json
python crisp.py native oesophageal ../CRISP_Base_Data/config/oesophageal_allocations.json
python crisp.py native oesophageal ../CRISP_Base_Data/config/oesophageal_scope.json
python crisp.py native oesophageal ../CRISP_Base_Data/config/oesophageal_freeze.json
python crisp.py native oesophageal ../CRISP_Base_Data/config/oesophageal_evaluate.json
```

Obtain the source count files and VCFs before this sequence. Preparation uses all
462 tumour and normal count columns and normalises on each native gene domain.
It intersects the 7,399-gene colorectal query domain with oesophageal support to
obtain the 7,388-gene comparison domain. Model identities, exclusions and category
assignments remain fixed across states. The allocation stage links numerical
scores to the model roster. Scope construction defines the functional source and
categories, and coefficient construction supplies the shared model-gene requests
for functional evaluation.

| Configuration | Role |
| --- | --- |
| `oesophageal_sources.json` | Source URLs, source paths and complete biological roster |
| `oesophageal_mutation_rule.json` | Mutation eligibility and exclusions |
| `oesophageal_inputs.json` | Tumour/normal counts, gene mapping and mutation matrix |
| `oesophageal_prepare.json` | Native normalisation, graph states and common domain |
| `oesophageal_compute.json` | Strict four-cell DawnRank calculation through `precision` |
| `oesophageal_allocations.json` | Allocation and degree records |
| `oesophageal_scope.json` | Functional-source and category definitions |
| `oesophageal_freeze.json` | Shared coefficient and label-request construction |
| `oesophageal_evaluate.json` | Functional contrasts and donor-composition summaries |

## Statistics

For supplied allocations, use the `base_statistics_*` configurations in
`CRISP_Base_Data/config`. For example:

```bash
python crisp.py native statistics ../CRISP_Base_Data/config/base_statistics_paired_crc.json
python crisp.py native statistics ../CRISP_Base_Data/config/base_statistics_endpoint_crc.json
```

The software templates below describe each operation:

| Configuration | Purpose |
| --- | --- |
| `statistics_paired_crc.json` | Joint and degree-relative updates, seed contributions |
| `statistics_draws_crc.json` | Ordered donors and 10,000 paired donor draws |
| `statistics_endpoint_crc.json` | Eight-coordinate colorectal update family |
| `statistics_categories_crc.json` | Colorectal category allocations and contributions |
| `statistics_categories_esca.json` | Oesophageal category allocations and contributions |
| `statistics_category_sensitivity_crc.json` | Colorectal quantity/fraction family |
| `statistics_category_sensitivity_esca.json` | Oesophageal quantity/fraction family |
| `statistics_signed_accounting.json` | Signed allocation and donor masses |

Colorectal statistics use the graph-wide `degree_allocations.json` from precision output. Colorectal draws consist of `donor_order.json` and `draw_indices.json`; oesophageal functional evaluation supplies donor order in `composition.json` and draws in `DONOR_DRAWS.npy`. Each contains 10,000 draws, aligned to lexicographic donor identifiers, with 83 or 58 sampled donors per draw.

Category ratios use pooled positive/returned mass and are recalculated in each draw. Zero-allocation donors stay in the cohort denominator. Unknown-label identification bounds and donor-composition envelopes are reported separately. Update and category families each use their specified calibration and endpoint set.

## Rank and candidate diagnostics

Use `python crisp.py native diagnostics` with the appropriate
`CRISP_Base_Data/config/base_diagnostics_*.json` configuration. Normalise the
combined allocation ledger before rank or influence calculations; construct
native menus before degree-menu comparisons, and menu records before the
opportunity comparison. The software templates below describe each operation:

| Configuration | Purpose |
| --- | --- |
| `diagnostics_prodigy.json` | Native PRODIGY carrier-to-allocation conversion |
| `diagnostics_normalise.json` | Combined model and allocation ledger |
| `diagnostics_rank.json` | Budgets, seed variation, turnover and label sensitivity |
| `diagnostics_influence.json` | Gene and donor contributions |
| `diagnostics_native_menus.json` | Complete emitted native menus |
| `diagnostics_degree_menu.json` | Emitted-menu and full-query degree comparisons |
| `diagnostics_menu_records.json` | Paired candidate-menu records |
| `diagnostics_opportunity.json` | Menu-conditioned candidate overlap |

Retain all declared budgets `1`, `5`, `10`, `20` and the complete emitted menus. The normalisation and rank calculations consume the corresponding four-budget records. Label files must cover every allocated model–gene identity, including degree allocations and all requested budgets. Full schemas are implemented in [statistics/workflow.py](../src/crisp/statistics/workflow.py) and [diagnostics/workflow.py](../src/crisp/diagnostics/workflow.py).

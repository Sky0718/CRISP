# Data and dependencies

The companion data deposit contains two input collections. `Supplementary_Data`
provides the model-gene records and summaries used by the `reproduce` and
`analysis` commands. `CRISP_Base_Data` provides native scores, allocations,
prepared inputs and dataset-specific configurations for the `native` commands.
Analysis code is distributed in this repository.

Extract the 16 Supplementary Data ZIP parts and `CRISP_Base_Data.zip` into a
common parent directory. Use this layout for the supplied configuration paths:

```text
CRISP/                 code checkout
CRISP_Base_Data/        base data and dataset-specific configurations
Supplementary_Data/    supplementary records and figure sources
```

Pass `Supplementary_Data` itself, rather than one ZIP part or a component
subdirectory, to `--data`. Commands accept either relative or absolute input paths.

## Supplementary inputs

The supplementary collection contains model-gene ledgers, rational allocation
weights, functional-screen outcomes, graph-control inputs, continuous depletion
records and numerical figure sources. `MANIFEST.json` records file checksums.
`FIGURE_SOURCE_MAP.json` maps the five main and four supplementary figures to
their numerical records. The deposit-level `FILE_INDEX.csv` identifies each
file's ZIP archive and checksum.

| Record group | Role |
| --- | --- |
| `Union_Sensitivity_Analysis` | Cancer-listed allocations, missing-outcome bounds, donor/role decomposition, budget and graph-control comparisons |
| `Union_Validation` | Continuous depletion summaries, measured-label sensitivity and four cancer-role groups |
| `Reader_Data/Data_10/Union_Analysis` | Colorectal and oesophageal union summaries and donor-composition envelopes |
| `Allocation_and_Graph_Controls` | Donor decomposition, candidate context and graph-control inputs |
| `Colorectal_Functional_Analysis` | Crossed graph/input category allocations and continuous readouts |
| `Variant_and_Network_Context` | Exact-allele annotation, relative connectivity and donor-composition figure inputs |
| `Main_Figure_Reproduction` | Fixed plotting inputs for Figures 1, 3 and 5 |
| `Figure_Source_Data` | Numerical carriers for every main and supplementary figure |

The primary cohort comprises 85 colorectal models from 83 donors. The oesophageal comparison comprises 59 models from 58 donors in the same source study. Models are averaged within donor before equal donor weighting. Functional outcomes are associated with literal model-gene identities; unmeasured outcomes remain separate from measured zeros.

## Base inputs and external sources

The configuration index `CRISP_Base_Data/CONFIGURATIONS.json` identifies
the command for each file under `CRISP_Base_Data/config`. These configurations
connect supplied records and dependent calculations. The files under this code
repository's `config` directory describe the software interfaces and shared
rules. Use the data collection's configurations for its corresponding analyses,
or adapt a configuration's paths to a new dataset.

Recomputing rankings from original molecular and network inputs requires the specified STRING physical-network releases, tumour and reference expression data, variant records and upstream method implementations. Functional-screen data supply the evaluation outcomes. Source identifiers and input contracts are given in the manuscript, Supplementary Information, supplied data records and [Native calculations](native.md). Obtain provider-controlled inputs through their original access routes.

## Paths and output directories

CLI data paths may be relative to the working directory. Paths inside native JSON configurations are relative to the configuration file. Executable names such as `Rscript` may resolve through `PATH`.

Use a fresh output directory outside the code repository. For supplementary
commands, it must also be outside the input-data tree. Native configurations
identify separate calculation destinations; check their `output` field before
invocation and choose a new destination if it already exists. If a configuration
is copied elsewhere, update its relative paths to preserve their targets.

When relocating inputs, retain complete model populations, candidate queries,
graph states and algorithm seeds. A missing functional outcome stays unknown,
separate from a measured zero.

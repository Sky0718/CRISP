# CRISP

CRISP (Cancer-gene Ranking Interpretation with Screen Profiles) links changes in cancer-gene rankings to functional evidence in the same tumour models. It follows which model-gene pairs gain or lose priority after a network update, how much of the selected allocation has measured fitness support, and how changes within donors and in donor composition contribute to pooled summaries.

This repository accompanies *Cancer-gene enrichment and fitness support diverge after a network update*. The study pairs network-based rankings with depletion screens in 85 colorectal tumour-derived organoid models from 83 donors, with a same-study comparison in 59 oesophageal models from 58 donors. Its central question is whether selecting more cancer-listed genes also concentrates model-matched fitness support.

## What the code provides

- Paired allocation summaries that keep measured-positive, measured-zero and unknown functional outcomes separate.
- Cancer-role and fitness-category comparisons, continuous depletion summaries, shared-label bounds and donor decompositions.
- Candidate-level variant context, graph comparisons and the numerical inputs for the manuscript figures.

Scientific definitions, interpretation and cohort-specific procedures are given in the manuscript and Supplementary Information. This repository describes how to access the corresponding calculations and figures.

## Installation

Use Python 3.10 or later in a dedicated environment. From the repository root, install the analysis dependencies and the CRISP package:

```bash
python -m pip install -r requirements.txt
python -m pip install --no-deps -e .
```

For figure generation, also install `requirements_figures.txt`. Native ranking calculations use a separate environment with `requirements_native.txt` and the R dependencies described in [Native calculations](docs/native.md).

## Prepare the data

The companion data deposit contains 17 ZIP archives: `CRISP_Base_Data.zip` and 16 Supplementary Data archives. Extract them into one directory, preserving their internal paths. The four Data 05 parts are independent ZIP archives and are extracted in the same way as the other files.

For both the analysis commands and the supplied native configurations, use this arrangement:

```text
workspace/
  CRISP/                 # this code checkout
  CRISP_Base_Data/        # base scores, allocations and native configurations
  Supplementary_Data/    # analysis records, tables and figure inputs
```

The commands below assume the working directory is `CRISP`. Research data and new results remain outside the code checkout. [Data and dependencies](docs/data.md) explains the input groups and the provider files needed for source-level molecular reconstruction.

## Start with the principal summaries

```bash
python crisp.py reproduce --data ../Supplementary_Data --output ../CRISP_results
python crisp.py reproduce --data ../Supplementary_Data --output ../CRISP_continuous --union-readouts
```

The first command derives cancer-listed allocation summaries, positive fractions, shared-label bounds and donor/role decompositions from the supplied records. The second includes continuous depletion summaries and measured-label sensitivity. Each command takes a new output directory and leaves its input records unchanged.

Use `python crisp.py analysis --help` to select a component. The [Analysis guide](docs/analyses.md) maps each component to its scientific question and input records. [Figure generation](docs/figures.md) describes the plotting dependencies and maps output files to the manuscript figures.

## Native ranking calculations

```bash
python crisp.py native statistics ../CRISP_Base_Data/config/base_statistics_paired_crc.json
python crisp.py native precision ../CRISP_Base_Data/config/precision_crc_first.json
```

The statistical example reads the supplied allocation records. The precision example requires the specified R environment and prepared graph input. Native JSON paths are relative to the configuration file, so the sibling-directory arrangement above matters. [Native calculations](docs/native.md) distinguishes analyses of supplied scores from molecular preparation and ranking calculations that require additional provider inputs.

## Tested hardware

The project has completed runs on the following systems:

| System | Processor | Graphics | Memory |
|---|---|---|---:|
| Ubuntu Server 24.04.4 LTS | Intel Core i9-14900K | NVIDIA GeForce GTX 750 Ti | 128 GB DDR5-4000 |
| Ubuntu Server 24.04.4 LTS | Intel Core i9-14900K | NVIDIA GeForce GTX 1050 Ti | 128 GB DDR5-4000 |
| Ubuntu Server 24.04.4 LTS | Intel Core i9-14900K | NVIDIA GeForce GTX 1050 Ti | 64 GB DDR5-4800 |
| Ubuntu Server 24.04.4 LTS | 2 × Intel Xeon Gold 6258R | Integrated graphics | 256 GB DDR4-2666 |
| Ubuntu Server 24.04.4 LTS | 2 × Intel Xeon Platinum 8269CY | Integrated graphics | 256 GB DDR4-2666 |
| Ubuntu Server 24.04.4 LTS | 2 × Intel Xeon Platinum 8272CL | Integrated graphics | 256 GB DDR4-2666 |
| Ubuntu Server 24.04.4 LTS | 2 × Intel Xeon Platinum 8375C | Integrated graphics | 512 GB DDR4-2666 |
| Windows 11 25H2 | Intel Core Ultra 9 285H | NVIDIA GeForce RTX 5070 | 32 GB LPDDR5X-7467 |

## Server cluster involved in the study

The project has completed runs on the server cluster involving the following servers:

| ID | System | Processor | Graphics | Memory |
|---|---|---|---|---:|
| 1 | Ubuntu Server 24.04.4 LTS | Intel Core i9-14900K | NVIDIA GeForce GTX 1050 Ti | 64 GB DDR5-4800 |
| 2 | Ubuntu Server 24.04.4 LTS | Intel Core i9-14900K | NVIDIA GeForce GTX 1050 Ti | 64 GB DDR5-4800 |
| 3 | Ubuntu Server 24.04.4 LTS | Intel Core i9-14900K | NVIDIA GeForce GTX 1050 Ti | 64 GB DDR5-4800 |
| 4 | Ubuntu Server 24.04.4 LTS | Intel Core i9-14900K | NVIDIA GeForce GTX 750 Ti | 128 GB DDR5-4000 |
| 5 | Ubuntu Server 24.04.4 LTS | Intel Core i9-14900K | NVIDIA GeForce GTX 1050 Ti | 128 GB DDR5-4000 |
| 6 | Ubuntu Server 24.04.4 LTS | Intel Core i9-14900K | NVIDIA GeForce GTX 1050 Ti | 128 GB DDR5-4000 |
| 7 | Ubuntu Server 24.04.4 LTS | Intel Core i9-14900K | NVIDIA GeForce GTX 1050 Ti | 128 GB DDR5-4000 |
| 8 | Ubuntu Server 24.04.4 LTS | Intel Core i9-14900K | NVIDIA GeForce GTX 1050 Ti | 128 GB DDR5-4000 |
| 9 | Ubuntu Server 24.04.4 LTS | Intel Core i9-14900K | NVIDIA GeForce GTX 1050 Ti | 128 GB DDR5-4000 |
| 10 | Ubuntu Server 24.04.4 LTS | 2 × Intel Xeon Gold 6258R | Integrated graphics | 256 GB DDR4-2666 |
| 11 | Ubuntu Server 24.04.4 LTS | 2 × Intel Xeon Gold 6258R | Integrated graphics | 256 GB DDR4-2666 |
| 12 | Ubuntu Server 24.04.4 LTS | 2 × Intel Xeon Platinum 8269CY | Integrated graphics | 256 GB DDR4-2666 |
| 13 | Ubuntu Server 24.04.4 LTS | 2 × Intel Xeon Platinum 8272CL | Integrated graphics | 256 GB DDR4-2666 |
| 14 | Ubuntu Server 24.04.4 LTS | 2 × Intel Xeon Platinum 8272CL | Integrated graphics | 256 GB DDR4-2666 |
| 15 | Ubuntu Server 24.04.4 LTS | 2 × Intel Xeon Platinum 8375C | Integrated graphics | 512 GB DDR4-2666 |
| 16 | Ubuntu Server 24.04.4 LTS | 2 × Intel Xeon Platinum 8375C | Integrated graphics | 512 GB DDR4-2666 |
| 17 | Ubuntu Server 24.04.4 LTS | 2 × Intel Xeon Platinum 8375C | Integrated graphics | 512 GB DDR4-2666 |
| 18 | Microsoft Windows 11 25H2 | Intel Core Ultra 9 285H | NVIDIA GeForce RTX 5070 | 32 GB LPDDR5X-7467 |

## Licence and citation

Original project code is distributed under the [Apache License 2.0](LICENSE). The bundled DawnRank-derived R routines are covered by [GNU GPL v3](licenses/GPL-3.0.txt), as detailed in [DawnRank attribution](licenses/DawnRank.md). External methods and data retain their own terms; see [Third-party components](docs/third-party.md) and [RBO attribution](licenses/RBO.txt).

Citation metadata are provided in [CITATION.cff](CITATION.cff).

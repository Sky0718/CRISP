# Cancer-listed union analysis

This analysis measures how the allocation and functional support of
cancer-listed genes change between the old and joint-update states. It uses
fixed model-gene allocations and is exploratory, with endpoints defined after
exposure to the allocation results.

## Populations and inputs

CRC includes 85 models from 83 donors and three policies: PRODIGY, DawnRank
and full-query degree. All five PRODIGY seeds are averaged within model.
ESCA includes 59 models from 58 donors, with DawnRank and degree. Models have
equal weight within donor and donors have equal weight. Method contrasts
compare G0H0 with G1H1; degree contrasts compare G0 with G1.

The cancer-listed union comprises ONCOGENE, TSG, DUAL_ONCOGENE_TSG and
OTHER_OR_UNSPECIFIED. These roles are disjoint. Selected BAGEL2 calls
distinguish unknown from measured zero; each unknown model-gene label is
shared across all states.

Paths in this protocol are relative to the data collection's
`Reader_Data/Data_10/Union_Analysis` directory. Primary records under
`inputs/archive` are:
- results/crc_precision/degree_allocations.json
- results/crc_prodigy_allocations.json
- results/crc_precision/allocations.json
- data/crc/selected_labels.json
- data/crc/categories.json
- results/esca_allocations/allocations.json
- results/esca_evaluation/selected_labels.json
- results/esca_scope/scope.json

The loaders are provided by crisp.statistics.inputs. Model and donor totals
reconcile with inputs/reference/Category_Sensitivity/model_category_levels.csv
and donor_category_levels.csv. CRC means and bounds also reconcile with
inputs/reference/Cancer_Allocation_Overview/Figure_Overview_Source.csv and
the inputs/reference/functional_content ledger.

## Endpoints and families

State endpoints are returned mass R, measured-positive mass H, measured-zero
mass Z, unknown mass U and pooled H/R. Masses are allocation per ten offered
slots. H/R includes unknown allocation in R. Changes comprise delta_R,
delta_H, delta_H_over_R and the lower/upper sharp shared-unknown bounds for
H change. Relative update-minus-degree bounds combine method and degree
coefficients before optimisation.

Each context has its own simultaneous family: five coordinates per policy
plus two relative-bound coordinates per ranking method. CRC therefore has
19 coordinates and ESCA has 12. Duplicate and constant coordinates remain
in their families, which are calibrated separately from other endpoint
families.

## Donor-composition sensitivity

All 10,000 paired whole-donor draws use the stored order in
inputs/reference/Category_Sensitivity/donor_index.csv and its corresponding
colorectal_donor_draws.csv and oesophageal_donor_draws.csv. The seeds are
20261004 for CRC and 20261005 for ESCA. Each draw samples D donors with
replacement, retaining nested models, seed averages and paired states.

Linear endpoint scales use the fixed original-donor standard error

sqrt(sum((x_d-mean(x))^2)/(D(D-1))).

For a pooled H/R difference, the donor influence is

(H_d1-p1 R_d1)/mean(R1) - (H_d0-p0 R_d0)/mean(R0).

Ratios are recomputed in every draw. Sharp bounds are also recomputed:
nonnegative resampling weights and identities distinct across donors make
weighted donor bounds the exact shared-label endpoints. Selected draws
receive an independent identity-level check.

The multiplier q is the linearly interpolated 95th percentile of the
maximum absolute fixed-scale deviation over the complete context family.
Coordinate envelopes are theta +/- q*SE. The outer bound envelope is
[lower_estimate-q*SE_lower, upper_estimate+q*SE_upper].
A constant coordinate contributes zero only when every draw equals its
original value. Undefined denominators or nonconstant zero-scale endpoints
disable the affected family while retaining all coordinates and draws.
Monte Carlo quantile uncertainty uses binomial order statistics.

These envelopes describe sensitivity to the composition of the recorded
cohort. Assay-error uncertainty and between-study variation are separate
quantities. No population sampling confidence intervals or P values are
assigned to these envelopes.

## Influence and numerical checks

Outputs include every donor contribution and sign, all leave-one-donor-out
estimates, bounds and pooled fractions, extremes with donor identities, and
contribution concentration. Checks cover SHA-256 bindings, exact rational
identities, conservation, population counts and five-seed completeness.
Complete numerical source tables, draw endpoint matrices, source/output
manifests and numerical comparisons accompany the summaries.

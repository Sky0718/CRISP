# Cancer-listed union sensitivity

This analysis evaluates completed positive-fraction bounds, donor and role
contributions, and sensitivity to budget and graph construction. It uses fixed
model-gene allocations and source labels in an explanatory post-exposure
comparison. CRC comprises 85 models from 83 donors, with the
five-seed PRODIGY average, DawnRank and degree. ESCA comprises 59 models
from 58 donors, with DawnRank and degree. The four mutually exclusive
cancer roles, unknown-call subtypes, fractional allocations, equal
model-within-donor weights and equal donor weights define all endpoints.

## Completed positive-fraction bounds

The completed positive-fraction change has fixed returned-mass denominators
R0 and R1. Each model-gene identity has coefficient

population_weight*(new/R1-old/R0).

Shared identities are combined before minimising or maximising unknown
binary calls. Explicit endpoint label assignments attain the sharp bounds.
The known component equals H1/R1-H0/R0, and the sharp interval is contained
in the conservative independent-arm interval. These are fixed-cohort
identification bounds, separate from the 19-coordinate CRC and 12-coordinate
ESCA donor-composition families for observed-fraction endpoints.

## Donor and role decomposition

The observed union H/R change is decomposed by donor and separately by
exclusive cancer role using the midpoint-centred identity. The components
are within-common-support change, common-support reweighting, entry and
exit. Zero-allocation units remain represented with undefined fractions.
Exact rational summation reconciles all components with the direct pooled
difference. Complete role levels and changes accompany the descriptive
decomposition.

## Budgets and graph controls

Budget sensitivity covers k=1, 5, 10 and 20 using the supplied allocations
or rankings. Graph sensitivity covers all 27 specified comparisons: 18
fixed-input DawnRank contrasts and nine degree contrasts across the original,
matched-edge-count and alternative-threshold graph pairs. Each comparison
retains the complete CRC population and contributes R, H, Z, U, H/R and sharp
completed-fraction change bounds. Missing required inputs stop the calculation.
These finite comparisons introduce no additional resampling family.

## Numerical checks and outputs

Checks reconcile k=10 levels and graph-control reference states with
their source records, enumerate small binary optimisation fixtures,
evaluate attaining labels on the data, and establish allocation
conservation and every decomposition identity. Outputs include exact
rational and decimal tables, input hashes and complete numerical
comparisons, including opposing changes and intervals spanning zero.

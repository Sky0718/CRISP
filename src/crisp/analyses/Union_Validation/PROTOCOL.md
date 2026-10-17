# Cancer-listed union functional and label-sensitivity analysis

The analysis characterises continuous depletion within the selected
cancer-listed allocation, the minimum measured-label changes needed to alter
its binary conclusion, and the contribution of each cancer role. It is an
exploratory post-exposure analysis. Complete outputs retain reversals,
incomplete coverage and envelopes spanning zero.

## Population and allocations

CRC comprises 85 models nested in 83 donors. PRODIGY allocations average
seeds 104729, 130363, 155921, 196613 and 228017 within model. DawnRank and
degree are deterministic. Donors receive weight 1/83 and models within a
donor receive equal weight. The cancer-listed union contains ONCOGENE, TSG,
DUAL_ONCOGENE_TSG and OTHER_OR_UNSPECIFIED. The old and joint-update states
share the ten-slot offered budget, model queries, annotations and source
calls. Their distinct graph/input states and resulting allocations are
supplied records and remain fixed throughout this analysis.

## Continuous functional readout

Every selected CRC union model-gene identity links to its original Table 6
continuous identity record, with model/sample, library and raw decimal
bindings. Depletion is D=-LFC using source-scaled LFC directly. For each
state and each of PRODIGY, DawnRank and degree, endpoints are returned
allocation R, finite-LFC allocation A, missing allocation R-A, coverage A/R,
observed depletion sum S=sum(alpha*w*D), and pooled measured depletion S/A.
Contrasts are new minus old. S covers measured allocation; complete-return
depletion is identified when A=R. Unmeasured continuous outcomes remain
missing.

The paired resampling design draws 83 whole donors with replacement 10,000
times using NumPy Generator(PCG64), seed 470016. Each draw retains nested
models, averaged seeds and both states, and recomputes sums and ratios.
The 15-coordinate family comprises three policies x changes in R, A, A/R,
S and S/A. Each coordinate's bootstrap sample SD (ddof=1) standardises
draw-minus-original differences. The common multiplier is the linearly
interpolated 95th percentile of the maximum absolute standardised deviation
over all coordinates. Envelopes are estimate +/- multiplier*SD. Constant
coordinates contribute zero standardised deviation. Complete donor draws,
endpoints, family membership and seed accompany model-gene, model and
donor records. The envelopes measure within-cohort composition sensitivity.

## Exact measured-label flips

The five context-policy contrasts are CRC PRODIGY, DawnRank and degree,
and ESCA DawnRank and degree. ESCA contains 59 models from 58 donors;
both contexts retain their complete populations and donor weights.

For each identity, shared normalised allocation coefficients combine as

c_i=alpha_i*(w1_i/R1-w0_i/R0).

R0 and R1 remain fixed. Unknown labels attain the joint upper bound with
y_i=1 where c_i>0 and y_i=0 otherwise. An admissible perturbation flips one
measured binary identity consistently across seeds and states. Its signed
upper-bound change is

g_i=c_i*(1-2*y_i).

Positive gains are sorted in descending exact rational order, with sample
ID and gene as deterministic tie-breakers. The largest k gains give the
maximum increase attainable with at most k measured-label flips.
Writing U for the initial completed-fraction upper bound in this section,
the minimal count is the smallest k satisfying

U+sum(g_(1:k))>=0.

The count is zero if U>=0 and unreachable if all admissible gains are
insufficient. Outputs include the complete gain curve, all gains, one
deterministic minimal witness, exact before/after bounds, role/donor/
label-direction composition, and the witness fraction among all measured
selected identities. The CRC conservative count lower bounds are 19 for
PRODIGY and 11 for DawnRank; the exact minimum counts attain these bounds.
This endpoint is a finite adversarial
perturbation count, not an estimated assay error rate or error probability.

## Four-role composition

All four roles, both states and all five context-policy contrasts contribute
R, measured-positive H, measured-zero Z, unknown U, H/R and union allocation
share. Exact midpoint within-role and role-reweighting contributions sum to
the observed union fraction change. Here U denotes unknown allocation,
distinct from the upper-bound notation above. A measured-zero source
depletion call describes the assay endpoint, not a gene's oncogenic role.
The role decomposition is descriptive accounting.

## Numerical checks and outputs

Binary coefficients, flips, role levels and decomposition use exact rational
arithmetic; LFC source values retain their decimal-to-rational binding.
Independent ledger reconstruction checks before/after fractions for each
minimal witness. Maximal k-1 gain establishes minimality, supported by
small enumerated fixtures.

Continuous checks independently aggregate source allocations and donor
records through direct matrix arithmetic, reconcile model/donor mass and
old/new states, and recompute stored draw endpoints and envelopes. Direct
input hashes bind each result to its source. Outputs use a separate
directory from their source records and contain all three analysis
components and their numerical comparisons.

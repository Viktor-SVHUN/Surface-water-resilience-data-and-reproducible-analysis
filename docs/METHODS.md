# Operational methods and interpretation

## 1. Source records and annual aggregation

The monitoring workbook contains annual summaries for monitoring sites: an annual mean, the reported contributing sample count and water-body/site identifiers. The baseline period is 2000–2022. It is not a sequence of independently dated within-year samples.

Records outside the selected physicochemical scope, without valid identifiers/values, outside the time window or with incompatible units are excluded and counted. Mass concentrations are harmonised to g/L for computation; selected figures display mg/L. Turbidity in NTU is not pooled with mass-based turbidity records.

Copies identical in the retained site/parameter/year/value/unit/sample-count fields are consolidated. Remaining versions for the same site/parameter/year/unit are averaged symmetrically, with conflicting values exported. A water body's annual value is the equal-weight mean of its available site means. Sample-count-weighted values are an alternative scenario; contributing counts are approximate coverage annotations, not audited unique sample totals.

## 2. Registry and recipient graph

VOR is the stable water-body/code identifier. KTJ identifies a monitoring site. Registry fields keep their source names; `VT_VOR_2015` is joined to monitoring `VOR`, and `BEF_VT_VOR` provides the designated recipient.

The river catalogue also contains lake-category (`WB_CAT=LW`) entries. Shared river/lake identifiers are consolidated, using the lake record where the explicit category and recipient agree. The final type is not determined by river-catalogue membership alone.

The computational graph's node set is the union of registered identifiers, observed monitoring codes and recipient endpoints. This is broader than the set of registered water bodies. Observed unregistered codes may therefore be present as isolated nodes. Flags `registered` and `observed` distinguish these populations; only registered River/Lake series with an applicable direction/threshold enter threshold-based local analysis. The graph is directed upstream-to-recipient. The release inventory reports graph and registry counts separately.

Network candidates are DCI, UPS, BC, DCO, DOWN, HC_IN, PRC and AC_alpha09, in that priority order. Outcome-blind pruning removes a descriptor when its absolute Spearman correlation with an already retained one reaches 0.90. Selection is applied to observed registered nodes. Additional diagnostic descriptors may be exported without entering this primary set. `metric_dictionary.csv` defines their algorithms and projections.

## 3. Local indicators and ranking

For value x and threshold T, the oriented deviation is d = s(x − T), with s = +1 when lower values are favourable and s = −1 when higher values are favourable. Dissolved oxygen and oxygen saturation use s = −1. Equality is compliant; the implementation uses a relative numerical tolerance of 10^−12 when classifying local annual values.

| Field | Definition | Unit |
| --- | --- | --- |
| `N_good`, `N_bad` | Counts of compliant/non-compliant observed annual values | Annual observations |
| `Amax_bad` | max(0, max d) | Parameter unit |
| `Amax_good` | max(0, max −d) | Parameter unit |
| `Amax_bad_relative`, `Amax_good_relative` | Corresponding magnitude divided by T | Dimensionless |
| `Lmax_bad_years`, `Lmax_good_years` | Longest adverse/compliant run | Years |
| `ARL_bad_years`, `ARL_good_years` | Mean run length for the specified state only | Years |
| `tR_first_observed_years` | First later compliant year minus the worst adverse year | Calendar years |

A run requires consecutive observed calendar years; an unobserved year splits it. An absent state has an unavailable average run length. If no failure occurs, recovery time is undefined (`no_adverse_year`). If a failure occurs but no subsequent compliant observation exists, recovery is right-censored. Transition bounds describe the observed annual transition; they do not rule out unobserved within-year returns.

The ranking score is the mean of `Lmax_bad_years / n_years` across parameters with at least eight observed years. A body needs at least five such parameters. Critical status is the upper tail at or above the empirical 90th percentile; ties use the minimum rank. This is not a statutory classification.

Sensitivity scenarios use direction-aware threshold shifts of ±5%, sample-count weighting, a common parameter set and omission of each calendar year. In the code, the common set requires a parameter's eligible-body count to reach 80% of the largest such count among parameters. Scenario agreement uses Spearman score correlation, Jaccard similarity of critical sets and top-five overlap. The all-scenario and year-omission Jaccard medians have different scopes and are reported separately. Omission envelopes are not sampling confidence intervals.

## 4. Horizontal visibility graphs

For an ordered annual series, z = s(x − mean x). A directed edge i→j, i<j, exists when every intermediate z[k] is strictly less than min(z[i],z[j]). Ties block visibility. Positive z means relatively adverse compared with the same body's mean, not necessarily regulatory non-compliance. Strictly increasing transformations leave the graph unchanged; a sign reversal need not do so.

Primary HVG series require eight years and a known orientation in a registered River/Lake body. Gap duration is not encoded as an edge length. The six primary metric families are DZO, LeM_U, Tra_U, GLE_D, AVPL_D and VAs_D; outcome-blind redundancy pruning follows the fixed priority order. D denotes directed metrics and U the undirected projection. The metric dictionary also identifies secondary descriptors, including explicit random/targeted vertex-deletion AUCs.

Sensitivity sampling selects up to 50 series per parameter with at least 12 years. Random and contiguous-block deletions remove 10%, 25% and 50% of values with 20 repetitions. The across-series correlation compares the baseline metric to the median perturbed metric. Additional outputs cover sample weighting, two-year bins and minimum-length association scenarios at 6, 8, 10 and 15 years. Shorter-series sensitivity computes the additional eligible series explicitly rather than filtering only the eight-year baseline table.

## 5. Associations and uncertainty

Each water body supplies one outcome per parameter. Descriptor and outcome ranks are separately residualised against the control matrix; their residual correlation is the controlled effect.

Controls include available catchment size, river length/lake area, depth, artificial/modified status, observed years, site count, type, altitude, geochemistry, slope and substrate. Numeric missingness is explicitly flagged, numeric values use within-analysis median replacement and ranks, categorical fields use dummy variables, and rank-deficient columns are removed by pivoted QR. The exported `controls` field records the actual fitted design; `typology_code` is also tested as a grouping in the separate profile analysis but is not an additional dummy block in this control-matrix function.

Tests require at least ten bodies and eight sub-basin clusters. Unknown sub-basin labels share a single cluster. Predictors with at most three distinct levels are inference-ineligible if any level has fewer than ten bodies or five supporting clusters. Effect-only results remain present with an explicit status.

The slope test uses a cluster sandwich variance with factor G/(G−1) × (n−1)/(n−k−1), where k is the number of columns in the control design, including the intercept. The reference distribution is t with G−1 degrees of freedom. Intervals use 299 paired sub-basin resamples of fixed residual ranks, conditional on the fitted ranks/controls; they are not fully refitted Spearman-bootstrap intervals. Dependence crossing sub-basin boundaries is not eliminated by this procedure.

BH correction is applied separately to local versus HVG and controlled versus unadjusted test families, pooling their parameter-specific tests. Statistical significance is q<0.05. An absolute coefficient of 0.3 is a pragmatic practical-interest threshold, not a universal scientific boundary.

## 6. Upstream-conditioned compliance and placebo

Direct predecessors in the recipient graph define upstream neighbours. Annual pairs require an observed downstream body and at least one observed upstream neighbour. Exposure means at least one observed upstream neighbour is non-compliant. TR is downstream compliance among exposed years; at least five exposed years are needed for the reported interval. The moving-block bootstrap uses three-year calendar blocks and 299 repetitions, retaining the calendar gap pattern.

The exported clean-upstream counterpart requires **complete upstream observation** with all direct neighbours compliant. The table separately records exposure coverage, complete-coverage exposed estimates and matched/exposed counts. Consequently, clean and exposed estimates may use different subsets of years. Their difference is descriptive, not a causal buffering estimate.

Mantel–Haenszel odds ratios compare downstream non-compliance with versus without observed upstream non-compliance, stratified by body. The code uses the Robins–Breslow–Greenland log-OR variance and a continuity-adjusted CMH statistic. The exact implementation is in `mantel_haenszel`.

For the placebo, each body is assigned the same number of pseudo-upstream bodies with no directed path to or from it. Candidates from the same sub-basin are preferred when enough are available. Assignment is seeded and fixed across years. The ratio of connected and placebo odds ratios is supplemented by an interval obtained by adding their log-OR variances and omitting covariance. This is an independence approximation; shared outcomes alone do not establish that the interval is conservative. It does not provide causal attribution or account for every temporal/spatial dependence mechanism. BH correction is applied across parameters.

## 7. Physicochemical profiles and exploratory display

Profiles use the temporal mean of s(x/T−1) for each body and parameter. Bodies need at least three observed parameters; parameters need observations in at least ten bodies. Retained missing cells are replaced by the parameter median, and each column is standardised to zero mean/unit variance.

PERMANOVA uses Euclidean centroid sums of squares with 999 label permutations. PERMDISP evaluates distances to group centroids with permutation-based comparison; dispersion differences complicate a centroid-only interpretation. Typological groups smaller than ten are excluded. Univariate River/Lake comparisons use non-imputed profile values, Mann–Whitney tests and BH adjustment.

t-SNE uses only the physicochemical feature matrix: two dimensions, perplexity 30, learning rate 200, 1000 iterations, PCA initialisation, Euclidean metric, Barnes–Hut method, angle 0.5 and seed 42. Type is added after embedding for colouring. Coordinates are exploratory and may vary with numerical-library versions; distances between displayed groups are not evidence of a typological effect.

## 8. Scope and reproducibility

These operational indicators concern annual physicochemical summaries, not individual disturbance/recovery events or complete ecological status. Thresholds are the supplied type-level operational values; legal equivalence or subtype-specific completeness is not asserted. No official-status input is bundled, so no official RBMP assessment benchmark is claimed.

The release preserves all four source input files byte-for-byte. A full run records configuration, dependency versions, input hashes, program hash and self-tests. The distributed results are accompanied by a validation report; numerical agreement does not validate all modelling assumptions. No external upload or data download is performed by the analysis program.

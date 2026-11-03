# Release validation

The distributed program was executed end-to-end on all four bundled source inputs using the full configuration recorded in `results/manifest.json`. The run completed without changing the source input bytes. Results in `results/` are from that execution, not reduced diagnostic settings.

## Checks performed

- All 24 built-in scientific self-tests passed, including independent ts2vg edge agreement on continuous and tied synthetic series.
- Direct-loop and running-maximum HVG edges agree on all eligible annual series; the disagreement table has no rows.
- 48 retained result tables were compared with the supplied full-run tables. Primary deterministic results and selected inferential metrics agree at relative tolerance 10^-8 and absolute tolerance 10^-10. `NUMERICAL_REGRESSION.csv` distinguishes agreement, secondary diagnostic differences and tables with no numeric columns.
- Input hashes match the source archive; `PROVENANCE.json` records the source and release identifiers.
- Output paths are portable; existing non-empty output directories are rejected.
- Packaging checks cover filename consistency, source compilation, table inventory and archive checksums.
- Six release-level unit tests passed, including distributed-ranking regression and working-directory-independent input discovery.
- All nine LaTeX table exports and the numerical macro file compile together without LaTeX warnings or errors in the validation document.

## Reproduced counts

| Quantity | Value |
| --- | ---: |
| Source monitoring workbook rows | 269,288 |
| Unique registered rivers after category reconciliation | 889 |
| Unique registered lakes | 189 |
| Computational graph nodes | 1,207 |
| Recipient arcs | 1,066 |
| Local body–parameter rows | 8,745 |
| HVG body–parameter series | 2,975 |
| Ranking-eligible bodies | 258 |
| Critical bodies | 26 |
| All-scenario median Jaccard similarity | 0.677419 |
| Single-year-omission median Jaccard similarity | 0.625000 |

The two Jaccard summaries use different scenario scopes; their difference is not a data discrepancy. The graph is broader than the registered population, as explained in METHODS.

## Secondary numerical variation

Two non-selected Infomap modularities differ between the supplied run and the validation environment (AEP703, orthophosphate; AIH047, conductivity). This secondary stochastic metric also changes the full HVG redundancy and PCA diagnostics slightly. Hydrological redundancy coefficients show small differences consistent with ranking numerical ties in diagnostic descriptors. These discrepancies are not labelled exact reproduction: `PLATFORM_VARIATION.csv` lists the affected cells and both values. The retained hydrological/HVG descriptor decisions, controlled associations, ranking scores, critical membership and upstream odds ratios are unchanged within the stated tolerance.

## Boundaries of validation

The comparison excludes degenerate, non-selected HITS diagnostics (`AUT`, `HBS`), explanatory/schema fields and exploratory t-SNE coordinates. Library/platform changes may affect such diagnostics and display coordinates. Scientific comparisons concern numerical results, not bytewise identity of PDF metadata or CSV serialisation.

Algorithmic checks and numerical reproduction do not establish causal validity, official ecological-status equivalence or all uncertainty assumptions. In particular, the odds-ratio ratio retains the documented covariance-omitting approximation. No official-status assessment dataset is included or represented as a completed benchmark. Refer to METHODS for these interpretation boundaries.

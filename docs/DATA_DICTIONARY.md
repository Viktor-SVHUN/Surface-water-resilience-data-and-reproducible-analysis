# Data dictionary and file conventions

## Source inputs

All four files are stored in `data/raw/`. Their source bytes and worksheet structures are unchanged; only filenames are standardised.

| File | Purpose | Primary identifiers |
| --- | --- | --- |
| `river_registry.csv` | River catalogue, including entries explicitly categorised as lakes | `VT_VOR_2015`, `BEF_VT_VOR` |
| `lake_registry.csv` | Lake catalogue and corresponding typological attributes | `VT_VOR_2015`, `BEF_VT_VOR` |
| `water_quality_annual_site_summaries.xlsx` | Annual site means, contributing sample counts and monitoring metadata | `VOR`, `KTJ`, `ANYAG_NEV_ANGOL`, `DATUM_YEAR` |
| `operational_thresholds.xlsx` | Type-level thresholds and favourable direction | `indicator`, `Lake`, `River`, `Boundary` |

VOR codes are strings, not row numbers. Numeric `ID` values in processed output are deterministic convenience indices; joins should use VOR. A recipient link runs from the row's source VOR to its designated `BEF_VT_VOR`.

### Monitoring workbook fields used by the program

| Source field | Meaning |
| --- | --- |
| `VOR` | Water-body or monitoring code |
| `KTJ` | Monitoring-site identifier |
| `ANYAG_NEV_ANGOL` | Source parameter name |
| `ATLAG` | Reported annual site mean |
| `DATUM_YEAR` | Calendar year of the summary |
| `MERTEKEGYSEG_NEV` | Source unit label |
| `MERESSZAM` | Reported sample count behind the annual summary |
| `FOLDRAJZI_SZELESSEG`, `FOLDRAJZI_HOSSZUSAG` | Site latitude/longitude |
| `VIZTEST_MEGNEVEZES` | Source water-body name |

### Registry fields used by the program

`VT_NEV_2015` provides the body name; `ALEGYSEG` the sub-basin label; `WB_CAT` the explicit category. `KOZV_VGY_MERET` and `TELJES_VGY_MERET` are direct/total catchment areas; `HOSSZ_MERT` is river length, `SUM_AREA_KM2` lake area and `ATLMELYS` mean depth. `ARTIFICIAL` and `MODIFIED` contain source Y/N flags. `ALT_CAT`, `SLOPE_CAT`, `GRAIN_CAT`, `GEOL_CAT`/`GEOL_CHEM_CAT` and `TIPUS_UJ`/`TIPUS` supply available typological classes. Missing or inapplicable fields remain explicit.

### Threshold workbook

`indicator` joins to the normalised parameter label. `Lake` and `River` contain thresholds. `Boundary=lower` means lower concentration is favourable (s=+1); `Boundary=higher` means higher concentration is favourable (s=−1). The supplied workbook is an operational input, not a complete legal-threshold catalogue.

## Processed data

- `monitoring_site_annual_summaries.csv`: one row per VOR/site/parameter/year/unit after consolidating source copies. `version_count` and `conflicting_values` retain consolidation diagnostics.
- `water_body_annual_summaries.csv`: one row per VOR/parameter/year/unit. `value` is the equal-site mean; `value_sample_weighted` is the alternative sample-count-weighted mean. `n_sites` counts available sites. `threshold`, `direction`, `adverse_deviation` and `eligible_local` specify threshold interpretation and eligibility.
- `hydrological_nodes.csv`: the computational graph's identifiers and metadata. `registered` and `observed` distinguish registry membership and monitoring availability. `Type` is the reconciled category; `catalogue_priority_type` is the category implied by river-first catalogue membership, retained only for reconciliation provenance.
- `recipient_links.csv`: directed `from`/`to` VOR pairs. The graph export is `hydrological_network.graphml`.
- `water_body_parameter_coverage.csv`: per-pair observed years, first/last years, internal gaps, site coverage and analysis-eligibility flags.

Mass parameters use g/L internally. Figures and LaTeX tables may display mg/L explicitly; conductivity is uS/cm, saturation is percent and accepted turbidity is NTU. Parameter names remain source labels so that joins are stable. An unavailable threshold is not replaced by the body's own mean.

## Main results

The inventory `RESULTS_INDEX.csv` lists every distributed result with its row/column count and purpose. `RESULTS_SCHEMA.csv` lists the actual columns and inferred storage types. Operational definitions, including units and censoring, are in `METHODS.md` and `results/tables/metric_dictionary.csv`.

### Associations

`rho` is a rank or residual-rank coefficient. `n`, `n_basins` and `n_controls` give the fitted sample/design dimensions. `p_cluster`, `q_BH`, `ci_low` and `ci_high` refer to the cluster procedure in METHODS. `p_asymptotic_raw` is supplementary. `controls` records the actual design columns; `status` explains effect-only, constant or insufficient-data cases. `FDR_significant` denotes q<0.05, not practical importance or causation. Missing inference is not represented as nonsignificance.

### Ranking and sensitivity

`score` is the average normalised longest adverse run. `eligible` identifies bodies meeting the parameter-count rule. `critical_cutoff`, `critical` and `rank` describe the baseline upper-tail rule. Scenario files distinguish threshold, weighting, common-parameter and single-year-omission alternatives. `rank_min`, `rank_max` and `critical_year_omission_fraction` are scenario summaries, not confidence limits.

### Upstream comparisons

Event tables identify `year`, `downstream_good` and `observed_upstream_adverse`. The matched table also gives observed/total direct upstream counts and complete-coverage flags. `conditional_compliance` is TR; `compliance_when_all_upstream_observed_clean` uses fully observed clean upstream years. `OR_MH`, `OR_MH_placebo` and `ROR_network_vs_placebo` denote the stratified connected, placebo and ratio estimates. Their interval and p/q fields follow the precise uncertainty definitions in METHODS.

### HVG metrics and perturbations

Metric suffixes `_D` and `_U` indicate directed and undirected projections. Deletion-replicate files retain VOR, parameter, deletion method/fraction, replicate index, metric, baseline/perturbed values and absolute changes. Concordance files aggregate across series, not across the nodes of a single graph. Series-length results specify `minimum_years` explicitly.

### Profiles

`quality_profile_input_matrix.csv` is the body-by-parameter matrix before median imputation/scaling; blank entries are the cells imputed for the profile analysis. `quality_profile_tsne.csv` contains VOR-linked exploratory coordinates and descriptive metadata. The parameter columns, filters and transformation are documented in METHODS and recorded in the manifest.

## Optional official-status input

No official-status dataset is included. The program optionally recognises `official_water_body_status.csv` in the selected input directory, with `VOR,status,assessment_year`, one row per VOR. If supplied, it creates a joined table only. Matching labels does not by itself establish temporal comparability or demonstrate incremental assessment value. No empty official-comparison placeholder is distributed.

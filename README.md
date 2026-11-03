# Surface-water resilience: data and reproducible analysis

Research software and data accompanying *Local and topological resilience of surface water bodies: how upstream non-compliance propagates through the Hungarian river network* by Viktor Sebestyén, Zsolt T. Kosztyán and János Abonyi.

The analysis combines annual physicochemical monitoring-site summaries, water-body registries, recipient links and operational thresholds. It computes local compliance/recovery indicators, persistence rankings, upstream-conditioned compliance, hydrological descriptors, horizontal visibility graph (HVG) metrics, controlled associations, typological profiles and sensitivity analyses.

## Contents

| Location | Contents |
| --- | --- |
| `data/raw/` | Four source inputs; their bytes are preserved |
| `src/surface_water_resilience.py` | Complete analysis program |
| `results/tables/` | Full machine-readable scientific results in CSV format |
| `results/figures/` | Individual vector-PDF and PNG figures |
| `results/latex/` | Numerical macros, LaTeX table fragments and figure exports |
| `results/manifest.json` | Configuration, software versions, input hashes and run metadata |
| `results/self_tests.csv` | Algorithmic self-test results |
| `docs/` | Methods, input schema, result inventory and provenance |
| `tests/` | Automated regression and packaging tests |
| `tools/verify_release.py` | Standard-library archive-integrity check |
| `MANIFEST.sha256` | SHA-256 checksums for all other distributed files |

CSV files are the canonical tabular outputs. They use UTF-8 with a byte-order mark, comma separators, decimal points and a header row. Empty numeric fields denote unavailable estimates, not zeros. Source parameter labels and Hungarian water-body names are preserved as data identifiers.

## Installation

The release was validated using Python 3.12. Use a dedicated environment. From the extracted `surface-water-resilience` directory:

```bash
python -m venv .venv
```

Activate it on Windows:

```powershell
.venv\Scripts\Activate.ps1
```

Or on Linux/macOS:

```bash
source .venv/bin/activate
```

Install the tested dependencies:

```bash
python -m pip install -r requirements.txt
```

The required Excel engine is openpyxl. The independent HVG reference package, ts2vg, is included in the tested requirements. No R installation, internet data download or external service is required to run the analysis after dependency installation.

## Reproduce the analysis

```bash
python src/surface_water_resilience.py --self-test-only
python src/surface_water_resilience.py
```

Default input paths are resolved relative to the program location, not the current working directory. A full run writes to a new `outputs/` directory, leaving the distributed `results/` unchanged. An existing non-empty output directory is rejected to prevent stale and newly generated files from being mixed. To run again, use another directory:

```bash
python src/surface_water_resilience.py --output-dir outputs_repeat
```

Alternative input directories must contain the same four canonical filenames:

```bash
python src/surface_water_resilience.py --input-dir /path/to/inputs --output-dir /path/to/new_outputs
```

In Spyder, open `src/surface_water_resilience.py` and run the file with no command-line arguments. The same portable defaults apply; editing machine-specific paths is unnecessary.

### Reduced diagnostic run

```bash
python src/surface_water_resilience.py --quick --output-dir outputs_quick
```

This deliberately reduces resampling and permutation settings and is **not** equivalent to the distributed full-configuration results. Other options are `--skip-tsne`, `--skip-sensitivity` and `--min-years N`. Run `--help` for the command-line interface. Configuration changes are recorded in each run's manifest.

## Verify the release

```bash
python tools/verify_release.py
python -m unittest discover -s tests -v
```

The integrity check requires only Python's standard library. The scientific tests require the installed analysis dependencies. `docs/VALIDATION.md` distinguishes algorithmic checks, full-run execution and agreement with the supplied numerical results. Hashes identify the exact inputs and program; they are not evidence that every scientific assumption is correct.

## Interpretation

- Monitoring records are **annual site summaries**, not individually dated samples. Annual body values give equal weight to sites after source-copy consolidation.
- The analysis covers physicochemical supporting elements, not whole-ecosystem resilience or full ecological status.
- Missing calendar years interrupt local runs. HVG visibility uses observation order, not the elapsed duration of a gap.
- Rankings are descriptive screening tools. Rank ranges are sensitivity envelopes, not confidence intervals.
- Upstream and network-position relationships are associational. The odds-ratio comparison uses the uncertainty approximation documented in `docs/METHODS.md`.
- The supplied inputs contain no official-status assessment table. No empirical comparison with an official RBMP status classification is represented as a result.
- t-SNE is exploratory; type is a colouring variable, not an embedding feature.

See `docs/METHODS.md` for full operational definitions and limitations, and `docs/DATA_DICTIONARY.md` for input fields, units and joins.

## LaTeX outputs

`results/latex/analysis_values.tex` contains numerical definitions named `\res...`; `table_*.tex` files contain tabular fragments. They can be used independently of any particular manuscript project. Tables require the LaTeX `booktabs` package; numbers and labels should be interpreted using the corresponding CSV table and methods description. Per-parameter association figures are also combined in the two supplementary heatmap PDFs.

## Citation and reuse

Use the author/title metadata in `CITATION.cff` when citing this research package. No unassigned DOI or publication identifier is inserted. Source-data provenance and reuse notices are described in `NOTICE.md`. This archive does not make an additional licence grant for the software or third-party data.

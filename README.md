# Paleoclimatic drivers of Acropora palmata populations

This repository contains:
- The compiled coral and climate dataset
- Supplementary Quarto document (R script) with exploratory analysis, KDE generation, bootstrapping, and GAM fitting steps

## Data Availability and Licensing
The **compiled dataset** integrates information from multiple published sources.  
Original datasets remain under the copyright and licensing terms of their respective authors and publishers.  
This compilation is released under the **Creative Commons Attribution 4.0 International (CC BY 4.0)** license, which permits reuse provided appropriate credit is given to both the compiler and the original sources.  
Please cite the original publications when reusing specific subsets of data.

## Code Availability
All R scripts and Quarto documents are released under the **MIT License**.  
This allows reuse, modification, and distribution with attribution.

## Which license covers which files
| License | Files | License text |
|---|---|---|
| **CC BY 4.0** | Data and outputs: `*.xlsx`, `*.csv`, `*.rds`, figures (`*.pdf`), tables (`*.docx`), the rendered report (`regional_analysis_Final.html`, `regional_analysis_Final_files/`) and the contents of `archive/`, `kaufman2020_CPS/`, `florida_local_SST_records/`, `M3_sensitivity/` and `sensitivity_round3/` that are not code | [`LICENSE`](LICENSE) |
| **MIT** | Code: `regional_analysis_Final.qmd`, `sensitivity_round3/round3_functions.R`, `M3_sensitivity/*.py` and `*.Rproj` | [`LICENSE-MIT`](LICENSE-MIT) |

Third-party data deposited here (e.g. the Marine20 calibration curve, the Kaufman et al. 2020 CPS ensemble, Florida SST records, Schmitt et al. 2025 event layers) remain under the terms of their original publications; cite the original sources when reusing them.

## Repository structure
- Repository root: current inputs read by `regional_analysis_Final.qmd`, rendered report, figures and tables.
- `archive/`: superseded versions of the input tables, kept for traceability and not read by the analysis (see `archive/README.md`).

## Citation
If you use this dataset or code, please cite:  
Medina-Valmaseda, A.E. et al. (2026). *Paleoclimatic drivers of Acropora palmata populations across Caribbean Reefs over the last 6000 Years*. Scientific Reports (in revision).

# archive/

Superseded versions of the input tables, kept for traceability. **None of these files is read by `regional_analysis_Final.qmd`**; the current inputs are in the repository root.

| File | What it is | Superseded by |
|---|---|---|
| `ages_kernel_reg_updt.xlsx` | First version of the fossil chronology (updt1). Baseline against which the changes listed in sheet `Trazabilidad` of the current chronology are documented, and input of `M3_sensitivity/build_updt2.py`. | `ages_kernel_reg_updt2.xlsx` |
| `df_GAM130326_regimes_v3.xlsx` | Climate table v3: hurricane windows counted from 2022 CE taken in order as BP bins (see Section 1.2 of the Quarto document). | `df_GAM130326_regimes_v5.xlsx` |
| `df_GAM130326_regimes_v4.xlsx` | Climate table v4: hurricane metrics recounted in BP bins; temperature from a single CPS ensemble member (V1). Used for the runs documented in `sensitivity_round3/README.md`. | `df_GAM130326_regimes_v5.xlsx` |
| `df_GAM130326_regimes_v4c.xlsx` | Variant of v4 with the youngest bin extended to 2022 CE. | `df_GAM130326_regimes_v5c.xlsx` |

The Quarto document mentions some of these files by name when it describes the corrections made between versions; they are kept here so that those descriptions can be checked.

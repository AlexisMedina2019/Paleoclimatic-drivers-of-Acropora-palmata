# CPS temperature reconstruction, 0–30°N (Kaufman et al. 2020)

`0to30500-6bands-100yr2kwindow.csv`: 500-member ensemble of the Composite Plus Scale (CPS) reconstruction for the 0–30°N latitude band (all longitudes). Columns: `binAges` (bin centres, yr BP, 1950 CE; 100-yr bins, 0–12,000 yr BP) and members `V1`–`V500` (°C anomalies).

Source: Kaufman, D. et al. (2020) Holocene global mean surface temperature, a multi-method reconstruction approach. *Scientific Data* 7, 201. File `ScientificDataAnalysis/CPS/0to30500-6bands-100yr2kwindow.csv` of the authors' code repository, github.com/nickmckay/Temperature12k (commit ea85457, 2020-06-10; MIT licence; archived release doi:10.5281/zenodo.3888590). Published reconstructions: NOAA Paleoclimatology, doi:10.25921/vzys-1280.

Use in this repository (round 3, 2026-10-30): `GMHT_CPS0_30N` in `df_GAM130326_regimes_v5.xlsx` is the ensemble median, re-referenced to 0 at 100 yr BP and interpolated at the centres of the model bins; `dT_raw` is the change from the preceding (older) century into each bin. Section 4.7.5 of `regional_analysis_Final.qmd` refits the models with individual members. The previous table (v4) used member V1 (see sheet `CPS_traceability`).

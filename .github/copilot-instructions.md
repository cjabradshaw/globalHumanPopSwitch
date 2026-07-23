# Repository guidance

## Workflow and validation

This is an R analysis repository, not an R package or application. There are no build, test, lint, or test-runner configurations, so there is no full-suite or single-test command.

The available lightweight validation is an R syntax check:

```sh
Rscript -e 'for (file in c("scripts/humpoptransitionR1.R", "scripts/new_lmer_AIC_tables3.r", "scripts/r.squared.R")) parse(file)'
```

`README.md` lists the analysis dependencies. The main script additionally imports `sandwich`:

```r
boot, dismo, gbm, ggplot2, gridExtra, lmtest, performance, plotrix,
sjPlot, sandwich, tmvnsim, truncnorm, wCorr
```

Do not present `Rscript scripts/humpoptransitionR1.R` as a supported clean-checkout command. The script uses bare relative paths for every input and helper, while tracked data and scripts reside in separate directories. It also expects the untracked country inputs `NER.csv`, `COD.csv`, `MLI.csv`, `TCD.csv`, `AGO.csv`, `NGA.csv`, `BDI.csv`, `BFA.csv`, `GMB.csv`, and `UGA.csv`. A full run therefore requires a prepared working directory containing the main script, both helper scripts, all input CSVs, and any missing country datasets; it writes generated CSVs into that same directory.

## Architecture

`scripts/humpoptransitionR1.R` is the sequential analysis entry point. It loads packages, sources the two local helpers, defines inline statistical helpers, then executes the complete analysis from top to bottom:

- It combines historical and UN population series, calculates log growth rates, and fits Ricker and Gompertz models for three global phases: 1800--1949, 1950--1961, and 1962 onwards.
- It repeats the Ricker/Gompertz workflow for UN regions, China, ten high-fertility countries, and East/South-East Asia excluding China. Each produces fit and input CSV artefacts.
- It joins demographic, regional classification, and PPP-GDP data for regional demographic summaries.
- It models temperature anomaly, ecological footprint, and CO2-e emissions against per-capita consumption and population using linear models, Newey-West errors, and boosted regression trees with 1,000 bootstrap iterations.
- It finishes with global age-structure calculations.
- It then reads published WID global and selected regional Gini series from `WID_DATA_DIR`, merges their 1980--2023 overlap with the global environmental indicators, and exports trend-controlled and first-difference Newey-West models.
- It fits six resampled WID boosted regression tree models for temperature anomaly, ecological footprint, and CO2-e emissions, with population, per-capita consumption, and one of two global Gini definitions as predictors.
- It downloads and caches regional energy use and ecological footprint for the WID-defined Sub-Saharan Africa, Latin America, and MENA aggregates in `data/WID_regional_energy_footprint.csv`.
- It sums Climate Watch country greenhouse-gas sectors, supplied through OWID, into cached regional CO2-e totals and fits six regional CO2-e BRTs. This source covers 1990--2023; regional climate BRTs remain out of scope until an aligned response is sourced.

`scripts/new_lmer_AIC_tables3.r` provides AICc/AIC/BIC model-comparison helpers (`aicW` and related functions). `scripts/r.squared.R` defines the `r.squared` S3 generic and methods for `lm`, `merMod`, and `lme` models. Input CSVs in `data/` are the authoritative datasets; `www/` only supplies README branding assets.

## Repository-specific conventions

- Preserve CSV headers, units, ordering, and names. The main script uses direct column names and positional column indexes (especially for boosted regression tree inputs), so a schema-only change can silently alter an analysis.
- Keep the three global phase boundaries and their uncertainty assumptions aligned: 5% for phase 1, 2% for phase 2, and 1% for phase 3. Phase 2 deliberately removes the 1950 transition for its Ricker fit.
- Treat generated CSVs as run artefacts. Names encode their analysis scope, such as `Phase3_Ricker_Gompertz_fits.csv`, `ESEA_Ricker_Gompertz_fits.csv`, and `NETA.BRT.boot.pred.med.csv`.
- The main script sources `new_lmer_AIC_tables3.R`, but the tracked helper is named `new_lmer_AIC_tables3.r`. Do not change either spelling or case without updating the source call; the mismatch breaks case-sensitive environments.
- Bootstrap and resampling sections use random sampling without a repository-wide seed. When comparing changed analytical results, control the R RNG state explicitly and document the seed used for that run.
- The WID block requires `WID_DATA_DIR` to name the external directory containing `WID_data_XX.csv` files. It uses WID's published World, Sub-Saharan Africa, Latin America, and MENA aggregates; do not construct the remaining project regions by averaging national Gini coefficients.
- The WID BRT block runs 1,000 bootstrap resamples per outcome/Gini combination by default. `WID_BRT_ITER` can reduce this only for exploratory runs; retain the default for reported analyses.
- Regional energy use is population-weighted from OWID's per-capita primary-energy series. Regional ecological footprint is summed from Global Footprint Network country accounts; set `WID_REFRESH_REGIONAL_ENVIRONMENT=true` to refresh the cached source-derived file.
 - Regional greenhouse-gas emissions are summed across Climate Watch sectors in OWID's `greenhouse-gas-emissions-by-sector.csv`. They are territorial tonnes of CO2-equivalents, cover 1990--2023, and are cached in `data/WID_regional_emissions.csv`; set `WID_REFRESH_REGIONAL_EMISSIONS=true` to refresh them.

# Long COVID Paxlovid Thesis Analysis Code

This folder contains cleaned R scripts for the thesis analysis.

## Scripts

- `01_imputation_iptw_regression_models.R`: multiple imputation, propensity-score/IPTW models, covariate balance, and regression summaries.
- `02_entropy_balancing_aipw_missingness.R`: entropy balancing, AIPW estimates, and observed predictors of missingness.
- `03_mediation_complete_case_mi.R`: complete-case and multiple-imputation causal mediation analyses for candidate mediators.

## Data

The scripts expect the original Excel input files to be available in the R working directory. Data files are not included here because they may contain sensitive clinical information.

## Reproducibility

Each script checks that required R packages are installed and stops with a clear message if any package is missing. Package installation is intentionally not performed inside the scripts.

## Suggested run order

1. Run `01_imputation_iptw_regression_models.R`.
2. Run `02_entropy_balancing_aipw_missingness.R`.
3. Run `03_mediation_complete_case_mi.R`.

## Privacy

Do not commit raw patient-level data, Excel files, or identifiable output to a public repository.

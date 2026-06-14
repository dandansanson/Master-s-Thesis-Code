# Multiple-imputation, IPTW, and regression analyses
# File: 01_imputation_iptw_regression_models.R
#
# Purpose:
#   Reproducible analysis script for the Long COVID/Paxlovid thesis.
#
# Notes:
#   - Input data files are expected in the working directory.
#   - This script does not install packages automatically.
#   - Run the scripts in numerical order unless using one script independently.

# Required packages ------------------------------------------------------------
required_packages <- c(
  "readxl",
  "mice",
  "dplyr",
  "tidyr",
  "cobalt",
  "broom",
  "openxlsx",
  "ggplot2",
  "naniar",
  "VIM",
  "WeightIt"
)

missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) {
  stop(
    "Install the following packages before running this script: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

invisible(lapply(required_packages, library, character.only = TRUE))

# Settings ------------------------------------------------------------
file_name <- "For Imputation.xlsx"
sheet_name <- "for imputation"

m_imputations <- 20
max_iterations <- 20
seed_num <- 123

outcomes <- c(
  "new_improve_grorcoop",
  "new_improved_grandcoop",
  "improve_clinician",
  "improve_patient"
)

outcome_labels <- c(
  "Delta Grade or COOP",
  "Delta Grade and COOP",
  "Clinician-reported improvement",
  "Patient-reported improvement"
)

treatment_var <- "paxlovid_treatment"

# Covariates only
# general_treatment, dusoi, covid_wave, and baseline_coop are excluded
covariates <- c(
  "age",
  "sex",
  "Number_comorb",
  "Vaccination_number"
)

vars_needed <- c(outcomes, treatment_var, covariates)

# Load data ------------------------------------------------------------
data_raw <- readxl::read_excel(file_name, sheet = sheet_name)

missing_vars <- setdiff(vars_needed, names(data_raw))

if (length(missing_vars) > 0) {
  stop(
    "These variables are missing from the Excel file: ",
    paste(missing_vars, collapse = ", ")
  )
}

df <- data_raw %>%
  dplyr::select(dplyr::all_of(vars_needed))

cat("\nData dimensions:\n")
print(dim(df))

cat("\nVariable names:\n")
print(names(df))

# Check missingness pattern first ------------------------------------------------------------
missing_table <- data.frame(
  Variable = names(df),
  Missing_n = colSums(is.na(df)),
  Missing_pct = round(colMeans(is.na(df)) * 100, 2)
) %>%
  dplyr::arrange(desc(Missing_pct))

cat("\nMissingness summary:\n")
print(missing_table)

missing_pattern <- mice::md.pattern(df, plot = FALSE)

cat("\nMissingness pattern:\n")
print(missing_pattern)

if (anyNA(df)) {
  VIM::aggr(
    df,
    col = c("navyblue", "red"),
    numbers = TRUE,
    sortVars = TRUE,
    labels = names(df),
    cex.axis = 0.7,
    gap = 3,
    ylab = c("Missing data", "Pattern")
  )
  
  print(naniar::vis_miss(df))
}

# Recode variable types ------------------------------------------------------------
factor_vars <- c(
  outcomes,
  "paxlovid_treatment",
  "sex"
)

numeric_vars <- c(
  "age",
  "Number_comorb",
  "Vaccination_number"
)

df <- df %>%
  dplyr::mutate(
    dplyr::across(dplyr::all_of(factor_vars), as.factor),
    dplyr::across(dplyr::all_of(numeric_vars), ~ as.numeric(.))
  )

cat("\nFactor levels:\n")
print(lapply(df[, factor_vars], levels))

cat("\nFinal data structure:\n")
str(df)

# Multiple imputation setup ------------------------------------------------------------
choose_factor_method <- function(x) {
  if (!any(is.na(x))) return("")
  
  nlev <- nlevels(droplevels(x))
  
  if (nlev == 2) return("logreg")
  if (nlev > 2) return("polyreg")
  
  return("")
}

if (anyNA(df)) {
  
  init <- mice::mice(df, maxit = 0, printFlag = FALSE)
  meth <- init$method
  pred <- init$predictorMatrix
  
  for (v in factor_vars) {
    meth[v] <- choose_factor_method(df[[v]])
  }
  
  for (v in numeric_vars) {
    meth[v] <- ifelse(any(is.na(df[[v]])), "pmm", "")
  }
  
  diag(pred) <- 0
  
  cat("\nImputation methods:\n")
  print(meth)
  
  imp <- mice::mice(
    df,
    m = m_imputations,
    maxit = max_iterations,
    method = meth,
    predictorMatrix = pred,
    seed = seed_num,
    printFlag = TRUE
  )
  
  plot(imp)
  
  imp_list <- mice::complete(imp, action = "all")
  
} else {
  
  cat("\nNo missing data detected. Skipping multiple imputation.\n")
  imp <- NULL
  imp_list <- rep(list(df), m_imputations)
}

# Convert Paxlovid treatment to 0/1 ------------------------------------------------------------
make_binary_numeric <- function(x) {
  x_chr <- tolower(trimws(as.character(x)))
  
  positive_values <- c(
    "1", "yes", "y", "true", "treated",
    "paxlovid", "paxlovid treatment", "paxlovid_treatment"
  )
  
  negative_values <- c(
    "0", "no", "n", "false", "untreated",
    "control", "standard", "usual care", "no paxlovid", "non-paxlovid"
  )
  
  dplyr::case_when(
    is.na(x_chr) ~ NA_real_,
    x_chr %in% positive_values ~ 1,
    x_chr %in% negative_values ~ 0,
    TRUE ~ NA_real_
  )
}

add_treatment_numeric <- function(dat) {
  dat$pax_num <- make_binary_numeric(dat[[treatment_var]])
  
  bad_values <- unique(as.character(
    dat[[treatment_var]][is.na(dat$pax_num) & !is.na(dat[[treatment_var]])]
  ))
  
  if (length(bad_values) > 0) {
    stop(
      "Some paxlovid_treatment values were not recognized: ",
      paste(bad_values, collapse = ", "),
      "\nEdit make_binary_numeric() to match your coding."
    )
  }
  
  dat
}

imp_list_model <- lapply(imp_list, add_treatment_numeric)

cat("\nTreatment coding check, first imputed dataset:\n")
print(table(
  imp_list_model[[1]][[treatment_var]],
  imp_list_model[[1]]$pax_num,
  useNA = "ifany"
))

# Helper functions for model fitting and pooling ------------------------------------------------------------
covariate_rhs <- paste(covariates, collapse = " + ")

adjusted_rhs <- paste(
  c("pax_num", covariates),
  collapse = " + "
)

fit_model_set <- function(dat_list, outcome_var, rhs, weights_col = NULL) {
  
  fit_list <- lapply(dat_list, function(dat) {
    
    model_formula <- as.formula(
      paste(outcome_var, "~", rhs)
    )
    
    if (is.null(weights_col)) {
      glm(
        model_formula,
        data = dat,
        family = binomial()
      )
    } else {
      glm(
        model_formula,
        data = dat,
        weights = dat[[weights_col]],
        family = quasibinomial()
      )
    }
  })
  
  pooled <- mice::pool(mice::as.mira(fit_list))
  pooled_summary <- summary(
    pooled,
    conf.int = TRUE,
    exponentiate = TRUE
  )
  
  list(
    fits = fit_list,
    pooled = pooled,
    summary = pooled_summary
  )
}

save_model_summary <- function(results_list, labels, method_name) {
  dplyr::bind_rows(lapply(seq_along(results_list), function(i) {
    results_list[[i]]$summary %>%
      dplyr::mutate(
        Outcome = labels[i],
        Method = method_name
      )
  }))
}

extract_pax_effects <- function(results_list, labels, method_name) {
  dplyr::bind_rows(lapply(seq_along(results_list), function(i) {
    results_list[[i]]$summary %>%
      dplyr::filter(term == "pax_num") %>%
      dplyr::transmute(
        Outcome = labels[i],
        Method = method_name,
        OR = round(estimate, 3),
        CI_low = round(`2.5 %`, 3),
        CI_high = round(`97.5 %`, 3),
        p_value = round(p.value, 3)
      )
  }))
}

# Model without weighting - Treatment + covariates, no IPTW ------------------------------------------------------------
unweighted_results <- lapply(outcomes, function(outcome) {
  fit_model_set(
    dat_list = imp_list_model,
    outcome_var = outcome,
    rhs = adjusted_rhs,
    weights_col = NULL
  )
})

full_unweighted <- save_model_summary(
  unweighted_results,
  outcome_labels,
  "Unweighted adjusted"
)

cat("\nUnweighted adjusted model summaries:\n")
print(full_unweighted)

# Check covariate balance before IPTW ------------------------------------------------------------
balance_formula <- as.formula(
  paste("pax_num ~", covariate_rhs)
)

balance_to_df <- function(bal_obj) {
  out <- as.data.frame(bal_obj$Balance)
  out$Covariate <- rownames(out)
  rownames(out) <- NULL
  out %>% dplyr::relocate(Covariate)
}

get_balance_df <- function(dat, weights = NULL, phase, imputation) {
  
  if (is.null(weights)) {
    bal <- cobalt::bal.tab(
      balance_formula,
      data = dat,
      binary = "std",
      continuous = "std",
      s.d.denom = "pooled"
    )
  } else {
    bal <- cobalt::bal.tab(
      balance_formula,
      data = dat,
      weights = weights,
      method = "weighting",
      un = TRUE,
      binary = "std",
      continuous = "std",
      s.d.denom = "pooled"
    )
  }
  
  balance_to_df(bal) %>%
    dplyr::mutate(
      Phase = phase,
      Imputation = imputation,
      .before = 1
    )
}

pre_balance_first <- cobalt::bal.tab(
  balance_formula,
  data = imp_list_model[[1]],
  binary = "std",
  continuous = "std",
  s.d.denom = "pooled"
)

cat("\nCovariate balance before IPTW, first imputed dataset:\n")
print(pre_balance_first)

pre_love <- cobalt::love.plot(
  pre_balance_first,
  stats = "mean.diffs",
  abs = TRUE,
  thresholds = c(m = 0.1)
)

print(pre_love)

pre_balance_all <- dplyr::bind_rows(lapply(seq_along(imp_list_model), function(i) {
  get_balance_df(
    dat = imp_list_model[[i]],
    weights = NULL,
    phase = "Before IPTW",
    imputation = i
  )
}))

# Estimate propensity scores and IPTW ------------------------------------------------------------
ps_formula <- as.formula(
  paste("pax_num ~", covariate_rhs)
)

add_ps_weights <- function(dat) {
  
  ps_model <- glm(
    ps_formula,
    data = dat,
    family = binomial()
  )
  
  dat$ps <- predict(ps_model, newdata = dat, type = "response")
  
  dat$ps <- pmin(pmax(dat$ps, 0.01), 0.99)
  
  p_treat <- mean(dat$pax_num == 1, na.rm = TRUE)
  
  dat$sw <- ifelse(
    dat$pax_num == 1,
    p_treat / dat$ps,
    (1 - p_treat) / (1 - dat$ps)
  )
  
  cap <- quantile(dat$sw, 0.99, na.rm = TRUE)
  dat$sw_trunc <- pmin(dat$sw, as.numeric(cap))
  
  list(
    data = dat,
    ps_model = ps_model
  )
}

ps_list <- lapply(imp_list_model, add_ps_weights)

imp_list_w <- lapply(ps_list, `[[`, "data")
ps_models <- lapply(ps_list, `[[`, "ps_model")

weight_summary <- dplyr::bind_rows(lapply(seq_along(imp_list_w), function(i) {
  dat <- imp_list_w[[i]]
  
  data.frame(
    Imputation = i,
    n = nrow(dat),
    treated_n = sum(dat$pax_num == 1, na.rm = TRUE),
    untreated_n = sum(dat$pax_num == 0, na.rm = TRUE),
    sw_min = min(dat$sw, na.rm = TRUE),
    sw_mean = mean(dat$sw, na.rm = TRUE),
    sw_median = median(dat$sw, na.rm = TRUE),
    sw_p95 = as.numeric(quantile(dat$sw, 0.95, na.rm = TRUE)),
    sw_p99 = as.numeric(quantile(dat$sw, 0.99, na.rm = TRUE)),
    sw_max = max(dat$sw, na.rm = TRUE),
    sw_trunc_max = max(dat$sw_trunc, na.rm = TRUE)
  )
}))

cat("\nWeight summary:\n")
print(weight_summary)

weight_plot <- ggplot(imp_list_w[[1]], aes(x = sw_trunc)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of Truncated Stabilized Weights",
    x = "Truncated stabilized weight",
    y = "Frequency"
  )

print(weight_plot)

ps_overlap_plot <- ggplot(
  imp_list_w[[1]],
  aes(x = ps, fill = factor(pax_num))
) +
  geom_density(alpha = 0.4) +
  theme_minimal() +
  labs(
    title = "Propensity Score Overlap",
    x = "Propensity score",
    fill = "Paxlovid"
  )

print(ps_overlap_plot)

# IPTW model without doubly robust estimator - Outcome model = treatment only, weighted by IPTW ------------------------------------------------------------
iptw_results <- lapply(outcomes, function(outcome) {
  fit_model_set(
    dat_list = imp_list_w,
    outcome_var = outcome,
    rhs = "pax_num",
    weights_col = "sw_trunc"
  )
})

full_iptw <- save_model_summary(
  iptw_results,
  outcome_labels,
  "IPTW only"
)

cat("\nIPTW-only model summaries:\n")
print(full_iptw)

# Check covariate balance after IPTW ------------------------------------------------------------
post_balance_first <- cobalt::bal.tab(
  balance_formula,
  data = imp_list_w[[1]],
  weights = imp_list_w[[1]]$sw_trunc,
  method = "weighting",
  un = TRUE,
  binary = "std",
  continuous = "std",
  s.d.denom = "pooled"
)

cat("\nCovariate balance after IPTW, first imputed dataset:\n")
print(post_balance_first)

post_love <- cobalt::love.plot(
  post_balance_first,
  stats = "mean.diffs",
  abs = TRUE,
  thresholds = c(m = 0.1)
)

print(post_love)

post_balance_all <- dplyr::bind_rows(lapply(seq_along(imp_list_w), function(i) {
  get_balance_df(
    dat = imp_list_w[[i]],
    weights = imp_list_w[[i]]$sw_trunc,
    phase = "After IPTW",
    imputation = i
  )
}))

balance_all <- dplyr::bind_rows(
  pre_balance_all,
  post_balance_all
)

# Doubly robust estimator - Outcome model = treatment + covariates, weighted by IPTW ------------------------------------------------------------
dr_results <- lapply(outcomes, function(outcome) {
  fit_model_set(
    dat_list = imp_list_w,
    outcome_var = outcome,
    rhs = adjusted_rhs,
    weights_col = "sw_trunc"
  )
})

full_dr <- save_model_summary(
  dr_results,
  outcome_labels,
  "Doubly robust"
)

cat("\nDoubly robust model summaries:\n")
print(full_dr)

# Compare Paxlovid effect across methods ------------------------------------------------------------
comparison_table <- dplyr::bind_rows(
  extract_pax_effects(
    unweighted_results,
    outcome_labels,
    "Unweighted adjusted"
  ),
  extract_pax_effects(
    iptw_results,
    outcome_labels,
    "IPTW only"
  ),
  extract_pax_effects(
    dr_results,
    outcome_labels,
    "Doubly robust"
  )
)

comparison_wide <- comparison_table %>%
  dplyr::mutate(
    CI = paste0("(", CI_low, ", ", CI_high, ")")
  ) %>%
  dplyr::select(Outcome, Method, OR, CI, p_value) %>%
  tidyr::pivot_wider(
    names_from = Method,
    values_from = c(OR, CI, p_value)
  )

cat("\nComparison of Paxlovid effect across methods:\n")
print(comparison_table)

cat("\nWide comparison table:\n")
print(comparison_wide)

# Save results ------------------------------------------------------------
missing_pattern_df <- as.data.frame(missing_pattern)
missing_pattern_df$Pattern <- rownames(missing_pattern_df)
rownames(missing_pattern_df) <- NULL

missing_pattern_df <- missing_pattern_df %>%
  dplyr::relocate(Pattern)

full_all <- dplyr::bind_rows(
  full_unweighted,
  full_iptw,
  full_dr
)

write.csv(missing_table, "missingness_summary.csv", row.names = FALSE)
write.csv(missing_pattern_df, "missingness_pattern.csv", row.names = FALSE)
write.csv(weight_summary, "iptw_weight_summary.csv", row.names = FALSE)
write.csv(balance_all, "covariate_balance_before_after_iptw.csv", row.names = FALSE)
write.csv(comparison_table, "paxlovid_effect_comparison_long.csv", row.names = FALSE)
write.csv(comparison_wide, "paxlovid_effect_comparison_wide.csv", row.names = FALSE)
write.csv(full_all, "full_model_summaries.csv", row.names = FALSE)

ggsave(
  "weights_distribution_imputation1.png",
  weight_plot,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  "ps_overlap_imputation1.png",
  ps_overlap_plot,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  "balance_before_iptw_imputation1.png",
  pre_love,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  "balance_after_iptw_imputation1.png",
  post_love,
  width = 8,
  height = 5,
  dpi = 300
)

wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb, "Missingness")
openxlsx::writeData(wb, "Missingness", missing_table)

openxlsx::addWorksheet(wb, "Missing Pattern")
openxlsx::writeData(wb, "Missing Pattern", missing_pattern_df)

openxlsx::addWorksheet(wb, "Weight Summary")
openxlsx::writeData(wb, "Weight Summary", weight_summary)

openxlsx::addWorksheet(wb, "Balance Before After")
openxlsx::writeData(wb, "Balance Before After", balance_all)

openxlsx::addWorksheet(wb, "Pax Comparison Long")
openxlsx::writeData(wb, "Pax Comparison Long", comparison_table)

openxlsx::addWorksheet(wb, "Pax Comparison Wide")
openxlsx::writeData(wb, "Pax Comparison Wide", comparison_wide)

openxlsx::addWorksheet(wb, "Full Model Summaries")
openxlsx::writeData(wb, "Full Model Summaries", full_all)

openxlsx::saveWorkbook(
  wb,
  "missingness_unweighted_iptw_dr_results_no_general_treatment.xlsx",
  overwrite = TRUE
)

# Session info ------------------------------------------------------------
sessionInfo()

# Additional WeightIt balancing weights ---------------------------------------

# Covariates only
# general_treatment, dusoi, covid_wave, and baseline_coop are excluded
covariates <- c(
  "age",
  "sex",
  "Number_comorb",
  "Vaccination_number"
)

covariate_rhs <- paste(covariates, collapse = " + ")

balance_formula <- as.formula(
  paste("pax_num ~", covariate_rhs)
)

add_balancing_weights <- function(dat) {
  
  # Method 1: propensity score weighting with logistic regression
  # Add squares to improve balance for numeric covariates
  w_obj <- WeightIt::weightit(
    pax_num ~ age + I(age^2) + sex +
      Number_comorb + I(Number_comorb^2) +
      Vaccination_number + I(Vaccination_number^2),
    data = dat,
    method = "glm",
    estimand = "ATE",
    stabilize = TRUE
  )
  
  dat$sw <- w_obj$weights
  
  # Optional truncation to avoid unstable extreme weights
  cap_low <- quantile(dat$sw, 0.01, na.rm = TRUE)
  cap_high <- quantile(dat$sw, 0.99, na.rm = TRUE)
  
  dat$sw_trunc <- pmin(
    pmax(dat$sw, as.numeric(cap_low)),
    as.numeric(cap_high)
  )
  
  list(
    data = dat,
    weightit_object = w_obj
  )
}

weightit_list <- lapply(imp_list_model, add_balancing_weights)

imp_list_w <- lapply(weightit_list, `[[`, "data")
weightit_objects <- lapply(weightit_list, `[[`, "weightit_object")

# Check covariate balance after weighting ------------------------------------------------------------

bal_before_after <- cobalt::bal.tab(
  balance_formula,
  data = imp_list_w[[1]],
  weights = imp_list_w[[1]]$sw_trunc,
  method = "weighting",
  un = TRUE,
  binary = "std",
  continuous = "std",
  s.d.denom = "pooled"
)

print(bal_before_after)

love_plot <- cobalt::love.plot(
  bal_before_after,
  stats = "mean.diffs",
  abs = TRUE,
  thresholds = c(m = 0.1),
  var.order = "unadjusted"
)

print(love_plot)

### Without weighting
# Adjusted model without weighting - Treatment + covariates - No general_treatment, dusoi, covid_wave, or baseline_coop ------------------------------------------------------------

covariates <- c(
  "age",
  "sex",
  "Number_comorb",
  "Vaccination_number"
)

outcomes <- c(
  "new_improve_grorcoop",
  "new_improved_grandcoop",
  "improve_clinician",
  "improve_patient"
)

outcome_labels <- c(
  "Delta Grade or COOP",
  "Delta Grade and COOP",
  "Clinician-reported improvement",
  "Patient-reported improvement"
)

adjusted_rhs <- paste(
  c("pax_num", covariates),
  collapse = " + "
)

fit_adjusted_model <- function(imp_list_model, outcome_var) {
  
  fit_list <- lapply(imp_list_model, function(dat) {
    
    model_formula <- as.formula(
      paste(outcome_var, "~", adjusted_rhs)
    )
    
    glm(
      model_formula,
      data = dat,
      family = binomial()
    )
  })
  
  pooled <- mice::pool(mice::as.mira(fit_list))
  
  pooled_summary <- summary(
    pooled,
    conf.int = TRUE,
    exponentiate = TRUE
  )
  
  list(
    fits = fit_list,
    pooled = pooled,
    summary = pooled_summary
  )
}

adjusted_results <- lapply(outcomes, function(outcome) {
  fit_adjusted_model(
    imp_list_model = imp_list_model,
    outcome_var = outcome
  )
})

names(adjusted_results) <- outcome_labels

# Print full adjusted model results ------------------------------------------------------------

for (i in seq_along(adjusted_results)) {
  cat("\n====================================================\n")
  cat("Adjusted model:", outcome_labels[i], "\n")
  cat("====================================================\n")
  print(adjusted_results[[i]]$summary)
}

# Extract Paxlovid effect only ------------------------------------------------------------

adjusted_pax_effect <- dplyr::bind_rows(lapply(seq_along(adjusted_results), function(i) {
  
  adjusted_results[[i]]$summary %>%
    dplyr::filter(term == "pax_num") %>%
    dplyr::transmute(
      Outcome = outcome_labels[i],
      OR = round(estimate, 3),
      CI_low = round(`2.5 %`, 3),
      CI_high = round(`97.5 %`, 3),
      p_value = round(p.value, 3)
    )
}))

print(adjusted_pax_effect)

# Save adjusted model results ------------------------------------------------------------

adjusted_full_results <- dplyr::bind_rows(lapply(seq_along(adjusted_results), function(i) {
  
  adjusted_results[[i]]$summary %>%
    dplyr::mutate(
      Outcome = outcome_labels[i],
      Method = "Adjusted model"
    )
}))

write.csv(
  adjusted_pax_effect,
  "adjusted_model_paxlovid_effect.csv",
  row.names = FALSE
)

write.csv(
  adjusted_full_results,
  "adjusted_model_full_results.csv",
  row.names = FALSE
)

# Weighted adjusted model - Outcome model = treatment + covariates - Weights = sw_trunc - No general_treatment, dusoi, covid_wave, or baseline_coop ------------------------------------------------------------

covariates <- c(
  "age",
  "sex",
  "Number_comorb",
  "Vaccination_number"
)

outcomes <- c(
  "new_improve_grorcoop",
  "new_improved_grandcoop",
  "improve_clinician",
  "improve_patient"
)

outcome_labels <- c(
  "Delta Grade or COOP",
  "Delta Grade and COOP",
  "Clinician-reported improvement",
  "Patient-reported improvement"
)

adjusted_rhs <- paste(
  c("pax_num", covariates),
  collapse = " + "
)

fit_weighted_adjusted_model <- function(imp_list_w, outcome_var) {
  
  fit_list <- lapply(imp_list_w, function(dat) {
    
    model_formula <- as.formula(
      paste(outcome_var, "~", adjusted_rhs)
    )
    
    glm(
      model_formula,
      data = dat,
      weights = sw_trunc,
      family = quasibinomial()
    )
  })
  
  pooled <- mice::pool(mice::as.mira(fit_list))
  
  pooled_summary <- summary(
    pooled,
    conf.int = TRUE,
    exponentiate = TRUE
  )
  
  list(
    fits = fit_list,
    pooled = pooled,
    summary = pooled_summary
  )
}

weighted_adjusted_results <- lapply(outcomes, function(outcome) {
  fit_weighted_adjusted_model(
    imp_list_w = imp_list_w,
    outcome_var = outcome
  )
})

names(weighted_adjusted_results) <- outcome_labels

# Print full weighted adjusted model results ------------------------------------------------------------

for (i in seq_along(weighted_adjusted_results)) {
  cat("\n====================================================\n")
  cat("Weighted adjusted model:", outcome_labels[i], "\n")
  cat("====================================================\n")
  print(weighted_adjusted_results[[i]]$summary)
}

# Extract Paxlovid effect only ------------------------------------------------------------

weighted_adjusted_pax_effect <- dplyr::bind_rows(lapply(seq_along(weighted_adjusted_results), function(i) {
  
  weighted_adjusted_results[[i]]$summary %>%
    dplyr::filter(term == "pax_num") %>%
    dplyr::transmute(
      Outcome = outcome_labels[i],
      OR = round(estimate, 3),
      CI_low = round(`2.5 %`, 3),
      CI_high = round(`97.5 %`, 3),
      p_value = round(p.value, 3)
    )
}))

print(weighted_adjusted_pax_effect)

# Save weighted adjusted model results ------------------------------------------------------------

weighted_adjusted_full_results <- dplyr::bind_rows(lapply(seq_along(weighted_adjusted_results), function(i) {
  
  weighted_adjusted_results[[i]]$summary %>%
    dplyr::mutate(
      Outcome = outcome_labels[i],
      Method = "Weighted adjusted model"
    )
}))

write.csv(
  weighted_adjusted_pax_effect,
  "weighted_adjusted_model_paxlovid_effect.csv",
  row.names = FALSE
)

write.csv(
  weighted_adjusted_full_results,
  "weighted_adjusted_model_full_results.csv",
  row.names = FALSE
)

# Models without IPTW and without doubly robust estimation ------------------------------------------------------------

library(mice)
library(dplyr)
library(openxlsx)

# Covariates only
# general_treatment, dusoi, covid_wave, and baseline_coop excluded
covariates <- c(
  "age",
  "sex",
  "Number_comorb",
  "Vaccination_number"
)

outcomes <- c(
  "new_improve_grorcoop",
  "new_improved_grandcoop",
  "improve_clinician",
  "improve_patient"
)

outcome_labels <- c(
  "Delta Grade or COOP",
  "Delta Grade and COOP",
  "Clinician-reported improvement",
  "Patient-reported improvement"
)

# Crude unweighted model - outcome ~ Paxlovid only ------------------------------------------------------------

fit_crude_model <- function(imp_list_model, outcome_var) {
  
  fit_list <- lapply(imp_list_model, function(dat) {
    
    model_formula <- as.formula(
      paste0(outcome_var, " ~ pax_num")
    )
    
    glm(
      model_formula,
      data = dat,
      family = binomial()
    )
  })
  
  pooled <- mice::pool(mice::as.mira(fit_list))
  
  pooled_summary <- summary(
    pooled,
    conf.int = TRUE,
    exponentiate = TRUE
  )
  
  list(
    fits = fit_list,
    pooled = pooled,
    summary = pooled_summary
  )
}

crude_results <- lapply(outcomes, function(outcome) {
  fit_crude_model(
    imp_list_model = imp_list_model,
    outcome_var = outcome
  )
})

names(crude_results) <- outcome_labels

# Adjusted unweighted model - outcome ~ Paxlovid + covariates ------------------------------------------------------------

adjusted_rhs <- paste(
  c("pax_num", covariates),
  collapse = " + "
)

fit_adjusted_model <- function(imp_list_model, outcome_var) {
  
  fit_list <- lapply(imp_list_model, function(dat) {
    
    model_formula <- as.formula(
      paste(outcome_var, "~", adjusted_rhs)
    )
    
    glm(
      model_formula,
      data = dat,
      family = binomial()
    )
  })
  
  pooled <- mice::pool(mice::as.mira(fit_list))
  
  pooled_summary <- summary(
    pooled,
    conf.int = TRUE,
    exponentiate = TRUE
  )
  
  list(
    fits = fit_list,
    pooled = pooled,
    summary = pooled_summary
  )
}

adjusted_results <- lapply(outcomes, function(outcome) {
  fit_adjusted_model(
    imp_list_model = imp_list_model,
    outcome_var = outcome
  )
})

names(adjusted_results) <- outcome_labels

# Print full results ------------------------------------------------------------

for (i in seq_along(outcomes)) {
  
  cat("\n====================================================\n")
  cat("Crude unweighted model:", outcome_labels[i], "\n")
  cat("====================================================\n")
  print(crude_results[[i]]$summary)
  
  cat("\n====================================================\n")
  cat("Adjusted unweighted model:", outcome_labels[i], "\n")
  cat("====================================================\n")
  print(adjusted_results[[i]]$summary)
}

# Extract Paxlovid effect only ------------------------------------------------------------

extract_pax_effect <- function(results_list, labels, method_name) {
  
  bind_rows(lapply(seq_along(results_list), function(i) {
    
    results_list[[i]]$summary %>%
      filter(term == "pax_num") %>%
      transmute(
        Outcome = labels[i],
        Method = method_name,
        OR = round(estimate, 3),
        CI_low = round(`2.5 %`, 3),
        CI_high = round(`97.5 %`, 3),
        p_value = round(p.value, 3)
      )
  }))
}

crude_pax_effect <- extract_pax_effect(
  crude_results,
  outcome_labels,
  "Crude unweighted"
)

adjusted_pax_effect <- extract_pax_effect(
  adjusted_results,
  outcome_labels,
  "Adjusted unweighted"
)

pax_effect_comparison <- bind_rows(
  crude_pax_effect,
  adjusted_pax_effect
)

print(pax_effect_comparison)

# Full model summaries ------------------------------------------------------------

save_model_summary <- function(results_list, labels, method_name) {
  
  bind_rows(lapply(seq_along(results_list), function(i) {
    
    results_list[[i]]$summary %>%
      mutate(
        Outcome = labels[i],
        Method = method_name
      )
  }))
}

full_crude <- save_model_summary(
  crude_results,
  outcome_labels,
  "Crude unweighted"
)

full_adjusted <- save_model_summary(
  adjusted_results,
  outcome_labels,
  "Adjusted unweighted"
)

full_all_no_weights <- bind_rows(
  full_crude,
  full_adjusted
)

# Save results ------------------------------------------------------------

write.csv(
  pax_effect_comparison,
  "no_weighting_paxlovid_effect_comparison.csv",
  row.names = FALSE
)

write.csv(
  full_all_no_weights,
  "no_weighting_full_model_results.csv",
  row.names = FALSE
)

wb <- createWorkbook()

addWorksheet(wb, "Paxlovid Effect")
writeData(wb, "Paxlovid Effect", pax_effect_comparison)

addWorksheet(wb, "Full Results")
writeData(wb, "Full Results", full_all_no_weights)

saveWorkbook(
  wb,
  "no_weighting_models_results.xlsx",
  overwrite = TRUE
)

# Unweighted crude models ------------------------------------------------------
# Crude model: outcome ~ paxlovid_treatment.
# No covariates, IPTW, or doubly robust adjustment.

# Load data ------------------------------------------------------------

data_raw <- read_excel("For Imputation.xlsx", sheet = "for imputation")

outcomes <- c(
  "new_improve_grorcoop",
  "new_improved_grandcoop",
  "improve_clinician",
  "improve_patient"
)

outcome_labels <- c(
  "Delta Grade or COOP",
  "Delta Grade and COOP",
  "Clinician-reported improvement",
  "Patient-reported improvement"
)

treatment_var <- "paxlovid_treatment"

vars_needed <- c(outcomes, treatment_var)

missing_vars <- setdiff(vars_needed, names(data_raw))

if (length(missing_vars) > 0) {
  stop(
    "These variables are missing from the Excel file: ",
    paste(missing_vars, collapse = ", ")
  )
}

df <- data_raw[, vars_needed, drop = FALSE]
# Convert variables to 0/1 numeric ------------------------------------------------------------

make_binary_numeric <- function(x) {
  x_chr <- tolower(trimws(as.character(x)))
  
  case_when(
    x_chr %in% c("1", "yes", "y", "true", "treated", "paxlovid") ~ 1,
    x_chr %in% c("0", "no", "n", "false", "untreated", "control", "no paxlovid") ~ 0,
    is.na(x_chr) ~ NA_real_,
    TRUE ~ suppressWarnings(as.numeric(x_chr))
  )
}

df$pax_num <- make_binary_numeric(df[[treatment_var]])

for (outcome in outcomes) {
  df[[paste0(outcome, "_num")]] <- make_binary_numeric(df[[outcome]])
}

cat("\nTreatment coding check:\n")
print(table(df[[treatment_var]], df$pax_num, useNA = "ifany"))

cat("\nOutcome coding checks:\n")
for (outcome in outcomes) {
  cat("\n", outcome, "\n")
  print(table(df[[outcome]], df[[paste0(outcome, "_num")]], useNA = "ifany"))
}

# Keep numeric treatment and numeric outcomes ------------------------------------------------------------

model_vars <- c(
  "pax_num",
  paste0(outcomes, "_num")
)

df_model <- dplyr::select(df, dplyr::all_of(model_vars))

# Multiple imputation ------------------------------------------------------------

init <- mice(df_model, maxit = 0, printFlag = FALSE)
meth <- init$method
pred <- init$predictorMatrix

meth[] <- "logreg"
meth["pax_num"] <- "logreg"

diag(pred) <- 0

imp <- mice(
  df_model,
  m = 20,
  maxit = 20,
  method = meth,
  predictorMatrix = pred,
  seed = 123,
  printFlag = TRUE
)

imp_list_model <- complete(imp, action = "all")

# Fit crude models - outcome ~ pax_num ------------------------------------------------------------

fit_crude_model <- function(imp_list_model, outcome_var_num) {
  
  fit_list <- lapply(imp_list_model, function(dat) {
    glm(
      as.formula(paste0(outcome_var_num, " ~ pax_num")),
      data = dat,
      family = binomial()
    )
  })
  
  pooled <- pool(as.mira(fit_list))
  
  summary(
    pooled,
    conf.int = TRUE,
    exponentiate = TRUE
  )
}

crude_results <- lapply(paste0(outcomes, "_num"), function(outcome_num) {
  fit_crude_model(imp_list_model, outcome_num)
})

names(crude_results) <- outcome_labels

# Print full crude model results ------------------------------------------------------------

for (i in seq_along(crude_results)) {
  cat("\n====================================================\n")
  cat("Unweighted unadjusted crude model:", outcome_labels[i], "\n")
  cat("Model:", paste0(outcomes[i], "_num"), "~ pax_num\n")
  cat("====================================================\n")
  print(crude_results[[i]])
}

# Combine full crude results ------------------------------------------------------------

full_crude_results <- bind_rows(lapply(seq_along(crude_results), function(i) {
  crude_results[[i]] %>%
    mutate(
      Outcome = outcome_labels[i],
      Outcome_variable = outcomes[i],
      Method = "Unweighted unadjusted",
      Model = paste0(outcomes[i], " ~ paxlovid_treatment")
    ) %>%
    select(
      Outcome,
      Outcome_variable,
      Method,
      Model,
      term,
      estimate,
      `2.5 %`,
      `97.5 %`,
      p.value,
      everything()
    )
}))

print(full_crude_results)

# Extract Paxlovid effect only ------------------------------------------------------------

crude_pax_effect <- full_crude_results %>%
  filter(term == "pax_num") %>%
  transmute(
    Outcome,
    Method,
    OR = round(estimate, 3),
    CI = paste0(
      "(",
      round(`2.5 %`, 3),
      ", ",
      round(`97.5 %`, 3),
      ")"
    ),
    p_value = round(p.value, 3)
  )

print(crude_pax_effect)

# Save results ------------------------------------------------------------

write.csv(
  full_crude_results,
  "unweighted_unadjusted_full_model_results.csv",
  row.names = FALSE
)

write.csv(
  crude_pax_effect,
  "unweighted_unadjusted_paxlovid_effect.csv",
  row.names = FALSE
)

wb <- createWorkbook()

addWorksheet(wb, "Full Crude Results")
writeData(wb, "Full Crude Results", full_crude_results)

addWorksheet(wb, "Paxlovid Effect")
writeData(wb, "Paxlovid Effect", crude_pax_effect)

saveWorkbook(
  wb,
  "unweighted_unadjusted_model_results.xlsx",
  overwrite = TRUE
)

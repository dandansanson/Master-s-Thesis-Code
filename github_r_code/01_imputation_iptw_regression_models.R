# Multiple-imputation, IPTW, and regression analyses
# Long COVID/Paxlovid thesis analysis
# Input data files should be in the working directory.

# Required packages ----
required_packages <- c("readxl","mice","dplyr", "tidyr","cobalt","broom","openxlsx","ggplot2","naniar","VIM","WeightIt")

missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) {
  stop(
    "Install the following packages before running this script: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

invisible(lapply(required_packages, library, character.only = TRUE))

# Settings ----
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
covariates <- c("age","sex","Number_comorb","Vaccination_number")

vars_needed <- c(outcomes, treatment_var, covariates)

# Load data ----
data_raw <- readxl::read_excel(file_name, sheet = sheet_name)

missing_vars <- setdiff(vars_needed, names(data_raw))

if (length(missing_vars) > 0) {
  stop("These variables are missing from the Excel file: ",paste(missing_vars, collapse = ", "))
}

lc_pat <- data_raw %>%
  dplyr::select(dplyr::all_of(vars_needed))

cat("\ndata dimensions:\n")
print(dim(lc_pat))

cat("\nVariable names:\n")
print(names(lc_pat))

# Check missingness pattern first ----
missing_table <- data.frame(
  Variable = names(lc_pat),
  Missing_n = colSums(is.na(lc_pat)),
  Missing_pct = round(colMeans(is.na(lc_pat)) * 100, 2)
) %>%
  dplyr::arrange(desc(Missing_pct))

cat("\nMissingness summary:\n")
print(missing_table)

missing_pattern <- mice::md.pattern(lc_pat, plot = FALSE)

cat("\nMissingness pattern:\n")
print(missing_pattern)

if (anyNA(lc_pat)) {
  VIM::aggr(lc_pat,
    col = c("navyblue", "red"),
    numbers = TRUE,
    sortVars = TRUE,
    labels = names(lc_pat),
    cex.axis = 0.7,
    gap = 3,
    ylab = c("Missing data", "Pattern")
  )

  print(naniar::vis_miss(lc_pat))
}

# Recode variable types ----
factor_vars <- c( outcomes,"paxlovid_treatment","sex")

numeric_vars <- c("age","Number_comorb","Vaccination_number")

lc_pat <- lc_pat %>%
  dplyr::mutate(
    dplyr::across(dplyr::all_of(factor_vars), as.factor),
    dplyr::across(dplyr::all_of(numeric_vars), ~ as.numeric(.))
  )


# Multiple imputation setup ----
choose_factor_method <- function(x) {
  if (!any(is.na(x))) return("")

  nlev <- nlevels(droplevels(x))

  if (nlev == 2) return("logreg")
  if (nlev > 2) return("polyreg")

  return("")
}

if (anyNA(lc_pat)) {
  init <- mice::mice(lc_pat, maxit = 0, printFlag = FALSE)
  meth <- init$method
  pred <- init$predictorMatrix

  for (v in factor_vars) {
    meth[v] <- choose_factor_method(lc_pat[[v]])
  }

  for (v in numeric_vars) {
    meth[v] <- ifelse(any(is.na(lc_pat[[v]])), "pmm", "")
  }

  diag(pred) <- 0

  cat("\nImputation methods:\n")
  print(meth)

  imp <- mice::mice(
    lc_pat,
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
  imp_list <- rep(list(lc_pat), m_imputations)
}

# Convert Paxlovid treatment to 0/1 ----
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

add_treatment_numeric <- function(lc_patient) {
  lc_patient$pax_num <- make_binary_numeric(lc_patient[[treatment_var]])

  bad_values <- unique(as.character(
    lc_patient[[treatment_var]][is.na(lc_patient$pax_num) & !is.na(lc_patient[[treatment_var]])]
  ))

  if (length(bad_values) > 0) {
    stop(
      "Some paxlovid_treatment values were not recognized: ",
      paste(bad_values, collapse = ", "),
      "\nEdit make_binary_numeric() to match your coding."
    )
  }

  lc_patient
}

imp_list_model <- lapply(imp_list, add_treatment_numeric)

cat("\nTreatment coding check, first imputed dataset:\n")
print(table(
  imp_list_model[[1]][[treatment_var]],
  imp_list_model[[1]]$pax_num,
  useNA = "ifany"
))

# Helper functions for model fitting and pooling ----
covariate_rhs <- paste(covariates, collapse = " + ")

adjusted_rhs <- paste(
  c("pax_num", covariates),
  collapse = " + "
)

fit_model_set <- function(lc_patient_list, outcome_var, rhs, weights_col = NULL) {
  fit_list <- lapply(lc_patient_list, function(lc_patient) {
    model_formula <- as.formula(
      paste(outcome_var, "~", rhs)
    )

    if (is.null(weights_col)) {
      glm(
        model_formula,
        data = lc_patient,
        family = binomial()
      )
    } else {
      glm(
        model_formula,
        data = lc_patient,
        weights = lc_patient[[weights_col]],
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

# Model without weighting - Treatment + covariates, no IPTW ----
unweighted_results <- lapply(outcomes, function(outcome) {
  fit_model_set(
    lc_patient_list = imp_list_model,
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

# Check covariate balance before IPTW ----
balance_formula <- as.formula(
  paste("pax_num ~", covariate_rhs)
)

balance_to_lc_pat <- function(bal_obj) {
  out <- as.data.frame(bal_obj$Balance)
  out$Covariate <- rownames(out)
  rownames(out) <- NULL
  out %>% dplyr::relocate(Covariate)
}

get_balance_lc_pat <- function(lc_patient, weights = NULL, phase, imputation) {
  if (is.null(weights)) {
    bal <- cobalt::bal.tab(
      balance_formula,
      data = lc_patient,
      binary = "std",
      continuous = "std",
      s.d.denom = "pooled"
    )
  } else {
    bal <- cobalt::bal.tab(
      balance_formula,
      data = lc_patient,
      weights = weights,
      method = "weighting",
      un = TRUE,
      binary = "std",
      continuous = "std",
      s.d.denom = "pooled"
    )
  }

  balance_to_lc_pat(bal) %>%
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
  get_balance_lc_pat(
    lc_patient = imp_list_model[[i]],
    weights = NULL,
    phase = "Before IPTW",
    imputation = i
  )
}))

# Estimate propensity scores and IPTW ----
ps_formula <- as.formula(
  paste("pax_num ~", covariate_rhs)
)

add_ps_weights <- function(lc_patient) {
  ps_model <- glm(
    ps_formula,
    data = lc_patient,
    family = binomial()
  )

  lc_patient$ps <- predict(ps_model, newdata = lc_patient, type = "response")

  lc_patient$ps <- pmin(pmax(lc_patient$ps, 0.01), 0.99)

  p_treat <- mean(lc_patient$pax_num == 1, na.rm = TRUE)

  lc_patient$sw <- ifelse(
    lc_patient$pax_num == 1,
    p_treat / lc_patient$ps,
    (1 - p_treat) / (1 - lc_patient$ps)
  )

  cap <- quantile(lc_patient$sw, 0.99, na.rm = TRUE)
  lc_patient$sw_trunc <- pmin(lc_patient$sw, as.numeric(cap))

  list(
    data = lc_patient,
    ps_model = ps_model
  )
}

ps_list <- lapply(imp_list_model, add_ps_weights)

imp_list_w <- lapply(ps_list, `[[`, "data")
ps_models <- lapply(ps_list, `[[`, "ps_model")

weight_summary <- dplyr::bind_rows(lapply(seq_along(imp_list_w), function(i) {
  lc_patient <- imp_list_w[[i]]

  data.frame(
    Imputation = i,
    n = nrow(lc_patient),
    treated_n = sum(lc_patient$pax_num == 1, na.rm = TRUE),
    untreated_n = sum(lc_patient$pax_num == 0, na.rm = TRUE),
    sw_min = min(lc_patient$sw, na.rm = TRUE),
    sw_mean = mean(lc_patient$sw, na.rm = TRUE),
    sw_median = median(lc_patient$sw, na.rm = TRUE),
    sw_p95 = as.numeric(quantile(lc_patient$sw, 0.95, na.rm = TRUE)),
    sw_p99 = as.numeric(quantile(lc_patient$sw, 0.99, na.rm = TRUE)),
    sw_max = max(lc_patient$sw, na.rm = TRUE),
    sw_trunc_max = max(lc_patient$sw_trunc, na.rm = TRUE)
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

# IPTW model without doubly robust estimator - Outcome model = treatment only, weighted by IPTW ----
iptw_results <- lapply(outcomes, function(outcome) {
  fit_model_set(
    lc_patient_list = imp_list_w,
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

# Check covariate balance after IPTW ----
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
  get_balance_lc_pat(
    lc_patient = imp_list_w[[i]],
    weights = imp_list_w[[i]]$sw_trunc,
    phase = "After IPTW",
    imputation = i
  )
}))

balance_all <- dplyr::bind_rows(
  pre_balance_all,
  post_balance_all
)

# Doubly robust estimator - Outcome model = treatment + covariates, weighted by IPTW ----
dr_results <- lapply(outcomes, function(outcome) {
  fit_model_set(
    lc_patient_list = imp_list_w,
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

# Compare Paxlovid effect across methods ----
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

# Save results ----
missing_pattern_lc_pat <- as.data.frame(missing_pattern)
missing_pattern_lc_pat$Pattern <- rownames(missing_pattern_lc_pat)
rownames(missing_pattern_lc_pat) <- NULL

missing_pattern_lc_pat <- missing_pattern_lc_pat %>%
  dplyr::relocate(Pattern)

full_all <- dplyr::bind_rows(
  full_unweighted,
  full_iptw,
  full_dr
)
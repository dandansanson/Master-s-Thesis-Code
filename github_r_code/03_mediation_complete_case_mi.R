# Complete-case and multiple-imputation mediation analyses
# File: 03_mediation_complete_case_mi.R
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
  "dplyr",
  "tidyr",
  "mice",
  "mediation",
  "openxlsx"
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

m_imp <- 20
maxit_imp <- 20
n_sims <- 1000
use_boot <- TRUE

mediators <- c(
  "VL",
  "PDZ1K1",
  "KALRN",
  "CETP"
)

covariates <- c(
  "age",
  "sex",
  "Number_comorb",
  "Vaccination_number"
)

outcomes <- c(
  "outcome_delta_or_coop",
  "outcome_delta_and_coop",
  "outcome_clinician",
  "outcome_patient"
)

outcome_labels <- c(
  "Patient or Clinician reported improvement",
  "Patient and Clinician reported improvement",
  "Clinician-reported improvement",
  "Patient-reported improvement"
)

# Load data ------------------------------------------------------------

data_raw <- readxl::read_excel("for mediation analysis with VL.xlsx")

vars_needed <- c(
  "Treatment",
  mediators,
  "new_improve_grorcoop",
  "new_improved_grandcoop",
  "improve_clinician",
  "improve_patient",
  covariates
)

missing_vars <- setdiff(vars_needed, names(data_raw))

if (length(missing_vars) > 0) {
  stop("Missing variables: ", paste(missing_vars, collapse = ", "))
}

df <- data_raw[, vars_needed, drop = FALSE]

# Recode variables ------------------------------------------------------------

make_binary_numeric <- function(x) {
  x_chr <- trimws(tolower(as.character(x)))
  
  dplyr::case_when(
    x_chr %in% c("1", "yes", "y", "improved", "treated", "true") ~ 1,
    x_chr %in% c("0", "no", "n", "not improved", "untreated", "false") ~ 0,
    is.na(x) ~ NA_real_,
    TRUE ~ NA_real_
  )
}

recode_treatment <- function(x) {
  x_num <- suppressWarnings(as.numeric(as.character(x)))
  x_chr <- trimws(tolower(as.character(x)))
  
  dplyr::case_when(
    x_num == 2 | x_chr %in% c("paxlovid") ~ 1,
    x_num == 0 | x_chr %in% c("no treatment", "untreated", "control", "none") ~ 0,
    x_num == 1 | x_chr %in% c("general", "general treatment") ~ NA_real_,
    is.na(x) ~ NA_real_,
    TRUE ~ NA_real_
  )
}

df <- df %>%
  dplyr::mutate(
    treatment_num = recode_treatment(Treatment),
    
    outcome_delta_or_coop = make_binary_numeric(new_improve_grorcoop),
    outcome_delta_and_coop = make_binary_numeric(new_improved_grandcoop),
    outcome_clinician = make_binary_numeric(improve_clinician),
    outcome_patient = make_binary_numeric(improve_patient),
    
    age = as.numeric(age),
    Number_comorb = as.numeric(Number_comorb),
    Vaccination_number = as.numeric(Vaccination_number),
    sex = as.factor(sex)
  )

for (m in mediators) {
  df[[m]] <- as.numeric(df[[m]])
}

# Exclude general treatment and missing treatment
df <- df[!is.na(df$treatment_num), , drop = FALSE]

cat("\nTreatment distribution after excluding general/missing treatment:\n")
print(table(df$treatment_num, useNA = "ifany"))

# Helper functions ------------------------------------------------------------

prepare_analysis_data <- function(dat) {
  
  binary_vars <- c(
    "treatment_num",
    "outcome_delta_or_coop",
    "outcome_delta_and_coop",
    "outcome_clinician",
    "outcome_patient"
  )
  
  for (v in binary_vars) {
    dat[[v]] <- as.numeric(as.character(dat[[v]]))
  }
  
  dat$sex <- as.factor(dat$sex)
  
  dat
}

extract_mediation <- function(med_result, mediator_name, outcome_name, analysis_name, imputation = NA) {
  
  s <- summary(med_result)
  
  data.frame(
    Mediator = mediator_name,
    Analysis = analysis_name,
    Imputation = imputation,
    Outcome = outcome_name,
    
    Effect = c("ACME", "ADE", "Total effect", "Proportion mediated"),
    
    Estimate = c(
      s$d.avg,
      s$z.avg,
      s$tau.coef,
      s$n.avg
    ),
    
    SE = c(
      sd(med_result$d.avg.sims, na.rm = TRUE),
      sd(med_result$z.avg.sims, na.rm = TRUE),
      sd(med_result$tau.sims, na.rm = TRUE),
      sd(med_result$n.avg.sims, na.rm = TRUE)
    ),
    
    CI_low = c(
      s$d.avg.ci[1],
      s$z.avg.ci[1],
      s$tau.ci[1],
      s$n.avg.ci[1]
    ),
    
    CI_high = c(
      s$d.avg.ci[2],
      s$z.avg.ci[2],
      s$tau.ci[2],
      s$n.avg.ci[2]
    ),
    
    p_value = c(
      s$d.avg.p,
      s$z.avg.p,
      s$tau.p,
      s$n.avg.p
    )
  )
}

run_mediation_one <- function(dat, mediator_name, outcome_var, outcome_name, analysis_name, imputation = NA) {
  
  dat <- prepare_analysis_data(dat)
  dat$Y <- dat[[outcome_var]]
  
  med_formula <- as.formula(
    paste(mediator_name, "~ treatment_num +", paste(covariates, collapse = " + "))
  )
  
  out_formula <- as.formula(
    paste("Y ~ treatment_num +", mediator_name, "+", paste(covariates, collapse = " + "))
  )
  
  med_model <- lm(
    med_formula,
    data = dat
  )
  
  out_model <- glm(
    out_formula,
    data = dat,
    family = binomial(link = "logit")
  )
  
  med_model$call$formula <- med_formula
  out_model$call$formula <- out_formula
  
  med_result <- mediation::mediate(
    model.m = med_model,
    model.y = out_model,
    treat = "treatment_num",
    mediator = mediator_name,
    treat.value = 1,
    control.value = 0,
    boot = use_boot,
    sims = n_sims
  )
  
  extract_mediation(
    med_result = med_result,
    mediator_name = mediator_name,
    outcome_name = outcome_name,
    analysis_name = analysis_name,
    imputation = imputation
  )
}

safe_run_mediation_one <- function(dat, mediator_name, outcome_var, outcome_name, analysis_name, imputation = NA) {
  
  tryCatch(
    {
      run_mediation_one(
        dat = dat,
        mediator_name = mediator_name,
        outcome_var = outcome_var,
        outcome_name = outcome_name,
        analysis_name = analysis_name,
        imputation = imputation
      )
    },
    error = function(e) {
      data.frame(
        Mediator = mediator_name,
        Analysis = analysis_name,
        Imputation = imputation,
        Outcome = outcome_name,
        Effect = c("ACME", "ADE", "Total effect", "Proportion mediated"),
        Estimate = NA_real_,
        SE = NA_real_,
        CI_low = NA_real_,
        CI_high = NA_real_,
        p_value = NA_real_,
        Error = e$message
      )
    }
  )
}

pool_rubin <- function(est, se) {
  
  ok <- !is.na(est) & !is.na(se)
  est <- est[ok]
  se <- se[ok]
  
  m <- length(est)
  
  if (m == 0) {
    return(data.frame(
      Estimate = NA_real_,
      SE = NA_real_,
      CI_low = NA_real_,
      CI_high = NA_real_,
      p_value = NA_real_,
      N_imputations_used = 0
    ))
  }
  
  q_bar <- mean(est)
  u_bar <- mean(se^2)
  b <- stats::var(est)
  
  if (m == 1 || is.na(b)) b <- 0
  
  total_var <- u_bar + (1 + 1 / m) * b
  total_se <- sqrt(total_var)
  
  if (b == 0 || u_bar == 0) {
    df <- Inf
  } else {
    r <- ((1 + 1 / m) * b) / u_bar
    df <- (m - 1) * (1 + 1 / r)^2
  }
  
  crit <- ifelse(is.infinite(df), 1.96, qt(0.975, df = df))
  p_value <- 2 * pt(abs(q_bar / total_se), df = df, lower.tail = FALSE)
  
  data.frame(
    Estimate = q_bar,
    SE = total_se,
    CI_low = q_bar - crit * total_se,
    CI_high = q_bar + crit * total_se,
    p_value = p_value,
    N_imputations_used = m
  )
}

# Run complete-case and MI mediation for each mediator ------------------------------------------------------------

all_results <- list()

for (mediator_var in mediators) {
  
  cat("\n========================================\n")
  cat("Running mediation for mediator:", mediator_var, "\n")
  cat("========================================\n")
  
  vars_for_analysis <- c(
    "treatment_num",
    mediator_var,
    outcomes,
    covariates
  )
  
  df_analysis <- df[, vars_for_analysis, drop = FALSE]
  
  # Missingness table for this mediator
  missing_table <- data.frame(
    Mediator = mediator_var,
    Variable = names(df_analysis),
    Missing_n = colSums(is.na(df_analysis)),
    Missing_pct = round(colMeans(is.na(df_analysis)) * 100, 2)
  )
  
  # Complete-case analysis
  dat_cc <- df_analysis[stats::complete.cases(df_analysis), , drop = FALSE]
  dat_cc <- prepare_analysis_data(dat_cc)
  
  cat("\nComplete-case sample size for", mediator_var, ":", nrow(dat_cc), "\n")
  
  cc_results <- dplyr::bind_rows(lapply(seq_along(outcomes), function(i) {
    safe_run_mediation_one(
      dat = dat_cc,
      mediator_name = mediator_var,
      outcome_var = outcomes[i],
      outcome_name = outcome_labels[i],
      analysis_name = "Complete case",
      imputation = NA
    )
  }))
  
  # Multiple imputation
  df_imp <- df_analysis
  
  binary_vars <- c(
    "treatment_num",
    "outcome_delta_or_coop",
    "outcome_delta_and_coop",
    "outcome_clinician",
    "outcome_patient"
  )
  
  numeric_vars <- c(
    mediator_var,
    "age",
    "Number_comorb",
    "Vaccination_number"
  )
  
  for (v in binary_vars) {
    df_imp[[v]] <- factor(df_imp[[v]], levels = c(0, 1))
  }
  
  df_imp$sex <- as.factor(df_imp$sex)
  
  init <- mice::mice(
    df_imp,
    maxit = 0,
    printFlag = FALSE
  )
  
  meth <- init$method
  pred <- init$predictorMatrix
  
  meth["treatment_num"] <- ""
  
  for (v in setdiff(binary_vars, "treatment_num")) {
    meth[v] <- ifelse(any(is.na(df_imp[[v]])), "logreg", "")
  }
  
  if (any(is.na(df_imp$sex))) {
    meth["sex"] <- ifelse(nlevels(droplevels(df_imp$sex)) == 2, "logreg", "polyreg")
  } else {
    meth["sex"] <- ""
  }
  
  for (v in numeric_vars) {
    meth[v] <- ifelse(any(is.na(df_imp[[v]])), "pmm", "")
  }
  
  diag(pred) <- 0
  
  imp <- mice::mice(
    df_imp,
    m = m_imp,
    maxit = maxit_imp,
    method = meth,
    predictorMatrix = pred,
    seed = 123,
    printFlag = TRUE
  )
  
  imp_list <- mice::complete(imp, action = "all")
  
  mi_results_each <- dplyr::bind_rows(lapply(seq_along(imp_list), function(j) {
    
    dat_j <- prepare_analysis_data(imp_list[[j]])
    
    dplyr::bind_rows(lapply(seq_along(outcomes), function(i) {
      safe_run_mediation_one(
        dat = dat_j,
        mediator_name = mediator_var,
        outcome_var = outcomes[i],
        outcome_name = outcome_labels[i],
        analysis_name = "Multiple imputation",
        imputation = j
      )
    }))
  }))
  
  mi_pooled <- mi_results_each %>%
    dplyr::group_by(Mediator, Outcome, Effect) %>%
    dplyr::group_modify(~ pool_rubin(.x$Estimate, .x$SE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      Analysis = "Multiple imputation pooled",
      Imputation = NA
    ) %>%
    dplyr::select(
      Mediator,
      Analysis,
      Imputation,
      Outcome,
      Effect,
      Estimate,
      SE,
      CI_low,
      CI_high,
      p_value,
      N_imputations_used
    )
  
  cc_results_for_compare <- cc_results %>%
    dplyr::mutate(N_imputations_used = NA_integer_) %>%
    dplyr::select(
      Mediator,
      Analysis,
      Imputation,
      Outcome,
      Effect,
      Estimate,
      SE,
      CI_low,
      CI_high,
      p_value,
      N_imputations_used
    )
  
  comparison_long <- dplyr::bind_rows(
    cc_results_for_compare,
    mi_pooled
  ) %>%
    dplyr::mutate(
      Estimate = round(Estimate, 3),
      SE = round(SE, 3),
      CI_low = round(CI_low, 3),
      CI_high = round(CI_high, 3),
      p_value = round(p_value, 3),
      CI = paste0("(", CI_low, ", ", CI_high, ")")
    )
  
  comparison_wide <- comparison_long %>%
    dplyr::select(Mediator, Outcome, Effect, Analysis, Estimate, CI, p_value) %>%
    tidyr::pivot_wider(
      names_from = Analysis,
      values_from = c(Estimate, CI, p_value)
    )
  
  all_results[[mediator_var]] <- list(
    missing_table = missing_table,
    cc_results = cc_results,
    mi_results_each = mi_results_each,
    mi_pooled = mi_pooled,
    comparison_long = comparison_long,
    comparison_wide = comparison_wide
  )
}

# Combine and save all mediator results ------------------------------------------------------------

missing_all <- dplyr::bind_rows(lapply(all_results, `[[`, "missing_table"))
cc_all <- dplyr::bind_rows(lapply(all_results, `[[`, "cc_results"))
mi_each_all <- dplyr::bind_rows(lapply(all_results, `[[`, "mi_results_each"))
mi_pooled_all <- dplyr::bind_rows(lapply(all_results, `[[`, "mi_pooled"))
comparison_long_all <- dplyr::bind_rows(lapply(all_results, `[[`, "comparison_long"))
comparison_wide_all <- dplyr::bind_rows(lapply(all_results, `[[`, "comparison_wide"))

write.csv(missing_all, "all_mediators_missingness.csv", row.names = FALSE)
write.csv(cc_all, "all_mediators_complete_case_results.csv", row.names = FALSE)
write.csv(mi_each_all, "all_mediators_each_imputed_dataset_results.csv", row.names = FALSE)
write.csv(mi_pooled_all, "all_mediators_MI_pooled_results.csv", row.names = FALSE)
write.csv(comparison_long_all, "all_mediators_CC_vs_MI_long.csv", row.names = FALSE)
write.csv(comparison_wide_all, "all_mediators_CC_vs_MI_wide.csv", row.names = FALSE)

wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb, "Missingness")
openxlsx::writeData(wb, "Missingness", missing_all)

openxlsx::addWorksheet(wb, "Complete Case")
openxlsx::writeData(wb, "Complete Case", cc_all)

openxlsx::addWorksheet(wb, "MI Each Dataset")
openxlsx::writeData(wb, "MI Each Dataset", mi_each_all)

openxlsx::addWorksheet(wb, "MI Pooled")
openxlsx::writeData(wb, "MI Pooled", mi_pooled_all)

openxlsx::addWorksheet(wb, "Comparison Long")
openxlsx::writeData(wb, "Comparison Long", comparison_long_all)

openxlsx::addWorksheet(wb, "Comparison Wide")
openxlsx::writeData(wb, "Comparison Wide", comparison_wide_all)

openxlsx::saveWorkbook(
  wb,
  "all_mediators_mediation_CC_vs_MI.xlsx",
  overwrite = TRUE
)

cat("\nAll mediation analyses complete.\n")

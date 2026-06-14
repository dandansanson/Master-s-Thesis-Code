
m_imp <- 20
maxit_imp <- 20
n_sims <- 1000
use_boot <- TRUE

mediators <- c("VL", "PDZ1K1", "KALRN", "CETP", "viral load", "orf1ab")

covariates <- c("age","sex","Number_comorb","Vaccination_number")

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

# Load data ----
data_raw <- readxl::read_excel("for mediation analysis.xlsx")

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

lc_pat <- data_raw[, vars_needed, drop = FALSE]

# Recode variables ----
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

lc_pat <- lc_pat %>%
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
  lc_pat[[m]] <- as.numeric(lc_pat[[m]])
}

# Exclude general treatment and missing treatment
lc_pat <- lc_pat[!is.na(lc_pat$treatment_num), , drop = FALSE]

cat("\nTreatment distribution after excluding general/missing treatment:\n")
print(table(lc_pat$treatment_num, useNA = "ifany"))

# Helper functions ----
prepare_analysis_data <- function(lc_patient) {
  binary_vars <- c(
    "treatment_num",
    "outcome_delta_or_coop",
    "outcome_delta_and_coop",
    "outcome_clinician",
    "outcome_patient"
  )

  for (v in binary_vars) {
    lc_patient[[v]] <- as.numeric(as.character(lc_patient[[v]]))
  }

  lc_patient$sex <- as.factor(lc_patient$sex)

  lc_patient
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

run_mediation_one <- function(lc_patient, mediator_name, outcome_var, outcome_name, analysis_name, imputation = NA) {
  lc_patient <- prepare_analysis_data(lc_patient)
  lc_patient$Y <- lc_patient[[outcome_var]]

  med_formula <- as.formula(
    paste(mediator_name, "~ treatment_num +", paste(covariates, collapse = " + "))
  )

  out_formula <- as.formula(
    paste("Y ~ treatment_num +", mediator_name, "+", paste(covariates, collapse = " + "))
  )

  med_model <- lm(
    med_formula,
    data = lc_patient
  )

  out_model <- glm(
    out_formula,
    data = lc_patient,
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

safe_run_mediation_one <- function(lc_patient, mediator_name, outcome_var, outcome_name, analysis_name, imputation = NA) {
  tryCatch(
    {
      run_mediation_one(
        lc_patient = lc_patient,
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
    lc_pat <- Inf
  } else {
    r <- ((1 + 1 / m) * b) / u_bar
    lc_pat <- (m - 1) * (1 + 1 / r)^2
  }

  crit <- ifelse(is.infinite(lc_pat), 1.96, qt(0.975, lc_pat = lc_pat))
  p_value <- 2 * pt(abs(q_bar / total_se), lc_pat = lc_pat, lower.tail = FALSE)

  data.frame(
    Estimate = q_bar,
    SE = total_se,
    CI_low = q_bar - crit * total_se,
    CI_high = q_bar + crit * total_se,
    p_value = p_value,
    N_imputations_used = m
  )
}

# Run complete-case and MI mediation for each mediator ----
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

  lc_pat_analysis <- lc_pat[, vars_for_analysis, drop = FALSE]

  # Missingness table for this mediator
  missing_table <- data.frame(
    Mediator = mediator_var,
    Variable = names(lc_pat_analysis),
    Missing_n = colSums(is.na(lc_pat_analysis)),
    Missing_pct = round(colMeans(is.na(lc_pat_analysis)) * 100, 2)
  )

  # Complete-case analysis
  lc_patient_cc <- lc_pat_analysis[stats::complete.cases(lc_pat_analysis), , drop = FALSE]
  lc_patient_cc <- prepare_analysis_data(lc_patient_cc)

  cat("\nComplete-case sample size for", mediator_var, ":", nrow(lc_patient_cc), "\n")

  cc_results <- dplyr::bind_rows(lapply(seq_along(outcomes), function(i) {
    safe_run_mediation_one(
      lc_patient = lc_patient_cc,
      mediator_name = mediator_var,
      outcome_var = outcomes[i],
      outcome_name = outcome_labels[i],
      analysis_name = "Complete case",
      imputation = NA
    )
  }))

  # Multiple imputation
  lc_pat_imp <- lc_pat_analysis

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
    lc_pat_imp[[v]] <- factor(lc_pat_imp[[v]], levels = c(0, 1))
  }

  lc_pat_imp$sex <- as.factor(lc_pat_imp$sex)

  init <- mice::mice(
    lc_pat_imp,
    maxit = 0,
    printFlag = FALSE
  )

  meth <- init$method
  pred <- init$predictorMatrix

  meth["treatment_num"] <- ""

  for (v in setdiff(binary_vars, "treatment_num")) {
    meth[v] <- ifelse(any(is.na(lc_pat_imp[[v]])), "logreg", "")
  }

  if (any(is.na(lc_pat_imp$sex))) {
    meth["sex"] <- ifelse(nlevels(droplevels(lc_pat_imp$sex)) == 2, "logreg", "polyreg")
  } else {
    meth["sex"] <- ""
  }

  for (v in numeric_vars) {
    meth[v] <- ifelse(any(is.na(lc_pat_imp[[v]])), "pmm", "")
  }

  diag(pred) <- 0

  imp <- mice::mice(
    lc_pat_imp,
    m = m_imp,
    maxit = maxit_imp,
    method = meth,
    predictorMatrix = pred,
    seed = 123,
    printFlag = TRUE
  )

  imp_list <- mice::complete(imp, action = "all")

  mi_results_each <- dplyr::bind_rows(lapply(seq_along(imp_list), function(j) {
    lc_patient_j <- prepare_analysis_data(imp_list[[j]])

    dplyr::bind_rows(lapply(seq_along(outcomes), function(i) {
      safe_run_mediation_one(
        lc_patient = lc_patient_j,
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

# Combine and save all mediator results ----
missing_all <- dplyr::bind_rows(lapply(all_results, `[[`, "missing_table"))
cc_all <- dplyr::bind_rows(lapply(all_results, `[[`, "cc_results"))
mi_each_all <- dplyr::bind_rows(lapply(all_results, `[[`, "mi_results_each"))
mi_pooled_all <- dplyr::bind_rows(lapply(all_results, `[[`, "mi_pooled"))
comparison_long_all <- dplyr::bind_rows(lapply(all_results, `[[`, "comparison_long"))
comparison_wide_all <- dplyr::bind_rows(lapply(all_results, `[[`, "comparison_wide"))

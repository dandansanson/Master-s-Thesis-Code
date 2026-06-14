# Entropy balancing, AIPW, and missingness diagnostics
# File: 02_entropy_balancing_aipw_missingness.R
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
  "WeightIt",
  "ggplot2",
  "openxlsx",
  "broom"
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

# Load data ------------------------------------------------------------
data_raw <- readxl::read_excel(
  "For Imputation.xlsx",
  sheet = "for imputation"
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

treatment_var <- "paxlovid_treatment"

covariates <- c(
  "age",
  "sex",
  "Number_comorb",
  "Vaccination_number"
)

vars_needed <- c(outcomes, treatment_var, covariates)

missing_vars <- setdiff(vars_needed, names(data_raw))

if (length(missing_vars) > 0) {
  stop(
    "These variables are missing from the Excel file: ",
    paste(missing_vars, collapse = ", ")
  )
}

df <- data_raw[, vars_needed, drop = FALSE]

# Recode variables ------------------------------------------------------------
df$sex <- as.factor(df$sex)
df$paxlovid_treatment <- as.factor(df$paxlovid_treatment)

for (v in outcomes) {
  df[[v]] <- as.factor(df[[v]])
}

df$age <- as.numeric(df$age)
df$Number_comorb <- as.numeric(df$Number_comorb)
df$Vaccination_number <- as.numeric(df$Vaccination_number)

# Missingness check ------------------------------------------------------------
missing_table <- data.frame(
  Variable = names(df),
  Missing_n = colSums(is.na(df)),
  Missing_pct = round(colMeans(is.na(df)) * 100, 2)
)

missing_table <- missing_table[order(-missing_table$Missing_pct), ]

print(missing_table)
print(mice::md.pattern(df, plot = FALSE))

# Multiple imputation ------------------------------------------------------------
choose_factor_method <- function(x) {
  if (!any(is.na(x))) return("")
  nlev <- nlevels(droplevels(x))
  if (nlev == 2) return("logreg")
  if (nlev > 2) return("polyreg")
  return("")
}

factor_vars <- c(outcomes, "paxlovid_treatment", "sex")
numeric_vars <- c("age", "Number_comorb", "Vaccination_number")

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
  
  imp <- mice::mice(
    df,
    m = 20,
    maxit = 20,
    method = meth,
    predictorMatrix = pred,
    seed = 123,
    printFlag = TRUE
  )
  
  imp_list <- mice::complete(imp, action = "all")
  
} else {
  imp <- NULL
  imp_list <- list(df)
}

# Convert treatment and outcomes to 0/1 numeric ------------------------------------------------------------
make_binary_numeric <- function(x) {
  x_chr <- tolower(trimws(as.character(x)))
  x_num <- suppressWarnings(as.numeric(x_chr))
  
  dplyr::case_when(
    x_chr %in% c("1", "yes", "y", "true", "treated", "paxlovid", "improved") ~ 1,
    x_chr %in% c("0", "no", "n", "false", "untreated", "control", "no paxlovid", "not improved") ~ 0,
    !is.na(x_num) & x_num %in% c(0, 1) ~ x_num,
    is.na(x_chr) ~ NA_real_,
    TRUE ~ NA_real_
  )
}

add_numeric_vars <- function(dat) {
  dat$pax_num <- make_binary_numeric(dat[[treatment_var]])
  dat$sex <- as.factor(dat$sex)
  
  for (outcome in outcomes) {
    dat[[paste0(outcome, "_num")]] <- make_binary_numeric(dat[[outcome]])
  }
  
  dat
}

imp_list_model <- lapply(imp_list, add_numeric_vars)

cat("\nTreatment coding check:\n")
print(table(
  imp_list_model[[1]]$paxlovid_treatment,
  imp_list_model[[1]]$pax_num,
  useNA = "ifany"
))

cat("\nOutcome coding checks:\n")
for (outcome in outcomes) {
  cat("\n", outcome, "\n")
  print(table(
    imp_list_model[[1]][[outcome]],
    imp_list_model[[1]][[paste0(outcome, "_num")]],
    useNA = "ifany"
  ))
}

# Entropy balancing ------------------------------------------------------------
balance_formula <- pax_num ~ age + sex + Number_comorb + Vaccination_number

add_entropy_weights <- function(dat) {
  
  dat$sex <- as.factor(dat$sex)
  
  w_obj <- WeightIt::weightit(
    balance_formula,
    data = dat,
    method = "ebal",
    estimand = "ATE"
  )
  
  dat$sw <- w_obj$weights
  
  list(
    data = dat,
    weightit_object = w_obj
  )
}

weight_list <- lapply(imp_list_model, add_entropy_weights)

imp_list_w <- lapply(weight_list, `[[`, "data")
weightit_objects <- lapply(weight_list, `[[`, "weightit_object")

# Covariate balance ------------------------------------------------------------
bal_after <- cobalt::bal.tab(
  balance_formula,
  data = imp_list_w[[1]],
  weights = imp_list_w[[1]]$sw,
  method = "weighting",
  un = TRUE,
  binary = "std",
  continuous = "std",
  s.d.denom = "pooled"
)

print(bal_after)

love_after <- cobalt::love.plot(
  bal_after,
  stats = "mean.diffs",
  abs = TRUE,
  thresholds = c(m = 0.1),
  var.order = "unadjusted"
)

print(love_after)

# Weight summary ------------------------------------------------------------
ess <- function(w) {
  sum(w, na.rm = TRUE)^2 / sum(w^2, na.rm = TRUE)
}

weight_summary <- dplyr::bind_rows(lapply(seq_along(imp_list_w), function(i) {
  dat <- imp_list_w[[i]]
  
  data.frame(
    Imputation = i,
    Min = min(dat$sw, na.rm = TRUE),
    Mean = mean(dat$sw, na.rm = TRUE),
    Median = median(dat$sw, na.rm = TRUE),
    P95 = as.numeric(quantile(dat$sw, 0.95, na.rm = TRUE)),
    P99 = as.numeric(quantile(dat$sw, 0.99, na.rm = TRUE)),
    Max = max(dat$sw, na.rm = TRUE),
    ESS = ess(dat$sw)
  )
}))

print(weight_summary)

# Entropy-weighted unadjusted model ------------------------------------------------------------
fit_weighted_unadjusted <- function(imp_list_w, outcome_var) {
  
  yvar <- paste0(outcome_var, "_num")
  
  fit_list <- lapply(imp_list_w, function(dat) {
    glm(
      as.formula(paste0(yvar, " ~ pax_num")),
      data = dat,
      weights = sw,
      family = quasibinomial()
    )
  })
  
  pooled <- mice::pool(mice::as.mira(fit_list))
  
  summary(
    pooled,
    conf.int = TRUE,
    exponentiate = TRUE
  )
}

weighted_unadjusted_results <- lapply(outcomes, function(outcome) {
  fit_weighted_unadjusted(imp_list_w, outcome)
})

names(weighted_unadjusted_results) <- outcome_labels

# Entropy-weighted adjusted model ------------------------------------------------------------
fit_weighted_adjusted <- function(imp_list_w, outcome_var) {
  
  yvar <- paste0(outcome_var, "_num")
  
  fit_list <- lapply(imp_list_w, function(dat) {
    glm(
      as.formula(paste0(
        yvar,
        " ~ pax_num + age + sex + Number_comorb + Vaccination_number"
      )),
      data = dat,
      weights = sw,
      family = quasibinomial()
    )
  })
  
  pooled <- mice::pool(mice::as.mira(fit_list))
  
  summary(
    pooled,
    conf.int = TRUE,
    exponentiate = TRUE
  )
}

weighted_adjusted_results <- lapply(outcomes, function(outcome) {
  fit_weighted_adjusted(imp_list_w, outcome)
})

names(weighted_adjusted_results) <- outcome_labels

# Combine weighted regression results ------------------------------------------------------------
combine_results <- function(results_list, labels, method_name) {
  dplyr::bind_rows(lapply(seq_along(results_list), function(i) {
    results_list[[i]] %>%
      dplyr::mutate(
        Outcome = labels[i],
        Method = method_name
      )
  }))
}

full_weighted_unadjusted <- combine_results(
  weighted_unadjusted_results,
  outcome_labels,
  "Entropy weighted unadjusted logistic regression"
)

full_weighted_adjusted <- combine_results(
  weighted_adjusted_results,
  outcome_labels,
  "Entropy weighted adjusted logistic regression"
)

full_weighted_models <- dplyr::bind_rows(
  full_weighted_unadjusted,
  full_weighted_adjusted
)

# Formal doubly robust estimator: AIPW ------------------------------------------------------------

pool_scalar <- function(q, u) {
  m <- length(q)
  qbar <- mean(q, na.rm = TRUE)
  ubar <- mean(u, na.rm = TRUE)
  
  if (m == 1) {
    total_var <- ubar
    df <- Inf
  } else {
    b <- stats::var(q, na.rm = TRUE)
    total_var <- ubar + (1 + 1 / m) * b
    
    if (is.na(b) || b == 0) {
      df <- Inf
    } else {
      df <- (m - 1) * (1 + ubar / ((1 + 1 / m) * b))^2
    }
  }
  
  se <- sqrt(total_var)
  crit <- ifelse(is.infinite(df), qnorm(0.975), qt(0.975, df))
  p_value <- ifelse(
    is.infinite(df),
    2 * pnorm(abs(qbar / se), lower.tail = FALSE),
    2 * pt(abs(qbar / se), df = df, lower.tail = FALSE)
  )
  
  list(
    estimate = qbar,
    se = se,
    ci_low = qbar - crit * se,
    ci_high = qbar + crit * se,
    p_value = p_value,
    df = df
  )
}

aipw_one_dataset <- function(dat, outcome_var) {
  
  yvar <- paste0(outcome_var, "_num")
  
  Y <- dat[[yvar]]
  A <- dat$pax_num
  n <- nrow(dat)
  
  ps_model <- glm(
    pax_num ~ age + sex + Number_comorb + Vaccination_number,
    data = dat,
    family = binomial()
  )
  
  e <- predict(ps_model, type = "response")
  e <- pmin(pmax(e, 0.01), 0.99)
  
  outcome_model <- glm(
    as.formula(paste0(
      yvar,
      " ~ pax_num + age + sex + Number_comorb + Vaccination_number"
    )),
    data = dat,
    family = binomial()
  )
  
  dat1 <- dat
  dat1$pax_num <- 1
  
  dat0 <- dat
  dat0$pax_num <- 0
  
  m1 <- predict(outcome_model, newdata = dat1, type = "response")
  m0 <- predict(outcome_model, newdata = dat0, type = "response")
  
  phi1 <- m1 + A / e * (Y - m1)
  phi0 <- m0 + (1 - A) / (1 - e) * (Y - m0)
  
  mu1 <- mean(phi1, na.rm = TRUE)
  mu0 <- mean(phi0, na.rm = TRUE)
  
  inf1 <- phi1 - mu1
  inf0 <- phi0 - mu0
  
  eps <- 1e-6
  mu1b <- pmin(pmax(mu1, eps), 1 - eps)
  mu0b <- pmin(pmax(mu0, eps), 1 - eps)
  
  rd <- mu1 - mu0
  inf_rd <- inf1 - inf0
  se_rd <- sqrt(stats::var(inf_rd, na.rm = TRUE) / n)
  
  log_rr <- log(mu1b / mu0b)
  inf_log_rr <- inf1 / mu1b - inf0 / mu0b
  se_log_rr <- sqrt(stats::var(inf_log_rr, na.rm = TRUE) / n)
  
  log_or <- qlogis(mu1b) - qlogis(mu0b)
  inf_log_or <- inf1 / (mu1b * (1 - mu1b)) -
    inf0 / (mu0b * (1 - mu0b))
  se_log_or <- sqrt(stats::var(inf_log_or, na.rm = TRUE) / n)
  
  data.frame(
    mu1 = mu1,
    mu0 = mu0,
    rd = rd,
    se_rd = se_rd,
    log_rr = log_rr,
    se_log_rr = se_log_rr,
    log_or = log_or,
    se_log_or = se_log_or
  )
}

aipw_for_outcome <- function(imp_list_model, outcome_var, outcome_label) {
  
  per_imp <- dplyr::bind_rows(lapply(imp_list_model, function(dat) {
    aipw_one_dataset(dat, outcome_var)
  }))
  
  rd_pool <- pool_scalar(
    q = per_imp$rd,
    u = per_imp$se_rd^2
  )
  
  rr_pool <- pool_scalar(
    q = per_imp$log_rr,
    u = per_imp$se_log_rr^2
  )
  
  or_pool <- pool_scalar(
    q = per_imp$log_or,
    u = per_imp$se_log_or^2
  )
  
  risk1_pool <- pool_scalar(
    q = per_imp$mu1,
    u = rep(0, nrow(per_imp))
  )
  
  risk0_pool <- pool_scalar(
    q = per_imp$mu0,
    u = rep(0, nrow(per_imp))
  )
  
  data.frame(
    Outcome = outcome_label,
    Outcome_variable = outcome_var,
    Method = "AIPW doubly robust estimator",
    Estimand = c(
      "Risk if treated",
      "Risk if untreated",
      "Risk difference",
      "Risk ratio",
      "Marginal odds ratio"
    ),
    Estimate = c(
      risk1_pool$estimate,
      risk0_pool$estimate,
      rd_pool$estimate,
      exp(rr_pool$estimate),
      exp(or_pool$estimate)
    ),
    CI_low = c(
      risk1_pool$ci_low,
      risk0_pool$ci_low,
      rd_pool$ci_low,
      exp(rr_pool$ci_low),
      exp(or_pool$ci_low)
    ),
    CI_high = c(
      risk1_pool$ci_high,
      risk0_pool$ci_high,
      rd_pool$ci_high,
      exp(rr_pool$ci_high),
      exp(or_pool$ci_high)
    ),
    p_value = c(
      NA,
      NA,
      rd_pool$p_value,
      rr_pool$p_value,
      or_pool$p_value
    )
  )
}

aipw_results <- dplyr::bind_rows(lapply(seq_along(outcomes), function(i) {
  aipw_for_outcome(
    imp_list_model = imp_list_model,
    outcome_var = outcomes[i],
    outcome_label = outcome_labels[i]
  )
}))

cat("\nFormal AIPW doubly robust results:\n")
print(aipw_results)

# Paxlovid effect from weighted regressions ------------------------------------------------------------
pax_effect_weighted <- full_weighted_models %>%
  dplyr::filter(term == "pax_num") %>%
  dplyr::transmute(
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

print(pax_effect_weighted)

# Save results ------------------------------------------------------------
write.csv(
  missing_table,
  "missingness_summary.csv",
  row.names = FALSE
)

write.csv(
  weight_summary,
  "entropy_weight_summary.csv",
  row.names = FALSE
)

write.csv(
  full_weighted_models,
  "entropy_weighted_regression_results.csv",
  row.names = FALSE
)

write.csv(
  pax_effect_weighted,
  "entropy_weighted_paxlovid_effect.csv",
  row.names = FALSE
)

write.csv(
  aipw_results,
  "aipw_doubly_robust_results.csv",
  row.names = FALSE
)

ggsave(
  "balance_after_entropy_weighting.png",
  love_after,
  width = 8,
  height = 5,
  dpi = 300
)

wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb, "Missingness")
openxlsx::writeData(wb, "Missingness", missing_table)

openxlsx::addWorksheet(wb, "Weight Summary")
openxlsx::writeData(wb, "Weight Summary", weight_summary)

openxlsx::addWorksheet(wb, "Weighted Regressions")
openxlsx::writeData(wb, "Weighted Regressions", full_weighted_models)

openxlsx::addWorksheet(wb, "Weighted Pax Effect")
openxlsx::writeData(wb, "Weighted Pax Effect", pax_effect_weighted)

openxlsx::addWorksheet(wb, "AIPW Doubly Robust")
openxlsx::writeData(wb, "AIPW Doubly Robust", aipw_results)

openxlsx::saveWorkbook(
  wb,
  "entropy_weighted_and_aipw_results.xlsx",
  overwrite = TRUE
)

sessionInfo()

# Missingness indicators
df$miss_patient <- ifelse(is.na(df$improve_patient), 1, 0)
df$miss_clinician <- ifelse(is.na(df$improve_clinician), 1, 0)
df$miss_grorcoop <- ifelse(is.na(df$new_improve_grorcoop), 1, 0)
df$miss_grandcoop <- ifelse(is.na(df$new_improved_grandcoop), 1, 0)

# Check missingness by treatment and sex
table(df$miss_patient, df$paxlovid_treatment, useNA = "ifany")
table(df$miss_patient, df$sex, useNA = "ifany")

# Compare age by missingness
t.test(age ~ miss_patient, data = df)

# Compare comorbidities by missingness
t.test(Number_comorb ~ miss_patient, data = df)

# Logistic model for missingness in patient-reported outcome
miss_model_patient <- glm(
  miss_patient ~ paxlovid_treatment + age + sex +
    Number_comorb + Vaccination_number,
  data = df,
  family = binomial()
)

summary(miss_model_patient)

# Optional: missingness in any outcome
df$miss_any_outcome <- ifelse(
  is.na(df$improve_patient) |
    is.na(df$improve_clinician) |
    is.na(df$new_improve_grorcoop) |
    is.na(df$new_improved_grandcoop),
  1,
  0
)

miss_model_any <- glm(
  miss_any_outcome ~ paxlovid_treatment + age + sex +
    Number_comorb + Vaccination_number,
  data = df,
  family = binomial()
)

summary(miss_model_any)

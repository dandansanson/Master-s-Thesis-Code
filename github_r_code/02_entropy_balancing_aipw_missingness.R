
# Load data ----
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

lc_pat <- data_raw[, vars_needed, drop = FALSE]

# Recode variables ----
lc_pat$sex <- as.factor(lc_pat$sex)
lc_pat$paxlovid_treatment <- as.factor(lc_pat$paxlovid_treatment)

for (v in outcomes) {
  lc_pat[[v]] <- as.factor(lc_pat[[v]])
}

lc_pat$age <- as.numeric(lc_pat$age)
lc_pat$Number_comorb <- as.numeric(lc_pat$Number_comorb)
lc_pat$Vaccination_number <- as.numeric(lc_pat$Vaccination_number)

# Missingness check ----
missing_table <- data.frame(
  Variable = names(lc_pat),
  Missing_n = colSums(is.na(lc_pat)),
  Missing_pct = round(colMeans(is.na(lc_pat)) * 100, 2)
)

missing_table <- missing_table[order(-missing_table$Missing_pct), ]

print(missing_table)
print(mice::md.pattern(lc_pat, plot = FALSE))

# Multiple imputation ----
choose_factor_method <- function(x) {
  if (!any(is.na(x))) return("")
  nlev <- nlevels(droplevels(x))
  if (nlev == 2) return("logreg")
  if (nlev > 2) return("polyreg")
  return("")
}

factor_vars <- c(outcomes, "paxlovid_treatment", "sex")
numeric_vars <- c("age", "Number_comorb", "Vaccination_number")

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

  imp <- mice::mice(
    lc_pat,
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
  imp_list <- list(lc_pat)
}

# Convert treatment and outcomes to 0/1 numeric ----
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

add_numeric_vars <- function(lc_patient) {
  lc_patient$pax_num <- make_binary_numeric(lc_patient[[treatment_var]])
  lc_patient$sex <- as.factor(lc_patient$sex)

  for (outcome in outcomes) {
    lc_patient[[paste0(outcome, "_num")]] <- make_binary_numeric(lc_patient[[outcome]])
  }

  lc_patient
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

# Entropy balancing ----
balance_formula <- pax_num ~ age + sex + Number_comorb + Vaccination_number

add_entropy_weights <- function(lc_patient) {
  lc_patient$sex <- as.factor(lc_patient$sex)

  w_obj <- WeightIt::weightit(
    balance_formula,
    data = lc_patient,
    method = "ebal",
    estimand = "ATE"
  )

  lc_patient$sw <- w_obj$weights

  list(
    data = lc_patient,
    weightit_object = w_obj
  )
}

weight_list <- lapply(imp_list_model, add_entropy_weights)

imp_list_w <- lapply(weight_list, `[[`, "data")
weightit_objects <- lapply(weight_list, `[[`, "weightit_object")

# Covariate balance ----
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

# Weight summary ----
ess <- function(w) {
  sum(w, na.rm = TRUE)^2 / sum(w^2, na.rm = TRUE)
}

weight_summary <- dplyr::bind_rows(lapply(seq_along(imp_list_w), function(i) {
  lc_patient <- imp_list_w[[i]]

  data.frame(
    Imputation = i,
    Min = min(lc_patient$sw, na.rm = TRUE),
    Mean = mean(lc_patient$sw, na.rm = TRUE),
    Median = median(lc_patient$sw, na.rm = TRUE),
    P95 = as.numeric(quantile(lc_patient$sw, 0.95, na.rm = TRUE)),
    P99 = as.numeric(quantile(lc_patient$sw, 0.99, na.rm = TRUE)),
    Max = max(lc_patient$sw, na.rm = TRUE),
    ESS = ess(lc_patient$sw)
  )
}))

print(weight_summary)

# Entropy-weighted unadjusted model ----
fit_weighted_unadjusted <- function(imp_list_w, outcome_var) {
  yvar <- paste0(outcome_var, "_num")

  fit_list <- lapply(imp_list_w, function(lc_patient) {
    glm(
      as.formula(paste0(yvar, " ~ pax_num")),
      data = lc_patient,
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

# Entropy-weighted adjusted model ----
fit_weighted_adjusted <- function(imp_list_w, outcome_var) {
  yvar <- paste0(outcome_var, "_num")

  fit_list <- lapply(imp_list_w, function(lc_patient) {
    glm(
      as.formula(paste0(
        yvar,
        " ~ pax_num + age + sex + Number_comorb + Vaccination_number"
      )),
      data = lc_patient,
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

# Combine weighted regression results ----
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

# Formal doubly robust estimator: AIPW ----
pool_scalar <- function(q, u) {
  m <- length(q)
  qbar <- mean(q, na.rm = TRUE)
  ubar <- mean(u, na.rm = TRUE)

  if (m == 1) {
    total_var <- ubar
    lc_pat <- Inf
  } else {
    b <- stats::var(q, na.rm = TRUE)
    total_var <- ubar + (1 + 1 / m) * b

    if (is.na(b) || b == 0) {
      lc_pat <- Inf
    } else {
      lc_pat <- (m - 1) * (1 + ubar / ((1 + 1 / m) * b))^2
    }
  }

  se <- sqrt(total_var)
  crit <- ifelse(is.infinite(lc_pat), qnorm(0.975), qt(0.975, lc_pat))
  p_value <- ifelse(
    is.infinite(lc_pat),
    2 * pnorm(abs(qbar / se), lower.tail = FALSE),
    2 * pt(abs(qbar / se), lc_pat = lc_pat, lower.tail = FALSE)
  )

  list(
    estimate = qbar,
    se = se,
    ci_low = qbar - crit * se,
    ci_high = qbar + crit * se,
    p_value = p_value,
    lc_pat = lc_pat
  )
}

aipw_one_dataset <- function(lc_patient, outcome_var) {
  yvar <- paste0(outcome_var, "_num")

  Y <- lc_patient[[yvar]]
  A <- lc_patient$pax_num
  n <- nrow(lc_patient)

  ps_model <- glm(
    pax_num ~ age + sex + Number_comorb + Vaccination_number,
    data = lc_patient,
    family = binomial()
  )

  e <- predict(ps_model, type = "response")
  e <- pmin(pmax(e, 0.01), 0.99)

  outcome_model <- glm(
    as.formula(paste0(
      yvar,
      " ~ pax_num + age + sex + Number_comorb + Vaccination_number"
    )),
    data = lc_patient,
    family = binomial()
  )

  lc_patient1 <- lc_patient
  lc_patient1$pax_num <- 1

  lc_patient0 <- lc_patient
  lc_patient0$pax_num <- 0

  m1 <- predict(outcome_model, newdata = lc_patient1, type = "response")
  m0 <- predict(outcome_model, newdata = lc_patient0, type = "response")

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
  per_imp <- dplyr::bind_rows(lapply(imp_list_model, function(lc_patient) {
    aipw_one_dataset(lc_patient, outcome_var)
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

# Paxlovid effect from weighted regressions ----
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


# Missingness indicators
lc_pat$miss_patient <- ifelse(is.na(lc_pat$improve_patient), 1, 0)
lc_pat$miss_clinician <- ifelse(is.na(lc_pat$improve_clinician), 1, 0)
lc_pat$miss_grorcoop <- ifelse(is.na(lc_pat$new_improve_grorcoop), 1, 0)
lc_pat$miss_grandcoop <- ifelse(is.na(lc_pat$new_improved_grandcoop), 1, 0)

# Check missingness by treatment and sex
table(lc_pat$miss_patient, lc_pat$paxlovid_treatment, useNA = "ifany")
table(lc_pat$miss_patient, lc_pat$sex, useNA = "ifany")

# Compare age by missingness
t.test(age ~ miss_patient, data = lc_pat)

# Compare comorbidities by missingness
t.test(Number_comorb ~ miss_patient, data = lc_pat)

# Logistic model for missingness in patient-reported outcome
miss_model_patient <- glm(
  miss_patient ~ paxlovid_treatment + age + sex +
    Number_comorb + Vaccination_number,
  data = lc_pat,
  family = binomial()
)

summary(miss_model_patient)

# Optional: missingness in any outcome
lc_pat$miss_any_outcome <- ifelse(
  is.na(lc_pat$improve_patient) |
    is.na(lc_pat$improve_clinician) |
    is.na(lc_pat$new_improve_grorcoop) |
    is.na(lc_pat$new_improved_grandcoop),
  1,
  0
)

miss_model_any <- glm(
  miss_any_outcome ~ paxlovid_treatment + age + sex +
    Number_comorb + Vaccination_number,
  data = lc_pat,
  family = binomial()
)

summary(miss_model_any)

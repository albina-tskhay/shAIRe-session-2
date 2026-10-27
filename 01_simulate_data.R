
# =============================================================================
# 01_simulate_data.R
# Simulate a realistic observational dataset for the causal inference workshop
#
# Scenario: comparing two biologics for severe asthma in routine care
#   A = 1: dupilumab
#   A = 0: mepolizumab
#   Y:     asthma exacerbation within 12 months (1 = yes, 0 = no)
#
# Confounding by indication: sicker patients are more likely to receive
# dupilumab AND more likely to exacerbate, creating a spurious association.
#
# # =============================================================================

library(tidyverse)

simulate_asthma <- function(n = 10000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # --- Step 1: baseline confounders ------------------------------------------
  age           <- rnorm(n, mean = 48, sd = 12)              # years
  sex           <- rbinom(n, 1, 0.55)                        # 1 = female
  bmi           <- rnorm(n, mean = 29, sd = 5)               # kg/m^2
  smoking       <- rbinom(n, 1, 0.25)                        # 1 = current/former
  prev_exac     <- pmin(rpois(n, lambda = 1.2), 5)           # exacerbations, prior year
  fev1_pct      <- pmax(pmin(rnorm(n, 68, 15), 120), 20)     # FEV1 % predicted
  eosinophils   <- round(pmin(rlnorm(n, 5.5, 0.6), 1500))    # cells/uL
  allergic_rhin <- rbinom(n, 1, 0.40)                        # 1 = yes
  
  # --- Step 2: treatment assignment (TRUE propensity score model) ------------
  # Confounding by indication: patients with more exacerbations, lower FEV1,
  # and allergic comorbidity are more often given dupilumab; patients with very
  # high eosinophils are more often given mepolizumab.
  # Two HIDDEN terms (age x smoking, quadratic BMI) also drive prescribing.
  lp_treatment <- -1.4 +
    0.008 * age +
    0.10  * sex +
    0.015 * bmi +
    0.25  * smoking +
    0.50  * prev_exac -
    0.025 * fev1_pct -
    0.002 * eosinophils +
    0.30  * allergic_rhin +
    0.008 * age * smoking +          # hidden interaction
    0.020 * (bmi - 29)^2             # hidden quadratic term
  treatment <- rbinom(n, 1, plogis(lp_treatment))
  
  # --- Step 3: outcome (TRUE outcome model) -----------------------------------
  # Dupilumab has a protective causal effect (log-OR = -0.40).
  # The same two hidden terms also affect exacerbation risk.
  outcome_lp <- function(a) {
    -1.8 -
      0.40  * a +                    # true causal effect (conditional log-OR)
      0.020 * age -
      0.15  * sex +
      0.025 * bmi +
      0.50  * smoking +
      0.45  * prev_exac -
      0.022 * fev1_pct +
      0.002 * eosinophils +
      0.20  * allergic_rhin +
      0.012 * age * smoking +        # hidden interaction
      0.004 * (bmi - 29)^2           # hidden quadratic term
  }
  exacerbation <- rbinom(n, 1, plogis(outcome_lp(treatment)))
  
  # --- Step 4: true causal effects (both potential outcomes for everyone) ----
  p_Y1 <- plogis(outcome_lp(1))      # risk if everyone received dupilumab
  p_Y0 <- plogis(outcome_lp(0))      # risk if everyone received mepolizumab
  truth <- c(
    rd             = mean(p_Y1 - p_Y0),
    rr             = mean(p_Y1) / mean(p_Y0),
    marginal_or    = (mean(p_Y1) / (1 - mean(p_Y1))) / (mean(p_Y0) / (1 - mean(p_Y0))),
    conditional_or = exp(-0.40)
  )
  
  data <- data.frame(
    id = seq_len(n), treatment, exacerbation,
    age = round(age, 1), sex, bmi = round(bmi, 1), smoking, prev_exac,
    fev1_pct = round(fev1_pct, 1), eosinophils, allergic_rhin
  )
  list(data = data, truth = truth)
}

sim <- simulate_asthma(n = 10000, seed = 42)

cat("\n========== TRUE CAUSAL EFFECTS ==========\n")
cat("True ATE (risk difference):", round(sim$truth["rd"], 4), "\n")
cat("True ATE (risk ratio):     ", round(sim$truth["rr"], 4), "\n")
cat("True marginal OR:          ", round(sim$truth["marginal_or"], 4), "\n")
cat("True conditional OR:       ", round(sim$truth["conditional_or"], 4), "\n")
cat("Exacerbation rate (overall):", round(mean(sim$data$exacerbation), 3), "\n")
cat("Proportion on dupilumab:   ", round(mean(sim$data$treatment), 3), "\n")
cat("==========================================\n\n")

# Add labeled factors for display
sim_data <- sim$data %>%
  mutate(
    treatment_lab     = factor(treatment, levels = c(0, 1),
                               labels = c("Mepolizumab", "Dupilumab")),
    exacerbation_lab  = factor(exacerbation, levels = c(0, 1), labels = c("No", "Yes")),
    sex_lab           = factor(sex, levels = c(0, 1), labels = c("Male", "Female")),
    smoking_lab       = factor(smoking, levels = c(0, 1), labels = c("Never", "Current/Former")),
    allergic_rhin_lab = factor(allergic_rhin, levels = c(0, 1), labels = c("No", "Yes"))
  )

write_csv(sim_data, "sim_asthma_data.csv")
cat("Dataset saved to sim_asthma_data.csv (n =", nrow(sim_data), ")\n")

true_effects <- tibble(
  estimand   = c("ATE (Risk Difference)", "ATE (Risk Ratio)",
                 "Marginal OR", "Conditional OR"),
  true_value = unname(sim$truth)
)
write_csv(true_effects, "true_effects.csv")
cat("True effects saved to true_effects.csv\n")
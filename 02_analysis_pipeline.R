# =============================================================================
# 02_analysis_pipeline.R
# Causal Inference Workshop: IPTW and G-Methods
#
# Question: what is the effect of dupilumab (A = 1) vs mepolizumab (A = 0)
#           on the 12-month risk of asthma exacerbation?
#
#   Part 1: Descriptive statistics and positivity
#   Part 2: Unadjusted (crude) estimation
#   Part 3: Conventional logistic regression
#   Part 4: Inverse probability of treatment weighting (IPTW)
#   Part 5: G-formula (parametric standardization)
#   Part 6: Doubly robust estimation (AIPW) [bonus]
#   Part 7: Compare all estimates to the truth
#  #
# Every causal method is run twice:
#   SIMPLE  = main effects only (what most analysts fit). Misspecified here,
#             because the true models contain an age x smoking interaction and
#             a quadratic BMI term.
#   CORRECT = includes those two terms (matches the true data-generating
#             models; only knowable because we simulated the data).
# =============================================================================

library(tidyverse)
library(tableone)
library(cobalt)
library(sandwich)
library(lmtest)
library(boot)
library(survey)

sim_data     <- read_csv("sim_asthma_data.csv", show_col_types = FALSE)
true_effects <- read_csv("true_effects.csv", show_col_types = FALSE)
true_rd  <- true_effects$true_value[true_effects$estimand == "ATE (Risk Difference)"]
true_mor <- true_effects$true_value[true_effects$estimand == "Marginal OR"]
true_cor <- true_effects$true_value[true_effects$estimand == "Conditional OR"]

sim_data <- sim_data %>%
  mutate(treatment_lab = factor(treatment_lab, levels = c("Mepolizumab", "Dupilumab")))

cat("True ATE (risk difference, dupilumab - mepolizumab):", round(true_rd, 4), "\n\n")

# Colors used in all plots
col_trt <- c("Mepolizumab" = "#4393C3", "Dupilumab" = "#D6604D")

# --- Model formulas used throughout ------------------------------------------
ps_simple   <- treatment ~ age + sex + bmi + smoking + prev_exac +
  fev1_pct + eosinophils + allergic_rhin
ps_correct  <- update(ps_simple, . ~ . + age:smoking + I((bmi - 29)^2))

out_simple  <- exacerbation ~ treatment + age + sex + bmi + smoking + prev_exac +
  fev1_pct + eosinophils + allergic_rhin
out_correct <- update(out_simple, . ~ . + age:smoking + I((bmi - 29)^2))


# =============================================================================
# PART 1: DESCRIPTIVE STATISTICS AND POSITIVITY
# =============================================================================

cat("\n========== PART 1: DESCRIPTIVE STATISTICS ==========\n\n")

# --- 1a. Table 1: baseline characteristics by treatment ---
covariates <- c("age", "sex", "bmi", "smoking", "prev_exac",
                "fev1_pct", "eosinophils", "allergic_rhin")

table1 <- CreateTableOne(vars = covariates, strata = "treatment_lab",
                         data = sim_data, test = FALSE)
print(table1, smd = TRUE)
# SMD > 0.1 indicates meaningful imbalance (confounding by indication)

# --- 1b. Estimate the propensity score (probability of receiving dupilumab) ---
ps_model_simple  <- glm(ps_simple,  data = sim_data, family = binomial)
ps_model_correct <- glm(ps_correct, data = sim_data, family = binomial)

sim_data$ps_simple  <- predict(ps_model_simple,  type = "response")
sim_data$ps_correct <- predict(ps_model_correct, type = "response")

# --- 1c. Positivity: do both groups cover the same PS range? ---
cat("\nPropensity score (correct model) by treatment group:\n")
sim_data %>%
  group_by(treatment_lab) %>%
  summarise(min = min(ps_correct), q25 = quantile(ps_correct, 0.25),
            median = median(ps_correct), q75 = quantile(ps_correct, 0.75),
            max = max(ps_correct), .groups = "drop") %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  print()

ps_plot <- ggplot(sim_data, aes(x = ps_correct, fill = treatment_lab)) +
  geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
  scale_fill_manual(values = col_trt) +
  labs(title = "Propensity Score Distribution by Treatment",
       x = "Estimated probability of receiving dupilumab", y = "Count", fill = NULL) +
  theme_minimal(base_size = 16) +
  theme(legend.position = "top")
print(ps_plot)

cat("\nObservations with extreme PS (< 0.02 or > 0.98):",
    sum(sim_data$ps_correct < 0.02 | sim_data$ps_correct > 0.98), "\n")


# =============================================================================
# PART 2: UNADJUSTED (CRUDE) ESTIMATION
# =============================================================================

cat("\n========== PART 2: UNADJUSTED ESTIMATION ==========\n\n")

risk_by_group <- sim_data %>%
  group_by(treatment_lab) %>%
  summarise(n = n(), events = sum(exacerbation), risk = mean(exacerbation), .groups = "drop")
print(risk_by_group)

risk_dupi <- risk_by_group$risk[risk_by_group$treatment_lab == "Dupilumab"]
risk_mepo <- risk_by_group$risk[risk_by_group$treatment_lab == "Mepolizumab"]
crude_rd  <- risk_dupi - risk_mepo
crude_or  <- (risk_dupi / (1 - risk_dupi)) / (risk_mepo / (1 - risk_mepo))

cat("\nCrude risk difference (dupilumab - mepolizumab):", round(crude_rd, 4), "\n")
cat("True ATE:", round(true_rd, 4), "   Bias:", round(crude_rd - true_rd, 4), "\n")
cat("--> Dupilumab looks worse than it is: it is given to sicker patients.\n")


# =============================================================================
# PART 3: CONVENTIONAL LOGISTIC REGRESSION
# =============================================================================

cat("\n========== PART 3: LOGISTIC REGRESSION ==========\n\n")

# The usual approach: outcome ~ treatment + confounders, report the adjusted OR
lr_simple  <- glm(out_simple,  data = sim_data, family = binomial)
lr_correct <- glm(out_correct, data = sim_data, family = binomial)

or_lr_simple  <- exp(coef(lr_simple)["treatment"])
or_lr_correct <- exp(coef(lr_correct)["treatment"])

cat("Adjusted OR, simple model (main effects):", round(or_lr_simple, 3), "\n")
cat("Adjusted OR, correct model:              ", round(or_lr_correct, 3), "\n")
cat("True CONDITIONAL OR:", round(true_cor, 3), "   True MARGINAL OR:", round(true_mor, 3), "\n")
# Logistic regression reports a CONDITIONAL odds ratio. Even when the
# model is correct, it estimates a different quantity than the marginal
#(population-average) effect, because the OR is non-collapsible.
#  It also gives no risk difference, which is what clinicians need.


# =============================================================================
# PART 4: INVERSE PROBABILITY OF TREATMENT WEIGHTING (IPTW)
# =============================================================================

cat("\n========== PART 4: IPTW ==========\n\n")

# --- 4a. Stabilized weights from each PS model ---
# Dupilumab:   P(A = 1) / PS        Mepolizumab: P(A = 0) / (1 - PS)
p_dupi <- mean(sim_data$treatment)
sim_data <- sim_data %>%
  mutate(
    sw_simple  = ifelse(treatment == 1, p_dupi / ps_simple,  (1 - p_dupi) / (1 - ps_simple)),
    sw_correct = ifelse(treatment == 1, p_dupi / ps_correct, (1 - p_dupi) / (1 - ps_correct))
  )

cat("Stabilized weights, simple PS:\n");  print(summary(sim_data$sw_simple))
cat("Stabilized weights, correct PS:\n"); print(summary(sim_data$sw_correct))

# --- 4c. Estimate the ATE ---
# Weighted risks in each arm (identical to a weighted GLM of Y on A)
iptw_risks <- function(y, a, w) c(p1 = weighted.mean(y[a == 1], w[a == 1]),
                                  p0 = weighted.mean(y[a == 0], w[a == 0]))

r_iptw_simple  <- iptw_risks(sim_data$exacerbation, sim_data$treatment, sim_data$sw_simple)
r_iptw_correct <- iptw_risks(sim_data$exacerbation, sim_data$treatment, sim_data$sw_correct)
iptw_rd_simple  <- unname(r_iptw_simple["p1"]  - r_iptw_simple["p0"])
iptw_rd_correct <- unname(r_iptw_correct["p1"] - r_iptw_correct["p0"])

# Same estimate via a weighted regression with robust standard errors
design_iptw <- svydesign(ids = ~1, weights = ~sw_correct, data = sim_data)
iptw_glm    <- svyglm(exacerbation ~ treatment, design = design_iptw, family = quasibinomial)
cat("\nWeighted GLM (correct PS), log-OR for dupilumab:\n")
print(summary(iptw_glm)$coefficients)

cat("\nIPTW risk difference, simple PS: ", round(iptw_rd_simple, 4),
    "  Bias:", round(iptw_rd_simple - true_rd, 4), "\n")
cat("IPTW risk difference, correct PS:", round(iptw_rd_correct, 4),
    "  Bias:", round(iptw_rd_correct - true_rd, 4), "\n")


# =============================================================================
# PART 5: G-FORMULA (PARAMETRIC STANDARDIZATION)
# =============================================================================

cat("\n========== PART 5: G-FORMULA ==========\n\n")

# Step 1: fit an outcome model. Step 2: predict everyone's risk under
# dupilumab and under mepolizumab. Step 3: average and contrast.
gformula_risks <- function(formula, data) {
  fit <- glm(formula, data = data, family = binomial)
  c(p1 = mean(predict(fit, newdata = mutate(data, treatment = 1), type = "response")),
    p0 = mean(predict(fit, newdata = mutate(data, treatment = 0), type = "response")))
}

r_gf_simple  <- gformula_risks(out_simple,  sim_data)
r_gf_correct <- gformula_risks(out_correct, sim_data)
gf_rd_simple  <- unname(r_gf_simple["p1"]  - r_gf_simple["p0"])
gf_rd_correct <- unname(r_gf_correct["p1"] - r_gf_correct["p0"])

cat("Risk if everyone got dupilumab (correct model): ", round(r_gf_correct["p1"], 4), "\n")
cat("Risk if everyone got mepolizumab (correct model):", round(r_gf_correct["p0"], 4), "\n\n")
cat("G-formula risk difference, simple model: ", round(gf_rd_simple, 4),
    "  Bias:", round(gf_rd_simple - true_rd, 4), "\n")
cat("G-formula risk difference, correct model:", round(gf_rd_correct, 4),
    "  Bias:", round(gf_rd_correct - true_rd, 4), "\n")


# =============================================================================
# PART 6: DOUBLY ROBUST ESTIMATION (AIPW) [BONUS]
# =============================================================================

cat("\n========== PART 6: DOUBLY ROBUST (AIPW) ==========\n\n")

# AIPW = g-formula prediction + IPTW-weighted average of the residuals.
# Consistent if EITHER the PS model OR the outcome model is correct.
aipw_risks <- function(ps_formula, out_formula, data) {
  ps  <- predict(glm(ps_formula, data = data, family = binomial), type = "response")
  fit <- glm(out_formula, data = data, family = binomial)
  mu1 <- predict(fit, newdata = mutate(data, treatment = 1), type = "response")
  mu0 <- predict(fit, newdata = mutate(data, treatment = 0), type = "response")
  a <- data$treatment; y <- data$exacerbation
  c(p1 = mean(mu1 + a / ps * (y - mu1)),
    p0 = mean(mu0 + (1 - a) / (1 - ps) * (y - mu0)))
}

r_aipw_simple  <- aipw_risks(ps_simple,  out_simple,  sim_data)   # both wrong
r_aipw_one     <- aipw_risks(ps_correct, out_simple,  sim_data)   # only PS correct
r_aipw_correct <- aipw_risks(ps_correct, out_correct, sim_data)   # both correct
aipw_rd_simple  <- unname(r_aipw_simple["p1"]  - r_aipw_simple["p0"])
aipw_rd_one     <- unname(r_aipw_one["p1"]     - r_aipw_one["p0"])
aipw_rd_correct <- unname(r_aipw_correct["p1"] - r_aipw_correct["p0"])

cat("AIPW, both models simple:        ", round(aipw_rd_simple, 4),  "  Bias:", round(aipw_rd_simple - true_rd, 4), "\n")
cat("AIPW, only the PS model correct: ", round(aipw_rd_one, 4),     "  Bias:", round(aipw_rd_one - true_rd, 4), "\n")
cat("AIPW, both models correct:       ", round(aipw_rd_correct, 4), "  Bias:", round(aipw_rd_correct - true_rd, 4), "\n")
cat("--> One correct model is enough: that is double robustness.\n")


# =============================================================================
# PART 7: COMPARE ALL ESTIMATES
# =============================================================================

cat("\n========== PART 7: SUMMARY ==========\n\n")

# --- 7a. Bootstrap 95% CIs (resample patients, refit every model) ---
# One bootstrap run computes all estimators at once. R = 200 takes about 1
# to 2 minutes; set run_bootstrap <- FALSE to skip it during the live session.
run_bootstrap <- TRUE

all_rd <- function(data, idx = seq_len(nrow(data))) {
  d  <- data[idx, ]
  pd <- mean(d$treatment)
  sw <- function(f) {
    ps <- predict(glm(f, data = d, family = binomial), type = "response")
    ifelse(d$treatment == 1, pd / ps, (1 - pd) / (1 - ps))
  }
  rd <- function(r) unname(r["p1"] - r["p0"])
  c(crude         = mean(d$exacerbation[d$treatment == 1]) - mean(d$exacerbation[d$treatment == 0]),
    iptw_simple   = rd(iptw_risks(d$exacerbation, d$treatment, sw(ps_simple))),
    iptw_correct  = rd(iptw_risks(d$exacerbation, d$treatment, sw(ps_correct))),
    gf_simple     = rd(gformula_risks(out_simple, d)),
    gf_correct    = rd(gformula_risks(out_correct, d)),
    aipw_simple   = rd(aipw_risks(ps_simple, out_simple, d)),
    aipw_one      = rd(aipw_risks(ps_correct, out_simple, d)),
    aipw_correct  = rd(aipw_risks(ps_correct, out_correct, d)))
}

method_labels <- c(
  crude        = "Unadjusted (crude)",
  iptw_simple  = "IPTW, simple PS",
  iptw_correct = "IPTW, correct PS",
  gf_simple    = "G-formula, simple model",
  gf_correct   = "G-formula, correct model",
  aipw_simple  = "AIPW, both simple",
  aipw_one     = "AIPW, only PS correct",
  aipw_correct = "AIPW, both correct")

estimates <- all_rd(sim_data)

if (run_bootstrap) {
  set.seed(123)
  boot_all <- boot(sim_data, all_rd, R = 200)
  ci <- t(sapply(seq_along(estimates), function(i)
    boot.ci(boot_all, type = "perc", index = i)$percent[4:5]))
} else {
  ci <- matrix(NA, nrow = length(estimates), ncol = 2)
}

results <- tibble(
  key    = names(estimates),
  Method = method_labels[names(estimates)],
  Model  = case_when(key == "crude" ~ "No adjustment",
                     grepl("correct", key) ~ "Correctly specified",
                     key == "aipw_one" ~ "One model correct",
                     TRUE ~ "Misspecified (main effects only)"),
  RD     = estimates,
  lower  = ci[, 1],
  upper  = ci[, 2],
  Bias   = estimates - true_rd,
  `% of true effect` = 100 * estimates / true_rd
)

results %>%
  select(Method, RD, lower, upper, Bias, `% of true effect`) %>%
  mutate(across(c(RD, lower, upper, Bias), ~ round(.x, 4)),
         `% of true effect` = round(`% of true effect`, 0)) %>%
  print(n = Inf)

# --- 7b. Odds ratio scale: where logistic regression fits in ---
marg_or <- function(r) unname((r["p1"] / (1 - r["p1"])) / (r["p0"] / (1 - r["p0"])))
or_table <- tibble(
  Method = c("TRUE marginal OR (target)", "TRUE conditional OR",
             "Unadjusted", "Logistic regression, simple (conditional)",
             "Logistic regression, correct (conditional)",
             "IPTW, simple PS", "IPTW, correct PS",
             "G-formula, simple model", "G-formula, correct model",
             "AIPW, both correct"),
  OR = c(true_mor, true_cor, crude_or, or_lr_simple, or_lr_correct,
         marg_or(r_iptw_simple), marg_or(r_iptw_correct),
         marg_or(r_gf_simple), marg_or(r_gf_correct), marg_or(r_aipw_correct))
) %>% mutate(OR = round(OR, 3))
cat("\nOdds ratio scale:\n")
print(or_table, n = Inf)

# --- 7c. Summary plot ---
plot_df <- results %>%
  mutate(Method = fct_rev(factor(Method, levels = method_labels)),
         Model  = factor(Model, levels = c("No adjustment", "Misspecified (main effects only)",
                                           "One model correct", "Correctly specified")))

comparison_plot <- ggplot(plot_df, aes(x = RD, y = Method, color = Model)) +
  geom_vline(xintercept = true_rd, linetype = "dashed", color = "black", linewidth = 1) +
  geom_vline(xintercept = 0, color = "grey70") +
  { if (run_bootstrap) geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.25, linewidth = 0.9, orientation = "y") } +
  geom_point(size = 4.5) +
  annotate("text", x = true_rd, y = 8.45, label = paste("True ATE =", round(true_rd, 3)),
           hjust = 1.05, size = 5) +
  scale_color_manual(values = c("No adjustment" = "#B2182B",
                                "Misspecified (main effects only)" = "#F4A582",
                                "One model correct" = "#92C5DE",
                                "Correctly specified" = "#2166AC")) +
  labs(title = "Dupilumab vs Mepolizumab: Estimated vs True ATE",
       x = "Risk difference (dupilumab - mepolizumab)", y = NULL, color = NULL) +
  theme_minimal(base_size = 16) +
  theme(legend.position = "bottom", legend.direction = "vertical")
print(comparison_plot)



cat("\n========== WORKSHOP COMPLETE ==========\n")

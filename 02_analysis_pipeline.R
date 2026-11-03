# =============================================================================
# 02_analysis_pipeline.R
# Causal inference workshop: IPTW, g-formula, and doubly robust standardization
#
# Question: effect of mepolizumab (treatment = 1) vs dupilumab (treatment = 0)
#           on the 12-month risk of exacerbation (exacerbation = 1)
#
# Each method is run twice:
#   "simple"  = main effects only (what we usually fit; misspecified here)
#   "correct" = adds age*smoking and I((bmi - 29)^2), as in the true models

# =============================================================================

library(tableone)
library(survey)
library(sandwich)
library(lmtest)
library(ggplot2)

ObsData <- read.csv("sim_asthma_data.csv")
TRUE_ATE <- read.csv("true_effects.csv")$TRUE_ATE
TRUE_ATE   # the number we want to recover: about -0.083


# =============================================================================
# PART 1: DESCRIPTIVE STATISTICS AND POSITIVITY
# =============================================================================

# Table 1 by treatment group, with standardized mean differences (SMD)
tab <- CreateTableOne(vars = c("age", "sex", "bmi", "smoking", "prev_exac", "fev1_pct", "eosinophils", "allergic_rhin"),
                      strata = "treatment", data = ObsData, test = FALSE)
print(tab, smd = TRUE)

# Propensity score: probability of receiving mepolizumab given the confounders
e <- glm(treatment ~ age + sex + bmi + smoking + prev_exac + fev1_pct + eosinophils + allergic_rhin
         + age*smoking + I((bmi - 29)^2),
         data = ObsData, family = binomial)$fitted.values

# Positivity: do the two groups overlap?
ggplot(ObsData, aes(x = e, fill = factor(treatment, labels = c("Dupilumab", "Mepolizumab")))) +
  geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
  labs(x = "Probability of receiving mepolizumab", y = "Count", fill = NULL) +
  theme_minimal(base_size = 16)

summary(e[ObsData$treatment == 1])   # mepolizumab patients
summary(e[ObsData$treatment == 0])   # dupilumab patients
sum(e < 0.02 | e > 0.98)             # patients with extreme propensity scores


# =============================================================================
# PART 2: UNADJUSTED (CRUDE) ESTIMATE
# =============================================================================

risk_mepo <- mean(ObsData$exacerbation[ObsData$treatment == 1])
risk_dupi <- mean(ObsData$exacerbation[ObsData$treatment == 0])
c(risk_mepo, risk_dupi) |> round(digits = 3)

RD_crude <- risk_mepo - risk_dupi
RD_crude |> round(digits = 3)   # compare with TRUE_ATE 
TRUE_ATE  


# =============================================================================
# PART 3: INVERSE PROBABILITY OF TREATMENT WEIGHTING (IPTW)
# =============================================================================

# Step 1: propensity scores

# Simple (misspecified) model
e_mis <- glm(treatment ~ age + sex + bmi + smoking + prev_exac + fev1_pct + eosinophils + allergic_rhin,
             data = ObsData, family = binomial)$fitted.values
# Correct model
e <- glm(treatment ~ age + sex + bmi + smoking + prev_exac + fev1_pct + eosinophils + allergic_rhin
         + age*smoking + I((bmi - 29)^2),
         data = ObsData, family = binomial)$fitted.values

# Step 2: weights = 1 / probability of the treatment actually received
omega_mis <- ifelse(ObsData$treatment == 1, 1/e_mis, 1/(1 - e_mis))
omega     <- ifelse(ObsData$treatment == 1, 1/e,     1/(1 - e))
summary(omega_mis)
summary(omega)

# Step 3: check balance in the weighted data (SMD < 0.1)
# We also add the two hidden terms to see whether they are balanced
ObsData$age_x_smoking <- ObsData$age * ObsData$smoking
ObsData$bmi_squared   <- (ObsData$bmi - 29)^2

pseudo_mis <- svydesign(ids = ~ 1, data = ObsData, weights = ~ omega_mis)
tab_mis <- svyCreateTableOne(vars = c("age", "sex", "bmi", "smoking", "prev_exac", "fev1_pct", "eosinophils", "allergic_rhin",
                                      "age_x_smoking", "bmi_squared"),
                             strata = "treatment", data = pseudo_mis, test = FALSE)
print(tab_mis, smd = TRUE)   # simple PS: main effects balanced, but bmi_squared is not

pseudo <- svydesign(ids = ~ 1, data = ObsData, weights = ~ omega)
tab_w <- svyCreateTableOne(vars = c("age", "sex", "bmi", "smoking", "prev_exac", "fev1_pct", "eosinophils", "allergic_rhin",
                                    "age_x_smoking", "bmi_squared"),
                           strata = "treatment", data = pseudo, test = FALSE)
print(tab_w, smd = TRUE)     # correct PS: everything balanced

# Step 4: compare outcomes in the weighted data
# Linear model gives the risk difference; logistic model gives the marginal OR
msm_RD_mis <- lm(exacerbation ~ treatment, weights = omega_mis, data = ObsData)
msm_RD     <- lm(exacerbation ~ treatment, weights = omega,     data = ObsData)

RD_IPW_mis <- msm_RD_mis$coef["treatment"]
RD_IPW     <- msm_RD$coef["treatment"]

c(RD_IPW_mis, RD_IPW) |> round(digits = 3)
TRUE_ATE 
# =============================================================================
# PART 4: G-FORMULA (G-COMPUTATION)
# =============================================================================

# Step 1: outcome models
# Simple (misspecified) model
Q_mis <- glm(exacerbation ~ treatment + age + sex + bmi + smoking + prev_exac + fev1_pct + eosinophils + allergic_rhin,
             data = ObsData, family = binomial)
# Correct model
Q <- glm(exacerbation ~ treatment + age + sex + bmi + smoking + prev_exac + fev1_pct + eosinophils + allergic_rhin
         + age*smoking + I((bmi - 29)^2),
         data = ObsData, family = binomial)

# Step 2: copy the data twice: everyone on mepolizumab, everyone on dupilumab
A1Data <- A0Data <- ObsData
A1Data$treatment <- 1
A0Data$treatment <- 0

# Step 3: predict each patient's risk in both worlds, then average
Y_A1_mis <- predict(Q_mis, A1Data, type = "response")
Y_A0_mis <- predict(Q_mis, A0Data, type = "response")
Y_A1     <- predict(Q,     A1Data, type = "response")
Y_A0     <- predict(Q,     A0Data, type = "response")

pred_A1_mis <- mean(Y_A1_mis); pred_A0_mis <- mean(Y_A0_mis)
pred_A1     <- mean(Y_A1);     pred_A0     <- mean(Y_A0)
c(pred_A1, pred_A0) |> round(digits = 3)   # risk if everyone got mepolizumab / dupilumab

RD_gcomp_mis <- pred_A1_mis - pred_A0_mis
RD_gcomp     <- pred_A1 - pred_A0
c(RD_gcomp_mis, RD_gcomp) |> round(digits = 3)

TRUE_ATE 
# =============================================================================
# PART 5: DOUBLY ROBUST STANDARDIZATION
# =============================================================================
# Same steps as the g-formula, but the outcome model is fitted on the
# IPTW-weighted data. Correct if EITHER the weights OR the outcome model is correct.

# Outcome models fitted with weights
# (R warns about "non-integer #successes": this is expected with weights and can be ignored)
Q_dr         <- glm(exacerbation ~ treatment + age + sex + bmi + smoking + prev_exac + fev1_pct + eosinophils + allergic_rhin
                    + age*smoking + I((bmi - 29)^2),
                    weights = omega, data = ObsData, family = binomial)        # both correct
Q_dr_Qmis    <- glm(exacerbation ~ treatment + age + sex + bmi + smoking + prev_exac + fev1_pct + eosinophils + allergic_rhin,
                    weights = omega, data = ObsData, family = binomial)        # outcome model wrong, weights correct
Q_dr_emis    <- glm(exacerbation ~ treatment + age + sex + bmi + smoking + prev_exac + fev1_pct + eosinophils + allergic_rhin
                    + age*smoking + I((bmi - 29)^2),
                    weights = omega_mis, data = ObsData, family = binomial)    # outcome model correct, weights wrong
Q_dr_mismis  <- glm(exacerbation ~ treatment + age + sex + bmi + smoking + prev_exac + fev1_pct + eosinophils + allergic_rhin,
                    weights = omega_mis, data = ObsData, family = binomial)    # both wrong

# Predict in the two worlds (A1Data and A0Data from Part 5) and average
RD_DRS        <- mean(predict(Q_dr,        A1Data, type = "response")) - mean(predict(Q_dr,        A0Data, type = "response"))
RD_DRS_Qmis   <- mean(predict(Q_dr_Qmis,   A1Data, type = "response")) - mean(predict(Q_dr_Qmis,   A0Data, type = "response"))
RD_DRS_emis   <- mean(predict(Q_dr_emis,   A1Data, type = "response")) - mean(predict(Q_dr_emis,   A0Data, type = "response"))
RD_DRS_mismis <- mean(predict(Q_dr_mismis, A1Data, type = "response")) - mean(predict(Q_dr_mismis, A0Data, type = "response"))

data.frame(both_correct = RD_DRS, Q_wrong = RD_DRS_Qmis, weights_wrong = RD_DRS_emis, both_wrong = RD_DRS_mismis) |> round(digits = 3)

TRUE_ATE 

# =============================================================================
# PART 6: BOOTSTRAP 95% CONFIDENCE INTERVALS (all methods)
# =============================================================================
# Resample patients with replacement and redo every step of Parts 2 to 6.
# Repeat B times; the spread of the B estimates gives the 95% CI.
# (R warns about "non-integer #successes" for the weighted models: expected, can be ignored)

boot_crude <- boot_IPW_mis <- boot_IPW <- boot_gcomp_mis <- boot_gcomp <- c()
boot_DRS_mismis <- boot_DRS_Qmis <- boot_DRS_emis <- boot_DRS <- c()
B <- 50   # use 500 or 1000 in a real analysis; 50 takes about 1 minute

for (i in 1:B) {
  # Draw the bootstrap sample
  db_boot <- ObsData[sample(1:nrow(ObsData), size = nrow(ObsData), replace = TRUE), ]
  
  # Unadjusted
  boot_crude[i] <- mean(db_boot$exacerbation[db_boot$treatment == 1]) - mean(db_boot$exacerbation[db_boot$treatment == 0])
  
  # IPTW: propensity scores, weights, weighted risk difference
  e_mis_b <- glm(treatment ~ age + sex + bmi + smoking + prev_exac + fev1_pct + eosinophils + allergic_rhin,
                 data = db_boot, family = binomial)$fitted.values
  e_b     <- glm(treatment ~ age + sex + bmi + smoking + prev_exac + fev1_pct + eosinophils + allergic_rhin
                 + age*smoking + I((bmi - 29)^2),
                 data = db_boot, family = binomial)$fitted.values
  db_boot$omega_mis_b <- ifelse(db_boot$treatment == 1, 1/e_mis_b, 1/(1 - e_mis_b))
  db_boot$omega_b     <- ifelse(db_boot$treatment == 1, 1/e_b,     1/(1 - e_b))
  
  boot_IPW_mis[i] <- lm(exacerbation ~ treatment, weights = omega_mis_b, data = db_boot)$coef["treatment"]
  boot_IPW[i]     <- lm(exacerbation ~ treatment, weights = omega_b,     data = db_boot)$coef["treatment"]
  
  # Copies of the bootstrap data: everyone on mepolizumab / everyone on dupilumab
  A1Data_b <- A0Data_b <- db_boot
  A1Data_b$treatment <- 1
  A0Data_b$treatment <- 0
  
  # G-formula
  Q_mis_b <- glm(exacerbation ~ treatment + age + sex + bmi + smoking + prev_exac + fev1_pct + eosinophils + allergic_rhin,
                 data = db_boot, family = binomial)
  Q_b     <- glm(exacerbation ~ treatment + age + sex + bmi + smoking + prev_exac + fev1_pct + eosinophils + allergic_rhin
                 + age*smoking + I((bmi - 29)^2),
                 data = db_boot, family = binomial)
  boot_gcomp_mis[i] <- mean(predict(Q_mis_b, A1Data_b, type = "response")) - mean(predict(Q_mis_b, A0Data_b, type = "response"))
  boot_gcomp[i]     <- mean(predict(Q_b,     A1Data_b, type = "response")) - mean(predict(Q_b,     A0Data_b, type = "response"))
  
  # Doubly robust standardization: outcome models fitted with weights
  Q_dr_mismis_b <- glm(exacerbation ~ treatment + age + sex + bmi + smoking + prev_exac + fev1_pct + eosinophils + allergic_rhin,
                       weights = omega_mis_b, data = db_boot, family = binomial)
  Q_dr_Qmis_b   <- glm(exacerbation ~ treatment + age + sex + bmi + smoking + prev_exac + fev1_pct + eosinophils + allergic_rhin,
                       weights = omega_b, data = db_boot, family = binomial)
  Q_dr_emis_b   <- glm(exacerbation ~ treatment + age + sex + bmi + smoking + prev_exac + fev1_pct + eosinophils + allergic_rhin
                       + age*smoking + I((bmi - 29)^2),
                       weights = omega_mis_b, data = db_boot, family = binomial)
  Q_dr_b        <- glm(exacerbation ~ treatment + age + sex + bmi + smoking + prev_exac + fev1_pct + eosinophils + allergic_rhin
                       + age*smoking + I((bmi - 29)^2),
                       weights = omega_b, data = db_boot, family = binomial)
  boot_DRS_mismis[i] <- mean(predict(Q_dr_mismis_b, A1Data_b, type = "response")) - mean(predict(Q_dr_mismis_b, A0Data_b, type = "response"))
  boot_DRS_Qmis[i]   <- mean(predict(Q_dr_Qmis_b,   A1Data_b, type = "response")) - mean(predict(Q_dr_Qmis_b,   A0Data_b, type = "response"))
  boot_DRS_emis[i]   <- mean(predict(Q_dr_emis_b,   A1Data_b, type = "response")) - mean(predict(Q_dr_emis_b,   A0Data_b, type = "response"))
  boot_DRS[i]        <- mean(predict(Q_dr_b,        A1Data_b, type = "response")) - mean(predict(Q_dr_b,        A0Data_b, type = "response"))
}

# Example: standard error and 95% CI for the g-formula (correct model)
sd(boot_gcomp)
quantile(boot_gcomp, probs = c(0.025, 0.975)) |> round(digits = 3)


# =============================================================================
# PART 7: SUMMARY
# =============================================================================

results <- data.frame(
  method   = c("Unadjusted", "IPTW, simple PS", "IPTW, correct PS",
               "G-formula, simple model", "G-formula, correct model",
               "DR, both simple", "DR, outcome model wrong", "DR, weights wrong", "DR, both correct"),
  estimate = c(RD_crude, RD_IPW_mis, RD_IPW, RD_gcomp_mis, RD_gcomp,
               RD_DRS_mismis, RD_DRS_Qmis, RD_DRS_emis, RD_DRS),
  lower    = c(quantile(boot_crude, 0.025), quantile(boot_IPW_mis, 0.025), quantile(boot_IPW, 0.025),
               quantile(boot_gcomp_mis, 0.025), quantile(boot_gcomp, 0.025),
               quantile(boot_DRS_mismis, 0.025), quantile(boot_DRS_Qmis, 0.025), quantile(boot_DRS_emis, 0.025), quantile(boot_DRS, 0.025)),
  upper    = c(quantile(boot_crude, 0.975), quantile(boot_IPW_mis, 0.975), quantile(boot_IPW, 0.975),
               quantile(boot_gcomp_mis, 0.975), quantile(boot_gcomp, 0.975),
               quantile(boot_DRS_mismis, 0.975), quantile(boot_DRS_Qmis, 0.975), quantile(boot_DRS_emis, 0.975), quantile(boot_DRS, 0.975)))
results$percent_of_truth <- round(100 * results$estimate / TRUE_ATE)
results$estimate <- round(results$estimate, 3)
results$lower    <- round(results$lower, 3)
results$upper    <- round(results$upper, 3)
results

ggplot(results, aes(x = estimate, y = factor(method, levels = rev(method)))) +
  geom_vline(xintercept = TRUE_ATE, linetype = "dashed") +
  geom_vline(xintercept = 0, color = "grey") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.25, color = "#2166AC") +
  geom_point(size = 4, color = "#2166AC") +
  labs(x = "Risk difference (mepolizumab - dupilumab), with 95% CI", y = NULL,
       title = "Estimates vs true ATE (dashed line)") +
  theme_minimal(base_size = 16)

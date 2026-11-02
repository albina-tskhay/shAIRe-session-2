# =============================================================================
# 01_simulate_data.R
# Simulate an observational dataset for the causal inference workshop
#
#   treatment    = 1: mepolizumab, 0: dupilumab
#   exacerbation = asthma exacerbation within 12 months (1 = yes, 0 = no)
#
# Because we simulate the data ourselves, we know the true treatment effect.
# =============================================================================

set.seed(6)    # for reproducibility.
n <- 10000     # number of patients

# -----------------------------------------------------------------------------
# Baseline confounders
# -----------------------------------------------------------------------------
age           <- rnorm(n, mean = 48, sd = 12)                     # years
sex           <- rbinom(n, size = 1, prob = 0.55)                 # 1 = female
bmi           <- rnorm(n, mean = 29, sd = 5)                      # kg/m2
smoking       <- rbinom(n, size = 1, prob = 0.25)                 # 1 = current/former smoker
prev_exac     <- pmin(rpois(n, lambda = 1.2), 5)                  # exacerbations in prior year (max 5)
fev1_pct      <- pmax(pmin(rnorm(n, mean = 68, sd = 15), 120), 20) # FEV1 % predicted (20 to 120)
eosinophils   <- round(pmin(rlnorm(n, meanlog = 5.5, sdlog = 0.6), 1500)) # blood eosinophils, cells/uL
allergic_rhin <- rbinom(n, size = 1, prob = 0.40)                 # 1 = allergic rhinitis

# -----------------------------------------------------------------------------
# Treatment: who receives mepolizumab?
# Patients with high eosinophils, more prior exacerbations, and lower FEV1 are
# more likely to receive mepolizumab; patients with allergic rhinitis are more
# likely to receive dupilumab (confounding by indication).
# Two "hidden" terms also play a role: age x smoking and (bmi - 29)^2.
# -----------------------------------------------------------------------------
treatment <- rbinom(n, size = 1, prob = plogis(-2.5 + 0.008*age + 0.10*sex + 0.015*bmi + 0.25*smoking
                                               + 0.50*prev_exac - 0.025*fev1_pct + 0.003*eosinophils
                                               - 0.30*allergic_rhin
                                               + 0.008*age*smoking + 0.020*(bmi - 29)^2))

# -----------------------------------------------------------------------------
# Potential outcomes
# Risk of exacerbation if the patient received mepolizumab (note the -0.40)
# and if the same patient received dupilumab (no -0.40).
# The same two hidden terms also affect the risk of exacerbation.
# -----------------------------------------------------------------------------
risk_Y.1 <- plogis(-1.8 - 0.40 + 0.020*age - 0.15*sex + 0.025*bmi + 0.50*smoking + 0.45*prev_exac
                   - 0.022*fev1_pct + 0.002*eosinophils + 0.20*allergic_rhin
                   + 0.012*age*smoking + 0.004*(bmi - 29)^2)
risk_Y.0 <- plogis(-1.8 - 0    + 0.020*age - 0.15*sex + 0.025*bmi + 0.50*smoking + 0.45*prev_exac
                   - 0.022*fev1_pct + 0.002*eosinophils + 0.20*allergic_rhin
                   + 0.012*age*smoking + 0.004*(bmi - 29)^2)

Y.1 <- rbinom(n, size = 1, prob = risk_Y.1)   # outcome if treated with mepolizumab
Y.0 <- rbinom(n, size = 1, prob = risk_Y.0)   # outcome if treated with dupilumab

# Observed outcome: the potential outcome under the treatment actually received
exacerbation <- Y.1*treatment + Y.0*(1 - treatment)

# -----------------------------------------------------------------------------
# True effects (only known because we simulated the data)
# -----------------------------------------------------------------------------
TRUE_EY.1 <- mean(risk_Y.1)   # risk if everyone received mepolizumab
TRUE_EY.0 <- mean(risk_Y.0)   # risk if everyone received dupilumab
TRUE_ATE  <- TRUE_EY.1 - TRUE_EY.0                                     # risk difference
TRUE_MOR  <- (TRUE_EY.1*(1 - TRUE_EY.0)) / ((1 - TRUE_EY.1)*TRUE_EY.0) # marginal odds ratio

round(c(TRUE_EY.1 = TRUE_EY.1, TRUE_EY.0 = TRUE_EY.0, TRUE_ATE = TRUE_ATE, TRUE_MOR = TRUE_MOR), 3)
mean(treatment)      # proportion on mepolizumab
mean(exacerbation)   # overall risk of exacerbation

# -----------------------------------------------------------------------------
# Save the observed data (without the potential outcomes) and the true effects
# -----------------------------------------------------------------------------
ObsData <- data.frame(id = 1:n, age = round(age, 1), sex, bmi = round(bmi, 1), smoking, prev_exac,
                      fev1_pct = round(fev1_pct, 1), eosinophils, allergic_rhin,
                      treatment, exacerbation)
write.csv(ObsData, "sim_asthma_data.csv", row.names = FALSE)

write.csv(data.frame(TRUE_ATE, TRUE_MOR), "true_effects.csv", row.names = FALSE)
Causal Inference Workshop: IPTW and G-Methods for Respiratory Health
Overview
This workshop introduces inverse probability of treatment weighting (IPTW), the g-formula, and doubly robust estimation through a hands-on R exercise with a simulated asthma dataset. The target audience is respiratory health clinicians with basic R experience. Duration: 30 to 40 minutes.
Files
File	Description
00_installation_guide	R and RStudio installation instructions
00_install_packages.R	Installs and checks all required packages
01_simulate_data.R	Generates the dataset and computes the true causal effects
02_analysis_pipeline.R	Full analysis: descriptive statistics, logistic regression, IPTW, g-formula, AIPW, repeated simulations
sim_asthma_data.csv	Generated dataset (created by 01_simulate_data.R)
true_effects.csv	True causal effects used to evaluate each method
Clinical Scenario
An observational comparison of two biologics for severe asthma:
•	A = 1: mepolizumab (anti-IL-5)
•	A = 0: dupilumab (anti-IL-4Rα)
•	Outcome: asthma exacerbation within 12 months (yes/no)
The dataset mimics confounding by indication. Patients with high blood eosinophils, more prior exacerbations, and lower lung function are more often prescribed mepolizumab, and these same patients are at higher risk of exacerbation. Patients with allergic comorbidity are more often prescribed dupilumab. Without adjustment, mepolizumab looks worse than it is, because it is given to sicker patients. This matches the lecture example, where patients with high eosinophils are more often given mepolizumab.
How the Data Are Generated
Confounders (n = 10,000)
Variable	Distribution	Clinical meaning
age	Normal(48, 12)	Age in years
sex	Bernoulli(0.55)	1 = female
bmi	Normal(29, 5)	Body mass index
smoking	Bernoulli(0.25)	1 = current or former smoker
prev_exac	Poisson(1.2), capped at 5	Exacerbations in the prior year
fev1_pct	Normal(68, 15), clipped to [20, 120]	FEV1 % predicted
eosinophils	Log-normal(5.5, 0.6), capped at 1,500	Blood eosinophils (cells/uL)
allergic_rhin	Bernoulli(0.40)	1 = comorbid allergic rhinitis
Treatment Assignment (true propensity score model)
logit P(mepolizumab | L) = -2.5
  + 0.008 * age
  + 0.10  * sex
  + 0.015 * bmi
  + 0.25  * smoking
  + 0.50  * prev_exac          <-- strong confounder
  - 0.025 * fev1_pct           <-- lower FEV1 -> more likely mepolizumab
  + 0.003 * eosinophils        <-- higher eos -> more likely mepolizumab
  - 0.30  * allergic_rhin      <-- allergic rhinitis -> more likely dupilumab
  + 0.008 * age * smoking      <-- hidden interaction
  + 0.020 * (bmi - 29)^2       <-- hidden quadratic term
Outcome (true outcome model)
logit P(exacerbation | A, L) = -1.8
  - 0.40  * mepolizumab        <-- true causal effect (conditional log-OR)
  + 0.020 * age
  - 0.15  * sex
  + 0.025 * bmi
  + 0.50  * smoking
  + 0.45  * prev_exac
  - 0.022 * fev1_pct
  + 0.002 * eosinophils
  + 0.20  * allergic_rhin
  + 0.012 * age * smoking      <-- hidden interaction
  + 0.004 * (bmi - 29)^2       <-- hidden quadratic term
The two hidden terms appear in both models. A main-effects model (the "simple" model that most analysts fit) leaves them out, so it is misspecified. The "correct" model includes them.
True Effects (seed = 42)
The true ATE is computed by applying the outcome model to every patient under both treatments and averaging: ATE = mean[ P(Y = 1 | A = 1, L) − P(Y = 1 | A = 0, L) ].
Estimand	True value
ATE, risk difference	−0.083 (mepolizumab lowers the 12-month risk by about 8 percentage points)
Risk ratio	0.82
Marginal odds ratio	0.71
Conditional odds ratio	0.67 (= exp(−0.40))
About 25% of patients receive mepolizumab, and the overall exacerbation rate is 45%.
What the Analysis Shows
Results from the main dataset (risk difference, with % of the true effect recovered):
Method	Model	Estimate	% of true effect
Unadjusted	none	+0.078	−94% (wrong direction)
IPTW	simple PS	−0.058	70%
IPTW	correct PS	−0.090	108%
G-formula	simple outcome model	−0.057	69%
G-formula	correct outcome model	−0.085	102%
AIPW	both simple	−0.058	70%
AIPW	only PS correct	−0.091	110%
AIPW	both correct	−0.092	111%
Averaged over 100 repeated simulations (Part 8 of the pipeline), the correctly specified IPTW, g-formula, and AIPW recover 100% of the true effect on average. The simple (misspecified) versions recover 71% to 74%, and the unadjusted comparison points in the wrong direction. In any single dataset, correctly specified estimates differ from the truth by a little because of random noise; the repeated simulations show that this difference averages out.
Key messages:
1.	Unadjusted comparison: confounding by indication makes mepolizumab look harmful.
2.	Logistic regression: the adjusted odds ratio it reports is a conditional effect, a different estimand from the marginal effect. A marginal risk difference can be obtained from the same model by standardization (predicting everyone under both treatments), but that calculation is the g-formula, so it equals the "G-formula, simple model" row. Its main-effects model is also misspecified.
3.	IPTW, g-formula, AIPW (simple models): remove most of the confounding bias, but some remains because the models miss the hidden terms.
4.	IPTW, g-formula, AIPW (correct models): recover the true effect on average.
5.	AIPW with only one correct model: still recovers the truth, which is the double robustness property.
6.	Balance checks: the simple PS balances all main effects (SMD < 0.1), yet BMI squared remains imbalanced. Balance should also be checked on plausible nonlinear terms.
How to Run
r
# Step 1: install packages (once)
source("00_install_packages.R")

# Step 2: generate the data
source("01_simulate_data.R")

# Step 3: run the analysis (about 2 to 3 minutes including the bootstrap and Part 8)
source("02_analysis_pipeline.R")


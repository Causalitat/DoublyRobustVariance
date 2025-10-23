# Program: ucla_v2.r
# Purpose: This program applies the geex ZINB and Hurdle functions to the UCLA NMES1988 dataset.

# Load required libraries
library(AER)
library(pscl)
library(geex)
library(MASS)

# Load the custom functions
source("example/R/02_ZINB_Functions.R")
source("example/R/03_Hurdle_Functions.R")

# Load the data in R from AER
data("NMES1988")

# Select variables used in the analysis
dat <- NMES1988[, c(1, 6:8, 13, 15, 18)]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Causal Question:
# What is the Average Treatment Effect (ATE) of private insurance on the number
# of doctor visits among people aged 66 and over?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Prepare the data for causal analysis
# Outcome (Y): visits
dat$Y <- dat$visits
# Exposure (X): insurance (1 = yes, 0 = no)
dat$X <- as.numeric(dat$insurance) - 1
# Covariates (Z): hospital, health, chronic, gender, school
dat <- dat[complete.cases(dat), ]

# Define the model formulas
propensity_formula <- X ~ hospital + health + chronic + gender + school
outcome_formula_zinb <- Y ~ X + hospital + health + chronic + gender + school | hospital + health + chronic + gender + school
outcome_formula_hurdle <- Y ~ X + hospital + health + chronic + gender + school | hospital + health + chronic + gender + school

# Apply the three geex ZINB estimators
wtd_results_zinb <- geex_WTD_ZINB(dat, propensity_formula, outcome_formula_zinb)
gf_results_zinb <- geex_GF_ZINB(dat, outcome_formula_zinb)
aipw_results_zinb <- geex_AIPW_ZINB(dat, propensity_formula, outcome_formula_zinb)

# Apply the three geex Hurdle estimators
wtd_results_hurdle <- geex_WTD_hurdle(dat, propensity_formula, outcome_formula_hurdle)
gf_results_hurdle <- geex_GF_hurdle(dat, outcome_formula_hurdle)
aipw_results_hurdle <- geex_AIPW_hurdle(dat, propensity_formula, outcome_formula_hurdle)

# Print the results
print("--- ZINB Estimators ---")
print("Weighted ZINB Estimator:")
print(wtd_results_zinb)
print("G-Formula ZINB Estimator:")
print(gf_results_zinb)
print("AIPW ZINB Estimator:")
print(aipw_results_zinb)

print("--- Hurdle Estimators ---")
print("Weighted Hurdle Estimator:")
print(wtd_results_hurdle)
print("G-Formula Hurdle Estimator:")
print(gf_results_hurdle)
print("AIPW Hurdle Estimator:")
print(aipw_results_hurdle)

# Program: 11_ZINB_Interaction_Check.R
# Purpose: Verify that geex_AIPW_ZINB_Int correctly handles interactions between Treatment (X) and Covariates.

rm(list=ls())

library(geex)
library(pscl)
library(numDeriv)
library(MASS)

# Source functions (Path relative to project root)
# Source functions (Path relative to project root)
source("example/R/02_ZINB_Functions.R")
source("example/R/04_Hurdle_NegBin_Functions_V2.R")

set.seed(42)
n <- 200

# Generate Data with Interaction Effect
# Truth: Treatment effect depends on Z
Z <- rnorm(n)
X_prob <- plogis(0.5 * Z)
X <- rbinom(n, 1, X_prob)

# Count part: log(mu) = 1 + 0.5*X + 0.5*Z + 1.0*(X*Z)
# Interaction term 1.0 means effect of X varies by Z.
mu_true <- exp(1 + 0.5*X + 0.5*Z + 1.0*X*Z) 

# Zero part: logit(phi) = -1 + 0.2*X + 0.2*Z (No interaction here for simplicity, but wrapper adds it)
phi_true <- plogis(-1 + 0.2*X + 0.2*Z)

# Simulate ZINB Outcome
Y <- numeric(n)
is_zero <- rbinom(n, 1, phi_true)
counts <- rnegbin(n, mu = mu_true, theta = 2)
Y <- ifelse(is_zero==1, 0, counts)

data <- data.frame(Y=Y, X=X, Z=Z)

cat("Data Generated.\n")
cat("----------------------------------------------------------------\n")
cat("1. Running ZINB Interaction Model (Doubly Robust w/ AIPW)...\n")
cat("----------------------------------------------------------------\n")

# Formulas
ps_formula <- X ~ Z
cov_formula <- ~ Z # Wrapper will expand to X*Z

# Run ZINB Wrapper
tryCatch({
  res <- geex_AIPW_ZINB_Int(data, ps_formula, cov_formula)
  print(res)
  cat("\nSuccess: ZINB Interaction model ran without error.\n")
}, error = function(e) {
  cat("\nError running interaction model:\n")
  print(e)
})

cat("\n----------------------------------------------------------------\n")
cat("2. Running Hurdle NB V2 Interaction Model (Weighted)...\n")
cat("   (Using 04_Hurdle_NegBin_Functions_V2.R)\n")
cat("----------------------------------------------------------------\n")

tryCatch({
  # geex_WTD_Hurdle_NB_Int is defined in V2 file and uses the V2 estimator (estimated theta)
  res_hurdle <- geex_WTD_Hurdle_NB_Int(data, ps_formula, cov_formula)
  print(res_hurdle$result)
  cat("Theta Est:", res_hurdle$theta_est, "\n")
  cat("\nSuccess: Hurdle V2 Interaction model ran without error.\n")
}, error = function(e) {
  cat("\nError running Hurdle V2 model:\n")
  print(e)
})

# Program: 06_Comparison_Analysis.R
# Purpose: Compare Manual ZTNB Est (04) vs Refactored PSCL ZTNB Est (03)

library(geex)
library(boot)
library(rootSolve)
library(MASS)
library(pscl)
library(numDeriv)

# Source both function files
# Note: Ensure no function name clashes. 
# 04: geex_Hurdle_NB
# 03: geex_AIPW_hurdle
source("04_Hurdle_NegBin_Functions.R")
source("03_Hurdle_Functions.R")

set.seed(12345)
n <- 1000

# ------------------------------------------------------------------------------
# 1. Generate Synthetic Data (Hurdle Negative Binomial)
# ------------------------------------------------------------------------------
Z1 <- rnorm(n)
Z2 <- rbinom(n, 1, 0.5)
p_X <- plogis(0.5 * Z1 - 0.5 * Z2)
X <- rbinom(n, 1, p_X)

# Outcome Models
# Zero part (Logistic)
logit_p_pos <- -0.5 + 1.0 * X + 0.5 * Z1 + 0.5 * Z2
p_pos <- plogis(logit_p_pos)
is_positive <- rbinom(n, 1, p_pos)

# Count part (ZTNB)
log_mu <- 1.0 + 0.5 * X + 0.2 * Z1
mu <- exp(log_mu)
theta_true <- 2.0 # Dispersion

Y_count <- rep(0, n)
for(i in 1:n) {
  if(is_positive[i] == 1) {
    val <- 0
    while(val == 0) {
      val <- rnbinom(1, size = theta_true, mu = mu[i])
    }
    Y_count[i] <- val
  }
}
Y <- Y_count
data <- data.frame(Y=Y, X=X, Z1=Z1, Z2=Z2)

# Formulas
# 04 style
ps_formula <- X ~ Z1 + Z2
zero_formula <- Y_bin ~ X + Z1 + Z2
count_formula <- Y ~ X + Z1

# 03 style (pscl style)
# outcome_formula: count_part | zero_part
# Note: pscl hurdle uses same predictors for count unless specified with |
# LHS is Y
outcome_formula_pscl <- Y ~ X + Z1 | X + Z1 + Z2

cat("Running Comparison...\n")

# ------------------------------------------------------------------------------
# 2. Run Manual Estimator (04)
# ------------------------------------------------------------------------------
cat("\n[Method 1] Manual ZTNB (04_Hurdle_NegBin_Functions.R)...\n")
tryCatch({
  res_04 <- geex_Hurdle_NB(data, ps_formula, zero_formula, count_formula)
  print(res_04$result)
  cat("Theta Used:", res_04$theta_used, "\n")
}, error = function(e) print(e))

# ------------------------------------------------------------------------------
# 3. Run Refactored PSCL Estimator (03)
# ------------------------------------------------------------------------------
cat("\n[Method 2] Refactored PSCL (03_Hurdle_Functions.R) - AIPW...\n")
tryCatch({
  res_03 <- geex_AIPW_hurdle(data, ps_formula, outcome_formula_pscl, dist = "negbin")
  print(res_03)
}, error = function(e) print(e))

cat("\nComparison Complete.\n")

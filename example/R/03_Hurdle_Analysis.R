# Program: 03_Hurdle_Analysis_NB.R
# Purpose: Compare Poisson Hurdle vs NB Hurdle on Overdispersed Data

library(geex)
library(boot)
library(rootSolve)
library(MASS)

# Source functions
source("02_Hurdle_Functions.R")
source("04_Hurdle_NegBin_Functions.R")

set.seed(12345)
n <- 2000

# ------------------------------------------------------------------------------
# 1. Generate Synthetic Data (Hurdle Negative Binomial)
# ------------------------------------------------------------------------------
# Covariates
Z1 <- rnorm(n)
Z2 <- rbinom(n, 1, 0.5)

p_X <- plogis(0.5 * Z1 - 0.5 * Z2)
X <- rbinom(n, 1, p_X)

# Outcomes
# Part 1: Zero vs Non-Zero (Binary)
logit_p_pos <- -0.5 + 1.0 * X + 0.5 * Z1 + 0.5 * Z2
p_pos <- plogis(logit_p_pos)
is_positive <- rbinom(n, 1, p_pos)

# Part 2: Count given Positive (Zero-Truncated NB)
log_mu <- 1.0 + 0.5 * X + 0.2 * Z1
mu <- exp(log_mu)
theta_true <- 1.5 # Dispersion parameter (Overdispersion)

Y_count <- rep(0, n)
for(i in 1:n) {
  if(is_positive[i] == 1) {
    val <- 0
    while(val == 0) {
      # Use base R rnbinom
      val <- rnbinom(1, size = theta_true, mu = mu[i])
    }
    Y_count[i] <- val
  }
}
Y <- Y_count
data <- data.frame(Y=Y, X=X, Z1=Z1, Z2=Z2)

cat("Data Summary:\n")
cat("Mean:", mean(data$Y), "Var:", var(data$Y), "Dispersion Index:", var(data$Y)/mean(data$Y), "\n")
cat("Proportion Zeros:", mean(data$Y==0), "\n\n")

# Formulas
ps_formula <- X ~ Z1 + Z2
zero_formula <- Y_bin ~ X + Z1 + Z2 
count_formula <- Y ~ X + Z1

# ------------------------------------------------------------------------------
# 2. Run Estimators
# ------------------------------------------------------------------------------

cat("------------------------------------------------------\n")
cat("1. Running Poisson Hurdle (ZTP) Estimator...\n")
# Ideally, this should have larger SE or be biased if overdispersion matters for variance?
# Usually, point estimates might be consistent, but SE underestimated.
tryCatch({
  res_pois <- geex_Hurdle(data, ps_formula, zero_formula, count_formula)
  print(res_pois)
}, error = function(e) print(e))

cat("\n------------------------------------------------------\n")
cat("2. Running Negative Binomial Hurdle (ZTNB) Estimator...\n")
tryCatch({
  res_nb <- geex_Hurdle_NB(data, ps_formula, zero_formula, count_formula)
  print(res_nb$result)
  cat("Theta Used:", res_nb$theta_used, "\n")
}, error = function(e) print(e))

cat("\n------------------------------------------------------\n")
cat("Comparison Complete.\n")

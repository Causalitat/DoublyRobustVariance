# Program: 06_Hurdle_NB_Comparison.R
# Purpose: Compare V1 (Fixed Theta) and V2 (Estimated Theta) Hurdle NB Estimators
#          to see the impact on Standard Errors.
#          Uses separate environments to avoid function name collisions.

library(geex)
library(boot)
library(rootSolve)
library(MASS)

# ------------------------------------------------------------------------------
# 1. Load Functions into Isolated Environments
# ------------------------------------------------------------------------------

# Create Environment for V1
env_v1 <- new.env()
# We need helper functions in there too if they are not global
sys.source("02_Hurdle_Functions.R", envir = env_v1)
sys.source("04_Hurdle_NegBin_Functions.R", envir = env_v1)

# Create Environment for V2
env_v2 <- new.env()
sys.source("02_Hurdle_Functions.R", envir = env_v2)
sys.source("04_Hurdle_NegBin_Functions_V2.R", envir = env_v2)

# ------------------------------------------------------------------------------
# 2. Generate Synthetic Data
# ------------------------------------------------------------------------------
set.seed(999) 
n <- 2000

# Covariates
Z1 <- rnorm(n)
Z2 <- rbinom(n, 1, 0.5)

p_X <- plogis(0.5 * Z1 - 0.5 * Z2)
X <- rbinom(n, 1, p_X)

# Outcomes
# Part 1: Zero vs Non-Zero (Binary)
logit_p_pos <- -0.5 + 1.2 * X + 0.5 * Z1 + 0.5 * Z2
p_pos <- plogis(logit_p_pos)
is_positive <- rbinom(n, 1, p_pos)

# Part 2: Count given Positive (Zero-Truncated NB)
log_mu <- 1.0 + 0.5 * X + 0.2 * Z1
mu <- exp(log_mu)
theta_true <- 0.8 

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

cat("Data Summary (N=", n, "):\n")
cat("Mean:", round(mean(data$Y),3), "Var:", round(var(data$Y),3), "\n")
cat("Dispersion Index (Var/Mean):", round(var(data$Y)/mean(data$Y),3), "\n")
cat("Proportion Zeros:", round(mean(data$Y==0),3), "\n\n")

# Formulas
ps_formula <- X ~ Z1 + Z2
zero_formula <- Y_bin ~ X + Z1 + Z2 
count_formula <- Y ~ X + Z1

# ------------------------------------------------------------------------------
# 3. Compare Estimators
# ------------------------------------------------------------------------------

cat("----------------------------------------------------------------\n")
cat("Comparison of Unweighted (Doubly Robust) Estimators\n")
cat("----------------------------------------------------------------\n")

# Run V1
cat("Running V1 (Fixed Theta)...\n")
res_v1 <- tryCatch({
  env_v1$geex_Hurdle_NB(data, ps_formula, zero_formula, count_formula)
}, error = function(e) { print(e); list(result=data.frame(ATE=NA, SE=NA, Type="V1-Error")) })

# Run V2
cat("Running V2 (Estimated Theta)...\n")
res_v2 <- tryCatch({
  env_v2$geex_Hurdle_NB(data, ps_formula, zero_formula, count_formula)
}, error = function(e) { print(e); list(result=data.frame(ATE=NA, SE=NA, Type="V2-Error")) })

# Display
df_dr <- rbind(res_v1$result, res_v2$result)
print(df_dr)

cat("\n----------------------------------------------------------------\n")
cat("Comparison of Weighted (IPW) Estimators\n")
cat("----------------------------------------------------------------\n")

# Run V1
cat("Running V1 Weighted...\n")
res_w_v1 <- tryCatch({
  env_v1$geex_WTD_Hurdle_NB(data, ps_formula, zero_formula, count_formula)
}, error = function(e) { print(e); list(result=data.frame(ATE=NA, SE=NA, Type="V1-Error")) })

# Run V2
cat("Running V2 Weighted...\n")
res_w_v2 <- tryCatch({
  env_v2$geex_WTD_Hurdle_NB(data, ps_formula, zero_formula, count_formula)
}, error = function(e) { print(e); list(result=data.frame(ATE=NA, SE=NA, Type="V2-Error")) })

df_w <- rbind(res_w_v1$result, res_w_v2$result)
print(df_w)

cat("\n----------------------------------------------------------------\n")
cat("Conclusion on Variance Inflation:\n")
se_v1 <- df_dr$SE[1]
se_v2 <- df_dr$SE[2]
if(!is.na(se_v1) && !is.na(se_v2)) {
  cat("DR SE Increase:", round((se_v2 - se_v1)/se_v1 * 100, 2), "%\n")
}

se_w1 <- df_w$SE[1]
se_w2 <- df_w$SE[2]
if(!is.na(se_w1) && !is.na(se_w2)) {
  cat("IPW SE Increase:", round((se_w2 - se_w1)/se_w1 * 100, 2), "%\n")
}

# Program: 08_Interaction_Check_V2.R
# Purpose: Verify that V2 estimators handle Interaction Terms correctly.

library(geex)
library(MASS)
library(boot)

# Load V2 functions into a clean environment to be sure
env_v2 <- new.env()
# We need helper functions
if(file.exists("02_Hurdle_Functions.R")) sys.source("02_Hurdle_Functions.R", envir = env_v2)
sys.source("04_Hurdle_NegBin_Functions_V2.R", envir = env_v2)

# ------------------------------------------------------------------------------
# 1. Generate Data with Strong Treatment Effect Heterogeneity
# ------------------------------------------------------------------------------
set.seed(2025)
n <- 2500

# Confounder / Modifier
Z <- rnorm(n, mean = 2, sd = 1) 

# Treatment Assignment
p <- plogis(-1 + 0.5 * Z)
X <- rbinom(n, 1, p)

# Outcome Generation
# Interaction: Treatment effect depends heavily on Z.
# Y ~ NegBin(mu, theta)
# log(mu) = 1 + 0.5*Z + 0.0*X + 0.8 * (X * Z)
# Interpretation: Main effect of X is 0, but interaction is 0.8.
# So for Z=2 (mean), effect of X is 1.6.

log_mu <- 1.0 + 0.5 * Z + 0.0 * X + 0.8 * (X * Z)
mu <- exp(log_mu)
theta_true <- 2.0

Y <- rnbinom(n, size = theta_true, mu = mu)
# Note: This is simpler than Hurdle (just NB), but Hurdle functions work on it too (outcomes > 0).
# Let's ensure no zeros to test the Count Part specifically, 
# OR let it naturally have zeros (standard NB has zeros).
# The Hurdle estimator treats 0s via the Zero part.
# If we generate from NB, we have "structural" zeros and "sampling" zeros.
# The Hurdle model separates them. It works fine.

data <- data.frame(Y=Y, X=X, Z=Z)

cat("Data Summary:\n")
cat("Avg Y|X=0:", mean(data$Y[data$X==0]), "\n")
cat("Avg Y|X=1:", mean(data$Y[data$X==1]), "\n\n")

# ------------------------------------------------------------------------------
# 2. Define Formulas WITH Interactions
# ------------------------------------------------------------------------------

# Propensity
ps_f <- X ~ Z

# Hurdle Parts with Interaction X*Z
# Note: We must include the interaction in the formula explicitly
# Count part: Y ~ X * Z  -> expands to X + Z + X:Z
count_f_int <- Y ~ X * Z
zero_f_int  <- as.numeric(Y > 0) ~ X * Z

# ------------------------------------------------------------------------------
# 3. Run V2 DR Estimator with Interactions
# ------------------------------------------------------------------------------

cat("Running V2 DR Estimator with Interactions...\n")

res_dr <- tryCatch({
  env_v2$geex_Hurdle_NB(data, ps_f, zero_f_int, count_f_int)
}, error = function(e) { print(e); NULL })

if(!is.null(res_dr)) {
  print(res_dr$result)
  cat("Theta Estimated:", res_dr$theta_est, "\n")
  cat("(Truth was 2.0)\n")
}

# ------------------------------------------------------------------------------
# 4. Check against Truth (Approximation)
# ------------------------------------------------------------------------------
# True Marginal Causal Effect
# E[Y(1)] = E[exp(1 + 0.5Z + 0.0 + 0.8Z)] = E[exp(1 + 1.3Z)]
# E[Y(0)] = E[exp(1 + 0.5Z)]
# Since Z ~ N(2, 1), E[exp(a + bZ)] = exp(a + 2b + 0.5*b^2)

mean_Y1 <- exp(1 + 1.3*2 + 0.5 * 1.3^2) # exp(1 + 2.6 + 0.845) = exp(4.445)
mean_Y0 <- exp(1 + 0.5*2 + 0.5 * 0.5^2) # exp(1 + 1 + 0.125) = exp(2.125)

true_ate <- mean_Y1 - mean_Y0
cat("\nTrue Analytical ATE (approx):", true_ate, "\n")

cat("\nAnalysis: If the model supports interaction, the ATE should be close to", round(true_ate,1), "\n")

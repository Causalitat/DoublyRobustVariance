
# Test M-estimation Skill
# This script implements a simple IPW estimator using geex to verify the skill instructions.

library(geex)
library(numDeriv)

# 1. Simulate Data
set.seed(123)
n <- 1000
W <- rnorm(n)
A <- rbinom(n, 1, plogis(0.5 * W))
Y <- 2 * A + W + rnorm(n)
data <- data.frame(W = W, A = A, Y = Y)

# 2. Fit Nuisance Model (Propensity Score)
ps_model <- glm(A ~ W, data = data, family = binomial)

# 3. Define Estimating Function (IPW for mu1 = E[Y^1])
# Target: mu1 = E[Y^1]
# Estimating Equation: (A * Y) / e(W) - mu1 = 0

ipw_estFUN <- function(data, models){
  A <- data$A
  Y <- data$Y
  W <- data$W
  
  # PS Model components
  Xe <- grab_design_matrix(data = data, rhs_formula = grab_fixed_formula(models$e))
  e_pos <- 1:ncol(Xe)
  e_scores <- grab_psiFUN(models$e, data)
  
  function(theta){
    # Recover parameters
    beta_ps <- theta[e_pos]
    mu_1 <- theta[max(e_pos) + 1]
    
    # Calculate scores
    e <- plogis(Xe %*% beta_ps)
    # Calculate scores
    e <- plogis(Xe %*% beta_ps)
    s_ps <- e_scores(beta_ps)
    s_ipw <- (A * Y) / e - mu_1
    
    c(s_ps, s_ipw)
  }
}

# 4. Initial Values
init_mu1 <- mean(data$Y[data$A==1] / fitted(ps_model)[data$A==1])
roots <- c(coef(ps_model), init_mu1)

# 5. Run m_estimate
results <- m_estimate(
  estFUN = ipw_estFUN,
  data = data,
  roots = roots,
  outer_args = list(models = list(e = ps_model)),
  compute_roots = FALSE
)

# 6. Output Results
print("M-estimation Results:")
print(results)

# Extract SE
n_params <- length(results@estimates)
se <- sqrt(results@vcov[n_params, n_params])
print(paste("Estimated SE for mu1:", se))

# Compare with naive calculation (ignoring PS uncertainty)
weights <- 1 / fitted(ps_model)
naive_se <- sd(data$Y[data$A==1] * weights[data$A==1]) / sqrt(sum(data$A==1))
# Note: This isn't exactly the naive SE but close enough for a sanity check that M-estimation doesn't explode.
print(paste("Naive check SE (approx):", naive_se))

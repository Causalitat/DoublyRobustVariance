# Program: 02_Hurdle_Functions.R
# Purpose: Extension of DR estimators for Zero-Inflated/Hurdle Count Outcomes
# Conceptual Model:
# E[Y|X,Z] = Pr(Y>0|X,Z) * E[Y|Y>0, X,Z]
# Part 1 (Zero): Pr(Y>0) ~ Logistic(X, Z)
# Part 2 (Count): Y|Y>0 ~ ZeroTruncatedPoisson(lambda(X,Z))

library(geex)
library(rootSolve) # For gradient helper if needed, but we do manual scores

# Helper to grab design matrix
gdm <- function(data, formula) {
  model.matrix(formula, data = data)
}

# ------------------------------------------------------------------------------
# 1. Component Score Functions
# ------------------------------------------------------------------------------

# Score for Zero-Truncated Poisson (ZTP)
# U(beta) = X * (Y - E[Y|Y>0])
score_ztp <- function(Y, X_mat, beta) {
  lambda <- exp(X_mat %*% beta)
  denom <- 1 - exp(-lambda)
  mu_ztp <- lambda / denom
  
  resid <- (Y - mu_ztp)
  resid[Y == 0] <- 0
  
  # Use sweep for safe element-wise multiplication
  # sweep(MARGIN=1) signals row-wise multiplication, which is what we want
  # (Each row of X_mat multiplied by corresponding resid element)
  return(sweep(X_mat, 1, as.vector(resid), "*"))
}

# Score for Standard Logistic (for the Zero/Hurdle part)
score_binary <- function(Y_binary, X_mat, beta) {
  probs <- plogis(X_mat %*% beta)
  resid <- Y_binary - probs
  return(sweep(X_mat, 1, as.vector(resid), "*"))
}

# ------------------------------------------------------------------------------
# 2. Estimating Function
# ------------------------------------------------------------------------------

estfun_Hurdle_PI <- function(data, models) {
  
  X <- data$X
  Y <- data$Y
  Y_bin <- as.numeric(Y > 0)
  
  # Extract formulas from the provided glm/dummy objects
  formula_e <- grab_fixed_formula(models$e)
  formula_z <- grab_fixed_formula(models$zero)
  formula_c <- grab_fixed_formula(models$count)
  
  # Design Matrices
  Xe <- gdm(data, formula_e)
  Xz <- gdm(data, formula_z)
  Xc <- gdm(data, formula_c)
  
  # Design Matrices for Counterfactuals (X=0, X=1)
  data0 <- data; data0$X <- 0
  data1 <- data; data1$X <- 1
  
  Xz0 <- gdm(data0, formula_z)
  Xz1 <- gdm(data1, formula_z)
  Xc0 <- gdm(data0, formula_c)
  Xc1 <- gdm(data1, formula_c)
  
  # Parameter indices
  p_e <- ncol(Xe)
  p_z <- ncol(Xz)
  p_c <- ncol(Xc)
  
  idx_e <- 1:p_e
  idx_z <- (p_e + 1):(p_e + p_z)
  idx_c <- (p_e + p_z + 1):(p_e + p_z + p_c)
  idx_cm1 <- p_e + p_z + p_c + 1
  idx_cm0 <- p_e + p_z + p_c + 2
  idx_ate <- p_e + p_z + p_c + 3
  
  # Pre-calculate propensity scores scores function if using geex helper
  # But here we can just write the score manually for consistency/speed
  # Score for PS (Logistic): X * (A - p)
  
  function(theta) {
    # 1. Unpack Parameters
    # theta is a single vector of length p_e + p_z + p_c + 3
    # Ensure dimensions match
    
    beta_e <- theta[idx_e]
    beta_z <- theta[idx_z]
    beta_c <- theta[idx_c]
    mu1    <- theta[idx_cm1]
    mu0    <- theta[idx_cm0]
    ate    <- theta[idx_ate]
    
    # 2. Compute Propensity Scores
    ps <- plogis(Xe %*% beta_e)
    # Numerical stability for PS
    ps <- pmax(pmin(ps, 1 - 1e-6), 1e-6)
    score_e <- sweep(Xe, 1, as.vector(X - ps), "*")
    
    # 3. Compute Hurdle Model Scores
    # Part 1: Binary (Zero vs Non-Zero)
    score_z <- score_binary(Y_bin, Xz, beta_z)
    
    # Part 2: Count (ZTP) with Lambda safeguards
    lambda_c_s <- exp(Xc %*% beta_c) # raw lambda
    lambda_c_s <- pmax(lambda_c_s, 1e-6)
    
    denom_c <- 1 - exp(-lambda_c_s)
    mu_ztp_c <- lambda_c_s / denom_c
    
    resid_c <- (Y - mu_ztp_c)
    resid_c[Y == 0] <- 0
    score_c <- sweep(Xc, 1, as.vector(resid_c), "*")
    
    # 4. Compute Predictions for DR Estimator
    
    # Under X=1
    prob_pos_1 <- plogis(Xz1 %*% beta_z)
    lambda_1   <- exp(Xc1 %*% beta_c)
    lambda_1   <- pmax(lambda_1, 1e-6)
    mu_ztp_1   <- lambda_1 / (1 - exp(-lambda_1))
    m1         <- prob_pos_1 * mu_ztp_1
    
    # Under X=0
    prob_pos_0 <- plogis(Xz0 %*% beta_z)
    lambda_0   <- exp(Xc0 %*% beta_c)
    lambda_0   <- pmax(lambda_0, 1e-6)
    mu_ztp_0   <- lambda_0 / (1 - exp(-lambda_0))
    m0         <- prob_pos_0 * mu_ztp_0
    
    # 5. Influence Function Parts
    psi_cm1 <- (X * Y - (X - ps) * m1) / ps - mu1
    psi_cm0 <- ((1 - X) * Y + (X - ps) * m0) / (1 - ps) - mu0
    psi_ate <- mu1 - mu0 - ate
    
    # Stack equations
    scores <- cbind(score_e, score_z, score_c, psi_cm1, psi_cm0, psi_ate)
    
    return(unname(scores))
  }
}

# ------------------------------------------------------------------------------
# 3. Wrapper Function
# ------------------------------------------------------------------------------

geex_Hurdle <- function(data, 
                        propensity_formula, 
                        zero_formula, 
                        count_formula,
                        initial_roots = NULL) {
  
  # Fit initial models to get roots if not provided
  
  # 1. Propensity
  m_ps <- glm(propensity_formula, data = data, family = binomial)
  
  # 2. Zero Part (Logistic on I(Y>0))
  data$Y_bin <- as.numeric(data$Y > 0)
  m_zero <- glm(zero_formula, data = data, family = binomial)
  
  # 3. Count Part (ZTP)
  # We use a standard Poisson on Y>0 subset as a starting guess.
  # Ideally, use a ZTP library, but standard Poisson is usually close enough for start values.
  # Or use simple optim to find MLE for ZTP.
  data_pos <- data[data$Y > 0, ]
  
  # Custom MLE for ZTP starter
  nll_ztp <- function(beta, y, X) {
    lambda <- exp(X %*% beta)
    loglik <- -lambda + y * log(lambda) - log(1 - exp(-lambda)) # dropping constants
    return(-sum(loglik))
  }
  
  start_beta_c <- coef(glm(count_formula, data = data_pos, family = poisson))
  X_pos <- model.matrix(count_formula, data_pos)
  Y_pos <- data_pos$Y
  
  # Optimize to get better ZTP roots
  opt_c <- optim(par = start_beta_c, fn = nll_ztp, y = Y_pos, X = X_pos, method = "BFGS")
  root_c <- opt_c$par
  
  # 4. Initial Claims for Means
  # Predict on full data
  Xz_full <- model.matrix(zero_formula, data)
  Xc_full <- model.matrix(count_formula, data)
  
  p_pos <- plogis(Xz_full %*% coef(m_zero))
  lam   <- exp(Xc_full %*% root_c)
  mu_ztp <- lam / (1 - exp(-lam))
  fitted_Y <- p_pos * mu_ztp
  
  # Simple DR initial guess (or just plug-in)
  # We'll just pass 0s or simple means if we trust geex to find root, 
  # but geex requires good roots.
  
  # Let's calculate the plugin estimates for roots
  # (Re-using logic from estfun would be cleaner, but just quick calc here)
  dat0 <- data; dat0$X <- 0
  dat1 <- data; dat1$X <- 1
  
  Xz0 <- model.matrix(zero_formula, dat0)
  Xc0 <- model.matrix(count_formula, dat0)
  p0 <- plogis(Xz0 %*% coef(m_zero))
  l0 <- exp(Xc0 %*% root_c)
  m0 <- mean(p0 * (l0 / (1 - exp(-l0))))
  
  Xz1 <- model.matrix(zero_formula, dat1)
  Xc1 <- model.matrix(count_formula, dat1)
  p1 <- plogis(Xz1 %*% coef(m_zero))
  l1 <- exp(Xc1 %*% root_c)
  m1 <- mean(p1 * (l1 / (1 - exp(-l1))))
  
  roots <- c(coef(m_ps), coef(m_zero), root_c, m1, m0, m1 - m0)
  roots <- unname(roots)
  
  # Setup models list for estfun (passed mainly for formula extraction)
  # We construct dummy glm objects for the count part just to carry the formula
  models <- list(e = m_ps, 
                 zero = m_zero, 
                 count = glm(count_formula, data=data, family=poisson))
                 
  geex_results <- m_estimate(
    estFUN = estfun_Hurdle_PI,
    data = data,
    roots = roots,
    compute_roots = FALSE, 
    outer_args = list(models = models)
  )
  
  # Extract results
  n_params <- length(geex_results@estimates)
  est <- geex_results@estimates[n_params]
  se  <- sqrt(geex_results@vcov[n_params, n_params])
  
  return(data.frame(ATE = est, SE = se, Type = "Hurdle-PI"))
}

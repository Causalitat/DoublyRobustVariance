# Program: 04_Hurdle_NegBin_Functions.R
# Purpose: Extension of DR estimators for Hurdle Models with Negative Binomial Counts
# Conceptual Model:
# Part 1 (Zero): Pr(Y>0) ~ Logistic(X, Z)
# Part 2 (Count): Y|Y>0 ~ ZeroTruncatedNegBin(mu, theta)

library(geex)
library(MASS) # For glm.nb to get initial theta

# Helper to grab design matrix (same as before)
gdm <- function(data, formula) {
  model.matrix(formula, data = data)
}

# ------------------------------------------------------------------------------
# 1. Component Score Functions
# ------------------------------------------------------------------------------

# Score for Zero-Truncated Negative Binomial (ZTNB)
# Derived Analytic Score for beta coefficients (assuming fixed theta)
# U(beta) = X^T * W * (Y - E[Y|Y>0])
# Where E[Y|Y>0] = mu / (1 - P0)
# P0 = (1 + mu/theta)^(-theta)
# W = 1 / (1 + mu/theta)   <-- weight from NB variance structure
score_ztnb <- function(Y, X_mat, beta, theta) {
  mu <- exp(X_mat %*% beta)
  
  # Numerical safeguards for mu
  mu <- pmax(mu, 1e-6)
  
  # Calculate P0 (Prob of 0 in standard NB)
  # P0 = (theta / (theta + mu))^theta
  # or exp( theta * (log(theta) - log(theta+mu)) )
  prob_zero <- (theta / (theta + mu))^theta
  
  # Safeguard P0 to avoid division by zero (if P0 ~ 1)
  prob_zero <- pmin(prob_zero, 1 - 1e-8)
  
  # Expected value of ZTNB
  mu_ztnb <- mu / (1 - prob_zero)
  
  # Residual
  resid <- Y - mu_ztnb
  resid[Y == 0] <- 0
  
  # Weight from GLM score derivation for NB
  # Score_NB = (Y-mu)/(1+mu/theta) * X
  weight <- 1 / (1 + mu / theta)
  
  # Combined residual for score
  # U = X * weight * (Y - mu_ztnb)  ?
  # Let's double check the correction term.
  # The derivation: U_ztnb = U_nb + term.
  # U_nb = (Y - mu) * weight * X
  # term = - (mu * weight * X) * (P0 / (1-P0))
  # Total = weight * X * [ (Y - mu) - mu * P0/(1-P0) ]
  #       = weight * X * [ Y - mu * (1 + P0/(1-P0)) ]
  #       = weight * X * [ Y - mu * (1/(1-P0)) ]
  #       = weight * X * (Y - mu_ztnb)
  # Correct.
  
  weighted_resid <- as.vector(weight * resid)
  
  # Safely sweep
  return(sweep(X_mat, 1, weighted_resid, "*"))
}

score_binary <- function(Y_binary, X_mat, beta) {
  probs <- plogis(X_mat %*% beta)
  probs <- pmax(pmin(probs, 1 - 1e-6), 1e-6)
  resid <- Y_binary - probs
  return(sweep(X_mat, 1, as.vector(resid), "*"))
}

# ------------------------------------------------------------------------------
# 2. Estimating Function
# ------------------------------------------------------------------------------

estfun_Hurdle_NB <- function(data, models) {
  
  X <- data$X
  Y <- data$Y
  Y_bin <- as.numeric(Y > 0)
  
  formula_e <- grab_fixed_formula(models$e)
  formula_z <- grab_fixed_formula(models$zero)
  formula_c <- grab_fixed_formula(models$count) # Dummy object
  
  # Retrieve fixed theta from models list
  theta_fixed <- models$theta_fixed
  if(is.null(theta_fixed)) stop("theta_fixed must be provided in models list")

  Xe <- gdm(data, formula_e)
  Xz <- gdm(data, formula_z)
  Xc <- gdm(data, formula_c)
  
  data0 <- data; data0$X <- 0
  data1 <- data; data1$X <- 1
  
  Xz0 <- gdm(data0, formula_z)
  Xz1 <- gdm(data1, formula_z)
  Xc0 <- gdm(data0, formula_c)
  Xc1 <- gdm(data1, formula_c)
  
  p_e <- ncol(Xe)
  p_z <- ncol(Xz)
  p_c <- ncol(Xc)
  
  idx_e <- 1:p_e
  idx_z <- (p_e + 1):(p_e + p_z)
  idx_c <- (p_e + p_z + 1):(p_e + p_z + p_c)
  idx_cm1 <- p_e + p_z + p_c + 1
  idx_cm0 <- p_e + p_z + p_c + 2
  idx_ate <- p_e + p_z + p_c + 3
  
  function(theta) {
    beta_e <- theta[idx_e]
    beta_z <- theta[idx_z]
    beta_c <- theta[idx_c]
    mu1    <- theta[idx_cm1]
    mu0    <- theta[idx_cm0]
    ate    <- theta[idx_ate]
    
    # 1. Propensity (Logistic)
    ps <- plogis(Xe %*% beta_e)
    ps <- pmax(pmin(ps, 1 - 1e-6), 1e-6)
    score_e <- sweep(Xe, 1, as.vector(X - ps), "*")
    
    # 2. Zero Hurdle (Logistic)
    score_z <- score_binary(Y_bin, Xz, beta_z)
    
    # 3. Count Part (ZTNB)
    score_c <- score_ztnb(Y, Xc, beta_c, theta_fixed)
    
    # 4. Predictions
    # Helper to get expected mean ZTNB
    mean_ztnb <- function(X_m, b_c, th) {
       mu <- exp(X_m %*% b_c)
       mu <- pmax(mu, 1e-6)
       p0 <- (th / (th + mu))^th
       p0 <- pmin(p0, 1 - 1e-8)
       return(mu / (1 - p0))
    }
    
    mu_ztnb_1 <- mean_ztnb(Xc1, beta_c, theta_fixed)
    mu_ztnb_0 <- mean_ztnb(Xc0, beta_c, theta_fixed)
    
    prob_pos_1 <- plogis(Xz1 %*% beta_z)
    prob_pos_0 <- plogis(Xz0 %*% beta_z)
    
    m1 <- prob_pos_1 * mu_ztnb_1
    m0 <- prob_pos_0 * mu_ztnb_0
    
    # 5. Influence Functions
    psi_cm1 <- (X * Y - (X - ps) * m1) / ps - mu1
    psi_cm0 <- ((1 - X) * Y + (X - ps) * m0) / (1 - ps) - mu0
    psi_ate <- mu1 - mu0 - ate
    
    scores <- cbind(score_e, score_z, score_c, psi_cm1, psi_cm0, psi_ate)
    return(unname(scores))
  }
}

# ------------------------------------------------------------------------------
# 3. Wrapper Function
# ------------------------------------------------------------------------------

geex_Hurdle_NB <- function(data, 
                           propensity_formula, 
                           zero_formula, 
                           count_formula) {
  
  # 1. Fit Propensity
  m_ps <- glm(propensity_formula, data = data, family = binomial)
  
  # 2. Fit Hurdle Zero
  data$Y_bin <- as.numeric(data$Y > 0)
  m_zero <- glm(zero_formula, data = data, family = binomial)
  
  # 3. Fit Count (Initial Estimate for theta and beta)
  # Use MASS::glm.nb on positive data to get theta and starting beta
  # Note: Standard NB is approximation of ZTNB
  data_pos <- data[data$Y > 0, ]
  
  # We use tryCatch because glm.nb might fail if data is not overdispersed
  m_count_nb <- tryCatch({
    glm.nb(count_formula, data = data_pos)
  }, error = function(e) {
    # Fallback to Poisson if NB fails (theta = Inf)
    warning("NB fit failed, falling back to Poisson (fixed large theta)")
    glm(count_formula, data = data_pos, family = poisson)
  })
  
  if(inherits(m_count_nb, "negbin")) {
    theta_est <- m_count_nb$theta
    start_beta_c <- coef(m_count_nb)
  } else {
    theta_est <- 10000 # Large value ~ Poisson
    start_beta_c <- coef(m_count_nb)
  }
  
  # Refine ZTNB means? 
  # We could do a custom MLE for ZTNB here if we strictly want good start roots,
  # but standard NB is usually very close for initialization.
  root_c <- start_beta_c
  
  # 4. Initial Means
  Xc_full <- model.matrix(count_formula, data)
  Xz_full <- model.matrix(zero_formula, data)
  
  # Helper for init
    mean_ztnb_init <- function(X_m, b_c, th) {
       mu <- exp(X_m %*% b_c) 
       p0 <- (th / (th + mu))^th
       # If theta is huge, p0 -> exp(-mu)
       return(mu / (1 - p0))
    }
  
  mu_ztnb_full <- mean_ztnb_init(Xc_full, root_c, theta_est)
  p_pos_full <- plogis(Xz_full %*% coef(m_zero))
  
  # Re-calculate m1/m0 roughly
  # (Proper way is to use counterfactual design matrices inside estfun, 
  # but here we just need a number for the root)
  dat1 <- data; dat1$X <- 1
  dat0 <- data; dat0$X <- 0
  
  Xc1 <- model.matrix(count_formula, dat1); Xz1 <- model.matrix(zero_formula, dat1)
  Xc0 <- model.matrix(count_formula, dat0); Xz0 <- model.matrix(zero_formula, dat0)
  
  m1_vec <- plogis(Xz1 %*% coef(m_zero)) * mean_ztnb_init(Xc1, root_c, theta_est)
  m0_vec <- plogis(Xz0 %*% coef(m_zero)) * mean_ztnb_init(Xc0, root_c, theta_est)
  
  m1_start <- mean(m1_vec)
  m0_start <- mean(m0_vec)
  
  roots <- c(coef(m_ps), coef(m_zero), root_c, m1_start, m0_start, m1_start - m0_start)
  roots <- unname(roots)
  
  # Models list
  # We construct a dummy glm object for count to carry formula
  # And PASS theta_fixed
  m_count_dummy <- glm(count_formula, data = data, family = poisson)
  
  models <- list(e = m_ps,
                 zero = m_zero,
                 count = m_count_dummy,
                 theta_fixed = theta_est)
  
  # GEEX
  geex_results <- m_estimate(
    estFUN = estfun_Hurdle_NB,
    data = data,
    roots = roots,
    compute_roots = FALSE,
    outer_args = list(models = models)
  )
  
  n_params <- length(geex_results@estimates)
  est <- geex_results@estimates[n_params]
  se <- sqrt(geex_results@vcov[n_params, n_params])
  
  return(list(result=data.frame(ATE=est, SE=se, Type="Hurdle-ZTNB"), 
              theta_used=theta_est))
}

# ------------------------------------------------------------------------------
# 4. Weighted (IPW) Hurdle NB Estimator
# ------------------------------------------------------------------------------

estfun_WTD_Hurdle_NB <- function(data, models) {
  
  X <- data$X # Treatment assignments
  # Note: Weights depend on Propensity Score e(X)
  
  Y <- data$Y
  Y_bin <- as.numeric(Y > 0)
  
  formula_e <- grab_fixed_formula(models$e)
  formula_z <- grab_fixed_formula(models$zero)
  formula_c <- grab_fixed_formula(models$count)
  
  theta_fixed <- models$theta_fixed
  if(is.null(theta_fixed)) stop("theta_fixed must be provided")

  Xe <- gdm(data, formula_e)
  Xz <- gdm(data, formula_z)
  Xc <- gdm(data, formula_c)
  
  # For predictions (G-computation step)
  data0 <- data; data0$X <- 0
  data1 <- data; data1$X <- 1
  
  Xz0 <- gdm(data0, formula_z)
  Xz1 <- gdm(data1, formula_z)
  Xc0 <- gdm(data0, formula_c)
  Xc1 <- gdm(data1, formula_c)
  
  p_e <- ncol(Xe)
  p_z <- ncol(Xz)
  p_c <- ncol(Xc)
  
  idx_e <- 1:p_e
  idx_z <- (p_e + 1):(p_e + p_z)
  idx_c <- (p_e + p_z + 1):(p_e + p_z + p_c)
  idx_cm1 <- p_e + p_z + p_c + 1
  idx_cm0 <- p_e + p_z + p_c + 2
  idx_ate <- p_e + p_z + p_c + 3
  
  function(theta) {
    beta_e <- theta[idx_e]
    beta_z <- theta[idx_z]
    beta_c <- theta[idx_c]
    mu1    <- theta[idx_cm1]
    mu0    <- theta[idx_cm0]
    ate    <- theta[idx_ate]
    
    # 1. Propensity Scores
    ps <- plogis(Xe %*% beta_e)
    ps <- pmax(pmin(ps, 1 - 1e-6), 1e-6)
    score_e <- sweep(Xe, 1, as.vector(X - ps), "*")
    
    # 2. Weights
    # IPW weight: A/e + (1-A)/(1-e)
    # W is used to weight the Hurdle Scores
    W <- (X / ps) + ((1 - X) / (1 - ps))
    
    # 3. Weighted Hurdle Scores
    # Zero Part
    score_z_raw <- score_binary(Y_bin, Xz, beta_z)
    score_z_wtd <- sweep(score_z_raw, 1, W, "*")
    
    # Count Part
    score_c_raw <- score_ztnb(Y, Xc, beta_c, theta_fixed)
    score_c_wtd <- sweep(score_c_raw, 1, W, "*")
    
    # 4. Predictions (G-Computation using Weighted Parameters)
    mean_ztnb <- function(X_m, b_c, th) {
       mu <- exp(X_m %*% b_c)
       mu <- pmax(mu, 1e-6)
       p0 <- (th / (th + mu))^th
       p0 <- pmin(p0, 1 - 1e-8)
       return(mu / (1 - p0))
    }
    
    mu_ztnb_1 <- mean_ztnb(Xc1, beta_c, theta_fixed)
    mu_ztnb_0 <- mean_ztnb(Xc0, beta_c, theta_fixed)
    
    prob_pos_1 <- plogis(Xz1 %*% beta_z)
    prob_pos_0 <- plogis(Xz0 %*% beta_z)
    
    m1 <- prob_pos_1 * mu_ztnb_1
    m0 <- prob_pos_0 * mu_ztnb_0
    
    # 5. Causal Means (Standardization)
    # psi = m1 - mu1
    psi_cm1 <- m1 - mu1
    psi_cm0 <- m0 - mu0
    psi_ate <- mu1 - mu0 - ate
    
    scores <- cbind(score_e, score_z_wtd, score_c_wtd, psi_cm1, psi_cm0, psi_ate)
    return(unname(scores))
  }
}

geex_WTD_Hurdle_NB <- function(data, 
                               propensity_formula, 
                               zero_formula, 
                               count_formula) {
  
  # 1. Fit Propensity
  m_ps <- glm(propensity_formula, data = data, family = binomial)
  ps_pred <- predict(m_ps, type="response")
  weights <- (data$X / ps_pred) + ((1 - data$X) / (1 - ps_pred))
  
  # 2. Fit Weighted Hurdle Zero
  data$Y_bin <- as.numeric(data$Y > 0)
  m_zero <- glm(zero_formula, data = data, family = binomial, weights = weights)
  
  # 3. Fit Weighted Count (Initial Estimate)
  data_pos <- data[data$Y > 0, ]
  # Subset weights for positive cases
  # Note: weights vector aligns with 'data'. Need to subset.
  weights_pos <- weights[data$Y > 0]
  
  m_count_nb <- tryCatch({
    glm.nb(count_formula, data = data_pos, weights = weights_pos)
  }, error = function(e) {
    warning("NB fit failed, falling back to Poisson")
    glm(count_formula, data = data_pos, family = poisson, weights = weights_pos)
  })
  
  if(inherits(m_count_nb, "negbin")) {
    theta_est <- m_count_nb$theta
    start_beta_c <- coef(m_count_nb)
  } else {
    theta_est <- 10000
    start_beta_c <- coef(m_count_nb)
  }
  
  root_c <- start_beta_c
  
  # 4. Roots for Means
  Xc_full <- model.matrix(count_formula, data)
  Xz_full <- model.matrix(zero_formula, data)
  
  mean_ztnb_init <- function(X_m, b_c, th) {
       mu <- exp(X_m %*% b_c) 
       p0 <- (th / (th + mu))^th
       return(mu / (1 - p0))
  }
  
  # G-Computation step for initial values using weighted model
  # Note: We need to predict for everyone as if treated/untreated
  dat1 <- data; dat1$X <- 1
  dat0 <- data; dat0$X <- 0
  
  Xc1 <- model.matrix(count_formula, dat1)
  Xz1 <- model.matrix(zero_formula, dat1)
  Xc0 <- model.matrix(count_formula, dat0)
  Xz0 <- model.matrix(zero_formula, dat0)
  
  m1_vec <- plogis(Xz1 %*% coef(m_zero)) * mean_ztnb_init(Xc1, root_c, theta_est)
  m0_vec <- plogis(Xz0 %*% coef(m_zero)) * mean_ztnb_init(Xc0, root_c, theta_est)
  
  m1_start <- mean(m1_vec)
  m0_start <- mean(m0_vec)
  
  roots <- c(coef(m_ps), coef(m_zero), root_c, m1_start, m0_start, m1_start - m0_start)
  roots <- unname(roots)
  
  # Models list
  m_count_dummy <- glm(count_formula, data = data, family = poisson)
  models <- list(e = m_ps,
                 zero = m_zero,
                 count = m_count_dummy,
                 theta_fixed = theta_est)
  
  # GEEX
  geex_results <- m_estimate(
    estFUN = estfun_WTD_Hurdle_NB,
    data = data,
    roots = roots,
    compute_roots = FALSE,
    outer_args = list(models = models)
  )
  
  n_params <- length(geex_results@estimates)
  est <- geex_results@estimates[n_params]
  se <- sqrt(geex_results@vcov[n_params, n_params])
  
  return(list(result=data.frame(ATE=est, SE=se, Type="WTD-Hurdle-ZTNB"), 
              theta_used=theta_est))
}

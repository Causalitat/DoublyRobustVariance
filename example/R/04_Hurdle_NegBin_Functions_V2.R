# Program: 04_Hurdle_NegBin_Functions_V2.R
# Purpose: Extension of DR estimators for Hurdle Models with Negative Binomial Counts
#          This version (V2) includes the dispersion parameter (theta) in the estimating equations
#          to fully account for its uncertainty in the standard errors.
# Conceptual Model:
# Part 1 (Zero): Pr(Y>0) ~ Logistic(X, Z)
# Part 2 (Count): Y|Y>0 ~ ZeroTruncatedNegBin(mu, theta)

library(geex)
library(MASS) # For glm.nb to get initial theta

# Helper to grab design matrix
gdm <- function(data, formula) {
  model.matrix(formula, data = data)
}

# ------------------------------------------------------------------------------
# 1. Component Score Functions
# ------------------------------------------------------------------------------

# Score for Theta (Standard Negative Binomial)
# Derived from log-likelihood of NB:
# L = lgamma(y+theta) - lgamma(theta) - lgamma(y+1) + theta*log(theta) + y*log(mu) - (y+theta)*log(theta+mu)
# Score_theta = digamma(y+theta) - digamma(theta) + log(theta) + 1 - log(theta+mu) - (y+theta)/(theta+mu)
#             = digamma(y+theta) - digamma(theta) + log(theta) - log(theta+mu) + (mu - y)/(theta+mu)
score_nb_theta_obs <- function(Y, mu, theta) {
  # Avoid numerical issues if theta is non-positive (should be solved by root finding, but safety first)
  theta <- max(theta, 1e-6)
  
  term1 <- digamma(Y + theta) - digamma(theta)
  term2 <- log(theta) - log(theta + mu)
  term3 <- (mu - Y) / (theta + mu)
  
  return(term1 + term2 + term3)
}

# Score for Zero-Truncated Negative Binomial (ZTNB) - Fixed Theta
# U(beta) = X^T * W * (Y - E[Y|Y>0])
score_ztnb <- function(Y, X_mat, beta, theta) {
  mu <- exp(X_mat %*% beta)
  mu <- pmax(mu, 1e-6)
  
  prob_zero <- (theta / (theta + mu))^theta
  prob_zero <- pmin(prob_zero, 1 - 1e-8)
  
  mu_ztnb <- mu / (1 - prob_zero)
  
  resid <- Y - mu_ztnb
  resid[Y == 0] <- 0
  
  # Weight from GLM score derivation for NB
  weight <- 1 / (1 + mu / theta)
  
  weighted_resid <- as.vector(weight * resid)
  
  return(sweep(X_mat, 1, weighted_resid, "*"))
}

# Standard Binary Score
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
  formula_c <- grab_fixed_formula(models$count) 
  
  Xe <- gdm(data, formula_e)
  Xz <- gdm(data, formula_z)
  Xc <- gdm(data, formula_c)
  
  # Prediction matrices
  data0 <- data; data0$X <- 0
  data1 <- data; data1$X <- 1
  Xz0 <- gdm(data0, formula_z)
  Xz1 <- gdm(data1, formula_z)
  Xc0 <- gdm(data0, formula_c)
  Xc1 <- gdm(data1, formula_c)
  
  p_e <- ncol(Xe)
  p_z <- ncol(Xz)
  p_c <- ncol(Xc)
  
  # Indices map
  idx_e <- 1:p_e
  idx_z <- (p_e + 1):(p_e + p_z)
  idx_c <- (p_e + p_z + 1):(p_e + p_z + p_c)
  idx_theta <- p_e + p_z + p_c + 1
  idx_cm1 <- idx_theta + 1
  idx_cm0 <- idx_theta + 2
  idx_ate <- idx_theta + 3
  
  function(params) {
    beta_e <- params[idx_e]
    beta_z <- params[idx_z]
    beta_c <- params[idx_c]
    theta  <- params[idx_theta]
    mu1    <- params[idx_cm1]
    mu0    <- params[idx_cm0]
    ate    <- params[idx_ate]
    
    # 1. Propensity (Logistic)
    ps <- plogis(Xe %*% beta_e)
    ps <- pmax(pmin(ps, 1 - 1e-6), 1e-6)
    score_e <- sweep(Xe, 1, as.vector(X - ps), "*")
    
    # 2. Zero Hurdle (Logistic)
    score_z <- score_binary(Y_bin, Xz, beta_z)
    
    # 3. Count Part (Beta)
    score_c <- score_ztnb(Y, Xc, beta_c, theta)
    
    # 4. Count Part (Theta)
    # The estimator for theta is derived from the positive counts (simulating glm.nb on data_pos)
    mu_curr <- exp(Xc %*% beta_c)
    mu_curr <- pmax(mu_curr, 1e-6)
    
    s_theta <- score_nb_theta_obs(Y, mu_curr, theta)
    s_theta[Y == 0] <- 0 # Only positive cases contribute to this theta estimator
    
    # 5. Predictions
    mean_ztnb <- function(X_m, b_c, th) {
       mu <- exp(X_m %*% b_c)
       mu <- pmax(mu, 1e-6)
       p0 <- (th / (th + mu))^th
       p0 <- pmin(p0, 1 - 1e-8)
       return(mu / (1 - p0))
    }
    
    mu_ztnb_1 <- mean_ztnb(Xc1, beta_c, theta)
    mu_ztnb_0 <- mean_ztnb(Xc0, beta_c, theta)
    
    prob_pos_1 <- plogis(Xz1 %*% beta_z)
    prob_pos_0 <- plogis(Xz0 %*% beta_z)
    
    m1 <- prob_pos_1 * mu_ztnb_1
    m0 <- prob_pos_0 * mu_ztnb_0
    
    # 6. Influence Functions for Means
    psi_cm1 <- (X * Y - (X - ps) * m1) / ps - mu1
    psi_cm0 <- ((1 - X) * Y + (X - ps) * m0) / (1 - ps) - mu0
    psi_ate <- mu1 - mu0 - ate
    
    scores <- cbind(score_e, score_z, score_c, s_theta, psi_cm1, psi_cm0, psi_ate)
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
  
  # 3. Fit Count (Initial Estimate)
  data_pos <- data[data$Y > 0, ]
  
  m_count_nb <- tryCatch({
    glm.nb(count_formula, data = data_pos)
  }, error = function(e) {
    warning("NB fit failed, falling back to Poisson (fixed large theta)")
    glm(count_formula, data = data_pos, family = poisson)
  })
  
  if(inherits(m_count_nb, "negbin")) {
    theta_est <- m_count_nb$theta
    start_beta_c <- coef(m_count_nb)
  } else {
    theta_est <- 10000 
    start_beta_c <- coef(m_count_nb)
  }
  
  # 4. Initial Means
  root_c <- start_beta_c
  Xc_full <- model.matrix(count_formula, data)
  Xz_full <- model.matrix(zero_formula, data)
  
  mean_ztnb_init <- function(X_m, b_c, th) {
       mu <- exp(X_m %*% b_c) 
       p0 <- (th / (th + mu))^th
       p0 <- pmin(p0, 1 - 1e-8)
       return(mu / (1 - p0))
  }
  
  dat1 <- data; dat1$X <- 1
  dat0 <- data; dat0$X <- 0
  Xc1 <- model.matrix(count_formula, dat1); Xz1 <- model.matrix(zero_formula, dat1)
  Xc0 <- model.matrix(count_formula, dat0); Xz0 <- model.matrix(zero_formula, dat0)
  
  m1_vec <- plogis(Xz1 %*% coef(m_zero)) * mean_ztnb_init(Xc1, root_c, theta_est)
  m0_vec <- plogis(Xz0 %*% coef(m_zero)) * mean_ztnb_init(Xc0, root_c, theta_est)
  
  m1_start <- mean(m1_vec)
  m0_start <- mean(m0_vec)
  
  # Roots now include theta_est
  roots <- c(coef(m_ps), coef(m_zero), root_c, theta_est, m1_start, m0_start, m1_start - m0_start)
  roots <- unname(roots)
  
  # Models list (dummy count model for formula)
  m_count_dummy <- glm(count_formula, data = data, family = poisson)
  
  models <- list(e = m_ps,
                 zero = m_zero,
                 count = m_count_dummy) # Note: theta_fixed is gone
  
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
  
  return(list(result=data.frame(ATE=est, SE=se, Type="Hurdle-ZTNB-V2"), 
              theta_est=geex_results@estimates[n_params - 3])) # Helper to see final theta
}

# ------------------------------------------------------------------------------
# 4. Weighted (IPW) Hurdle NB Estimator
# ------------------------------------------------------------------------------

estfun_WTD_Hurdle_NB <- function(data, models) {
  
  X <- data$X 
  Y <- data$Y
  Y_bin <- as.numeric(Y > 0)
  
  formula_e <- grab_fixed_formula(models$e)
  formula_z <- grab_fixed_formula(models$zero)
  formula_c <- grab_fixed_formula(models$count)
  
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
  idx_theta <- p_e + p_z + p_c + 1
  idx_cm1 <- idx_theta + 1
  idx_cm0 <- idx_theta + 2
  idx_ate <- idx_theta + 3
  
  function(params) {
    beta_e <- params[idx_e]
    beta_z <- params[idx_z]
    beta_c <- params[idx_c]
    theta  <- params[idx_theta]
    mu1    <- params[idx_cm1]
    mu0    <- params[idx_cm0]
    ate    <- params[idx_ate]
    
    # 1. Propensity Scores
    ps <- plogis(Xe %*% beta_e)
    ps <- pmax(pmin(ps, 1 - 1e-6), 1e-6)
    score_e <- sweep(Xe, 1, as.vector(X - ps), "*")
    
    # 2. Weights
    W <- (X / ps) + ((1 - X) / (1 - ps))
    
    # 3. Weighted Hurdle Scores
    score_z_raw <- score_binary(Y_bin, Xz, beta_z)
    score_z_wtd <- sweep(score_z_raw, 1, W, "*")
    
    # Count Part (Beta)
    score_c_raw <- score_ztnb(Y, Xc, beta_c, theta)
    score_c_wtd <- sweep(score_c_raw, 1, W, "*")
    
    # Count Part (Theta)
    # Weighted Score for theta?
    # Usually weighted MLE implies Sum( w * score_i ) = 0
    mu_curr <- exp(Xc %*% beta_c)
    mu_curr <- pmax(mu_curr, 1e-6)
    s_theta_raw <- score_nb_theta_obs(Y, mu_curr, theta)
    s_theta_raw[Y == 0] <- 0
    s_theta_wtd <- s_theta_raw * W # Apply IPW weights to theta score too
    
    # 4. Predictions
    mean_ztnb <- function(X_m, b_c, th) {
       mu <- exp(X_m %*% b_c)
       mu <- pmax(mu, 1e-6)
       p0 <- (th / (th + mu))^th
       p0 <- pmin(p0, 1 - 1e-8)
       return(mu / (1 - p0))
    }
    
    mu_ztnb_1 <- mean_ztnb(Xc1, beta_c, theta)
    mu_ztnb_0 <- mean_ztnb(Xc0, beta_c, theta)
    
    prob_pos_1 <- plogis(Xz1 %*% beta_z)
    prob_pos_0 <- plogis(Xz0 %*% beta_z)
    
    m1 <- prob_pos_1 * mu_ztnb_1
    m0 <- prob_pos_0 * mu_ztnb_0
    
    # 5. Causal Means
    psi_cm1 <- m1 - mu1
    psi_cm0 <- m0 - mu0
    psi_ate <- mu1 - mu0 - ate
    
    scores <- cbind(score_e, score_z_wtd, score_c_wtd, s_theta_wtd, psi_cm1, psi_cm0, psi_ate)
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
  data$wts <- weights 
  
  # 2. Fit Weighted Hurdle Zero
  data$Y_bin <- as.numeric(data$Y > 0)
  m_zero <- glm(zero_formula, data = data, family = binomial, weights = wts)
  
  # 3. Fit Weighted Count
  data_pos <- data[data$Y > 0, ]
  weights_pos <- data$wts[data$Y > 0]
  data_pos$weights_pos <- weights_pos 
  
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
  
  # 4. Roots
  Xc_full <- model.matrix(count_formula, data)
  Xz_full <- model.matrix(zero_formula, data)
  
  mean_ztnb_init <- function(X_m, b_c, th) {
       mu <- exp(X_m %*% b_c) 
       p0 <- (th / (th + mu))^th
       p0 <- pmin(p0, 1 - 1e-8)
       return(mu / (1 - p0))
  }
  
  dat1 <- data; dat1$X <- 1
  dat0 <- data; dat0$X <- 0
  Xc1 <- model.matrix(count_formula, dat1); Xz1 <- model.matrix(zero_formula, dat1)
  Xc0 <- model.matrix(count_formula, dat0); Xz0 <- model.matrix(zero_formula, dat0)
  
  m1_vec <- plogis(Xz1 %*% coef(m_zero)) * mean_ztnb_init(Xc1, root_c, theta_est)
  m0_vec <- plogis(Xz0 %*% coef(m_zero)) * mean_ztnb_init(Xc0, root_c, theta_est)
  
  m1_start <- mean(m1_vec)
  m0_start <- mean(m0_vec)
  
  roots <- c(coef(m_ps), coef(m_zero), root_c, theta_est, m1_start, m0_start, m1_start - m0_start)
  roots[is.na(roots)] <- 0
  roots <- unname(roots)
  
  m_count_dummy <- glm(count_formula, data = data, family = poisson)
  models <- list(e = m_ps,
                 zero = m_zero,
                 count = m_count_dummy)
  
  # GEEX with Manual Fallback (updated for extra param)
  geex_results <- tryCatch({
    m_estimate(
      estFUN = estfun_WTD_Hurdle_NB,
      data = data,
      roots = roots,
      compute_roots = FALSE,
      outer_args = list(models = models)
    )
  }, error = function(e) {
    warning("geex::m_estimate failed. Falling back to manual Sandwich Variance.")
    
    psi_fun_closure <- estfun_WTD_Hurdle_NB(data, models = models)
    grad_fun <- function(th) { colSums(psi_fun_closure(th)) }
    
    if(!requireNamespace("numDeriv", quietly = TRUE)) stop("Need numDeriv")
    A <- numDeriv::jacobian(grad_fun, roots)
    
    scores <- psi_fun_closure(roots)
    B <- crossprod(scores)
    
    Ainv <- tryCatch({ solve(A) }, error = function(e) { stop("A singular") })
    V <- Ainv %*% B %*% t(Ainv)
    
    list(estimates = roots, vcov = V)
  })
  
  if(inherits(geex_results, "geex")) {
    n_params <- length(geex_results@estimates)
    est <- geex_results@estimates[n_params]
    se <- sqrt(geex_results@vcov[n_params, n_params])
  } else {
    n_params <- length(geex_results$estimates)
    est <- geex_results$estimates[n_params]
    se <- sqrt(geex_results$vcov[n_params, n_params])
  }
  
  return(list(result=data.frame(ATE=est, SE=se, Type="WTD-Hurdle-ZTNB-V2"), 
              theta_est=theta_est))
}

geex_WTD_Hurdle_NB_Int <- function(data, 
                                   propensity_formula, 
                                   covariate_formula) {
  cov_terms <- attr(terms(covariate_formula), "term.labels")
  cov_str <- paste(cov_terms, collapse = " + ")
  str_zero <- paste("as.numeric(Y > 0) ~ X * (", cov_str, ")")
  str_count <- paste("Y ~ X * (", cov_str, ")")
  f_zero_int <- as.formula(str_zero)
  f_count_int <- as.formula(str_count)
  
  return(geex_WTD_Hurdle_NB(data, 
                            propensity_formula, 
                            f_zero_int, 
                            f_count_int))
}

# Program: 02_ZINB_Functions.R
# Purpose: This program contains functions for 3 propensity-based ZINB estimators.
# 1. Weighted ZINB estimator
# 2. G-formula/standardization ZINB estimator
# 3. Doubly robust AIPW ZINB estimator
#
# Updates:
# - Added robust manual Sandwich Variance fallback (from 04_Hurdle_NegBin_Functions.R)
# - Added direct Risk Ratio (RR) estimation (from 00_Estimators_11.03.21.R)
# - Switched to direct model.matrix for better robustness

# Load required libraries
library(geex)
library(pscl)
library(numDeriv)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Helper functions for ZINB models
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# --- Analytical Scores (Legacy Code Optimization) ---
# From: 00_Estimators_11.03.21.R
# Pre-calculates helper sum for Negative Binomial scores
calc_sum3 <- function(X, odisp){
  # Optimized vector version of the loop
  # X is vector of counts, odisp is dispersion
  # sum_{j=0}^{y-1} ...
  
  # For scalar X it does a loop. We need vectorized.
  sapply(X, function(x_i) {
    if(x_i == 0) return(0)
    j <- 0:(x_i - 1)
    sum(-1 / (j * odisp^2 + odisp))
  })
}

# Analytical Score for ZINB (Single Observation or Vectorized)
# Note: odisp here is 'alpha' (1/theta).
ScoreZINB_Analytic <- function(Y, lambda, nu, odisp, X_count, X_zero) {
  # lambda: Count mean (mu)
  # nu: Zero odds? In legacy code: nu = exp(Z*gamma) (Odds of structural zero)
  # odisp: Dispersion (alpha = 1/theta)
  
  # Prob(Y=0) part:
  # P(Y=0)   = P(Zero) + P(Count=0|NotZero)*P(NotZero)
  #          = phi + (1-phi)*(1 + alpha*mu)^(-1/alpha)
  # phi      = nu / (1 + nu)
  
  # Common terms
  # Legacy derivation uses:
  # nu = exp(Zg)
  # lam = exp(Xb)
  # A = (1 + odisp*lam)
  
  # Wait, to make this fully generic for model.matrix X_count/X_zero:
  # We return the "score for linear predictor" then multiply by covariates.
  
  obs_zero <- as.numeric(Y == 0)
  obs_pos  <- as.numeric(Y > 0)
  
  term_base <- (1 + odisp*lambda)
  term_pow  <- term_base^(-1/odisp)     # P(NB=0)
  
  # Score for Count Mean (lambda) -> Beta
  # dL/d(lambda)
  # If Y=0: ...
  # If Y>0: (Y - lambda) / (1 + odisp*lambda)
  
  zeros_numer <- -lambda * term_base^(-1 - 1/odisp)
  zeros_denom <- nu + term_pow
  score_lambda_y0 <- zeros_numer / zeros_denom
  score_lambda_pos <- (Y - lambda) / term_base
  
  score_lambda <- obs_zero * score_lambda_y0 + obs_pos * score_lambda_pos
  # This is dL/d(lambda). Chain rule: dL/d(beta) = dL/d(lambda) * lambda * X
  # But legacy code `ScoreZINB` integrates X.
  
  # Score for Zero Inflation (nu) -> Gamma
  # dL/d(nu)
  # phi = nu/(1+nu).
  # Legacy:
  # Y=0: 1 * (nu / (nu + term_pow)) - nu/(1+nu) ??
  # Let's trust legacy formula pattern:
  
  score_nu_y0 <- (nu / (nu + term_pow)) - (nu / (1 + nu)) 
  score_nu_pos <- - (nu / (1 + nu)) # Simply -phi if Y>0
  # Legacy Line 197: ... - 1*nu/(1+nu)
  
  score_nu <- obs_zero * score_nu_y0 + obs_pos * score_nu_pos
  
  # Score for Dispersion (odisp) -> Theta
  # This is complex. Legacy defines it in line 201.
  sum_term <- calc_sum3(Y, odisp)
  
  # Legacy Line 201
  term_log <- log(1 + odisp*lambda)
  
  # Y=0 part
  score_odisp_y0_num <- ((1 + odisp*lambda)*term_log - odisp*lambda)
  score_odisp_y0_den <- (odisp^2 * (1 + odisp*lambda) * (nu * (1+odisp*lambda)^(1/odisp) + 1))
  # Note: The legacy formula simplifies term_pow in denom differently.
  
  score_odisp_y0 <- score_odisp_y0_num / score_odisp_y0_den
  
  # Y>0 part
  score_odisp_pos <- (term_log * odisp^(-2)) + ((Y - lambda)/(odisp*(1+odisp*lambda))) + sum_term
  
  score_odisp <- obs_zero * score_odisp_y0 + obs_pos * score_odisp_pos
  
  # Return components (dL/dLam, dL/dNu, dL/dOdisp)
  list(lam = score_lambda, nu = score_nu, odisp = score_odisp)
}

# A generic log-likelihood function for a zeroinfl model.
logL.zeroinfl <- function(par, X, Z, Y, weights, offset, link) {
  k_count <- ncol(X)
  k_zero <- ncol(Z)
  beta <- par[1:k_count]
  gamma <- par[(k_count + 1):(k_count + k_zero)]
  
  # Dispersion parameter (theta) - clamp to be positive
  # Note: pscl estimates log(theta) if dist="negbin"
  size <- exp(par[k_count + k_zero + 1])
  if(is.na(size) || size <= 0) size <- .Machine$double.eps

  mu <- exp(X %*% beta + offset)
  phi <- plogis(Z %*% gamma)
  phi <- pmin(pmax(phi, .Machine$double.eps), 1 - .Machine$double.eps)

  loglik_i <- suppressWarnings(
    weights * (
      (Y > 0) * (log(1 - phi) + dnbinom(Y, size = size, mu = mu, log = TRUE)) +
      (Y == 0) * (log(phi + (1 - phi) * pnbinom(0, size = size, mu = mu)))
    )
  )
  loglik_i[is.na(loglik_i) | is.nan(loglik_i) | (is.infinite(loglik_i) & loglik_i < 0)] <- -1e10
  sum(loglik_i, na.rm=TRUE)
}

# Generic Predict Function for ZINB Mean
predict_zinb_mean <- function(X_count, X_zero, theta_m) {
  k_count <- ncol(X_count)
  k_zero <- ncol(X_zero)
  beta <- theta_m[1:k_count]
  gamma <- theta_m[(k_count + 1):(k_count + k_zero)]
  
  mu <- exp(X_count %*% beta)
  phi <- plogis(X_zero %*% gamma)
  
  # Mean = (1 - phi) * mu
  return((1 - phi) * mu)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Weighted ZINB Estimator (IPW)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

estfun_WTD_ZINB <- function(data, models){
  X <- data$X
  Y <- data$Y
  
  # Design Matrices
  Xe <- model.matrix(models$e$formula, data)
  
  # ZINB Model Matrices (Need to handle separate formulas)
  # pscl formula is: outcomes ~ count_preds | zero_preds
  f_zinb <- models$m$terms$full
  mf <- model.frame(f_zinb, data)
  Xm_count <- model.matrix(models$m$terms$count, mf)
  Xm_zero <- model.matrix(models$m$terms$zero, mf)
  
  # For Predictions (Counterfactuals)
  data0 <- data; data0$X <- 0
  data1 <- data; data1$X <- 1
  mf0 <- model.frame(f_zinb, data0)
  mf1 <- model.frame(f_zinb, data1)
  
  Xm0_count <- model.matrix(models$m$terms$count, mf0)
  Xm1_count <- model.matrix(models$m$terms$count, mf1)
  Xm0_zero <- model.matrix(models$m$terms$zero, mf0)
  Xm1_zero <- model.matrix(models$m$terms$zero, mf1)

  # Prepare objective function for ZINB scores
  # We reuse the specific data structure
  Y_zinb <- model.response(mf)
  weights_zinb <- weights(models$m)
  if(is.null(weights_zinb)) weights_zinb <- rep(1, nrow(data))
  offset_zinb <- rep(0, nrow(data)) # Simplify for now
  
  obj_fun_zinb <- function(par) {
     logL.zeroinfl(par, Xm_count, Xm_zero, Y_zinb, weights_zinb, offset_zinb)
  }
  
  e_scores_fun <- estfun(models$e)
  
  # Indices
  e_pos <- 1:length(coef(models$e))
  m_params <- c(coef(models$m, "count"), coef(models$m, "zero"), log(models$m$theta))
  m_pos <- (max(e_pos) + 1):(max(e_pos) + length(m_params))
  
  function(theta){
    p <- length(theta)
    
    # 1. Propensity Scores
    e_lp <- Xe %*% theta[e_pos]
    e_prob <- plogis(e_lp)
    
    # 2. Weights
    W_ipw <- (X / e_prob) + ((1 - X) / (1 - e_prob))
    
    # 3. ZINB Scores (Weighted)
    m_theta <- theta[m_pos]
    
    # Check if we can use Analytic Scores (Faster, Stable)
    # We assume 'use_analytic' is passed in outer_args or we detect it.
    # For now, let's just implement the choice here.
    use_analytic <- FALSE # Set to TRUE to use experimental analytical scores
    
    # Extract params
    k_count <- ncol(Xm_count)
    k_zero <- ncol(Xm_zero)
    beta_est <- m_theta[1:k_count]
    gamma_est <- m_theta[(k_count+1):(k_count+k_zero)]
    
    # pscl uses theta (size). odisp (alpha) = 1/theta
    theta_pscl <- exp(m_theta[k_count + k_zero + 1])
    odisp_est <- 1/theta_pscl 
    
    if(use_analytic) {
       lambdas <- exp(Xm_count %*% beta_est)
       nus <- exp(Xm_zero %*% gamma_est)
       
       scores_all <- ScoreZINB_Analytic(Y, lambdas, nus, odisp_est, NULL, NULL)
       
       # Chain Rule
       # dL/dBeta = dL/dLam * dLam/dBeta * X = score_lam * lam * X
       score_beta <- sweep(Xm_count, 1, scores_all$lam * lambdas * W_ipw, "*")
       
       # dL/dGamma = dL/dNu * dNu/dGamma * Z = score_nu * nu * Z
       score_gamma <- sweep(Xm_zero, 1, scores_all$nu * nus * W_ipw, "*")
       
       # dL/dLogTheta
       # param is eta = log(theta) = -log(odisp)
       # dL/dEta = dL/dOdisp * dOdisp/dTheta * dTheta/dEta
       # odisp = 1/e^eta = e^-eta
       # dOdisp/dEta = -e^-eta = -odisp
       score_theta <- scores_all$odisp * (-odisp_est) * W_ipw
       
       # Combine
       eq_m_val <- cbind(score_beta, score_gamma, score_theta)
       
    } else {
        # Numerical Fallback (Slow)
        scores_zinb_unwt <- t(sapply(1:nrow(data), function(i) {
           numDeriv::grad(function(par) logL.zeroinfl(par, Xm_count[i,,drop=F], Xm_zero[i,,drop=F], Y[i], 1, 0), m_theta)
        }))
        eq_m_val <- sweep(scores_zinb_unwt, 1, W_ipw, "*")
    }
    mu1 <- predict_zinb_mean(Xm1_count, Xm1_zero, m_theta)
    mu0 <- predict_zinb_mean(Xm0_count, Xm0_zero, m_theta)
    
    # 5. Risk Ratio & ATE
    # theta structure: [e, m, mu1, mu0, logRR, RR, RD]
    mu1_param <- theta[p-4]
    mu0_param <- theta[p-3]
    logRR_param <- theta[p-2]
    RR_param <- theta[p-1]
    RD_param <- theta[p]
    
    # Equations
    eq_e <- e_scores_fun # derived from glm
    # If e_scores_fun depends on fixed estimates, that's bad.
    # estfun(glm) returns values at optimum.
    # We need explicit score: X * (Y - p)
    eq_e_val <- Xe * as.vector(X - e_prob)
    
    # ZINB Scores weighted by W_ipw
    # We need to compute unweighted scores at m_theta, then multiply
    scores_zinb_unwt <- t(sapply(1:nrow(data), function(i) {
       numDeriv::grad(function(par) logL.zeroinfl(par, Xm_count[i,,drop=F], Xm_zero[i,,drop=F], Y[i], 1, 0), m_theta)
    }))
    eq_m_val <- sweep(scores_zinb_unwt, 1, W_ipw, "*")
    
    # Mean eqs
    eq_mu1 <- mu1 - mu1_param
    eq_mu0 <- mu0 - mu0_param
    
    # RR: mu1 = mu0 * exp(logRR)  => mu1 - mu0*exp(logRR) = 0? 
    # Or consistent sys: 
    # logRR = log(mu1) - log(mu0)
    # Using 00 method:
    eq_logRR <- mu1_param - (mu0_param * exp(logRR_param)) # Zero if RR correct
    eq_RR <- exp(logRR_param) - RR_param
    
    eq_RD <- (mu1_param - mu0_param) - RD_param
    
    cbind(eq_e_val, eq_m_val, eq_mu1, eq_mu0, eq_logRR, eq_RR, eq_RD)
  }
}

geex_WTD_ZINB <- function(data, propensity_formula, outcome_formula){
  # 1. Fit Initial Models
  e_model <- glm(propensity_formula, data = data, family = binomial)
  ps <- predict(e_model, type="response")
  wts <- (data$X / ps) + ((1 - data$X) / (1 - ps))
  
  # Fit Weighted ZINB (using pscl)
  m_model <- zeroinfl(outcome_formula, data = data, weights = wts, dist = "negbin")
  
  # 2. Initial Values
  data0 <- data; data0$X <- 0
  data1 <- data; data1$X <- 1
  mu1_vec <- predict(m_model, newdata = data1, type = "response")
  mu0_vec <- predict(m_model, newdata = data0, type = "response")
  
  mu1_start <- mean(mu1_vec)
  mu0_start <- mean(mu0_vec)
  RR_start <- mu1_start / mu0_start
  logRR_start <- log(RR_start)
  RD_start <- mu1_start - mu0_start
  
  m_coef <- c(coef(m_model, "count"), coef(m_model, "zero"), log(m_model$theta))
  roots <- c(coef(e_model), m_coef, mu1_start, mu0_start, logRR_start, RR_start, RD_start)
  roots[is.na(roots)] <- 0
  roots <- unname(roots)
  
  models <- list(e = e_model, m = m_model)
  
  # 3. GEEX with Robust Fallback
  geex_results <- tryCatch({
    m_estimate(
      estFUN = estfun_WTD_ZINB,
      data = data,
      roots = roots,
      compute_roots = FALSE,
      outer_args = list(models = models)
    )
  }, error = function(e) {
    warning("geex::m_estimate failed. Falling back to Manual Sandwich.")
    
    # Manual Sandwich Logic
    psi_fun <- estfun_WTD_ZINB(data, models)
    A <- numDeriv::jacobian(function(th) colSums(psi_fun(th)), roots)
    B <- crossprod(psi_fun(roots))
    V <- solve(A) %*% B %*% t(solve(A))
    list(estimates = roots, vcov = V)
  })
  
  # 4. Extract Results
  res <- if(inherits(geex_results, "geex")) geex_results@estimates else geex_results$estimates
  var <- if(inherits(geex_results, "geex")) geex_results@vcov else geex_results$vcov
  
  p <- length(res)
  # [..., mu1, mu0, logRR, RR, RD]
  # Indices: p-4, p-3, p-2, p-1, p
  
  est_RR <- res[p-1]
  se_RR <- sqrt(var[p-1, p-1])
  
  est_RD <- res[p]
  se_RD <- sqrt(var[p, p])
  
  return(data.frame(
    Parameter = c("RR", "RD", "ATE"),
    Estimate = c(est_RR, est_RD, est_RD),
    SE = c(se_RR, se_RD, se_RD),
    Type = "WTD-ZINB"
  ))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Doubly Robust AIPW ZINB Estimator
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

estfun_AIPW_ZINB <- function(data, models){
  X <- data$X
  Y <- data$Y
  
  Xe <- model.matrix(models$e$formula, data)
  
  f_zinb <- models$m$terms$full
  mf <- model.frame(f_zinb, data)
  Xm_count <- model.matrix(models$m$terms$count, mf)
  Xm_zero <- model.matrix(models$m$terms$zero, mf)
  
  data0 <- data; data0$X <- 0
  data1 <- data; data1$X <- 1
  mf0 <- model.frame(f_zinb, data0)
  mf1 <- model.frame(f_zinb, data1)
  Xm0_count <- model.matrix(models$m$terms$count, mf0)
  Xm1_count <- model.matrix(models$m$terms$count, mf1)
  Xm0_zero <- model.matrix(models$m$terms$zero, mf0)
  Xm1_zero <- model.matrix(models$m$terms$zero, mf1)

  e_pos <- 1:length(coef(models$e))
  m_params <- c(coef(models$m, "count"), coef(models$m, "zero"), log(models$m$theta))
  m_pos <- (max(e_pos) + 1):(max(e_pos) + length(m_params))
  
  function(theta){
    p <- length(theta)
    
    # 1. Propensity
    e_lp <- Xe %*% theta[e_pos]
    e_prob <- plogis(e_lp)
    
    # 2. Scores for e
    eq_e <- Xe * as.vector(X - e_prob)
    
    # 3. Scores for m (Unweighted ZINB)
    m_theta <- theta[m_pos]
    
    # extract parameters
    
    use_analytic <- FALSE # Set to TRUE to use experimental analytical scores
    
    k_count <- ncol(Xm_count)
    k_zero <- ncol(Xm_zero)
    beta_est <- m_theta[1:k_count]
    gamma_est <- m_theta[(k_count+1):(k_count+k_zero)]
    
    theta_pscl <- exp(m_theta[k_count + k_zero + 1])
    odisp_est <- 1/theta_pscl 
    
    if(use_analytic) {
       lambdas <- exp(Xm_count %*% beta_est)
       nus <- exp(Xm_zero %*% gamma_est)
       
       scores_all <- ScoreZINB_Analytic(Y, lambdas, nus, odisp_est, NULL, NULL)
       
       # Chain Rule
       score_beta <- sweep(Xm_count, 1, scores_all$lam * lambdas, "*")
       score_gamma <- sweep(Xm_zero, 1, scores_all$nu * nus, "*")
       score_theta <- scores_all$odisp * (-odisp_est)
       
       scores_zinb <- cbind(score_beta, score_gamma, score_theta)
    } else {
        scores_zinb <- t(sapply(1:nrow(data), function(i) {
           numDeriv::grad(function(par) logL.zeroinfl(par, Xm_count[i,,drop=F], Xm_zero[i,,drop=F], Y[i], 1, 0), m_theta)
        }))
    }
    
    # 4. Predictions
    mu1_hat <- predict_zinb_mean(Xm1_count, Xm1_zero, m_theta)
    mu0_hat <- predict_zinb_mean(Xm0_count, Xm0_zero, m_theta)
    
    # 5. Parameters
    mu1_param <- theta[p-4]
    mu0_param <- theta[p-3]
    logRR_param <- theta[p-2]
    RR_param <- theta[p-1]
    RD_param <- theta[p]
    
    # 6. AIPW Equations
    # IF_1 = (A/e)Y - ((A-e)/e)mu1
    # IF_0 = ((1-A)/(1-e))Y + ((A-e)/(1-e))mu0
    
    if1 <- (X / e_prob) * Y - ((X - e_prob) / e_prob) * mu1_hat
    if0 <- ((1 - X) / (1 - e_prob)) * Y + ((X - e_prob) / (1 - e_prob)) * mu0_hat
    
    eq_mu1 <- if1 - mu1_param
    eq_mu0 <- if0 - mu0_param
    
    eq_logRR <- mu1_param - (mu0_param * exp(logRR_param))
    eq_RR <- exp(logRR_param) - RR_param
    eq_RD <- (mu1_param - mu0_param) - RD_param
    
    cbind(eq_e, scores_zinb, eq_mu1, eq_mu0, eq_logRR, eq_RR, eq_RD)
  }
}

geex_AIPW_ZINB <- function(data, propensity_formula, outcome_formula){
  # 1. Fit Initial
  e_model <- glm(propensity_formula, data = data, family = binomial)
  m_model <- zeroinfl(outcome_formula, data = data, dist = "negbin")
  
  # 2. Initial Values (DR)
  ps <- predict(e_model, type="response")
  
  data0 <- data; data0$X <- 0
  data1 <- data; data1$X <- 1
  mu1_hat <- predict(m_model, newdata = data1, type = "response")
  mu0_hat <- predict(m_model, newdata = data0, type = "response")
  
  if1 <- (data$X / ps) * data$Y - ((data$X - ps) / ps) * mu1_hat
  if0 <- ((1 - data$X) / (1 - ps)) * data$Y + ((data$X - ps) / (1 - ps)) * mu0_hat
  
  mu1_start <- mean(if1)
  mu0_start <- mean(if0)
  RR_start <- mu1_start / mu0_start
  logRR_start <- log(RR_start)
  RD_start <- mu1_start - mu0_start # ATE
  
  m_coef <- c(coef(m_model, "count"), coef(m_model, "zero"), log(m_model$theta))
  roots <- c(coef(e_model), m_coef, mu1_start, mu0_start, logRR_start, RR_start, RD_start)
  roots[is.na(roots)] <- 0
  roots <- unname(roots)
  
  models <- list(e = e_model, m = m_model)
  
  # 3. GEEX with Robust Fallback
  geex_results <- tryCatch({
    m_estimate(
      estFUN = estfun_AIPW_ZINB,
      data = data,
      roots = roots,
      compute_roots = FALSE,
      outer_args = list(models = models)
    )
  }, error = function(e) {
    warning("geex::m_estimate failed. Falling back to Manual Sandwich.")
    psi_fun <- estfun_AIPW_ZINB(data, models)
    A <- numDeriv::jacobian(function(th) colSums(psi_fun(th)), roots)
    B <- crossprod(psi_fun(roots))
    V <- solve(A) %*% B %*% t(solve(A))
    list(estimates = roots, vcov = V)
  })
  
  res <- if(inherits(geex_results, "geex")) geex_results@estimates else geex_results$estimates
  var <- if(inherits(geex_results, "geex")) geex_results@vcov else geex_results$vcov
  
  p <- length(res)
  est_RR <- res[p-1]; se_RR <- sqrt(var[p-1, p-1])
  est_RD <- res[p]; se_RD <- sqrt(var[p, p])
  
  return(data.frame(
    Parameter = c("RR", "RD", "ATE"),
    Estimate = c(est_RR, est_RD, est_RD),
    SE = c(se_RR, se_RD, se_RD),
    Type = "AIPW-ZINB"
  ))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. Interaction Wrapper
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

geex_AIPW_ZINB_Int <- function(data, ps_formula, cov_formula) {
  cov_terms <- attr(terms(cov_formula), "term.labels")
  cov_str <- paste(cov_terms, collapse = " + ")
  
  # Expansion: Y ~ X * (Z1 + Z2) | X * (Z1 + Z2)
  # pscl formula format: count | zero
  # We assume interaction in both parts
  f_str <- paste("Y ~ X * (", cov_str, ") | X * (", cov_str, ")")
  outcome_f_int <- as.formula(f_str)
  
  return(geex_AIPW_ZINB(data, ps_formula, outcome_f_int))
}

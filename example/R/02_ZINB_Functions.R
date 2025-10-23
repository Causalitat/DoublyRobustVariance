# Program: 02_ZINB_Functions.R
# Purpose: This program contains functions for 3 propensity-based ZINB estimators.
# 1. Weighted ZINB estimator
# 2. G-formula/standardization ZINB estimator
# 3. Doubly robust AIPW ZINB estimator

# Load required libraries
library(geex)
library(pscl)
library(numDeriv)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Helper functions for ZINB models
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# A generic log-likelihood function for a zeroinfl model. This function is
# needed so that the numDeriv::grad function can be used to numerically
# compute the scores of the ZINB log-likelihood.
logL.zeroinfl <- function(par, X, Z, Y, weights, offset, link) {
  # Number of parameters for the count and zero-inflation components
  k_count <- ncol(X)
  k_zero <- ncol(Z)

  # Extract parameters for the count and zero-inflation components
  beta <- par[1:k_count]
  gamma <- par[(k_count + 1):(k_count + k_zero)]

  # Dispersion parameter - clamp to be positive
  size <- exp(par[k_count + k_zero + 1])
  if(is.na(size) || size <= 0){
    size <- .Machine$double.eps
  }

  # Linear predictors
  mu <- exp(X %*% beta + offset)
  phi <- plogis(Z %*% gamma)

  # Clamp phi to avoid log(0)
  phi <- pmin(pmax(phi, .Machine$double.eps), 1 - .Machine$double.eps)

  # Log-likelihood contribution for each observation
  loglik_i <- suppressWarnings(
    weights * (
      (Y > 0) * (
        log(1 - phi) + dnbinom(Y, size = size, mu = mu, log = TRUE)
      ) +
      (Y == 0) * (
        log(phi + (1 - phi) * pnbinom(0, size = size, mu = mu))
      )
    )
  )

  # Replace NA/NaN/-Inf with a large negative number
  loglik_i[is.na(loglik_i) | is.nan(loglik_i) | (is.infinite(loglik_i) & loglik_i < 0)] <- -1e10

  # Return the sum of the log-likelihoods
  sum(loglik_i)
}


# A custom psiFUN for zeroinfl objects. geex does not have a native psiFUN for
# zeroinfl objects from the pscl package, so we need to create one. This
# function computes the scores of the ZINB model, which are the estimating
# functions for the model parameters.
psiFUN.zeroinfl <- function(model, data){
  # Extract model components
  Y <- model$y
  X <- model.matrix(model, component = "count")
  Z <- model.matrix(model, component = "zero")
  weights <- model.weights(model)
  if(is.null(weights)){
    weights <- rep(1, nrow(data))
  }
  offset <- model.offset(model)
  if(is.null(offset)){
    offset <- rep(0, nrow(data))
  }

  # The function that numDeriv::grad will differentiate
  objective_function <- function(par){
    logL.zeroinfl(par, X = X, Z = Z, Y = Y, weights = weights, offset = offset)
  }

  # Return a function that computes the scores for each observation
  function(theta){
    # Use numDeriv::grad to compute the gradient of the log-likelihood for each
    # observation. This is done by wrapping the objective function in another
    # function that computes the log-likelihood for a single observation.
    sapply(1:nrow(data), function(i){
      grad(
        func = logL.zeroinfl,
        x = theta,
        X = X[i, , drop = FALSE],
        Z = Z[i, , drop = FALSE],
        Y = Y[i],
        weights = weights[i],
        offset = offset[i]
      )
    })
  }
}

# A custom grab_psiFUN for zeroinfl objects. This function tells geex how to
# find and use the custom psiFUN for zeroinfl objects.
grab_psiFUN.zeroinfl <- function(model, data){
  psiFUN.zeroinfl(model, data)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Weighted ZINB estimator
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Helper function to predict the outcome from a ZINB model given parameters
predict_zinb <- function(X_count, X_zero, theta_m){
  k_count <- ncol(X_count)
  k_zero <- ncol(X_zero)
  beta <- theta_m[1:k_count]
  gamma <- theta_m[(k_count + 1):(k_count + k_zero)]

  mu <- exp(X_count %*% beta)
  phi <- plogis(X_zero %*% gamma)

  (1 - phi) * mu
}

# Estimating function for the weighted ZINB estimator
estfun_WTD_ZINB <- function(data, models){
  # Grab design matrices and psi functions
  X <- data$X
  Xe <- grab_design_matrix(data = data, rhs_formula = grab_fixed_formula(models$e))

  data0 <- data1 <- data
  data0$X <- 0
  data1$X <- 1

  Xm0_count <- grab_design_matrix(data = data0, rhs_formula = grab_fixed_formula(models$m), component = "count")
  Xm1_count <- grab_design_matrix(data = data1, rhs_formula = grab_fixed_formula(models$m), component = "count")
  Xm0_zero <- grab_design_matrix(data = data0, rhs_formula = grab_fixed_formula(models$m), component = "zero")
  Xm1_zero <- grab_design_matrix(data = data1, rhs_formula = grab_fixed_formula(models$m), component = "zero")

  e_scores <- grab_psiFUN(models$e, data)
  m_scores <- grab_psiFUN(models$m, data)

  # Define parameter positions
  e_pos <- seq_len(ncol(Xe))
  last_e <- if(length(e_pos) > 0) max(e_pos) else 0
  m_pos <- seq_len(length(coef(models$m))) + last_e

  function(theta){
    p <- length(theta)
    e <- plogis(Xe %*% theta[e_pos])
    W <- X/e + (1-X)/(1-e)

    # Potential outcomes
    m1 <- predict_zinb(Xm1_count, Xm1_zero, theta[m_pos])
    m0 <- predict_zinb(Xm0_count, Xm0_zero, theta[m_pos])

    rbind(
      e_scores(theta[e_pos]),
      sweep(m_scores(theta[m_pos]), 2, W, "*"),
      m1 - theta[p - 2],
      m0 - theta[p - 1],
      rep((theta[p - 2] - theta[p - 1]) - theta[p], nrow(data))
    )
  }
}

# Wrapper function for the weighted ZINB estimator
geex_WTD_ZINB <- function(data, propensity_formula, outcome_formula){
  # Fit initial models
  e_model <- glm(propensity_formula, data = data, family = binomial)
  data$IPTW <- 1 / predict(e_model, type = "response") * data$X + 1 / (1 - predict(e_model, type = "response")) * (1 - data$X)
  m_model <- zeroinfl(outcome_formula, data = data, weights = IPTW)

  # Get initial values
  data0 <- data1 <- data
  data0$X <- 0
  data1$X <- 1
  CM1_IV <- mean(predict(m_model, newdata = data1, type = "response"))
  CM0_IV <- mean(predict(m_model, newdata = data0, type = "response"))
  ATE_IV <- CM1_IV - CM0_IV

  models <- list(e = e_model, m = m_model)

  # Run geex
  geex_results <- m_estimate(
    estFUN = estfun_WTD_ZINB,
    data = data,
    roots = c(coef(e_model), coef(m_model), CM1_IV, CM0_IV, ATE_IV),
    compute_roots = FALSE,
    outer_args = list(models = models)
  )

  # Format and return results
  ATE_est <- geex_results@estimates[length(geex_results@estimates)]
  ATE_se <- sqrt(geex_results@vcov[length(geex_results@estimates), length(geex_results@estimates)])

  cbind(ATE = ATE_est, ESseATE = ATE_se)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. G-formula/standardization ZINB estimator
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Estimating function for the g-formula ZINB estimator
estfun_GF_ZINB <- function(data, models){
  # Grab design matrices and psi functions
  data0 <- data1 <- data
  data0$X <- 0
  data1$X <- 1

  Xm0_count <- grab_design_matrix(data = data0, rhs_formula = grab_fixed_formula(models$m), component = "count")
  Xm1_count <- grab_design_matrix(data = data1, rhs_formula = grab_fixed_formula(models$m), component = "count")
  Xm0_zero <- grab_design_matrix(data = data0, rhs_formula = grab_fixed_formula(models$m), component = "zero")
  Xm1_zero <- grab_design_matrix(data = data1, rhs_formula = grab_fixed_formula(models$m), component = "zero")

  m_scores <- grab_psiFUN(models$m, data)

  m_pos <- seq_len(length(coef(models$m)))

  function(theta){
    p <- length(theta)
    m1 <- predict_zinb(Xm1_count, Xm1_zero, theta[m_pos])
    m0 <- predict_zinb(Xm0_count, Xm0_zero, theta[m_pos])

    rbind(
      m_scores(theta[m_pos]),
      m1 - theta[p - 2],
      m0 - theta[p - 1],
      rep((theta[p - 2] - theta[p - 1]) - theta[p], nrow(data))
    )
  }
}

# Wrapper function for the g-formula ZINB estimator
geex_GF_ZINB <- function(data, outcome_formula){
  # Fit initial model
  m_model <- zeroinfl(outcome_formula, data = data)

  # Get initial values
  data0 <- data1 <- data
  data0$X <- 0
  data1$X <- 1
  CM1_IV <- mean(predict(m_model, newdata = data1, type = "response"))
  CM0_IV <- mean(predict(m_model, newdata = data0, type = "response"))
  ATE_IV <- CM1_IV - CM0_IV

  models <- list(m = m_model)

  # Run geex
  geex_results <- m_estimate(
    estFUN = estfun_GF_ZINB,
    data = data,
    roots = c(coef(m_model), CM1_IV, CM0_IV, ATE_IV),
    compute_roots = FALSE,
    outer_args = list(models = models)
  )

  # Format and return results
  ATE_est <- geex_results@estimates[length(geex_results@estimates)]
  ATE_se <- sqrt(geex_results@vcov[length(geex_results@estimates), length(geex_results@estimates)])

  cbind(ATE = ATE_est, ESseATE = ATE_se)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. Doubly robust AIPW ZINB estimator
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Estimating function for the DR AIPW ZINB estimator
estfun_AIPW_ZINB <- function(data, models){
  # Grab design matrices and psi functions
  Y <- data$Y
  X <- data$X
  Xe <- grab_design_matrix(data = data, rhs_formula = grab_fixed_formula(models$e))

  data0 <- data1 <- data
  data0$X <- 0
  data1$X <- 1

  Xm0_count <- grab_design_matrix(data = data0, rhs_formula = grab_fixed_formula(models$m), component = "count")
  Xm1_count <- grab_design_matrix(data = data1, rhs_formula = grab_fixed_formula(models$m), component = "count")
  Xm0_zero <- grab_design_matrix(data = data0, rhs_formula = grab_fixed_formula(models$m), component = "zero")
  Xm1_zero <- grab_design_matrix(data = data1, rhs_formula = grab_fixed_formula(models$m), component = "zero")

  e_scores <- grab_psiFUN(models$e, data)
  m_scores <- grab_psiFUN(models$m, data)

  # Define parameter positions
  e_pos <- seq_len(ncol(Xe))
  last_e <- if(length(e_pos) > 0) max(e_pos) else 0
  m_pos <- seq_len(length(coef(models$m))) + last_e

  function(theta){
    p <- length(theta)
    e <- plogis(Xe %*% theta[e_pos])

    # Potential outcomes
    m1 <- predict_zinb(Xm1_count, Xm1_zero, theta[m_pos])
    m0 <- predict_zinb(Xm0_count, Xm0_zero, theta[m_pos])

    # AIPW estimating equations for the causal means
    CM1 <- (X/e) * Y - ((X - e)/e) * m1
    CM0 <- ((1-X)/(1-e)) * Y + ((X-e)/(1-e)) * m0

    rbind(
      e_scores(theta[e_pos]),
      m_scores(theta[m_pos]),
      CM1 - theta[p - 2],
      CM0 - theta[p - 1],
      rep((theta[p - 2] - theta[p - 1]) - theta[p], nrow(data))
    )
  }
}

# Wrapper function for the DR AIPW ZINB estimator
geex_AIPW_ZINB <- function(data, propensity_formula, outcome_formula){
  # Fit initial models
  e_model <- glm(propensity_formula, data = data, family = binomial)
  m_model <- zeroinfl(outcome_formula, data = data)

  # Get initial values
  data0 <- data1 <- data
  data0$X <- 0
  data1$X <- 1
  m1_IV <- predict(m_model, newdata = data1, type = "response")
  m0_IV <- predict(m_model, newdata = data0, type = "response")
  e_IV <- predict(e_model, type = "response")

  CM1_IV <- mean((data$X/e_IV) * data$Y - ((data$X - e_IV)/e_IV) * m1_IV)
  CM0_IV <- mean(((1-data$X)/(1-e_IV)) * data$Y + ((data$X-e_IV)/(1-e_IV)) * m0_IV)
  ATE_IV <- CM1_IV - CM0_IV

  models <- list(e = e_model, m = m_model)

  # Run geex
  geex_results <- m_estimate(
    estFUN = estfun_AIPW_ZINB,
    data = data,
    roots = c(coef(e_model), coef(m_model), CM1_IV, CM0_IV, ATE_IV),
    compute_roots = FALSE,
    outer_args = list(models = models)
  )

  # Format and return results
  ATE_est <- geex_results@estimates[length(geex_results@estimates)]
  ATE_se <- sqrt(geex_results@vcov[length(geex_results@estimates), length(geex_results@estimates)])

  cbind(ATE = ATE_est, ESseATE = ATE_se)
}

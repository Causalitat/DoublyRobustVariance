# Program: 07_Double_Robustness_Fast.R
# Purpose: Direct Calculation of Bias for DR and WTD Estimators (No variance/GEEX)

library(MASS)
library(pscl)

set.seed(42)
n <- 2000 

# ------------------------------------------------------------------------------
# Data Generation
# ------------------------------------------------------------------------------
generate_data <- function(n) {
  Z1 <- rnorm(n)
  Z2 <- rnorm(n)
  Z2_sq <- Z2^2
  
  # Truth
  p_X <- plogis(-0.5 + 0.8 * Z1 - 0.5 * Z2_sq)
  X <- rbinom(n, 1, p_X)
  
  logit_p_pos <- -0.5 + 0.5 * X + 0.5 * Z1 + 0.3 * Z2_sq
  p_pos <- plogis(logit_p_pos)
  is_positive <- rbinom(n, 1, p_pos)
  
  log_mu <- 1.0 + 0.5 * X + 0.5 * Z1 + 0.2 * Z2_sq
  mu <- exp(log_mu)
  theta_true <- 2.0
  
  Y <- rep(0, n)
  for(i in 1:n) {
    if(is_positive[i] == 1) {
      val <- 0
      while(val == 0) val <- rnbinom(1, size = theta_true, mu = mu[i])
      Y[i] <- val
    }
  }
  return(data.frame(Y=Y, X=X, Z1=Z1, Z2=Z2, Z2_sq=Z2_sq))
}

TRUE_ATE <- 2.12

# ------------------------------------------------------------------------------
# Estimation Helper
# ------------------------------------------------------------------------------
calc_estimates <- function(dat, f_ps, f_zero, f_count) {
  
  # 1. Fit Propensity
  m_ps <- glm(f_ps, data = dat, family = binomial)
  e_hat <- predict(m_ps, type = "response")
  e_hat <- pmax(pmin(e_hat, 1 - 1e-6), 1e-6)
  
  # 2. Fit Outcome Models (Standard / Unweighted)
  dat$Y_bin <- as.numeric(dat$Y > 0)
  m_zero <- glm(f_zero, data = dat, family = binomial)
  
  dat_pos <- dat[dat$Y > 0, ]
  m_count <- tryCatch({
    glm.nb(f_count, data = dat_pos)
  }, error = function(e) glm(f_count, data = dat_pos, family = poisson))
  
  if(inherits(m_count, "negbin")) theta_est <- m_count$theta else theta_est <- 10000
  
  # 3. Fit Weighted Outcome Models
  dat$wts <- (dat$X / e_hat) + ((1 - dat$X) / (1 - e_hat))
  m_zero_wtd <- glm(f_zero, data = dat, family = binomial, weights = wts)
  
  weights_pos <- dat$wts[dat$Y > 0]
  dat_pos$weights_pos <- weights_pos
  m_count_wtd <- tryCatch({
    glm.nb(f_count, data = dat_pos, weights = weights_pos)
  }, error = function(e) glm(f_count, data = dat_pos, family = poisson, weights = weights_pos))
  
  if(inherits(m_count_wtd, "negbin")) theta_wtd <- m_count_wtd$theta else theta_wtd <- 10000
  
  # 4. Predictions Helper
  mean_ztnb <- function(X_mat, beta, th) {
     mu <- exp(X_mat %*% beta)
     p0 <- (th/(th+mu))^th
     return(mu/(1-p0))
  }
  
  dat1 <- dat; dat1$X <- 1
  dat0 <- dat; dat0$X <- 0
  
  # Matrices
  Xz1 <- model.matrix(f_zero, dat1); Xz0 <- model.matrix(f_zero, dat0)
  Xc1 <- model.matrix(f_count, dat1); Xc0 <- model.matrix(f_count, dat0)
  
  # --- DR Predictions (Unweighted Parameters) ---
  m1_z_dr <- plogis(Xz1 %*% coef(m_zero))
  m0_z_dr <- plogis(Xz0 %*% coef(m_zero))
  m1_c_dr <- mean_ztnb(Xc1, coef(m_count), theta_est)
  m0_c_dr <- mean_ztnb(Xc0, coef(m_count), theta_est)
  
  m1_dr_hat <- m1_z_dr * m1_c_dr
  m0_dr_hat <- m0_z_dr * m0_c_dr
  
  # --- WTD Predictions (Weighted Parameters) ---
  m1_z_wtd <- plogis(Xz1 %*% coef(m_zero_wtd))
  m0_z_wtd <- plogis(Xz0 %*% coef(m_zero_wtd))
  m1_c_wtd <- mean_ztnb(Xc1, coef(m_count_wtd), theta_wtd)
  m0_c_wtd <- mean_ztnb(Xc0, coef(m_count_wtd), theta_wtd)
  
  m1_wtd_hat <- m1_z_wtd * m1_c_wtd
  m0_wtd_hat <- m0_z_wtd * m0_c_wtd
  
  # 5. Calculation
  # DR (AIPW)
  # mu1 = mean( (A*Y)/e - ((A-e)/e)*m1 )
  dr_mu1 <- mean( (dat$X * dat$Y)/e_hat - ((dat$X - e_hat)/e_hat)*m1_dr_hat )
  dr_mu0 <- mean( ((1 - dat$X) * dat$Y)/(1 - e_hat) + ((dat$X - e_hat)/(1 - e_hat))*m0_dr_hat )
  dr_ate <- dr_mu1 - dr_mu0
  
  # WTD (Standardization on Weighted Model)
  wtd_mu1 <- mean(m1_wtd_hat)
  wtd_mu0 <- mean(m0_wtd_hat)
  wtd_ate <- wtd_mu1 - wtd_mu0
  
  return(c(dr_ate, wtd_ate))
}


cat("Running Fast Check (Direct Calculation)...\n")
dat <- generate_data(n)

f_ps_corr <- X ~ Z1 + I(Z2^2)
f_ps_wrong <- X ~ Z1 + Z2
f_zero_corr <- Y_bin ~ X + Z1 + I(Z2^2)
f_zero_wrong <- Y_bin ~ X + Z1 + Z2
f_count_corr <- Y ~ X + Z1 + I(Z2^2)
f_count_wrong <- Y ~ X + Z1 + Z2

# 1. Both Correct
res1 <- calc_estimates(dat, f_ps_corr, f_zero_corr, f_count_corr)

# 2. Outcome Wrong
res2 <- calc_estimates(dat, f_ps_corr, f_zero_wrong, f_count_wrong)

# 3. Propensity Wrong
res3 <- calc_estimates(dat, f_ps_wrong, f_zero_corr, f_count_corr)

# 4. Both Wrong
res4 <- calc_estimates(dat, f_ps_wrong, f_zero_wrong, f_count_wrong)

cat("\nRESULTS (Bias Check)\n")
cat("Truth:", TRUE_ATE, "\n")
cat(sprintf("%-20s %-10s %-10s %-10s %-10s\n", "Scenario", "DR Est", "WTD Est", "DR Bias", "WTD Bias"))
cat(sprintf("%-20s %-10.4f %-10.4f %-10.4f %-10.4f\n", "1. Both Correct", res1[1], res1[2], res1[1]-TRUE_ATE, res1[2]-TRUE_ATE))
cat(sprintf("%-20s %-10.4f %-10.4f %-10.4f %-10.4f\n", "2. Outcome Wrong", res2[1], res2[2], res2[1]-TRUE_ATE, res2[2]-TRUE_ATE))
cat(sprintf("%-20s %-10.4f %-10.4f %-10.4f %-10.4f\n", "3. Propensity Wrong", res3[1], res3[2], res3[1]-TRUE_ATE, res3[2]-TRUE_ATE))
cat(sprintf("%-20s %-10.4f %-10.4f %-10.4f %-10.4f\n", "4. Both Wrong", res4[1], res4[2], res4[1]-TRUE_ATE, res4[2]-TRUE_ATE))

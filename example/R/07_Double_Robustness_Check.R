# Program: 07_Double_Robustness_Check.R
# Purpose: Verify Double Robustness of Hurdle Estimators
# We test 3 scenarios:
# 1. Both Models Correct
# 2. Outcome Model Misspecified (Propensity Correct)
# 3. Propensity Model Misspecified (Outcome Correct)

library(geex)
library(MASS)
library(pscl)
library(numDeriv)

source("04_Hurdle_NegBin_Functions.R")

set.seed(42)
n <- 500
n_sim <- 50 # Number of simulations

# ------------------------------------------------------------------------------
# Data Generation
# ------------------------------------------------------------------------------
generate_data <- function(n) {
  # Covariates
  Z1 <- rnorm(n) # Linear relevant
  Z2 <- rnorm(n) # Non-linear relevant (Z2^2)
  Z2_sq <- Z2^2
  
  # Treatment Assignment (Truth depends on Z1 and Z2_sq)
  # logit(e) = -0.5 + 0.8*Z1 - 0.5*Z2_sq
  p_X <- plogis(-0.5 + 0.8 * Z1 - 0.5 * Z2_sq)
  X <- rbinom(n, 1, p_X)
  
  # Outcome Model (Truth depends on Z1 and Z2_sq)
  # Zero part
  logit_p_pos <- -0.5 + 0.5 * X + 0.5 * Z1 + 0.3 * Z2_sq
  p_pos <- plogis(logit_p_pos)
  is_positive <- rbinom(n, 1, p_pos)
  
  # Count part
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
  
  # True ATE (Approximate via super-population)
  # We can calculate this by setting X=1 and X=0 for everyone
  # For speed, we just use the realized sample truth or a large N approx
  return(data.frame(Y=Y, X=X, Z1=Z1, Z2=Z2, Z2_sq=Z2_sq))
}

# Calculate True ATE via large sample
calc_true_ate <- function() {
  d_large <- generate_data(20000)
  
  # Function to get mean
  get_mean <- function(data_in, x_val) {
    mu <- exp(1.0 + 0.5 * x_val + 0.5 * data_in$Z1 + 0.2 * data_in$Z2_sq)
    p_pos <- plogis(-0.5 + 0.5 * x_val + 0.5 * data_in$Z1 + 0.3 * data_in$Z2_sq)
    theta <- 2.0
    p0 <- (theta / (theta + mu))^theta
    mu_zt <- mu / (1 - p0)
    return(mean(p_pos * mu_zt))
  }
  
  m1 <- get_mean(d_large, 1)
  m0 <- get_mean(d_large, 0)
  return(m1 - m0)
}
TRUE_ATE <- calc_true_ate()
cat("True ATE:", TRUE_ATE, "\n")


# ------------------------------------------------------------------------------
# Simulation Loop
# ------------------------------------------------------------------------------
run_sim <- function(n_sim, n) {
  results <- data.frame()
  
  # Formulas
  # Correct: Includes Z2_sq (or I(Z2^2))
  # Incorrect: Only includes Z2 (Linear)
  
  f_ps_corr <- X ~ Z1 + I(Z2^2)
  f_ps_wrong <- X ~ Z1 + Z2
  
  f_zero_corr <- as.numeric(Y>0) ~ X + Z1 + I(Z2^2)
  f_zero_wrong <- as.numeric(Y>0) ~ X + Z1 + Z2
  
  f_count_corr <- Y ~ X + Z1 + I(Z2^2)
  f_count_wrong <- Y ~ X + Z1 + Z2
  
  for(i in 1:n_sim) {
    dat <- generate_data(n)
    
    # helper wrapper to run and catch
    run_est <- function(func, ps_f, z_f, c_f) {
       tryCatch({
         r <- func(dat, ps_f, z_f, c_f)
         return(r$result$ATE)
       }, error = function(e) NA)
    }

    # 1. Both Correct
    dr_corr_corr <- run_est(geex_Hurdle_NB, f_ps_corr, f_zero_corr, f_count_corr)
    wtd_corr_corr <- run_est(geex_WTD_Hurdle_NB, f_ps_corr, f_zero_corr, f_count_corr)
    
    # 2. Outcome Wrong (Propensity Correct)
    dr_corr_wrong <- run_est(geex_Hurdle_NB, f_ps_corr, f_zero_wrong, f_count_wrong)
    wtd_corr_wrong <- run_est(geex_WTD_Hurdle_NB, f_ps_corr, f_zero_wrong, f_count_wrong)
    
    # 3. Propensity Wrong (Outcome Correct)
    dr_wrong_corr <- run_est(geex_Hurdle_NB, f_ps_wrong, f_zero_corr, f_count_corr)
    wtd_wrong_corr <- run_est(geex_WTD_Hurdle_NB, f_ps_wrong, f_zero_corr, f_count_corr)
    
    # 4. Both Wrong
    dr_wrong_wrong <- run_est(geex_Hurdle_NB, f_ps_wrong, f_zero_wrong, f_count_wrong)
    wtd_wrong_wrong <- run_est(geex_WTD_Hurdle_NB, f_ps_wrong, f_zero_wrong, f_count_wrong)
    
    results <- rbind(results, data.frame(
      sim=i,
      dr_cc=dr_corr_corr, wtd_cc=wtd_corr_corr,
      dr_cw=dr_corr_wrong, wtd_cw=wtd_corr_wrong,
      dr_wc=dr_wrong_corr, wtd_wc=wtd_wrong_corr,
      dr_ww=dr_wrong_wrong, wtd_ww=wtd_wrong_wrong
    ))
    if(i %% 10 == 0) cat(".")
  }
  cat("\n")
  return(results)
}

cat("Running Simulation...\n")
res_df <- run_sim(n_sim, n)

# ------------------------------------------------------------------------------
# Analysis
# ------------------------------------------------------------------------------
analyze_bias <- function(vec, truth) {
  mean(vec - truth, na.rm=TRUE)
}
analyze_rmse <- function(vec, truth) {
  sqrt(mean((vec - truth)^2, na.rm=TRUE))
}

cat("\nSimulation Results (N =", n, ", Reps =", n_sim, ")\n")
cat("True ATE:", round(TRUE_ATE, 4), "\n\n")

cat(sprintf("%-25s %-10s %-10s\n", "Scenario", "DR Bias", "WTD Bias"))
cat(sprintf("%-25s %-10.4f %-10.4f\n", "1. Both Correct", 
            analyze_bias(res_df$dr_cc, TRUE_ATE), analyze_bias(res_df$wtd_cc, TRUE_ATE)))
cat(sprintf("%-25s %-10.4f %-10.4f\n", "2. Outcome Wrong", 
            analyze_bias(res_df$dr_cw, TRUE_ATE), analyze_bias(res_df$wtd_cw, TRUE_ATE)))
cat(sprintf("%-25s %-10.4f %-10.4f\n", "3. Propensity Wrong", 
            analyze_bias(res_df$dr_wc, TRUE_ATE), analyze_bias(res_df$wtd_wc, TRUE_ATE)))
cat(sprintf("%-25s %-10.4f %-10.4f\n", "4. Both Wrong", 
            analyze_bias(res_df$dr_ww, TRUE_ATE), analyze_bias(res_df$wtd_ww, TRUE_ATE)))

cat("\nInterpretations:\n")
cat("- 'Outcome Wrong': DR should have low bias (Robust). WTD should have low bias (depends on weights).\n")
cat("- 'Propensity Wrong': DR should have low bias (Robust). WTD should have HIGH bias.\n")

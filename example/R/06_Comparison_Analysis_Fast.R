# Program: 06_Comparison_Analysis_Fast.R
# Purpose: FAST Compare Manual ZTNB Est (04) vs Refactored PSCL ZTNB Est (03)
# N = 100 for speed

library(geex)
library(boot)
library(rootSolve)
library(MASS)
library(pscl)
library(numDeriv)

source("04_Hurdle_NegBin_Functions.R")
source("03_Hurdle_Functions.R")

set.seed(12345)
n <- 100

# ------------------------------------------------------------------------------
# 1. Generate Synthetic Data
# ------------------------------------------------------------------------------
Z1 <- rnorm(n)
Z2 <- rbinom(n, 1, 0.5)
p_X <- plogis(0.5 * Z1 - 0.5 * Z2)
X <- rbinom(n, 1, p_X)

logit_p_pos <- -0.5 + 1.0 * X + 0.5 * Z1 + 0.5 * Z2
p_pos <- plogis(logit_p_pos)
is_positive <- rbinom(n, 1, p_pos)

log_mu <- 1.0 + 0.5 * X + 0.2 * Z1
mu <- exp(log_mu)
theta_true <- 2.0 

Y_count <- rep(0, n)
for(i in 1:n) {
  if(is_positive[i] == 1) {
    val <- 0
    while(val == 0) val <- rnbinom(1, size = theta_true, mu = mu[i])
    Y_count[i] <- val
  }
}
Y <- Y_count
data <- data.frame(Y=Y, X=X, Z1=Z1, Z2=Z2)
data$Y_bin <- as.numeric(data$Y > 0) # Pre-calculate for fallback

ps_formula <- X ~ Z1 + Z2
zero_formula <- Y_bin ~ X + Z1 + Z2
count_formula <- Y ~ X + Z1
outcome_formula_pscl <- Y ~ X + Z1 | X + Z1 + Z2

cat("Running Comparison (N=100)...\n")

# ------------------------------------------------------------------------------
# 2. Run Manual Estimator (04)
# ------------------------------------------------------------------------------
res_04_est <- NA
res_04_se <- NA

tryCatch({
  res_04 <- geex_Hurdle_NB(data, ps_formula, zero_formula, count_formula)
  res_04_est <- res_04$result$ATE
  res_04_se <- res_04$result$SE
}, error = function(e) {
  # Manual Fallback
  e_model <- glm(ps_formula, data = data, family = binomial)
  m_zero <- glm(zero_formula, data = data, family = binomial)
  data_pos <- data[data$Y > 0, ]
  m_count_nb <- tryCatch({
    glm.nb(count_formula, data = data_pos)
  }, error = function(e) glm(count_formula, data = data_pos, family = poisson))
  
  theta_est <- if(inherits(m_count_nb, "negbin")) m_count_nb$theta else 10000
  root_c <- coef(m_count_nb)
  
  mean_ztnb_init <- function(X_m, b_c, th) {
       mu <- exp(X_m %*% b_c); p0 <- (th/(th+mu))^th
       return(mu/(1-p0))
  }
  
  dat1 <- data; dat1$X <- 1; dat0 <- data; dat0$X <- 0
  Xc1 <- model.matrix(count_formula, dat1); Xz1 <- model.matrix(zero_formula, dat1)
  Xc0 <- model.matrix(count_formula, dat0); Xz0 <- model.matrix(zero_formula, dat0)
  
  m1_vec <- plogis(Xz1 %*% coef(m_zero)) * mean_ztnb_init(Xc1, root_c, theta_est)
  m0_vec <- plogis(Xz0 %*% coef(m_zero)) * mean_ztnb_init(Xc0, root_c, theta_est)
  m1_start <- mean(m1_vec); m0_start <- mean(m0_vec)
  
  roots <- c(coef(e_model), coef(m_zero), root_c, m1_start, m0_start, m1_start - m0_start)
  roots <- unname(roots)
  m_count_dummy <- glm(count_formula, data = data, family = poisson)
  models <- list(e = e_model, zero = m_zero, count = m_count_dummy, theta_fixed = theta_est)
  
  psi_fun <- estfun_Hurdle_NB(data, models = models)
  grad_fun <- function(theta) colSums(psi_fun(theta))
  A <- jacobian(grad_fun, roots)
  B <- crossprod(psi_fun(roots))
  Sigma <- solve(A) %*% B %*% solve(t(A))
  n_p <- length(roots)
  res_04_est <<- roots[n_p]
  res_04_se <<- sqrt(Sigma[n_p, n_p])
})

# ------------------------------------------------------------------------------
# 3. Run Refactored PSCL Estimator (03)
# ------------------------------------------------------------------------------
res_03_est <- NA
res_03_se <- NA

tryCatch({
  res_03 <- geex_AIPW_hurdle(data, ps_formula, outcome_formula_pscl, dist = "negbin")
  res_03_est <- res_03[1,1]
  res_03_se <- res_03[1,2]
}, error = function(e) print(e))

# ------------------------------------------------------------------------------
# 4. Run Weighted Manual Estimator (04 IPW)
# ------------------------------------------------------------------------------
res_04_wtd_est <- NA
res_04_wtd_se <- NA

tryCatch({
  res_04_wtd <- geex_WTD_Hurdle_NB(data, ps_formula, zero_formula, count_formula)
  res_04_wtd_est <- res_04_wtd$result$ATE
  res_04_wtd_se <- res_04_wtd$result$SE
}, error = function(e) {
  # Manual Fallback for WTD
  # 1. Fit Propensity
  m_ps <- glm(ps_formula, data = data, family = binomial)
  ps_pred <- predict(m_ps, type="response")
  W <- (data$X / ps_pred) + ((1 - data$X) / (1 - ps_pred))
  
  # 2. Fit Weighted Models
  data$W <- W
  m_zero <- glm(zero_formula, data = data, family = binomial, weights = W)
  
  data_pos <- data[data$Y > 0, ]
  # glm.nb looks for weights in data if provided
  m_count_nb <- tryCatch({
    glm.nb(count_formula, data = data_pos, weights = W)
  }, error = function(e) glm(count_formula, data = data_pos, family = poisson, weights = W))
  
  theta_est <- if(inherits(m_count_nb, "negbin")) m_count_nb$theta else 10000
  root_c <- coef(m_count_nb)
  
  # Means
  mean_ztnb_init <- function(X_m, b_c, th) {
       mu <- exp(X_m %*% b_c); p0 <- (th/(th+mu))^th
       return(mu/(1-p0))
  }
  
  dat1 <- data; dat1$X <- 1; dat0 <- data; dat0$X <- 0
  Xc1 <- model.matrix(count_formula, dat1); Xz1 <- model.matrix(zero_formula, dat1)
  Xc0 <- model.matrix(count_formula, dat0); Xz0 <- model.matrix(zero_formula, dat0)
  
  m1_vec <- plogis(Xz1 %*% coef(m_zero)) * mean_ztnb_init(Xc1, root_c, theta_est)
  m0_vec <- plogis(Xz0 %*% coef(m_zero)) * mean_ztnb_init(Xc0, root_c, theta_est)
  m1_start <- mean(m1_vec); m0_start <- mean(m0_vec)
  
  roots <- c(coef(m_ps), coef(m_zero), root_c, m1_start, m0_start, m1_start - m0_start)
  roots <- unname(roots)
  m_count_dummy <- glm(count_formula, data = data, family = poisson)
  models <- list(e = m_ps, zero = m_zero, count = m_count_dummy, theta_fixed = theta_est)
  
  psi_fun <- estfun_WTD_Hurdle_NB(data, models = models)
  grad_fun <- function(theta) colSums(psi_fun(theta))
  A <- jacobian(grad_fun, roots)
  B <- crossprod(psi_fun(roots))
  Sigma <- solve(A) %*% B %*% solve(t(A))
  n_p <- length(roots)
  res_04_wtd_est <<- roots[n_p]
  res_04_wtd_se <<- sqrt(Sigma[n_p, n_p])
})


# ------------------------------------------------------------------------------
# 5. Final Comparison
# ------------------------------------------------------------------------------
cat(sprintf("%-20s %-10s %-10s\n", "Method", "ATE", "SE"))
cat(sprintf("%-20s %-10.4f %-10.4f\n", "1. Manual DR ZTNB", res_04_est, res_04_se))
cat(sprintf("%-20s %-10.4f %-10.4f\n", "2. PSCL DR ZTNB", res_03_est, res_03_se))
cat(sprintf("%-20s %-10.4f %-10.4f\n", "3. Manual WTD ZTNB", res_04_wtd_est, res_04_wtd_se))

cat("Comparison Complete.\n")

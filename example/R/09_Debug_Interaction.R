# Program: 09_Debug_Interaction.R
# Purpose: Manually calculate Bread (A) and Meat (B) to debug dimensions

library(geex)
library(MASS)
library(pscl)
library(numDeriv)

source("04_Hurdle_NegBin_Functions.R")

set.seed(123)
n <- 200
Z1 <- rnorm(n)
X <- rbinom(n, 1, 0.5)
mu <- exp(1 + 0.5*X + 0.5*Z1 + 0.5*X*Z1)
Y <- rnbinom(n, size=2, mu=mu) 
Y[rbinom(n, 1, 0.3) == 1] <- 0
data <- data.frame(Y=Y, X=X, Z1=Z1)

# Formulas
ps_f <- X ~ Z1
cov_f <- ~ Z1
# Interaction Formulas
str_zero <- paste("as.numeric(Y > 0) ~ X * (Z1)")
str_count <- paste("Y ~ X * (Z1)")
zero_f <- as.formula(str_zero)
count_f <- as.formula(str_count)

cat("Running Debug setup...\n")

# 1. Setup Models
m_ps <- glm(ps_f, data = data, family = binomial)
ps_pred <- predict(m_ps, type="response")
data$wts <- (data$X / ps_pred) + ((1 - data$X) / (1 - ps_pred))

data$Y_bin <- as.numeric(data$Y > 0)
m_zero <- glm(zero_f, data = data, family = binomial, weights = data$wts)

data_pos <- data[data$Y > 0, ]
weights_pos <- data_pos$wts
data_pos$weights_pos <- weights_pos
m_count_nb <- glm.nb(count_f, data = data_pos, weights = weights_pos)
theta_est <- m_count_nb$theta
root_c <- coef(m_count_nb)

# Roots
Xc_full <- model.matrix(count_f, data)
Xz_full <- model.matrix(zero_f, data)

mean_ztnb_init <- function(X_m, b_c, th) {
   mu <- exp(X_m %*% b_c) 
   p0 <- (th / (th + mu))^th
   return(mu / (1 - p0))
}

dat1 <- data; dat1$X <- 1
dat0 <- data; dat0$X <- 0
Xc1 <- model.matrix(count_f, dat1); Xz1 <- model.matrix(zero_f, dat1)
Xc0 <- model.matrix(count_f, dat0); Xz0 <- model.matrix(zero_f, dat0)

m1_vec <- plogis(Xz1 %*% coef(m_zero)) * mean_ztnb_init(Xc1, root_c, theta_est)
m0_vec <- plogis(Xz0 %*% coef(m_zero)) * mean_ztnb_init(Xc0, root_c, theta_est)
m1_start <- mean(m1_vec); m0_start <- mean(m0_vec)

roots <- c(coef(m_ps), coef(m_zero), root_c, m1_start, m0_start, m1_start - m0_start)
roots[is.na(roots)] <- 0
roots <- unname(roots)

# Models list
m_count_dummy <- glm(count_f, data = data, family = poisson)
models <- list(e = m_ps, zero = m_zero, count = m_count_dummy, theta_fixed = theta_est)

cat("Calculating ESTFUN...\n")
# Get the function
psi_fun_closure <- estfun_WTD_Hurdle_NB(data, models = models)

# 2. Calculate B (Meat)
cat("Calculating B (Crossprod)...\n")
scores <- psi_fun_closure(roots)
cat("Scores Dim:", dim(scores), "\n") # Should be 200 x 13

B <- crossprod(scores)
cat("B Dim:", dim(B), "\n")

# 3. Calculate A (Bread / Jacobian)
cat("Calculating A (Jacobian)...\n")
# geex sums columns before jacobian
grad_fun <- function(theta) {
  S <- psi_fun_closure(theta)
  return(colSums(S))
}

A <- jacobian(grad_fun, roots)
cat("A Dim:", dim(A), "\n")

# 4. Invert A
cat("Inverting A...\n")
Ainv <- tryCatch({ solve(A) }, error = function(e) { print(e); return(NULL) })

if(!is.null(Ainv)) {
  cat("Ainv Dim:", dim(Ainv), "\n")
  
  cat("Calculating V = Ainv %*% B...\n")
  V_part <- Ainv %*% B
  cat("Success.\n")
} else {
  cat("A inversion failed.\n")
}

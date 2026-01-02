# Program: 08_Interaction_Check.R
# Purpose: Verify the Interaction Wrapper geex_WTD_Hurdle_NB_Int works

library(geex)
library(MASS)
library(pscl)
library(numDeriv)

source("04_Hurdle_NegBin_Functions.R")

set.seed(123)
n <- 200

# Generate simple data with effect modification
# Treatment effect is larger when Z1 is high
Z1 <- rnorm(n)
X <- rbinom(n, 1, 0.5)
mu <- exp(1 + 0.5*X + 0.5*Z1 + 0.5*X*Z1) # Interaction in Truth
Y <- rnbinom(n, size=2, mu=mu) 
# Make some zeros
Y[rbinom(n, 1, 0.3) == 1] <- 0

data <- data.frame(Y=Y, X=X, Z1=Z1)

# Formulas
ps_f <- X ~ Z1
cov_f <- ~ Z1

cat("Running Interaction Model Check...\n")
res_int <- geex_WTD_Hurdle_NB_Int(data, ps_f, cov_f)

print(res_int$result)

cat("\nDone.\n")

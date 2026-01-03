# Program: 10_ZINB_Comparison.R
# Purpose: Compare legacy (00_Estimators) vs New (02_ZINB_Functions) ZINB estimators
# Replication of logic from noHeapingCount/05_ZINB.R

rm(list=ls())

# libraries
library(numDeriv)
library(boot)
library(geex)
library(MASS)
library(pscl)

# Source Legacy and New Functions
# Note: We must be careful about function name conflicts.
# Legacy file path needs to be correct relative to this script or absolute.
# We are in example/R. legacy is in noHeapingCount/
# Let's assume we run this from example/R

# Source NEW functions
source("example/R/02_ZINB_Functions.R")

# Source LEGACY functions (adjust path)
# We need to temporarily rename or isolate if there are conflicts.
# Checking names:
# Legacy: estfun_dr_ZINB, geex_DR
# New: estfun_AIPW_ZINB, geex_AIPW_ZINB
# No conflict in export names.
source("noHeapingCount/00_Estimators_11.03.21.R")

set.seed(123)
num.part <- 500
sim <- 1
log.RR <- 0.5

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Data Generation (from 05_ZINB.R)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# generate age of each woman at start of study
agevis <- as.vector(runif(num.part, min = 20, max = 40))
part.id <- seq(1:num.part)
simdat1 <- as.data.frame(cbind(part.id, agevis))

# create a drug use and SE rate for each woman 
simdat1$drugrate <- inv.logit(0.08 - (simdat1$agevis/100)) 
simdat1$druguse <- rbinom(num.part, 1, simdat1$drugrate) 
simdat1$sxexrate <- inv.logit(-2.9 - simdat1$agevis/100 + 1.2*simdat1$druguse)
simdat1$sxex <- rbinom(num.part, 1, simdat1$sxexrate) 

# incarceration (Treatment A)
# Mapping: inc -> X
simdat1$incrate <- inv.logit(-1.5 + (1 - simdat1$agevis/100) + 0.5*simdat1$druguse + 0.5*simdat1$sxex)
simdat1$inc <- rbinom(num.part, 1, simdat1$incrate)
simdat1$X <- simdat1$inc # Map for new functions

# create prob of excess zero
simdat1$exczero0p <- inv.logit(-2.5 + (simdat1$agevis/100) - 0.3*simdat1$druguse - 2*simdat1$sxex)
simdat1$exczero1p <- inv.logit(-2.5 + (simdat1$agevis/100) - 0.3*simdat1$druguse - 2*simdat1$sxex)
simdat1$exczero0 <- rbinom(num.part, 1, simdat1$exczero0p)
simdat1$exczero1 <- rbinom(num.part, 1, simdat1$exczero1p)

# create means for susceptible pop
# Legacy code uses log.RR to shift lambda1
simdat1$lambda0 <- exp(-1 - 0.005*simdat1$agevis + 0.7*simdat1$druguse + 3.5*simdat1$sxex)
simdat1$lambda1 <- exp(-1 - 0.005*simdat1$agevis + 0.7*simdat1$druguse + 3.5*simdat1$sxex + log.RR)

simdat1$SP0 <- ifelse(simdat1$exczero0==1, 0, rnegbin(num.part, simdat1$lambda0, theta=2))
simdat1$SP1 <- ifelse(simdat1$exczero1==1, 0, rnegbin(num.part, simdat1$lambda1, theta=2))

# assign observed outcome
simdat1$SP <- ifelse(simdat1$inc==0, simdat1$SP0, simdat1$SP1)
simdat1$Y <- simdat1$SP # Map for new functions

# age/10 transformation
simdat1$agevis10 <- simdat1$agevis/10

cat("Data Generated. Simulating Estimators...\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Legacy Estimator (geex_DR)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Need to setup arguments exactly like 05_ZINB.R

# Weights setup (Legacy requires this pre-calc for initial values)
denom <- glm(inc ~ druguse + agevis10 + sxex, data = simdat1, family = binomial("logit"))
simdat2 <- simdat1 # Copy
simdat2$denom.probs <- denom$fitted.values

# Initial model fit for ZINB params (mod.iv)
thetamod <- zeroinfl(SP ~ 1 + inc + druguse + sxex +  agevis10 | 1 + druguse + sxex +  agevis10, 
                     data = simdat2, dist = "negbin", link="logit")
init.ZINB <- c(coef(thetamod), log(thetamod$theta)) # Note: Legacy uses log(1/theta) sometimes? 
# In 00: log(1/thetamod$theta) on line 132.
# Let's check ML.ZINB function in 00 to see what it optims.
# 00 line 271: logodisp<-par[10]; exp(logodisp)*lam / (1 + exp(logodisp)*lam)
# This implies odisp is dispersion param 'alpha' where Var = mu + alpha*mu^2 ?
# pscl theta is usually 1/alpha (the 'size' param).
# So 1/theta = alpha.
# So log(1/theta) is correct for the manual ML in 00.

# Legacy init code:
init.ZINB.legacy <- c(coef(thetamod), log(1/thetamod$theta))

# We use the legacy ML function to find roots? 
# Actually 05_ZINB calls ML.ZINB(simdat2, init.ZINB) to optimize.
# We will skip the manual optim step if possible and just use thetamod pars if close enough,
# OR we try to run the legacy optim if it works.
mod.iv <- tryCatch(ML.ZINB(simdat2, init.ZINB.legacy), error=function(e) init.ZINB.legacy)

# Predict Counterfactuals (Legacy Logic)
# Legacy manually predicts:
CM0.iv.dr <- mean(exp(mod.iv[1]+mod.iv[3]*simdat2$druguse+mod.iv[4]*simdat2$sxex+mod.iv[5]*simdat2$agevis10)*(1-inv.logit(mod.iv[6]+mod.iv[7]*simdat2$druguse+mod.iv[8]*simdat2$sxex+mod.iv[9]*simdat2$agevis10)))
# No, that's PG logic.
# DR logic uses the AIPW formula with predictions:
pred0.dr <- exp(mod.iv[1]+mod.iv[3]*simdat2$druguse+mod.iv[4]*simdat2$sxex+mod.iv[5]*simdat2$agevis10)*(1-inv.logit(mod.iv[6]+mod.iv[7]*simdat2$druguse+mod.iv[8]*simdat2$sxex+mod.iv[9]*simdat2$agevis10))
pred1.dr <- exp(mod.iv[1]+mod.iv[2]+mod.iv[3]*simdat2$druguse+mod.iv[4]*simdat2$sxex+mod.iv[5]*simdat2$agevis10)*(1-inv.logit(mod.iv[6]+mod.iv[7]*simdat2$druguse+mod.iv[8]*simdat2$sxex+mod.iv[9]*simdat2$agevis10))

# AIPW Means
CM1.iv.dr <- mean(((simdat2$inc*simdat2$SP - (simdat2$inc-simdat2$denom.probs)*pred1.dr))/simdat2$denom.probs)
CM0.iv.dr <- mean((((1-simdat2$inc)*simdat2$SP + (simdat2$inc-simdat2$denom.probs)*pred0.dr))/(1-simdat2$denom.probs))
RR.iv.dr <- CM1.iv.dr/CM0.iv.dr

cat("Legacy Init Values: RR =", RR.iv.dr, "\n")

# Run Legacy Estimator
cat("Running Legacy geex_DR...\n")
res_legacy <- tryCatch({
  geex_DR(simdat2, estfun_dr_ZINB, inc ~ druguse + agevis10 + sxex, mod.iv, CM0.iv.dr, CM1.iv.dr, RR.iv.dr)
}, error = function(e) {
  cat("Legacy failed:", e$message, "\n")
  return(NULL)
})


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# New Estimator (geex_AIPW_ZINB)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\nRunning New geex_AIPW_ZINB...\n")

# Formulas
# Count part: ~ druguse + sxex + agevis10 (Legacy includes 'inc' as first term usually)
# Legacy count coefs: Intercept, inc, druguse, sxex, agevis10 (5 params)
# Legacy zero coefs: Intercept, druguse, sxex, agevis10 (4 params)
# Wait, legacy 'thetamod' formula:
# SP ~ 1 + inc + druguse + sxex + agevis10 | 1 + druguse + sxex + agevis10
# Note: In pscl, the first part is Count, second is Zero.
# So Count depends on Trt+Covs. Zero depends on Covs only.

ps_formula <- X ~ druguse + agevis10 + sxex
# New functions use 'X' and 'Y' columns we created.
# Outcome formula for pscl:
out_formula <- Y ~ X + druguse + sxex + agevis10 | druguse + sxex + agevis10

res_new <- geex_AIPW_ZINB(simdat1, ps_formula, out_formula)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Comparison
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("\n--- RESULTS COMPARISON ---\n")
if(!is.null(res_legacy)) {
  # Legacy returns S4 object
  est_leg <- res_legacy@estimates
  vcov_leg <- res_legacy@vcov
  p <- length(est_leg)
  cat("Legacy (RR): Est =", est_leg[p], " SE =", sqrt(vcov_leg[p, p]), "\n")
} else {
  cat("Legacy: Failed\n")
}

cat("New (RR):    Est =", res_new$Estimate[1], " SE =", res_new$SE[1], "\n")
cat("New (ATE):   Est =", res_new$Estimate[3], " SE =", res_new$SE[3], "\n")

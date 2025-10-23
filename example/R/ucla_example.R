# Program: ucla_example.R
# Purpose: This program contains the example from the UCLA statistics website.

# Loading required packages for the workshop
library(AER)        # Applied Econometrics with R
library(pscl)       # Political Science Computational Laboratory
library(MASS)       # Modern Applied Statistics with S
library(sjPlot)     # For visualizing regression models
library(lmtest)     # For hypothesis testing in linear models
library(emmeans)    # For estimated marginal means
library(performance) # For model performance evaluation
library(sandwich)   #Robust Covariance Matrix Estimators
library(ggplot2)
library(dplyr)
library(tidyr)

# Load the data in R from AER
data("NMES1988")

# Select variables used in the workshop
dat <- NMES1988[, c(1, 6:8, 13, 15, 18)]

# Poisson Regression model for number of visits
formula <- visits ~ hospital + health + chronic + gender + school + insurance
m.pois <- glm(formula = formula,
                    family  = poisson(link = "log"),
                    data    = dat)

# Negative binomial Model for number of visits
m.nb <- glm.nb(formula, data = dat, link = log)

# Zero-Inflated Model for the number of visits
m.zeroinf.pois <- zeroinfl(
  visits ~ hospital + health + chronic + gender + school + insurance | #the count component
                             chronic + insurance + school + gender, #the zero component
                           data = dat, dist = "poisson")
m.zeroinf.nb <- zeroinfl(visits ~ hospital + health + chronic + gender + school + insurance |
                           chronic + insurance + school + gender,
                           data = dat, dist = "negbin")

# Hurdle Model for the number of visits
m.hurdle.pois <- hurdle(
  visits ~ hospital + health + chronic + gender + school + insurance,
  data = dat, dist = "poisson")
m.hurdle.nb <- hurdle(
  visits ~ hospital + health + chronic + gender + school + insurance,
  data = dat, dist = "negbin")

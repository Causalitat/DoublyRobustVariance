# Program: fish_data_exercise.R
# Purpose: This program contains the fish data exercise from the UCLA statistics website.

# Load required libraries
library(pscl)
library(MASS)
library(sjPlot)

# Read the data from idre website
ex.data <- read.csv("https://stats.idre.ucla.edu/stat/data/fish.csv")

# Make the variable 'camper' a factor since it is binary
ex.data$camper <- factor(ex.data$camper, labels = c("no", "yes"))

# Fit the negative binomial model
ex.m.nb <- glm.nb(count ~ child + camper + persons, data = ex.data)

# Fit the zero-inflated negative binomial model
ex.m.zeroinf.nb <- zeroinfl(count ~ child + camper | persons,
                            dist = "negbin", data = ex.data)

# Fit the hurdle negative binomial model
ex.m.hurdle.nb <- hurdle(count ~ child + camper + persons |
                         child + camper + persons,
                         dist = "negbin", data = ex.data)

# Use tab_model to display the models' summaries
tab_model(ex.m.zeroinf.nb, ex.m.hurdle.nb, ex.m.nb)

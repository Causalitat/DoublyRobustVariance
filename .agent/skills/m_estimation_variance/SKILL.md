---
name: M-Estimation Variance
description: Learn how to implement M-estimation for variance (sandwich estimator) using the geex package, specifically for causal inference estimators.
---

# M-Estimation Variance Skill

## 1. Overview
M-estimation (estimating equations) provides a unified framework for asymptotic variance estimation. It is particularly useful for constructing "sandwich" variance estimators that account for uncertainty in nuisance parameters (e.g., propensity scores, outcome models) in causal inference.

This skill guides you through implementing M-estimation using the `geex` R package.

## 2. Prerequisites
Ensure you have the following packages installed:
```r
library(geex)
library(numDeriv) # Often needed for derivatives
library(rootSolve) # Useful for finding roots if needed manually
```

## 3. The `geex` Workflow
The core of `geex` is the `m_estimate()` function. To use it, you need to define an `estFUN` (estimating function) that returns a function for calculating the score equations (psi functions) for a single observation.

### Structure of an `estFUN`
The `estFUN` takes `data` (and optionally `models`) as input and returns a function of `theta` (the parameters).

```r
my_estFUN <- function(data) {
    # Extract data variables for the current observation (or all if vectorized)
    # geex usually passes the full data frame to the outer function,
    # and expects the inner function to return a matrix of scores (n x p)
    # OR a list of functions if doing manual splitting.
    
    # Typically, for standard use, we define it to work on the full dataset:
    A <- data$A
    Y <- data$Y
    # ... extraction ...

    function(theta) {
        # theta is the vector of ALL parameters (nuisance + target)
        # Calculate scores for each parameter
        
        # ... logic ...
        
        # Return a matrix of scores: n rows (observations), p columns (parameters)
        return(cbind(score_param1, score_param2, ...)) 
    }
}
```

## 4. Common Score Equations (Templates)

### Logistic Regression (Propensity Score)
For a logistic model $logit(P(A=1|X)) = X\beta$:
$$ \psi(O_i, \beta) = (A_i - \text{plogis}(X_i^T \beta)) X_i $$

```r
# Helper to grab scores from a fitted glm object
# Usage in estFUN:
# e_scores <- grab_psiFUN(models$e, data)
# ... inside inner function ...
# current_score <- e_scores(theta[indices])
```

### Linear Regression (Outcome Model)
For a linear model $Y = X\beta + \epsilon$:
$$ \psi(O_i, \beta) = (Y_i - X_i^T \beta) X_i $$

## 5. Step-by-Step Implementation Guide

### Step 1: Fit Nuisance Models (Standard Way)
First, fit your nuisance models (e.g., propensity score, outcome model) using standard R functions like `glm()` or `lm()`. This gives you good starting values for the parameters.

```r
ps_model <- glm(A ~ W1 + W2, data = data, family = binomial)
out_model <- lm(Y ~ A + W1 + W2, data = data)
```

### Step 2: Define the Stacked Estimating Function
Create the function that stacks the scores for the nuisance parameters and your target parameter(s).

**Example: Inverse Probability Weighting (IPW)**
Target: $\mu_1 = E[Y^1]$
Estimating Equation: $\frac{A Y}{e(X)} - \mu_1 = 0$

```r
ipw_estFUN <- function(data, models){
  A <- data$A
  Y <- data$Y
  
  # Get design matrix and score function for propensity model
  Xe <- grab_design_matrix(data = data, rhs_formula = grab_fixed_formula(models$e))
  e_pos <- 1:ncol(Xe) # Indices for PS parameters
  e_scores <- grab_psiFUN(models$e, data)
  
  function(theta){
    # 1. Recover Parameters
    beta_ps <- theta[e_pos]
    mu_1 <- theta[max(e_pos) + 1] # Target parameter is last
    
    # 2. Calculate Nuisance Values
    e <- plogis(Xe %*% beta_ps)
    
    # 3. Calculate Scores
    # PS Scores
    s_ps <- e_scores(beta_ps)
    
    # IPW Score: (AY)/e - mu
    s_ipw <- (A * Y) / e - mu_1
    
    # 4. Return Combined Scores
    # geex expects a vector of scores for the current unit
    c(s_ps, s_ipw)
  }
}
```

### Step 3: Run `m_estimate`
Call `m_estimate` with the `estFUN`, the data, the roots (initial values), and any external models.

```r
# Initial values
# Nuisance params from fitted models
# Target param (mu_1) can be calculated simply first
init_mu1 <- mean(data$Y[data$A==1] / fitted(ps_model)[data$A==1]) 
roots <- c(coef(ps_model), init_mu1)

results <- m_estimate(
  estFUN = ipw_estFUN,
  data = data,
  roots = roots,
  outer_args = list(models = list(e = ps_model)),
  compute_roots = FALSE # We provide good roots, so no need to solve for them
)
```

### Step 4: Extract Variance
The variance-covariance matrix is stored in `results@vcov`. The standard error for your target parameter (usually the last one) is the square root of the corresponding diagonal element.

```r
n_params <- length(results@estimates)
variance <- results@vcov[n_params, n_params]
se <- sqrt(variance)
print(paste("Estimate:", results@estimates[n_params]))
print(paste("SE:", se))
```

## 6. Tips for Doubly Robust Estimators (AIPW/DR)
For AIPW, your estimating equation for $\mu_1$ involves both the propensity score $e(X)$ and the outcome model $m_1(X)$.

$$ \psi_{DR}(O_i) = \frac{A_i Y_i}{e_i} - \frac{A_i - e_i}{e_i} m_{1,i} - \mu_1 = 0 $$

You will need to:
1.  Pass both `ps_model` and `out_model` in `outer_args`.
2.  Extract design matrices and score functions for **both**.
3.  In the inner function, recover $\beta_{ps}$ and $\beta_{out}$ from `theta`.
4.  Compute $e_i$ and $m_{1,i}$.
5.  Compute the DR score.
6.  `c()` all scores: PS scores, Outcome scores, DR score.

## 7. Production Grade Patterns (Best Practices)
Based on high-quality statistical code (e.g., from `geex` papers), follow these patterns for robust implementation.

### 7.1. Modular Design
Separate the Estimating Function (`estfun`) from the Wrapper Function (`geex_wrapper`).

*   **`estfun_NAME`**: Pure function argument for `m_estimate`. Contains the math.
*   **`geex_NAME`**: User-facing function. Fits initial `glm` models, prepares roots, and calls `m_estimate`.

### 7.2. Explicit Design Matrices for Counterfactuals
When calculating potential outcomes ($m_1 = E[Y|A=1, W]$), do not rely on `predict()`. Instead:
1.  Use `grab_design_matrix` to get $X_m$ for the observed data.
2.  Create copies of data with $A=1$ (`data1`) and $A=0$ (`data0`).
3.  Use `grab_design_matrix` on these copies to get $X_{m1}$ and $X_{m0}$.
4.  Matrix-multiply inside the `estfun`: `m1 <- Xm1 %*% theta[m_pos]`.

**Example Pattern:**
```r
estfun_AIPW <- function(data, models){
  # 1. Setup Design Matrices (Observed and Counterfactual)
  Xe  <- grab_design_matrix(data = data, rhs_formula = grab_fixed_formula(models$e))
  Xm  <- grab_design_matrix(data = data, rhs_formula = grab_fixed_formula(models$m))
  
  # Counterfactuals
  data1 <- data; data1$A <- 1
  data0 <- data; data0$A <- 0
  Xm1 <- grab_design_matrix(data = data1, rhs_formula = grab_fixed_formula(models$m))
  Xm0 <- grab_design_matrix(data = data0, rhs_formula = grab_fixed_formula(models$m))
  
  # 2. Parameter Indices
  e_pos <- 1:ncol(Xe)
  m_pos <- (max(e_pos) + 1):(max(e_pos) + ncol(Xm))
  
  # 3. Score Functions
  e_scores <- grab_psiFUN(models$e, data)
  m_scores <- grab_psiFUN(models$m, data)
  
  function(theta){
    # Recover Parameters & Values
    e  <- plogis(Xe %*% theta[e_pos])
    m1 <- Xm1 %*% theta[m_pos] # E[Y|A=1]
    m0 <- Xm0 %*% theta[m_pos] # E[Y|A=0]
    
    # ... Score Calculation ...
  }
}
```

### 7.3. Handling Initial Values (Roots)
Always compute good initial values (roots) before calling `m_estimate`.
*   Nuisance parameters: Use `coef(model)`.
*   Target parameters: Use the plug-in estimate from your `glm` models (e.g., `mean(predict(m_model, data1))`).
*   Pass `compute_roots = FALSE` to `m_estimate` to save computation time and avoid convergence issues.

### 7.4. Complex Likelihoods (e.g., Zero-Inflated)
For complex models like ZINB (Zero-Inflated Negative Binomial), you cannot use standard `glm`. You must manually define the score function component.

See `noHeapingCount/00_Estimators_11.03.21.R` for examples of `ScoreZINB` functions that define the gradient of the log-likelihood manually.
```r
ScoreZINB <- function(A, L, ... parameters ...) {
   # Return vector of derivatives for each parameter
   c( score_param1, score_param2, ... )
}
```

### 7.5. Handling Many Regressors & Interactions (The Matrix Way)
Avoid manually writing out linear predictors like `theta[1] + theta[2]*A + theta[3]*L1 ...`. This is error-prone for large models or interactions.

Instead, rely on `grab_design_matrix` and matrix multiplication (`%*%`).

**Bad (Manual):**
```r
lam = exp(theta[1] + theta[2]*A + theta[3]*W1 + theta[4]*A*W1)
```

**Good (Automatic):**
1.  Define the model formula with interactions (e.g., `Y ~ A * W1 + .`).
2.  Extract the design matrix $X$.
3.  Compute $X \beta$ directly.

```r
# 1. Get Design Matrix (Handle A*W1, splines, etc. automatically)
X_lam <- grab_design_matrix(data = data, rhs_formula = grab_fixed_formula(models$zi_model))
lam_pos <- 1:ncol(X_lam)

function(theta) {
    # 2. Compute Linear Predictor
    lin_pred <- X_lam %*% theta[lam_pos]
    lam <- exp(lin_pred)
    
    # ...
}
```
```

## 8. Debugging Tips
If your estimator fails to converge or produces weird results, inspect the internal "bread" and "meat" matrices using `m_estimation_basis`.

```r
# 1. Create the basis object instead of running m_estimate
basis <- new("m_estimation_basis",
             .estFUN = my_estFUN,
             .data   = data)

# 2. Check roots (are they finding zeros?)
roots <- estimate_GFUN_roots(basis)
print(roots)

# 3. Check Derivatives (Bread) & Variance (Meat)
# Compute matrices at the found roots
mats <- estimate_sandwich_matrices(basis, .theta = roots$root)
bread <- grab_bread(mats)
meat  <- grab_meat(mats)

print(bread) # Should be invertible
print(meat)  # Should be positive semi-definite
```

## 9. Advanced Controls

### 9.1. Root Finding Control
If the default root finder fails, you can customize it or use `uniroot` for single-parameter problems.

```r
m_estimate(...,
  root_control = setup_root_control(start = c(0.5, 0.5), # Explicit start values
                                    method = "multiroot") # or "uniroot"
)
```

### 9.2. Bias Correction
`geex` supports bias corrections (e.g., Fay's correction for small samples).

```r
# Example from geex documentation
bias_correction <- function(components, b){ ... } # Define correction logic

results <- m_estimate(...,
  corrections = list(
    fay = correction(bias_correction, b = 0.1)
  )
)
```



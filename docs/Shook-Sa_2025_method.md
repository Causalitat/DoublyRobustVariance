# Double Robust Variance Estimation with Parametric Working Models

**Based on:** Shook-Sa, B. E., et al. (2025). *Double robust variance estimation with parametric working models*. Biometrics.

## 1. Overview and Estimand

**Objective:** Estimate the Average Causal Effect (ACE) of a binary exposure $X$ on an outcome $Y$ and provide doubly robust variance estimators.

**Notation:**
*   $X \in \{0, 1\}$: Binary exposure.
*   $Y$: Outcome.
*   $Z$: Vector of baseline covariates.
*   $Y^x$: Potential outcome under exposure $x$.
*   $O_i = (X_i, Y_i, Z_i)$: Observed data for subject $i$.

**Estimand:** Average Causal Effect (ACE)
$$ \text{ACE} = \mu^1 - \mu^0, \quad \text{where } \mu^x = E(Y^x) $$

**Key Concept:**
While doubly robust (DR) estimators for the point estimate (ACE) are well-known, the commonly used influence-function-based variance estimator is **not** doubly robust (it requires both outcome and propensity models to be correct). This paper highlights that the **empirical sandwich variance estimator** and the **nonparametric bootstrap** *are* doubly robust variance estimators.

---

## 2. Doubly Robust Point Estimators

Three DR estimators are considered. Each solves a set of estimating equations.

### 2.1. Classic AIPW Estimator ($\widehat{DR}_C$)
$$ \widehat{DR}_C = \hat{\mu}_C^1 - \hat{\mu}_C^0 $$
where
$$ \hat{\mu}_C^1 = \frac{1}{n} \sum_{i=1}^n \frac{X_i Y_i}{\hat{e}_i} - \frac{X_i - \hat{e}_i}{\hat{e}_i} \hat{a}_1(Z_i, \hat{\gamma}) $$
$$ \hat{\mu}_C^0 = \frac{1}{n} \sum_{i=1}^n \frac{(1-X_i) Y_i}{1-\hat{e}_i} + \frac{X_i - \hat{e}_i}{1-\hat{e}_i} \hat{a}_0(Z_i, \hat{\gamma}) $$
*   $\hat{e}_i = e(Z_i, \hat{\alpha})$: Propensity score (e.g., from logistic regression).
*   $\hat{a}_x(Z_i, \hat{\gamma})$: Predicted outcome $E(Y|X=x, Z)$ from an outcome model.

### 2.2. Weighted Regression AIPW Estimator ($\widehat{DR}_{WR}$)
$$ \widehat{DR}_{WR} = \hat{\mu}_{WR}^1 - \hat{\mu}_{WR}^0 $$
$$ \hat{\mu}_{WR}^x = \frac{1}{n} \sum_{i=1}^n \hat{b}_x(Z_i, \hat{\beta}) $$
*   $\hat{b}_x(Z_i, \hat{\beta})$: Predicted outcome from a weighted outcome model, where weights are the Inverse Probability of Treatment Weights (IPTW).

### 2.3. Targeted Maximum Likelihood Estimation ($\widehat{DR}_{TMLE}$)
$$ \widehat{DR}_{TMLE} = \hat{\mu}_{TMLE}^1 - \hat{\mu}_{TMLE}^0 $$
$$ \hat{\mu}_{TMLE}^x = \frac{1}{n} \sum_{i=1}^n \hat{c}_x(O_i, \hat{\gamma}, \hat{\eta}_x) $$
*   Uses a "targeting" step to update initial outcome predictions using the propensity scores.

---

## 3. Variance Estimation Methods

### 3.1. Influence Function Based Variance Estimator (Not Doubly Robust)
The standard variance estimator based on the efficient influence function:
$$ \hat{V}(\widehat{DR})_{IF} = \frac{1}{n^2} \sum_{i=1}^n \hat{I}_i^2 $$
*   **Limitation:** This estimator is **only consistent if both** the propensity score and outcome models are correctly specified. If one is misspecified, this variance estimator can be biased (conservative or anti-conservative), leading to incorrect confidence intervals.

### 3.2. Empirical Sandwich Variance Estimator (Doubly Robust)
**Mechanism:** Treat the DR estimator and the nuisance parameters (propensity and outcome model coefficients) as a joint system of M-estimation equations.
Let $\theta$ be the full vector of parameters (nuisance parameters + ACE). The estimator $\hat{\theta}$ solves $\sum \psi(O_i, \theta) = 0$.

**Variance Formula:**
$$ V(\theta) = A(\theta)^{-1} B(\theta) \{ A(\theta)^{-1} \}^T $$
where:
*   $A(\theta) = E\{ -\partial \psi(O_i; \theta) / \partial \theta \}$
*   $B(\theta) = E\{ \psi(O_i; \theta) \psi(O_i; \theta)^T \}$

**Implementation:**
This estimator is available in software packages like `geex` (R) or `delicatessen` (Python). It automatically accounts for uncertainty in the nuisance parameter estimation.
*   **Advantage:** It provides valid variance estimates (and thus valid CIs) if **either** the propensity or outcome model is correctly specified.

### 3.3. Nonparametric Bootstrap (Doubly Robust)
**Mechanism:**
1.  Draw $B$ resamples of size $n$ with replacement from the original data.
2.  Re-estimate the nuisance models (propensity, outcome) and the DR point estimate ($\widehat{DR}_b$) in each resample.
3.  Compute the variance of the $B$ bootstrap estimates.

**Variance Formula:**
$$ \hat{V}(\widehat{DR})_{NB} = \frac{1}{B-1} \sum_{b=1}^B (\widehat{DR}_b - \widehat{DR}^*)^2 $$

*   **Advantage:** Also doubly robust for variance.
*   **Disadvantage:** Computationally intensive as it requires refitting models $B$ times.

---

## 4. Summary of Steps for Practitioners

1.  **Choose a DR Point Estimator:** (e.g., AIPW, Weighted Regression AIPW, or TMLE).
2.  **Specify Parametric Working Models:**
    *   Propensity score model: $P(X=1|Z)$.
    *   Outcome model: $E(Y|X, Z)$.
3.  **Estimate Variance:**
    *   **Do NOT** rely solely on the standard influence-function-based SE (often the default in packages) if you suspect model misspecification.
    *   **Use Empirical Sandwich Variance:** Define the full stack of estimating equations (propensity scores + outcome regression + DR estimator) and use a sandwich estimator.
    *   **OR Use Nonparametric Bootstrap:** Resample data, refit *all* models, calculate point estimate, repeat 500+ times, and compute the standard deviation of the estimates.

## 5. Software
*   **R:** `geex` package for M-estimation (sandwich variance).
*   **Python:** `delicatessen` package for M-estimation.

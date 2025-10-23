#Program: FUN_v2.r
#Purpose: This program contains revised functions for DR estimators that can handle categorical covariates.

library(numDeriv)
library(boot)
library(geex)
library(resample)

# Helper function to sanitize formula for a given dataset
sanitize_formula_for_data <- function(formula, data) {
    ter <- terms(formula)
    vars <- all.vars(formula)

    response <- as.character(formula[[2]])
    predictors <- attr(ter, "term.labels")

    single_vars <- vars[vars %in% names(data) & vars != response]

    bad_vars <- c()
    for (v in single_vars) {
        is_bad <- FALSE
        if (all(is.na(data[[v]]))) {
            is_bad <- TRUE
        } else if (is.factor(data[[v]])) {
            if (nlevels(droplevels(data[[v]])) < 2) {
                is_bad <- TRUE
            }
        } else if (is.numeric(data[[v]]) || is.integer(data[[v]])) {
            if (var(data[[v]], na.rm = TRUE) == 0) {
                is_bad <- TRUE
            }
        } else if (is.character(data[[v]])) {
            if (length(unique(data[[v]][!is.na(data[[v]])])) < 2) {
                is_bad <- TRUE
            }
        }
        if (is_bad) {
            bad_vars <- c(bad_vars, v)
        }
    }

    if (length(bad_vars) == 0) {
        return(formula)
    }

    final_predictors <- predictors
    for (bad_v in bad_vars) {
        final_predictors <- final_predictors[!grepl(paste0("\\b", bad_v, "\\b"), final_predictors)]
    }

    if (length(final_predictors) == 0) {
        new_formula_str <- paste(response, "~ 1")
    } else {
        new_formula_str <- paste(response, "~", paste(final_predictors, collapse = " + "))
    }

    return(as.formula(new_formula_str))
}

##########################################################################################################
#### IF Variance Estimator ####
##########################################################################################################

IF_Var <- function(exposure,outcome,prop.score,Yhat0,Yhat1,est.DR){

  X<-exposure
  Y<-outcome
  ehat<-prop.score
  IF = (X*Y/ehat - (1-X)*Y/(1-ehat)) - ((X-ehat)/ehat/(1-ehat))*(((1-ehat)*Yhat1)+(ehat*Yhat0))-est.DR
  sdATE=sqrt(sum(IF^2)/(length(IF)^2))

  return(sdATE)
}


##########################################################################################################
#### M-Estimators ####
##########################################################################################################

#### Classic (Plug-in) AIPW Estimation ####

estfun_PI <- function(data,models){

  X<-data$X
  Y<-data$Y

  Xe <- grab_design_matrix(data = data, rhs_formula = grab_fixed_formula(models$e))
  Xm <- grab_design_matrix(data = data, rhs_formula = grab_fixed_formula(models$m))

  data0 <- data
  data0$X <- rep(0,nrow(data))
  Xm0 <- grab_design_matrix(data = data0,rhs_formula = grab_fixed_formula(models$m))

  data1 <- data
  data1$X <- rep(1,nrow(data))
  Xm1 <- grab_design_matrix(data = data1,rhs_formula = grab_fixed_formula(models$m))

  e_pos <- seq_len(ncol(Xe))
  last_e <- if (length(e_pos) > 0) max(e_pos) else 0
  m_pos <- seq_len(ncol(Xm)) + last_e

  e_scores <- grab_psiFUN(models$e, data)
  m_scores <- grab_psiFUN(models$m, data)

  function(theta){
    p<-length(theta)
    e <- plogis(Xe %*% theta[e_pos])
    m0 <- Xm0 %*% theta[m_pos]
    m1 <- Xm1 %*% theta[m_pos]

    CM1<-((X*Y - (X-e)*m1)/e)
    CM0<-(((1-X)*Y + (X-e)*m0)/(1-e))
    c(e_scores(theta[e_pos]),
      m_scores(theta[m_pos]),
      CM1-theta[p-2],
      CM0-theta[p-1],
      theta[p-2]- theta[p-1]-theta[p])
  }
}

geex_PI <- function(data, propensity_formula, outcome_formula, CM0.IV, CM1.IV, ATE.IV){

  safe_prop_formula <- sanitize_formula_for_data(propensity_formula, data)
  safe_out_formula <- sanitize_formula_for_data(outcome_formula, data)

  e_model  <- glm(safe_prop_formula, data = data, family =binomial)
  m_model  <- glm(safe_out_formula, data = data)
  models <- list(e = e_model, m=m_model)

  geex_resultsdrPI <-m_estimate(
    estFUN = estfun_PI,
    data   = data,
    roots = c(coef(e_model),coef(m_model), CM1.IV, CM0.IV, ATE.IV),
    compute_roots = FALSE,
    outer_args = list(models = models))

    DR.PI<-geex_resultsdrPI@estimates[length(geex_resultsdrPI@estimates)]
    seDR.PI <- sqrt(geex_resultsdrPI@vcov[length(geex_resultsdrPI@estimates),length(geex_resultsdrPI@estimates)])
    DR.PI.all<-cbind(DR.PI,seDR.PI)

  return(DR.PI.all)
}


#### Weighted Regression AIPW Estimation ####

gdm <- \(data, model) {
  grab_design_matrix(data = data, rhs_formula = grab_fixed_formula(model))
}

estfun_WTD <- \(data, models) {
  X <- data$X
  Xe <- gdm(data, models$e)
  Xm <- gdm(data, models$m)
  data0 <- data1 <- data
  data0$X <- 0
  data1$X <- 1
  Xm0 <- gdm(data0, models$m)
  Xm1 <- gdm(data1, models$m)
  e_scores <- grab_psiFUN(models$e, data)
  m_scores <- grab_psiFUN(models$m, data)

  e_pos <- seq_len(ncol(Xe))
  last_e <- if(length(e_pos) > 0) max(e_pos) else 0
  m_pos <- seq_len(ncol(Xm)) + last_e


  \(theta){

    p <- length(theta)
    e <- plogis(Xe %*% theta[e_pos])
    m0 <- Xm0 %*% theta[m_pos]
    m1 <- Xm1 %*% theta[m_pos]
    W  <- 1/dbinom(X, size = 1, prob = e)

    c(e_scores(theta[e_pos]),
      W * m_scores(theta[m_pos]),
      m1 - theta[p - 2],
      m0 - theta[p - 1],
      theta[p - 2] - theta[p - 1] - theta[p]
    )
  }
}

geex_WTD <- function(data, propensity_formula, outcome_formula,coef.WTD, CM0.IV, CM1.IV, ATE.IV){

  safe_prop_formula <- sanitize_formula_for_data(propensity_formula, data)
  safe_out_formula <- sanitize_formula_for_data(outcome_formula, data)

  e_model  <- glm(safe_prop_formula, data = data, family =binomial)
  m_model  <- glm(safe_out_formula, data = data)
  models <- list(e = e_model, m=m_model)

  geex_resultsdrWTD <-m_estimate(
    estFUN = estfun_WTD,
    data   = data,
    roots = c(coef(e_model),coef.WTD, CM1.IV, CM0.IV, ATE.IV),
    compute_roots = FALSE,
    outer_args = list(models = models))

  DR.WTD<-geex_resultsdrWTD@estimates[length(geex_resultsdrWTD@estimates)]
  seDR.WTD <- sqrt(geex_resultsdrWTD@vcov[length(geex_resultsdrWTD@estimates),length(geex_resultsdrWTD@estimates)])
  DR.WTD.all<-cbind(DR.WTD,seDR.WTD)

  return(DR.WTD.all)
}


#### TMLE ####

ScoreWTD_tar<- function(Y_sc,mu,W){
  c(W*(Y_sc-mu))
}

estfun_TMLE <- function(data,models){

  X<-data$X
  Y<-data$Y
  Y_sc<-data$Y_scaled
  a<-data$a
  b<-data$b

  Xe <- grab_design_matrix(data = data, rhs_formula = grab_fixed_formula(models$e))
  Xm <- grab_design_matrix(data = data, rhs_formula = grab_fixed_formula(models$m))
  data0 <- data
  data0$X <- rep(0,nrow(data))
  Xm0 <- grab_design_matrix(data = data0,rhs_formula = grab_fixed_formula(models$m))

  data1 <- data
  data1$X <- rep(1,nrow(data))
  Xm1 <- grab_design_matrix(data = data1,rhs_formula = grab_fixed_formula(models$m))

  e_pos <- seq_len(ncol(Xe))
  last_e <- if (length(e_pos) > 0) max(e_pos) else 0
  m_pos <- seq_len(ncol(Xm)) + last_e

  e_scores <- grab_psiFUN(models$e, data)
  m_scores <- grab_psiFUN(models$m, data)

  function(theta){
    p<-length(theta)
    e <- plogis(Xe %*% theta[e_pos])
    m0 <- Xm0 %*% theta[m_pos]
    m1 <- Xm1 %*% theta[m_pos]

    last_m <- if (length(m_pos) > 0) max(m_pos) else last_e

    mu0 <- inv.logit(theta[last_m+1]+logit(m0))
    mu1 <- inv.logit(theta[last_m+2]+logit(m1))
    w0<- (1-X)/(1-e)
    w1<- X/e

    CM0<-inv.logit(logit(m0)+theta[last_m+1])*(b-a)+a
    CM1<-inv.logit(logit(m1)+theta[last_m+2])*(b-a)+a

    c(e_scores(theta[e_pos]),
      m_scores(theta[m_pos]),
      ScoreWTD_tar(Y_sc,mu0,w0),
      ScoreWTD_tar(Y_sc,mu1,w1),
      CM1-theta[p-2],
      CM0-theta[p-1],
      theta[p-2]- theta[p-1]-theta[p])
  }
}

geex_TMLE <- function(data, propensity_formula, outcome_formula, tar0.IV, tar1.IV, CM0.IV, CM1.IV, ATE.IV){

  safe_prop_formula <- sanitize_formula_for_data(propensity_formula, data)
  safe_out_formula <- sanitize_formula_for_data(outcome_formula, data)

  e_model  <- glm(safe_prop_formula, data = data, family =binomial)
  m_model  <- glm(safe_out_formula, data = data)
  models <- list(e = e_model, m=m_model)

  geex_resultsdrTMLE <-m_estimate(
    estFUN = estfun_TMLE,
    data   = data,
    roots = c(coef(e_model),coef(m_model), tar0.IV, tar1.IV, CM1.IV, CM0.IV, ATE.IV),
    compute_roots = FALSE,
    outer_args = list(models = models))

  DR.TMLE<-geex_resultsdrTMLE@estimates[length(geex_resultsdrTMLE@estimates)]
  seDR.TMLE <- sqrt(geex_resultsdrTMLE@vcov[length(geex_resultsdrTMLE@estimates),length(geex_resultsdrTMLE@estimates)])
  DR.TMLE.all<-cbind(DR.TMLE,seDR.TMLE)

  return(DR.TMLE.all)
}


##########################################################################################################
#### Bootstrap Estimators ####
##########################################################################################################

### Plug-in AIPW

bootstrap_PIAIPW <- function(B, bootdata, propensity_formula, outcome_formula){
  if(B == 0) return(0)
  if(B>0){
    boot.est <- matrix(NaN, nrow = B, ncol = 1)
    datbi<-samp.bootstrap(nrow(bootdata), B)
    for(i in 1:B){
      dati <- bootdata[datbi[,i],]
      skip_to_next=FALSE

      safe_prop_formula <- sanitize_formula_for_data(propensity_formula, dati)
      safe_out_formula <- sanitize_formula_for_data(outcome_formula, dati)

      tryCatch(prop.model  <- glm(safe_prop_formula, data = dati, family = binomial), error = function(e) { skip_to_next <<- TRUE})
      if(skip_to_next==FALSE) {
        dati$PS<-predict(prop.model,dati, type="response")
        dati$IPTW<-ifelse(dati$X==1,dati$PS^(-1),(1-dati$PS)^(-1))}

      alltrt.boot<-alluntrt.boot<-dati
      alltrt.boot$X<-1
      alluntrt.boot$X<-0

      tryCatch(out.model.unwt  <- glm(safe_out_formula, data = dati), error = function(e) { skip_to_next <<- TRUE})
      if(skip_to_next==FALSE) {
        dati$Y0.unwt<-(predict(out.model.unwt, alluntrt.boot, type = "response"))
        dati$Y1.unwt<-(predict(out.model.unwt, alltrt.boot, type = "response"))

        CM1.PI<-mean(((dati$X*dati$Y-(dati$X-dati$PS)*dati$Y1.unwt))/dati$PS)
        CM0.PI<-mean((((1-dati$X)*dati$Y+(dati$X-dati$PS)*dati$Y0.unwt))/(1-dati$PS))
        boot.est [i]<-CM1.PI - CM0.PI}
      if(skip_to_next==TRUE) {
        boot.est [i]<-NA}
    }
    boot.se<-sd(boot.est,na.rm=TRUE)
    return(boot.se)
  }}


### Weighted Regression AIPW

bootstrap_WRAIPW <- function(B, bootdata, propensity_formula, outcome_formula){
  if(B == 0) return(0)
  if(B>0){
    boot.est <- matrix(NaN, nrow = B, ncol = 1)
    datbi<-samp.bootstrap(nrow(bootdata), B)
    for(i in 1:B){
      dati <- bootdata[datbi[,i],]
      skip_to_next=FALSE

      safe_prop_formula <- sanitize_formula_for_data(propensity_formula, dati)
      safe_out_formula <- sanitize_formula_for_data(outcome_formula, dati)

      tryCatch(prop.model  <- glm(safe_prop_formula, data = dati, family =binomial), error = function(e) { skip_to_next <<- TRUE})
      if(skip_to_next==FALSE) {
        dati$PS<-predict(prop.model,dati, type="response")
        dati$IPTW<-ifelse(dati$X==1,dati$PS^(-1),(1-dati$PS)^(-1))}

      alltrt.boot<-alluntrt.boot<-dati
      alltrt.boot$X<-1
      alluntrt.boot$X<-0

      tryCatch(out.model.wt  <- glm(safe_out_formula, data = dati,weights=IPTW), error = function(e) { skip_to_next <<- TRUE})
      if(skip_to_next==FALSE) {
        dati$Y0.wt<-(predict(out.model.wt, alluntrt.boot, type = "response"))
        dati$Y1.wt<-(predict(out.model.wt, alltrt.boot, type = "response"))

        CM0.WTD<-mean(dati$Y0.wt)
        CM1.WTD<-mean(dati$Y1.wt)
        boot.est [i]<-CM1.WTD - CM0.WTD}

      if(skip_to_next==TRUE) {
        boot.est [i]<-NA}
    }
    boot.se<-sd(boot.est,na.rm=TRUE)
    return(boot.se)
  }}


### TMLE

bootstrap_TMLE <- function(B, bootdata, propensity_formula, outcome_formula){
  if(B == 0) return(0)
  if(B>0){
    boot.est <- matrix(NaN, nrow = B, ncol = 1)
    datbi<-samp.bootstrap(nrow(bootdata), B)
    for(i in 1:B){
      dati <- bootdata[datbi[,i],]
      skip_to_next=FALSE

      safe_prop_formula <- sanitize_formula_for_data(propensity_formula, dati)
      safe_out_formula <- sanitize_formula_for_data(outcome_formula, dati)

      tryCatch(prop.model  <- glm(safe_prop_formula, data = dati, family =binomial), error = function(e) { skip_to_next <<- TRUE})
      if(skip_to_next==FALSE) {
        dati$PS<-predict(prop.model,dati, type="response")}

      alltrt.boot<-alluntrt.boot<-dati
      alltrt.boot$X<-1
      alluntrt.boot$X<-0

      a.boot<-min(dati$Y)
      b.boot<-max(dati$Y)
      dati$a.boot<-a.boot
      dati$b.boot<-b.boot
      dati$Y_scaled<-(dati$Y-a.boot)/(b.boot-a.boot)

      tryCatch(out.model.sc.unwt <- glm(safe_out_formula, data = dati), error = function(e) { skip_to_next <<- TRUE})
      if(skip_to_next==FALSE) {
        dati$Y0.sc.unwt<-(predict(out.model.sc.unwt, alluntrt.boot, type = "response"))
        dati$Y1.sc.unwt<-(predict(out.model.sc.unwt, alltrt.boot, type = "response"))

        dati$Y0.sc.unwt<-ifelse(dati$Y0.sc.unwt<min(dati$Y_scaled),min(dati$Y_scaled),dati$Y0.sc.unwt)
        dati$Y0.sc.unwt<-ifelse(dati$Y0.sc.unwt>max(dati$Y_scaled),max(dati$Y_scaled),dati$Y0.sc.unwt)

        dati$Y1.sc.unwt<-ifelse(dati$Y1.sc.unwt<min(dati$Y_scaled),min(dati$Y_scaled),dati$Y1.sc.unwt)
        dati$Y1.sc.unwt<-ifelse(dati$Y1.sc.unwt>max(dati$Y_scaled),max(dati$Y_scaled),dati$Y1.sc.unwt)}

      tryCatch(target.model.Y0.boot <- glm(Y_scaled ~ 1, offset=logit(Y0.sc.unwt), weights = (1-X)/(1-PS),
                                           family = binomial(), data = dati), error = function(e) { skip_to_next <<- TRUE})
      tryCatch(target.model.Y1.boot <- glm(Y_scaled ~ 1, offset=logit(Y1.sc.unwt), weights = X/PS,
                                           family = binomial(), data = dati), error = function(e) { skip_to_next <<- TRUE})
      if(skip_to_next==FALSE) {
        dati$Y0.TMLE<-inv.logit(logit(dati$Y0.sc.unwt)+coef(target.model.Y0.boot)[1])*(b.boot-a.boot)+a.boot
        dati$Y1.TMLE<-inv.logit(logit(dati$Y1.sc.unwt)+coef(target.model.Y1.boot)[1])*(b.boot-a.boot)+a.boot
        CM0.TMLE<-mean(dati$Y0.TMLE)
        CM1.TMLE<-mean(dati$Y1.TMLE)
        boot.est [i]<-CM1.TMLE-CM0.TMLE}

      if(skip_to_next==TRUE) {
        boot.est [i]<-NA}

    }
    boot.se<-sd(boot.est,na.rm=TRUE)
    return(boot.se)
  }}

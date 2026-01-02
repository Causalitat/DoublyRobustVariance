#Program: 01_Analysis.R
#Developed by: BES (Thanks to Bradley Saul for streamlining weighted regression AIPW functions)
#Purpose: Analyze example dataset ('exampledata.csv') using 3 DR estimators
#Created: last updated 04.24.24

#load required libraries
library(geex)
library(boot)
library(stats)

#specify number of bootstrap resamples
num.boots<-1000

#set working directory
setwd("/Users/user/Documents/Replicate_Gemini/DoublyRobustVariance/example/R")


#bring in example dataset: 
### outcome=Y, exposure=X, covariates=(Z1, Z2, Z3)
### Estimate the effect of X on Y, assuming covariates Z1, Z2, Z3 provide conditional exchangeability. 
### The correct functional forms of the propensity and outcome models are specified in Section 4 of the manuscript
dat<-read.csv('exampledata.csv')

# load functions from source file. 
source('00_Functions.R')

  
  ############################################################################################################################################
  ##### Compute each of the three DR estimators on observed data, with 3 variance estimates (influence function, sandwich, and bootstrap)#####
  ############################################################################################################################################
  
  #########################################################
  ##### fit weight model (needed for all estimators) #####
  #########################################################
  
  prop.model <- glm(X ~ Z1 + Z2 + Z1:Z2 + Z3 + Z1:Z3, family = binomial(), data = dat)
  dat$PS<-predict(prop.model,dat, type="response")
  dat$IPTW<-ifelse(dat$X==1,dat$PS^(-1),(1-dat$PS)^(-1))
  
  #create copies of data sets for pseudo-potential outcomes
  alltrt<-dat
  alluntrt<-dat
  alltrt$X<-1
  alluntrt$X<-0
  
 
  #########################
  #classic (plug-in) AIPW#
  #########################
  
  #fit outcome model (unweighted) and get pseudo-potential outcomes
  out.model.unwt <- glm(Y ~ Z1*Z2*X, data = dat)
  dat$Y0.unwt<-(predict(out.model.unwt, alluntrt, type = "response"))
  dat$Y1.unwt<-(predict(out.model.unwt, alltrt, type = "response"))
  
  #estimate causal means and ACE
  CM1.PI<-mean(((dat$X*dat$Y-(dat$X-dat$PS)*dat$Y1.unwt))/dat$PS)
  CM0.PI<-mean((((1-dat$X)*dat$Y+(dat$X-dat$PS)*dat$Y0.unwt))/(1-dat$PS))
  ATE.PI<-CM1.PI - CM0.PI
  
  #estimate variance using M-estimation, bootstrap, and the IF
  #M-est: input dataset (with exposure X, outcome Y), propensity model specification, outcome model specification, and initial values for each causal mean and the ACE
  DR.PI<-geex_PI(dat, X ~ Z1 + Z2 + Z1:Z2 + Z3 + Z1:Z3, Y ~ X*Z1*Z2, CM0.PI, CM1.PI, ATE.PI)

  #bootstrap: number of bootstraps, input dataset (with exposure X, outcome Y),propensity model specification, outcome model specification
  bootseATE.PI<-bootstrap_PIAIPW(num.boots, dat, X ~ Z1 + Z2 + Z1:Z2 + Z3 + Z1:Z3, Y ~ X*Z1*Z2)
  
  #IF: input the exposure and outcome variables, propensity score, pseudo-potential outcomes under no exposure and exposure, and the estimate ACE
  IFseATE.PI<-IF_Var(dat$X,dat$Y,dat$PS,dat$Y0.unwt,dat$Y1.unwt,ATE.PI)
  
  #format output
  DR.PI2<-as.data.frame(cbind(DR.PI,bootseATE.PI,IFseATE.PI))
  colnames(DR.PI2)<-c('ATE', 'ESseATE','bootseATE','IFseATE')
  DR.PI2$type<-'Classic'
  
  
  ##########################
  #weighted regression AIPW#
  ##########################
  
  #fit outcome model (weighted) and get pseudo-potential outcomes
  out.model.wt <- glm(Y ~ Z1*Z2*X, weights=IPTW, data = dat)
  dat$Y0.wt<-(predict(out.model.wt, alluntrt, type = "response"))
  dat$Y1.wt<-(predict(out.model.wt, alltrt, type = "response"))
  
  #estimate causal means and ACE
  CM0.WTD<-mean(dat$Y0.wt)
  CM1.WTD<-mean(dat$Y1.wt)
  ATE.WTD<-CM1.WTD - CM0.WTD
  
  #estimate variance using M-estimation and the IF
  #M-est: input dataset (with exposure X, outcome Y), propensity model specification, outcome model specification, and initial values for weighted outcome model, each causal mean, and the ACE
  DR.WTD<-geex_WTD(dat, X ~ Z1 + Z2 + Z1:Z2 + Z3 + Z1:Z3, Y ~ Z1*Z2*X, coef(out.model.wt), CM0.WTD, CM1.WTD, ATE.WTD)

  #bootstrap: number of bootstrap resamples, input dataset (with exposure X, outcome Y), propensity model specification, outcome model specification
  bootseATE.WTD<-bootstrap_WRAIPW(num.boots, dat, X ~ Z1 + Z2 + Z1:Z2 + Z3 + Z1:Z3, Y ~ Z1*Z2*X)
  
  #IF: input the exposure and outcome variables, propensity score, pseudo-potential outcomes under no exposure and exposure, and the estimate ACE
  IFseATE.WTD<-IF_Var(dat$X,dat$Y,dat$PS, dat$Y0.wt, dat$Y1.wt,ATE.WTD)
  
  #format output
  DR.WTD2<-as.data.frame(cbind(DR.WTD,bootseATE.WTD,IFseATE.WTD))
  colnames(DR.WTD2) <-c('ATE', 'ESseATE','bootseATE','IFseATE')
  DR.WTD2$type<-'WTD'
  

  ######
  #TMLE#
  ######
  
  #scale the outcome for TMLE
  a<-min(dat$Y)
  b<-max(dat$Y)
  dat$a<-a
  dat$b<-b
  dat$Y_scaled<-(dat$Y-a)/(b-a)

  #fit scaled outcome model (unweighted) and get stage 1 pseudo-potential outcomes
  out.model.sc.unwt <- glm(Y_scaled ~ Z1*Z2*X, data = dat)
  dat$Y0.sc.unwt<-(predict(out.model.sc.unwt, alluntrt, type = "response"))
  dat$Y1.sc.unwt<-(predict(out.model.sc.unwt, alltrt, type = "response"))
  
  #truncate predictions to within the range of the data
  dat$Y0.sc.unwt<-ifelse(dat$Y0.sc.unwt<min(dat$Y_scaled),min(dat$Y_scaled),dat$Y0.sc.unwt)
  dat$Y0.sc.unwt<-ifelse(dat$Y0.sc.unwt>max(dat$Y_scaled),max(dat$Y_scaled),dat$Y0.sc.unwt)
  
  dat$Y1.sc.unwt<-ifelse(dat$Y1.sc.unwt<min(dat$Y_scaled),min(dat$Y_scaled),dat$Y1.sc.unwt)
  dat$Y1.sc.unwt<-ifelse(dat$Y1.sc.unwt>max(dat$Y_scaled),max(dat$Y_scaled),dat$Y1.sc.unwt)
  
  #fit TMLE targeting models
  target.model.Y0 <- glm(Y_scaled ~ 1, offset=logit(Y0.sc.unwt), weights = (1-X)/(1-PS), 
                         family = binomial(), data = dat)
  target.model.Y1 <- glm(Y_scaled ~ 1, offset=logit(Y1.sc.unwt), weights = X/PS, 
                         family = binomial(), data = dat)
  
  #get stage 2 pseudo-potential outcomes and estimate the causal means and the ACE
  dat$Y0.TMLE<-inv.logit(logit(dat$Y0.sc.unwt)+coef(target.model.Y0)[1])*(b-a)+a
  dat$Y1.TMLE<-inv.logit(logit(dat$Y1.sc.unwt)+coef(target.model.Y1)[1])*(b-a)+a
  CM0.TMLE<-mean(dat$Y0.TMLE)
  CM1.TMLE<-mean(dat$Y1.TMLE)
  ATE.TMLE<-CM1.TMLE-CM0.TMLE
  
  #estimate variance using M-estimation and the IF
  #M-est: input dataset (with exposure X, outcome Y), propensity model specification, outcome model specification (for scaled outcome), and initial values for targeting models, each causal mean, and the ACE
  DR.TMLE<-geex_TMLE(dat, X ~ Z1 + Z2 + Z1:Z2 + Z3 + Z1:Z3, Y_scaled ~ Z1*Z2*X, coef(target.model.Y0), coef(target.model.Y1), CM0.TMLE, CM1.TMLE, ATE.TMLE)

  #bootstrap: number of bootstraps, input dataset (with exposure X, outcome Y), propensity model specification, outcome model specification (for scaled outcome)
  bootseATE.TMLE<-bootstrap_TMLE(num.boots, dat, X ~ Z1 + Z2 + Z1:Z2 + Z3 + Z1:Z3, Y_scaled ~ Z1*Z2*X)
  
  #IF: input the exposure and outcome variables, propensity score, stage 2 pseudo-potential outcomes under no exposure and exposure, and the estimate ACE
  IFseATE.TMLE<-IF_Var(dat$X,dat$Y,dat$PS, dat$Y0.TMLE, dat$Y1.TMLE,ATE.TMLE)
  
  #format output
  DR.TMLE2<-as.data.frame(cbind(DR.TMLE,bootseATE.TMLE,IFseATE.TMLE))
  colnames(DR.TMLE2)<-c('ATE', 'ESseATE','bootseATE','IFseATE')
  DR.TMLE2$type<-'TMLE'
  
  
 
  ######### combine all results and output #############
  DR.all<-rbind(DR.PI2,DR.WTD2,DR.TMLE2)
  colnames(DR.all)<-c('ATE', 'ESseATE','bootseATE','IFseATE')
  write.csv(DR.all, "example_estimates.csv")
  
rm(list=ls())

library(dplyr)
library(survival)
library(survMisc)
library(MASS)

data <- read.csv("../data/survDataCleaned.csv")
data.raceCleaned <- subset(data, select = -race_collapsed)

## re-level so the models use the following reference 
data.raceCleaned$race_cleaned <- relevel(factor(data.raceCleaned$race_cleaned), ref = "white")
data.raceCleaned$FIGO <- relevel(factor(data.raceCleaned$FIGO), ref = "Stage I")

## Full Models
fullCoxMod <- coxph(Surv(time, delta) ~ size.intermediate + factor(race_cleaned) + factor(FIGO)+age_at_diagnosis, data=data.raceCleaned, ties = "breslow")
full.expAFT <- survreg(Surv(time, delta) ~ size.intermediate + factor(race_cleaned) + 
                         factor(FIGO) + age_at_diagnosis, data = data.raceCleaned, dist = "exponential")
full.weibullAFT <- survreg(Surv(time, delta) ~ size.intermediate + factor(race_cleaned) + 
                             factor(FIGO) + age_at_diagnosis, data = data.raceCleaned, dist = "weibull")
full.loglogisticAFT <- survreg(Surv(time, delta) ~ size.intermediate + factor(race_cleaned) + 
                                 factor(FIGO) + age_at_diagnosis, data = data.raceCleaned, dist = "loglogistic")
full.lognormalAFT <- survreg(Surv(time, delta) ~ size.intermediate + factor(race_cleaned) + 
                               factor(FIGO) + age_at_diagnosis, data = data.raceCleaned, dist = "lognormal")

## Selected Models via AIC and BIC 
n <- nrow(data.raceCleaned) # number of rows, use this for BIC
CoxPH.aic <- stepAIC(fullCoxMod, direction="both", 
                     k=2, trace=1)
CoxPH.bic <- stepAIC(fullCoxMod, direction="both", 
                     k=log(n), trace=1)
expAFT.aic <- stepAIC(full.expAFT, direction = "both", k=2, trace=1)
expAFT.bic <- stepAIC(full.expAFT, direction = "both", k=log(n), trace=1)
weibullAFT.aic <- stepAIC(full.weibullAFT, direction = "both", k=2, trace=1)
weibullAFT.bic <- stepAIC(full.weibullAFT, direction = "both", k=log(n), trace=1)
loglogisticAFT.aic <- stepAIC(full.loglogisticAFT, direction = "both", k=2, trace=1)
loglogisticAFT.bic <- stepAIC(full.loglogisticAFT, direction = "both", k=log(n), trace=1)
lognormalAFT.aic <- stepAIC(full.lognormalAFT, direction = "both", k=2, trace=1)
lognormalAFT.bic <- stepAIC(full.lognormalAFT, direction = "both", k=log(n), trace=1)

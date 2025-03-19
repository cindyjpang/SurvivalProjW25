#########################################################
#### Check PH assumptions of full Cox model

setwd("~/Desktop/coursework/biostat215/SurvivalProjW25/r")

library(survival)
mypal <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728")
mypal2 <- c("#17BECF", "#BCBD22", "#FF6347", "#6A5ACD", "darkred")


data <- (read.csv("data/survDataCleaned.csv")
         %>% mutate(race_cleaned = factor(race_cleaned),
                    FIGO = factor(FIGO)))
data$race_cleaned <- relevel(data$race_cleaned, ref = "white")
data$FIGO <- relevel(data$FIGO, ref = "Stage I")

fullCoxMod <- coxph(Surv(time, delta) ~ size.intermediate + factor(race_cleaned) 
                    + factor(FIGO)+age_at_diagnosis, data=data, ties = "breslow")

summary(fullCoxMod)

############# First test continuous variables using Schoenfeld residuals
# Test the proportional hazards assumption
schoenfeld_test <- cox.zph(fullCoxMod)

par(mfrow=c(1,2))
plot(schoenfeld_test, var = "age_at_diagnosis", 
     main = "Schoenfeld plot for age at diagnosis")

plot(schoenfeld_test, var = "size.intermediate", 
     main = "Schoenfeld plot for biopsy size")

############# Andersen plot for factors
reptime <- function(l, t){
  x <- numeric(max(t))
  for(i in min(t):max(t)){
    diff <- i - t
    diff <- diff[diff >= 0]
    x[i] <- l[which.min(diff)]
  }
  return(x)
}

getH <- function(s, fit) {
  H <- fit$hazard[fit$strata == s]
  tt <- fit$time[fit$strata == s]
  return(reptime(H, tt))
} 

#### Race
coxMod_race <- coxph(Surv(time, delta) ~ size.intermediate + factor(FIGO)
                     + age_at_diagnosis + strata(race_cleaned), data=data, ties = "breslow")
levels(data$race_cleaned)
fit3 <- basehaz(coxMod_race)
Hbase <- getH("white", fit3)

par(mfrow=c(1,2))
plot(getH("asian", fit3)[1:length(Hbase)] ~ Hbase, 
     col = mypal2[1],
     main = 'Anderson plot for race', ylab = 'H(g)', xlab = 'H(white)', 
     type = 's', xlim = c(0, 0.35), ylim = c(0, 0.5))
lines(getH("black", fit3)[1:length(Hbase)] ~ Hbase, col = mypal2[2])
lines(getH("hispanic", fit3)[1:length(Hbase)] ~ Hbase, col = mypal2[3])
lines(getH("unreported/other", fit3)[1:length(Hbase)] ~ Hbase, col = mypal2[4])
abline(a = 0, b = 1)
legend("topleft", col = mypal2, lty = 1, lwd = 2, bty = "n", cex=1.2,
       legend = levels(data$race_cleaned)[-1])

### FIGO
coxMod_figo <- coxph(Surv(time, delta) ~ size.intermediate + factor(race_cleaned)
                     + age_at_diagnosis + strata(FIGO), data=data, ties = "breslow")
levels(data$FIGO)
fit2 <- basehaz(coxMod_figo)
Hbase <- getH("Stage I", fit2)

plot(getH("Stage II", fit2)[1:length(Hbase)] ~ Hbase, 
     col = mypal[1],
     main = 'Anderson plot for FIGO stage', ylab = 'H(g)', xlab = 'H(Stage I)', 
     type = 's', xlim = c(0, 0.2), ylim = c(0, 1))
lines(getH("Stage III", fit2)[1:length(Hbase)] ~ Hbase, col = mypal[2])
lines(pmin(1, getH("Stage IV", fit2)[1:length(Hbase)]) ~ Hbase, col = mypal[3])
abline(a = 0, b = 1)
legend("topleft", col = mypal, lty = 1, lwd = 2, bty = "n", cex=1.2,
       legend = levels(data$FIGO)[-1])



########## Plot log H vs t

## Race
coxMod_race <- coxph(Surv(time, delta) ~ size.intermediate + factor(FIGO)
                     + age_at_diagnosis + strata(race_cleaned), data=data, ties = "breslow")
levels(data$race_cleaned)
fit3 <- basehaz(coxMod_race)

par(mfrow=c(1,2))
plot(log(fit3$hazard[fit3$strata == "white"]) ~ fit3$time[fit3$strata == "white"], 
     col = mypal2[1], type = 'l', ylim = c(-5, -1),
     main = 'Log Cumulative Hazard, race',  ylab = 'log H(t)', xlab = 'time')
lines(log(fit3$hazard[fit3$strata == "asian"]) ~ fit3$time[fit3$strata == "asian"],
      col = mypal2[2])
lines(log(fit3$hazard[fit3$strata == "black"]) ~ fit3$time[fit3$strata == "black"],
      col = mypal2[3])
lines(log(fit3$hazard[fit3$strata == "hispanic"]) ~ fit3$time[fit3$strata == "hispanic"],
      col = mypal2[4])
lines(log(fit3$hazard[fit3$strata == "unreported/other"]) ~ fit3$time[fit3$strata == "unreported/other"],
      col = mypal2[5])
legend("bottomright", lty = 1, lwd = 2, col = mypal2, bty = "n",
       legend = levels(data$race_cleaned))


## FIGO
coxMod_figo <- coxph(Surv(time, delta) ~ size.intermediate + factor(race_cleaned)
                     + age_at_diagnosis + strata(FIGO), data=data, ties = "breslow")
levels(data$FIGO)
fit2 <- basehaz(coxMod_figo)

plot(log(fit2$hazard[fit2$strata == "Stage I"]) ~ fit2$time[fit2$strata == "Stage I"], 
     col = mypal[1], type = 'l', ylim = c(-5, 2),
     main = 'Log Cumulative Hazard, FIGO stage', ylab = 'log H(t)', xlab = 'time')
lines(log(fit2$hazard[fit2$strata == "Stage II"]) ~ fit2$time[fit2$strata == "Stage II"], 
      col = mypal[2])
lines(log(fit2$hazard[fit2$strata == "Stage III"]) ~ fit2$time[fit2$strata == "Stage III"], 
      col = mypal[3])
lines(log(fit2$hazard[fit2$strata == "Stage IV"]) ~ fit2$time[fit2$strata == "Stage IV"], 
      col = mypal[4])
legend("bottomright", lty = 1, lwd = 2, col = mypal, bty = "n",
       legend = levels(data$FIGO))



















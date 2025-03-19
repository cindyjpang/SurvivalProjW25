
setwd("~/Desktop/coursework/biostat215/SurvivalProjW25/r")
library(survMisc)
library(gtsummary)
library(gt)

mypal <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728")
mypal2 <- c("#17BECF", "#BCBD22", "#FF6347", "#6A5ACD", "darkred")

data <- (read.csv("data/survDataCleaned.csv")
         %>% mutate(time_to_death_years = time / 365,
                    race = factor(race_cleaned), 
                    FIGO = factor(FIGO)))
head(data)

###### Compare survival by FIGO stage
par(mfrow=c(1,2))
fit <- survfit(Surv(time/365, delta) ~ FIGO, data = data, conf.type = 'none')
plot(fit, col = mypal, xlab = "Time (years)", ylab = "S(t)", cex.lab = 1.2)
legend("bottomleft", col = mypal, lty = 1, lwd = 2, 
       bty = "n", cex = 1.2, legend = levels(data$FIGO))

logranks <- ten(Surv(time, delta) ~ FIGO, data = data)
comp(logranks, p = c(1, 0), q = c(0, 1))
tbl <- data.frame(attr(logranks, "tft")$tft)
tbl <- tbl[,c(1,2,6,8)]
colnames(tbl) <- c("Test", "Z(t)", "Chi-square", "p-value")
tbl[,-1] <- round(tbl[,-1], 2)
tbl<- cbind(" " = c("Log-rank", "Gehan", "Tarone-Ware", "Peto-Peto", 
                    "Modified PP","FH p=1,q=0", "FH p=0,q=1"), tbl[,-1])
gt(tbl) %>% tab_style(
  style = cell_text(weight = "bold"),  
  locations = cells_column_labels()) %>% tab_style(
    style = cell_text(weight = "bold"),  
    locations = cells_body(columns = 1))



###### Compare survival by race
fit <- survfit(Surv(time/365, delta) ~ race, data = data, conf.type = 'none')
plot(fit, col = mypal2, xlab = "Time (years)", ylab = "S(t)", cex.lab = 1.2)
legend("bottomleft", col = mypal2, lty = 1, lwd =2,  bty = "n", cex = 1.2,
       legend = levels(data$race))


logranks <- ten(Surv(time, delta) ~ race, data = data %>% filter(!grepl("unreported", race)))
comp(logranks, p = c(1, 0), q = c(0, 1))
tbl <- data.frame(attr(logranks, "tft")$tft)
tbl <- tbl[,c(1,2,6,8)]
colnames(tbl) <- c("Test", "Z(t)", "Chi-square", "p-value")
tbl[,-1] <- round(tbl[,-1], 2)
tbl<- cbind(" " = c("Log-rank", "Gehan", "Tarone-Ware", "Peto-Peto", 
              "Modified PP","FH p=1,q=0", "FH p=0,q=1"), tbl[,-1])
gt(tbl) %>% tab_style(
  style = cell_text(weight = "bold"),  
  locations = cells_column_labels()) %>% tab_style(
  style = cell_text(weight = "bold"),  
  locations = cells_body(columns = 1))



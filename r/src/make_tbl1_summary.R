
setwd("~/Desktop/coursework/biostat215/SurvivalProjW25/r")
library(gtsummary)


data <- (read.csv("data/survDataCleaned.csv")
         %>% mutate(time_to_death_years = time / 365,
                    Race = factor(race_cleaned), 
                    FIGO = factor(FIGO),
                    delta = ifelse(delta == 1, "dead", "alive"))
         %>% rename(`Survival time (years)` = time_to_death_years,
                    `Size of biopsy (cm)` =  size.intermediate,
                    `Age at diagnosis` = age_at_diagnosis,
                    `Vital status` = delta,
                    `FIGO stage` = FIGO))

head(data)


vars <- c("Survival time (years)", "Size of biopsy (cm)", "Age at diagnosis", "Race", "FIGO stage")
factorVars <- c( "Race", "FIGO stage")

table1 <- (data[,c(vars, "Vital status")] %>% tbl_summary(by = `Vital status`)
           %>% add_p()  # Add p-values
           %>% modify_header(label = "**Variable**") %>% bold_labels() 
           %>% add_overall())

table1





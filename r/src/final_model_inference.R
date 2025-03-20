
setwd("~/Desktop/coursework/biostat215/SurvivalProjW25/r")

data <- (read.csv("data/survDataCleaned.csv")
         %>% mutate(race_cleaned = factor(race_cleaned),
                    FIGO = factor(FIGO)))
data$race_cleaned <- relevel(data$race_cleaned, ref = "white")

data$FIGO_collapsed <- ifelse(data$FIGO %in% c("Stage I", "Stage II"), "Stage I/II", 
                              ifelse(data$FIGO == "Stage III", "Stage III", "Stage IV"))
data$FIGO_collapsed <- relevel(as.factor(data$FIGO_collapsed), ref = "Stage I/II")

data$race_collapsed <- ifelse(data$race_cleaned %in% c("white", "black", "asian"), 
                              as.character(data$race_cleaned), "unreported/other")
data$race_collapsed <- relevel(factor(data$race_collapsed), ref = "white")

data2 <- survSplit(Surv(time, delta) ~ size.intermediate + factor(race_collapsed) + factor(FIGO_collapsed) + age_at_diagnosis,
                   data = data, cut = c(opt_tau), episode = "tgroup")
data2$tgroup <- as.integer(data2$tgroup == 2)
colnames(data2) <- gsub("factor", "", colnames(data2))
colnames(data2) <- gsub("[()]", "", colnames(data2))

newCoxMod <- coxph(Surv(time, delta) ~ size.intermediate + factor(race_collapsed) 
                   + factor(FIGO_collapsed) + age_at_diagnosis +age_at_diagnosis:strata(tgroup),
                   data=data2, ties = "breslow") 
summary(newCoxMod)

cox_results <- tidy(newCoxMod, exponentiate = TRUE)  # Converts HR to exp(coef)
cox_results$conf.low <- cox_results$estimate - 1.96 * cox_results$std.error
cox_results$conf.high <- cox_results$estimate + 1.96 * cox_results$std.error
cox_results$term <-  c("Biopsy size cm", "Race/asian", "Race/black", "Race/unreported.other",
                       "Stage III", "Stage IV", "Age at diagnosis before 893", 
                       "Age at diagnosis after 893")
cox_results %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  mutate(
    HR = sprintf("%.2f (%.2f - %.2f)", estimate, conf.low, conf.high),
    p.value = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
  ) %>%
  select(term, HR, p.value) %>%
  gt() %>%
  cols_label(term = " ", HR = "Hazard Ratio (95% CI)", p.value = "p-Value") %>%
  fmt_missing(columns = everything(), missing_text = "N/A") %>%
  tab_options(table.font.size = px(14)) %>% tab_style(
    style = cell_text(weight = "bold"),  
    locations = cells_column_labels()) %>% tab_style(
      style = cell_text(weight = "bold"),  
      locations = cells_body(columns = 1))


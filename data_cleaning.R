library(dplyr)

#### data/cleaning
dat <- (read.csv("data/survProjData.csv") 
        %>% dplyr::select(sample, id, 
                          race = race.demographic, 
                          ethicity = ethnicity.demographic, 
                          vital_status = vital_status.demographic,
                          age_at_index = age_at_index.demographic, 
                          days_to_birth = days_to_birth.demographic,
                          year_of_birth = year_of_birth.demographic,
                          days_to_death = days_to_death.demographic,
                          year_of_death = year_of_death.demographic,
                          figo_stage.diagnoses, 
                          days_to_last_follow_up.diagnoses,
                          age_at_diagnosis = age_at_earliest_diagnosis_in_years.diagnoses.xena_derived,
                          primary_diagnosis.diagnoses,
                          size.shortest = shortest_dimension.samples, 
                          size.intermediate = intermediate_dimension.samples,
                          size.longest = longest_dimension.samples,
                          sample_type.samples))

# Convert to days since time to death is reported in days 
dat$age_at_index_days <- dat$age_at_index * 365 

######### Determine age at end of study. This will be the censored times for patients who are still alive
END_YEAR <- 2020
age_at_end <- END_YEAR - dat$year_of_birth 
days_in_study <- (age_at_end - dat$age_at_index) * 365
dat$days_to_death_censored <- ifelse(is.na(dat$days_to_death), days_in_study, dat$days_to_death)

######### Collapse FIGO stage categories into Stage 1-4
dat$FIGO <- gsub("[ABC]+$", "", dat$figo_stage.diagnoses)
# Four subjects have missing stages. Just randomly assign them to a stage
dat$FIGO[dat$FIGO == ""] <- paste("Stage", c("I", "II", "III", "IV"))
table(dat$FIGO)

######## Impute the 20 missing biopsy sizes with the mean
dat$size.intermediate[is.na(dat$size.intermediate)] <- mean(dat$size.intermediate, na.rm = TRUE)
dat$size.shortest[is.na(dat$size.shortest)] <- mean(dat$size.shortest, na.rm = TRUE)
dat$size.longest[is.na(dat$size.longest)] <- mean(dat$size.longest, na.rm = TRUE)

######## Collapse racial/ethnic categories 
dat <- dat %>% mutate(race_cleaned = case_when(ethicity == "hispanic or latino" ~ "hispanic",
                                                race == "white" ~ "white",
                                               race == "black or african american" ~ "black",
                                               race == "asian" ~ "asian", 
                                               TRUE ~ "unreported/other"))

# the counts are really imbalanced. we may need to collapse categories even further
dat$race_collapsed <- case_when(dat$race_cleaned == "white" ~ "white",
                                dat$race_cleaned %in% c("black", "hispanic") ~ "black/hisp.",
                                dat$race_cleaned == "asian" ~ "asian",
                                TRUE ~ "other/unreported")

table(dat$race_cleaned)
table(dat$race_collapsed)

####### Use clean_dat for analyses
dat$dead_indicator <- as.integer(dat$vital_status == "Dead")
clean_dat <- dat %>% dplyr::select(sample, 
                                   age_at_index, days_to_death, days_to_death_censored, dead_indicator,
                                   year_of_birth, year_of_death, age_at_diagnosis,
                                   size.shortest, size.intermediate, size.longest,
                                   FIGO, race_cleaned, race_collapsed)

summary(clean_dat)



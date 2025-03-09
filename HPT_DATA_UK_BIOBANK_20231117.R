#######  Osei Asibey - 2023.11.17  ######
###########################################################
# Getting Familiar with GP_Clinical Data in Shapter

# data loading and manipulation
install.packages("tidyverse")
library(tidyverse)
library(haven)
library(lukesRlib)
library(openxlsx)
library(vroom)
library(dplyr)
library(stringr)
library(nlme)
library(lme4)

# set default ggplot theme
theme_set(theme_bw())

# load packages for survival analysis and plotting
library(survival)
library(survivalROC)
library(ggsurvfit)
library(tidycmprsk)

# Load package Graphs
library(GGally)
library(ggplot2)
library(reshape2)

# load UK Biobank baseline data --------------------------------------------------------------------
# external script renames/formats common variables
#Setting Directory                      
setwd("I:/Projects/Osei")           

hpt_ukb_v = read.csv("H:/Projects/Biobank/Osei/Datafiles/hpt_ukb_v.csv")

hpt_ukb_comb = read.csv("H:/Projects/Biobank/Osei/Datafiles/hpt_ukb_v2v3.csv")

# Data Cleaning (Outlyers)
################################################################################
1# Rearranging Systolic and Diastolic BPs
hpt_ukb_v1 = hpt_ukb_v|>
  mutate(
    syst_bp = case_when(
      !is.na(systolic_bp) & is.na(diastolic_bp) ~ systolic_bp,
      is.na(systolic_bp) & !is.na(diastolic_bp) ~ NA,
      !is.na(systolic_bp) & !is.na(diastolic_bp) ~ if_else(systolic_bp > diastolic_bp, systolic_bp, diastolic_bp),
      TRUE ~ NA_real_),
    diast_bp = case_when(
      !is.na(systolic_bp) & is.na(diastolic_bp) ~ NA,
      is.na(systolic_bp) & !is.na(diastolic_bp) ~ diastolic_bp,
      !is.na(systolic_bp) & !is.na(diastolic_bp) ~ if_else(systolic_bp < diastolic_bp, systolic_bp, diastolic_bp),
      TRUE ~ NA_real_)
  ) %>%
  select(n_eid_14631, event_dt, READV2V3_CODE, BP_CATEGORIES.x, BP_CATEGORIES.y, BP_CATEGORIES, syst_bp, diast_bp)

hpt_ukb_v1 = hpt_ukb_v1|>
  select(n_eid_14631, event_dt, READV2V3_CODE, BP_CATEGORIES.x, BP_CATEGORIES.y, BP_CATEGORIES, syst_bp, diast_bp)


2# Setting Outlyers
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
1# Rename syst_bp and diast_bp to systolic_bp and Diastolic BP
hpt_ukb_v1 = hpt_ukb_v1|>
  rename("systolic_bp"="syst_bp", "diastolic_bp"= "diast_bp")


2# Rounding Values
hpt_ukb_v1$systolic_bp <- as.numeric(hpt_ukb_v1$systolic_bp)
hpt_ukb_v1$systolic_bp <- round(hpt_ukb_v1$systolic_bp)

hpt_ukb_v1$diastolic_bp <- as.numeric(hpt_ukb_v1$diastolic_bp)
hpt_ukb_v1$diastolic_bp <- round(hpt_ukb_v1$diastolic_bp)

3# Systolic Range 350mmHg to 24mmHg 
hpt_ukb_v1$systolic_bp[which(hpt_ukb_v1$systolic_bp >= 350)] <- NA
hpt_ukb_v1$systolic_bp[which(hpt_ukb_v1$systolic_bp <= 24)] <- NA

4# Diastolic Range 150 to 24mmHg
hpt_ukb_v1$diastolic_bp[which(hpt_ukb_v1$diastolic_bp >= 300)] <- NA
hpt_ukb_v1$diastolic_bp[which(hpt_ukb_v1$diastolic_bp <= 24)] <- NA

3# Data Cleaning, Removing all NAs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sorting out NA in both Systolic BP and Diastolic BP

#Remove rows with NAs in both Systolic_BP and Diastolic_BP
hpt_ukb_v2 <- hpt_ukb_v1|>
  filter(!is.na(systolic_bp) | !is.na(diastolic_bp))

# Check if the removed rows are only NAs for both columns
#2 to see which values are NA or missing values
hpt_ukb_v2Check <- hpt_ukb_v1|>
  filter(is.na(as.numeric(systolic_bp)) | is.na(as.numeric(diastolic_bp)))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
4# Remove Duplicates 
hpt_ukb_v3 <- hpt_ukb_v2 %>%
  distinct(n_eid_14631, event_dt, systolic_bp, diastolic_bp, .keep_all = TRUE)

hpt_ukb_v3 <- hpt_ukb_v2 %>%
  distinct(n_eid_14631, event_dt, .keep_all = TRUE)

# Check
duplicates <- hpt_ukb_v2[duplicated(hpt_ukb_v2[c("n_eid_14631", "event_dt", "systolic_bp", "diastolic_bp")]) | 
                           duplicated(hpt_ukb_v2[c("n_eid_14631", "event_dt", "systolic_bp", "diastolic_bp")], fromLast = TRUE), ]

duplicates2 <- hpt_ukb_v2[duplicated(hpt_ukb_v2[c("n_eid_14631", "event_dt")]) | 
                           duplicated(hpt_ukb_v2[c("n_eid_14631", "event_dt")], fromLast = TRUE), ]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Setting Pulse Pressure 

hpt_ukb_v4 <- hpt_ukb_v3 %>%
  mutate(Pulse_Pressure = systolic_bp - diastolic_bp)

# Filter if Pulse pressure below 20
hpt_ukb_v4b = hpt_ukb_v2|>
  filter(n_eid_14631=="3924070" & event_dt=="2008-05-14")

rm(hpt_ukb_V2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Check for duplicates
duplicates <- hpt_ukb_v2[duplicated(hpt_ukb_v2[c("n_eid_14631", "event_dt")]) | 
                           duplicated(hpt_ukb_v2[c("n_eid_14631", "event_dt")], fromLast = TRUE), ]

# Check for multiple entr with systolic>diastolic and NAs
hpt_ukb_v3mult1 <- hpt_ukb_v3 %>%
  group_by(n_eid_14631, event_dt) %>%
  filter(
    n() > 1 & # Checks if there are more than one row in the group.
      sum(!is.na(systolic_bp)) > 1 & #  Checks if there are more than one non-NA value in the systolic_bp
      any(is.na(diastolic_bp)))  # Checks if there is at least one NA value in the diastolic_bp

hpt_ukb_v3multi2 <- hpt_ukb_v3 %>%
  group_by(n_eid_14631, event_dt) %>%
  filter(
    n() > 1 &
      sum(!is.na(diastolic_bp)) > 1 &
      any(is.na(systolic_bp)))

rm(hpt_ukb_v3multi2)

# View the duplicate rows
View(duplicates)

hpt_ukb_comb = hpt_ukb_comb|>
  select(n_eid_14631, event_dt, READV2V3_CODE,READV2_DESC.x, BP_CATEGORIES, value1, value2)

hpt_ukb_v3 = hpt_ukb_v3|>
  mutate(pulse_pressure = systolic_bp-diastolic_bp)|>
  select(n_eid_14631, event_dt, READV2V3_CODE, READV2_DESC.x, BP_CATEGORIES.x, BP_CATEGORIES.y, BP_CATEGORIES, systolic_bp, diastolic_bp, pulse_pressure)

## Removing Multiple entry with NA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
hpt_ukb_v4 <- anti_join(hpt_ukb_v4, hpt_ukb_v3mult1, by = c("n_eid_14631", "event_dt"))

hpt_ukb_v4 <- anti_join(hpt_ukb_v4, hpt_ukb_v3multi2, by = c("n_eid_14631", "event_dt"))


# Save Files
file_path <- "I:/Projects/Osei/Datafiles/hpt_ukb_v3.csv"
write.csv(hpt_ukb_v3, file = file_path, row.names = FALSE)





######## Osei Asibey - 2023.11.17  #######
###########################################################
# Loading Baseline UK Biobank Data and Executing Script
# Shapter


## Loading External Script for Baseline UK Biobank Data ------------------------
souce("I:/Projects/Osei/_baseline_ukb.20240110.R")


## Joining baseline data........................................................
hpt_ukb = left_join(hpt_ukb_v3, ukb, by="n_eid_14631")


## Loading Health Outcomes......................................................

# Load recent death data
ukb_d = read_dta("H:/Projects/BioBank/14631_ageing-well/Death data/ukb14631_death_20230412.dta")
hpt_ukb = left_join(hpt_ukb, ukb_d, by="n_eid_14631")
rm(ukb_d)

# load HES up to Oct 2022
ukb_hes2022 = read_dta("H:/Projects/BioBank/14631_ageing-well/HES up to 2022 Oct/ukb14631_HES_20230929.dta")

## Dropping Extra Columns
# Function to drop extra column not needed now
drop_extra_ukb_cols2 = function(x)  {
  keepers = grepl("^n_|^hes2022_|^chd_|^cardiomyopathy_|^hf_|^mi_|^stroke_|^cerevd_|^pvd_|^t2d_hes_|^dementia_f01_|^dementia_hes_",colnames(x))
  keepers[1] = TRUE
  x = x[,keepers]
  x
}
# Drop extra col
ukb_hes2022 = drop_extra_ukb_cols2(ukb_hes2022)
hpt_ukb = left_join(hpt_ukb, ukb_hes2022, by="n_eid_14631")


# Other functions that did not work
 # ukb_hes2022 = ukb_hes2022 |> select(n_eid_14631, matches(names(ukb_hes2022),"hes2022|chd|cardiomyopathy|hf|mi|stroke|cerevd|pvd|t2d_hes|dmemntia_f01|dementia_hes"))
 # hpt_ukb = left_join(hpt_ukb, ukb_d, by="n_eid_14631")

 # ukb_hes2022 <- ukb_hes2022 %>% select(n_eid_14631, str_subset(names(ukb_hes2022), "hes2022|chd|cardiomyopathy|hf|mi|stroke|cerevd|pvd|t2d_hes|dementia_f01|dementia_hes"))

 # ukb_hes22 <- ukb_h %>%
  # select(n_eid_14631, grep("hes2022|chd|cardiomyopathy|hf|mi|stroke|cerevd|pvd|t2d_hes|dementia_f01|dementia_hes", names(ukb_hes2022), value = TRUE))

 # columns_to_keep <- c("hes2022", "chd", "cardiomyopathy", "hf", "mi", "stroke", "cerevd", "pvd", "t2d_hes", "dementia_f01", "dementia_hes")
 # ukb_hes22 <- ukb_hes22 %>%
  #  select(n_eid_14631, matches(paste0(columns_to_keep, collapse = "|")))

#
# Other variables to load later
# load genetic covariates
# load APOE genotype
# load Hypertension polygenic scores
#



#
# Data Analysis
#

# Variables available 
names(hpt_ukb)

# Getting the number of Unique Participants 
hpt_ukb_uniq = hpt_ukb |>
  distinct(n_eid_14631, .keep_all = TRUE)

# Number of participants
length(hpt_ukb_uniq$n_eid_14631)

length(unique((hpt_ukb$n_eid_14631)))

hpt_ukb |> select(n_eid_14631) |> pull() |> unique() |> length()
length(unique((hpt_ukb$n_eid_14631)))

# Tabulation of Variables
# Summary of participant age
hpt_ukb_uniq |>select(age) |> summary()

hpt_ukb_uniq |> select(sex) |> table()

hpt_ukb_uniq |> select(chd_hes2022) |> table()
table(hpt_ukb_uniq$chd_hes2022, useNA = "always")

hpt_ukb_uniq |> select(hf_hes2022) |> table()

hpt_ukb_uniq |> select(cardiomyopathy_hes2022) |> table()

hpt_ukb_uniq |> select(pvd_hes2022) |> table()

hpt_ukb_uniq |> select(cerevd_hes2022) |> table()

hpt_ukb_uniq |> select(stroke_ischemic_hes2022) |> table()

hpt_ukb_uniq |> select(dementia_f01_hes2022) |> table()

hpt_ukb_uniq |> select(bmi_cats) |> table()
table(hpt_ukb_uniq$bmi_catg, useNA = "always")

hpt_ukb_uniq |> select(education) |> table()

hpt_ukb_uniq |> select(ethnicity) |> table()

hpt_ukb_uniq |> select(smoking_status) |> table()

hpt_ukb |> select(stroke_ischemic_hes2022) |> table()
table(hpt_ukb$stroke_ischemic_hes2022, useNA = "always")



#
# Setting Data for Mixed Effect Analysis
#


1.# Sub-setting Data for Inclusion 

#~ Specify the Start Date of participants
start_date <- as.Date("2017-01-01")

  # Calculate age at start date by subtracting date of birth from the specific date
hpt_ukb <- hpt_ukb |>
  mutate(age_startdate = round(as.numeric(difftime(start_date, date_of_birth, units = "days")/365.25)))

# hpt_ukb_stroke3 <- hpt_ukb_stroke3 |>
  # mutate(age_startdate = round(as.numeric(difftime(start_date, date_of_birth, units = "days")/365.25)))

  # Data stroke (creating variable for those with stroke before and "after 70")
hpt_ukb_stroke <- hpt_ukb |>
  mutate(stroke_after70 = round(as.numeric(difftime(stroke_ischemic_hes2022_df, start_date, units = "days")/365.25)))

  
2. # Filter age by the age of recruitment
hpt_ukb_stroke = hpt_ukb_stroke |> 
  filter(age_startdate > 70,                      # Include Participant who were 70 at start Date
        #stroke_ischemic_hes2022 == 1 & !is.na(stroke_ischemic_hes2022), # Include Participant who Who have had stroke before or Missing variable
         stroke_after70 > 0)                      # Include Participant who were 70 at start Date

  # Those with event outcome of stroke
hpt_ukb_stroke1 = hpt_ukb_stroke |> 
  filter(age_startdate > 70,                      # Include Participant who were 70 at start Date
         stroke_ischemic_hes2022 == 1 & !is.na(stroke_ischemic_hes2022), # Include Participant who Who have had stroke before or Missing variable
         stroke_after70 > 0)  

3. # Filter for patients with a minimum of 2 systolic_bp measurements and event_dt within 10 years
hpt_ukb_stroke2 <- hpt_ukb_stroke %>%
  mutate(event_dt = as.Date(event_dt))

hpt_ukb_stroke2 <- hpt_ukb_stroke2 %>%
  group_by(n_eid_14631) %>%
  filter(
    n() >= 2,
    any(!is.na(event_dt)),
    all(event_dt >= as.Date("2007-01-01")),
    all(event_dt <= as.Date("2017-01-01"))
  ) %>%
  arrange(event_dt)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
hpt_ukb_stroke2_check <- hpt_ukb_stroke2 %>%
  group_by(n_eid_14631) %>%
  filter(
    !(
      n() >= 2 &
        any(!is.na(event_dt)) &
        all(event_dt >= as.Date("2007-01-01")) &
        all(event_dt <= as.Date("2017-01-01"))
    )
  ) %>%
  arrange(event_dt)
rm(hpt_ukb_stroke_check, hpt_ukb_stroke2_check, hpt_ukb_stroke2_2, hpt_ukb_stroke2_3)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Print unique values of event_dt for each group
filtered_data %>%
  group_by(n_eid_14631) %>%
  summarise(unique_dates = toString(unique(event_dt)))


4. # Calculate the average systolic_bp at the head and tail of event_dt
hpt_ukb_stroke_avg <- hpt_ukb_stroke2 %>%
  group_by(n_eid_14631) %>%
  summarize(
    avg_systolic_bp_head = mean(head(systolic_bp, n = 3, na.rm = TRUE)),
    avg_systolic_bp_tail = mean(tail(systolic_bp, n = 3, na.rm = TRUE))
  )

hpt_ukb_stroke_avg2 <- hpt_ukb_stroke2 %>%
  group_by(n_eid_14631) %>%
  mutate(
    avg_systolic_bp_head = mean(head(systolic_bp, n = 3, na.rm = TRUE)),
    avg_systolic_bp_tail = mean(tail(systolic_bp, n = 3, na.rm = TRUE))
  )

5. # Calculate Decline_SBP as the difference between the averages
hpt_ukb_stroke_avg <- hpt_ukb_stroke_avg %>%
  mutate(Decline_SBP = avg_systolic_bp_head - avg_systolic_bp_tail)

6.# Determining the Declining SBP Categories
hpt_ukb_stroke_avg = hpt_ukb_stroke_avg |> mutate(Decline_SBP_Cat = as_factor(case_when(
  Decline_SBP >= -5 & Decline_SBP < 5 ~ 1,
  Decline_SBP >= 5 ~ 2,
  Decline_SBP < -5 ~ 0)),
  Decline_SBP_Cat = labelled(Decline_SBP_Cat, c(Normal = 1, Rising_SBP = 2, Declining_SBP = 3), label="Decline_SBP (categories)")
)

hpt_ukb_stroke_avg |> select(Decline_SBP_Cat) |> table()
hpt_ukb_stroke_avg |> select(Decline_SBP_Cat) |> summary()




hpt_ukb_stroke3 = hpt_ukb_stroke_avg
rm(highest_systolic, hpt_ukb_stroke2_check, hpt_ukb_stroke2_1, less_than_six_measurements, lowest_systolic, more_than_six_measurements)

## Joining baseline data........................................................
hpt_ukb_stroke3 = left_join(hpt_ukb_stroke3, ukb, by="n_eid_14631")


## Loading Health Outcomes......................................................

# Load recent death data
ukb_d = read_dta("H:/Projects/BioBank/14631_ageing-well/Death data/ukb14631_death_20230412.dta")
hpt_ukb_stroke3 = left_join(hpt_ukb_stroke3, ukb_d, by="n_eid_14631")
rm(ukb_d)

# load HES up to Oct 2022
ukb_hes2022 = read_dta("H:/Projects/BioBank/14631_ageing-well/HES up to 2022 Oct/ukb14631_HES_20230929.dta")

## Dropping Extra Columns
# Function to drop extra column not needed now
drop_extra_ukb_cols2 = function(x)  {
  keepers = grepl("^n_|^hes2022_|^chd_|^cardiomyopathy_|^hf_|^mi_|^stroke_|^cerevd_|^pvd_|^t2d_hes_|^dementia_f01_|^dementia_hes_",colnames(x))
  keepers[1] = TRUE
  x = x[,keepers]
  x
}
# Drop extra col
ukb_hes2022 = drop_extra_ukb_cols2(ukb_hes2022)
hpt_ukb_stroke3 = left_join(hpt_ukb_stroke3, ukb_hes2022, by="n_eid_14631")




7. # Selecting relevant columns
hpt_ukb_stroke3 <- hpt_ukb_stroke3 |>
  select(
    n_eid_14631, 
    avg_systolic_bp_head, 
    avg_systolic_bp_tail, 
    Decline_SBP, 
    Decline_SBP_Cat, 
    stroke_ischemic_hes2022, 
    stroke_ischemic_hes2022_df, 
    stroke_ischemic_hes2022_prev, 
    stroke_ischemic_hes2022_prev, 
    education, 
    smoking_status, 
    weight,
    height,
    bmi,
    bmi_cats, 
    waist,
    hip,
    whr,
    cholesterol,
    hdl,
    ldl,
    hba1c,
    sex, 
    ethnicity, 
    date_of_birth, 
    date_of_death, 
    dead)

8. # Regression Analysis
 
# Linear Regression
model <- lm(stroke_ischemic_hes2022 ~ Decline_SBP, data = hpt_ukb_stroke3)
summary(model)
confint(model, level = 0.95)

summary(lm(stroke_ischemic_hes2022 ~ Decline_SBP_Cat, data = hpt_ukb_stroke3))
confint(lm(stroke_ischemic_hes2022 ~ Decline_SBP_Cat, data = hpt_ukb_stroke3), level = 0.95)

# age_startdate
summary(lm(stroke_ischemic_hes2022 ~ Decline_SBP + sex, data = hpt_ukb_stroke3))
confint(lm(stroke_ischemic_hes2022 ~ Decline_SBP + sex, data = hpt_ukb_stroke3), level = 0.95)

# sex
summary(lm(stroke_ischemic_hes2022 ~ Decline_SBP + sex + age_startdate, data = hpt_ukb_stroke3))
confint(lm(stroke_ischemic_hes2022 ~ Decline_SBP + sex + age_startdate, data = hpt_ukb_stroke3), level = 0.95)

# Education
summary(lm(stroke_ischemic_hes2022 ~ Decline_SBP + education, data = hpt_ukb_stroke3))
confint(lm(stroke_ischemic_hes2022 ~ Decline_SBP + education, data = hpt_ukb_stroke3), level = 0.95)

# Ethnicity
summary(lm(stroke_ischemic_hes2022 ~ Decline_SBP + ethnicity, data = hpt_ukb_stroke3))
confint(lm(stroke_ischemic_hes2022 ~ Decline_SBP + ethnicity, data = hpt_ukb_stroke3), level = 0.95)

# Smoking Status
summary(lm(stroke_ischemic_hes2022 ~ Decline_SBP + smoking_status, data = hpt_ukb_stroke3))
confint(lm(stroke_ischemic_hes2022 ~ Decline_SBP + smoking_status, data = hpt_ukb_stroke3), level = 0.95)

# bmi
summary(lm(stroke_ischemic_hes2022 ~ Decline_SBP + bmi, data = hpt_ukb_stroke3))
confint(lm(stroke_ischemic_hes2022 ~ Decline_SBP + bmi, data = hpt_ukb_stroke3), level = 0.95)

# bmi_cats
summary(lm(stroke_ischemic_hes2022 ~ Decline_SBP + bmi_cats, data = hpt_ukb_stroke3))
confint(lm(stroke_ischemic_hes2022 ~ Decline_SBP + bmi_cats, data = hpt_ukb_stroke3), level = 0.95)

# cholesterol
summary(lm(stroke_ischemic_hes2022 ~ Decline_SBP + cholesterol, data = hpt_ukb_stroke3))
confint(lm(stroke_ischemic_hes2022 ~ Decline_SBP + cholesterol, data = hpt_ukb_stroke3), level = 0.95)

# hdl
summary(lm(stroke_ischemic_hes2022 ~ Decline_SBP + hdl, data = hpt_ukb_stroke3))
confint(lm(stroke_ischemic_hes2022 ~ Decline_SBP + hdl, data = hpt_ukb_stroke3), level = 0.95)

# ldl
summary(lm(stroke_ischemic_hes2022 ~ Decline_SBP + ldl, data = hpt_ukb_stroke3))
confint(lm(stroke_ischemic_hes2022 ~ Decline_SBP + ldl, data = hpt_ukb_stroke3), level = 0.95)

# whr
summary(lm(stroke_ischemic_hes2022 ~ Decline_SBP + whr, data = hpt_ukb_stroke3))
confint(lm(stroke_ischemic_hes2022 ~ Decline_SBP + whr, data = hpt_ukb_stroke3), level = 0.95)

# dead
summary(lm(stroke_ischemic_hes2022 ~ Decline_SBP + dead, data = hpt_ukb_stroke3))
confint(lm(stroke_ischemic_hes2022 ~ Decline_SBP + dead, data = hpt_ukb_stroke3), level = 0.95)

# Dead
summary(lm(stroke_ischemic_hes2022 ~ Decline_SBP + dead, data = hpt_ukb_stroke3))
confint(lm(stroke_ischemic_hes2022 ~ Decline_SBP + dead, data = hpt_ukb_stroke3), level = 0.95)

# hba1c
summary(lm(stroke_ischemic_hes2022 ~ Decline_SBP + hba1c, data = hpt_ukb_stroke3))
confint(lm(stroke_ischemic_hes2022 ~ Decline_SBP + hba1c, data = hpt_ukb_stroke3), level = 0.95)


# Multivariate Regression
summary(lm(stroke_ischemic_hes2022 ~ Decline_SBP + cholesterol + ldl + whr + dead , data = hpt_ukb_stroke3))
confint(lm(stroke_ischemic_hes2022 ~ Decline_SBP + dead, data = hpt_ukb_stroke3), level = 0.95)

summary(lm(stroke_ischemic_hes2022 ~ Decline_SBP + sex + ethnicity + cholesterol + ldl + whr + dead, data = hpt_ukb_stroke3))
confint(lm(stroke_ischemic_hes2022 ~ Decline_SBP + sex + ethnicity + cholesterol + ldl + whr + dead, data = hpt_ukb_stroke3), level = 0.95)

model = lm(stroke_ischemic_hes2022 ~ Decline_SBP + sex + ethnicity + cholesterol + ldl + whr + dead, data = hpt_ukb_stroke3)


# Log Regression
model.log <- glm(stroke_ischemic_hes2022 ~ Decline_SBP, data = hpt_ukb_stroke3, family = "binomial")
summary(model.log)
confint(model.log, level = 0.95)
exp(cbind(OR=coef(model.log), confint(model.log, level = 0.95)))


model.log2 <- glm(stroke_ischemic_hes2022 ~ Decline_SBP_Cat, data = hpt_ukb_stroke3, family = "binomial")
summary(model.log2)
confint(model.log2, level = 0.95)
exp(cbind(OR=coef(model.log2), confint(model.log2, level = 0.95)))



# Plotting the coefficients
plot(model, which = 1)

# Plotting the residuals
plot(model, which = 2)

# Plotting the leverage
plot(model, which = 5)



  # Dealing with Missing values
hpt_isch_stroke1 <- hpt_isch_stroke |>
  subset(complete.cases(Decline_SBP,Time_since_Start, bmi_cats, cholesterol, sex, ethnicity), )
hpt_isch_stroke1$n_eid_14631 <- as.factor(hpt_isch_stroke1$n_eid_14631)
 
  # Model Formula
model <- lme(stroke_ischemic_hes2022 ~ Decline_SBP + Time_since_Start + Decline_SBP:Time_since_Start + sex + ethnicity + bmi_cats + cholesterol + (1|n_eid_14631), data = hpt_isch_stroke1)

model <- lme(stroke_ischemic_hes2022 ~ Decline_SBP + Time_since_Start + Decline_SBP:Time_since_Start + sex + ethnicity + bmi_cats + cholesterol, random = ~1|n_eid_14631, data = hpt_isch_stroke1)

lme.ridge(stroke_ischemic_hes2022 ~ Decline_SBP + Time_since_Start + Decline_SBP:Time_since_Start + sex + ethnicity + bmi_cats + cholesterol + (1|n_eid_14631), data = hpt_isch_stroke1)

summary(model)

install.packages("car")
library(car)
vif(model)
install.packages("MASS")
library(MASS)


cor(hpt_isch_stroke1[, c("Decline_SBP", "Time_since_Start", "bmi_cats", "cholesterol", "sex", "ethnicity")])





lme(stroke_ischemic_hes2022 ~ Decline_SBP + Time_since_Start + Decline_SBP:Time_since_Start + sex + ethnicity + bmi_cats + cholesterol + (1|n_eid_14631), data = hpt_uk_IS)

rm(hpt_isch_stroke2)
length(unique(hpt_isch_stroke$n_eid_14631))




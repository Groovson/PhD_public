## Osei Asibey - 2023.11.17  ##
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Shapter
  
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
install.packages("survival")
library(survival)
library(survivalROC)
library(ggsurvfit)
library(tidycmprsk)

# Load package Graphs
library(GGally)
library(ggplot2)
library(reshape2)


# loading Clean Data on Systolic BPs
# Loading Baseline UK Biobank Data and Executing Script

# load sytolic data (hpt_ukb_v3)
hpt_ukb_v3 = read.csv("I:/Projects/Osei/Datafiles/hpt_ukb_v3.csv")

# load baseline data (ukb)
hpt_ukb_v3 = read.csv("I:/Projects/Osei/Datafiles/hpt_ukb_v3.csv")


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


#
# Setting Data for Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


1.# Sub-setting Data for Inclusion (hpt_ukb_stroke)

#~ Specify the Study Start Date, Age at Start Date, and those without Stroke at startdate of participants
start_date <- as.Date("2017-01-01")

# Calculate age at start date by subtracting date of birth from the Study Start Date
hpt_ukb_stroke <- hpt_ukb |>
  mutate(
    age_startdate = round(as.numeric(difftime(start_date, date_of_birth, units = "days")/365.25)),
    stroke_after70 = round(as.numeric(difftime(stroke_ischemic_hes2022_df, start_date, units = "days")/365.25)))

# The hospital inpatient data (from Hospital Episode Statistics, HES) is summarized in 4 variables for each condition
# For example, for CHD:
#   chd_hes2022       binary variable, ever diagnosed (up to Oct 2022). [1]=yes, [0]=no
#   chd_hes2022_df    date variable, date of first diagnosis. If `chd_hes2022`==1 date of diagnosis... if `chd_hes2022`==0, date of HES censoring (31 Oct 2022)
#   chd_hes2022_prev  binary variable, diagnosed before UK Biobank baseline assessment. [1]=yes, [0]=no
#   chd_hes2022_exp   binary variable, diagnosed after UK Biobank baseline assessment (excludes prevalent). [1]=yes, [0]=no

2.#  Filter for inclusion Criteria
hpt_ukb_stroke = hpt_ukb_stroke |> 
  filter(age_startdate > 70,                      # Include Participant who are more than 70 at start Date
         stroke_after70 > 0)                      # Include paient stroke free at study start date

3.#  Filter for patients with a minimum of 2 systolic_bp measurements and event_dt within 10 years
hpt_ukb_stroke2 <- hpt_ukb_stroke |>
  mutate(event_dt = as.Date(event_dt) |>      # Convert to Date
  group_by(n_eid_14631) |>                    # Group by study id
  filter(
    n() >= 2,                                 # Include those with more than 2 recors of systolic bp
    any(!is.na(event_dt)),                    # Accoun for missing data
    all(event_dt >= as.Date("2007-01-01")),   # Include all patient from "2007-01-01" to "2017-01-01"
    all(event_dt <= as.Date("2017-01-01"))
  ) %>%
  arrange(event_dt)

hpt_ukb_stroke2 <- hpt_ukb_stroke %>%
  mutate(event_dt = as.Date(event_dt))

hpt_ukb_stroke2 <- hpt_ukb_stroke2 %>%
  group_by(n_eid_14631) %>%
  filter(
    n() >= 2,
    any(!is.na(event_dt)),
    all(event_dt >= as.Date("2007-01-01")),
    all(event_dt <= as.Date("2017-01-01"))) %>%
  arrange(event_dt)

4.#  Calculate the average systolic_bp at the head and tail of event_dt
  hpt_ukb_stroke_avg <- hpt_ukb_stroke2 |>
    group_by(n_eid_14631) |>
    summarize(
      avg_systolic_bp_head = mean(head(systolic_bp, n = 3, na.rm = TRUE)),
      avg_systolic_bp_tail = mean(tail(systolic_bp, n = 3, na.rm = TRUE)))
  
5.# Calculate Decline_SBP as the difference between the averages
hpt_ukb_stroke_avg <- hpt_ukb_stroke_avg |>
  mutate(Decline_SBP =avg_systolic_bp_tail - avg_systolic_bp_head)
  
6.# Determining the Declining SBP Categories
hpt_ukb_stroke_avg = hpt_ukb_stroke_avg |> mutate(Decline_SBP_Cat = as_factor(case_when(
  Decline_SBP >= -5 & Decline_SBP < 5 ~ 1,
  Decline_SBP >= 5 ~ 2,
  Decline_SBP < -5 ~ 3)),
  Decline_SBP_Cat = labelled(Decline_SBP_Cat, c(Normal = 1, Rising_SBP = 2, Declining_SBP = 3), label="Decline_SBP (categories)")
)

hpt_ukb_stroke_avg |> select(Decline_SBP_Cat) |> table()
hpt_ukb_stroke_avg |> select(Decline_SBP_Cat) |> summary()




hpt_ukb_stroke3 = hpt_ukb_stroke_avg

## Joining baseline data........................................................
hpt_ukb_stroke3 = left_join(hpt_ukb_stroke3, ukb, by="n_eid_14631")


## Loading Health Outcomes......................................................

# Load recent death data
ukb_d = read_dta("H:/Projects/BioBank/14631_ageing-well/Death data/ukb14631_death_20230412.dta")
hpt_ukb_stroke3 = left_join(hpt_ukb_stroke3, ukb_d, by="n_eid_14631")

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

# Calculate age at start date by subtracting date of birth from the Study Start Date
hpt_ukb_stroke3 <- hpt_ukb_stroke3 |>
  mutate(
    age_startdate = round(as.numeric(difftime(start_date, date_of_birth, units = "days")/365.25)),
    stroke_after70 = round(as.numeric(difftime(stroke_ischemic_hes2022_df, start_date, units = "days")/365.25)))


 # Selecting relevant columns
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
    age_startdate,
    stroke_after70,
    date_of_birth, 
    date_of_death, 
    dead)


8.# Analysis

#~ Regression Analysis

# Log Regression
model.log <- glm(stroke_ischemic_hes2022 ~ Decline_SBP_Cat, data = hpt_ukb_stroke3, family = "binomial")
summary(model.log)
confint(model.log, level = 0.95)
exp(cbind(OR=coef(model.log), confint(model.log, level = 0.95)))

# Adjust for age and sex
model.log2 <- glm(stroke_ischemic_hes2022 ~ Decline_SBP_Cat + sex + age_startdate, data = hpt_ukb_stroke3, family = "binomial")
summary(model.log2)
confint(model.log2, level = 0.95)
exp(cbind(OR=coef(model.log2), confint(model.log2, level = 0.95)))




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
library(purrr)

# set default ggplot theme
theme_set(theme_bw())

# load packages for survival analysis and plotting
library(survival)
library(survivalROC)
library(ggsurvfit)
library(tidycmprsk)

# Load package Graphs
install.packages("GGally")
library(GGally)
library(ggplot2)
library(reshape2)

# Package for Mixed Effect methods
install.packages("Matrix")  # Install or update the Matrix package
install.packages("lme4", dependencies = TRUE)  # Reinstall lme4 with dependencies
library(lme4)
install.packages("mixedup")
library(mixedup)
install.packages("merTools")
library(merTools)
library(dplyr)

# Set Directory
setwd("H:/Projects/Biobank/Osei")   


# loading Clean Data on Systolic BPs
# Loading Baseline UK Biobank Data and Executing Script

# load Data
hpt_ukb_sub = read.csv("H:/Projects/Biobank/Osei/Datafiles/hpt_ukb_sub1.csv")
#ukb <- read.csv("H:/Projects/Biobank/Osei/Datafiles/ukb.csv")
cmri_ukb <- read.csv("H:/Projects/Biobank/Osei/Datafiles/cmri_ukb1.csv")


hpt_ukb_sub <- left_join(hpt_ukb_sub , cmri_ukb, by = "n_eid_14631")




## Creating Parameters for inclusion and exclusion ##################################################

# Study Start date
hpt_ukb_sub <- hpt_ukb_sub |>
  arrange(n_eid_14631,desc(event_dt))|>
  group_by(n_eid_14631)|>
  mutate(study_start_date = max(event_dt))

hpt_ukb_sub$event_dt <- as.Date(hpt_ukb_sub$event_dt)
hpt_ukb_sub$study_start_date <- as.Date(hpt_ukb_sub$study_start_date)

# Age at study start date
hpt_ukb_sub <- hpt_ukb_sub |>
  mutate(age_SSD= round(as.numeric(difftime(study_start_date, date_of_birth, units = "days")/365.25)))

# Those with CMRI measure +/- 2yrs when their last SBP was measured
hpt_ukb_sub$date_of_cmri <- as.Date(hpt_ukb_sub$ts_53_2_0)
hpt_ukb_sub1 <- hpt_ukb_sub |>
  mutate(range_cmri= round(as.numeric(difftime(study_start_date, date_of_cmri, units = "days")/365.25)))

# Determining Year for each event_dt from the last SBP measurement
hpt_ukb_sub1 <- hpt_ukb_sub1 |> mutate(year = as.numeric(floor((study_start_date - event_dt)/365.25 ))) 
hpt_ukb_sub1 <- hpt_ukb_sub1 |> mutate(year2=year*-1)

# Determining the number of SBP measured per each participants
hpt_ukb_sub1 <- hpt_ukb_sub1|>
  arrange(n_eid_14631, desc(event_dt))|>
  group_by(n_eid_14631)|>
  mutate(n_count=n())

# You Can skip the next step for the general population without the subclassification of BSP

### Setting initial systolic BP ranges #####
hpt_ukb_sub2 <- hpt_ukb_sub1|>
  arrange(n_eid_14631,desc(event_dt))|>
  group_by(n_eid_14631)|>
  mutate(start_sbp = tail(systolic_bp, n=1, na.rm = TRUE),
         avg_start_sbp=mean(tail(systolic_bp, n=2, na.rm = TRUE)))

hpt_ukb_sub2 <- hpt_ukb_sub2|> 
  mutate(
    hpt_cat1=case_when(
      start_sbp < 140 ~ "NormalBP", 
      start_sbp >= 140 ~ "HighBP"),
    hpt_cat2=case_when(
      avg_start_sbp < 140 ~ "NormalBP",
      avg_start_sbp >= 140 ~ "HighBP"))


hpt_ukb_sub2 <- hpt_ukb_sub2 |> filter(hpt_cat1=="HighBP")
## Filters #######################################################################################
hpt_ukb_sub2 <- hpt_ukb_sub2 |> filter(!is.na(date_of_cmri))           # Those with no CMRI Data
hpt_ukb_sub2 <- hpt_ukb_sub2 |> filter(n_count>2)                      # Those with 3 or more data
hpt_ukb_sub2 <- hpt_ukb_sub2 |> filter(age_SSD>=70)                    # Those age 70+ when their last SBP was measured
hpt_ukb_sub2 <- hpt_ukb_sub2 |> filter(range_cmri>=-2 & range_cmri<=2) # Those with CMRI measure +/- 2yrs when their last SBP was measured
hpt_ukb_sub2 <- hpt_ukb_sub2 |> filter(year2>=-10)



## Extracting Coefficients for Slopes ############################################################

# Mixed Model
model_mx <- lmer(systolic_bp ~ year2 + (year2|n_eid_14631), hpt_ukb_sub2)    # Gener
summary(model_mx)


# Extract random effects
random_effects <- ranef(model_mx)$n_eid_14631 

#random_coef <- coef(model_mx)$n_eid_14631

# Convert to data frame
rand_df <- data.frame(n_eid_14631 = rownames(random_effects), 
                      coef_slope = random_effects$`year2`)

rand_df <- rand_df |> mutate(Decline_SBP5 = coef_slope*5, Decline_SBP7=coef_slope*7, Decline_SBP10 = coef_slope*10)

# Determining the Declining SBP Categories and converting it into categorical variables
rand_df <- rand_df |> mutate(Decline_SBP5_Cat = as.factor(case_when(
  Decline_SBP5 < -5 ~ "Decline",
  Decline_SBP5 >= -5 & Decline_SBP5 < 5 ~ "Stable",
  Decline_SBP5 >= 5 ~ "Increase")))

rand_df <- rand_df |> mutate(Decline_SBP7_Cat = as.factor(case_when(
  Decline_SBP7 < -5 ~ "Decline",
  Decline_SBP7 >= -5 & Decline_SBP7 < 5 ~ "Stable",
  Decline_SBP7 >= 5 ~ "Increase")))

rand_df <- rand_df |> mutate(Decline_SBP10_Cat = as.factor(case_when(
  Decline_SBP10 < -5 ~ "Decline",
  Decline_SBP10 >= -5 & Decline_SBP10 < 5 ~ "Stable",
  Decline_SBP10 >= 5 ~ "Increase")))

# Loading other covariates
rand_df$n_eid_14631 <- as.numeric(rand_df$n_eid_14631)
hpt_ukb_sub3 <- left_join(rand_df, hpt_ukb_sub1, by = "n_eid_14631")



## setting up for log regression ###################################################################
hpt_ukb_sub3$Decline_SBP5_Cat <- relevel(hpt_ukb_sub3$Decline_SBP5_Cat, ref = "Stable")
hpt_ukb_sub3$Decline_SBP7_Cat <- relevel(hpt_ukb_sub3$Decline_SBP7_Cat, ref = "Stable")
hpt_ukb_sub3$Decline_SBP10_Cat <- relevel(hpt_ukb_sub3$Decline_SBP10_Cat, ref = "Stable")

# Convert hpt_dx_hx to a factor
hpt_ukb_sub3$hpt_dx_hx <- as.factor(hpt_ukb_sub3$hpt_dx_hx)

# Relevel the factor to set '0' as the reference level
hpt_ukb_sub3$hpt_dx_hx <- relevel(hpt_ukb_sub3$hpt_dx_hx, ref = "0")




model_lm <- lm(ASI ~ Decline_SBP5_Cat +age_SSD + sex + hpt_dx_hx, data = hpt_ukb_sub3)
summary(model_lm)
confint(model_lm, level = 0.95)

model_lm <- lm(LVESV ~ Decline_SBP5_Cat +age_SSD + sex + education + hpt_dx_hx + bmi_cats + ldl_Cat + smoking_status + t2d_hes2022, data = hpt_ukb_sub3)
summary(model_lm)
confint(model_lm, level = 0.95)


# Obtain the coefficients
coefficients <- coef(summary(model_lm))

# Obtain the confidence intervals
conf_intervals <- confint(model_lm)

# Combine coefficients and confidence intervals
summary_with_ci <- cbind(coefficients, conf_intervals)

# Convert to data frame
(summary_df <- as.data.frame(summary_with_ci))

# Write to CSV file
write.csv(summary_df, "summary_with_ci.csv", row.names = FALSE)




file_path <- "H:/Projects/Biobank/Osei/Datafiles/rand_cmr.csv"
write.csv(rand_df, file = file_path, row.names = FALSE)

install.packages("tidyverse")
library(tidyverse)
library(haven)
library(lukesRlib)
require(openxlsx)
library(vroom)
library(dplyr)

# set default ggplot theme
theme_set(theme_bw())

# load packages for survival analysis and plotting
library(survival)
library(survivalROC)
library(ggsurvfit)
library(tidycmprsk)

# Package for Mixed Effect methods
install.packages("Matrix")  # Install or update the Matrix package
install.packages("lme4", dependencies = TRUE)  # Reinstall lme4 with dependencies
library(lme4)
library(dplyr)


install.packages("mixedup")
library(mixedup)
install.packages("merTools")
library(merTools)


setwd("H:/Projects/BioBank/Osei/")

## Load data 
# Baseline Data
# Systolic BPs
cmri_ukb = read.csv("H:/Projects/BioBank/Osei/Datafiles/cmri_ukb1.csv")

sbp_decline <- read.csv("H:/Projects/Biobank/Osei/Datafiles/sbp_decline.csv")

hpt_ukb_sub <- read.csv("H:/Projects/Biobank/Osei/Datafiles/hpt_ukb_sub1.csv")

## Select Cardiac MRI Parameters
cmri_ukb1 = cmri_ukb|>
  select(
    n_eid_14631,
    ts_53_2_0,        # Date of attending assessment centre (Instance 2)
    n_93_2_0,         # Systolic BP at assessment centre for CMRI (iNSTANCE 2)
    n_93_2_1,
    n_21003_2_0,      # Age at attending assessment centre (instance 2)
    n_22421_2_0,      # Cardiac MRI Chanber Funciotns
    n_22422_2_0,
    n_22423_2_0,
    n_22420_2_0,
    n_22424_2_0,
    n_4195_2_0,
    n_4196_2_0,
    n_4207_2_0,
    n_22421_3_0,
    n_22422_3_0,
    n_22423_3_0,
    n_22420_3_0,
    n_22424_3_0,
    n_4195_3_0,
    n_4196_3_0,
    n_4207_3_0,
    n_21021_0_0,
    n_21021_1_0,
    n_21021_2_0,
    n_21021_3_0,
    n_12681_2_0,     # Cardiac Stiffness (PWA Augmented Index)
    n_12681_2_1,
    n_12681_2_2,
    n_12681_2_3,
    n_12681_3_0,
    n_12681_3_1,
    n_12681_3_2,
    n_12681_3_3,
    n_12681_3_4)




## Cardiac MRI
# Arterial Stifness
cmri_ukb1$ASI <- cmri_ukb1$n_21021_2_0
cmri_ukb1$ASI <- ifelse(is.na(cmri_ukb1$ASI), cmri_ukb1$n_21021_3_0, cmri_ukb1$ASI)  

# Pulse Wave Analysis for Arteria Stiffness
cmri_ukb1$pwv_aix <- cmri_ukb1$n_12681_2_0
cmri_ukb1$pwv_aix <- ifelse(is.na(cmri_ukb1$pwv_aix), cmri_ukb1$n_12681_3_0, cmri_ukb1$n_12681_2_0)

# Left Ventricular End
cmri_ukb1$LVSV <- cmri_ukb1$n_22423_2_0
cmri_ukb1$LVSV <- ifelse(is.na(cmri_ukb1$LVSV), cmri_ukb1$n_22423_3_0, cmri_ukb1$LVSV) 

# LV_ejection_fraction=n_22420_3_0,
cmri_ukb1$LVEF <- cmri_ukb1$n_22420_2_0
cmri_ukb1$LVEF <- ifelse(is.na(cmri_ukb1$LVEF), cmri_ukb1$n_22420_3_0,cmri_ukb1$LVEF)

# LV_cardiac_output=n_22424_3_0
cmri_ukb1$LVCO <- cmri_ukb1$n_22424_2_0
cmri_ukb1$LVCO <- ifelse(is.na(cmri_ukb1$LVCO), cmri_ukb1$n_22424_3_0, cmri_ukb1$LVCO)

# LV_end_diastolic_volume=n_22421_3_0,
cmri_ukb1$LVEDV <- cmri_ukb1$n_22421_2_0
cmri_ukb1$LVEDV <- ifelse(is.na(cmri_ukb1$LVEDV), cmri_ukb1$n_22421_3_0, cmri_ukb1$LVEDV)

# LV_end_systolic_volume=n_22422_3_0,
cmri_ukb1$LVESV <- cmri_ukb1$n_22422_2_0
cmri_ukb1$LVESV <- ifelse(is.na(cmri_ukb1$LVESV), cmri_ukb1$n_22422_3_0, cmri_ukb1$LVESV)


# systolic_bp at CMRI
cmri_ukb1 <- cmri_ukb
hpt_ukb_cmri <- cmri_ukb1
hpt_ukb_cmri$systolic_bp <- hpt_ukb_cmri$n_93_2_0
hpt_ukb_cmri$systolic_bp <- ifelse(is.na(hpt_ukb_cmri$systolic_bp), hpt_ukb_cmri$n_93_2_1, hpt_ukb_cmri$systolic_bp)

hpt_ukb_cmri <-hpt_ukb_cmri|> filter(!is.na(systolic_bp))
hpt_ukb_cmri <- hpt_ukb_cmri|> select(n_eid_14631, ts_53_2_0, systolic_bp, n_21003_2_0)
hpt_ukb_cmri <- hpt_ukb_cmri|> mutate(study_start_date = ts_53_2_0) 
hpt_ukb_cmri <- hpt_ukb_cmri|> rename("event_dt" = "ts_53_2_0") 



# Join Cardiac MRI with Baseline Data
hpt_ukb_sub1  <- left_join(hpt_ukb_sub , cmri_ukb1, by = "n_eid_14631")
hpt_ukb_sub1 <-hpt_ukb_sub1 |> mutate(study_start_date = ts_53_2_0)             # Date of attending assessment centre (Instance 2)
hpt_ukb_sub1 <-hpt_ukb_sub1 |> filter(!is.na(study_start_date))
#hpt_ukb_sub1 <-hpt_ukb_sub1 |> select(n_eid_14631, event_dt, systolic_bp, n_21003_2_0, study_start_date)
hpt_ukb_sub1$event_dt <- as.Date(hpt_ukb_sub1$event_dt)
hpt_ukb_sub1$study_start_date <- as.Date(hpt_ukb_sub1$study_start_date)


## Version A hpt_ukb_Sub1a (Include SBP at CMRI date)
hpt_ukb_sub1a <-rbind(hpt_ukb_sub1, hpt_ukb_cmri)
hpt_ukb_sub1a <-hpt_ukb_sub1a |> mutate(year = as.numeric(floor((study_start_date-event_dt)/365.25 )))         # Floor runs into the number down to the nearest integer 
hpt_ukb_sub1a <- hpt_ukb_sub1a|> mutate(year2=year*-1)
hpt_ukb_sub1a$age_SSD <- hpt_ukb_sub1a$n_21003_2_0

hpt_ukb_sub1a <- hpt_ukb_sub1a|>
  arrange(n_eid_14631, desc(event_dt))|>
  group_by(n_eid_14631)|>
  mutate(n_count=n())


## Version B hpt_ukb_Sub1a (Include SBP at CMRI date)
hpt_ukb_sub1b <- hpt_ukb_sub1 |> mutate(year = as.numeric(floor((study_start_date-event_dt)/365.25 )))         # Floor runs into the number down to the nearest integer 
hpt_ukb_sub1b <- hpt_ukb_sub1b |> mutate(year2=year*-1)
hpt_ukb_sub1b$age_SSD <- hpt_ukb_sub1b$n_21003_2_0

hpt_ukb_sub1b <- hpt_ukb_sub1b|>
  arrange(n_eid_14631, desc(event_dt))|>
  group_by(n_eid_14631)|>
  mutate(n_count=n())


hpt_ukb_sub1c <- hpt_ukb_sub1|>
  arrange(n_eid_14631, desc(event_dt))|>
  group_by(n_eid_14631)|>
  mutate(n_count=n())




#hpt_ukb_sub1a <- hpt_ukb_sub1 |> select(n_eid_14631, event_dt, date_of_birth, study_start_date, year, year2, systolic_bp, sex, age_SSD, n_21003_2_0, n_count)

## Participant age 70+ with 10 years Trajectory before SSD
# Filters hpt_ukb_sub1a
hpt_ukb_sub1a <- hpt_ukb_sub1a|> filter(age_SSD>=70)                            # Accounting for age +70
hpt_ukb_sub1a <- hpt_ukb_sub1a|> filter(n_count>2)                              # Participants more than 1 SBP 
hpt_ukb_sub1a <- hpt_ukb_sub1a|> filter(year2>=-10 & year2<=0)   


# Filters hpt_ukb_sub1b
hpt_ukb_sub1b <- hpt_ukb_sub1b|> filter(age_SSD>=70)                            # Accounting for age +70
hpt_ukb_sub1b <- hpt_ukb_sub1b|> filter(n_count>2)                              # Participants more than 1 SBP 
hpt_ukb_sub1b <- hpt_ukb_sub1b|> filter(year2>=-10 & year2<=0)                  


# Mixed Model
summary(lmer(systolic_bp ~ year2 + (year2|n_eid_14631), hpt_ukb_sub1a))       

model_mxa <- lmer(systolic_bp ~ year2  + (year2|n_eid_14631), hpt_ukb_sub1a)        # 10 years trajectory
model_mxb <- lmer(systolic_bp ~ year2  +  (year2|n_eid_14631), hpt_ukb_sub1b)

#model_mx1 <- lmer(systolic_bp ~ year2 + sex + age_SSD + (year2|n_eid_14631), hpt_ukb_sub1a)  # Normal BP

summary(model_mxb)
 

## Extracting Coefficients for Slopes ###########################################
# Extract random effects
random_effects <- ranef(model_mxa)$n_eid_14631 
random_effects2 <- ranef(model_mxb)$n_eid_14631 


random_coef <- coef(model_mx)$n_eid_14631

# Convert to data frame
rand_dfa <- data.frame(n_eid_14631 = rownames(random_effects), 
                      coef_slope = random_effects$`year2`)

rand_dfb <- data.frame(n_eid_14631 = rownames(random_effects2), 
                       coef_slope = random_effects2$`year2`)

rand_dfa <- rand_dfa |> mutate(Decline_SBP5 = coef_slope*5, Decline_SBP10 = coef_slope*10)
rand_dfb <- rand_dfb |> mutate(Decline_SBP5 = coef_slope*5, Decline_SBP10 = coef_slope*10)

# Determining the Declining SBP Categories and converting it into categorical variables
rand_dfa <- rand_dfa |> mutate(Decline_SBP5_Cat = as_factor(case_when(
  Decline_SBP5 < -5 ~ "Decline",
  Decline_SBP5 >= -5 & Decline_SBP5 < 5 ~ "Stable",
  Decline_SBP5 >= 5 ~ "Increase")))

rand_dfa <- rand_dfa |> mutate(Decline_SBP10_Cat = as_factor(case_when(
  Decline_SBP10 < -5 ~ "Decline",
  Decline_SBP10 >= -5 & Decline_SBP10 < 5 ~ "Stable",
  Decline_SBP10 >= 5 ~ "Increase")))

rand_dfb <- rand_dfb |> mutate(Decline_SBP5_Cat = as_factor(case_when(
  Decline_SBP5 < -5 ~ "Decline",
  Decline_SBP5 >= -5 & Decline_SBP5 < 5 ~ "Stable",
  Decline_SBP5 >= 5 ~ "Increase")))

rand_dfb <- rand_dfb |> mutate(Decline_SBP10_Cat = as_factor(case_when(
  Decline_SBP10 < -5 ~ "Decline",
  Decline_SBP10 >= -5 & Decline_SBP10 < 5 ~ "Stable",
  Decline_SBP10 >= 5 ~ "Increase")))

# Loading other baseline variables
rand_dfa$n_eid_14631 <- as.numeric(rand_dfa$n_eid_14631)
hpt_ukb_sub2 <- left_join(rand_dfa, hpt_ukb_sub1 %>% distinct(n_eid_14631, .keep_all = TRUE), by = "n_eid_14631")

rand_dfb$n_eid_14631 <- as.numeric(rand_dfb$n_eid_14631)
hpt_ukb_sub3 <- left_join(rand_dfb, hpt_ukb_sub1 %>% distinct(n_eid_14631, .keep_all = TRUE), by = "n_eid_14631")


## setting up for log regression ###################################################

hpt_ukb_sub2$Decline_SBP5_Cat <- relevel(hpt_ukb_sub2$Decline_SBP5_Cat, ref = "Stable")
hpt_ukb_sub2$Decline_SBP10_Cat <- relevel(hpt_ukb_sub2$Decline_SBP10_Cat, ref = "Stable")

hpt_ukb_sub3$Decline_SBP5_Cat <- relevel(hpt_ukb_sub3$Decline_SBP5_Cat, ref = "Stable")
hpt_ukb_sub3$Decline_SBP10_Cat <- relevel(hpt_ukb_sub3$Decline_SBP10_Cat, ref = "Stable")

# Linear Regression model
model_lm <- lm(!is.na(LVCO) ~ Decline_SBP5_Cat, data = hpt_ukb_sub2a)
summary(lm(!is.na(LVCO) ~ Decline_SBP5_Cat + age_SSD + sex, data = hpt_ukb_sub2))
confint(lm(!is.na(LVCO) ~ Decline_SBP5_Cat + age_SSD + sex, data = hpt_ukb_sub2), level = 0.95)



model_lm <- lm(LVCO ~ Decline_SBP5_Cat + age_SSD + sex, data = hpt_ukb_sub2)

summary(model_lm)
confint(model_lm, level = 0.95)



## Extracting Output in excel formart
# Obtain the coefficients
coefficients <- coef(summary(model_lm))

# Obtain the confidence intervals
conf_intervals <- confint(model_lm)

# Combine coefficients and confidence intervals
summary_with_ci <- cbind(coefficients, conf_intervals)

# Convert to data frame
summary_df <- as.data.frame(summary_with_ci)

# Write to CSV file
write.csv(summary_df, "summary_with_ci.csv", row.names = FALSE)









## Saving File 
file_path <- "H:/Projects/Biobank/Osei/Datafiles/cmri_ukb1.csv"
write.csv(cmri_ukb1, file = file_path, row.names = FALSE)










hpt_ukb_sub1a <- hpt_ukb_sub1a|>
  arrange(n_eid_14631, desc(event_dt))|>
  group_by(n_eid_14631)|>
  mutate(study_SD = max(event_dt))

hpt_ukb_sub1a <- hpt_ukb_sub1a|>
  mutate(age_event = round(as.numeric(difftime(event_dt, date_of_birth, units = "days")/365.25)))

hpt_ukb_sub1b <- hpt_ukb_sub1a|> filter(year2==-10)


# Calculating Mean sd
(mean_age <- mean(hpt_ukb_sub1b$age_event, na.rm=TRUE))
[1] 63.4717
> (mean_age <- sd(hpt_ukb_sub1b$age_event, na.rm=TRUE))
[1] 2.780362
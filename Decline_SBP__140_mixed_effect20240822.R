#  Osei Asibey - 2024.04.23  
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Snow
  
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
library(dplyr)
install.packages("Matrix")
install.packages("Matrix")  # Install or update the Matrix package
install.packages("lme4", dependencies = TRUE)  # Reinstall lme4 with dependencies
library(lme4)
install.packages("mixedup")
library(mixedup)
install.packages("merTools")
library(merTools)


# Set Directory
setwd("H:/Projects/Biobank/Osei")   


# loading Clean Data on Systolic BPs
# Loading Baseline UK Biobank Data and Executing Script

# load Data
hpt_ukb_sub = read.csv("H:/Projects/Biobank/Osei/Datafiles/hpt_ukb_sub1.csv")
#ukb <- read.csv("H:/Projects/Biobank/Osei/Datafiles/ukb.csv")
cmri_ukb <- read.csv("H:/Projects/Biobank/Osei/Datafiles/cmri_ukb1.csv")



## Setting Data for Mixed Effect Model###########################################
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

# Filter Age at 70years
hpt_ukb_mx <- hpt_ukb_sub|>
  filter(age_SSD>=70) 

# Setting initial systolic BP ranges
hpt_ukb_mx <-hpt_ukb_mx|>
  arrange(n_eid_14631,desc(event_dt))|>
  group_by(n_eid_14631)|>
  mutate(start_sbp = tail(systolic_bp, n=1, na.rm = TRUE),
         avg_start_sbp=mean(tail(systolic_bp, n=2, na.rm = TRUE)))

hpt_ukb_mx<- hpt_ukb_mx|> 
  mutate(
    hpt_cat1=case_when(
      start_sbp < 140 ~ "NormalBP", 
      start_sbp >= 140 ~ "HighBP"),
    hpt_cat2=case_when(
      avg_start_sbp < 140 ~ "NormalBP",
      avg_start_sbp >= 140 ~ "HighBP"))




hpt_ukb_mx$event_dt<-as.Date(hpt_ukb_mx$event_dt)
hpt_ukb_mx$study_start_date <- as.Date(hpt_ukb_mx$study_start_date)

hpt_ukb_mx <-hpt_ukb_mx |> mutate(year = as.numeric(floor((study_start_date-event_dt)/365.25 ))) # Floor runs into the number down to the nearest integer 

#hpt_ukb_mx <-hpt_ukb_mx |> select(n_eid_14631, event_dt, study_start_date, year, systolic_bp, age_SSD, sex, hpt_cat1, hpt_cat2)

hpt_ukb_mx <-hpt_ukb_mx |> mutate(year2=year*-1)
hpt_ukb_mx <-hpt_ukb_mx |> filter(year2>=-10)


hpt_ukb_mx1 <- hpt_ukb_mx|> filter(hpt_cat1=="HighBP")
hpt_ukb_mx2 <- hpt_ukb_mx|> filter(hpt_cat2=="HighBP")


## Regressional Models ##########################################################
# Linear Regression
stroke_lm = lm(systolic_bp ~ year_vs + age_SSD + sex, data = hpt_ukb_stroke2mx1)
summary(stroke_lm)
confint(stroke_lm, level = 0.95)
exp(cbind(OR=coef(stroke_lm), confint(stroke_lm, level = 0.95)))

# Mixed Model

#lm_mixed2 = lmer(systolic_bp ~ year2 + (1|n_eid_14631), data = hpt_ukb_stroke2mx2)
#lm_mixed = lmer(systolic_bp ~ year2 + age_SSD + sex + (year2 | n_eid_14631), data = hpt_ukb_stroke2mx2)
model_mx <- lmer(systolic_bp ~ year2 + (year2|n_eid_14631), hpt_ukb_mx)          # General

model_mx1 <- lmer(systolic_bp ~ year2 + (year2|n_eid_14631), hpt_ukb_mx1)   
model_mx2 <- lmer(systolic_bp ~ year2 + (year2|n_eid_14631), hpt_ukb_mx2)        # High BP

summary(model_mx2)
summary(model_mx)
summary(model_mx1)


## Extracting Coefficients for Slopes ###########################################
# Extract random effects
random_effmx <- ranef(model_mx1)$n_eid_14631 

random_coef <- coef(model_mx)$n_eid_14631

# Convert to data frame
rand_dfmx <- data.frame(n_eid_14631 = rownames(random_effmx), 
                        coef_slope = random_effmx$`year2`)

rand_dfmx <- rand_dfmx |> mutate(Decline_SBP5 = coef_slope*5, Decline_SBP7=coef_slope*7, Decline_SBP10 = coef_slope*10)


# Determining the Declining SBP Categories and converting it into categorical variables
rand_dfmx <- rand_dfmx |> mutate(Decline_SBP5_Cat = as.factor(case_when(
  Decline_SBP5 < -5 ~ "Decline",
  Decline_SBP5 >= -5 & Decline_SBP5 < 5 ~ "Stable",
  Decline_SBP5 >= 5 ~ "Increase")))

rand_dfmx <- rand_dfmx |> mutate(Decline_SBP7_Cat = as.factor(case_when(
  Decline_SBP7 < -5 ~ "Decline",
  Decline_SBP7 >= -5 & Decline_SBP7 < 5 ~ "Stable",
  Decline_SBP7 >= 5 ~ "Increase")))

rand_dfmx <- rand_dfmx |> mutate(Decline_SBP10_Cat = as.factor(case_when(
  Decline_SBP10 < -5 ~ "Decline",
  Decline_SBP10 >= -5 & Decline_SBP10 < 5 ~ "Stable",
  Decline_SBP10 >= 5 ~ "Increase")))

# Loading other covariates
rand_dfmx$n_eid_14631 <- as.numeric(rand_dfmx$n_eid_14631)
hpt_ukb_submx <- left_join(rand_dfmx, hpt_ukb_sub %>% distinct(n_eid_14631, .keep_all = TRUE), by = "n_eid_14631")




## setting up for log regression ####################################################
hpt_ukb_submx$Decline_SBP5_Cat <- relevel(hpt_ukb_submx$Decline_SBP5_Cat, ref = "Stable")
hpt_ukb_submx$Decline_SBP7_Cat <- relevel(hpt_ukb_submx$Decline_SBP7_Cat, ref = "Stable")
hpt_ukb_submx$Decline_SBP10_Cat <- relevel(hpt_ukb_submx$Decline_SBP10_Cat, ref = "Stable")


# Logistic Regression
model.log <- glm(stroke_ischemic_hes2022 ~ Decline_SBP5_Cat, data = hpt_ukb_submx, family = "binomial")
summary(model.log)
confint(model.log, level = 0.95)
exp(cbind(OR=coef(model.log), confint(model.log, level = 0.95)))

model.log1 <- glm(mi_hes2022 ~ Decline_SBP10_Cat, data = hpt_ukb_submx, family = "binomial")
summary(model.log1)
confint(model.log1, level = 0.95)
exp(cbind(OR=coef(model.log1), confint(model.log1, level = 0.95)))

# Fitting a Survival Analysis
# Conditions
survival_obj <- Surv(time = as.numeric(difftime(hpt_ukb_submx$stroke_ischemic_hes2022_df, hpt_ukb_submx$event_dt, units = "days")), 
                     event = hpt_ukb_submx$stroke_ischemic_hes2022)

# Fit a Cox Proportional Hazards model
cox_model <- coxph(survival_obj ~ Decline_SBP5_Cat + age_SSD + sex, data = hpt_ukb_submx)
summary(cox_model)



## All Cause of Mortality###################################################################################
hpt_ukb_submx$mort_ed <- as.Date("2022-10-31")
hpt_ukb_submx$dod <- hpt_ukb_submx$date_of_death
hpt_ukb_submx$dod <- as.Date(hpt_ukb_submx$dod)
hpt_ukb_submx$dod <- if_else(is.na(hpt_ukb_submx$dod), hpt_ukb_submx$mort_ed, hpt_ukb_submx$dod)



# Create the survival object
surv_obj <- Surv(time = as.numeric(difftime(hpt_ukb_submx$dod, hpt_ukb_submx$event_dt, units = "days")), 
                 event = hpt_ukb_submx$dead)


surv_obj <- Surv(time = as.numeric(difftime(hpt_ukb_submx$date_of_death, hpt_ukb_submx$event_dt, units = "days")), 
                 event = hpt_ukb_submx$dead)


# Fit the Cox proportional hazards model
cox_model <- coxph(surv_obj ~ Decline_SBP10_Cat + age_SSD + sex, data = hpt_ukb_submx)

# Print the summary of the model
summary(cox_model)


any(is.na(hpt_ukb_mort$dod))
any(is.na(hpt_ukb_mort$event_dt))

class(hpt_ukb_mort$event_dt)
class(hpt_ukb_mort$dod)

length(hpt_ukb_mort$dod)
length(hpt_ukb_mort$event_dt)









## Plot Coefficience
plot((model.log)$coef[-1, 1], col = ifelse(significant, "red", "black"),  pch = 19)



## Save File
file_path <- "H:/Projects/Biobank/Osei/Datafiles/hpt_ukb_sub1.csv"
write.csv(hpt_ukb_sub1, file = file_path, row.names = FALSE)


# Predict Model

## Predict how low systolic BP will be in 5years
new_data <- data.frame(year2 = 5, sex=0, age_SSD=0)

predict(stroke_mixed2, newdata = new_data, re.form = NA)


predictInterval(stroke_mixed2)   # for various model predictions, possibly with new data

REsim(gpa_mixed)             # mean, median and sd of the random effect estimates

plotREsim(REsim(stroke_mixed2))  # plot the interval estimates



## Plot

# Plot the data with a best fit line
library(merTools)

# Generate prediction intervals
pred_intervals <- predictInterval(model_mx, newdata = hpt_ukb_sub2, level = 0.9)

# Combine the original data with the intervals
hpt_ukb_sub3 <- cbind(hpt_ukb_sub2, pred_intervals)

# Plot with prediction intervals
ggplot(hpt_ukb_sub3, aes(x = year2, y = fit)) +
  geom_line(aes(y = fit), color = "blue", size = 1) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  labs(
    title = "Linear Mixed-Effects Model with Prediction Intervals",
    x = "Years from Start Date (year2)",
    y = "Predicted Systolic Blood Pressure (systolic_bp)"
  ) +
  theme_minimal()

library(ggplot2)
library(ggeffects)

# Generate predictions with a smooth line across the entire range of years
predictions <- ggpredict(model_mx, terms = "year2")

# Filter the predictions to only include the range from -5.0 to 0.0
predictions <- predictions %>%
  filter(x >= -5.0 & x <= 0.0)

# Plotting with ggplot
bf=ggplot(predictions, aes(x = x, y = predicted)) +
  geom_line(color = "#FF00FF", linewidth = 1) +  # Smooth line representing the fixed effect
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +  # Confidence interval
  labs(
    title = "SBP Decline Over 5years",
    x = "Years",
    y = "Systolic Blood Pressure (mmHg)"
  ) +
  theme_minimal()

ggsave("best_fit_lmer_models_sbp_decline_cmri5y_140+_2024.08.19.png", bf, width=15, height=13, units="cm", dpi=200, bg='white')


##OThers
# Boxplot
ggplot(hpt_ukb_stroke2mx2, aes(x = year_vs, y = systolic_bp)) +
  geom_boxplot() +
  labs(x = "Year", y = "Systolic BP", title = "Boxplot of Systolic BP over Years")

# Histogram
ggplot(hpt_ukb_stroke2mx2, aes(x = systolic_bp)) +
  geom_histogram(binwidth = 5, fill = "skyblue", color = "black") +
  labs(x = "Systolic BP", y = "Frequency", title = "Histogram of Systolic BP")

# Density
ggplot(hpt_ukb_stroke2mx2, aes(x = systolic_bp)) +
  geom_density(fill = "skyblue", color = "black") +
  labs(x = "Systolic BP", y = "Density", title = "Density Plot of Systolic BP")

# Scatter Plot and Smooth Line
ggplot(hpt_ukb_stroke2mx2, aes(x = year_vs, y = systolic_bp)) +
  geom_point() +  # Add points for the data
  geom_smooth(method = "lm", se = FALSE) +  # Add a best fit line
  labs(x = "Year", y = "Systolic BP", title = "Scatter Plot with Smoothed Line")









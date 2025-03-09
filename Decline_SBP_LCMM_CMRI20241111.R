#######  Osei Asibey - 2024.11.11  ######
###########################################################
# Latent Class Mixture Module on GP_Clinical Data in Snow

# Set Directory
setwd("H:/Projects/Biobank/Osei")   

# data loading and manipulation
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

# Plot Packages
library(ggplot2)
install.packages("catspec")
#library(catspec)
install.packages("gridExtra")
library(gridExtra)


# Load Package for LCMM
if (!require("lcmm")) install.packages("lcmm")
library(lcmm)



###### loading Clean Data on Systolic BPs #####
# Loading Baseline UK Biobank Data and Executing Script

# load Data
hpt_ukb_sub = read.csv("H:/Projects/Biobank/Osei/Datafiles/hpt_ukb_sub1.csv")
#ukb <- read.csv("H:/Projects/Biobank/Osei/Datafiles/ukb.csv")
cmri_ukb <- read.csv("H:/Projects/Biobank/Osei/Datafiles/cmri_ukb1.csv")


hpt_ukb_sub <- left_join(hpt_ukb_sub , cmri_ukb, by = "n_eid_14631")

## Diastolic BP 
hpt_ukb_sub <- hpt_ukb_sub |>
  mutate(pp = as.numeric(systolic_bp - diastolic_bp))

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


###################################################################################################################

## Extracting Coefficients for Slopes ############################################################

# Linear
hpt_ukb_lcmm <- hpt_ukb_sub1|>
  select(n_eid_14631, event_dt, systolic_bp, diastolic_bp, pp, age_SSD, sex, year, year2, n_count)



str(hpt_ukb_lcmm)

# Explore
p1 <- ggplot(hpt_ukb_lcmm, aes(year2, systolic_bp, group=n_eid_14631)) + 
  geom_line() + 
  scale_y_continuous(limits = c(80,240)) + 
  labs(x = "years", y = "Systolic BP", title = "One Line per person, all subjects")
p1

p2 <- ggplot(hpt_ukb_lcmm[hpt_ukb_lcmm$n_eid_14631 %in% unique(hpt_ukb_lcmm$n_eid_14631)[1:300] & hpt_ukb_lcmm$n_count>5,], aes(year2, systolic_bp, group=n_eid_14631)) + 
  geom_line() + 
  scale_y_continuous(limits = c(80,240)) + 
  labs(x = "years", y = "Systolic BP", title = "Just 300 random subjects with >5 total obs")

p3 <- ggplot(hpt_ukb_lcmm[hpt_ukb_lcmm$n_eid_14631 %in% unique(hpt_ukb_lcmm$n_eid_14631)[1:300] & hpt_ukb_lcmm$n_count>5,], aes(year2, systolic_bp, group=n_eid_14631)) + 
  geom_smooth(aes(group=n_eid_14631), method="lm", se=F) + 
  scale_y_continuous(limits = c(80,240)) + 
  labs(x = "years", y = "Systolic BP", title = "Just 300 random subjects with >5 total obs\n straight line per person")
p3

p4 <- ggplot(hpt_ukb_lcmm[hpt_ukb_lcmm$n_eid_14631 %in% unique(hpt_ukb_lcmm$n_eid_14631)[1:300] & hpt_ukb_lcmm$n_count>20,], aes(year2, systolic_bp, group=n_eid_14631)) + 
  geom_smooth(aes(group=n_eid_14631), method="loess", se=F) + 
  scale_y_continuous(limits = c(80,240)) + 
  labs(x = "years", y = "Systolic BP", title = "Just 300 random subjects with >5 total obs\n smoothed line per person")
p4

grid.arrange(p1, p2, p3, p4, nrow=2)


# Latent Class Model
## Filters
hpt_lcmm <- hpt_ukb_sub1 |> filter(!is.na(date_of_cmri))           # Those with no CMRI Data
hpt_lcmm <-hpt_lcmm|> filter(year2>=-10)
hpt_lcmm <-hpt_lcmm |> filter(systolic_bp >=60 & systolic_bp<=240)
hpt_lcmm <-hpt_lcmm |> filter(n_count>4)
hpt_lcmm <-hpt_lcmm |> filter(age_SSD>=70)

# Filter participants with less than 4 BP measurements
hpt_lcmm2 <- hpt_lcmm[!hpt_lcmm$n_eid_14631 %in% names(which(table(hpt_lcmm$n_eid_14631)<4)),]


# Adding categorical variable to the ID, it changes it from specific ID to increading numerals
hpt_lcmm2$n_eid_14631 <- as.numeric(as.factor(hpt_lcmm2$n_eid_14631))

summary(as.vector(table(hpt_lcmm2$n_eid_14631)))
(n_distinct(hpt_lcmm2$n_eid_14631))

hpt_lcmm2 <- hpt_lcmm2|>
  select(n_eid_14631, event_dt, systolic_bp, diastolic_bp, pp, age_SSD, sex, year, year2, education, smoking_status, bmi_cats,ldl_Cat, stroke_ischemic_hes2022, hf_hes2022, mi_hes2022, pvd_hes2022, t2d_hes2022, dead, n_count)


# 3 Groups
d1 <- lcmm(systolic_bp ~ year2, random = ~ year2, subject = 'n_eid_14631', idiag = TRUE, data = hpt_lcmm2, link = "linear")

d2 <- lcmm(systolic_bp ~ year2, random = ~ year2, subject = 'n_eid_14631', mixture = ~year2, ng = 2, idiag = TRUE, data = hpt_lcmm2, link = "linear", B = d1)
d3 <- lcmm(systolic_bp ~ year2, random = ~ year2, subject = 'n_eid_14631', mixture = ~year2, ng = 3, idiag = TRUE, data = hpt_lcmm2, link = "linear", B = d1)
d4 <- lcmm(systolic_bp ~ year2, random = ~ year2, subject = 'n_eid_14631', mixture = ~year2, ng = 4, idiag = TRUE, data = hpt_lcmm2, link = "linear", B = d1)
d5 <- lcmm(systolic_bp ~ year2, random = ~ year2, subject = 'n_eid_14631', mixture = ~year2, ng = 5, idiag = TRUE, data = hpt_lcmm2, link = "linear", B = d1)
d6 <- lcmm(systolic_bp ~ year2, random = ~ year2, subject = 'n_eid_14631', mixture = ~year2, ng = 6, idiag = TRUE, data = hpt_lcmm2, link = "linear", B = d1)
d7 <- lcmm(systolic_bp ~ year2, random = ~ year2, subject = 'n_eid_14631', mixture = ~year2, ng = 7, idiag = TRUE, data = hpt_lcmm2, link = "linear", B = d1)

summary(d4)
postprob(d4) 

## Sline
s1 <- lcmm(systolic_bp ~ year2, random = ~ year2, subject = 'n_eid_14631', idiag = TRUE, data = hpt_lcmm2, link = "linear")
s4 <- lcmm(systolic_bp ~ year2, random = ~ year2, subject = 'n_eid_14631', mixture = ~year2, ng = 4, idiag = TRUE, data = hpt_lcmm2, link = "splines", B = d1)

ss1 <- lcmm(systolic_bp ~ year2, random = ~ year2, subject = 'n_eid_14631', idiag = TRUE, data = hpt_lcmm2, link = "splines")
ss4 <- lcmm(systolic_bp ~ year2, random = ~ year2, subject = 'n_eid_14631', mixture = ~year2, ng = 4, idiag = TRUE, data = hpt_lcmm2, link = "splines", B = ss1)

summary(ss4)
postprob(ss4) 


# look at the post probs closer
round(summary(as.numeric(d2$pprob[d2$pprob[,"class"] == 1, "prob1"])), 2)
round(summary(as.numeric(d2$pprob[d2$pprob[,"class"] == 2, "prob2"])), 2)
round(summary(as.numeric(d2$pprob[d2$pprob[,"class"] == 3, "prob2"])), 2)


hpt_tab <- as.data.frame(d2$pprob[,1:2])
hpt_lcmm2$group2 <-  factor(hpt_tab$class[match(hpt_lcmm2$n_eid_14631, hpt_tab$n_eid_14631)])

hpt_tab <- as.data.frame(d3$pprob[,1:2])
hpt_lcmm2$group3 <-  factor(hpt_tab$class[match(hpt_lcmm2$n_eid_14631, hpt_tab$n_eid_14631)])

hpt_tab <- as.data.frame(d4$pprob[,1:2])
hpt_lcmm2$class <-  factor(hpt_tab$class[match(hpt_lcmm2$n_eid_14631, hpt_tab$n_eid_14631)])

hpt_tab <- as.data.frame(d5$pprob[,1:2])
hpt_lcmm2$group5 <-  factor(hpt_tab$class[match(hpt_lcmm2$n_eid_14631, hpt_tab$n_eid_14631)])

hpt_tab <- as.data.frame(d6$pprob[,1:2])
hpt_lcmm2$group6 <-  factor(hpt_tab$class[match(hpt_lcmm2$n_eid_14631, hpt_tab$n_eid_14631)])

hpt_tab <- as.data.frame(d7$pprob[,1:2])
hpt_lcmm2$group7 <-  factor(hpt_tab$class[match(hpt_lcmm2$n_eid_14631, hpt_tab$n_eid_14631)])



# Ensure hpt_tab$class is named with the corresponding IDs for easy matching
names(hpt_tab$class) <- hpt_tab$n_eid_14631

# Map 'class' to hpt_lcmm2$group4 based on n_eid_14631
#hpt_lcmm2$group2 <- factor(hpt_tab$class[as.character(hpt_lcmm2$n_eid_14631)])


summary_table <- hpt_lcmm2 %>%
  group_by(class) %>%
  summarise(
    Count = n_distinct(n_eid_14631)
  )

# Print the summary table
print(summary_table)


##### LCMMM for outcome
# Convert hpt_dx_hx to a factor
hpt_lcmm2$class <- as.factor(hpt_lcmm2$class)

# Relevel the factor to set '0' as the reference level
hpt_lcmm2$class <- relevel(hpt_lcmm2$class, ref = "1")

# Linear Model
model_lm <- lm(pwv_aix ~ class +age_SSD + sex, data = hpt_lcmm2)
summary(model_lm)
confint(model_lm, level = 0.95)


###### Creating Table for Graph
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





# plot check
p1 <- ggplot(hpt_lcmm2, aes(year2, systolic_bp, group=n_eid_14631)) + 
  geom_line() + 
  geom_smooth(aes(group=group3), method="loess", linewidth=2, se=F)  + 
  scale_y_continuous(limits = c(80,240)) + 
  labs(x = "years", y = "Systolic BP", colour="Latent Class", title="Raw")
p1

p2 <- ggplot(hpt_lcmm2, aes(year2, systolic_bp, group=n_eid_14631)) +
  geom_smooth(aes(group=n_eid_14631, colour=group3),linewidth=0.5, se=F) +
  geom_smooth(aes(group=group3), method="loess", size=2, se=T) + 
  scale_y_continuous(limits = c(80,240)) + 
  labs(x = "years", y = "Systolic BP", colour="Latent Class", title="Smoothed") +
  theme(legend.position="none")
p2

##########################################
p2 <- ggplot(hpt_lcmm2, aes(year2, systolic_bp, group = n_eid_14631)) +
  geom_smooth(aes(group = n_eid_14631, colour = group3), linewidth = 0.5, se = FALSE) + 
  geom_smooth(aes(group = group2), method = "gam", formula = y ~ s(x, bs = "cs"), linewidth = 2, se = TRUE) +
  scale_y_continuous(limits = c(80, 240)) + 
  labs(x = "Years", y = "Systolic BP", colour = "Latent Class", title = "Smoothed") +
  theme(legend.position = "none")
p2

p3 <- ggplot(hpt_lcmm2, aes(year2, systolic_bp, group = n_eid_14631)) +
  geom_smooth(aes(group = n_eid_14631, colour = group3), linewidth = 0.5, se = FALSE) + 
  geom_smooth(aes(group = group3), method = "gam", formula = y ~ s(x, bs = "cs"), linewidth = 2, se = TRUE) +
  scale_y_continuous(limits = c(80, 240)) + 
  labs(x = "Years", y = "Systolic BP", colour = "Latent Class", title = "Smoothed") +
  theme(legend.position = "none")
p3

p4 <- ggplot(hpt_lcmm2, aes(year2, systolic_bp, group = n_eid_14631)) +
  geom_smooth(aes(group = n_eid_14631, colour = group3), linewidth = 0.5, se = FALSE) + 
  geom_smooth(aes(group = group4), method = "gam", formula = y ~ s(x, bs = "cs"), linewidth = 2, se = TRUE) +
  scale_y_continuous(limits = c(80, 240)) + 
  labs(x = "Years", y = "Systolic BP", colour = "Latent Class", title = "Smoothed") +
  theme(legend.position = "none")
p4

p5 <- ggplot(hpt_lcmm2, aes(year2, systolic_bp, group = n_eid_14631)) +
  geom_smooth(aes(group = n_eid_14631, colour = group3), linewidth = 0.5, se = FALSE) + 
  geom_smooth(aes(group = group5), method = "gam", formula = y ~ s(x, bs = "cs"), linewidth = 2, se = TRUE) +
  scale_y_continuous(limits = c(80, 240)) + 
  labs(x = "Years", y = "Systolic BP", colour = "Latent Class", title = "Smoothed") +
  theme(legend.position = "none")
p5

grid.arrange(p2, p3, p4, p5,  ncol = 2, top = "3 Latent Classes")

#########


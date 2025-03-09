
############################ Uni-variate, and Cross Tabulation Analysis####

# Summary Statistics and Cross Tabulation
summary_stats <- hpt_ukb_submx %>%
  group_by(Decline_SBP_Cat, cerevd_hes2022) %>%
  summarize(count = n()) %>%
  mutate(percent = count / sum(count) * 100)

summary_stats2 <- hpt_ukb_stroke2a %>%
  group_by(cerevd_hes2022) %>%
  summarize(count = n()) %>%
  mutate(percent = count / sum(count) * 100)

summary_stats <- hpt_ukb_stroke2a %>%
  group_by(Decline_SBP_Cat) %>%
  summarize(
    count = n(),
    mean_age = mean(systolic_bp, na.rm = TRUE),
    median_age = median(systolic_bp, na.rm = TRUE),
    sd_age = sd(systolic_bp, na.rm = TRUE)
  ) %>%
  mutate(percent = count / sum(count) * 100)

summary_stats2 <- hpt_ukb_stroke2a %>%
  #group_by(Decline_SBP_Cat) %>%
  summarize(
    count = n(),
    mean_age = mean(systolic_bp, na.rm = TRUE),
    median_age = median(systolic_bp, na.rm = TRUE),
    sd_age = sd(systolic_bp, na.rm = TRUE)
  ) %>%
  mutate(percent = count / sum(count) * 100)


####################### Test For Univariate Association (P-Values)######
## ANOVA summary
anova_result <- aov(age_SSD ~ Decline_SBP, data = hpt_ukb_stroke2a)
anova_result <- aov(systolic_bp ~ Decline_SBP, data = hpt_ukb_stroke2a)

# Age and blood pressure changes (Decline_SBP_Cat) - Analysis of Variance (ANOVA)
summary(aov(systolic_bp ~ Decline_SBP_Cat, data = hpt_ukb_stroke2a))
#  The dependent(Outcome) variable age is a continuous variable
#  THe independent (Predictor/Exposure) variable is a categorical variable 
#  So ANOVA is the best statistical methods used for modeling relationships between variables to compare means score in the grps


## Chi-square test (Categorical to Categorical)

#  BMI
(chi_sq_result <- chisq.test( table(hpt_ukb_stroke2a$Decline_SBP_Cat, hpt_ukb_stroke2a$bmi_cats)))
#  Smoking Status
(chi_sq_result <- chisq.test( table(hpt_ukb_stroke2a$Decline_SBP_Cat, hpt_ukb_stroke2a$smoking_status)))
#  Education 
(chi_sq_result <- chisq.test( table(hpt_ukb_stroke2a$Decline_SBP_Cat, hpt_ukb_stroke2a$education)))
#  Ethnicity
(chi_sq_result <- chisq.test( table(hpt_ukb_stroke2a$Decline_SBP_Cat, hpt_ukb_stroke2a$ethnicity)))
#  Sex
(chi_sq_result <- chisq.test( table(hpt_ukb_stroke2a$Decline_SBP_Cat, hpt_ukb_stroke2a$sex)))
#  LDL_Cat
(chi_sq_result <- chisq.test( table(hpt_ukb_stroke2a$Decline_SBP_Cat, hpt_ukb_stroke2a$ldl_Cat)))


## Logistic Regression

#  Stroke
summary(glm(stroke_ischemic_hes2022 ~ Decline_SBP_Cat, data = hpt_ukb_stroke2a, family = "binomial"))
#  Heart failure
summary(glm(hf_hes2022 ~ Decline_SBP_Cat, data = hpt_ukb_stroke2a, family = "binomial"))
#  MI
summary(glm(mi_hes2022 ~ Decline_SBP_Cat, data = hpt_ukb_stroke2a, family = "binomial"))
# Perepheral Vascular Dx
summary(glm(pvd_hes2022 ~ Decline_SBP_Cat, data = hpt_ukb_stroke2a, family = "binomial"))
#  Mortality
summary(glm(dead ~ Decline_SBP_Cat, data = hpt_ukb_stroke2a, family = "binomial"))



3. # Calculting Summary values
summary_stats <- hpt_ukb_stroke2a %>%
  group_by(Decline_SBP_Cat, sex) %>%
  summarize(count = n()) %>%
  mutate(percent = count / sum(count) * 100)

summary_stats2 <- hpt_ukb_stroke2a %>%
  group_by(sex) %>%
  summarize(count = n()) %>%
  mutate(percent = count / sum(count) * 100)

# Print summary statistics
print(summary_stats)
print(summary_stats2)

########













## Osei Asibey Owusu -- 2024.01.10
## Exploring UK Biobank Baseline Assessment 


# data loading and manipulation
install.packages("tidyverse")
library(tidyverse)
library(haven)
library(lukesRlib)
library(openxlsx)
library(vroom)
library(dplyr)

# set default ggplot theme
theme_set(theme_bw())

# load packages for survival analysis and plotting
library(survival)
library(survivalROC)
library(ggsurvfit)
library(tidycmprsk)

load("H:/Projects/BioBank/14631_ageing-well/data downloaded 2021 October - update all fields/ukb48975_0_sub.RDat")
load("H:/Projects/BioBank/14631_ageing-well/data downloaded 2021 October - update all fields/ukb48975_0.RDat")

## rename age/sex and other commonly used
ukb_0_sub = ukb_0_sub |> rename(
  assessment_date = ts_53_0_0,
  month_of_birth = n_52_0_0,
  year_of_birth = n_34_0_0,
  age = n_21022_0_0,
  sex = n_31_0_0,
  height = n_50_0_0,
  weight = n_21002_0_0,
  bmi = n_21001_0_0,
  waist = n_48_0_0,
  hip = n_49_0_0,
  bmd_heel = n_3148_0_0,
  grip_left = n_46_0_0,
  grip_right = n_47_0_0
)

## rename haematology
ukb_0_sub = ukb_0_sub |> rename(
  hbc = n_30020_0_0,
  mcv = n_30040_0_0,
  rdw = n_30070_0_0,
  mrv = n_30260_0_0,
  ret_imm_frac = n_30280_0_0,
  mpv = n_30100_0_0,
  pdw = n_30110_0_0,
  plt_crit = n_30090_0_0,
  mscv = n_30270_0_0,
  mch = n_30050_0_0,
  mch_c = n_30060_0_0,
  rbc_n = n_30010_0_0,
  rbc_p = n_30230_0_0,
  ret_p = n_30240_0_0,
  ret_n = n_30250_0_0,
  plt_n = n_30080_0_0,
  rbcn_n = n_30170_0_0,
  hct_p = n_30030_0_0,
  wbc_n = n_30000_0_0,
  neu_n = n_30140_0_0,
  neu_p = n_30200_0_0,
  lym_n = n_30120_0_0,
  lym_p = n_30180_0_0,
  mon_n = n_30130_0_0,
  mon_p = n_30190_0_0,
  eos_n = n_30150_0_0,
  eos_p = n_30210_0_0,
  bas_n = n_30160_0_0,
  bas_p = n_30220_0_0
)

## rename blood biomarkers
ukb_0_sub = ukb_0_sub |> rename(
  albumin = n_30600_0_0,
  alkaline_phos = n_30610_0_0,
  alanine_amino = n_30620_0_0,
  apolipoprot_a = n_30630_0_0,
  apolipoprot_b = n_30640_0_0,
  aspartate_ami = n_30650_0_0,
  bilirubin_tot = n_30660_0_0,
  urea = n_30670_0_0,
  calcium = n_30680_0_0,
  cholesterol = n_30690_0_0,
  creatinine = n_30700_0_0,
  crp = n_30710_0_0,
  cystatin_c = n_30720_0_0,
  glutamyltrans = n_30730_0_0,
  glucose = n_30740_0_0,
  hba1c = n_30750_0_0,
  hdl = n_30760_0_0,
  igf1 = n_30770_0_0,
  ldl = n_30780_0_0,
  lipoprotein_a = n_30790_0_0,
  oestradiol = n_30800_0_0,
  phosphate = n_30810_0_0,
  rheumatoid_fa = n_30820_0_0,
  shbg = n_30830_0_0,
  bilirubin_dir = n_30840_0_0,
  testosterone = n_30850_0_0,
  total_protein = n_30860_0_0,
  triglycerides = n_30870_0_0,
  urate = n_30880_0_0,
  vitamin_d = n_30890_0_0
)

## generate date of birth, waist:hip, etc
ukb_0_sub = ukb_0_sub |> mutate(
  date_of_birth = dmy(paste0(15,"-",month_of_birth,"-",year_of_birth)),
  whr = waist/hip,
  gripmax = grip_left,
  gripmax = if_else(!is.na(grip_right) & grip_right > grip_left, grip_right, gripmax),
  dynapenia_ewgsop = case_when(
    sex == 0 & gripmax >= 20 ~ 0,
    sex == 0 & gripmax <  20 ~ 1,
    sex == 1 & gripmax >= 30 ~ 0,
    sex == 1 & gripmax <  30 ~ 1,
    TRUE ~ NA_real_)
)

## BMI categories
ukb_0_sub = ukb_0_sub |> mutate(bmi_cats = as_factor(case_when(
  bmi >= 18.5 & bmi < 25 ~ 1,
  bmi >= 25 & bmi < 30 ~ 2,
  bmi >= 30 ~ 3,
  bmi < 18.5 ~ 4)),
  bmi_cats = labelled(bmi_cats, c(Normal = 1, Overweight = 2, Obese = 3, Underweight = 4), label="BMI (categories)")
)

## ethnicity
ukb_0_sub = ukb_0_sub |> mutate(ethnicity = as_factor(case_when(
  n_21000_0_0 %in% c(1,1001,1002,1003) ~ 1,
  n_21000_0_0 %in% c(2,2001,2002,2003,2004) ~ 2,
  n_21000_0_0 %in% c(3,3001,3002,3003,3004) ~ 3,
  n_21000_0_0 %in% c(4,4001,4002,4003,4004) ~ 4,
  n_21000_0_0 %in% c(5) ~ 5,
  n_21000_0_0 %in% c(6) ~ 6)),
  ethnicity = labelled(ethnicity, c(White = 1, Mixed = 2, Asian = 3, Black = 4, Chinese = 5, Other = 6), label="Ethnicity (self-report)")
)

## education
ukb_0_sub = ukb_0_sub |> mutate(education = case_when(
  n_6138_0_0==-7 ~ 0,
  n_6138_0_0==4 | n_6138_0_1==4 | n_6138_0_2==4 | n_6138_0_3==4 | n_6138_0_4==4 | n_6138_0_5==4 ~ 1,
  n_6138_0_0==3 | n_6138_0_1==3 | n_6138_0_2==3 | n_6138_0_3==3 | n_6138_0_4==3 | n_6138_0_5==3 ~ 2,
  n_6138_0_0==2 | n_6138_0_1==2 | n_6138_0_2==2 | n_6138_0_3==2 | n_6138_0_4==2 | n_6138_0_5==2 ~ 3,
  n_6138_0_0==5 | n_6138_0_1==5 | n_6138_0_2==5 | n_6138_0_3==5 | n_6138_0_4==5 | n_6138_0_5==5 ~ 3,
  n_6138_0_0==6 | n_6138_0_1==6 | n_6138_0_2==6 | n_6138_0_3==6 | n_6138_0_4==6 | n_6138_0_5==6 ~ 4,
  n_6138_0_0==1 | n_6138_0_1==1 | n_6138_0_2==1 | n_6138_0_3==1 | n_6138_0_4==1 | n_6138_0_5==1 ~ 5),
  education = labelled(education, c(None=0, CSEs=1, GCSEs=2, Alevels=3, ProfQual=4, Degree=5), label="Highest education level attained")
)

## smoking status
ukb_0_sub = ukb_0_sub |> mutate(smoking_status = if_else(
  n_20116_0_0 %in% c(0,1,2), n_20116_0_0, NA_real_),
  smoking_status = labelled(smoking_status, c(Never=1, Previous=2, Current=3), label="Smoking Status (Self-Report")
)

## blood pressure
ukb_0_sub = ukb_0_sub |> mutate(
  sbp_0_avg = case_when(
    !is.na(n_4080_0_0) & !is.na(n_4080_0_1) ~ (n_4080_0_0 + n_4080_0_1)/2,
    !is.na(n_93_0_0) & !is.na(n_93_0_1) ~ (n_93_0_0 + n_93_0_1)/2),
  dbp_0_avg = case_when(
    !is.na(n_4079_0_0) & !is.na(n_4079_0_1) ~ (n_4079_0_0 + n_4079_0_1)/2,
    !is.na(n_94_0_0) & !is.na(n_94_0_1) ~ (n_94_0_0 + n_94_0_1)/2),
  
  sbp_0_min = pmin(n_4080_0_0, n_4080_0_1, n_93_0_0, n_93_0_1, na.rm=TRUE),
  dbp_0_min = pmin(n_4079_0_0, n_4079_0_1, n_94_0_0, n_94_0_1, na.rm=TRUE)
)

## had menopause
ukb_0_sub = ukb_0_sub |> mutate(
  had_menopause=case_when(n_2724_0_0==0~0,  n_2724_0_0==1~1),
  had_menopause_with_hyst=case_when(n_2724_0_0==0~0,  n_2724_0_0==1~1,  n_2724_0_0==2~1),
  age_at_menopause=case_when(n_3581_0_0>0~n_3581_0_0, TRUE~NA_real_),
  years_since_menopause=case_when(age>age_at_menopause~age-age_at_menopause, TRUE~NA_real_)
)

## parents lifespan
ukb_0_sub = ukb_0_sub |> mutate(
  fathers_age_death = if_else(n_1807_0_0>10, n_1807_0_0, NA_real_),
  mothers_age_death = if_else(n_3526_0_0>10, n_1807_0_0, NA_real_),
  fathers_age_death_z = lukesRlib::z_trans(fathers_age_death),
  mothers_age_death_z = lukesRlib::z_trans(mothers_age_death),
  parents_age_death_u = if_else(!is.na(fathers_age_death) & !is.na(mothers_age_death), (fathers_age_death+mothers_age_death)/2, NA_real_),
  parents_age_death_zu = if_else(!is.na(fathers_age_death_z) & !is.na(mothers_age_death_z), (fathers_age_death_z+mothers_age_death_z)/2, NA_real_)
)


##
## updated/recoded assessment info 
ukb_t = read_dta("H:/Projects/BioBank/14631_ageing-well/data downloaded 2021 October - update all fields/ukb48975.assessment_info.recode_0.dta")
ukb_0_sub = left_join(ukb_0_sub, ukb_t, by="n_eid_14631")
rm(ukb_t)

##
## function to drop the non-renamed cols...
drop_extra_ukb_cols = function(x)  {
  keepers = ! grepl("^n_|^ts_",colnames(x))
  keepers[1] = TRUE
  x = x[,keepers]
  x
}

## Dropping Extra Columns
ukb = drop_extra_ukb_cols(ukb_0_sub)
#ukb = ukb_0_sub
rm(ukb_0_sub)


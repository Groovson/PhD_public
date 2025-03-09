#######  Osei Asibey - 2025.03.01  ######
###########################################################
# Blood Pressure Variables in GP_Clinical Data

# data loading and manipulation
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

# load UK Biobank baseline data --------------------------------------------------------------------
# external script renames/formats common variables
#Setting Directory                      
setwd("I:/Projects/Osei")           
#Obtaining data from source
ukb_hpt_20231026 = read_dta("H:/Projects/BioBank/14631_ageing-well/GP up to 2017 June/gp_clinical_subset.dta")
ukb_hpt = read_dta("H:/Projects/BioBank/14631_ageing-well/GP up to 2017 June/gp_clinical.dta")
hpt_ukb = ukb_hpt

## Importing Source Codes
source("H:/Projects/BioBank/Olivia/Code_lists/_get_diagnostic_codes.20220711_3.R")
ctv3_codelist<-get_ctv3_matches(your_read_codelist$h) 

rm(lkp_read2, map_ctv3_icd10, map_ctv3_icd9, map_icd9_icd10, map_read2_icd10, map_read2_icd9)

## Importing read Codes
read_2_codelist <- read.xlsx("H:/Projects/BioBank/Osei/Codes/read_2_codingOAO.xlsx")
read_2_code_sbp <- read.xlsx("I:/Projects/Osei/Codes/BP_read_codes2_SBP.xlsx")
read_2_code_GenBP <- read.xlsx("I:/Projects/Osei/Codes/BP_read_codes2_GeneralBP.xlsx")
read_2_code_dbp <- read.xlsx("I:/Projects/Osei/Codes/BP_read_codes2_DBP.xlsx")

### Filtering read_2 BP Codes from (map_read2_ctv_3)
map_read2_ctv3_bpSBP <- merge(map_read2_ctv3, read_2_code_sbp, by.x="READV2_CODE",by.y="READV2_CODES")
map_read2_ctv3_bpDBP <- merge(map_read2_ctv3, read_2_code_dbp, by.x="READV2_CODE",by.y="READV2_CODES")
map_read2_ctv3_bpGen <- merge(map_read2_ctv3, read_2_code_GenBP, by.x="READV2_CODE",by.y="READV2_CODES")

map_read2_ctv3_bpGen<- map_read2_ctv3_bpGen[-c(1, 3), ] 

################################################################################
0# Working with Both Systolic_bp, Diastolic BP and Gen_BP
# remove excess readcodes
rm(lkp_ctv3, lkp_icd10, lkp_icd9, lkp_read2, map_ctv3_icd10, map_ctv3_icd9, map_icd9_icd10, 
   map_read2_icd10, map_read2_icd9)

### Filtering read_2 BP Codes from (map_read2_ctv_3)
map_read2_ctv3_bp <- merge(map_read2_ctv3, read_2_codelist, by.x="READV2_CODE",by.y="READV2_CODES")

map_read2_ctv3_bp <- map_read2_ctv3_bp[-c(1, 3), ] 

# Merging Read codes with Data
hpt_ukbV2<-merge(hpt_ukb,map_read2_ctv3_bp,by.x="read_2",by.y="READV2_CODE")   
hpt_ukbV3<-merge(hpt_ukb,map_read2_ctv3_bp,by.x="read_3",by.y="READV3_CODE") 

colnames(hpt_ukbV2)         
colnames(hpt_ukbV3)

names(hpt_ukbV2) %in% names(hpt_ukbV3)
names(hpt_ukbV3) %in% names(hpt_ukbV2)

#Renaming READV2_CODE AND READV3_CODE AS READV2V3_CODES
hpt_ukbV2 <- hpt_ukbV2 |> 
  rename("READV2V3_CODE"="READV3_CODE")
hpt_ukbV3 = hpt_ukbV3 |> 
  rename("READV2V3_CODE"="READV2_CODE")

#Combine from both READV2 AND V3 CODES
hpt_ukb_v <- rbind(hpt_ukbV2, hpt_ukbV3)

col_remove <- c("value3", "CHAPTER", "READV2_DESC.x", "TERMV2_ORDER", "TERMV2_TYPE", "TERMV3_CODE", 
                "TERMV3_DESC", "IS_ASSURED", "READV2_DESC.y", "TERMV3_TYPE", "KEYS", "X5")
hpt_ukb <- hpt_ukb_v |> 
  select(- one_of(col_remove))

rm(map_ctv3_read2, map_read2_ctv3, map_read2_ctv3_bp, map_read2_ctv3_bpDBP, map_read2_ctv3_bpGen, map_read2_ctv3_bpSBP, read_2_code_dbp, read_2_code_GenBP, read_2_code_sbp, read_2_codelist)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Creating a column systolic_bp to contain all values of Systolic BP

# Converting into whole values
hpt_ukb$value1 <- as.numeric(hpt_ukb$value1)
hpt_ukb$value2 <- as.numeric(hpt_ukb$value2)

# Filling Missing variables with NA
hpt_ukb$value1[hpt_ukb$value1==""]<-NA
hpt_ukb$value2[hpt_ukb$value2==""]<-NA

hpt_ukb$value1[hpt_ukb$value1==0]<-NA
hpt_ukb$value2[hpt_ukb$value2==0]<-NA



#hpt_ukbV2V3sbp2$systolic_bp <-as.numeric(hpt_ukbV2V3sbp2$systolic_bp)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Cleaning Data

# Filtering out rows with missing values in both value1 and value2 
hpt_ukb <- hpt_ukb %>%
  filter(!is.na(as.numeric(value1)) | !is.na(as.numeric(value2)))

1# Creating a Data for Systolic_BP
################################################################################

# Filter Only systolic BPs
hpt_ukb_sbp <- hpt_ukb |>
filter(BP_CATEGORIES==1)

# Creating a systolic_bp column
hpt_ukb_sbp$systolic_bp <- hpt_ukb_sbp$value1
hpt_ukb_sbp$systolic_bp <- if_else(is.na(hpt_ukb_sbp$systolic_bp), hpt_ukb_sbp$value2, hpt_ukb_sbp$systolic_bp)

1b# Working on others designated as Systolic BPs (99)
hpt_ukb_sbp2 <- hpt_ukb |>
  filter(READV2V3_CODE=="246K."| READV2V3_CODE=="XaI9f") |>
  filter(is.na(as.numeric(value1)) | is.na(as.numeric(value2)))

hpt_ukb_sbp2$systolic_bp <- hpt_ukb_sbp2$value1
hpt_ukb_sbp2$systolic_bp <- if_else(is.na(hpt_ukb_sbp2$systolic_bp), hpt_ukb_sbp2$value2, hpt_ukb_sbp2$systolic_bp)

hpt_ukb_sbp <- rbind(hpt_ukb_sbp, hpt_ukb_sbp2)



2# Creating a Data for Systolic_BP
################################################################################

# Filter Only systolic BPs
hpt_ukb_dbp <- hpt_ukb |>
  filter(BP_CATEGORIES==2)

# Creating a systolic_bp column
hpt_ukb_dbp$diastolic_bp <- hpt_ukb_dbp$value1
hpt_ukb_dbp$diastolic_bp <- if_else(is.na(hpt_ukb_dbp$diastolic_bp), hpt_ukb_dbp$value2, hpt_ukb_dbp$diastolic_bp)

2b# Working on others designated as Systolic BPs (99)
hpt_ukb_dbp2 <- hpt_ukb |>
  filter(READV2V3_CODE=="246L."| READV2V3_CODE=="XaI9g") |>
  filter(is.na(as.numeric(value1)) | is.na(as.numeric(value2)))

hpt_ukb_dbp2$diastolic_bp <- hpt_ukb_dbp2$value1
hpt_ukb_dbp2$diastolic_bp <- if_else(is.na(hpt_ukb_dbp2$diastolic_bp), hpt_ukb_dbp2$value2, hpt_ukb_dbp2$diastolic_bp)

hpt_ukb_dbp <- rbind(hpt_ukb_dbp, hpt_ukb_dbp2)


3# "Full_Join of both systolic BP data and Diastolic BP data to get systolic_bp along side)
################################################################################

# Creating a unique identifier for each measurement within each dataset
hpt_ukb_sbp <- hpt_ukb_sbp|>
  arrange(n_eid_14631, event_dt, desc(systolic_bp))|>
  mutate(measurement_id_sbp = row_number())

hpt_ukb_dbp <- hpt_ukb_dbp|>
  arrange(n_eid_14631, event_dt, desc(diastolic_bp))|>
  mutate(measurement_id_dbp = row_number())

# Performing the join
hpt_ukb_j <- full_join(hpt_ukb_sbp, hpt_ukb_dbp, by = c("n_eid_14631", "event_dt", "measurement_id_sbp" = "measurement_id_dbp"))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Removing excess column
hpt_ukb_j <- hpt_ukb_j|>
  rename("READV2V3_CODE" = "READV2V3_CODE.x","BP_CATEGORIES" = "BP_CATEGORIES.x")

hpt_ukb_j <- hpt_ukb_j|>
  select(n_eid_14631, event_dt, READV2V3_CODE, BP_CATEGORIES,systolic_bp, diastolic_bp)

4# Working the GeneralBP records with values in both value1 and 2
################################################################################

hpt_ukb_g <- hpt_ukb|>
  filter(BP_CATEGORIES==3 | BP_CATEGORIES==99)

# Creating new variable systolic_bp (Maximum value of value1and2) and diastolic_bp(Minimum value of value1and2) 
hpt_ukb_g$systolic_bp <- if_else(hpt_ukb_g$value1 > hpt_ukb_g$value2, hpt_ukb_g$value1, hpt_ukb_g$value2)
hpt_ukb_g$diastolic_bp <- if_else(hpt_ukb_g$value1 < hpt_ukb_g$value2, hpt_ukb_g$value1, hpt_ukb_g$value2)

hpt_ukb_g <- hpt_ukb_g |> 
  select(n_eid_14631, event_dt, READV2V3_CODE, BP_CATEGORIES, systolic_bp, diastolic_bp)|>
  filter(!is.na(systolic_bp) & !is.na(diastolic_bp))

5# Combine the two Data, 1. Joint systolic and diastolic (hpt_ukbV2V3_joint) and General ()
################################################################################
hpt_ukb1 <- rbind(hpt_ukb_j, hpt_ukb_g)

rm(hpt_ukb_dbp, hpt_ukb_dbp2, hpt_ukb_g, hpt_ukb_j, hpt_ukb_sbp, hpt_ukb_sbp2, hpt_ukb_v, hpt_ukbGen_check, hpt_ukbGen_check2, hpt_ukbV2, hpt_ukbV3)





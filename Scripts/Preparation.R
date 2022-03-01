# -------------------------------------------------------
# Preparation--------------------------------------------
# -------------------------------------------------------


# 1 Preparation------------------------------------------

# 1.1 Setup----------------------------------------------

## Load Packages
library("tidyverse")
library("readxl")


## Set Slurm directory
setwd("/home/nkappelmann/OPTIMA/ImmunePrediction")

## Load psychometric data
# Slurmgate
#load("/home/nkappelmann/OPTIMA/OPTIMA_Analyses/Data/Processed/02_PsychometricData.RData")

# MPI local
load("./Data/02_PsychometricData.RData")

## Load cytokine data and related IDs
# Slurmgate
#cyto_dat = readRDS("/binder/mgp/datasets/2020_ImmuneDepression/cytokine/02_data_for_analysis_imputed.rds")
#cyto_IDs = read_excel("/binder/mgp/datasets/2020_ImmuneDepression/cytokine/01_become_optima_randomization_STD_Controls_Pheno_1219.xlsx")

# MPI local
cyto_dat = readRDS("./Data/02_data_for_analysis_imputed.rds")
cyto_IDs = read_excel("./Data/01_become_optima_randomization_STD_Controls_Pheno_1219.xlsx")


## Load BMI data
bmi = read.csv(file = "./Data/BMI.csv", header = TRUE, sep = ";")


# Save biobank study IDs
cyto_dat$ID.bio = cyto_dat$Row.names

# Merge OPTIMA IDs
cyto_dat = merge(cyto_dat, cyto_IDs[, c("ID", "Lagerung_BC.Tube")], 
                 by.x = "ID.bio", by.y = "Lagerung_BC.Tube", all.x = TRUE)

# Delete duplicates from merging
cyto_dat = cyto_dat[!duplicated(cyto_dat$ID),]

# Save variable names in cytokine reference
cyto_ref = data.frame(vars = NA,
                      labels = colnames(cyto_dat)[10:ncol(cyto_dat)],
                      stringsAsFactors = FALSE)

# Change column names by excluding hyphens and slashes
colnames(cyto_dat) = gsub("/", "_", colnames(cyto_dat), fixed = TRUE)
colnames(cyto_dat) = gsub("-", ".", colnames(cyto_dat), fixed = TRUE)
colnames(cyto_dat) = gsub(" ", "_", colnames(cyto_dat), fixed = TRUE)
cyto_ref$vars = colnames(cyto_dat)[10:ncol(cyto_dat)]
cyto_ref = cyto_ref[cyto_ref$vars != "ID",]

## Remove BeCOME data
cyto_dat = cyto_dat[grepl("PTP", cyto_dat$ID, fixed = TRUE),]

## Remove OPTIMA participants without cytokine data
# Note: PTP0114 is excluded as this person clicked through all of their outcome data.
dat = dat[dat$ID %in% cyto_dat$ID,]

## Remove age and sex, so these don't get duplicated
dat[, c("sex", "age")] = NULL

## Merge psychometric and cytokine data
dat = merge(dat, cyto_dat, by = "ID", all.x = TRUE)


## Merge BMI data
dat = merge(dat, bmi, by = "ID", all.x = TRUE)


# 1.2 Coding---------------------------------------------

## Individual diagnoses
dat = dat %>%
   mutate(
      t0_cidi_depression_episode = case_when(
         t0_DF320 == 1 | t0_DF321 == 1 | t0_DF322 == 1 | t0_DF323 == 1 ~ 1,
         is.na(t0_cidi_diagnsum) ~ NA_real_,
         TRUE ~ 0),
      t0_cidi_depression_recurrent = case_when(
         t0_DF330 == 1 | t0_DF331 == 1 | t0_DF332 == 1 | t0_DF333 == 1 ~ 1,
         is.na(t0_cidi_diagnsum) ~ NA_real_,
         TRUE ~ 0),
      t0_cidi_anxiety = case_when(
         t0_DF4000 == 1 | t0_DF4001 == 1 | t0_DF401 == 1 | t0_DF4021 == 1 |
            t0_DF4022 == 1 | t0_DF4023 == 1 | t0_DF4024 == 1 | t0_DF4025 == 1 |
            t0_DF409 == 1 | t0_DF410 == 1 | t0_DF411 == 1 | t0_DF428 == 1 |
            t0_DF431 == 1 | t0_DF451 == 1 | t0_DF452 == 1 | t0_DF454 == 1 |
            t0_DF459 == 1 ~ 1,
         is.na(t0_cidi_diagnsum) ~ NA_real_,
         TRUE ~ 0),
      t0_cidi_substance = case_when(
         t0_DF101 == 1 | t0_DF102 == 1 | t0_DF111 == 1 | t0_DF112 == 1 |
            t0_DF121 == 1 | t0_DF122 == 1 | t0_DF131 == 1 | t0_DF132 == 1 |
            t0_DF141 == 1 | t0_DF151 == 1 | t0_DF152 == 1 | t0_DF161 == 1 |
            t0_DF172 == 1 ~ 1,
         is.na(t0_cidi_diagnsum) ~ NA_real_,
         TRUE ~ 0)
   )

## Create number of previous diagnoses variable divided by age and standardise
dat$t0_diagn_by_age = dat$t0_cidi_diagnsum / dat$age


## PID domains
dat$t0_pid_antago = rowMeans(dat[, paste0("t0_pid_", c("mani", "dece", "gran"))])
dat$t0_pid_negaff = rowMeans(dat[, paste0("t0_pid_", c("anxi", "emot", "sepa"))])
dat$t0_pid_disinh = rowMeans(dat[, paste0("t0_pid_", c("irre", "impu", "dist"))])
dat$t0_pid_detach = rowMeans(dat[, paste0("t0_pid_", c("anhe", "inti", "with"))])
dat$t0_pid_psycho = rowMeans(dat[, paste0("t0_pid_", c("ecce", "perc", "unus"))])

## Get Top 27% CRP cut-off
dat$hsCRP_inflamed = factor(ifelse(dat$CRP > quantile(dat$hsCRP, prob=1-27/100), 
                                 "Inflamed", "Non-Inflamed"),
                          levels = c("Non-Inflamed", "Inflamed"))


## Ward type
dat$ward_type = ifelse(dat$ward %in% c("ST 1", "ST 3", "ST 4"), "Ward", "Day-clinic") %>%
   factor(levels = c("Day-clinic", "Ward"))

## Ethnicity
dat$ethnicity = recode(dat$ethnicity,
                       '0' = "German",
                       '1' = "Other")



## Set binary variable indicating if cytokine data are complete (following imputation)
cyto_ref[, "na.sum"] = NA
for(i in cyto_ref$vars) {
      cyto_ref[cyto_ref$vars == i, "na.sum"] = sum(is.na(dat[, i]))
}
cyto_ref$no.na = ifelse(cyto_ref$na.sum >= 1, 0, 1)


## Somatic BDI variable
# Note: This is the sum of 4=lack of pleasure, 15=loss of energy, 16=sleeping problems, 
#       18=changes in appetite, 19=concentration difficulty, 20=tiredness or fatigue,
#       21=loss of interest in sex)
for(i in 0:7)  {
   dat[, paste0("t", i, "_bdi_som")] = rowSums(dat[, paste0("t", i, "_bdi_", c(4, 15, 16, 18:21))])
}


# 2 Preprocessing----------------------------------------

# 2.1 Reduce data to relevant variables------------------

## Data is reduced to relevant variables

dat = dat[, c("ID", "sex", "age", "ward_type", "ethnicity",
              "BMI", "hsCRP_inflamed", "t0_diagn_by_age",
              paste0("t0_cidi_", c("depression_episode", "depression_recurrent", 
                                   "anxiety", "substance")), 
              paste0("t", 0:7, "_bdi"), paste0("t", 0:7, "_bdi_som"), "t0_madrs",
              paste0("t0_pid_", c("negaff", "detach", "psycho", "antago", "disinh")),
              paste0("t0_bsi_", c("soma", "zwan", "unsi", "depr", "angs", "aggr", 
                                  "phob", "para", "psyc", "PSDI")),
              cyto_ref$vars)]


# 2.2 Interpolation of depressive symptom outcome--------

# 2.2.1 BDI----------------------------------------------

## Create locf variable
dat$t7_bdi_locf = NA

## Run loop to carry forward the last observed bdi-value
for(i in 1:nrow(dat))   {
      for(j in 7:0)   {
            if(is.na(dat[i, "t7_bdi_locf"]) & !is.na(dat[i, paste0("t", j, "_bdi")]))     {
                  dat[i, "t7_bdi_locf"] = dat[i, paste0("t", j, "_bdi")]
                  }
      }
}

## Create change variable by subtracting the t0 from the t7_locf score.
dat$bdi_locf_improve = dat$t0_bdi - dat$t7_bdi_locf


# 2.2.2 BDI Somatic--------------------------------------

dat$t7_bdi_som_locf = NA

## Run loop to carry forward the last observed bdi-value
for(i in 1:nrow(dat))   {
   for(j in 7:0)   {
      if(is.na(dat[i, "t7_bdi_som_locf"]) & !is.na(dat[i, paste0("t", j, "_bdi_som")]))     {
         dat[i, "t7_bdi_som_locf"] = dat[i, paste0("t", j, "_bdi_som")]
      }
   }
}

## Create change variable by subtracting the t0 from the t7_locf score.
dat$bdi_som_locf_improve = dat$t0_bdi_som - dat$t7_bdi_som_locf



# 3 Merge genetic data-----------------------------------

## Load corrected PRS data
load("./Data/PRScs/PRScs_corrected.RData")

## Define prs.vars
prs.vars = colnames(prs_dat)[grepl("_PRS", colnames(prs_dat), fixed = TRUE)]

## Merge data
dat = merge(dat, prs_dat[, c("IID", prs.vars)], by.x = "ID", by.y = "IID", 
            all.x = TRUE)


# 3 Merge transcriptomic data----------------------------

## Load data
load("./Data/TranscriptomicData_processed_PConly.RData")

## Rename df
rna = rna_proc

## Define rna.vars
rna.vars = colnames(rna)[!grepl("ID", colnames(rna))]

## Reduce RNA data to relevant individuals
rna = rna[rna$ID %in% dat$ID,]

## Remove duplicate IDs
rna = rna[!duplicated(rna$ID),]

## Merge data
dat = merge(dat, rna[, c("ID", rna.vars)], by = "ID", all.x = FALSE)



# 6 Save Data--------------------------------------------

## Save main data

save(dat, file = "./Data/OPTIMA_Cytokine_PreprocessedData.RData")


## Save info of inflammatory variables
save(cyto_ref, file = "./Data/Cytokine_Reference.RData")



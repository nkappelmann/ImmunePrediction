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

## Create binarised and standardised sex variable
dat$sex_std = ifelse(dat$sex == "female", 0.5, -0.5)

## Create standardised age variable
dat$age_std = scale(dat$age)

## Create standardised bmi variable
dat$BMI_std = scale(dat$BMI)

## Create standardised t0_bdi variable
dat$t0_bdi_std = scale(dat$t0_bdi)

## Get Top 27% CRP cut-off
dat$hsCRP_inflamed = factor(ifelse(dat$CRP > quantile(dat$hsCRP, prob=1-27/100), 
                                 "Inflamed", "Non-Inflamed"),
                          levels = c("Non-Inflamed", "Inflamed"))


## Ward
dat$ward_type = ifelse(grepl("ST", dat$ward, fixed = TRUE), "Ward", "Day-clinic")


## Other baseline variables

# Country of origin
dat$countryorigin = ifelse(dat$countryorigin == 0, "Germany", "Other")

# Ethnicity
dat$ethnicity = ifelse(dat$ethnicity == 0, "German", "Other")

# Current employment (code retired as extra category)
dat$employed = ifelse(dat$employed == 1, "Employed", "Unemployed")
dat[!is.na(dat$unemploymentreason) & dat$unemploymentreason == "retired", "employed"] = "Retired" 


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

# 2.1 Interpolation of depressive symptom outcome--------

# 2.1.1 BDI----------------------------------------------

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


# 2.1.2 BDI Somatic--------------------------------------

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



# 2.1.3 MADRS--------------------------------------------

## Create empty bdi_locf variable
dat$t7_madrs_locf = NA

## Run loop to carry forward the last observed bdi-value
for(i in 1:nrow(dat))   {
      for(j in c(7, 4, 0))   {
            if(is.na(dat[i, "t7_madrs_locf"]) & !is.na(dat[i, paste0("t", j, "_madrs")]))     {
                  dat[i, "t7_madrs_locf"] = dat[i, paste0("t", j, "_madrs")]
            }
      }
}

## Create change variable by subtracting the t0 from the t7_locf score.
dat$madrs_locf_improve = dat$t0_madrs - dat$t7_madrs_locf


# 3 Merge genetic data-----------------------------------

## Load corrected PRS data
load("/home/nkappelmann/OPTIMA/ImmunePrediction/Data/PRScs/PRScs_corrected.RData")

## Define prs.vars
prs.vars = paste0(list.files("/home/nkappelmann/PRScs/ImmunePrediction/"), 
                  "_PRS")


## Scale all PRSs
prs_dat[,prs.vars] = scale(prs_dat[,prs.vars])

## Merge data
dat = merge(dat, prs_dat[, c("IID", prs.vars)], by.x = "ID", by.y = "IID", 
            all.x = TRUE)


# 3 Merge transcriptomic data----------------------------

## Load data
load("/home/nkappelmann/OPTIMA/ImmunePrediction/Data/TranscriptomicData_processed_PConly.RData")

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



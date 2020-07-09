# -------------------------------------------------------
# Preparation--------------------------------------------
# -------------------------------------------------------


# 1 Preparation------------------------------------------

# 1.1 Setup----------------------------------------------

## Load Packages
library("tidyverse")
library("readxl")

## Load psychometric data
load("/home/nkappelmann/OPTIMA/OPTIMA_Analyses/Data/Processed/02_PsychometricData.RData")

## Load cytokine data and related IDs
cyto_dat = readRDS("/binder/mgp/datasets/2020_ImmuneDepression/cytokine/02_data_for_analysis.rds")
cyto_IDs = read_excel("/binder/mgp/datasets/2020_ImmuneDepression/cytokine/01_become_optima_randomization_STD_Controls_Pheno_1219.xlsx")

# Save biobank study IDs
cyto_dat$ID.bio = row.names(cyto_dat)

# Merge OPTIMA IDs
cyto_dat = merge(cyto_dat, cyto_IDs[, c("ID", "Lagerung_BC.Tube")], 
                 by.x = "ID.bio", by.y = "Lagerung_BC.Tube", all.x = TRUE)

# Delete duplicates from merging
cyto_dat = cyto_dat[!duplicated(cyto_dat$ID),]

# Save variable names in cytokine reference
cyto_ref = data.frame(vars = NA,
                      labels = colnames(cyto_dat)[9:ncol(cyto_dat)],
                      stringsAsFactors = FALSE)

# Change column names by excluding hyphens and slashes
colnames(cyto_dat) = gsub("/", "_", colnames(cyto_dat), fixed = TRUE)
colnames(cyto_dat) = gsub("-", ".", colnames(cyto_dat), fixed = TRUE)
colnames(cyto_dat) = gsub(" ", "_", colnames(cyto_dat), fixed = TRUE)
cyto_ref$vars = colnames(cyto_dat)[9:ncol(cyto_dat)]
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


# 1.2 Coding---------------------------------------------

## Create binarised and standardised sex variable
dat$sex_std = ifelse(dat$sex == "female", 0.5, -0.5)

## Create standardised age variable
dat$age_std = scale(dat$age)

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


# 2 Preprocessing----------------------------------------

# 2.1 LOCF of depressive symptom outcome-----------------

# 2.1.1 BDI----------------------------------------------

## Create empty bdi_locf variable
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


# 2.1.2 MADRS--------------------------------------------

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


# 3 Age, Sex & Non-normality Corrections-----------------

# 3.1 Log-transformation---------------------------------

# 3.2 Age & Sex correction-------------------------------

# 3.3 Scaling--------------------------------------------


# 4 Missing Data Imputation------------------------------

## Applying Nearest-Neighbour Imputation of Inflammatory Cytokines




# 5 Save Data--------------------------------------------

## Save main data
save(dat, file = "/binder/mgp/datasets/2020_ImmuneDepression/cytokine/Nils_Preprocessed/OPTIMA_Cytokine_PreprocessedData.RData")

## Save info of inflammatory variables
save(cyto_ref, file = "/binder/mgp/datasets/2020_ImmuneDepression/cytokine/Nils_Preprocessed/Cytokine_Reference.RData")

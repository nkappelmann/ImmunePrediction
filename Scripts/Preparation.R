# -------------------------------------------------------
# Preparation--------------------------------------------
# -------------------------------------------------------


# 1 Preparation------------------------------------------

# 1.1 Setup----------------------------------------------

## Load Packages
library("tidyverse")

## Load psychometric data
load("/home/nkappelmann/OPTIMA/OPTIMA_Analyses/Data/Processed/02_PsychometricData.RData")

## Load cytokine data
cyto_dat = readRDS("/binder/mgp/datasets/2020_ImmuneDepression/cytokine/01_data_for_analysis.rds")

# Save Study IDs
cyto_dat$ID = row.names(cyto_dat)

# Save variable names in cytokine reference
cyto_ref = data.frame(vars = NA,
                      labels = colnames(cyto_dat)[8:ncol(cyto_dat)])

# Change column names by excluding hyphens and slashes
colnames(cyto_dat) = gsub("/", "_", colnames(cyto_dat), fixed = TRUE)
colnames(cyto_dat) = gsub("-", ".", colnames(cyto_dat), fixed = TRUE)
colnames(cyto_dat) = gsub(" ", "_", colnames(cyto_dat), fixed = TRUE)
cyto_ref$vars = colnames(cyto_dat)[8:ncol(cyto_dat)]


## Merge psychometric and cytokine data
dat = merge(dat, cyto_dat, by = "ID", all.x = TRUE)


# 1.2 Coding---------------------------------------------


## CRP cut-off >3mg/L
dat$CRP_greater3 = factor(ifelse(dat$hsCRP > 3, ">3mg/L", "<=3mg/L"),
                          levels = c("<=3mg/L", ">3mg/L"))


## Sex as binary number and category
dat$sex_num = ifelse(dat$sex2 == "weiblich", 1, 0)
dat$sex_cat = ifelse(dat$sex2 == "weiblich", "Female", "Male")


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



# 2 Preprocessing----------------------------------------

# 2.1 LOCF of depressive symptom outcome-----------------

# 2.1.1 BDI----------------------------------------------

## Create empty bdi_locf variable
dat$t7_bdi_locf = NA

## Run loop to carry forward the last observed bdi-value
for(i in 1:nrow(dat))   {
      for(j in 7:1)   {
            if(is.na(dat[i, "t7_bdi_locf"]) & !is.na(dat[i, paste0("t", j, "_bdi")]))     {
                  dat[i, "t7_bdi_locf"] = dat[i, paste0("t", j, "_bdi")]
                  }
      }
}


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


# 3 Age, Sex & Non-normality Corrections-----------------

# 3.1 Log-transformation---------------------------------

# 3.2 Age & Sex correction-------------------------------

# 3.3 Scaling--------------------------------------------


# 4 Missing Data Imputation------------------------------

## Applying Nearest-Neighbour Imputation of Inflammatory Cytokines




# 5 Save Data--------------------------------------------

## Save main data
save(dat, file = "")

## Save info of inflammatory variables
save(cyto_ref, file = "")

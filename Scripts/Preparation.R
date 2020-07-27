# -------------------------------------------------------
# Preparation--------------------------------------------
# -------------------------------------------------------


# 1 Preparation------------------------------------------

# 1.1 Setup----------------------------------------------

## Load Packages
library("tidyverse")
library("readxl")
library("lme4")

## Load psychometric data
load("/home/nkappelmann/OPTIMA/OPTIMA_Analyses/Data/Processed/02_PsychometricData.RData")

## Load cytokine data and related IDs
cyto_dat = readRDS("/binder/mgp/datasets/2020_ImmuneDepression/cytokine/02_data_for_analysis_imputed.rds")
cyto_IDs = read_excel("/binder/mgp/datasets/2020_ImmuneDepression/cytokine/01_become_optima_randomization_STD_Controls_Pheno_1219.xlsx")

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

# 2.1 Interpolation of depressive symptom outcome--------

# 2.1.1 BDI----------------------------------------------

## Create longitudinal BDI data
dat_l = with(dat,
             data.frame(ID = rep(ID, 8),
                        time = rep(0:7, each = nrow(dat)),
                        bdi = c(t0_bdi, t1_bdi, t2_bdi, t3_bdi, t4_bdi, t5_bdi, t6_bdi, t7_bdi)))

## Obtain IDs with minimum of two BDI values
ids_with_bdi = dat[, paste0("t", 0:7, "_bdi")] %>% is.na() %>% rowSums < 7

## Create prediction model based on individuals with complete BDI data
model = lmer(bdi ~ time + I(time^2) + (1 + time | ID), 
             data = dat_l[rep(ids_with_bdi, 8),], REML = TRUE)

## Extract predicted BDI for individuals without t7_bdi
predict_dat = dat_l[dat_l$ID %in% dat[is.na(dat$t7_bdi), "ID"], ]



## Extract person-specific random intercept and random slope for time
model_ran_dat = ranef(model)$ID
model_ran_dat$ID = row.names(model_ran_dat)
colnames(model_ran_dat) = c("bdi_ran_int", "bdi_ran_slope", "ID")

dat = merge(dat, model_ran_dat, by = "ID", all.x = TRUE)




## Create empty bdi_intpol variable
dat$t7_bdi_intpol = NA

## Obtain rows of individuals with to-be-interpolated data
rows_to_intpol = which(is.na(dat$t7_bdi))


## Use the appox function for interpolation
for(i in rows_to_intpol)   {
   # Create subset of to-be-interpolated participant
   id_dat = dat[i, c("ID", paste0("t", 0:7, "_bdi"))]
   
   # Perform interpolation only if at least 2 datapoints are available
   if(sum(!is.na(id_dat[1, 2:9])) >= 2)  {
      
   # Use the appox function to interpolate bdi values
   intpol_dat = data.frame(approx(x = 1:8, y = as.numeric(id_dat[1, paste0("t", 0:7, "_bdi")]), 
                                  xout = 1:8, rule = 2,
                                  na.rm = TRUE, method = "linear"))
   
   # TEMP: Print result
   print(id_dat)
   print(t(intpol_dat))
   
   # Save value for T7
   dat[i, "t7_bdi_intpol"] = intpol_dat[8, "y"]
   }
}

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

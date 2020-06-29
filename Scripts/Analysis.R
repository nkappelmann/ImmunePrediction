# -------------------------------------------------------
# Analysis-----------------------------------------------
# -------------------------------------------------------


# 1 Preparation------------------------------------------

## Load packages
library("tidyverse")
library("RColorBrewer")
library("caret")
library("glmnet")
library("PubHelper")


## Load data
load("")

## Load cytokine reference
load("")



# 2 Descriptive statistics-------------------------------

# 2.1 Baseline Tables------------------------------------

## Baseline Table of clinical and sociodemographic characteristics
baselineTable(dat, round_dec = 2, placeholder = "   ",
              vars = c("sex_cat", "age", "socialclass_cat", 
                       "countryorigin", "ethnicity", "employed",
                       "t0_bdi", "t0_madrs", "t0_bsi_depr",
                       "ward_type",
                       "t0_cidi_diagnsum"),
              labels = c("Sex", "Age", "Social class",
                         "Country of Origin", "Ethnicity", "Employment status",
                         "BDI", "MADRS", "BSI: Depression",
                         "Hospitalisation in",
                         "Number of psychiatric diagnoses"))

## Baseline Tables of inflammatory parameters
baselineTable(dat, round_dec = 2, placeholder = "   ", grouping.var = "crp_greater3",
              vars = cyto_ref$vars,
              labels = cyto_ref$labels)


# 3 Model Training---------------------------------------

# 3.1 Elastic Net Regression-----------------------------


#trainControl(method = "LOOCV")
#trainControl(method = "boot")

# Zou and Hastie’s(2005) recommended default alpha = 0.5.

train(x = dat[,-y], y = dat[,y], 
      method = "glmnet", metric = "RMSE",
      trControl = trainControl(method = "LOOCV"))


# 3.2 Random Forest--------------------------------------

# Surpassing Variable importance threshold?





# 4 Content validity analysis-----------------------------


# 4.1 MADRS----------------------------------------------



# 4.2 BSI-Depression-------------------------------------


# 4.3 WHODAS---------------------------------------------


# 4.4 CIDI Depression------------------------------------


# 4.5 Dropout--------------------------------------------




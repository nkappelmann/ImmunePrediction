# -------------------------------------------------------
# Analysis-----------------------------------------------
# -------------------------------------------------------


# 1 Preparation------------------------------------------

## Load packages
library("tidyverse")
library("ggcorrplot")
library("RColorBrewer")
library("gghalves")
library("caret")
library("glmnet")
library("PubHelper")


## Load data
# Slurmgate
#load("/binder/mgp/datasets/2020_ImmuneDepression/cytokine/Nils_Preprocessed/OPTIMA_Cytokine_PreprocessedData.RData")

# MPI local
load("./Data/OPTIMA_Cytokine_PreprocessedData.RData")


# Index participants with minimum of two BDI observations
ids_with_bdi = dat[, paste0("t", 0:7, "_bdi")] %>% is.na() %>% rowSums < 7

## Load cytokine reference
# Slurmgate
#load("/binder/mgp/datasets/2020_ImmuneDepression/cytokine/Nils_Preprocessed/Cytokine_Reference.RData")

# MPI local
load("./Data/Cytokine_Reference.RData")


## Source nested cross-validation function
source("./Scripts/functions.R")




# 2 Descriptive statistics-------------------------------

# 2.1 Baseline Tables------------------------------------

## Baseline Table of clinical and sociodemographic characteristics
Table1 = baselineTable(dat, round_dec = 2, placeholder = "   ", grouping.var = "hsCRP_inflamed",
                       vars = c("sex", "age", "BMI", "socialclass_cat", 
                                "countryorigin", "ethnicity", "employed",
                                "t0_bdi", "t0_madrs", "t0_bsi_depr",
                                "ward_type",
                                "t0_cidi_diagnsum"),
                       labels = c("Sex", "Age", "BMI", "Social class",
                                  "Country of Origin", "Ethnicity", "Employment status",
                                  "BDI", "MADRS", "BSI: Depression",
                                  "Hospitalisation in",
                                  "Number of psychiatric diagnoses"))

## Baseline Tables of inflammatory parameters
TableS1 = baselineTable(dat, round_dec = 2, placeholder = "   ", 
                        grouping.var = "hsCRP_inflamed",
                        vars = cyto_ref$vars,
                        labels = cyto_ref$labels)

exportPubHelpercsv(Table1, file = "./Results/Table1.csv")
exportPubHelpercsv(TableS1, file = "./Results/TableS1.csv")


# 3 ML Analysis Pipeline---------------------------------

# 3.1 Baseline Covariates--------------------------------

## Define parameters
x = c("t0_bdi_std", "sex_std", "age_std", "BMI_std")

## Run analysis
set.seed(8)
covariates.base.output = nested.cv(data = dat[ids_with_bdi,], 
                                 x = x,
                                 y = "t7_bdi_locf",
                                 k.outer = 5, 
                                 k.inner = 5, 
                                 num_repeats = 100, 
                                 perm.test = FALSE,
                                 runGLMnet = TRUE,
                                 runRF = TRUE,
                                 runSVM = TRUE,
                                 runBART = FALSE)

# Save results
save(covariates.base.output, file = "./Results/covariates.base.output.RData")

## Run permutation analysis
set.seed(9)
covariates.base.perm = nested.cv(data = dat[ids_with_bdi,], 
                               x = x,
                               y = "t7_bdi_locf",
                               k.outer = 5, 
                               k.inner = 5, 
                               num_repeats = 100, 
                               perm.test = TRUE,
                               runGLMnet = TRUE,
                               runRF = TRUE,
                               runBART = FALSE)

save(covariates.base.perm, file = "./Results/covariates.base.perm.RData")




# 3.2 Cytokines------------------------------------------


# 3.2.1 Without Covariates-------------------------------

## Define parameters
x = cyto_ref$vars

## Run analysis
set.seed(10)
cytokine.base.output = nested.cv(data = dat[ids_with_bdi,], 
                                 x = x,
                                 y = "t7_bdi_locf",
                                 k.outer = 5, 
                                 k.inner = 5, 
                                 num_repeats = 100, 
                                 perm.test = FALSE,
                                 runGLMnet = TRUE,
                                 runRF = TRUE,
                                 runSVM = TRUE,
                                 runBART = FALSE)

# Save results
save(cytokine.base.output, file = "./Results/cytokine.base.output.RData")

## Run permutation analysis
set.seed(11)
cytokine.base.perm = nested.cv(data = dat[ids_with_bdi,], 
                               x = x,
                               y = "t7_bdi_locf",
                               k.outer = 5, 
                               k.inner = 5, 
                               num_repeats = 100, 
                               perm.test = TRUE,
                               runGLMnet = TRUE,
                               runRF = TRUE,
                               runBART = FALSE)

save(cytokine.base.perm, file = "./Results/cytokine.base.perm.RData")


## Summary statistics
with(cytokine.base.output$fit, by(RMSE, model, summary))
with(ML_results$fit[ML_results$fit$model == "glmnet",], hist(RMSE))

with(ML_perm$fit, by(RMSE, model, summary))


t.test(ML_perm$fit[ML_perm$fit$model == "glmnet", "RMSE"], 
       ML_results$fit[ML_results$fit$model == "glmnet", "RMSE"], var.equal = TRUE)
t.test(ML_perm$fit[ML_perm$fit$model == "glmnet", "Rsquared"], 
       ML_results$fit[ML_results$fit$model == "glmnet", "Rsquared"], var.equal = TRUE)
t.test(ML_perm$fit[ML_perm$fit$model == "rf", "RMSE"], 
       ML_results$fit[ML_results$fit$model == "rf", "RMSE"], var.equal = TRUE)
t.test(ML_perm$fit[ML_perm$fit$model == "rf", "Rsquared"], 
       ML_results$fit[ML_results$fit$model == "rf", "Rsquared"], var.equal = TRUE)





# 3.2.2 With Covariates----------------------------------


## Define parameters
x = c(cyto_ref$vars, "t0_bdi_std", "sex_std", "age_std", "BMI_std")

## Run analysis
set.seed(12)
cytokine.comb.output = nested.cv(data = dat[ids_with_bdi,], 
                                 x = x,
                                 y = "t7_bdi_locf",
                                 k.outer = 5, 
                                 k.inner = 5, 
                                 num_repeats = 100, 
                                 perm.test = FALSE,
                                 runGLMnet = TRUE,
                                 runRF = TRUE,
                                 runSVM = TRUE,
                                 runBART = FALSE)

# Save results
save(cytokine.comb.output, file = "./Results/cytokine.comb.output.RData")

## Run permutation analysis
set.seed(13)
cytokine.comb.perm = nested.cv(data = dat[ids_with_bdi,], 
                               x = x,
                               y = "t7_bdi_locf",
                               k.outer = 5, 
                               k.inner = 5, 
                               num_repeats = 100, 
                               perm.test = TRUE,
                               runGLMnet = TRUE,
                               runRF = TRUE,
                               runBART = FALSE)

save(cytokine.comb.perm, file = "./Results/cytokine.comb.perm.RData")





# 3.3 Gene-expression------------------------------------


# 3.2.1 Without Covariates-------------------------------


# 3.2.2 With Covariates----------------------------------



# 3.4 Polygenic Risk Scores------------------------------


# 3.2.1 Without Covariates-------------------------------


# 3.2.2 With Covariates----------------------------------



# 3.5 Combined Immunophenotyping-------------------------


# 3.2.1 Without Covariates-------------------------------


# 3.2.2 With Covariates----------------------------------



# 4 Model Evaluation-------------------------------------

## Create evaluation data.frame
model.comparison = expand.grid(model = c("Covariates only", "Cytokines", "Gene-expression", 
                                         "Polygenic risk scores", "Combined immunophenotyping"),
                               covariates = c("Without covariates", "With covariates"))
model.comparison$tval = NA
model.comparison$pval = NA


# 4.1 With Covariates------------------------------------

## Covariate model



## Cytokine model
t.test.output = t.test(cytokine.base.output$fit$RMSE, cytokine.base.perm$fit$RMSE)
model.comparison[model.comparison$model == "Cytokines" & 
                    model.comparison$covariates == "Without covariates", c("tval", "pval")] = 
c(t.test.output$statistic, t.test.output$p.value)


## Gene-expression model


## Polygenic risk score model



# 4.2 Without Covariates---------------------------------

## Covariate model



## Cytokine model


## Gene-expression model


## Polygenic risk score model



# 4 Content validity analysis----------------------------


# 4.1 MADRS----------------------------------------------



# 4.2 BSI-Depression-------------------------------------


# 4.3 WHODAS---------------------------------------------


# 4.4 CIDI Depression------------------------------------


# 4.5 Dropout--------------------------------------------


# 5 Visualisation----------------------------------------


# 5.1 Preparation----------------------------------------

## Collate fit statistics
# Add meta-data
covariates.base.output$fit$pred = "Covariates"
covariates.base.output$fit$pred.type = "Model Prediction"
covariates.base.output$fit$covariates = "With covariates"
covariates.base.perm$fit$pred = "Covariates"
covariates.base.perm$fit$pred.type = "Null Model"
covariates.base.perm$fit$covariates = "With covariates"

cytokine.base.output$fit$pred = "Cytokines"
cytokine.base.output$fit$pred.type = "Model Prediction"
cytokine.base.output$fit$covariates = "Without covariates"
cytokine.base.perm$fit$pred = "Cytokines"
cytokine.base.perm$fit$pred.type = "Null Model"
cytokine.base.perm$fit$covariates = "Without covariates"

cytokine.comb.output$fit$pred = "Cytokines"
cytokine.comb.output$fit$pred.type = "Model Prediction"
cytokine.comb.output$fit$covariates = "With covariates"
cytokine.comb.perm$fit$pred = "Cytokines"
cytokine.comb.perm$fit$pred.type = "Null Model"
cytokine.comb.perm$fit$covariates = "With covariates"

# Rowbind data
fit.stats = rbind.data.frame(covariates.base.output$fit,
                             covariates.base.perm$fit,
                             cytokine.base.output$fit,
                             cytokine.base.perm$fit,
                             cytokine.comb.output$fit,
                             cytokine.comb.perm$fit)

## Recode model
fit.stats$algorithm = recode(fit.stats$model,
                             'glmnet' = "Elastic net regression",
                             'rf' = "Random forest",
                             'bart' = "BARTmachine")

## Set factor levels
fit.stats$covariates = factor(fit.stats$covariates, 
                              levels = c("Without covariates", "With covariates"))

vis_dat = with(dat, 
               data.frame(t7_bdi_locf = rep(t7_bdi_locf, 6),
                          t7_bdi_pred = c(pred.base.enet, pred.unadj.enet, pred.adj.enet,
                                          pred.base.rf, pred.unadj.rf, pred.adj.rf),
                          method = c(rep("Elastic Net", 3*nrow(dat)), 
                                     rep("Random Forest", 3*nrow(dat))),
                          model = rep(rep(c("Baseline", "Cytokines", "Baseline+Cytokines"), 
                                          each = nrow(dat)), 2),
                          r2 = rep(c(summary(lm(t7_bdi_locf ~ pred.base.enet))$r.squared, 
                                     summary(lm(t7_bdi_locf ~ pred.unadj.enet))$r.squared, 
                                     summary(lm(t7_bdi_locf ~ pred.adj.enet))$r.squared,
                                     summary(lm(t7_bdi_locf ~ pred.base.rf))$r.squared, 
                                     summary(lm(t7_bdi_locf ~ pred.unadj.rf))$r.squared, 
                                     summary(lm(t7_bdi_locf ~ pred.adj.rf))$r.squared),
                                   each = nrow(dat))
               )
)

## Set factor leves
vis_dat$model = factor(vis_dat$model, levels = c("Baseline", "Cytokines", "Baseline+Cytokines"))



# 5.2 Cytokine Descriptive Statistics--------------------

## Get complete data proportion per cytokine
cyto_ref$perc.complete = (nrow(dat) - cyto_ref$na.sum) / nrow(dat) * 100

## Get complete variable
cyto_ref$complete = ifelse(cyto_ref$no.na == 1, "Complete", "Not Complete")

## Create barplot
ggplot(cyto_ref, aes(x = reorder(labels, perc.complete), y = perc.complete, fill = complete)) +
   geom_bar(aes(fill = complete), col = "black", stat = "identity") +
   scale_fill_brewer(palette = "Dark2") +
   labs(x = "", y = "Complete Data (%)") +
   coord_flip() +
   theme_bw() +
   theme(legend.position = "none",
         legend.title = element_blank())



# 5.3 Cytokine Correlations------------------------------


## Save correlation matrix and p-value matrix
cyto_cor = cor(dat[, cyto_ref$vars], use = "pairwise.complete.obs", method = "pearson")
cyto_cor_pmat = cor_pmat(dat[, cyto_ref$vars])

## Set labels 
colnames(cyto_cor) = cyto_ref$labels
rownames(cyto_cor) = cyto_ref$labels
colnames(cyto_cor_pmat) = cyto_ref$labels
rownames(cyto_cor_pmat) = cyto_ref$labels

# Reduce to those cytokines with no missing values (i.e., used in analyses)
cyto_cor = cyto_cor[row.names(cyto_cor) %in% cyto_ref[cyto_ref$no.na == 1, "labels"],
                    colnames(cyto_cor) %in% cyto_ref[cyto_ref$no.na == 1, "labels"]]
cyto_cor_pmat = cyto_cor_pmat[row.names(cyto_cor_pmat) %in% cyto_ref[cyto_ref$no.na == 1, "labels"],
                              row.names(cyto_cor_pmat) %in% cyto_ref[cyto_ref$no.na == 1, "labels"]]

## Export with width=900; height=612
ggcorrplot(cyto_cor, 
           hc.order = TRUE, type = "lower", 
           legend.title = "", show.legend = TRUE,
           #p.mat = cyto_cor_pmat,
           outline.col = "white", lab = TRUE, digits = 1, lab_size = 3,
           #ggtheme = ggplot2::theme_minimal(),
           colors = c(brewer.pal(3, "Dark2")[1], "white", brewer.pal(3, "Dark2")[2])) 



# 5.4 ML Prediction Scatter Plots------------------------

## Prediction (export with: width = 800, height = standard)
ggplot(vis_dat, aes(x = t7_bdi_pred, y = t7_bdi_locf)) +
      geom_smooth(aes(col = model, fill = model), method = "lm") +
      geom_point(aes(fill = model), col = "black", shape = 21, size = 2.5) +
      facet_grid(model~method) +
      #geom_text(x = 10, y = 40, label = expression("Spearman's"~rho == 0.34), parse = TRUE) +
      labs(x = "BDI: Predicted Symptom Score", y = "BDI: Actual Symptom Score",
           fill = "Model: ", col = "Model: ") +
      scale_fill_brewer(palette = "Dark2") +
      scale_color_brewer(palette = "Dark2") +
      theme_bw() +
      theme(legend.position = "top",
            #strip.text.x = element_text(colour = "black"),
            #strip.background = element_rect(colour = "black", fill = "white")
            )



# 5.5 Absolute Residual Boxplots-------------------------

## Calculate absolute residuals
vis_dat$abs_residuals = with(vis_dat, abs(t7_bdi_locf - t7_bdi_pred))

## Plot
ggplot(vis_dat, aes(x = method, y = abs_residuals)) +
      geom_boxplot(aes(fill = model), alpha = 0.6, position = position_dodge(width = 0.8),
                   outlier.alpha = 0) +
      geom_jitter(aes(fill = model), shape = 21, col = "black", 
                  position = position_jitterdodge(dodge.width = 0.8)) +
      scale_fill_brewer(palette = "Dark2") +
      labs(x = "Model", y = "Absolute Residuals") +
      theme_bw() +
      theme(legend.position = "top",
            legend.title = element_blank())



# 5.6 Fit Statistics-------------------------------------

## RMSE (width=700; height=473)
ggplot(fit.stats[fit.stats$model != "bart",], aes(x = pred, y = RMSE)) +
   geom_half_violin(aes(fill = pred.type), col = "black", alpha = 0.8) +
   geom_half_boxplot(aes(fill = pred.type), outlier.alpha = 0, alpha = 0.8, col = "black",
                     errorbar.draw = TRUE, side = "r") +
   geom_half_point(aes(fill = pred.type, col = pred.type), size = 0.2, alpha = 0.3) +
   facet_grid(algorithm~covariates, scales = "free_x", space = "free_x") +
   scale_fill_brewer(palette = "Dark2") +
   scale_color_brewer(palette = "Dark2") +
   #scale_y_continuous(limits = c(0, 20)) +
   coord_cartesian(expand = FALSE) +
   labs(x = "", y = "RMSE") +
   theme_bw() +
   theme(legend.position = "top",
         legend.title = element_blank())



## R2  (width=700; height=473)
ggplot(fit.stats[fit.stats$model != "bart",], aes(x = pred, y = Rsquared)) +
   geom_half_violin(aes(fill = pred.type), col = "black", alpha = 0.8) +
   geom_half_boxplot(aes(fill = pred.type), outlier.alpha = 0, alpha = 0.8, col = "black",
                     errorbar.draw = TRUE, side = "r") +
   geom_half_point(aes(fill = pred.type, col = pred.type), size = 0.2, alpha = 0.3) +
   facet_grid(algorithm~covariates, scales = "free_x", space = "free_x") +
   scale_fill_brewer(palette = "Dark2") +
   scale_color_brewer(palette = "Dark2") +
   #scale_y_continuous(limits = c(0, 0.5)) +
   coord_cartesian(expand = FALSE) +
   labs(x = "", y = expression(R^2)) +
   theme_bw() +
   theme(legend.position = "top",
         legend.title = element_blank())




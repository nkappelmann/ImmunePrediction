# -------------------------------------------------------
# Evaluation---------------------------------------------
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

## Set Slurm directory
setwd("/home/nkappelmann/OPTIMA/ImmunePrediction")


## Load data
load("./Data/OPTIMA_Cytokine_PreprocessedData.RData")


## Load ML model output
#load("./Results/covariates.output.RData")
load("./Results/clinical.output.RData")
load("./Results/cytokine.output.RData")
load("./Results/rna.output.RData")
load("./Results/prs.output.RData")
load("./Results/omics.output.RData")
load("./Results/omicsplusclin.output.RData")


## Index participants used in analysis
# Index participants with minimum of two BDI observations
ids_with_bdi = dat[, paste0("t", 0:7, "_bdi")] %>% is.na() %>% rowSums < 7

# Index participants with genetic data
ids_with_geneticdata = !is.na(dat$BMI_PRS)

# Index participants with transcriptomic data
ids_with_rnadata = !is.na(dat$rna.PC1) # Present in all individuals

# Index participants with all relevant data
ids_for_analysis = ids_with_bdi & ids_with_geneticdata & ids_with_rnadata

## Load cytokine reference
load("./Data/Cytokine_Reference.RData")

## Define vars
prs.vars = colnames(dat)[grepl("PRS", colnames(dat))]
rna.vars = colnames(dat)[grepl("rna.PC", colnames(dat))]
clin.vars = c("t0_bdi_std", "t0_madrs_std", "sex_std", "age_std", "BMI_std", "t0_diagn_by_age_std",
              paste0("t0_bsi_", c("soma", "zwan", "unsi", "depr", "angs", "aggr", 
                                  "phob", "para", "psyc"), "_std"))


## Source nested cross-validation function
source("./Scripts/functions.R")


# 2 Preparation------------------------------------------

# 2.1 Fit statistics-------------------------------------

## Fit statistics are combined
fit.stats = rbind.data.frame(clinical.output$fit,
                             cytokine.output$fit,
                             prs.output$fit,
                             rna.output$fit,
                             omics.output$fit,
                             omicsplusclin.output$fit)

## Meta-data is added
fit.stats$feature.set = rep(c("Clinical", "Cytokines", "PRS", 
                          "Gene expression", "Multi-omics",
                          "Multi-omics\n+\nClinical"), each = nrow(clinical.output$fit))




# 3 Descriptive statistics-------------------------------

# 3.1 Baseline Tables------------------------------------

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





# 4 Model Evaluation-------------------------------------

# 4.1 Paired t-tests-------------------------------------


## Conduct paired t-test by group
t.test.results = t.test_loop(fit.stats)

## Recode model
t.test.results$algorithm = recode(t.test.results$model,
                        'glmnet' = "Elastic net regression",
                        'rf' = "Random forest",
                        'knn' = "k-nearest neighbour")


## Plot results
# RMSE_delta mean + sd
ggplot(t.test.results, aes(x = feature.set, y = RMSE_delta_mean)) +
   geom_hline(yintercept = 0, col = "lightgrey", size = 2) +
   geom_errorbar(aes(ymin = RMSE_delta_mean - RMSE_delta_sd, 
                     ymax = RMSE_delta_mean + RMSE_delta_sd, group = model), 
                 stat = "identity", position = position_dodge(width = 0.3), width = 0.2) +
   geom_point(aes(fill = algorithm), position = position_dodge(width = 0.3), shape = 21, 
              col = "black", size = 4) +
   labs(x = "", y = bquote(Delta*"RMSE: Median (95% CI)")) +
   scale_fill_brewer(palette = "Dark2") +
   theme_bw() +
   theme(legend.position = "top",
         legend.title = element_blank())

# RMSE_delta median 95% quantile
ggplot(t.test.results, aes(x = feature.set, y = RMSE_delta_median)) +
   geom_hline(yintercept = 0, col = "lightgrey", size = 2) +
   geom_errorbar(aes(ymin = RMSE_delta_5quantile, ymax = RMSE_delta_95quantile, group = model), 
                 stat = "identity", position = position_dodge(width = 0.4), width = 0.2) +
   geom_point(aes(fill = algorithm), position = position_dodge(width = 0.4), shape = 21, 
              col = "black", size = 3) +
   labs(x = "", y = bquote(Delta*"RMSE: Median (95% CI)")) +
   scale_fill_brewer(palette = "Dark2") +
   theme_bw() +
   theme(legend.position = "top",
         legend.title = element_blank(),
         plot.background = element_rect(fill = "transparent", color = NA),
         legend.background = element_rect(fill = "transparent"))
ggsave(filename = paste0("./Plots/", Sys.Date(), "_PredPerformance_RMSE.png"), 
       device = "png", width = 6, height = 4, dpi = 300, bg = "transparent")



# 4.2 XY-------------------------------------------------

## Create evaluation data.frame
model.comparison = expand.grid(model = c("Covariates only", "Cytokines", "Gene-expression", 
                                         "Polygenic risk scores", "Combined immunophenotyping"),
                               covariates = c("Without covariates", "With covariates"))
model.comparison$tval = NA
model.comparison$pval = NA

# Remove non-sensival row
model.comparison = model.comparison[2:nrow(model.comparison),]


# 4.1 With Covariates------------------------------------

## Covariate model
t.test.output = t.test(covariates.base.output$fit$RMSE, covariates.base.perm$fit$RMSE)
model.comparison[model.comparison$model == "Covariates only" & 
                       model.comparison$covariates == "With covariates", c("tval", "pval")] = 
      c(t.test.output$statistic, t.test.output$p.value)


## Cytokine model
# Without covariates
t.test.output = t.test(cytokine.base.output$fit$RMSE, cytokine.base.perm$fit$RMSE)
model.comparison[model.comparison$model == "Cytokines" & 
                       model.comparison$covariates == "Without covariates", c("tval", "pval")] = 
      c(t.test.output$statistic, t.test.output$p.value)

# With covariates
t.test.output = t.test(cytokine.comb.output$fit$RMSE, cytokine.comb.perm$fit$RMSE)
model.comparison[model.comparison$model == "Cytokines" & 
                       model.comparison$covariates == "With covariates", c("tval", "pval")] = 
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





# 5.1.1 Fit statistics-----------------------------------

## Collate fit statistics
# Add meta-data
covariates.base.output$fit[, c("pred", "pred.type", "covariates")] =
      rep(c("Clinical", "Model Prediction", "With covariates"),
          each = nrow(covariates.base.output$fit))
covariates.base.perm$fit[, c("pred", "pred.type", "covariates")] =
      rep(c("Clinical", "Null Model", "With covariates"),
          each = nrow(covariates.base.perm$fit))

cytokine.base.output$fit[, c("pred", "pred.type", "covariates")] =
      rep(c("Cytokines", "Model Prediction", "Without covariates"),
          each = nrow(cytokine.base.output$fit))
cytokine.base.perm$fit[, c("pred", "pred.type", "covariates")] =
      rep(c("Cytokines", "Null Model", "Without covariates"),
          each = nrow(cytokine.base.perm$fit))

cytokine.comb.output$fit[, c("pred", "pred.type", "covariates")] =
      rep(c("Cytokines", "Model Prediction", "With covariates"),
          each = nrow(cytokine.comb.output$fit))
cytokine.comb.perm$fit[, c("pred", "pred.type", "covariates")] =
      rep(c("Cytokines", "Null Model", "With covariates"),
          each = nrow(cytokine.comb.perm$fit))

prs.base.output$fit[, c("pred", "pred.type", "covariates")] =
      rep(c("PRS", "Model Prediction", "Without covariates"),
          each = nrow(prs.base.output$fit))
prs.base.perm$fit[, c("pred", "pred.type", "covariates")] =
      rep(c("PRS", "Null Model", "Without covariates"),
          each = nrow(prs.base.perm$fit))

prs.comb.output$fit[, c("pred", "pred.type", "covariates")] =
      rep(c("PRS", "Model Prediction", "With covariates"),
          each = nrow(prs.base.output$fit))
prs.comb.perm$fit[, c("pred", "pred.type", "covariates")] =
      rep(c("PRS", "Null Model", "With covariates"),
          each = nrow(prs.base.perm$fit))

rna.base.output$fit[, c("pred", "pred.type", "covariates")] =
      rep(c("Gene expression", "Model Prediction", "Without covariates"),
          each = nrow(prs.base.output$fit))
rna.base.perm$fit[, c("pred", "pred.type", "covariates")] =
      rep(c("Gene expression", "Null Model", "Without covariates"),
          each = nrow(prs.base.perm$fit))

rna.comb.output$fit[, c("pred", "pred.type", "covariates")] =
      rep(c("Gene expression", "Model Prediction", "With covariates"),
          each = nrow(prs.base.output$fit))
rna.comb.perm$fit[, c("pred", "pred.type", "covariates")] =
      rep(c("Gene expression", "Null Model", "With covariates"),
          each = nrow(prs.base.perm$fit))

omics.base.output$fit[, c("pred", "pred.type", "covariates")] =
      rep(c("Multi-omics", "Model Prediction", "Without covariates"), 
          each = nrow(omics.base.output$fit))
omics.base.perm$fit[, c("pred", "pred.type", "covariates")] =
      rep(c("Multi-omics", "Null Model", "Without covariates"), 
          each = nrow(omics.base.perm$fit))

omics.comb.output$fit[, c("pred", "pred.type", "covariates")] =
      rep(c("Multi-omics\n+\nClinical", "Model Prediction", "With covariates"), 
          each = nrow(omics.comb.output$fit))
omics.comb.perm$fit[, c("pred", "pred.type", "covariates")] =
      rep(c("Multi-omics\n+\nClinical", "Null Model", "With covariates"), 
          each = nrow(omics.comb.perm$fit))





# 5.1.2 Variable importance------------------------------

## Transform variable importance data to long data.frame
varImp_long = omicsplusclin.output$varImp %>%
   pivot_longer(cols = age_std:VEGF.D, names_to = "vars", values_to = "varImp") %>%
   filter(type == "pred") %>%
   group_by(vars) %>% 
   mutate(mean_varImp = mean(varImp)) %>%
   ungroup() %>%
   arrange(desc(mean_varImp))

# Define algorithm
varImp_long$algorithm = varImp_long$model %>%
   recode('glmnet' = "Elastic net regression",
          'rf' = "Random forest",
          'knn' = "k-nearest neighbour")




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


# 5.4 PRS Correlations-----------------------------------


## Save correlation matrix and p-value matrix
prs_cor = cor(dat[, prs.vars], use = "pairwise.complete.obs", method = "pearson")
prs_cor_pmat = cor_pmat(dat[, prs.vars])

## Set labels 
colnames(prs_cor) = gsub("_", " ", prs.vars)
rownames(prs_cor) = gsub("_", " ", prs.vars)
colnames(prs_cor_pmat) = gsub("_", " ", prs.vars)
rownames(prs_cor_pmat) = gsub("_", " ", prs.vars)

## Export with width=1000; height=700
ggcorrplot(prs_cor, 
           hc.order = TRUE, type = "lower", 
           legend.title = "", show.legend = TRUE,
           #p.mat = cyto_cor_pmat,
           outline.col = "white", lab = TRUE, digits = 1, lab_size = 3,
           #ggtheme = ggplot2::theme_minimal(),
           colors = c(brewer.pal(3, "Dark2")[1], "white", brewer.pal(3, "Dark2")[2])) 



# 5.5 ML Prediction Scatter Plots------------------------

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



# 5.6 Absolute Residual Boxplots-------------------------

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



# 5.7 Fit Statistics-------------------------------------

## RMSE
ggplot(fit.stats, aes(x = pred, y = RMSE)) +
      geom_half_violin(aes(fill = pred.type), col = "black", alpha = 0.8) +
      geom_half_boxplot(aes(fill = pred.type), outlier.alpha = 0, alpha = 0.8, col = "black",
                        errorbar.draw = TRUE, side = "r") +
      geom_half_point(aes(fill = pred.type, col = pred.type), size = 0.2, alpha = 0.3) +
      facet_grid(algorithm~., scales = "free_x", space = "free_x") +
      scale_fill_brewer(palette = "Dark2") +
      scale_color_brewer(palette = "Dark2") +
      #scale_y_continuous(limits = c(0, 20)) +
      coord_cartesian(expand = FALSE) +
      labs(x = "", y = "RMSE") +
      theme_bw() +
      theme(legend.position = "top",
            legend.title = element_blank(),
            plot.background = element_rect(fill = "transparent", colour = NA),
            panel.background = element_rect(fill = "transparent", colour = NA),
            legend.background = element_rect(fill = "transparent", colour = NA))
ggsave(filename = paste0("./Plots/", Sys.Date(), "_PredPerformance_RMSE.png"), 
       bg = "transparent", dpi = 600, width = 6, height = 5)


## R2  (width=750; height=507)
ggplot(fit.stats, aes(x = pred, y = Rsquared)) +
      geom_half_violin(aes(fill = pred.type), col = "black", alpha = 0.8) +
      geom_half_boxplot(aes(fill = pred.type), outlier.alpha = 0, alpha = 0.8, col = "black",
                        errorbar.draw = TRUE, side = "r") +
      geom_half_point(aes(fill = pred.type, col = pred.type), size = 0.2, alpha = 0.3) +
      facet_grid(algorithm~., scales = "free_x", space = "free_x") +
      scale_fill_brewer(palette = "Dark2") +
      scale_color_brewer(palette = "Dark2") +
      #scale_y_continuous(limits = c(0, 0.5)) +
      coord_cartesian(expand = FALSE) +
      labs(x = "", y = expression(R^2)) +
      theme_bw() +
      theme(legend.position = "top",
            legend.title = element_blank(),
            plot.background = element_rect(fill = "transparent", colour = NA),
            panel.background = element_rect(fill = "transparent", colour = NA),
            legend.background = element_rect(fill = "transparent", colour = NA))
ggsave(filename = paste0("./Plots/", Sys.Date(), "_PredPerformance_R2.png"), 
       bg = "transparent", dpi = 600, width = 6, height = 5)



## Bar charts
# Get summary data
fit.stats.summary = fit.stats %>%
      group_by(algorithm, pred.type, pred) %>%
      summarise(RMSE_conf.low = quantile(RMSE, probs = 0.025),
                RMSE_conf.high = quantile(RMSE, probs = 0.975),
                Rsquared_conf.low = quantile(Rsquared, probs = 0.025, na.rm = TRUE),
                Rsquared_conf.high = quantile(Rsquared, probs = 0.975, na.rm = TRUE),
                RMSE_se = sqrt(var(RMSE) * (1/5 + (98*(1/5))/(98*(4/5)))),
                Rsquared_se = sqrt(var(Rsquared, na.rm = TRUE) * 
                                         (1/5 + (98*(1/5))/(98*(4/5)))),
                RMSE = mean(RMSE), 
                Rsquared = mean(Rsquared, na.rm = TRUE),
                Rsquared.median = median(Rsquared, na.rm = TRUE)
      )

# RMSE
ggplot(fit.stats.summary, aes(x = pred, y = RMSE, group = pred.type)) +
      geom_bar(aes(fill = pred.type), col = "black", alpha = 0.8,
               position = position_dodge(width = 1), stat = "identity") +
      geom_errorbar(aes(group = pred.type, ymin = RMSE - RMSE_se, ymax = RMSE + RMSE_se), 
                    col = "black", width = 0.2,
                    position = position_dodge(width = 1), stat = "identity") +
      facet_grid(algorithm~., scales = "free_x", space = "free_x") +
      scale_fill_brewer(palette = "Dark2") +
      scale_color_brewer(palette = "Dark2") +
      #scale_y_continuous(limits = c(0, 20)) +
      coord_cartesian(ylim = c(10, 18)) +
      labs(x = "", y = "RMSE"
           #y = "RMSE (Robust SE)"
      ) +
      theme_bw() +
      theme(legend.position = "top",
            legend.title = element_blank(),
            plot.background = element_rect(fill = "transparent", colour = NA),
            panel.background = element_rect(fill = "transparent", colour = NA),
            legend.background = element_rect(fill = "transparent", colour = NA))
ggsave(filename = paste0("./Plots/", Sys.Date(), "_PredPerformance_RMSE_bar.png"), 
       bg = "transparent", dpi = 600, width = 6, height = 5)

# R2
ggplot(fit.stats.summary, aes(x = pred, y = Rsquared.median, group = pred.type)) +
      geom_bar(aes(fill = pred.type), col = "black", alpha = 0.8,
               position = position_dodge(width = 1), stat = "identity") +
      geom_errorbar(aes(group = pred.type, ymin = Rsquared - Rsquared_se, 
                        ymax = Rsquared + Rsquared_se), col = "black", width = 0.2,
                    position = position_dodge(width = 1), stat = "identity") +
      facet_grid(algorithm~., scales = "free_x", space = "free_x") +
      scale_fill_brewer(palette = "Dark2") +
      scale_color_brewer(palette = "Dark2") +
      #scale_y_continuous(limits = c(0, 20)) +
      labs(x = "", y = expression(paste("Median ", R^2))
           #y = expression(paste(R^2, " (Robust SE)"))
      ) +
      theme_bw() +
      theme(legend.position = "top",
            legend.title = element_blank(),
            plot.background = element_rect(fill = "transparent", colour = NA),
            panel.background = element_rect(fill = "transparent", colour = NA),
            legend.background = element_rect(fill = "transparent", colour = NA))
ggsave(filename = paste0("./Plots/", Sys.Date(), "_PredPerformance_R2_bar.png"), 
       bg = "transparent", dpi = 600, width = 6, height = 5)



## Comparison test
model = lm(RMSE ~ pred.type + algorithm, 
           data = fit.stats[fit.stats$pred == "Clinical",])
summary(model)
model = lm(RMSE ~ pred.type + algorithm, 
           data = fit.stats[fit.stats$pred == "Cytokines",])
summary(model)
model = lm(RMSE ~ pred.type + algorithm, 
           data = fit.stats[fit.stats$pred == "PRS",])
summary(model)
model = lm(RMSE ~ pred.type + algorithm, 
           data = fit.stats[fit.stats$pred == "Gene expression",])
summary(model)
model = lm(RMSE ~ pred.type + algorithm, 
           data = fit.stats[fit.stats$pred == "Multi-omics",])
summary(model)
model = lm(RMSE ~ pred.type + algorithm, 
           data = fit.stats[fit.stats$pred == "Multi-omics\n+\nClinical",])
summary(model)
model = lm(RMSE ~ pred + algorithm, 
           data = fit.stats[fit.stats$pred.type == "Model Prediction",])
summary(model)



# Cytokines
with(fit.stats[fit.stats$algorithm == "Elastic net regression" & 
                     fit.stats$covariates == "Without covariates" & 
                     fit.stats$pred == "Cytokines", ], 
     t.test(RMSE ~ pred.type))

# Covariates
with(fit.stats[fit.stats$algorithm == "Elastic net regression" & 
                     fit.stats$covariates == "With covariates" & 
                     fit.stats$pred == "Covariates", ], 
     t.test(RMSE ~ pred.type))



# 5.7 Variable importance--------------------------------


## Reduce data to most important variables
vis_varImp_long = varImp_long %>%
      filter(vars %in% head(unique(varImp_long$vars), 20))

# Set factor levels
vis_varImp_long$vars = factor(vis_varImp_long$vars, levels = unique(vis_varImp_long$vars))



## Visualise top 20 variables for model prediction only
ggplot(vis_varImp_long, aes(x = vars, y = varImp)) +
      stat_boxplot(aes(fill = algorithm), col = "black", outlier.alpha = 0, alpha = 0.8) +
      geom_jitter(aes(fill = algorithm, col = algorithm), size = 0.1, alpha = 0.2) +
      scale_fill_brewer(palette = "Dark2") +
      scale_color_brewer(palette = "Dark2") +
      #scale_x_discrete(limits = rev(levels(vis_varImp.stats.comb$var))) +
      facet_grid(algorithm~.) +
      #coord_flip() +
      labs(x = "", y = "Variable Importance") +
      theme_bw() +
      theme(legend.position = "top",
            legend.title = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.background = element_rect(fill = "transparent", colour = NA),
            #panel.background = element_rect(fill = "transparent", colour = NA),
            legend.background = element_rect(fill = "transparent", colour = NA))
ggsave(filename = paste0("./Plots/", Sys.Date(), "_varImp.png"), 
       bg = "transparent", dpi = 600, width = 8, height = 6)

## Visualise all variables and compared to Null Model
ggplot(na.omit(arrange(varImp.stats, var)), 
       aes(x = var, y = varImp)) +
      stat_boxplot(aes(fill = pred.type), col = "black", outlier.alpha = 0, alpha = 0.8,
                   position = position_dodge(width = 1)) +
      geom_jitter(aes(fill = pred.type, col = pred.type), size = 0.08, alpha = 0.3,
                  position = position_jitterdodge(dodge.width = 1)) +
      scale_fill_brewer(palette = "Dark2") +
      scale_color_brewer(palette = "Dark2") +
      scale_x_discrete(limits = rev(levels(varImp.stats.comb$ind$var))) +
      facet_grid(.~algorithm) +
      coord_flip() +
      labs(x = "", y = "Variable Importance") +
      theme_bw() +
      theme(legend.position = "top",
            legend.title = element_blank())


# Second attempt with geom_halves
ggplot(na.omit(arrange(varImp.stats, var)), 
       aes(x = var, y = varImp)) +
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




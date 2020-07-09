# -------------------------------------------------------
# Analysis-----------------------------------------------
# -------------------------------------------------------


# 1 Preparation------------------------------------------

## Load packages
library("tidyverse")
library("ggcorrplot")
library("RColorBrewer")
library("caret")
library("glmnet")
library("PubHelper")


## Load data
load("/binder/mgp/datasets/2020_ImmuneDepression/cytokine/Nils_Preprocessed/OPTIMA_Cytokine_PreprocessedData.RData")

## Load cytokine reference
load("/binder/mgp/datasets/2020_ImmuneDepression/cytokine/Nils_Preprocessed/Cytokine_Reference.RData")



# 1.1 Define Parameters----------------------------------

## Save predictors, covariates, and outcome variables
x = cyto_ref[cyto_ref$no.na == 1, "vars"]
y = "bdi_locf_improve"
covariates = c("age_std", "sex_std")

## Define trainControl
fitControl = trainControl(method = "repeatedcv",
                          number = 10,
                          repeats = 10,
                          savePredictions = TRUE)





# 2 Descriptive statistics-------------------------------

# 2.1 Baseline Tables------------------------------------

## Baseline Table of clinical and sociodemographic characteristics
Table1 = baselineTable(dat, round_dec = 2, placeholder = "   ", grouping.var = "hsCRP_inflamed",
                       vars = c("sex", "age", "socialclass_cat", 
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
TableS1 = baselineTable(dat, round_dec = 2, placeholder = "   ", 
                        grouping.var = "hsCRP_inflamed",
                        vars = cyto_ref$vars,
                        labels = cyto_ref$labels)

exportPubHelpercsv(Table1, file = "./Results/Table1.csv")
exportPubHelpercsv(TableS1, file = "./Results/TableS1.csv")


# 3 Model Training---------------------------------------

## Create empty data.frame to extract fit statistics
fit.stats = data.frame(method = character(),
                       model = character(),
                       r2 = numeric(),
                       rmse = numeric(),
                       stringsAsFactors = FALSE)


# 3.1 Elastic Net Regression-----------------------------


#trainControl(method = "boot")

# Zou and Hastie’s(2005) recommended default alpha = 0.5.



# 3.1.1 NULL Model---------------------------------------

set.seed(1)
model.null.enet = train(x = dat[,c(x, covariates)], y = sample(dat[,y], replace = FALSE), 
                        method = "glmnet", metric = "RMSE", #alpha = 0.5,
                        trControl = fitControl)
print(model.null.enet)

## Save fit statistics
fit.stats[nrow(fit.stats) + 1, ] = NA
fit.stats[nrow(fit.stats), c("method", "model")] = c("Elastic Net", "Null")
fit.stats[nrow(fit.stats), c("rmse", "r2")] = model.null.enet$results %>%
      filter(RMSE == min(RMSE)) %>% 
      select(RMSE, Rsquared)




# 3.1.2 Baseline Model-----------------------------------

set.seed(1)
model.base.enet = train(x = dat[,covariates], y = dat[,y], 
                         method = "glmnet", metric = "RMSE", #alpha = 0.5,
                         trControl = fitControl)
print(model.base.enet)

## Save fit statistics
fit.stats[nrow(fit.stats) + 1, ] = NA
fit.stats[nrow(fit.stats), c("method", "model")] = c("Elastic Net", "Baseline")
fit.stats[nrow(fit.stats), c("rmse", "r2")] = model.base.enet$results %>%
      filter(RMSE == min(RMSE)) %>% 
      select(RMSE, Rsquared)

## Save LOOCF predictions
tuningvalues = as.numeric(model.base.enet$finalModel$tuneValue)
dat$pred.base.enet = model.base.enet$pred %>% 
      filter(alpha == tuningvalues[1] & lambda == tuningvalues[2]) %>% 
      arrange(rowIndex) %>% pull(pred)



# 3.1.3 Cytokine Model-----------------------------------


set.seed(1)
model.unadj.enet = train(x = dat[,x], y = dat[,y], 
                         method = "glmnet", metric = "RMSE", #alpha = 0.5,
                         trControl = fitControl)
print(model.unadj.enet)

## Save fit statistics
fit.stats[nrow(fit.stats) + 1, ] = NA
fit.stats[nrow(fit.stats), c("method", "model")] = c("Elastic Net", "Cytokines")
fit.stats[nrow(fit.stats), c("rmse", "r2")] = model.unadj.enet$results %>%
      filter(RMSE == min(RMSE)) %>% 
      select(RMSE, Rsquared)

## Save LOOCF predictions
tuningvalues = as.numeric(model.unadj.enet$finalModel$tuneValue)
dat$pred.unadj.enet = model.unadj.enet$pred %>% 
      filter(alpha == tuningvalues[1] & lambda == tuningvalues[2]) %>% 
      arrange(rowIndex) %>% pull(pred)


# 3.1.4 Combined-----------------------------------------


set.seed(1)
model.adj.enet = train(x = dat[,c(x, covariates)], y = dat[,y], 
                         method = "glmnet", metric = "RMSE", #alpha = 0.5,
                         trControl = fitControl)
print(model.adj.enet)

## Save fit statistics
fit.stats[nrow(fit.stats) + 1, ] = NA
fit.stats[nrow(fit.stats), c("method", "model")] = c("Elastic Net", "Baseline+Cytokines")
fit.stats[nrow(fit.stats), c("rmse", "r2")] = model.adj.enet$results %>%
      filter(RMSE == min(RMSE)) %>% 
      select(RMSE, Rsquared)


## Save LOOCF predictions
tuningvalues = as.numeric(model.adj.enet$finalModel$tuneValue)
dat$pred.adj.enet = model.adj.enet$pred %>% 
      filter(alpha == tuningvalues[1] & lambda == tuningvalues[2]) %>% 
      arrange(rowIndex) %>% pull(pred)



# 3.2 Random Forest--------------------------------------

# Surpassing Variable importance threshold?


# 3.2.1 Null Model---------------------------------------

set.seed(1)
model.null.rf = train(x = dat[,c(x, covariates)], y = sample(dat[,y], replace = FALSE), 
                      method = "rf", metric = "RMSE", #alpha = 0.5,
                      trControl = fitControl)
print(model.null.rf)

## Save fit statistics
fit.stats[nrow(fit.stats) + 1, ] = NA
fit.stats[nrow(fit.stats), c("method", "model")] = c("Random Forest", "Null")
fit.stats[nrow(fit.stats), c("rmse", "r2")] = model.null.rf$results %>%
      filter(RMSE == min(RMSE)) %>% 
      select(RMSE, Rsquared)



# 3.2.2 Baseline Model-----------------------------------

set.seed(1)
model.base.rf = train(x = dat[,covariates], y = dat[,y], 
                        method = "rf", metric = "RMSE", #alpha = 0.5,
                        trControl = fitControl)
print(model.base.rf)

## Save fit statistics
fit.stats[nrow(fit.stats) + 1, ] = NA
fit.stats[nrow(fit.stats), c("method", "model")] = c("Random Forest", "Baseline")
fit.stats[nrow(fit.stats), c("rmse", "r2")] = model.base.rf$results %>%
      filter(RMSE == min(RMSE)) %>% 
      select(RMSE, Rsquared)


## Save LOOCF predictions
dat$pred.base.rf = model.base.rf$pred %>% 
      filter(mtry == as.numeric(model.base.rf$bestTune)) %>% 
      arrange(rowIndex) %>% pull(pred)


# 3.1.3 Cytokine Model-----------------------------------


set.seed(1)
model.unadj.rf = train(x = dat[,x], y = dat[,y], 
                         method = "rf", metric = "RMSE", #alpha = 0.5,
                         trControl = fitControl)
print(model.unadj.rf)

## Save fit statistics
fit.stats[nrow(fit.stats) + 1, ] = NA
fit.stats[nrow(fit.stats), c("method", "model")] = c("Random Forest", "Cytokines")
fit.stats[nrow(fit.stats), c("rmse", "r2")] = model.unadj.rf$results %>%
      filter(RMSE == min(RMSE)) %>% 
      select(RMSE, Rsquared)


## Save LOOCF predictions
dat$pred.unadj.rf = model.unadj.rf$pred %>% 
      filter(mtry == as.numeric(model.unadj.rf$bestTune)) %>% 
      arrange(rowIndex) %>% pull(pred)


# 3.1.4 Combined-----------------------------------------


set.seed(1)
model.adj.rf = train(x = dat[,c(x, covariates)], y = dat[,y], 
                     method = "rf", metric = "RMSE", #alpha = 0.5,
                     trControl = fitControl)
print(model.adj.rf)

## Save fit statistics
fit.stats[nrow(fit.stats) + 1, ] = NA
fit.stats[nrow(fit.stats), c("method", "model")] = c("Random Forest", "Baseline+Cytokines")
fit.stats[nrow(fit.stats), c("rmse", "r2")] = model.adj.rf$results %>%
      filter(RMSE == min(RMSE)) %>% 
      select(RMSE, Rsquared)


## Save LOOCF predictions
dat$pred.adj.rf = model.adj.rf$pred %>% 
      filter(mtry == as.numeric(model.adj.rf$bestTune)) %>% 
      arrange(rowIndex) %>% pull(pred)



# 4 Content validity analysis----------------------------


# 4.1 MADRS----------------------------------------------



# 4.2 BSI-Depression-------------------------------------


# 4.3 WHODAS---------------------------------------------


# 4.4 CIDI Depression------------------------------------


# 4.5 Dropout--------------------------------------------


# 5 Visualisation----------------------------------------


# 5.1 Preparation----------------------------------------

## Create data
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





# 5.2 Cytokine Correlations------------------------------


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



# 5.3 ML Prediction Scatter Plots------------------------

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



# 5.3 Absolute Residual Boxplots-------------------------

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



# 5.5 Fit Statistics-------------------------------------

## Set factor levels
fit.stats$model = factor(fit.stats$model, levels = c("Null", "Baseline", "Cytokines", 
                                                     "Baseline+Cytokines"))

ggplot(fit.stats, aes(x = method, y = r2)) +
      geom_bar(aes(fill = model), stat = "identity", position = position_dodge(), col = "black") +
      scale_fill_brewer(palette = "Dark2") +
      scale_y_continuous(limits = c(0, 0.15)) +
      coord_cartesian(expand = FALSE) +
      labs(x = "Model", y = expression(R^2)) +
      theme_bw() +
      theme(legend.position = "top",
            legend.title = element_blank())


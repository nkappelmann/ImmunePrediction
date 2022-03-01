# -------------------------------------------------------
# Evaluation---------------------------------------------
# -------------------------------------------------------


# 1 Preparation------------------------------------------

## Load packages
library("tidyverse")
library("tableone")
library("ggcorrplot")
library("RColorBrewer")

## Set Slurm directory
setwd("/home/nkappelmann/OPTIMA/ImmunePrediction")


## Load data
load("./Data/OPTIMA_Cytokine_PreprocessedData.RData")


## Load ML model output
load("./Results/clinical.output.RData")
load("./Results/cytokine.output.RData")
load("./Results/rna.output.RData")
load("./Results/prs.output.RData")
load("./Results/omics.output.RData")
load("./Results/omicsplusclin.output.RData")

load("./Results/clinical.output.som.RData")
load("./Results/cytokine.output.som.RData")
load("./Results/rna.output.som.RData")
load("./Results/prs.output.som.RData")
load("./Results/omics.output.som.RData")
load("./Results/omicsplusclin.output.som.RData")


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
clin.vars = c("t0_madrs", "sex", "age", "BMI", "t0_diagn_by_age",
              paste0("t0_cidi_", c("depression_episode", "depression_recurrent", 
                                   "anxiety", "substance")), 
              paste0("t0_pid_", c("negaff", "detach", "psycho", "antago", "disinh")),
              paste0("t0_bsi_", c("soma", "zwan", "unsi", "depr", "angs", "aggr", 
                                  "phob", "para", "psyc")))


## Source nested cross-validation function
source("./Scripts/functions.R")


# 2 Preparation------------------------------------------

# 2.1 Fit statistics-------------------------------------

# 2.1.1 bdi_locf_improve---------------------------------

## Fit statistics are combined
fit_stats = collate_fit_results(list(
   clinical.output$fit, cytokine.output$fit, prs.output$fit, 
   rna.output$fit, omics.output$fit, omicsplusclin.output$fit
   ))


# A results table is created
fit_stats_table = get_fit_results_table(fit_stats)

write_csv(fit_stats_table, file = "./Results/Fit_Stats_Table.csv")


# 2.1.1 bdi_locf_som_improve-----------------------------

fit_stats_som = collate_fit_results(list(
   clinical.output.som$fit, cytokine.output.som$fit, prs.output.som$fit, 
   rna.output.som$fit, omics.output.som$fit, omicsplusclin.output.som$fit
))


# A results table is created
fit_stats_som_table = get_fit_results_table(fit_stats_som)

write_csv(fit_stats_som_table, file = "./Results/Fit_Stats_Som_Table.csv")


# 3 Descriptive statistics-------------------------------

# 3.1 Baseline Tables------------------------------------

# Prepare data
table1_dat = dat %>%
   select(sex, age, ward_type, ethnicity, BMI, t0_bdi, t0_madrs, t0_bsi_depr)
table1_dat$included = ifelse(ids_for_analysis, "Included", "Excluded") %>%
   factor(levels = c("Included", "Excluded"))

table1_vars = colnames(table1_dat)[colnames(table1_dat) != "included"]
table1_catvars = c("sex", "ward_type", "ethnicity", "included")



## Create and export table
tab1 = CreateTableOne(vars = table1_vars, strata = "included", 
                       data = table1_dat, factorVars = table1_catvars)


tab1Mat = print(tab1, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE,
                showAllLevels = TRUE)
write.table(tab1Mat, file = "./Results/Table1.csv", quote = FALSE, sep = ";")



# 3.2 Treatment improvements-----------------------------

# 3.2.1 BDI change---------------------------------------

# Overall
by(dat$bdi_locf_improve, ids_for_analysis, summary)
by(dat$bdi_locf_improve, ids_for_analysis, sd)

# SMD
with(dat[ids_for_analysis], summary(bdi_locf_improve / sd(bdi_locf_improve)))

# 3.2.2 BDI somatic change-------------------------------

# Overall
by(dat$bdi_som_locf_improve, ids_for_analysis, summary)
by(dat$bdi_som_locf_improve, ids_for_analysis, sd)

# SMD
with(dat[ids_for_analysis], summary(bdi_som_locf_improve / sd(bdi_locf_improve)))


# 4 Model Evaluation-------------------------------------



## Plot results
plot_fit_results(fit_stats)
ggsave(filename = paste0("./Plots/", Sys.Date(), "_PredPerformance_RMSE.png"), 
       device = "png", width = 8, height = 8, dpi = 300, bg = "transparent")

plot_fit_results(fit_stats_som)
ggsave(filename = paste0("./Plots/", Sys.Date(), "_PredPerformanceSom_RMSE.png"), 
       device = "png", width = 8, height = 8, dpi = 300, bg = "transparent")


# RMSE delta
ggplot(fit_stats, aes(x = feature.set, y = RMSE_delta)) +
   geom_hline(yintercept = 0, col = "lightgrey", size = 2) +
   geom_violin(aes(fill = algorithm), alpha = 0.5, position = position_dodge(0)) + # outlier.alpha = 0, 
   #geom_jitter(aes(col = feature.set), size = 0.8, alpha = 0.5) +
   stat_summary(fun = "mean", geom = "crossbar", width = 0.3, col = "black") +
   stat_summary(fun = "mean", geom = "point", size = 2, col = "black") +
   scale_fill_brewer(palette = "Dark2") +
   scale_color_brewer(palette = "Dark2") +
   labs(x = "", y = bquote(Delta*"RMSE (Model-based error reduction)")) +
   facet_grid(algorithm ~ .) +
   theme_bw() +
   theme(legend.position = "top",
         legend.title = element_blank(),
         strip.background = element_rect(colour="black", fill="white"),
         strip.text = element_text(face = "bold", size = 15))




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




# 5.1.2 Variable importance------------------------------


plot_variable_importance(omicsplusclin.output)
ggsave(filename = paste0("./Plots/", Sys.Date(), "_VarImp.png"), 
       device = "png", width = 6, height = 6, dpi = 300, bg = "transparent")

plot_variable_importance(omicsplusclin.output.som)
ggsave(filename = paste0("./Plots/", Sys.Date(), "_VarImp_Som.png"), 
       device = "png", width = 6, height = 6, dpi = 300, bg = "transparent")




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



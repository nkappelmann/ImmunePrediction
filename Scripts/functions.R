# -------------------------------------------------------
# Functions----------------------------------------------
# -------------------------------------------------------

# Load packages------------------------------------------

library("tidyverse")
library("tidymodels")
library("furrr")
library("vip")
library("progressr")
library("RColorBrewer")


# ml_pipeline--------------------------------------------

ml_pipeline <- function(
   df,
   x,
   y,
   k_outer = 10, # number of outer CV folds
   k_inner = 10, # number of inner CV folds
   num_repeats = 100, # number of repeats
   parallel = TRUE, # parallelise over repetitions
   seed = 1, # set Seed for analyses
   show_progress = TRUE # Show progress bar for repetitions
)  {
   
   
   ## Subset data to relevant variables only
   df = df[, c(x, y)] %>%
      mutate_if(is.character, as.factor)
   colnames(df)[colnames(df) == y] = "y"
   
   ## Set running seed, which is incremented by 1 every time it is used
   runningSeed = seed
   
   ## Define parameters
   obs.all = nrow(df)
   fold_vector = rep_len(1:k_outer, length.out = obs.all)
   
   ## Define models
   # Random forest
   rf_spec = rand_forest(mtry = tune(), trees = 1000, min_n = tune()) %>%
      set_engine("ranger", importance = "impurity") %>% 
      set_mode("regression")
   
   rf_grid = grid_latin_hypercube(finalize(mtry(), df[, x]), min_n(), size = 10)
   
   # Elastic net
   glmnet_spec = linear_reg(penalty = tune(), mixture = tune()) %>%
      set_engine("glmnet") %>%
      set_mode("regression")
   
   glmnet_grid = grid_latin_hypercube(penalty(), mixture(), size = 10)
   
   # SVM
   # svm_spec = svm_poly(cost = tune(), degree = tune()) %>%
   #    set_engine("kernlab") %>%
   #    set_mode("regression")
   # 
   # svm_grid = grid_latin_hypercube(cost(), degree(), size = 10)
   
   # k-nearest neighbour
   knn_spec =
      nearest_neighbor(neighbors = tune(), weight_func = tune()) %>%
      set_engine("kknn") %>%
      set_mode("regression")

   knn_grid = grid_latin_hypercube(neighbors(), weight_func(), size = 10)
   
   
   # neural net
   # nnet_spec = mlp(penalty = tune(), hidden_units = tune()) %>%
   #    set_engine("nnet") %>%
   #    set_mode("regression")
   # 
   # nnet_grid = grid_latin_hypercube(penalty(), hidden_units(), size = 10)
   

   
   # Create list with repetitions data for mapping function
   df_list = list()
   
   for(i in 1:num_repeats) {
      df_list[[i]] = df
      df_list[[i]]$seed = runningSeed + i
      
      set.seed(df_list[[i]]$seed)
      df_list[[i]]$fold_outer = sample(fold_vector)
      set.seed(df_list[[i]]$seed)
      df_list[[i]]$y_perm = sample(df$y)
   }
   
   
   # Inner CV-----------------------------
   
   inner_cv = function(
      wflow, wflow_perm, model_grid, inner_cv_splits, train_df, test_df, algorithm, outer
   )   {
      
      # Define unique seed
      inner_seed = as.numeric(paste0(train_df$seed[1], outer))
      
      # Fit and select best model
      set.seed(inner_seed)
      model_fit = wflow %>% 
         tune_grid(resamples = inner_cv_splits, grid = model_grid, metrics = metric_set(rmse))
      set.seed(inner_seed)
      model_fit_perm = wflow_perm %>% 
         tune_grid(resamples = inner_cv_splits, grid = model_grid, metrics = metric_set(rmse))
      
      best_model = model_fit %>%
         select_best("rmse")
      best_model_perm = model_fit_perm %>%
         select_best("rmse")
      
      # Refit to full training data
      set.seed(inner_seed)
      final_wflow = wflow %>% 
         finalize_workflow(best_model)
      set.seed(inner_seed)
      final_wflow_perm = wflow_perm %>% 
         finalize_workflow(best_model_perm)
      
      set.seed(inner_seed)
      final_model_fit = final_wflow %>%
         fit(train_df)
      set.seed(inner_seed)
      final_model_fit_perm = final_wflow_perm %>%
         fit(train_df)
      
      test_pred = predict(final_model_fit, test_df)
      test_pred_perm = predict(final_model_fit_perm, test_df)
      
      # Get fit indices
      fit_output = tibble(
         k = outer,
         algorithm = algorithm
      )
      fit_output[, c("RMSE", "R2")] = t(get_fit_indices(test_df, test_pred, y = "actual"))
      fit_output[, c("RMSE_perm", "R2_perm")] = t(get_fit_indices(test_df, test_pred_perm, y = "perm"))
      
      # Get hyperparameters
      fit_output[, c("mtry", "min_n", "penalty", "mixture", 
                     "cost", "degree", "neighbors", "weight_func", "hidden_units")] = NA
      if(algorithm == "rf")   {
         fit_output[, c("mtry", "min_n")] = best_model[, c("mtry", "min_n")]
      } else if(algorithm == "glmnet") {
         fit_output[, c("penalty", "mixture")] = best_model[, c("penalty", "mixture")]
      } else if(algorithm == "svm") {
         fit_output[, c("cost", "degree")] = best_model[, c("cost", "degree")]
      } else if(algorithm == "knn") {
         fit_output[, c("neighbors", "weight_func")] = best_model[, c("neighbors", "weight_func")]
      } else if(algorithm == "nnet") {
         fit_output[, c("penalty", "hidden_units")] = best_model[, c("penalty", "hidden_units")]
      } else {stop("Please specify the algorithm correctly")}
      
      # Get variable importances
      if(algorithm %in% c("rf", "glmnet", "nnet"))   {
         importance_output = final_model_fit %>%   
            extract_fit_parsnip() %>% 
            vi() %>%
            select(Variable, Importance) %>%
            pivot_wider(names_from = Variable, values_from = Importance) %>%
            mutate(k = outer, algorithm = algorithm)
      } else   {
         importance_output = tibble(Variable = x, Importance = NA) %>%
            pivot_wider(names_from = Variable, values_from = Importance) %>%
            mutate(k = outer, algorithm = algorithm)
      }
      
      
      output = list(fit = fit_output, importance = importance_output)
      
      return(output)
      
   }
   
   # Outer CV-----------------------------
      
   outer_cv = function(object)   {
      
      for(outer in 1:k_outer) {
         
         if(outer == 1) {cat("\tFold no. 1----")
         } else {cat(outer, "----", sep = "")}
         
         train_df = object[object$fold_outer != k_outer,]
         test_df = object[object$fold_outer == k_outer,]
         
         outer_cv_seed = object$seed
         
         set.seed(outer_cv_seed)
         inner_cv_splits = vfold_cv(train_df, v = k_inner)
         outer_cv_seed = outer_cv_seed + 1
         
         # Create Recipe
         rec = create_recipe(train_df, y = "actual")
         rec_perm = create_recipe(train_df, y = "permuted")
         
         
         ## Random forest
         rf_wflow = workflow() %>% 
            add_model(rf_spec)
         rf_wflow_perm = rf_wflow %>% add_recipe(rec_perm)
         rf_wflow = rf_wflow %>% add_recipe(rec)
         
         rf_output = inner_cv(wflow = rf_wflow, wflow_perm = rf_wflow_perm, 
                              model_grid = rf_grid, inner_cv_splits = inner_cv_splits, 
                              train_df = train_df, test_df = test_df, 
                              algorithm = "rf", outer = outer)
         
         
         ## Elastic net regression
         glmnet_wflow = workflow() %>% 
            add_model(glmnet_spec)
         glmnet_wflow_perm = glmnet_wflow %>% add_recipe(rec_perm)
         glmnet_wflow = glmnet_wflow %>% add_recipe(rec)
         
         glmnet_output = inner_cv(wflow = glmnet_wflow, wflow_perm = glmnet_wflow_perm, 
                                  model_grid = glmnet_grid, inner_cv_splits = inner_cv_splits, 
                                  train_df = train_df, test_df = test_df, 
                                  algorithm = "glmnet", outer = outer)
         
         
         ## k-Nearest Neighbour
         knn_wflow = workflow() %>%
            add_model(knn_spec)
         knn_wflow_perm = knn_wflow %>% add_recipe(rec_perm)
         knn_wflow = knn_wflow %>% add_recipe(rec)

         knn_output = inner_cv(wflow = knn_wflow, wflow_perm = knn_wflow_perm,
                               model_grid = knn_grid, inner_cv_splits = inner_cv_splits, 
                               train_df = train_df, test_df = test_df, 
                               algorithm = "knn", outer = outer)
         
         
         ## SVM
         # svm_wflow = workflow() %>%
         #    add_model(svm_spec)
         # svm_wflow_perm = svm_wflow %>% add_recipe(rec_perm)
         # svm_wflow = svm_wflow %>% add_recipe(rec)
         # 
         # svm_output = inner_cv(wflow = svm_wflow, wflow_perm = svm_wflow_perm,
         #                       model_grid = svm_grid, algorithm = "svm",
         #                       inner_cv_splits = inner_cv_splits, 
         #                       train_df = train_df, test_df = test_df, outer = outer)   
         
         ## Neural network
         # nnet_wflow = workflow() %>%
         #    add_model(nnet_spec)
         # nnet_wflow_perm = nnet_wflow %>% add_recipe(rec_perm)
         # nnet_wflow = nnet_wflow %>% add_recipe(rec)
         # 
         # nnet_output = inner_cv(wflow = nnet_wflow, wflow_perm = nnet_wflow_perm, 
         #                        model_grid = nnet_grid, inner_cv_splits = inner_cv_splits, 
         #                        train_df = train_df, test_df = test_df, 
         #                        algorithm = "nnet", outer = outer)   
         
         if(exists("temp_output")) {rm(temp_output)}
         temp_output = purrr::map2(rf_output, glmnet_output, bind_rows)
         temp_output = purrr::map2(temp_output, knn_output, bind_rows)
         # temp_output = purrr::map2(temp_output, nnet_output, bind_rows)
         
         ## Merge results
         if(outer == 1)   {
            output = temp_output
         } else   {
            output = purrr::map2(output, temp_output, bind_rows)
         }
         
      }
      
      return(output)
      
   } 
   
   # Set up parallelisation
   if(parallel == TRUE) {
      if(future::supportsMulticore())  {
         plan(multicore)
      } else{
         plan(multisession)
      }
      
      outer_cv_mapping = function(x) {
         p = progressor(steps = length(x))
         
         result = future_map_dfr(x, ~{
            p()
            outer_cv(.x)
         }, .options = furrr_options(seed = TRUE))
         
         return(result)
      }
      
      # Run processing with/without progress bar
      if(show_progress) {
         with_progress({
            nested_cv_results = outer_cv_mapping(df_list)
         })
         
      } else   {nested_cv_results = outer_cv_mapping(df_list)}
   } else   {nested_cv_results = map_dfr(df_list, outer_cv)}
   
   
   
   return(nested_cv_results)
}



# create_recipe------------------------------------------

create_recipe = function(train, y_type = "actual")  {
   
   if(y_type == "actual") {
      rec = recipe(y ~ ., data = train) %>%
         update_role(y_perm, new_role = "id variable")
   } else if(y_type == "permuted") {
      rec = recipe(y_perm ~ ., data = train) %>%
         update_role(y, new_role = "id variable")
   } else   {stop("Please define y as 'actual' or 'permuted'")}
   rec = rec %>%
      update_role(seed, new_role = "id variable") %>%
      update_role(fold_outer, new_role = "id variable") %>%
      step_dummy(all_nominal_predictors()) %>%
      step_zv(all_predictors()) %>%
      step_normalize(all_numeric_predictors()) %>%
      step_impute_knn(all_predictors())
   
   return(rec)
}


# get_fit_indices----------------------------------------

get_fit_indices = function(df, pred, y_type = "actual")   {
   output = c(
      RMSE = rmse_vec(truth = df[, ifelse(y_type == "actual", "y", "y_perm")], 
                      estimate = pred$.pred),
      R2 = rsq_vec(truth = df[, ifelse(y_type == "actual", "y", "y_perm")], 
                        estimate = pred$.pred)
   )
   return(output)
}


# collate_fit_results------------------------------------

collate_fit_results = function(
   fit_results # list of fit results
)  {
   
   
   # Fit statistics are combined
   fit_stats = bind_rows(fit_results)
   
   # Meta-data is added
   fit_stats$feature.set = rep(c("Clinical", "Cytokines", "PRS", 
                                 "Gene expression", "Multi-omics",
                                 "Multi-omics\n+\nClinical"), each = nrow(fit_results[[1]])) %>%
      factor(levels = c("Clinical", "Cytokines", "PRS", "Gene expression", "Multi-omics", 
                        "Multi-omics\n+\nClinical"))
   
   # Algorithm is recoded
   fit_stats$algorithm = recode(fit_stats$algorithm,
                                'rf' = "Random forest",
                                'glmnet' = "Elastic net",
                                'knn' = "k-Nearest neighbour") %>%
      factor(levels = c("Elastic net", "Random forest", "k-Nearest neighbour"))
   
   # RMSE delta is calculated
   fit_stats$RMSE_delta = fit_stats$RMSE - fit_stats$RMSE_perm
   
   # A binary variable to quantify model prediction over permutation
   fit_stats$model_better = ifelse(fit_stats$RMSE_delta < 0, "Better than chance", "Worse than chance")
   
   
   return(fit_stats)
   
}


# get_fit_results_table----------------------------------

get_fit_results_table = function(
   fit_results # collate_fit_results output
)  {
   
   if(sum(is.na(fit_results$R2) > 0))  {
      warning(paste0("There are ", sum(is.na(fit_results$R2)), " NAs in fit_results$R2."))
   }
   
   fit_stats_table = fit_results %>%
      group_by(feature.set, algorithm) %>%
      summarise(RMSE_MeanSD = paste0(round(mean(RMSE), 2), " (", round(sd(RMSE), 2), ")"),
                R2_MedianIQR = paste0(round(median(R2, na.rm = TRUE), 2), " (", round(quantile(R2, 0.25, na.rm = TRUE), 2), ", ",
                                      round(quantile(R2, 0.75, na.rm = TRUE), 2), ")"),
                RMSE_delta_MeanSD = paste0(round(mean(RMSE_delta), 2), " (", round(sd(RMSE_delta), 2), ")"),
                RMSE_delta_pval = 1 - sum(RMSE_delta < 0) / n())
   
   fit_stats_table$feature.set = gsub("\n", "", fit_stats_table$feature.set)
   
   
   return(fit_stats_table)
}



# plot_fit_results---------------------------------------

plot_fit_results = function(
   dat # collate_fit_results output
)  {
   
   # Note: Code adapted from: 
   # https://stackoverflow.com/questions/36203195/fill-specific-regions-in-geom-violin-plot
   
   p = ggplot() + 
      geom_violin(data = dat, aes(x = feature.set, y = RMSE_delta)) +
      facet_grid(algorithm ~ .)
   p_build = ggplot2::ggplot_build(p)$data[[1]]
   
   #This comes directly from the source of geom_violin
   p_build = transform(p_build,
                       xminv = x - violinwidth * (x - xmin),
                       xmaxv = x + violinwidth * (xmax - x))
   
   p_build = rbind(plyr::arrange(transform(p_build, x = xminv), y),
                   plyr::arrange(transform(p_build, x = xmaxv), -y))
   
   #Add our fill variable
   p_build$fill_group = ifelse(p_build$y >= 0,'Null model','Actual model')
   p_build$group1 = with(p_build,interaction(factor(group), factor(fill_group)))
   
   p_build$algorithm = p_build$PANEL %>%
      recode('1' = levels(dat$algorithm)[1],
             '2' = levels(dat$algorithm)[2],
             '3' = levels(dat$algorithm)[3])
   
   #Note the use of the group aesthetic here with our computed version,
   # group1
   p_fill = ggplot() + 
      geom_hline(yintercept = 0, col = "darkgrey", size = 0.5) +
      geom_polygon(data = p_build,
                   aes(x = x, y = y, group = group1, fill = fill_group),
                   alpha = 0.8, col = "black") +
      scale_fill_brewer(palette = "Dark2") +
      scale_color_brewer(palette = "Dark2") +
      scale_x_continuous(breaks = 1:length(levels(dat$feature.set)),
                         labels = levels(dat$feature.set)) +
      labs(x = "", y = bquote(Delta*"RMSE (Model-based error reduction)"),
           fill = "Better prediction by: ") +
      facet_grid(algorithm ~ .) +
      theme_bw() +
      theme(legend.position = "top",
            strip.background = element_rect(colour="black", fill="white"),
            strip.text = element_text(face = "bold", size = 15))
   
   return(p_fill)
   
}


# plot_variable_importance-------------------------------

plot_variable_importance = function(
   ml_output, # output from ml_pipeline function
   top_n_predictors = 20, # number of top predictors to be plotted
   precision_type = "IQR" # Alternative: CI
)  {
   
   var_imp = ml_output$importance %>%
      select(!c(sex, k)) %>%
      filter(algorithm != "knn") %>%
      group_by(algorithm) %>%
      pivot_longer(!contains("algorithm"), names_to = "variable", values_to = "importance") %>%
      group_by(algorithm, variable) %>%
      summarise(mean = mean(importance), 
                sd = sd(importance), 
                ci.lb = quantile(importance, probs = 0.025),
                ci.ub = quantile(importance, probs = 0.975),
                median = median(importance), 
                iqr.lb = quantile(importance, probs = 0.25),
                iqr.ub = quantile(importance, probs = 0.75)) %>%
      group_by(algorithm) %>%
      arrange(algorithm, desc(median)) %>%
      slice_head(n = top_n_predictors) %>%
      arrange(algorithm, median) %>%
      ungroup() %>%
      mutate(order = row_number())
   
   
   var_imp$algorithm = algorithm = recode(var_imp$algorithm,
                                          'rf' = "Random forest",
                                          'glmnet' = "Elastic net",
                                          'knn' = "k-Nearest neighbour") %>%
      factor(levels = c("Elastic net", "Random forest", "k-Nearest neighbour")) %>%
      droplevels()
   
   var_imp$variable_group = ifelse(var_imp$variable %in% clin.vars, "Clinical",
                                   ifelse(var_imp$variable %in% cyto_ref$vars, "Cytokines", 
                                          ifelse(var_imp$variable %in% rna.vars, "Gene expression",
                                                 "PRS")))
   
   
   var_imp$variable = gsub("rna.", "", var_imp$variable, fixed = TRUE)
   var_imp$variable = gsub("_PRS", "", var_imp$variable, fixed = TRUE)
   var_imp$variable = gsub("t0_", "", var_imp$variable, fixed = TRUE)
   var_imp$variable = gsub("VEGF.A", "VEGF-A", var_imp$variable, fixed = TRUE)
   var_imp$variable = gsub("cidi", "CIDI", var_imp$variable, fixed = TRUE)
   var_imp$variable = gsub(".", "", var_imp$variable, fixed = TRUE)
   var_imp$variable = gsub("alpha", "a", var_imp$variable, fixed = TRUE)
   var_imp$variable = gsub("beta", "b", var_imp$variable, fixed = TRUE)
   var_imp$variable = gsub("diagn_by_age", "Age-corrected\n#diagnoses", var_imp$variable, fixed = TRUE)
   var_imp$variable = gsub("bsi", "BSI", var_imp$variable, fixed = TRUE)
   var_imp$variable = gsub("depr", "Depression", var_imp$variable, fixed = TRUE)
   var_imp$variable = gsub("Depressionession_episode", "depressive\nepisode", var_imp$variable, fixed = TRUE)
   var_imp$variable = gsub("Depressionession_recurrent", "recurrent\ndepression", var_imp$variable, fixed = TRUE)
   var_imp$variable = gsub("unsi", "Interper-\nsonal Sensitivity", var_imp$variable, fixed = TRUE)
   var_imp$variable = gsub("angs", "Anxiety", var_imp$variable, fixed = TRUE)
   var_imp$variable = gsub("phob", "Phobic\nAnxiety", var_imp$variable, fixed = TRUE)
   var_imp$variable = gsub("pid", "PID", var_imp$variable, fixed = TRUE)
   var_imp$variable = gsub("antago", "Antagonism", var_imp$variable, fixed = TRUE)
   var_imp$variable = gsub("disinh", "Disinhibition", var_imp$variable, fixed = TRUE)
   var_imp$variable = gsub("_", " ", var_imp$variable, fixed = TRUE)
   
   if(precision_type == "IQR")   {
      var_imp$lb = var_imp$iqr.lb
      var_imp$ub = var_imp$iqr.ub
   } else if(precision_type == "CI")   {
      var_imp$lb = var_imp$ci.lb
      var_imp$ub = var_imp$ci.ub
   } else   {stop("Please set precision_type to 'IQR' or 'CI'.")}
   
   varImp_plot = ggplot(var_imp, aes(x = order, y = median)) +
      geom_errorbar(aes(ymin = lb, ymax = ub)) +
      geom_point(aes(fill = variable_group), shape = 21, size = 3) +
      labs(x = "", y = paste0("Variable Importance: Median (", 
                              ifelse(precision_type == "IQR", "IQR", "95% CI"), ")")) +
      scale_fill_brewer(palette = "Set2") +
      facet_wrap(algorithm~., scales = "free") +
      scale_x_continuous(breaks = var_imp$order,
                         labels = var_imp$variable,
                         expand = c(0, 0)) +
      coord_flip() + 
      theme_bw() +
      theme(legend.position = "top",
            legend.title = element_blank(),
            strip.background = element_rect(colour="black", fill="white"),
            strip.text = element_text(face = "bold", size = 15))
   
   return(varImp_plot)
   
}

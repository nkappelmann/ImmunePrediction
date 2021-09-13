# -------------------------------------------------------
# Functions----------------------------------------------
# -------------------------------------------------------

# Load packages------------------------------------------

library("tidyverse")
library("tidymodels")
library("furrr")
library("vip")
library("baguette")



# ml_pipeline--------------------------------------------

ml_pipeline <- function(
   df,
   x,
   y,
   k_outer = 10, # number of outer CV folds
   k_inner = 10, # number of inner CV folds
   seed = 1 # set Seed for analyses
)  {
   
   
   ## Subset data to relevant variables only
   df = df[, c(x, y)] %>%
      mutate_if(is.character, as.factor)
   colnames(df)[colnames(df) == y] = "y"
   
   ## Set running seed, which is incremented by 1 every time it is used
   runningSeed = seed
   
   ## Define parameters
   obs.all = nrow(df)
   cores = parallel::detectCores() # to do: remove
   fold_vector = rep_len(1:k_outer, length.out = obs.all)
   
   ## Define models
   # Random forest
   rf_spec = rand_forest(mtry = tune(), trees = 1000, min_n = tune()) %>%
      set_engine("ranger") %>%
      set_mode("regression")
   
   rf_grid = grid_regular(finalize(mtry(), df[, x]), min_n(), levels = 5)
   
   # Elastic net
   glmnet_spec = linear_reg(penalty = tune(), mixture = tune()) %>%
      set_engine("glmnet") %>%
      set_mode("regression")
   
   glmnet_grid = grid_regular(penalty(), mixture(), levels = 5)
   
   # SVM
   svm_spec = svm_poly(cost = tune(), degree = tune()) %>%
      set_engine("kernlab") %>%
      set_mode("regression")
   
   svm_grid = grid_regular(cost(), degree(), levels = 5)
   
   
   ## Start nested CV
   for(num_repeat in 1:num_repeats) {
      
      cat("\nRepeat no.", num_repeat, "\n")
      
      # Assign folds and create permuted Y
      set.seed(runningSeed)
      fold_outer = sample(fold_vector)
      set.seed(runningSeed)
      df$y_perm = sample(df$y)
      runningSeed = runningSeed + 1
      
      
      for(outer in 1:k_outer) {
         
         if(outer == 1) {cat("\tFold no. 1----")
         } else {cat(outer, "----", sep = "")}
         
         train_df = df[fold_outer != num_repeat,]
         test_df = df[fold_outer == num_repeat,]
         
         set.seed(runningSeed)
         inner_cv_splits = vfold_cv(train_df, v = k_inner)
         runningSeed = runningSeed + 1
         
         
         # Create Recipe
         rec = create_recipe(train_df, y = "actual")
         rec_perm = create_recipe(train_df, y = "permuted")
         
         
         # Run inner CV----------------------
         
         inner_cv = function(wflow, wflow_perm, model_grid, algorithm)   {
            
            # Fit and select best model
            model_fit = wflow %>% 
               tune_grid(resamples = inner_cv_splits, grid = model_grid, metrics = metric_set(rmse))
            model_fit_perm = wflow_perm %>% 
               tune_grid(resamples = inner_cv_splits, grid = model_grid, metrics = metric_set(rmse))
            
            best_model = model_fit %>%
               select_best("rmse")
            best_model_perm = model_fit_perm %>%
               select_best("rmse")
            
            # Refit to full training data
            final_wflow = wflow %>% 
               finalize_workflow(best_model)
            final_wflow_perm = wflow_perm %>% 
               finalize_workflow(best_model_perm)
            
            final_model_fit = final_wflow %>%
               fit(train_df)
            final_model_fit_perm = final_wflow_perm %>%
               fit(train_df)
            
            test_pred = predict(final_model_fit, test_df)
            test_pred_perm = predict(final_model_fit_perm, test_df)
            
            ## Get fit indices
            fit_output = tibble(
               num_repeat = num_repeat,
               k = k_outer,
               algorithm = algorithm,
               hyperparameters = best_model
            )
            fit_output[, c("RMSE", "R2")] = t(get_fit_indices(test_df, test_pred, y = "actual"))
            fit_output[, c("RMSE_perm", "R2_perm")] = t(get_fit_indices(test_df, test_pred_perm, y = "perm"))
            
            fit_output[, c("mtry", "min_n", "penalty", "mixture", "cost", "degree")] = NA
            if(algorithm == "rf")   {
               fit_output[, c("mtry", "min_n")] = best_model[, c("mtry", "min_n")]
            } else if(algorithm == "glmnet") {
               fit_output[, c("penalty", "mixture")] = best_model[, c("penalty", "mixture")]
            } else if(algorithm == "svm") {
               fit_output[, c("cost", "degree")] = best_model[, c("cost", "degree")]
            }
            else {stop("Please specify the algorithm correctly")}
            
            
            ## Get variable importance
            #var_imp(final_model_fit)
            
            
            return(fit_output)
            
         }
         
         
         ## Random forest
         rf_wflow = workflow() %>% 
            add_model(rf_spec)
         rf_wflow_perm = rf_wflow %>% add_recipe(rec_perm)
         rf_wflow = rf_wflow %>% add_recipe(rec)
         
         rf_output = inner_cv(wflow = rf_wflow, wflow_perm = rf_wflow_perm, 
                              model_grid = rf_grid, algorithm = "rf")
         
         
         ## Elastic net regression
         glmnet_wflow = workflow() %>% 
            add_model(glmnet_spec)
         glmnet_wflow_perm = glmnet_wflow %>% add_recipe(rec_perm)
         glmnet_wflow = glmnet_wflow %>% add_recipe(rec)
         
         glmnet_output = inner_cv(wflow = glmnet_wflow, wflow_perm = glmnet_wflow_perm, 
                                  model_grid = glmnet_grid, algorithm = "glmnet")
         
         
         ## k-Nearest Neighbour
         
         ## SVM
         svm_wflow = workflow() %>%
            add_model(svm_spec)
         svm_wflow_perm = svm_wflow %>% add_recipe(rec_perm)
         svm_wflow = svm_wflow %>% add_recipe(rec)
         
         svm_output = inner_cv(wflow = svm_wflow, wflow_perm = svm_wflow_perm, 
                                  model_grid = svm_grid, algorithm = "svm")   
         
         ## Random forest
         
         # Create workflow
         rf_wflow = workflow() %>% 
            add_model(rf_spec)
         rf_wflow_perm = rf_wflow %>% add_recipe(rec_perm)
         rf_wflow = rf_wflow %>% add_recipe(rec)
         
         # Fit
         rf_fit = rf_wflow %>% 
            tune_grid(resamples = inner_cv_splits, grid = rf_grid, metrics = metric_set(rmse))
         rf_fit_perm = rf_wflow_perm %>% 
            tune_grid(resamples = inner_cv_splits, grid = rf_grid, metrics = metric_set(rmse))
         
         # Select best model
         best_rf = rf_fit %>%
            select_best("rmse")
         best_rf_perm = rf_fit_perm %>%
            select_best("rmse")
         
         # Refit to full training data
         final_rf_wflow = rf_wflow %>% 
            finalize_workflow(best_rf)
         final_rf_wflow_perm = rf_wflow_perm %>% 
            finalize_workflow(best_rf_perm)
         
         final_rf_fit = final_rf_wflow %>%
            fit(train_df) 
         final_rf_fit_perm = final_rf_wflow_perm %>%
            fit(train_df)
         
         rf_test_pred = predict(final_rf_fit, test_df)
         rf_test_pred_perm = predict(final_rf_fit_perm, test_df)
         
         
         get_fit_indices(test_df, rf_test_pred, y = "actual")
         get_fit_indices(test_df, rf_test_pred_perm, y = "perm")
         
         temp_output = tibble(
            num_repeat = num_repeat,
            k = k_outer,
            algorithm = "rf",
            fit_indices = get_fit_indices(test_df, rf_test_pred, y = "actual"),
            fit_indices_perm = get_fit_indices(test_df, rf_test_pred_perm, y = "perm"),
            )
         
         
         
         
         ## Elastic net regression
         
         
         
         
      }
      
   }
   
   
   # Define nested CV cycle
   set.seed(runningSeed)
   nested_cv_splits = nested_cv(df, 
                                outside = vfold_cv(v = k_outer, repeats = 5), 
                                inside = vfold_cv(v = k_inner))
   runningSeed = runningSeed + 1

   
   output = tibble()
   
   object = nested_cv_splits$splits[[1]]
   
   

   
   ## Create workflow
   rf_wflow = 
      workflow() %>% 
      add_model(rf_spec) %>% 
      add_recipe(rec)
   
   
   
   rf_fit = rf_wflow %>% 
      tune_grid(
         resamples = train_folds,
         grid = rf_grid,
         metrics = metric_set(rmse)
      )
   
   ## Predict
   best_rf <- rf_fit %>%
      select_best("rmse")
   
   final_rf_wflow = rf_wflow %>% 
      finalize_workflow(best_rf)
   
   final_rf_fit = 
      final_rf_wflow %>%
      fit(train_df) 
   
   rf_testing_pred = predict(final_rf_fit, test_df)
   
   rmse_vec(truth = test_df[, y], estimate = rf_testing_pred$.pred)
   rsq_trad_vec(truth = test_df[, y], estimate = rf_testing_pred$.pred)
}

# create_recipe------------------------------------------

create_recipe = function(train, y = "actual")  {
   
   if(y == "actual") {
      rec = recipe(y ~ ., data = train) %>%
         update_role(y_perm, new_role = "id variable")
   } else if(y == "permuted") {
      y_name_old = y
      rec = recipe(y_perm ~ ., data = train) %>%
         update_role(y, new_role = "id variable")
   } else   {stop("Please define y as 'actual' or 'permuted'")}
   rec = rec %>%
      step_dummy(all_nominal_predictors()) %>%
      step_zv(all_predictors()) %>%
      step_normalize(all_numeric_predictors()) %>%
      step_impute_knn(all_predictors())
   
   return(rec)
}


# get_fit_indices----------------------------------------

get_fit_indices = function(df, pred, y = "actual")   {
   output = c(
      RMSE = rmse_vec(truth = df[, ifelse(y == "actual", "y", "y_perm")], 
                      estimate = pred$.pred),
      R2 = rsq_vec(truth = df[, ifelse(y == "actual", "y", "y_perm")], 
                        estimate = pred$.pred)
   )
   return(output)
}



# rf_rmse------------------------------------------------

rf_rmse <- function(object, cost = 1) {
   y_col = ncol(object$data)
   
   # Create model specifications
   rf_spec = rand_forest(mtry = tune(), trees = 1000, min_n = tune()) %>%
      set_engine("ranger") %>%
      set_mode("regression")
   
   rf_grid = grid_regular(
      levels = 5,
      finalize(mtry(), analysis(object) %>% dplyr::select(-y)), 
      min_n())
   
   # Create Recipe
   rec = recipe(as.formula(paste0(y, " ~ .")), data = analysis(object)) %>%
      step_dummy(all_nominal_predictors()) %>%
      step_zv(all_predictors()) %>%
      step_normalize(all_numeric_predictors()) %>%
      step_impute_knn(all_predictors())
   
   # Create workflow
   rf_wflow = 
      workflow() %>% 
      add_model(rf_spec) %>% 
      add_recipe(rec)
   
   model = 
      rand_forest(mode = "regression", cost = cost) %>% 
      set_engine("ranger") %>% 
      set_mode("regression") %>%
      fit(y ~ ., data = analysis(object))
   
   holdout_pred = 
      predict(model, assessment(object) %>% dplyr::select(-y)) %>% 
      bind_cols(assessment(object) %>% dplyr::select(y))
   
   rmse(holdout_pred, truth = y, estimate = .pred)$.estimate
}

# In some case, we want to parameterize the function over the tuning parameter:
rmse_wrapper <- function(cost, object) svm_rmse(object, cost)

# `object` will be an `rsplit` object for the bootstrap samples
tune_over_cost <- function(object) {
   tibble(cost = 2 ^ seq(-2, 8, by = 1)) %>% 
      mutate(RMSE = map_dbl(cost, rmse_wrapper, object = object))
}

# `object` is an `rsplit` object in `results$inner_resamples` 
summarize_tune_results <- function(object) {
   # Return row-bound tibble that has the 25 bootstrap results
   map_df(object$splits, tune_over_cost) %>%
      # For each value of the tuning parameter, compute the 
      # average RMSE which is the inner bootstrap estimate. 
      group_by(cost) %>%
      summarize(mean_RMSE = mean(RMSE, na.rm = TRUE),
                n = length(RMSE),
                .groups = "drop")
}


plan(multisession)

tuning_results <- future_map(results$inner_resamples, summarize_tune_results) 


# fit_inner_cv_wflow-------------------------------------

fit_inner_cv_wflow <- function(
   train,
   recipe_,
   model,
   tune_grid
) {
   ## Create workflow
   wflow = 
      workflow() %>% 
      add_model(model) %>% 
      add_recipe(recipe_)
   
   # Fit
   fit = wflow %>% 
      tune_grid(
         resamples = split,
         grid = tune_grid,
         metrics = metric_set(rmse)
      )
   
   best_rf <- rf_fit %>%
      select_best("rmse")
   
   final_rf_wflow = rf_wflow %>% 
      finalize_workflow(best_rf)
   
   final_rf_fit = 
      final_rf_wflow %>%
      fit(train_df) 
   
   
}


# nested.cv----------------------------------------------

nested.cv <- function(
   data, 
   x,
   y,
   k_outer = 10, # number of outer CV folds
   k_inner = 10, # number of inner CV folds
   num_repeats = 1, # number of CV repeats
   runGLMnet = TRUE,
   runRF = TRUE,
   runKNN = TRUE,
   seed = 1 # set Seed for analyses
)  {
   
   ## Subset data to relevant variables only
   data = data[, c(x, y)]
   
   ## Set running seed, which is incremented by 1 every time it is used
   runningSeed = seed
   
   ## Create folds
   obs.all = nrow(data)
   set.seed(runningSeed)
   fold.pool.outer = rep_len(sample(1:k_outer), length.out = obs.all)
   runningSeed = runningSeed + 1
   
   ## Define fitControl object for caret
   fitControl = trainControl(method = "cv",
                             number = k_inner)
   
   ## Define used models
   used.models = c()
   if(runGLMnet)  {used.models = c(used.models, "glmnet")}
   if(runRF)  {used.models = c(used.models, "rf")}
   if(runKNN)  {used.models = c(used.models, "knn")}
   
   
   ## Define tuneGrids
   glmnet.tuneGrid = expand.grid(alpha = seq(from = 0, to = 1, by = 0.1),
                                 lambda = seq(from = 0, to = 1, by = 0.1))
   rf.tuneGrid = data.frame(mtry = unique(round(seq(from = 2, to = length(x), 
                                                    length.out = 5))))
   knn.tuneGrid = expand.grid(k = 1:25)
   
   
   ## Create output data.frame for fit statistics
   fit.stats = expand.grid(num_repeat = 1:num_repeats,
                           k = 1:k_outer,
                           model = used.models,
                           type = c("pred", "perm"))
   fit.stats[, c("RMSE", "Rsquared", "MAE", "alpha", "lambda", "mtry", "k.knn")] = NA
   
   
   ## Create output data.frame for variable importance
   varImp.stats = fit.stats[, c("num_repeat", "k", "model", "type")]
   varImp.stats[, sort(x)] = NA
   
   
   ## Create output data.frame for individual predictions
   # Save rowID
   data$rowID = paste0("r", 1:nrow(data))
   
   # Create data.frame
   pred.stats = expand.grid(num_repeat = 1:num_repeats,
                            k = 1:k_outer,
                            model = used.models,
                            type = c("pred", "perm"))
   pred.stats[, unique(data$rowID)] = NA
   
   
   
   ## Run repeats
   for(repeats in 1:num_repeats)  {
      
      ## Print repeat number to index progress
      cat(paste0("\nRepeat no.:\t\t", repeats))
      
      ## Y is permuted for permutation test
      set.seed(runningSeed)
      data[, "y_perm"] = sample(data[, y], replace = FALSE)
      runningSeed = runningSeed + 1
      
      ## Assign folds
      set.seed(runningSeed)
      data$fold.outer = sample(fold.pool.outer)
      runningSeed = runningSeed + 1
      
      ## Run outer CV
      for(outer in 1:k_outer) {
         
         ## Print outer fold number to index progress
         cat(paste0("\n---Outer CV-Fold no.:\t", outer))
         
         ## Divide data into train and test sets
         train = data[data$fold.outer != outer,]
         test = data[data$fold.outer == outer,]
         
         # 1 glmnet-------------------------------------------
         
         if(runGLMnet)  {
            
            ## Run inner CV
            set.seed(runningSeed)
            glmnet.fit.pred = train(x = train[,x], y = train[,y], 
                               method = "glmnet", metric = "RMSE", 
                               trControl = fitControl,
                               tuneGrid = glmnet.tuneGrid)
            runningSeed = runningSeed + 1
            
            set.seed(runningSeed)
            glmnet.fit.perm = train(x = train[,x], y = train[,"y_perm"], 
                                    method = "glmnet", metric = "RMSE", 
                                    trControl = fitControl,
                                    tuneGrid = glmnet.tuneGrid)
            runningSeed = runningSeed + 1
            
            ## Save index rows to save fit.stats output
            fit.stats.index.pred = which(fit.stats$num_repeat == repeats & 
                                            fit.stats$k == outer & fit.stats$type == "pred" &
                                            fit.stats$model == "glmnet")
            fit.stats.index.perm = which(fit.stats$num_repeat == repeats & 
                                            fit.stats$k == outer & fit.stats$type == "perm" &
                                            fit.stats$model == "glmnet")
            
            ## Save final tuning parameters
            fit.stats[fit.stats.index.pred, c("alpha", "lambda")] = glmnet.fit.pred$bestTune
            fit.stats[fit.stats.index.perm, c("alpha", "lambda")] = glmnet.fit.perm$bestTune
            
            
            ## Predict in independent test set
            glmnet.preds.pred = predict(glmnet.fit.pred, newdata = test)
            glmnet.preds.perm = predict(glmnet.fit.perm, newdata = test)
            
            # Define pred.stats indices to save prediction results
            pred.stats.index.pred = which(pred.stats$num_repeat == repeats & pred.stats$model == "glmnet" &
                                             pred.stats$k == outer & pred.stats$type == "pred")
            pred.stats.index.perm = which(pred.stats$num_repeat == repeats & pred.stats$model == "glmnet" &
                                             pred.stats$k == outer & pred.stats$type == "perm")
            
            # Save predictions
            pred.stats[pred.stats.index.pred, test$rowID] = glmnet.preds.pred
            pred.stats[pred.stats.index.perm, test$rowID] = glmnet.preds.perm
            
            # Save fit statistics
            fit.stats[fit.stats.index.pred, c("RMSE", "Rsquared", "MAE")] = 
               postResample(test[, y], glmnet.preds.pred)
            fit.stats[fit.stats.index.perm, c("RMSE", "Rsquared", "MAE")] = 
               postResample(test[, y], glmnet.preds.perm)
            
            ## Save variable importance
            glmnet.varImp.pred = varImp(glmnet.fit.pred)$importance
            glmnet.varImp.perm = varImp(glmnet.fit.perm)$importance
            
            # save variable label and order by label
            glmnet.varImp.pred$vars = row.names(glmnet.varImp.pred)
            glmnet.varImp.pred = arrange(glmnet.varImp.pred, vars)
            glmnet.varImp.perm$vars = row.names(glmnet.varImp.perm)
            glmnet.varImp.perm = arrange(glmnet.varImp.perm, vars)
            
            # Save importance values
            varImp.stats[fit.stats.index.pred, sort(x)] = glmnet.varImp.pred$Overall
            varImp.stats[fit.stats.index.perm, sort(x)] = glmnet.varImp.perm$Overall
            
         }
         
         
         # 2 rf-----------------------------------------------
         
         if(runRF == TRUE) {
            
            ## Run inner CV
            set.seed(runningSeed)
            rf.fit.pred = train(x = train[,x], y = train[,y], 
                                method = "rf", metric = "RMSE", 
                                trControl = fitControl,
                                tuneGrid = rf.tuneGrid)
            runningSeed = runningSeed + 1
            
            set.seed(runningSeed)
            rf.fit.perm = train(x = train[,x], y = train[,"y_perm"], 
                                method = "rf", metric = "RMSE", 
                                trControl = fitControl,
                                tuneGrid = rf.tuneGrid)
            runningSeed = runningSeed + 1
            
            ## Save index row to save fit.stats output
            fit.stats.index.pred = which(fit.stats$num_repeat == repeats & 
                                            fit.stats$k == outer & fit.stats$type == "pred" &
                                            fit.stats$model == "rf")
            fit.stats.index.perm = which(fit.stats$num_repeat == repeats & 
                                            fit.stats$k == outer & fit.stats$type == "perm" &
                                            fit.stats$model == "rf")
            
            ## Save final tuning parameters
            fit.stats[fit.stats.index.pred, "mtry"] = rf.fit.pred$finalModel$tuneValue
            fit.stats[fit.stats.index.perm, "mtry"] = rf.fit.perm$finalModel$tuneValue
            
            ## Predict in independent test set
            rf.preds.pred = predict(rf.fit.pred, newdata = test)
            rf.preds.perm = predict(rf.fit.perm, newdata = test)
            
            # Define pred.stats indices to save prediction results
            pred.stats.index.pred = which(pred.stats$num_repeat == repeats & pred.stats$model == "rf" &
                                             pred.stats$k == outer & pred.stats$type == "pred")
            pred.stats.index.perm = which(pred.stats$num_repeat == repeats & pred.stats$model == "rf" &
                                             pred.stats$k == outer & pred.stats$type == "perm")
            
            
            # Save predictions
            pred.stats[pred.stats.index.pred, test$rowID] = rf.preds.pred
            pred.stats[pred.stats.index.perm, test$rowID] = rf.preds.perm
            
            
            # Save fit statistics
            fit.stats[fit.stats.index.pred, c("RMSE", "Rsquared", "MAE")] = 
               postResample(test[, y], rf.preds.pred)
            fit.stats[fit.stats.index.perm, c("RMSE", "Rsquared", "MAE")] = 
               postResample(test[, y], rf.preds.perm)
            
            ## Save variable importance
            rf.varImp.pred = varImp(rf.fit.pred)$importance
            rf.varImp.perm = varImp(rf.fit.perm)$importance
            
            # save variable label and order by label
            rf.varImp.pred$vars = row.names(rf.varImp.pred)
            rf.varImp.pred = arrange(rf.varImp.pred, vars)
            
            rf.varImp.perm$vars = row.names(rf.varImp.perm)
            rf.varImp.perm = arrange(rf.varImp.perm, vars)
            
            # Save importance values
            varImp.stats[fit.stats.index.pred, sort(x)] = rf.varImp.pred$Overall
            varImp.stats[fit.stats.index.perm, sort(x)] = rf.varImp.perm$Overall
            
         }
         
         # 3 knn----------------------------------------------
         
         if(runKNN == TRUE) {
            
            ## Run inner CV
            set.seed(runningSeed)
            knn.fit.pred = train(x = train[,x], y = train[,y], 
                                 method = "knn", metric = "RMSE", 
                                 trControl = fitControl,
                                 tuneGrid = knn.tuneGrid)
            runningSeed = runningSeed + 1
            
            set.seed(runningSeed)
            knn.fit.perm = train(x = train[,x], y = train[,"y_perm"], 
                                 method = "knn", metric = "RMSE", 
                                 trControl = fitControl,
                                 tuneGrid = knn.tuneGrid)
            runningSeed = runningSeed + 1
            
            ## Save index row to save fit.stats output
            fit.stats.index.pred = which(fit.stats$num_repeat == repeats & 
                                            fit.stats$k == outer & fit.stats$type == "pred" &
                                            fit.stats$model == "knn")
            fit.stats.index.perm = which(fit.stats$num_repeat == repeats & 
                                            fit.stats$k == outer & fit.stats$type == "perm" &
                                            fit.stats$model == "knn")
            
            ## Save final tuning parameters
            fit.stats[fit.stats.index.pred, "k.knn"] = knn.fit.pred$finalModel$tuneValue
            fit.stats[fit.stats.index.perm, "k.knn"] = knn.fit.perm$finalModel$tuneValue
            
            ## Predict in independent test set
            knn.preds.pred = predict(knn.fit.pred, newdata = test)
            knn.preds.perm = predict(knn.fit.perm, newdata = test)
            
            # Define pred.stats indices to save prediction results
            pred.stats.index.pred = which(pred.stats$num_repeat == repeats & pred.stats$model == "knn" &
                                             pred.stats$k == outer & pred.stats$type == "pred")
            pred.stats.index.perm = which(pred.stats$num_repeat == repeats & pred.stats$model == "knn" &
                                             pred.stats$k == outer & pred.stats$type == "perm")
            
            # Save predictions
            pred.stats[pred.stats.index.pred, test$rowID] = knn.preds.pred
            pred.stats[pred.stats.index.perm, test$rowID] = knn.preds.perm
            
            
            # Save fit statistics
            fit.stats[fit.stats.index.pred, c("RMSE", "Rsquared", "MAE")] = 
               postResample(test[, y], knn.preds.pred)
            fit.stats[fit.stats.index.perm, c("RMSE", "Rsquared", "MAE")] = 
               postResample(test[, y], knn.preds.perm)
            
            ## Save variable importance
            knn.varImp.pred = varImp(knn.fit.pred)$importance
            knn.varImp.perm = varImp(knn.fit.perm)$importance
            
            # save variable label and order by label
            knn.varImp.pred$vars = row.names(knn.varImp.pred)
            knn.varImp.pred = arrange(knn.varImp.pred, vars)
            
            knn.varImp.perm$vars = row.names(knn.varImp.perm)
            knn.varImp.perm = arrange(knn.varImp.perm, vars)
            
            # Save importance values
            varImp.stats[fit.stats.index.pred, sort(x)] = knn.varImp.pred$Overall
            varImp.stats[fit.stats.index.perm, sort(x)] = knn.varImp.perm$Overall
            
         }
         
         
      }
      
   }
   
   ## Subset fit.stats to fit.stats with available data
   for(i in colnames(fit.stats)[7:ncol(fit.stats)])   {
      if(sum(is.na(fit.stats[, i])) == nrow(fit.stats))  {fit.stats[, i] = NULL}
   }
   

   ## Collate output for fit statistics and variable importance in list
   output = list(fit = fit.stats, 
                 varImp = varImp.stats, 
                 pred = pred.stats)
   
   return(output)
   
}



# getmode------------------------------------------------

getmode <- function(v) {
   
   uniqv = unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
   
}




# t.test_loop--------------------------------------------

t.test_loop = function(
   x, # fit statistics dataset
   round_dec = 3
)  {
   
   ## Create empty output data.frame
   output = expand.grid(feature.set = unique(x$feature.set),
                        model = unique(x$model))
   output[, c("RMSE_mean", "RMSE_sd", "R2_median", "R2_25quantile", "R2_75quantile",
              "RMSE_delta_mean", "RMSE_delta_sd", "RMSE_delta_median", "RMSE_delta_5quantile", 
              "RMSE_delta_95quantile", "RMSE", "R2", "RMSE_delta", "tval", "pval", "quant_pval")] = NA
   
   ## Loop over combinations
   for(i in unique(output$feature.set))   {
      for(j in unique(output$model))   {
         
         ## Subset data
         x.pred = filter(x, feature.set == i & model == j & type == "pred")
         x.perm = filter(x, feature.set == i & model == j & type == "perm")
         
         ## Conduct paired t-test
         t.test.fit = t.test(x.pred$RMSE, x.perm$RMSE, paired = TRUE)
         
         ## Write output data
         output.index = which(output$feature.set == i & output$model == j)
         
         output[output.index, "RMSE_mean"] = mean(x.pred$RMSE)
         output[output.index, "RMSE_sd"] = sd(x.pred$RMSE)
         output[output.index, "R2_median"] = median(x.pred$Rsquared, na.rm = TRUE)
         output[output.index, "R2_25quantile"] = quantile(x.pred$Rsquared, probs = 0.25, na.rm = TRUE)
         output[output.index, "R2_75quantile"] = quantile(x.pred$Rsquared, probs = 0.75, na.rm = TRUE)
         
         RMSE_delta = x.pred$RMSE - x.perm$RMSE
         output[output.index, "RMSE_delta_mean"] = mean(RMSE_delta)
         output[output.index, "RMSE_delta_sd"] = sd(RMSE_delta)
         output[output.index, "RMSE_delta_median"] = median(RMSE_delta)
         output[output.index, c("RMSE_delta_5quantile", "RMSE_delta_95quantile")] = 
            quantile(RMSE_delta, probs = c(0.05, 0.95))
         output[output.index, "quant_pval"] = (1 - ecdf(RMSE_delta)(0))
         
         output[output.index, "tval"] = t.test.fit$statistic
         output[output.index, "pval"] = t.test.fit$p.value
         
      }
   }
   
   
   output[, "RMSE"] = with(output, paste0(round(RMSE_mean, round_dec), " (", 
                                          round(RMSE_sd, round_dec), ")"))
   output[, "R2"] = with(output, paste0(round(R2_median, round_dec), " (", 
                                        round(R2_25quantile, round_dec), "-", 
                                        round(R2_75quantile, round_dec), ")"))
   
   output[, "RMSE_delta"] = with(output, paste0(round(RMSE_delta_mean, round_dec), " (", 
                                          round(RMSE_delta_sd, round_dec), ")"))
   
   return(output)
   
}



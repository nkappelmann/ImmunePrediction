# -------------------------------------------------------
# Functions----------------------------------------------
# -------------------------------------------------------


# nested.cv----------------------------------------------

nested.cv <- function(
   data, 
   x,
   y,
   k.outer = 10, # number of outer CV folds
   k.inner = 10, # number of inner CV folds
   num_repeats = 1, # number of CV repeats
   perm.test = FALSE, # Should y be shuffled to obtain a permutation distribution?
   runGLMnet = TRUE,
   runRF = TRUE,
   runSVM = TRUE,
   runBART = TRUE
)  {
   
   ## Subset data to relevant variables only
   data = data[, c(x, y)]
   
   ## Create folds
   obs.all = nrow(data)
   fold.pool.outer = rep_len(sample(1:k.outer), length.out = obs.all)
   
   ## Define fitControl object for caret
   fitControl = trainControl(method = "cv",
                             number = k.inner)
   
   ## Define tuneGrid
   glmnet.tuneGrid = expand.grid(alpha = seq(from = 0, to = 1, by = 0.2),
                                 lambda = seq(from = 0, to = 1, by = 0.2))
   rf.tuneGrid = data.frame(mtry = unique(round(seq(from = 2, to = length(x), length.out = 5))))
   bart.tuneGrid = expand.grid(num_trees = c(10, 15, 20, 50), 
                               k = 2, 
                               alpha = 0.95, 
                               beta = 2, 
                               nu = 3)
   
   ## Create output data.frame for fit statistics
   fit.stats = expand.grid(num_repeat = 1:num_repeats,
                           k = 1:k.outer,
                           model = c("glmnet", "rf", "bart"))
   fit.stats[, c("RMSE", "Rsquared", "MAE", "alpha", "lambda", "mtry", 
                 "num_trees", "k.bart", "beta", "nu")] = NA
   
   
   ## Create output data.frame for variable importance
   varImp.stats = fit.stats[, c("num_repeat", "k", "model")]
   varImp.stats[, sort(x)] = NA
   
   
   ## Create output data.frame for individual predictions
   # Save rowID
   data$rowID = paste0("r", 1:nrow(data))
   
   # Create data.frame
   pred.stats = expand.grid(num_repeat = 1:num_repeats,
                            model = c("glmnet", "rf", "bart"))
   pred.stats[, unique(data$rowID)] = NA
   
   
   
   ## Run repeats
   for(repeats in 1:num_repeats)  {
      
      ## Print repeat number to index progress
      cat(paste0("\nRepeat no.:\t\t", repeats))
      
      ## Permute Y if a permutation test is indicated
      if(perm.test == TRUE)   {
         data[, y] = sample(data[, y], replace = FALSE)
      }
      
      ## Assign folds
      data$fold.outer = sample(fold.pool.outer)
      
      ## Run outer CV
      for(outer in 1:k.outer) {
         
         ## Print outer fold number to index progress
         cat(paste0("\n---Outer CV-Fold no.:\t", outer))
         
         ## Divide data into train and test sets
         train = data[data$fold.outer != outer,]
         test = data[data$fold.outer == outer,]
         
         # 1 glmnet-------------------------------------------
         
         if(runGLMnet)  {
            
            ## Run inner CV
            glmnet.fit = train(x = train[,x], y = train[,y], 
                               method = "glmnet", metric = "RMSE", 
                               trControl = fitControl,
                               tuneGrid = glmnet.tuneGrid)
            
            ## Save index row to save fit.stats output
            fit.stats.index = which(fit.stats$num_repeat == repeats & 
                                       fit.stats$k == outer & 
                                       fit.stats$model == "glmnet")
            
            ## Save final tuning parameters
            fit.stats[fit.stats.index, c("alpha", "lambda")] = glmnet.fit$bestTune
            
            
            
            ## Predict in independent test set
            glmnet.preds = predict(glmnet.fit, newdata = test)
            
            ## Create MANUAL PREDICTIONS
            glmnet.coefs = coef(glmnet.fit$finalModel, s = glmnet.fit$bestTune$lambda)
            manual.preds = as.vector(glmnet.coefs[1,] + glmnet.coefs[2,]*test$t0_bdi_std + 
                                        glmnet.coefs[3,]*test$sex_std + 
                                        glmnet.coefs[4,]*test$age_std + 
                                        glmnet.coefs[5,]*test$BMI_std)
            
            # Save predictions
            pred.stats[pred.stats$num_repeat == repeats & pred.stats$model == "glmnet", 
                       test$rowID] = glmnet.preds
            
            # Save fit statistics
            fit.stats[fit.stats.index, c("RMSE", "Rsquared", "MAE")] = 
               postResample(test[, y], glmnet.preds)
            
            ## Save variable importance
            glmnet.varImp = varImp(glmnet.fit)$importance
            
            # save variable label and order by label
            glmnet.varImp$vars = row.names(glmnet.varImp)
            glmnet.varImp = arrange(glmnet.varImp, vars)
            
            # Save importance values
            varImp.stats[fit.stats.index, sort(x)] = glmnet.varImp$Overall
            
         }
         
         
         # 2 rf-----------------------------------------------
         
         if(runRF == TRUE) {
            
            ## Run inner CV
            rf.fit = train(x = train[,x], y = train[,y], 
                           method = "rf", metric = "RMSE", 
                           trControl = fitControl,
                           tuneGrid = rf.tuneGrid
            )
            
            ## Save index row to save fit.stats output
            fit.stats.index = which(fit.stats$num_repeat == repeats & 
                                       fit.stats$k == outer & 
                                       fit.stats$model == "rf")
            
            ## Save final tuning parameters
            fit.stats[fit.stats.index, "mtry"] = rf.fit$finalModel$tuneValue
            
            ## Predict in independent test set
            rf.preds = predict(rf.fit, newdata = test)
            
            # Save predictions
            pred.stats[pred.stats$num_repeat == repeats & pred.stats$model == "rf", 
                       test$rowID] = rf.preds
            
            
            # Save fit statistics
            fit.stats[fit.stats.index, c("RMSE", "Rsquared", "MAE")] = 
               postResample(test[, y], rf.preds)
            
            ## Save variable importance
            rf.varImp = varImp(rf.fit)$importance
            
            # save variable label and order by label
            rf.varImp$vars = row.names(rf.varImp)
            rf.varImp = arrange(rf.varImp, vars)
            
            # Save importance values
            varImp.stats[fit.stats.index, sort(x)] = rf.varImp$Overall
            
         }
         
         # 3 BARTmachine--------------------------------------
         
         if(runBART == TRUE)  {
            
            ## Run inner CV
            bart.fit = train(x = train[,x], y = train[,y], 
                             method = "bartMachine", metric = "RMSE", 
                             trControl = fitControl,
                             tuneGrid = bart.tuneGrid
            )
            
            ## Save index row to save fit.stats output
            fit.stats.index = which(fit.stats$num_repeat == repeats & 
                                       fit.stats$k == outer & 
                                       fit.stats$model == "bart")
            
            ## Save final tuning parameters
            fit.stats[fit.stats.index, c("num_trees", "k.bart", "alpha", "beta", "nu")] = 
               bart.fit$finalModel$tuneValue
            
            ## Predict in independent test set
            bart.preds = predict(bart.fit, newdata = test)
            
            # Save predictions
            pred.stats[pred.stats$num_repeat == repeats & pred.stats$model == "bart", 
                       test$rowID] = bart.preds
            
            # Save fit statistics
            fit.stats[fit.stats.index, c("RMSE", "Rsquared", "MAE")] = 
               postResample(test[, y], bart.preds)
            
            ## Save variable importance
            bart.varImp = varImp(bart.fit)$importance
            
            # save variable label and order by label
            bart.varImp$vars = row.names(bart.varImp)
            bart.varImp = arrange(bart.varImp, vars)
            
            # Save importance values
            varImp.stats[fit.stats.index, sort(x)] = bart.varImp$Overall
            
         }
         
      }
      
   }
   
   ## Aggregate pred.stats across models using a voting system
   pred.aggregate.stats = pred.aggregate(preds = pred.stats)
   
   
   ## Collate output for fit statistics and variable importance in list
   output = list(fit = fit.stats, 
                 varImp = varImp.stats, 
                 pred = pred.stats,
                 pred.aggregate = pred.aggregate.stats)
   
   return(output)
   
}



# pred.aggregate-----------------------------------------

pred.aggregate <- function(
   preds
   ) {
   
   ## Get rowindex
   rowindex = paste0("r", 1:(ncol(preds) - 2))
   
   ## Remove prediction output for models that were not run
   preds = na.omit(preds)
   
   ## Round predictions
   preds[, rowindex] = round(preds[, rowindex])
   
   ## transpose
   tpreds = t(preds[, rowindex])
   
   ## Get mode
   pred.aggregate = apply(tpreds, 1, getmode) %>% as.vector()
   
   return(pred.aggregate)
}


# getmode------------------------------------------------

getmode <- function(v) {
   
   uniqv = unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
   
}

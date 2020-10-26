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
   runKNN = TRUE,
   runNNET = FALSE,
   runSVM = FALSE,
   runBART = FALSE
)  {
   
   ## Subset data to relevant variables only
   data = data[, c(x, y)]
   
   ## Create folds
   obs.all = nrow(data)
   fold.pool.outer = rep_len(sample(1:k.outer), length.out = obs.all)
   
   ## Define fitControl object for caret
   fitControl = trainControl(method = "cv",
                             number = k.inner)
   
   ## Define used models
   used.models = c()
   if(runGLMnet)  {used.models = c(used.models, "glmnet")}
   if(runRF)  {used.models = c(used.models, "rf")}
   if(runKNN)  {used.models = c(used.models, "knn")}
   if(runSVM)  {used.models = c(used.models, "svm")}
   if(runNNET)  {used.models = c(used.models, "nnet")}
   if(runBART)  {used.models = c(used.models, "bart")}
   
   ## Define tuneGrids
   glmnet.tuneGrid = expand.grid(alpha = seq(from = 0, to = 1, by = 0.2),
                                 lambda = seq(from = 0, to = 1, by = 0.2))
   rf.tuneGrid = data.frame(mtry = unique(round(seq(from = 2, to = length(x), length.out = 5))))
   knn.tuneGrid = expand.grid(k = 1:10)
   svm.tuneGrid = expand.grid(degree = 1:3, 
                              scale = c(0.001, 0.01, 0.1), 
                              C = c(0.25, 0.5, 1))
   nnet.tuneGrid = expand.grid(size = seq(from = 1, to = 10, by = 2),
                               decay = seq(from = 0.1, to = 0.5, by = 0.1))
   bart.tuneGrid = expand.grid(num_trees = c(10, 15, 20, 50), 
                               k = 2, 
                               alpha = 0.95, 
                               beta = 2, 
                               nu = 3)
   
   ## Create output data.frame for fit statistics
   fit.stats = expand.grid(num_repeat = 1:num_repeats,
                           k = 1:k.outer,
                           model = used.models)
   fit.stats[, c("RMSE", "Rsquared", "MAE", "alpha", "lambda", "mtry", "k.knn", 
                 "degree", "scale", "C", "size", "decay",
                 "num_trees", "k.bart", "beta", "nu")] = NA
   
   
   ## Create output data.frame for variable importance
   varImp.stats = fit.stats[, c("num_repeat", "k", "model")]
   varImp.stats[, sort(x)] = NA
   
   
   ## Create output data.frame for individual predictions
   # Save rowID
   data$rowID = paste0("r", 1:nrow(data))
   
   # Create data.frame
   pred.stats = expand.grid(num_repeat = 1:num_repeats,
                            model = used.models)
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
         
         # 3 knn----------------------------------------------
         
         if(runKNN == TRUE) {
            
            ## Run inner CV
            knn.fit = train(x = train[,x], y = train[,y], 
                           method = "knn", metric = "RMSE", 
                           trControl = fitControl,
                           tuneGrid = knn.tuneGrid
            )
            
            ## Save index row to save fit.stats output
            fit.stats.index = which(fit.stats$num_repeat == repeats & 
                                       fit.stats$k == outer & 
                                       fit.stats$model == "knn")
            
            ## Save final tuning parameters
            fit.stats[fit.stats.index, "k.knn"] = knn.fit$finalModel$tuneValue
            
            ## Predict in independent test set
            knn.preds = predict(knn.fit, newdata = test)
            
            # Save predictions
            pred.stats[pred.stats$num_repeat == repeats & pred.stats$model == "knn", 
                       test$rowID] = knn.preds
            
            
            # Save fit statistics
            fit.stats[fit.stats.index, c("RMSE", "Rsquared", "MAE")] = 
               postResample(test[, y], knn.preds)
            
            ## Save variable importance
            knn.varImp = varImp(knn.fit)$importance
            
            # save variable label and order by label
            knn.varImp$vars = row.names(knn.varImp)
            knn.varImp = arrange(knn.varImp, vars)
            
            # Save importance values
            varImp.stats[fit.stats.index, sort(x)] = knn.varImp$Overall
            
         }
         
         # 4 svm----------------------------------------------
         
         if(runSVM)  {
            
            ## Run inner CV
            svm.fit = train(x = train[,x], y = train[,y], 
                               method = "svmPoly", metric = "RMSE", 
                               trControl = fitControl,
                               tuneGrid = svm.tuneGrid)
            
            ## Save index row to save fit.stats output
            fit.stats.index = which(fit.stats$num_repeat == repeats & 
                                       fit.stats$k == outer & 
                                       fit.stats$model == "svm")
            
            ## Save final tuning parameters
            fit.stats[fit.stats.index, c("degree", "scale", "C")] = svm.fit$bestTune
            
            
            ## Predict in independent test set
            svm.preds = predict(svm.fit, newdata = test)
            
            
            # Save predictions
            pred.stats[pred.stats$num_repeat == repeats & pred.stats$model == "svm", 
                       test$rowID] = svm.preds
            
            # Save fit statistics
            fit.stats[fit.stats.index, c("RMSE", "Rsquared", "MAE")] = 
               postResample(test[, y], svm.preds)
            
            ## Save variable importance
            svm.varImp = varImp(svm.fit)$importance
            
            # save variable label and order by label
            svm.varImp$vars = row.names(svm.varImp)
            svm.varImp = arrange(svm.varImp, vars)
            
            # Save importance values
            varImp.stats[fit.stats.index, sort(x)] = svm.varImp$Overall
            
         }
         
         
         
         # 5 nnet---------------------------------------------
         
         if(runNNET)  {
            
            ## Run inner CV
            nnet.fit = train(x = train[,x], y = train[,y], 
                            method = "nnet", metric = "RMSE", 
                            trControl = fitControl,
                            tuneGrid = nnet.tuneGrid,
                            trace = FALSE)
            
            ## Save index row to save fit.stats output
            fit.stats.index = which(fit.stats$num_repeat == repeats & 
                                       fit.stats$k == outer & 
                                       fit.stats$model == "nnet")
            
            ## Save final tuning parameters
            fit.stats[fit.stats.index, c("size", "decay")] = nnet.fit$bestTune
            
            
            ## Predict in independent test set
            nnet.preds = predict(nnet.fit, newdata = test)
            
            
            # Save predictions
            pred.stats[pred.stats$num_repeat == repeats & pred.stats$model == "nnet", 
                       test$rowID] = nnet.preds
            
            # Save fit statistics
            fit.stats[fit.stats.index, c("RMSE", "Rsquared", "MAE")] = 
               postResample(test[, y], nnet.preds)
            
            ## Save variable importance
            nnet.varImp = varImp(nnet.fit)$importance
            
            # save variable label and order by label
            nnet.varImp$vars = row.names(nnet.varImp)
            nnet.varImp = arrange(nnet.varImp, vars)
            
            # Save importance values
            varImp.stats[fit.stats.index, sort(x)] = nnet.varImp$Overall
            
         }
         
         
         # 6 BARTmachine--------------------------------------
         
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
   
   ## Subset fit.stats to fit.stats with available data
   for(i in colnames(fit.stats)[7:ncol(fit.stats)])   {
      if(sum(is.na(fit.stats[, i])) == nrow(fit.stats))  {fit.stats[, i] = NULL}
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





# aggregate.varImp---------------------------------------


aggregate.varImp <- function(
   obj # nested.cv output object
)   {
   
   ## Define input and output
   input = obj$varImp
   temp = data.frame(num_repeat = numeric(),
                       k = numeric(),
                       model = character(),
                       var = character(),
                       varImp = numeric())
   
   
   ## Fill varImp data.frame
   for(i in 4:ncol(input))   {
      for(j in 1:nrow(input))   {
         
         ## Create empty new row
         temp[(nrow(temp) + 1), ] = NA
         
         ## Fill general info
         temp[nrow(temp), c("num_repeat", "k")] = 
            input[j, c("num_repeat", "k")]
         temp[nrow(temp), "model"] = 
            as.character(input[j, "model"])
         
         ## Fill var and varImp info
         temp[nrow(temp), "var"] = 
            colnames(input)[i]
         temp[nrow(temp), "varImp"] = 
            input[j, i]
      }
      
   }
   
   ## Reduce data and get median info
   output = rank.aggregate.varImp(x = temp)
   
   
   return(output)
}


# rank.aggregate.varImp----------------------------------

rank.aggregate.varImp = function(
   x # input should be varImp data in long format as provided by aggregate.varImp()
)  {
   
   # Create ranking on variable importance
   varImp.ranking = with(x, by(varImp, var, median, na.rm = TRUE))
   varImp.ranking = data.frame(var = names(varImp.ranking),
                               varImp = as.vector(varImp.ranking)) %>%
      arrange(desc(varImp))
   
   
   # Create levels of individual data
   x$var = factor(x$var, levels = varImp.ranking$var)
   
   # Arrange data by varImp
   x = arrange(x, var)
   
   # Create output list
   output = list(ind = x,
                 ranking = varImp.ranking)
   
   return(output)
}


# ------------------------------------
# Pipeline Test-----------------------
# Housing Data------------------------
# ------------------------------------


# 1 Preparation-----------------------

## Load packages
library("tidyverse")
library("caret")
library("caretEnsemble")
library("AmesHousing")


## Load data
house = make_ames()
house = na.omit(house)

## Preprocess data
# Change to data.frame
data = as.data.frame(house)

# Remove columns with low variance
data = data[, -nearZeroVar(data)]

# Remove factor and character variables
data = data[, sapply(data, class) %in% c("integer", "numeric")]

# Standardise variables
data = scale(data) %>% as.data.frame()

## Set y & x
y = "Sale_Price"
x = colnames(data)[colnames(data) %in% y == FALSE]

## Reduce rows to subsample
data = data[sample(1:nrow(data), size = 120, replace = FALSE),]



# 2 Run ML pipeline-------------------

# Define function parameters to run code outside function

k.outer = 10 # number of outer CV folds
k.inner = 10 # number of inner CV folds
num_repeats = 1 # number of CV repeats
repeats = 1


## Create folds
obs.all = nrow(data)
fold.pool.outer = rep_len(sample(1:k.outer), length.out = obs.all)

## Define used models
used.models = c("glmnet", "rf", "knn")

## Define tuneGrids
glmnet.tuneGrid = expand.grid(alpha = seq(from = 0, to = 1, by = 0.2),
                              lambda = seq(from = 0, to = 1, by = 0.2))
rf.tuneGrid = data.frame(mtry = unique(round(seq(from = 2, to = length(x), 
                                                 length.out = 5))))
knn.tuneGrid = expand.grid(k = 1:25)


## Create output data.frame for individual predictions
# Save rowID
data$rowID = paste0("r", 1:nrow(data))


## Assign folds
data$fold.outer = sample(fold.pool.outer)


## Divide data into train and test sets
train = data[data$fold.outer != outer,]
test = data[data$fold.outer == outer,]

## Define fitControl object for caret
fitControl = trainControl(method = "cv",
                          number = k.inner,
                          savePredictions = "final",
                          index = createResample(train[,y], 
                                                 k.inner))


## Run individual algorithms in inner CV

# glmnet
glmnet.fit = train(x = train[,x], y = train[,y], 
                   method = "glmnet", metric = "RMSE", 
                   trControl = fitControl,
                   tuneGrid = glmnet.tuneGrid)


# random forest
rf.fit = train(x = train[,x], y = train[,y], 
               method = "rf", metric = "RMSE", 
               trControl = fitControl,
               tuneGrid = rf.tuneGrid)

# k-nearest neighbours
knn.fit = train(x = train[,x], y = train[,y], 
                method = "knn", metric = "RMSE", 
                trControl = fitControl,
                tuneGrid = knn.tuneGrid)


## Run ensemble model
model_list = caretList(
      x = train[,x], y = train[,y], 
      trControl = fitControl,
      methodList = c("glmnet", "rf", "knn"))


# Apply greedy algorithm 
greedy_ensemble = caretEnsemble(
      model_list, 
      metric = "RMSE",
      trControl = fitControl)
summary(greedy_ensemble)


## Extract predictions on outer test set
model_preds = lapply(model_list, predict, newdata = test)
ensemble_preds = predict(greedy_ensemble, newdata = test)

## Get intercorrelations between model predictions
modelCor(resamples(model_list))

## Get R^2 with test set values
cor(model_preds$glmnet, test$Sale_Price)^2
cor(model_preds$rf, test$Sale_Price)^2
cor(model_preds$knn, test$Sale_Price)^2
cor(ensemble_preds, test$Sale_Price)^2


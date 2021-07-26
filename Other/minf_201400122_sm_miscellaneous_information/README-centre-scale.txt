Solubility Challenge

Dataset was taken verbatim, with descriptors, from Neetika and James. 
After the ensemble methodology was finalised for this test, the dataset was modified slightly:
	NULL, ZeroVar, and NearZeroVar columns (Descriptors) removed. Leaving 122 molecules with 114(?) descriptors.
	See "SC_All_Descriptors.csv" in this directory for final used dataset. 

R:
	source("folders.r)
	input<-read.csv("SC_All_Descriptors.csv",row.names=1,header=T)
	source("method.r")
	
Methodology:
	Define string of 5 seed numbers. 
	Repeat following for each of those seeds:
		Establish 1 in 10 k-fold CV loop:
			Split dataset into training and testing (testing held for external validation). Split defined by seed & fold number. 
				Train each model with caret package, with 10 fold CV for training
					Trained with 90:10 split of training data into training/internal validation.
					Trained with various tuning parameters
					Repeated across 10 folds
					Averaged. Best tuning parameters used for final model. 
				Trained on full training set. 
				Ensemble Predictors (greedy, linear) generated from models. Weightings based on internal validation RMSE.
				Output of training and internal validation > output_fold_k.txt, error_fold_k.csv < includes the finalised parameters for each model, weightings for the greedy ensemble.
			All models inc. ensemble models tested against external validation set (i.e. data they've never seen)
			Predictions vs actual data outputted > Predictions_fold_k.csv
			Repeated for 10 kfold CV
		Repeated for 5 seed numbers
	End
	
The predictions for each fold of each seed are then compiled in a seed_summary.xlsx file. 
These are compared to actual_value
RMSE and RSQ are calculated
	RMSE by;
		=SQRT(SUMSQ(A1:A10)/COUNTA(A1:A10))  	###for example, where A1:A10 is an array for the method, difference between predicted and actual value
	RSQ by;
		=RSQ(A1:A10,B1:B10)  					###for example, A1:A10 predicted value, B1:B10 actual value

Overal Summary compiled, see "Summary of Results.xlsx"
Contains the RMSE and R2 values for each fold of each seed. This used to produce
	rmse-summary.csv
	r2-summary.csv
which in turn are used to produce RMSE-Summary.png, R2-Summary.png (etc). 
		rmse<-read.csv("rmse-summary.csv")
		palette(rainbow(15))
		boxplot(rmse, main="TITLE", xlab="x axis label", ylab="etc", ylim=c(0,2.5), col=(paste(palette())))



########################
Notes:
	5 Seed Numbers used as 10 seems excessive with no discernible gain in quality. See a folder in "R Script Dev" for data comparisons. 
	Linear models(glm, lm, and linear kernals for rvm,svm,gauss) were removed due to their consistently low performance. 
	
	Final models used: 
		rvmRadial, rvmPoly, svmRadial, svmPoly, gaussprRadial, gaussPoly, nnet, knn, treebad, rf, cubist, pls, pcr, ens_greedy (greedy ensemble with caretEnsemble), ens_linear (glm ensemble stack) ## 15 models total
	
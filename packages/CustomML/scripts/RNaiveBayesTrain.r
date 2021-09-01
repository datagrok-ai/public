#name: RNaiveBayesTrain
#meta.mlname: RNaiveBayes
#meta.mlrole: train
#description: Custom R Train func for Naive Bayes
#language: r
#input: dataframe df
#input: string predict_column
#input: string laplace {category: Parameters}
#output: blob model

#>>> Initialize local environment and install/require necessary packages
library(caret)
library(e1071)
 
#if (!require(naiveBayes)) install.packages('naiveBayes')  #install package if not already there
#library(naiveBayes)
#install.packages('fastDummies')
#renv::install("fastDummies")
#install.packages("scales") 
#require(fastDummies)

#>>> Extract/Prepare train features and target variable
trainY <- df[[predict_column]]
trainX <- df[, -grep(predict_column, colnames(df))]

#>>> Convert categorical features into dummy variables
dummies <- dummyVars(" ~ .", data = trainX)
trainX = data.frame(predict(dummies, newdata = trainX))

#>>> Train model
trained_model <- naiveBayes(x=trainX, y=trainY, laplace = laplace)

#>>> Save trained model 
saveRDS(trained_model, file=model)
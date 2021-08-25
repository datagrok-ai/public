#name: RKNNApply
#meta.mlname: RKNN
#meta.mlrole: apply
#description: Custom R Apply func for KNN 
#language: r
#input: blob model
#input: dataframe df
#input: string namesKeys [original features' names]
#input: string namesValues [new features' names]
#output: dataframe data_out

#>>> Load necessary packages
library(caret)
library(e1071)
renv::init()
if (!require(kknn)) install.packages("kknn")  #install package if not already there
install.packages('kknn')
library(kknn)
if (!require(MLmetrics)) install.packages("MLmetrics")  #install package if not already there
library(MLmetrics)

#>>> If original(namesKeys) and new(namesValues) passed, map original names to new
namesKeys = strsplit(namesKeys, ',')[[1]]
namesValues = strsplit(namesValues, ',')[[1]]
if (length(namesKeys) > 0) {
  featuresNames = names(df)
  for (n in 0:length(namesKeys))
    names(df)[featuresNames == namesValues[n]] = namesKeys[n]
}

#>>> Retrieve saved/trained model
trained_model = readRDS(model)

#>>> Predict using trained model
predY = predict(trained_model, df)
data_out = data.frame(predY)
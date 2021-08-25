#name: RCaretNBApply
#meta.mlname: RCaretNB
#meta.mlrole: apply
#description: Custom R Apply func for Naive Bayes using Caret
#language: r
#input: blob model
#input: dataframe df
#input: string namesKeys
#input: string namesValues
#output: dataframe data_out

renv::init()
install.packages("naivebayes")
library(naivebayes)

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
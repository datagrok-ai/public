#name: imputeFAMDImpl
#description: Impute the missing values of a mixed dataset using the PCA method "factorial analysis for mixed data"
#language: r
#tags: template, demo
#input: dataframe data {slice:false} [Input data table]
#input: column_list catCols [list of categorical columns]
#input: int ncp = 9 [number of components used to predict the missing entries]
#input: string method = 'Regularized' [reconstruction formulae]
#input: double regCoeff = 1 [regularization coeff]
#output: dataframe imputedDF [imputed dataset]

require(missMDA)

data[ , unlist(catCols)] <- lapply(data[ , unlist(catCols)] , factor)
imputedDF <- imputeFAMD(data, ncp = ncp, method = method, coeff.ridge = regCoeff)
imputedDF <- as.data.frame(sapply(imputedDF$completeObs, as.numeric))
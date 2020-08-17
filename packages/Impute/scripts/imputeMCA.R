#name: imputePCAImpl
#description: Impute the missing values of a mixed dataset using the Multiple Correspondence Analysis method
#language: r
#tags: template, demo
#input: dataframe data [Input data table]
#input: int ncp = 9 [number of components used to predict the missing entries]
#input: string method = 'Regularized' [reconstruction formulae]
#input: double regCoeff = 1 [regularization coeff]
#output: dataframe outputDF [imputed dataset]
require(missMDA)


data <- lapply(data , factor)
imputedDF <- imputePCA(data, ncp = ncp, method = method, coeff.ridge = regCoeff)
imputedDF <- as.data.frame(sapply(imputedDF$completeObs, as.numeric))
#name: imputePCAImpl
#description: Impute the missing values of a continuous dataset using the PCA method
#language: r
#tags: template, demo
#input: dataframe data [Input data table]
#input: column_list catCols [list of categorical columns]
#input: string method = 'Regularized' [reconstruction formulae]
#input: double regCoeff = 1 [regularization coeff]
#output: dataframe imputedDF [imputed dataset]
require(missMDA)

ncomp <- estim_ncpPCA(X, ncp.max = ncol(data) - 2)
imputedDF <- imputePCA(data, ncp = ncomp$ncp, method = method, coeff.ridge = regCoeff)
imputedDF <- as.data.frame(sapply(imputedDF$completeObs, as.numeric))

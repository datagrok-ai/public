#name: imputePCAImpl
#description: Impute the missing values of a continuous dataset using the PCA method
#language: r
#tags: template, demo
#input: dataframe data [Input data table]
#input: column_list columns
#input: column_list catCols {slice:false} [list of categorical columns]
#input: string method = 'Regularized' [reconstruction formulae]
#input: double regCoeff = 1 [regularization coeff]
#output: dataframe imputedDF [imputed dataset]

require(missMDA)
require(gdata)

# convert all variables to numeric
vars_non_num <- names(data)[!sapply(data, is.numeric)]
bigMap <- mapLevels(data[,c(vars_non_num)])
if (length(vars_non_num) != 0) {
  data <- as.data.frame(sapply(data, as.integer)) }

ncomp <- estim_ncpPCA(X, ncp.max = ncol(data) - 2)
imputedDF <- imputePCA(data, ncp = ncomp$ncp, method = method, coeff.ridge = regCoeff)
imputedDF <- as.data.frame(sapply(imputedDF$completeObs, as.numeric))

mapLevels(imputedDF[,c(vars_non_num)]) <- bigMap
imputedDF <- imputedDF[,columns]

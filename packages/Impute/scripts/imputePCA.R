#name: imputePCAImpl
#description: Impute the missing values of a continuous dataset using the PCA method
#language: r
#tags: template, demo
#input: dataframe data {slice:false} [Input data table]
#input: column_list columns
#input: string method {choices: ["regularized", "em"]} [reconstruction formulae]
#input: double regCoeff = 1 [regularization coeff]
#output: dataframe imputedDF [imputed dataset]

require(missMDA)
require(dplyr)

vars_int <- names(data)[sapply(data, is.integer)]
data<-data[rowSums(is.na(data)) != ncol(data), ]

ncomp <- estim_ncpPCA(data, ncp.max = ncol(data) - 2)
imputedDF <- imputePCA(data, ncp = ncomp$ncp, method = method, coeff.ridge = regCoeff)
imputedDF <- as.data.frame(imputedDF$completeObs)

imputedDF <- imputedDF %>% mutate_at(c(vars_int), as.integer)
imputedDF <- imputedDF[,columns]

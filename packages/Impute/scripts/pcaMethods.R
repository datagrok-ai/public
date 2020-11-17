#name: pcaMethodsImpl
#description: Impute the missing values of a mixed dataset using nipals/ppca/bpca/nlpca methods
#language: r
#tags: template, demo
#input: dataframe data [Input data table]
#input: column_list columns
#input: string scaling {choices: ["uv", "vector", 'pareto"]}
#input: string method {choices: ["nipals", "ppca", "bpca", "nlpca"]}
#output: dataframe imputedDF [imputed dataset]

require(pcaMethods)
require(missMDA)
require(dplyr)

data<-data[rowSums(is.na(data)) != ncol(data), ]
vars_int <- names(data)[sapply(data, is.integer)]

method <- gsub(".*:","",method)
res <- prep(data, scale=scaling, center=TRUE, simple=FALSE)
ncomp <- missMDA::estim_ncpPCA(res$data, ncp.max = ncol(data) - 2)
if (ncomp$ncp > 1) {
  imputedDF <- pcaMethods::pca(res$data, method = method, center = FALSE, nPcs = ncomp$ncp)
} else {
  imputedDF <- pcaMethods::pca(res$data, method = method, center = FALSE, nPcs = 2)
}
imputedDF <- prep(imputedDF@completeObs, scale=res$scale, center=res$center, rev=TRUE)

imputedDF <- as.data.frame(imputedDF) %>% mutate_at(c(vars_int), as.integer)



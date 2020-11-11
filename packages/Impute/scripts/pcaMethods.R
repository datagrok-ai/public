#name: pcaMethodsImpl
#description: Impute the missing values of a mixed dataset using nipals/ppca/bpca/nlpca methods
#language: r
#tags: template, demo
#input: dataframe data [Input data table]
#input: column_list columns
#input: string scaling {choices: ["none", "uv", "vector", 'pareto"]}
#input: string method {choices: ["nipals", "ppca", "bpca", "nlpca"]}
#output: dataframe imputedDF [imputed dataset]

require(pcaMethods)
require(missMDA)
require(gdata)

# convert all variables to numeric
vars_non_num <- names(data)[!sapply(data, is.numeric)]
bigMap <- mapLevels(data[,c(vars_non_num)])
if (length(vars_non_num) != 0) {
  data <- as.data.frame(sapply(data, as.integer)) }

method <- gsub(".*:","",method)
if (scaling != 'none') {

  res <- prep(data, scale=scaling, center=TRUE, simple=FALSE)
  ncomp <- missMDA::estim_ncpPCA(res$data, ncp.max = ncol(data) - 2)
  if (ncomp$ncp > 1) {
    imputedDF <- pcaMethods::pca(res$data, method = method, center = FALSE, nPcs = ncomp$ncp)
  } else {
    imputedDF <- pca(res$data, method = method, center = FALSE, nPcs = 2)
  }
  imputedDF <- prep(imputedDF@completeObs, scale=res$scale, center=res$center, rev=TRUE)

} else {

  ncomp <- missMDA::estim_ncpPCA(data, ncp.max = ncol(data) - 2)
  if (ncomp$ncp > 0) {
    imputedDF <- pcaMethods::pca(data, method = method, center = FALSE, nPcs = ncomp$ncp)
  } else {
    imputedDF <- pca(data, method = method, center = FALSE, nPcs = 2)
  }
  imputedDF <- imputedDF@completeObs

}

mapLevels(imputedDF[,c(vars_non_num)]) <- bigMap
imputedDF <- imputedDF[,columns]



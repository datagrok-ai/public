#name: VsnNormalize
#description: Variance-stabilizing normalization via vsn
#language: r
#environment: channels: [conda-forge, bioconda], dependencies: [bioconductor-vsn]
#input: dataframe exprDf
#output: dataframe result

# exprDf columns are named s1, s2, ... containing raw (non-log2) intensities
# VSN performs a generalized log2 transformation that stabilizes variance

library(vsn)

exprMat <- as.matrix(exprDf)

# Run VSN normalization (returns glog2-transformed values)
vsnFit <- vsn::justvsn(exprMat)

result <- as.data.frame(vsnFit)
colnames(result) <- colnames(exprDf)

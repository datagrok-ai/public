#name: limmaDE
#description: Differential expression via limma moderated t-test
#language: r
#environment: channels: [conda-forge, bioconda], dependencies: [bioconductor-limma]
#input: dataframe exprDf
#input: int nGroup1
#input: double fcThreshold = 1.0
#input: double pThreshold = 0.05
#output: dataframe result

# exprDf columns are named s1, s2, ... with group1 first, then group2
nTotal <- ncol(exprDf)
exprMat <- as.matrix(exprDf)
nRows <- nrow(exprMat)
# 1-based input row index, carried through to the result so the client aligns
# stats to proteins by key — never by row position (some R outputs reorder).
rownames(exprMat) <- seq_len(nRows)
g1Idx <- 1:nGroup1
g2Idx <- (nGroup1 + 1):nTotal

hasLimma <- suppressWarnings(require(limma, quietly = TRUE))

if (hasLimma) {
  # Limma empirical Bayes moderated t-test
  group <- factor(c(rep("ctrl", nGroup1), rep("treat", nTotal - nGroup1)),
                  levels = c("ctrl", "treat"))
  design <- model.matrix(~ group)
  fit <- lmFit(exprMat, design)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")

  result <- data.frame(
    row = as.integer(rownames(tt)),
    log2FC = tt$logFC,
    p.value = tt$P.Value,
    adj.p.value = tt$adj.P.Val,
    significant = (abs(tt$logFC) >= fcThreshold) & (tt$adj.P.Val <= pThreshold),
    check.names = FALSE
  )
} else {
  stop("limma package is not available")
}

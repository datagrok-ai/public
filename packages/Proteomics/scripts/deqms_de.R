#name: deqmsDE
#description: DEqMS peptide-count-weighted differential expression
#language: r
#environment: channels: [conda-forge, bioconda], dependencies: [bioconductor-limma, bioconductor-deqms]
#input: dataframe exprDf
#input: int nGroup1
#input: dataframe peptideDf
#input: double fcThreshold = 1.0
#input: double pThreshold = 0.05
#output: dataframe result

# exprDf columns are named s1, s2, ... with group1 first, then group2
nTotal <- ncol(exprDf)
exprMat <- as.matrix(exprDf)
rownames(exprMat) <- seq_len(nrow(exprMat))
counts <- as.numeric(peptideDf[[1]])

# Design matrix (same as limma_de.R)
group <- factor(c(rep("ctrl", nGroup1), rep("treat", nTotal - nGroup1)),
                levels = c("ctrl", "treat"))
design <- model.matrix(~ group)

# Limma pipeline (required by DEqMS as base)
library(limma)
fit <- lmFit(exprMat, design)
fit <- eBayes(fit)

hasDeqms <- suppressWarnings(require(DEqMS, quietly = TRUE))

if (hasDeqms) {
  # DEqMS: peptide-count-weighted variance estimation
  fit$count <- counts
  fit <- spectraCounteBayes(fit)
  # outputResult() returns rows SORTED by significance, NOT input order. The
  # original 1-based input index is preserved in rownames; carry it through so
  # the client aligns stats to proteins by key, never by row position.
  tt <- outputResult(fit, coef = 2)

  result <- data.frame(
    row = as.integer(rownames(tt)),
    log2FC = tt$logFC,
    p.value = tt$sca.P.Value,
    adj.p.value = tt$sca.adj.pval,
    significant = (abs(tt$logFC) >= fcThreshold) & (tt$sca.adj.pval <= pThreshold),
    check.names = FALSE
  )
} else {
  # Fallback: use limma topTable (DEqMS not installed)
  warning("DEqMS not available, using limma instead")
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")

  result <- data.frame(
    row = as.integer(rownames(tt)),
    log2FC = tt$logFC,
    p.value = tt$P.Value,
    adj.p.value = tt$adj.P.Val,
    significant = (abs(tt$logFC) >= fcThreshold) & (tt$adj.P.Val <= pThreshold),
    check.names = FALSE
  )
}

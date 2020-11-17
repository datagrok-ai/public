#name: cleanMetaImpl
#description: rids the data of columns/rows with NA % above the threshold
#language: r
#tags: template, demo
#input: dataframe data [Input data table]
#output: graphics naPlot
#output: graphics matrixPlot
#output: graphics clusterPlot

require(dplyr)
require(mice)
require(stats)
require(ltm)
require(ggplot2)
require(ggdendro)

# METADATA function
get_data <- function(X) {

  # basic metadata
  comp <- sum(stats::complete.cases(X))
  rows <- nrow(X)
  cols <- ncol(X)
  mat <- stats::cor(X, use = "pairwise.complete.obs", method = "pearson")
  missfrac_per_df <- sum(is.na(X))/(nrow(X) * ncol(X))
  missfrac_per_var <- colMeans(is.na(X))
  na_per_df <- sum(is.na(X))
  na_per_var <- sapply(X, function(x) sum(length(which(is.na(x)))))
  mdpat <- mice::md.pattern(X, plot = F)
  data_names <- colnames(X)
  mdpat <- mdpat[, data_names]

  # checking min_PDM thresholds
  mdpat_count <- mdpat[-c(1, nrow(mdpat)), ]
  min_PDM <- c(5,10,20,50,100,200,500,1000)
  min_PDM_obs <- c()
  for (i in 1:length(min_PDM)) {
    index <- as.numeric(rownames(mdpat_count)) > min_PDM[i]
    mdpat_count_simple <- mdpat_count[index, ]
    min_PDM_obs[i] <- round(sum(100 * as.numeric(rownames(mdpat_count_simple)))/sum(as.numeric(rownames(mdpat_count))))
  }

  min_PDM_df <- cbind(min_PDM, min_PDM_obs)
  row.names(min_PDM_df) <- c(1:8)
  colnames(min_PDM_df) <- c("min_PDM_threshold", "perc_obs_retained")

  # NA correlation
  na_cor <- matrix(nrow = cols, ncol = cols)
  colnames(na_cor) <- colnames(X)
  row.names(na_cor) <- paste0(colnames(X), "_is.na")
  for (i in 1:cols) {
    for (j in 1:cols) {
      na_cor[i,j] <- ltm::biserial.cor(X[,j], is.na(X)[,i], use="complete.obs")
    }
  }

  na_cor <- data.table::as.data.table(na_cor, keep.rownames = TRUE)
  melted_cormat <- data.table::melt(na_cor, id = 'rn')
  p_cor <- ggplot(data = melted_cormat, aes(x=factor(rn, levels=unique(rn)), y=variable, fill=value)) +
    geom_tile() + labs(x= "", y= "") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("Variable - Variable NA Correlation Matrix") +
    theme(plot.title = element_text(hjust = 0.5))


  X_update <- as.data.frame(scale(X))
  nm1 <- names(X_update)[colSums(is.na(X_update)) > 0]
  mydescend <- function(x) {
    dplyr::desc(is.na(x))
  }
  arr_X <- X_update %>% dplyr::arrange_at(dplyr::vars(nm1), funs(mydescend))

  # matrix plot
  df_miss_id <- cbind(c(1:rows), arr_X)
  colnames(df_miss_id) <- c("Observations", colnames(arr_X))
  df_miss_id <- data.table::as.data.table(df_miss_id)
  df_melt <- data.table::melt(df_miss_id, id = c('Observations'))
  matrix_plot <- ggplot(df_melt, aes(x = variable, y = Observations)) + geom_tile(aes(fill = value)) +
    scale_fill_gradient(low = "white", high = "blue") + theme(panel.background = element_blank()) +
    ggtitle("Matrix plot of missing data") + theme(plot.title = element_text(hjust = 0.5))


  # cluster plot
  any_miss <- X_update[, which(!colSums(is.na(X_update)) == 0)]

  yesno <- any_miss %>% is.na
  d <- stats::dist(t(yesno), method = "binary")
  hc <- stats::hclust(d, method = "ward.D")
  hcdata <- ggdendro::dendro_data(hc)
  cluster_plot <- ggdendro::ggdendrogram(hcdata, theme_dendro = FALSE) + ggtitle("Cluster plot of missing data") +
    theme(plot.title = element_text(hjust = 0.5)) + labs(x = "variable", y = "Height")

  # output
  list(Complete_cases = comp, Rows = rows, Columns = cols, Corr_matrix = mat, Fraction_missingness = missfrac_per_df,
       Fraction_missingness_per_variable = missfrac_per_var, Total_NA = na_per_df,
       NA_per_variable = na_per_var, MD_Pattern = mdpat, NA_Correlations = na_cor, NA_Correlation_plot = p_cor,
       min_PDM_thresholds = min_PDM_df, Matrix_plot = matrix_plot, Cluster_plot = cluster_plot)
}


# convert all variables to numeric
vars_non_num <- names(data)[!sapply(data, is.numeric)]
if (length(vars_non_num) != 0) {
  data <- as.data.frame(sapply(data, as.numeric)) }

data <- data[rowSums(is.na(data)) != ncol(data), ]

metadata <- get_data(data)
naPlot <- print(metadata$NA_Correlation_plot)
matrixPlot <- print(metadata$Matrix_plot)
clusterPlot <- print(metadata$Cluster_plot)


#naPlot  <- print(ggplot() + theme_void())
#matrixPlot <- print(ggplot() + theme_void())
#clusterPlot <- print(ggplot() + theme_void())
#cleanDf <- X


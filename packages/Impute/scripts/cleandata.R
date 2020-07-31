#name: CleanDataImpl
#description: rids the data of columns/rows with NA % above the threshold
#language: r
#tags: template, demo
#input: dataframe data [Input data table]
#input: column_list columns [list of all columns of interest]
#input: double VarRemovalThreshold = 0.5 [% missingness per column]
#input: double IndRemovalThreshold = 0.5 [% missingness per row]
#output: dataframe X [processed dataframe]

# FUNCTION
clean <- function(X, var_remove = NULL, var_removal_threshold, ind_removal_threshold,
                  missingness_coding = NA) {
  
  # remove undesired variables
  if (!is.null(var_remove)) {X[var_remove] <- NULL }
  
  # give warning when strings are present variables
  strings_present <- sum(sapply(X, is.character)) > 0
  
  if (strings_present == TRUE) {
    stop("Warning! Your data contains string variables. Please inspect your data and either remove these variables using the
var_remove argument or convert them into type factor/numeric where applicable.") }
  
  # convert all variables to numeric
  vars_non_num <- names(X)[!sapply(X, is.numeric)]
  
  if (length(vars_non_num) != 0) {
    X <- as.data.frame(sapply(X, as.numeric)) }
  if (length(vars_non_num) != 0) {
    message(paste("Variable(s) ", (paste(vars_non_num, collapse = ", ")), " converted to numeric.",
                  sep = "")) }
  
  # convert to NA
  X <- as.data.frame(lapply(X, function(x) replace(x, x %in% missingness_coding,
                                                   NA)))
  
  # remove variables above missingness threshold
  missfrac_per_var <- colMeans(is.na(X))
  vars_above_thres <- colnames(X)[missfrac_per_var >= var_removal_threshold]
  if (length(vars_above_thres) != 0)
    new_df <- X[, -which(missfrac_per_var >= var_removal_threshold)] else new_df <- X
  
  if (length(vars_above_thres) != 0) {
    message(paste("Variable(s) ", (paste(vars_above_thres, collapse = ", ")), " removed due to exceeding the pre-defined removal threshold (>",
                  var_removal_threshold * 100, "%) for missingness.", sep = "")) }
  
  # remove individuals above missingness threshold
  missfrac_per_ind <- rowMeans(is.na(new_df))
  inds_above_thres <- rownames(X)[missfrac_per_ind >= ind_removal_threshold]
  if (length(inds_above_thres) != 0) {
    clean_df <- new_df[-which(missfrac_per_ind >= ind_removal_threshold), ] } else { clean_df <- new_df }
  
  if (length(inds_above_thres) != 0) {
    message(paste(length(inds_above_thres), " individual(s) removed due to exceeding the pre-defined removal threshold (>",
                  ind_removal_threshold * 100, "%) for missingness.", sep = "")) }
  
  return(clean_df)
  
}

X <- data[,columns]
X <- clean(X, 
           var_removal_threshold = VarRemovalThreshold, 
           ind_removal_threshold = IndRemovalThreshold)
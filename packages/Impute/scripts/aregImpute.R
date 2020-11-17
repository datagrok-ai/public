#name: aregImputeImpl
#description: uses mice alrorithm to perform multiple imputation
#language: r
#tags: template, demo
#input: dataframe data [Input data table]
#input: column_list columns
#input: int nk = 0
#input: bool tlinear = True
#input: string type {choices : ["pmm", "regression", "normpmm"]}
#input: int pmmtype = 1 {range: 1-4}
#input: string match = "weighted" {choices : ["weighted", "closest", "kclosest"]}
#input: string btMethod = "simple" {choices : ["simple", "approximate bayesian"]}
#input: int burnin = 10 [discard initial iterations]
#output: dataframe imputedDF [imputed dataset]

require(Hmisc)
require(gdata)
require(dplyr)

type <- gsub(".*:","",type)
if (type == 'regression') {
  btMethod <- 'approximate bayesian'
} else if (type == 'normpmm') {
  tlinear <- TRUE
}

data<-data[rowSums(is.na(data)) != ncol(data), ]

# convert all variables to numeric
vars_non_num <- names(data)[!sapply(data, is.numeric)]
vars_int <- names(data)[sapply(data, is.integer)]

if (length(vars_non_num) != 0) {
  bigMap <- mapLevels(data[,c(vars_non_num)])
  data <- data %>% mutate_at(c(vars_non_num), as.integer)
}

Xcolnames <- colnames(data)
Xformula <- stats::as.formula(paste("~", paste(Xcolnames, collapse = "+")))

hmisc_algo <- Hmisc::aregImpute(formula = Xformula,
                                data = data,
                                n.impute = 1,
                                nk = nk,
                                tlinear = tlinear,
                                type = type,
                                pmmtype = pmmtype,
                                match = match,
                                boot.method = btMethod,
                                burnin = burnin)

imputedDF <- as.data.frame(Hmisc::impute.transcan(hmisc_algo,
                                                 imputation = 1,
                                                 data = data,
                                                 list.out = TRUE,
                                                 pr = FALSE,
                                                 check = FALSE))

if (length(vars_non_num) != 0) {
  imputedDF <- imputedDF %>% mutate_at(c(vars_non_num), as.integer)
  mapLevels(imputedDF[,c(vars_non_num)]) <- bigMap
}

imputedDF <- imputedDF %>% mutate_at(c(vars_int), as.integer)
imputedDF <- imputedDF[,columns]
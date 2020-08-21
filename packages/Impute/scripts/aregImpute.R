#name: aregImputeImpl
#description: uses mice alrorithm to perform multiple imputation
#language: r
#tags: template, demo
#input: dataframe data [Input data table]
#input: int nk = 0
#input: bool tlinear = True
#input: string type {choices : ["pmm", "regression", "normpmm"]}
#input: int pmmtype = 1 {range: 1-4}
#input: string match = "weighted" {choices : ["weighted", "closest", "kclosest"]}
#input: string btMethod = "simple" {choices : ["simple", "approximate bayesian"]}
#input: int burnin = 10 [discard initial iterations]
#output: dataframe imputedDF [imputed dataset]

require(Hmisc)

type <- gsub(".*:","",type)
if (type == 'regression') {
  btMethod <- 'approximate bayesian'
} else if (type == 'normpmm') {
  tlinear <- TRUE
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
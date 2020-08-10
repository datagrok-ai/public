#name: aregImputeImpl
#description: uses mice alrorithm to perform multiple imputation
#language: r
#tags: template, demo
#input: dataframe data [Input data table]
#input: int burnin = 5 [discard initial iterations]
#output: dataframe outputDF [imputed dataset]

require(Hmisc)

Xcolnames <- colnames(data)
Xformula <- stats::as.formula(paste("~", paste(Xcolnames, collapse = "+")))

hmisc_algo <- Hmisc::aregImpute(formula = Xformula,
                                data = data,
                                n.impute = 1,
                                burnin = 5,
                                nk = 0,
                                type = "pmm",
                                pmmtype = 2)

outputDF <- as.data.frame(Hmisc::impute.transcan(hmisc_algo,
                                                 imputation = 1,
                                                 data = data,
                                                 list.out = TRUE,
                                                 pr = FALSE,
                                                 check = FALSE))
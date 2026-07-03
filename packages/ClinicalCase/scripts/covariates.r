#name: covariates
#friendlyName: Covariates Analysis
#description: Aalen additive regression showing how covariates affect survival over time
#language: r
#input: string coVariates [Covariate terms for the survival model (e.g. age + sex)]
#input: dataframe covariatesDf [Survival data with 'time' 'status' and covariate columns]
#output: graphics plot
#output: dataframe diagnostics

require("ggfortify")
require("survival")
require("ggplot2")

diagnostics = as.data.frame(installed.packages()[,c(1,3:4)])

modform <- as.formula(paste("Surv(time, status)", coVariates, sep = " ~ "))

aa_fit <- aareg(modform, data = covariatesDf)

autoplot(aa_fit) +
  theme_bw()


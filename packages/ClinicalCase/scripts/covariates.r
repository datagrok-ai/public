#name: covariates
#language: r
#input: string coVariates
#input: dataframe covariatesDf
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


#name: survivalPlot
#language: r
#input: string inputStrata = '1'
#input: string confInt = '0.95'
#input: dataframe survivalDf
#output: graphics plot
#output: dataframe diagnostics

require("ggfortify")
require("survival")
require("ggplot2")

diagnostics = as.data.frame(installed.packages()[,c(1,3:4)])

if (inputStrata == '') {
  inputStrata = '1'
}

modform <- as.formula(paste("Surv(time, status)", inputStrata, sep = " ~ "))

km_trt_fit <- survfit(modform, data=survivalDf, conf.int = as.numeric(confInt))

autoplot(km_trt_fit) +
  theme_bw()


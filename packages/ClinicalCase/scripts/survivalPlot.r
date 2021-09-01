#name: survivalPlot
#language: r
#input: string inputStrata = '1'
#input: string confInt = 0.95
#output: graphics

renv::init()
renv::install("ggfortify")

require("ggfortify")
require("tidyverse")
require("survival")
require("ggplot2")
require("dplyr")
require("ggplot")

data(veteran)

if (inputStrata == '') {
  inputStrata = '1'
}

modform <- as.formula(paste("Surv(time, status)", inputStrata, sep = " ~ "))

km_trt_fit <- survfit(modform, data=veteran, conf.int = as.numeric(confInt))
autoplot(km_trt_fit)  + 
  theme_bw()


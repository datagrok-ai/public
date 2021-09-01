#name: covariates
#language: r
#input: string coVariates
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

modform <- as.formula(paste("Surv(time, status)", coVariates, sep = " ~ "))

aa_fit <- aareg(modform, data = veteran)
autoplot(aa_fit) + 
  theme_bw()

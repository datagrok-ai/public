#name: Vmax
#description: Predict the minimum filter size required for the separation of effluent from bioreactors, given the constraints of batch size and total batch time based on a training dataset of time vs filtrate volume for a given filter type, filter area and pressure
#language: r
#input: dataframe test_data {editor: Compute:manualOutlierSelectionDialog; editor-button: Outliers...}
#input: double test_area = 3.5 {caption: Filter Area; units: cm²} [Filter area]
#input: double vbatch = 25 {caption: VBatch; units: L} [Desired batch size to process]
#input: double tbatch = 0.5 {caption: TBatch; units: hr} [Desired process time]
#input: double sf = 1.5 {caption: Safety Factor} [Safety Factor]
#output: dataframe processedData {viewer: Scatter Plot(x: "time (hr)", y: "t/V (hr/(L/m2))", showRegressionLine: "true"); category: Raw Plot}
#output: double vmax { category: Lin Reg Param }
#output: double q0 { category: Lin Reg Param }
#output: double slope_est { category: Lin Reg Param }
#output: double slope_err { category: Lin Reg Param }
#output: double icept_est { category: Lin Reg Param }
#output: double icept_err { category: Lin Reg Param }
#output: double rsq { category: Lin Reg Param }
#output: double f_stat { category: Lin Reg Param }
#output: double df { category: Lin Reg Param }
#output: double sigma { category: Lin Reg Param }
#output: double resid_ss { category: Lin Reg Param }
#output: double expl_ss { category: Lin Reg Param }
#output: double ss_tot { category: Lin Reg Param }
#output: double amin { category: Recommendation }

a <- as.data.frame(test_data)
a <- a[!a$isOutlier=='true',]
processedData <- a
processedData$time..min. <- processedData$time..min. / 60
names(processedData)[names(processedData) == "time..min."] <- "time (hr)"
processedData$filtrate.volume..mL. <- (processedData$filtrate.volume..mL./1000)/(test_area/10000)
names(processedData)[names(processedData) == "filtrate.volume..mL."] <- "V (L/m2)"
processedData["t/V (hr/(L/m2))"] <- processedData$`time (hr)` / processedData$`V (L/m2)`

lin_mod <- lm(`t/V (hr/(L/m2))` ~ `time (hr)`, data=processedData)
lin_mod_sum <- summary(lin_mod)
coeff <- lin_mod_sum$coefficients
slope_est <- coeff["`time (hr)`", "Estimate"]
slope_err <- coeff["`time (hr)`", "Std. Error"]
icept_est <- coeff["(Intercept)", "Estimate"]
icept_err <- coeff["(Intercept)", "Std. Error"]
rsq <- lin_mod_sum$r.squared
ad_rsq <- lin_mod_sum$adj.r.squared
f_stat <- lin_mod_sum$fstatistic[1]
df <- lin_mod_sum$df[2]
sigma <- lin_mod_sum$sigma
resid_ss <- sum(lin_mod_sum$residuals^2)
expl_ss <- sum((fitted(lin_mod) - mean(processedData$`t/V (hr/(L/m2))`))^2)
vmax <- round(1 / slope_est)
q0 <- round(1 / icept_est)
amin <- (vbatch / vmax) + (vbatch / tbatch / q0)
amin <- round(amin, digits = 1)
ss_tot <- resid_ss + expl_ss
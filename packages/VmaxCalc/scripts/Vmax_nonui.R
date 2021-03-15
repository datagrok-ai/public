#name: VmaxR
#language: r
#input: dataframe test_data
#input: double test_area
#input: double test_pressure
#input: double vbatch
#input: double tbatch
#output: dataframe proc_data
#output: dataframe reg_params

proc_data <- test_data
proc_data$time..min. <- proc_data$time..min./60
names(proc_data)[names(proc_data) == "time..min."] <- "time (hr)"
proc_data$filtrate.volume..mL. <- (proc_data$filtrate.volume..mL./1000)/(test_area/10000)
names(proc_data)[names(proc_data) == "filtrate.volume..mL."] <- "V (L/m2)"
proc_data["t/V (hr/(L/m2))"] <- proc_data$`time (hr)`/proc_data$`V (L/m2)`

lin_mod <- lm(`t/V (hr/(L/m2))` ~ `time (hr)`, data=proc_data)
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
expl_ss <- sum((fitted(lin_mod) - mean(proc_data$`t/V (hr/(L/m2))`))^2)
vmax <- round(1 / slope_est)
q0 <- round(1 / icept_est)
amin <- (vbatch / vmax) + (vbatch / tbatch / q0)

vars <- c('A min', 'V max', 'Q0', 'Slope', 'Slope Err', 'Intercept', 'Intercept Err',
          'R^2', 'Adj. R^2', 'F', 'df', 'sigma', 'SS res', 'SS exp', 'SS tot')
vals <- c(round(amin, digits = 1), vmax, q0, slope_est, slope_err, icept_est, icept_err, rsq, ad_rsq, f_stat,
          df, sigma, resid_ss, expl_ss, resid_ss+expl_ss)
reg_params <- data.frame(vars, vals)






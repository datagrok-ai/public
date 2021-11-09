#name: Maximal Volume
#language: r
#tags: model
#input: dataframe test_data {editor: Compute:manualOutlierSelectionDialog; editor-button: Outliers...}
#input: double test_area = 2.5 {caption: Filter Area; units: cm²} [Filter area]
#input: double test_pressure = 21 {caption: Pressure; units: psi} [Pressure]
#input: double vbatch = 11 {caption: VBatch; units: L} [Desired batch size to process]
#input: double tbatch = 0.8 {caption: TBatch; units: hr} [Desired process time]
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
#output: double amin { category: Min Filter Area }

vmax <- 3
q0 <- 6
slope_est <- vbatch / test_area
slope_err <- 3
icept_est <- slope_err / test_area
icept_err <- 2
rsq <- 2
f_stat <- test_pressure * test_area
df <- 5
sigma <- 4
resid_ss <- 2
expl_ss <- 7
ss_tot <- 5
amin <- vbatch / tbatch
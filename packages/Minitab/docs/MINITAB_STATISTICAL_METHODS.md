# Minitab Statistical Methods — Implementation Reference

> **Purpose:** Algorithm and implementation reference for each statistical analysis the Datagrok Minitab plugin must support. Includes: method description, inputs/outputs, R implementation path, and pharma regulatory context.

---

## 1. Statistical Process Control (SPC)

### 1.1 Control Charts — General Architecture

All control charts share a common structure:
- **Center line (CL):** Process average
- **Control limits (UCL/LCL):** ±3σ from center line (by default; ±2σ and ±1σ optional warning limits)
- **Data points:** Individual measurements or subgroup statistics plotted over time/batch sequence
- **Out-of-control signals:** Points or patterns that violate detection rules

**R package:** `qcc` — the primary package for SPC in R

```r
library(qcc)
# Example: Xbar-R chart
data <- matrix(c(...), ncol = 5)  # 5 observations per subgroup
q <- qcc(data, type = "xbar", nsigmas = 3)
```

### 1.2 Individuals and Moving Range Chart (I-MR)

**When used:** Individual measurements (subgroup size = 1). Most common in pharma: assay results, one result per batch.

**Inputs:**
- Column of numeric measurements
- Optional: subgroup labels (for display)

**Calculations:**
- Individual chart: CL = x̄, UCL = x̄ + 3(MR̄/d2), LCL = x̄ - 3(MR̄/d2)
  - d2 = 1.128 (for n=2 moving range)
  - MR̄ = mean of |xi - xi-1|
- Moving Range chart: CL = MR̄, UCL = D4 × MR̄, LCL = 0 (D3 × MR̄, but D3=0 for n=2)
  - D4 = 3.267, D3 = 0 (for n=2)

**R implementation:**
```r
library(qcc)
i_chart <- qcc(data, type = "xbar.one")
mr_chart <- qcc(data, type = "R")  # Or manually compute moving ranges
```

**Pharma context:** Default chart for release testing results, one value per batch.

### 1.3 Xbar-R Chart (Mean and Range)

**When used:** Subgroups of size 2–8 (small subgroups).

**Inputs:**
- Data matrix: rows = subgroups, columns = observations within subgroup
- Or: column of data + subgroup size

**Calculations:**
- x̄ chart: CL = x̄̄, UCL = x̄̄ + A2 × R̄, LCL = x̄̄ - A2 × R̄
- R chart: CL = R̄, UCL = D4 × R̄, LCL = D3 × R̄
- Control chart constants (A2, D3, D4) by subgroup size n — use standard table

### 1.4 P Chart and Laney P' Chart

**When used:** Proportion nonconforming (pass/fail data); varying sample sizes.

**P chart:** Standard binomial-based; assumes independence.

**Laney P' chart:** Corrected for overdispersion — **required for pharma when sample sizes are large** (analytical testing, where n per subgroup can be thousands of tablets).

```r
# Standard P chart
qcc(defectives, sizes = sample_sizes, type = "p")

# Laney adjustment for overdispersion
# Minitab's Laney P' chart uses sigma_z adjustment
# Implemented manually or via qcc extensions
```

**ICH Q7/GMP context:** P charts are used for visual inspection defect rates, sterility test positives, and out-of-specification rates.

### 1.5 Western Electric (WECO) Rules

Standard set of rules for detecting non-random patterns in control charts. Minitab numbers them 1–8 (same as Western Electric Statistical Quality Control handbook):

| Rule | Description |
|---|---|
| 1 | 1 point more than 3σ from center line |
| 2 | 9 points in a row on same side of center line |
| 3 | 6 points in a row, all increasing or all decreasing |
| 4 | 14 points in a row, alternating up and down |
| 5 | 2 of 3 consecutive points more than 2σ from CL (same side) |
| 6 | 4 of 5 consecutive points more than 1σ from CL (same side) |
| 7 | 15 points in a row within 1σ of CL (either side) |
| 8 | 8 points in a row more than 1σ from CL (either side) |

**Nelson Rules** are nearly identical; Minitab uses them as an alternative. Implement both; let user select.

**R implementation:**
```r
# qcc automatically applies WECO rules (select which ones)
q <- qcc(data, type = "xbar", rules = c(1, 2, 3, 4, 5, 6, 7, 8))
# violations() gives indices of out-of-control points
```

---

## 2. Process Capability Analysis

### 2.1 Normal Capability Indices

**Inputs:** Column of measurements, USL, LSL, optional target T

**Outputs:** Cp, Cpk, Pp, Ppk, Cpm (optional), % out of spec, Sigma level

**Calculations:**
```
# Within-subgroup (short-term) estimates:
Cp  = (USL - LSL) / (6 × σ_within)
CPU = (USL - x̄) / (3 × σ_within)
CPL = (x̄ - LSL) / (3 × σ_within)
Cpk = min(CPU, CPL)

# Overall (long-term) estimates:
Pp  = (USL - LSL) / (6 × σ_overall)
PPU = (USL - x̄) / (3 × σ_overall)
PPL = (x̄ - LSL) / (3 × σ_overall)
Ppk = min(PPU, PPL)

# σ_within (estimated from control chart):
σ_within = MR̄ / d2  (for I-MR chart)
σ_within = R̄ / d2   (for Xbar-R chart)

# σ_overall = sample standard deviation of all data
```

**Confidence intervals** for Cpk (Clopper-Pearson or chi-square based):
```
CI_lower = Cpk × sqrt(1 - z_α/2 × sqrt(1/(9n × Cpk²) + 1/(2(n-1))))
```

**R implementation:**
```r
library(SixSigma)
# Or compute manually; SixSigma has ss.ca.cp() function
library(qcc)
process.capability(qcc_object, spec.limits = c(LSL, USL), target = T)
```

**Pharma standards:**
- Cpk ≥ 1.33 (4σ, minimum for validated processes in most pharma guidelines)
- Cpk ≥ 1.67 (5σ, world-class / six-sigma-aligned)
- ICH Q6A specifies acceptance criteria but not capability indices directly

### 2.2 Non-Normal Capability

**When used:** Data does not follow normal distribution (e.g., dissolution data, impurities near zero).

**Options:**
1. **Transformation:** Box-Cox or Johnson transformation → fit normal → compute Cp/Cpk on transformed data
2. **Distribution fitting:** Fit a non-normal distribution (Weibull, lognormal, gamma, etc.) → compute capability using percentiles: Cp_nonnormal = (USL - LSL) / (P99.865 - P0.135)

**R implementation:**
```r
library(MASS)
# Box-Cox: find optimal lambda
bc <- boxcox(data ~ 1)
lambda <- bc$x[which.max(bc$y)]
data_transformed <- (data^lambda - 1) / lambda

# Johnson transformation: use SuppDists package
library(SuppDists)
fit_johnson <- JohnsonFit(data)

# Non-normal capability
library(capability)
capability(data, spec_lim = c(LSL, USL), distribution = "weibull")
```

### 2.3 Automated Capability Distribution (Minitab's Latest Feature)

Minitab's "Automated Capability" auto-selects the best distribution:
1. Test normality (Anderson-Darling)
2. If normal: use normal capability
3. If not normal: fit Weibull, lognormal, gamma; select best by AD statistic
4. Also test Box-Cox and Johnson transformations
5. Present results for best-fitting distribution

**Implementation:** Implement as a pipeline: normality test → distribution fitting with AIC/BIC comparison → capability calculation.

### 2.4 Capability Sixpack

The most-used Minitab output in pharma QC. Composite panel of six plots:

1. **I-MR chart or Xbar-R chart** (top left) — verify process stability
2. **Run chart** (top center) — check for trends/patterns
3. **Normal probability plot** (top right) — verify normality assumption
4. **Capability histogram** (bottom left) — distribution with spec limits overlay
5. **Last 25 values plot** (bottom center) — recent process behavior
6. **Capability summary** (bottom right) — Cp, Cpk, Pp, Ppk, % out of spec

**Implementation note:** Implement as a Datagrok composite viewer with 6 sub-panels. The layout must be exactly 2×3 (3 top, 3 bottom) — this is what pharma users recognize and expect.

---

## 3. Measurement System Analysis (MSA / Gage R&R)

### 3.1 Crossed Gage R&R (ANOVA Method)

**Inputs:**
- `measurement` column — the response variable
- `part` column — the parts being measured (categorical)
- `operator` column — the operators doing the measuring (categorical)
- `replicate` — implicit from data structure (each operator measures each part multiple times)

**Model:** Two-factor crossed ANOVA with interaction:
```
Y_ijk = μ + Part_i + Operator_j + (Part × Operator)_ij + ε_ijk
```

**ANOVA table:**

| Source | DF | SS | MS | EMS (Expected Mean Squares) |
|---|---|---|---|---|
| Parts | p-1 | SS_P | MS_P | σ²_e + nσ²_PO + noσ²_P |
| Operators | o-1 | SS_O | MS_O | σ²_e + nσ²_PO + npσ²_O |
| Parts×Operators | (p-1)(o-1) | SS_PO | MS_PO | σ²_e + nσ²_PO |
| Repeatability (Error) | po(n-1) | SS_e | MS_e | σ²_e |

**Variance components:**
```
σ²_Repeatability = MS_e
σ²_Reproducibility = (MS_O - MS_PO) / (p × n)
σ²_GRR = σ²_Repeatability + σ²_Reproducibility
σ²_Parts = (MS_P - MS_PO) / (o × n)
σ²_Total = σ²_GRR + σ²_Parts
```

**Key metrics:**
```
%GRR = 100 × σ_GRR / σ_Total           # Should be < 10% for acceptable
%StudyVariation = 100 × (6σ_GRR) / (6σ_Total)
%Tolerance = 100 × (6σ_GRR) / Tolerance_Range   # Where Tolerance = USL - LSL
NDC = 1.41 × (σ_Parts / σ_GRR)         # Number of Distinct Categories; should be ≥ 5
```

**R implementation:**
```r
library(qcc)
grd <- qcc.groups(measurement, interaction(part, operator))

# Or using lme4 for the ANOVA
library(lme4)
model <- lmer(measurement ~ (1|part) + (1|operator) + (1|part:operator), data = df)
vc <- VarCorr(model)
```

**Pharma context:** Required before validating any analytical method (ICH Q2). %GRR < 10% is typically required; 10–30% is marginal and requires justification.

### 3.2 EMP Method (Wheeler's Criteria)

Minitab recently added the EMP (Evaluating the Measurement Process) method as an alternative to the traditional %GRR approach. Wheeler argues that %GRR relative to process variation is more meaningful than relative to tolerance.

**EMP criteria:**
- First-class monitor: Test-retest error < 0.03 × process variation (σ_e/σ_P < 0.03)
- Second-class: 0.03–0.10
- Third-class: 0.10–0.30
- Fourth-class: > 0.30 (inadequate)

---

## 4. Hypothesis Tests

### 4.1 Equivalence Tests (TOST — Two One-Sided Tests)

**Critical for pharma** — used to demonstrate that two analytical methods, two batches, or two processes are equivalent.

**Inputs:**
- Two groups of data (or one group vs. a reference value)
- Equivalence bounds (δ): typically ±20% for bioequivalence; varies for other applications

**TOST approach:**
1. Test H0: |μ1 - μ2| ≥ δ (not equivalent) vs. Ha: |μ1 - μ2| < δ (equivalent)
2. Decomposed into two one-sided t-tests:
   - H01: μ1 - μ2 ≥ δ  (upper bound)
   - H02: μ1 - μ2 ≤ -δ (lower bound)
3. Equivalence declared if BOTH null hypotheses rejected at significance level α

**R implementation:**
```r
library(TOSTER)
TOSTtwo(m1, m2, sd1, sd2, n1, n2, low_eqbound = -delta, high_eqbound = delta, alpha = 0.05)
```

**Regulatory context:** ICH E9 statistical principles; FDA bioequivalence guidance; USP <1010> analytical data interpretation.

### 4.2 Normality Tests

Multiple tests available; Minitab offers all of them:

| Test | R Function | Best For |
|---|---|---|
| Anderson-Darling | `ad.test()` in `nortest` | General; most sensitive to tail behavior |
| Ryan-Joiner | Similar to Shapiro-Wilk | Shapiro-Wilk is the standard equivalent |
| Kolmogorov-Smirnov | `ks.test()` | Two-sample comparison |
| Shapiro-Wilk | `shapiro.test()` | Small samples (n < 50) |

**Pharma guidance:** Always test for normality before capability analysis. Anderson-Darling is the most commonly required test in pharma SOPs.

---

## 5. Design of Experiments (DOE)

### 5.1 Two-Level Factorial Designs

**R packages:** `FrF2` (primary), `DoE.base`

```r
library(FrF2)
# Full 2^3 factorial, 2 replicates, 2 centerpoints
design <- FrF2(8, 3, factor.names = c("Temperature", "Pressure", "Time"),
               replications = 2, ncenter = 2)

# Fractional 2^(6-2) design
design <- FrF2(16, 6, generators = c("ABCE", "BCDF"))
```

### 5.2 Response Surface Designs

**R package:** `rsm`

```r
library(rsm)
# Central Composite Design
ccd_design <- ccd(3, n0 = 3, alpha = "orthogonal",
                  coding = list(x1 ~ (Temp - 50)/10,
                                x2 ~ (Pressure - 100)/20,
                                x3 ~ (Time - 30)/5))
```

### 5.3 Design Analysis — Effects and ANOVA

After running experiments and collecting responses:

```r
# Fit the model
lm_model <- lm(Response ~ (A + B + C)^2 + I(A^2) + I(B^2) + I(C^2), data = design_data)

# ANOVA
anova(lm_model)

# Effects plot (half-normal plot of effects)
library(FrF2)
MEPlot(lm_model)   # Main effects plot
IAPlot(lm_model)   # Interaction plot

# Response optimizer: find factor settings that optimize response
# Use optim() or response surface optimization
```

---

## 6. Stability Studies (ICH Q1E)

### 6.1 Shelf-Life Estimation

**Regulatory reference:** ICH Q1E "Evaluation for Stability Data"

**Model (for assay decline):**
```
Y = β0 + β1 × Time + ε
```

**Algorithm:**
1. Test whether batches can be pooled (ANCOVA):
   - Test for common slope: F-test for batch × time interaction
   - Test for common intercept: F-test for batch effect after accounting for time
   - If pooled: use all batches in one model; otherwise: analyze batches separately and report minimum shelf life

2. Fit linear regression (pooled or per-batch)

3. Calculate 95% lower confidence bound (one-sided):
   ```
   ŷ(t) - t_{α, n-2} × SE(ŷ(t))
   ```

4. Shelf life = time when lower 95% confidence bound intersects acceptance criterion (e.g., 90% of label claim for assay)

**R implementation:**
```r
# Per-batch with test for poolability
library(nlme)

# ANCOVA: test for common slope
model_full <- lm(Assay ~ Time * Batch, data = stability_data)
model_common_slope <- lm(Assay ~ Time + Batch, data = stability_data)
anova(model_common_slope, model_full)  # Test H0: common slope

# If slopes are parallel, test common intercept
model_common_int <- lm(Assay ~ Time, data = stability_data)
anova(model_common_int, model_common_slope)

# Predict with confidence interval
predict_df <- data.frame(Time = seq(0, 36, by = 0.1))
pred <- predict(model_common_int, predict_df, interval = "confidence", level = 0.95)

# Find shelf life: where pred[,"lwr"] < acceptance_limit
shelf_life <- predict_df$Time[which(pred[,"lwr"] < 90)[1]]
```

---

## 7. Regression and ANOVA

### 7.1 Multiple Linear Regression

**R implementation (standard):**
```r
model <- lm(Y ~ X1 + X2 + X3, data = df)
summary(model)
anova(model)
plot(model)  # Residual diagnostic plots
```

**Required outputs matching Minitab:**
- Regression equation
- Coefficient table: Predictor, Coef, SE Coef, T, P-value
- Summary of fit: S (residual SE), R-sq, R-sq(adj), R-sq(pred)
- ANOVA table: Source, DF, SS, MS, F, P
- Residual plots: Normal probability plot of residuals, Residuals vs. Fitted, Histogram of residuals, Residuals vs. Order

### 7.2 General Linear Model (GLM) / ANOVA

```r
# One-way ANOVA
model <- aov(Y ~ Group, data = df)
summary(model)
TukeyHSD(model)  # Post-hoc comparisons

# Two-way ANOVA with interaction
model <- aov(Y ~ A * B, data = df)

# General linear model with random effects
library(lme4)
model <- lmer(Y ~ fixed_factor + (1 | random_factor), data = df)
```

---

## 8. Reliability and Distribution Analysis

### 8.1 Distribution Identification

**When used:** Before capability analysis; to identify best-fitting distribution for non-normal data.

**Distributions to test:** Normal, lognormal, Weibull, exponential, gamma, logistic, log-logistic, smallest extreme value, largest extreme value.

**Selection criteria:** Anderson-Darling statistic (lower = better fit); also AIC/BIC

```r
library(MASS)
library(fitdistrplus)

# Fit multiple distributions and compare
fit_norm <- fitdist(data, "norm")
fit_lnorm <- fitdist(data, "lnorm")
fit_weibull <- fitdist(data, "weibull")
fit_gamma <- fitdist(data, "gamma")

# Compare by AIC
gofstat(list(fit_norm, fit_lnorm, fit_weibull, fit_gamma))
```

### 8.2 Probability Plots

Normal probability plot: plot sorted data vs. expected normal quantiles (Q-Q plot).
```r
qqnorm(data); qqline(data)
# Or use ggplot2 for better formatting
```

Weibull probability plot: plot on Weibull probability paper (log-log scale).

---

## 9. Output Formatting Standards

### Minitab Output Convention

Minitab produces structured output in a consistent format. The plugin should replicate this style for user familiarity:

```
One-Sample T: Dissolution_%

Test of μ = 85 vs ≠ 85

Variable        N    Mean  StDev  SE Mean  95% CI             T      P
Dissolution_%  24  88.342  3.215    0.656  (86.989, 89.695)  5.09  0.000
```

### Capability Analysis Output Convention

```
Process Capability Report for Dissolution_%

LSL     Target     USL
80       90        100

Overall Capability
Pp     0.89
PPL    1.07
PPU    0.70
Ppk    0.70

Potential (Within) Capability
Cp     1.02
CPL    1.24
CPU    0.81
Cpk    0.81

Observed Performance    Expected Overall Perf   Expected Within Perf
% < LSL    0.00         % < LSL    0.51          % < LSL    0.06
% > USL    4.17         % > USL    5.36          % > USL    2.96
% Total    4.17         % Total    5.87          % Total    3.02
```

### Significance Level Reporting

Minitab reports p-values exactly but displays `0.000` for p < 0.0005. Follow this convention for compatibility. Always report exact p-values in the output data table even if display is rounded.

---

## 10. R Package Dependency Summary

| Minitab Feature | R Package | Notes |
|---|---|---|
| SPC control charts | `qcc` | Core SPC package; widely used |
| Process capability | `qcc`, `SixSigma` | `capability` package also useful |
| Distribution fitting | `fitdistrplus`, `MASS` | Comprehensive distribution fitting |
| Gage R&R / MSA | `lme4`, `qcc` | `lme4` for ANOVA; `qcc` for formatted output |
| DOE — factorial | `FrF2`, `DoE.base` | Industry-standard |
| DOE — response surface | `rsm` | Lenth's method for effects |
| Stability (ICH Q1E) | `stability` (custom) | May need custom implementation |
| Regression | base R `lm()` | No package needed |
| Equivalence tests | `TOSTER` | TOST implementation |
| Normality tests | `nortest` | A-D, Lilliefors, etc. |
| Power & sample size | `pwr`, base R | `power.t.test()` etc. |
| Survival / reliability | `survival`, `flexsurv` | |
| Multivariate | base R, `FactoMineR` | PCA, cluster |

---

*See also: `MINITAB_PLUGIN.md` for workflow documentation, `MINITAB_FILE_FORMATS.md` for format specifications.*

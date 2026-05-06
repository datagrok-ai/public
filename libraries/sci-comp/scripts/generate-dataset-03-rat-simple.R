# Generate synthetic rat PK dataset (03_rat_simple.csv) for NCA validation.
#
# Scope (Phase 0): 8 rats, single PO dose, single occasion, 1-compartment model
# with first-order absorption + multiplicative log-normal residual error and
# additive Gaussian noise. Fixed seed for reproducibility.
#
# Full multi-route / multi-occasion / multi-dose rat dataset is deferred to
# Phase 2 (datasets 04 and 06) — see Risk Register R7 in
# docs/nca_development_plan_v2.md.
#
# Output: src/nca/__tests__/datasets/03_rat_simple.csv
# Columns: Subject, Dose, Time, conc
#
# Run from libraries/sci-comp/ working directory:
#   Rscript scripts/generate-dataset-03-rat-simple.R

set.seed(42)

# ------------------------------------------------------------------
# Parameters
# ------------------------------------------------------------------
nRats <- 8L
times <- c(0.25, 0.5, 1, 2, 4, 6, 8, 12)  # hours post-dose
dose  <- 2.5                                # mg per rat (10 mg/kg, ~250 g body weight)

# Population mean PK parameters (1-compartment PO)
mu_ka <- 1.5      # absorption rate constant [1/hr]
mu_ke <- 0.3      # elimination rate constant [1/hr]  -> half-life ~ 2.31 hr
mu_V  <- 0.375    # volume of distribution [L]         (1.5 L/kg * 0.25 kg)
mu_F  <- 0.6      # bioavailability fraction [-]

# Inter-subject variability (log-normal SD on log scale)
sd_ka <- 0.20
sd_ke <- 0.20
sd_V  <- 0.20
sd_F  <- 0.10

# Residual error (conservative for Phase 0 — keeps profiles smooth, no BLQ-like
# zeros and no non-monotonic tails that would complicate the lambda_z best-fit
# validation gate. Edge cases with noisy/BLQ tails are covered by dedicated
# synthetic unit tests in Tasks 1.7 / 1.9, not by this dataset.)
sigma_prop <- 0.10   # proportional (CV)
sigma_add  <- 0.01   # additive [mg/L]

# ------------------------------------------------------------------
# Per-subject parameter sampling
# ------------------------------------------------------------------
ka <- mu_ka * exp(rnorm(nRats, 0, sd_ka))
ke <- mu_ke * exp(rnorm(nRats, 0, sd_ke))
V  <- mu_V  * exp(rnorm(nRats, 0, sd_V))
F_ <- mu_F  * exp(rnorm(nRats, 0, sd_F))

# ------------------------------------------------------------------
# Concentration generation
# ------------------------------------------------------------------
# 1-compartment first-order absorption:
#   C(t) = F*D/V * ka/(ka-ke) * (exp(-ke*t) - exp(-ka*t))
rows <- list()
for (i in seq_len(nRats)) {
  conc_true <- F_[i] * dose / V[i] *
               ka[i] / (ka[i] - ke[i]) *
               (exp(-ke[i] * times) - exp(-ka[i] * times))

  # Mixed proportional + additive residual error
  conc_obs <- conc_true * exp(rnorm(length(times), 0, sigma_prop)) +
              rnorm(length(times), 0, sigma_add)
  conc_obs <- pmax(conc_obs, 0)  # guard against tiny negatives from additive noise

  rows[[i]] <- data.frame(
    Subject = i,
    Dose    = dose,
    Time    = times,
    conc    = round(conc_obs, 4)
  )
}
df <- do.call(rbind, rows)

# ------------------------------------------------------------------
# Guards — fail loudly if generated data violates Phase 0 assumptions
# ------------------------------------------------------------------
if (any(df$conc <= 0)) {
  bad <- df[df$conc <= 0, ]
  stop("Generated data contains non-positive concentrations (",
       nrow(bad), " rows). Reduce noise parameters. Affected:\n",
       paste(capture.output(print(bad)), collapse = "\n"))
}

# Each subject's tail must be monotonically decreasing for the last 3 points
# (so that lambda_z best-fit has a clean target). If not — fail.
for (i in unique(df$Subject)) {
  sub <- df[df$Subject == i, ]
  sub <- sub[order(sub$Time), ]
  tail3 <- tail(sub$conc, 3)
  if (any(diff(tail3) >= 0)) {
    stop("Subject ", i, " has non-monotonic tail: ",
         paste(round(tail3, 4), collapse = " -> "),
         ". Reduce noise parameters or change seed.")
  }
}

# ------------------------------------------------------------------
# Write
# ------------------------------------------------------------------
out <- "src/nca/__tests__/datasets/03_rat_simple.csv"
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
write.csv(df, out, row.names = FALSE)

# ------------------------------------------------------------------
# Sanity output
# ------------------------------------------------------------------
cat("Wrote:", out, "\n")
cat("Rows:", nrow(df), " (expected:", nRats * length(times), ")\n")
cat("Subjects:", paste(unique(df$Subject), collapse = ", "), "\n")
cat("Time range: [", min(df$Time), ",", max(df$Time), "] hr\n")
cat("Conc range: [", min(df$conc), ",", max(df$conc), "] mg/L\n\n")

cat("Per-subject Cmax (observed):\n")
for (i in seq_len(nRats)) {
  sub <- df[df$Subject == i, ]
  cat(sprintf("  Subject %d: Cmax = %.3f at Time = %.2f hr\n",
              i, max(sub$conc), sub$Time[which.max(sub$conc)]))
}

#!/usr/bin/env Rscript
# =============================================================================
# regen-sparse-fixture.R — PKNCA + hand-derived Holder reference for UC-04 sparse
# NCA (sci-comp core/sparse.ts). Rule-18 oracle.
#
# Produces fixtures/05_mouse_sparse.json with three reference blocks:
#
#   1. DESTRUCTIVE (1 sample/animal) — the canonical sparse case. Each dose group
#      of the mouse dataset is run through PKNCA 0.12.1 pk.calc.sparse_auclast,
#      which yields sparse_auclast / sparse_auc_se / sparse_auc_df. An independent
#      hand-coded Holder A1+A3 implementation (below) is cross-checked to match
#      PKNCA to < 1e-6 — proving the formula reduces to Bailer when r_ij = 0.
#      Also emits the composite mean/sd/n/%cv/%blq per nominal time for
#      buildCompositeProfile validation.
#
#   2. BATCH (multiple samples/animal, partial overlap) — PKNCA returns df = NA
#      for batch designs ("Cannot yet calculate sparse degrees of freedom for
#      multiple samples per subject"), so the batch SE + df come from the
#      hand-coded Holder A1 covariance variance + A3 unbiased estimator +
#      Nedelman-Jia 1998 correlated Satterthwaite df. This is the independent
#      oracle for the batch path (synthesis RG-2 / AC-U2).
#
#   3. BATCH-REDUCES-TO-BAILER — the same hand-Holder run on a destructive
#      (r_ij = 0) design must equal the pure-diagonal Bailer variance, asserted
#      against the PKNCA destructive numbers above.
#
# PKNCA's exact algorithm (reverse-engineered from pk.calc.sparse_auc source):
#   weights      = c(0, diff(t)/2) + c(diff(t)/2, 0)        # over ALL design t
#   AUClast      = pk.calc.auc(mean, t, "AUClast", "linear") # to last measurable
#   var, df      = var_sparse_auc(weighted sparse means)     # Holder quadratic
# The t = 0 anchor (conc = 0) is required (PKNCA errors without it for
# extravascular profiles). Trailing all-BLQ timepoints are inert (mean = var = 0).
#
# Usage:  Rscript src/nca/__tests__/regen-sparse-fixture.R
# =============================================================================

suppressPackageStartupMessages({
  library(PKNCA)
  library(jsonlite)
  library(dplyr)
})
stopifnot(as.character(packageVersion("PKNCA")) == "0.12.1")

here <- function(...) file.path("src", "nca", "__tests__", ...)

# ---- Independent Holder A1 + A3 + Nedelman-Jia df implementation -------------
# samples: list over timepoints (ascending, EXCLUDING the t0=0 anchor); each
# element is a named numeric vector (names = animal id, values = concentration,
# BLQ already substituted to 0). times: the nominal timepoints (ascending, no t0).
# nblq: count of BLQ samples per timepoint (for the ">50% BLQ -> mean 0" rule).
# The t0=0 anchor (conc=0, fixed) is prepended internally for the weights.
holder_sparse <- function(times, samples, nblq = NULL) {
  t_all <- c(0, times)
  half  <- diff(t_all) / 2
  w_all <- c(0, half) + c(half, 0)   # PKNCA sparse_auc_weight_linear
  w     <- w_all[-1]                  # weight on each measured timepoint mean
  k     <- length(times)

  means <- vapply(samples, mean, numeric(1))
  r     <- vapply(samples, length, integer(1))
  v     <- vapply(samples, function(x) if (length(x) > 1) var(x) else 0, numeric(1))
  if (is.null(nblq)) nblq <- rep(0L, k)
  # PKNCA "arithmetic mean, <=50% BLQ": zero the mean when >50% of samples BLQ.
  means_auc <- ifelse(nblq / r > 0.5, 0, means)

  # Sigma_hat: covariance matrix of the timepoint means (Holder A1 + A3).
  Sig <- matrix(0, k, k)
  diag(Sig) <- v / r
  rij_mat <- matrix(0L, k, k)
  for (i in seq_len(k - 1)) for (j in (i + 1):k) {
    both <- intersect(names(samples[[i]]), names(samples[[j]]))
    rij  <- length(both)
    rij_mat[i, j] <- rij; rij_mat[j, i] <- rij
    if (rij >= 1) {
      xi <- samples[[i]][both] - means[i]
      xj <- samples[[j]][both] - means[j]
      # A3 unbiased covariance estimator (Holder 2001 eq A3).
      denom <- (rij - 1) + (1 - rij / r[i]) * (1 - rij / r[j])
      sij   <- sum(xi * xj) / denom
      Sig[i, j] <- rij * sij / (r[i] * r[j])   # cov of the means
      Sig[j, i] <- Sig[i, j]
    }
  }
  psi <- as.numeric(t(w) %*% Sig %*% w)         # Var(AUC) = w' Sigma w (A1)

  # Nedelman-Jia 1998 correlated Satterthwaite df (matrix form):
  #   df = (w' Sigma w)^2 / sum_{i,j} (w_i w_j Sigma_ij)^2 / nu_ij
  # nu_ii = r_i - 1 ; nu_ij = r_ij - 1 (off-diagonal, when r_ij >= 2).
  den <- 0
  for (i in seq_len(k)) {
    nu <- max(r[i] - 1, 1)
    den <- den + (w[i]^2 * Sig[i, i])^2 / nu
  }
  for (i in seq_len(k - 1)) for (j in (i + 1):k) if (Sig[i, j] != 0) {
    nu <- max(rij_mat[i, j] - 1, 1)
    den <- den + 2 * (w[i] * w[j] * Sig[i, j])^2 / nu
  }
  df <- if (den > 0) psi^2 / den else NA_real_

  # AUClast on the mean profile (anchored t0=0), standard linear trapezoidal.
  # Uses the >50%-BLQ-zeroed means so trailing BLQ-majority timepoints end AUClast.
  auc <- as.numeric(pk.calc.auc(conc = c(0, means_auc), time = t_all,
                                auc.type = "AUClast", method = "linear"))
  list(auc = auc, se = sqrt(psi), df = df, var = psi,
       weights = w, means = means, n = r, var_t = v, rij = rij_mat)
}

# ---- 1. DESTRUCTIVE: mouse dataset, per dose group --------------------------
raw <- read.csv(here("datasets", "04_mouse_sparse_destructive.csv"),
                stringsAsFactors = FALSE)
df <- raw %>% mutate(Conc_mg_L = ifelse(BLQ == "True", 0, Conc_ng_mL / 1000))

destructive <- list()
max_err <- 0
for (d in sort(unique(df$Dose_mg_kg))) {
  sub <- df %>% filter(Dose_mg_kg == d) %>% arrange(NominalTime_hr)
  # PKNCA oracle (per-animal t0=0 anchor).
  anchor <- sub %>% distinct(MouseID) %>% mutate(NominalTime_hr = 0, Conc_mg_L = 0)
  sub0 <- bind_rows(sub %>% select(MouseID, NominalTime_hr, Conc_mg_L), anchor) %>%
    arrange(MouseID, NominalTime_hr)
  pk <- unlist(pk.calc.sparse_auclast(conc = sub0$Conc_mg_L,
                                      time = sub0$NominalTime_hr,
                                      subject = sub0$MouseID))
  # Hand Holder cross-check.
  times <- sort(unique(sub$NominalTime_hr))
  samples <- lapply(times, function(tt) {
    s <- sub %>% filter(NominalTime_hr == tt); setNames(s$Conc_mg_L, s$MouseID)
  })
  nblq <- vapply(times, function(tt)
    sum(sub$BLQ[sub$NominalTime_hr == tt] == "True"), integer(1))
  h <- holder_sparse(times, samples, nblq)
  err <- max(abs(c(h$auc - pk["sparse_auclast"],
                   h$se  - pk["sparse_auc_se"],
                   h$df  - pk["sparse_auc_df"])))
  max_err <- max(max_err, err)
  cat(sprintf("dose %g: PKNCA auc=%.6f se=%.6f df=%.6f | hand auc=%.6f se=%.6f df=%.6f | err=%.2e\n",
              d, pk["sparse_auclast"], pk["sparse_auc_se"], pk["sparse_auc_df"],
              h$auc, h$se, h$df, err))

  comp <- sub %>% group_by(NominalTime_hr) %>%
    summarize(n = n(), nblq = sum(BLQ == "True"),
              mean = mean(Conc_mg_L), sd = sd(Conc_mg_L),
              cv = ifelse(mean > 0, 100 * sd / mean, NA_real_), .groups = "drop") %>%
    arrange(NominalTime_hr)

  destructive[[as.character(d)]] <- list(
    dose_mg_kg = d,
    sparse_auclast = unname(pk["sparse_auclast"]),
    sparse_auc_se  = unname(pk["sparse_auc_se"]),
    sparse_auc_df  = unname(pk["sparse_auc_df"]),
    composite = list(
      nominal_time = comp$NominalTime_hr, n = comp$n, n_blq = comp$nblq,
      mean = comp$mean, sd = comp$sd, cv_pct = comp$cv)
  )
}
cat(sprintf("max |PKNCA - hand| over destructive doses = %.2e (gate 1e-6)\n", max_err))
stopifnot(max_err < 1e-6)

# ---- 2. BATCH: synthetic partial-overlap design (hand-derived oracle) -------
# 6 animals, each sampled at exactly 2 of 3 timepoints; r_i = 4, r_ij = 2 for all
# pairs. PKNCA cannot produce df here (df = NA), so the hand-Holder values ARE
# the reference. A clean monotone-decay profile.
batch_samples <- list(
  "1" = c(A = 10.0, B = 12.0, C = 11.0, F = 9.5),  # t = 1
  "2" = c(A = 6.0,  B = 7.0,  D = 8.0,  E = 9.0),   # t = 2
  "4" = c(C = 2.0,  D = 3.0,  E = 2.5,  F = 1.5)    # t = 4
)
batch_times <- c(1, 2, 4)
hb <- holder_sparse(batch_times, batch_samples)
cat(sprintf("\nbatch: auc=%.8f se=%.8f df=%.8f (PKNCA df=NA for batch)\n",
            hb$auc, hb$se, hb$df))
# r_ij sanity: every off-diagonal overlap is 2.
stopifnot(hb$rij[1, 2] == 2, hb$rij[1, 3] == 2, hb$rij[2, 3] == 2)

# Verify PKNCA agrees on the batch AUC + SE (only df is unavailable for batch).
batch_long <- do.call(rbind, lapply(seq_along(batch_times), function(i) {
  s <- batch_samples[[i]]
  data.frame(subject = names(s), time = batch_times[i], conc = as.numeric(s))
}))
anchor_b <- data.frame(subject = unique(batch_long$subject), time = 0, conc = 0)
batch0 <- rbind(batch_long, anchor_b)
batch0 <- batch0[order(batch0$subject, batch0$time), ]
pkb <- unlist(pk.calc.sparse_auclast(conc = batch0$conc, time = batch0$time,
                                     subject = batch0$subject))
cat(sprintf("batch PKNCA: auc=%.8f se=%.8f df=%s\n",
            pkb["sparse_auclast"], pkb["sparse_auc_se"],
            as.character(pkb["sparse_auc_df"])))
batch_err <- max(abs(c(hb$auc - pkb["sparse_auclast"], hb$se - pkb["sparse_auc_se"])))
cat(sprintf("batch |PKNCA - hand| (auc,se) = %.2e (gate 1e-6)\n", batch_err))
stopifnot(batch_err < 1e-6)

# ---- 3. Emit fixture --------------------------------------------------------
fixture <- list(
  dataset = "05_mouse_sparse",
  pknca_version = "0.12.1",
  description = paste(
    "Sparse / destructive-sampling NCA oracle (UC-04). Destructive blocks are",
    "PKNCA pk.calc.sparse_auclast; the batch block is hand-derived Holder",
    "A1+A3 + Nedelman-Jia df (PKNCA returns df=NA for batch), cross-checked",
    "against PKNCA on auc+se."),
  config = list(
    auc_method = "linear",   # Bailer/Holder SE is linear-trapezoidal only (AD-3)
    blq_rule = "set-zero",
    sparse_mean_method = "arithmetic mean, <=50% BLQ",
    t0_anchor = "conc=0 at t=0 (required by PKNCA for extravascular profiles)",
    weight = "c(0, diff(t)/2) + c(diff(t)/2, 0) over all design timepoints"),
  destructive = unname(destructive),
  batch = list(
    description = "6 animals x 3 timepoints, each animal sampled twice; r_i=4, r_ij=2 all pairs",
    times = batch_times,
    samples = lapply(seq_along(batch_times), function(i) {
      list(time = batch_times[i],
           animal = names(batch_samples[[i]]),
           conc = as.numeric(batch_samples[[i]]))
    }),
    sparse_auclast = hb$auc,
    sparse_auc_se = hb$se,
    sparse_auc_df = hb$df,
    weights = hb$weights,
    means = unname(hb$means),
    n = unname(hb$n)
  )
)

out <- here("fixtures", "05_mouse_sparse.json")
write_json(fixture, out, auto_unbox = TRUE, digits = 12, pretty = TRUE)
cat(sprintf("\nwrote %s\n", out))
cat("done\n")

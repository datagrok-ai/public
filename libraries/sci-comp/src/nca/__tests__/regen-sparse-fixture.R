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
# w_extra (block 4, IV bolus): added to the weights — the C0 anchor's lever-arm weight t1/2
# times dC0/dmean_i (delta method); NULL everywhere else. It moves se / df only, not auc.
holder_sparse <- function(times, samples, nblq = NULL, w_extra = NULL) {
  t_all <- c(0, times)
  half  <- diff(t_all) / 2
  w_all <- c(0, half) + c(half, 0)   # PKNCA sparse_auc_weight_linear
  w     <- w_all[-1]                  # weight on each measured timepoint mean
  if (!is.null(w_extra)) w <- w + w_extra
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

# ---- 4. IV BOLUS: simulated destructive composites (sparse-iv-bolus-anchor F1) --
# No real destructive IV-bolus dataset exists in the repo (TODO DATA-GAP-SPARSE-IV-ORACLE),
# so this block SIMULATES one with a seeded generator and computes the oracle on the draw.
# Generator (deterministic: set.seed below, R's default Mersenne-Twister + Inversion):
#   C(t) = dose * (A e^(-alpha t) + B e^(-beta t)) * e^(0.1 z)          mg/L, 10 % prop. residual
#   every macro-constant lognormal per animal: p * exp(0.3 z - 0.3^2/2)   (omega 0.3)
#   destructive: one sample per animal, n = 3 per nominal time, 5 significant digits, no BLQ
#   (LLOQ 0.001 ng/mL — BLQ edges are unit-tested in sparse.test.ts / AC-3.1b, not here)
# Composites (keyed by Dose_mg_kg):
#   1 mg/kg  2-cmt A 8, alpha 4, B 2, beta 0.3; first sample 5 min, wide last gap (12 -> 24 h)
#   2 mg/kg  same model; first sample 15 min, 8 -> 24 h last gap (large back-extrapolated share)
#   3 mg/kg  slow 1-cmt (k 0.2/h), 15 / 30 min: redrawn until the composite rises at the
#            first pair (mean(0.5 h) >= mean(0.25 h)) — the chain's 'c1' path
#   4 mg/kg  2-cmt as 1 mg/kg with a dense tail (last gap 2 h) — the terminal note's silent case
# Oracle per composite (PKNCA 0.12.1, independent of sci-comp):
#   c0, c0_method       pk.calc.c0 on the composite mean (the method = the first one of the
#                       default chain that answers alone)
#   sparse_auclast_c0   pk.calc.auciv on the composite mean (linear), cross-checked against
#                       pk.calc.sparse_auclast with every animal's t = 0 anchor set to c0
#                       (equal to < 1e-10 or the script stops)
#   sparse_auclast      the (0, 0)-anchored value (today's sci-comp), with its SE / df
#   sparse_auc_se_const_c0  PKNCA's SE with c0 as a constant (the Phoenix convention)
#   delta_se / delta_df     hand-coded delta-method SE / df (holder_sparse with the C0
#                           sensitivity added to the weights): a TRANSCRIPTION check of
#                           sci-comp, not the independent SE oracle (that is the jest
#                           finite-difference gradient + the parametric calibration).
set.seed(20261110)
iv_models <- list(
  list(dose = 1, A = 8, alpha = 4, B = 2, beta = 0.3,
       times = c(0.083, 0.25, 0.5, 1, 2, 4, 8, 12, 24), path = "logslope"),
  list(dose = 2, A = 8, alpha = 4, B = 2, beta = 0.3,
       times = c(0.25, 0.5, 1, 2, 4, 8, 24), path = "logslope"),
  list(dose = 3, A = 0, alpha = 1, B = 1, beta = 0.2,
       times = c(0.25, 0.5, 1, 2, 4, 8, 12), path = "c1"),
  list(dose = 4, A = 8, alpha = 4, B = 2, beta = 0.3,
       times = c(0.25, 0.5, 1, 2, 3, 4, 6, 8, 10, 12), path = "logslope")
)
IV_N <- 3; IV_OMEGA <- 0.3; IV_SIGMA <- 0.1; IV_LLOQ_NG <- 0.001
iv_draw <- function(m) {
  e <- function() exp(IV_OMEGA * rnorm(1) - IV_OMEGA^2 / 2)
  do.call(rbind, lapply(m$times, function(tt) do.call(rbind, lapply(seq_len(IV_N), function(k) {
    A <- m$A * e(); al <- m$alpha * e(); B <- m$B * e(); be <- m$beta * e()
    c_mg_l <- m$dose * (A * exp(-al * tt) + B * exp(-be * tt)) * exp(IV_SIGMA * rnorm(1))
    data.frame(NominalTime_hr = tt, Conc_ng_mL = signif(1000 * c_mg_l, 5))
  }))))
}
iv_rows <- list(); animal <- 100
for (m in iv_models) {
  repeat {
    d <- iv_draw(m)
    mn <- tapply(d$Conc_ng_mL, d$NominalTime_hr, mean)
    rising <- mn[2] >= mn[1]
    if ((m$path == "c1") == rising) break  # c1 composite must rise, the others must fall
  }
  stopifnot(all(d$Conc_ng_mL > IV_LLOQ_NG))
  n <- nrow(d)
  ids <- sprintf("M%03d", animal + seq_len(n)); animal <- animal + n
  iv_rows[[length(iv_rows) + 1]] <- data.frame(
    StudyID = "MOUSE-IVB-007", MouseID = ids, Species = "Mouse", Strain = "CD-1",
    Sex = rep(c("F", "M"), length.out = n), BodyWeight_g = round(runif(n, 20, 30), 1),
    Route = "IV-bolus", Dose_mg_kg = m$dose, NominalTime_hr = d$NominalTime_hr,
    ActualTime_hr = d$NominalTime_hr, Conc_ng_mL = d$Conc_ng_mL, BLQ = "False",
    LLOQ_ng_mL = IV_LLOQ_NG, Matrix = "Plasma", SamplingType = "Destructive",
    Comment = "Terminal bleed (simulated)", stringsAsFactors = FALSE)
}
iv_csv <- do.call(rbind, iv_rows)
iv_csv_path <- here("datasets", "07_iv_bolus_sparse.csv")
write.csv(iv_csv, iv_csv_path, row.names = FALSE, quote = FALSE)

# The oracle reads the CSV back, so it runs on exactly the committed (rounded) values.
iv <- read.csv(iv_csv_path, stringsAsFactors = FALSE) %>% mutate(Conc_mg_L = Conc_ng_mL / 1000)
c0_method_of <- function(conc, time) {
  for (mth in c("c0", "logslope", "c1", "cmin", "set0")) {
    v <- suppressWarnings(pk.calc.c0(conc, time, method = mth))
    if (!is.na(v)) return(list(value = as.numeric(v), method = mth))
  }
  stop("no c0 method answered")
}
iv_out <- list()
for (m in iv_models) {
  sub <- iv %>% filter(Dose_mg_kg == m$dose) %>% arrange(NominalTime_hr)
  times <- sort(unique(sub$NominalTime_hr))
  samples <- lapply(times, function(tt) {
    s <- sub %>% filter(NominalTime_hr == tt); setNames(s$Conc_mg_L, s$MouseID)
  })
  means <- vapply(samples, mean, numeric(1))
  c0r <- c0_method_of(means, times)
  c0 <- c0r$value
  stopifnot(c0r$method == m$path)
  stopifnot(isTRUE(all.equal(c0, as.numeric(pk.calc.c0(means, times)), tolerance = 0)))

  # (0, 0) baseline and the auciv point estimate on the composite mean.
  auc0 <- as.numeric(pk.calc.auc.last(conc = c(0, means), time = c(0, times), method = "linear"))
  auciv <- as.numeric(pk.calc.auciv(conc = c(0, means), time = c(0, times), c0 = c0, auc = auc0,
                                    method = "linear"))
  # Cross-check: sparse_auclast with every animal's t = 0 anchor at c0.
  anc <- function(v) sub %>% distinct(MouseID) %>% mutate(NominalTime_hr = 0, Conc_mg_L = v)
  run_sparse <- function(v) {
    s0 <- bind_rows(sub %>% select(MouseID, NominalTime_hr, Conc_mg_L), anc(v)) %>%
      arrange(MouseID, NominalTime_hr)
    unlist(pk.calc.sparse_auclast(conc = s0$Conc_mg_L, time = s0$NominalTime_hr,
                                  subject = s0$MouseID))
  }
  pk_c0 <- run_sparse(c0)
  pk_00 <- run_sparse(0)
  stopifnot(abs(pk_c0["sparse_auclast"] - auciv) < 1e-10 * auciv)
  stopifnot(abs(pk_00["sparse_auclast"] - auc0) < 1e-10 * auc0)

  # Hand-coded delta SE / df: the C0 sensitivity distributed onto the weights.
  t1 <- times[1]; t2 <- times[2]
  grad <- rep(0, length(times))
  if (c0r$method == "logslope") {
    grad[1] <-  c0 * t2 / ((t2 - t1) * means[1])
    grad[2] <- -c0 * t1 / ((t2 - t1) * means[2])
  } else if (c0r$method %in% c("c1", "cmin")) {
    grad[which(means == c0)[1]] <- 1
  }
  h0 <- holder_sparse(times, samples)
  stopifnot(abs(h0$auc - auc0) < 1e-10 * auc0, abs(h0$se - pk_00["sparse_auc_se"]) < 1e-9)
  hd <- holder_sparse(times, samples, w_extra = t1 / 2 * grad)
  cat(sprintf(paste0("IV %g mg/kg: c0=%.6f (%s) auc(0,0)=%.6f auc(c0)=%.6f se(0,0)=%.6f ",
                     "se(const c0)=%.6f se(delta)=%.6f df(delta)=%.4f\n"),
              m$dose, c0, c0r$method, auc0, auciv, pk_00["sparse_auc_se"],
              pk_c0["sparse_auc_se"], hd$se, hd$df))
  iv_out[[length(iv_out) + 1]] <- list(
    dose_mg_kg = m$dose, times = times, path = m$path,
    c0 = c0, c0_method = c0r$method,
    sparse_auclast_c0 = auciv,
    sparse_auclast = unname(pk_00["sparse_auclast"]),
    sparse_auc_se = unname(pk_00["sparse_auc_se"]),
    sparse_auc_df = unname(pk_00["sparse_auc_df"]),
    sparse_auc_se_const_c0 = unname(pk_c0["sparse_auc_se"]),
    delta_se = hd$se, delta_df = hd$df,
    composite = list(nominal_time = times, n = unname(h0$n), mean = unname(means)))
}
iv_fixture <- list(
  dataset = "07_iv_bolus_sparse",
  pknca_version = "0.12.1",
  description = paste(
    "Simulated destructive IV-bolus composites (sparse-iv-bolus-anchor F1; generator in",
    "regen-sparse-fixture.R block 4). c0 / c0_method = PKNCA pk.calc.c0 on the composite",
    "mean; sparse_auclast_c0 = PKNCA pk.calc.auciv (linear) = pk.calc.sparse_auclast with",
    "the t=0 anchor at c0; sparse_auclast / _se / _df = the (0, 0)-anchored PKNCA run;",
    "sparse_auc_se_const_c0 = PKNCA's SE with c0 a constant (Phoenix); delta_se / delta_df",
    "= hand-coded delta-method transcription check."),
  generator = list(seed = 20261110, n_per_time = IV_N, omega = IV_OMEGA, sigma = IV_SIGMA,
                   lloq_ng_ml = IV_LLOQ_NG),
  config = list(auc_method = "linear", blq_rule = "set-zero", route = "IV-bolus"),
  composites = iv_out)
for (out_iv in c(here("fixtures", "07_iv_bolus_sparse.json"),
                 file.path("..", "..", "packages", "NCA", "files", "07_iv_bolus_sparse-reference.json"))) {
  write_json(iv_fixture, out_iv, auto_unbox = TRUE, digits = 12, pretty = TRUE)
  cat(sprintf("wrote %s\n", out_iv))
}
stopifnot(file.copy(iv_csv_path, file.path("..", "..", "packages", "NCA", "files", "07_iv_bolus_sparse.csv"),
                    overwrite = TRUE))
cat("done\n")

# =============================================================================
# regen-fixtures.R — PKNCA reference generator for the moment/lag parameters.
#
# Rule-18 oracle. Produces PKNCA 0.12.1 reference values for the FR-200 derived
# parameters (AUMClast, AUMCinf, MRT, Vss, Tlag, %AUMCextrap) on the three
# committed datasets, and authors a new IV-infusion fixture (04) from a known
# one-compartment model so the −T_inf/2 correction can be validated.
#
# Outputs:
#   __tests__/fixtures/_new_params.json   (merged into 01/02/03 by merge-fixtures.mjs)
#   __tests__/datasets/04_iv_infusion.csv (new)
#   __tests__/fixtures/04_iv_infusion.json (new, complete)
#
# It also re-derives the EIGHT already-committed parameters and prints the max
# relative error vs the committed fixtures — a non-zero error means the PKNCA
# configuration here does not reproduce the original run, so the NEW values
# would not be trustworthy. See REGEN.md for the exact invocation.
# =============================================================================

suppressMessages({
  library(PKNCA)
  library(jsonlite)
})

stopifnot(as.character(packageVersion("PKNCA")) == "0.12.1")

here <- function(...) file.path("src", "nca", "__tests__", ...)

# Datasets use either `Time` or `time`; normalise to lowercase.
read_ds <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  names(d)[names(d) == "Time"] <- "time"
  d$Subject <- as.character(d$Subject)
  d
}

# PKNCA options matching the committed fixtures' `config` block.
PKNCA.options(
  auc.method = "lin up/log down",
  min.hl.points = 3,
  min.hl.r.squared = 0.85,
  allow.tmax.in.half.life = FALSE,
  min.span.ratio = 2
)

# Parameters requested for every interval [0, Inf).
interval_template <- function(extra = character(0)) {
  cols <- c(
    "cmax", "tmax", "auclast", "aucinf.obs", "aucpext.obs",
    "half.life", "cl.obs", "vz.obs", "lambda.z",
    "aumclast", "aumcinf.obs", "mrt.obs", "tlag"
  )
  cols <- c(cols, extra)
  iv <- data.frame(start = 0, end = Inf)
  for (c in cols) iv[[c]] <- TRUE
  iv
}

# Extravascular convention (matches computeNca): when a subject has no t=0
# observation, the predose concentration is 0. PKNCA needs that record to
# integrate AUC/AUMC over [0, Inf); IV bolus uses its own c0 back-extrapolation
# and must NOT be zero-prepended.
ensure_t0_zero <- function(conc_df) {
  parts <- lapply(split(conc_df, conc_df$Subject), function(g) {
    if (any(g$time == 0)) return(g)
    rbind(data.frame(Subject = g$Subject[1], time = 0, conc = 0), g)
  })
  do.call(rbind, parts)
}

run_nca <- function(conc_df, dose_df, route, duration = 0, extra = character(0)) {
  if (route == "extravascular") conc_df <- ensure_t0_zero(conc_df)
  o_conc <- PKNCAconc(conc_df, conc ~ time | Subject)
  dose_df$duration <- duration
  o_dose <- PKNCAdose(dose_df, dose ~ time | Subject, route = route,
                      duration = "duration")
  d <- PKNCAdata(o_conc, o_dose, intervals = interval_template(extra))
  res <- suppressWarnings(pk.nca(d))
  as.data.frame(res)
}

# Pull one PPTESTCD value for a subject as a plain numeric (NA -> NA).
getp <- function(df, subj, code) {
  v <- df$PPORRES[df$Subject == subj & df$PPTESTCD == code]
  if (length(v) == 0) return(NA_real_) else return(as.numeric(v[1]))
}

# ---------------------------------------------------------------------------
# Datasets 01–03: compute new params + validate old params.
# ---------------------------------------------------------------------------
new_params <- list()

validate_old <- function(label, df, fixture_path, skip = character(0)) {
  fx <- fromJSON(fixture_path, simplifyVector = FALSE)
  maxerr <- 0; worst <- ""
  for (p in fx$profiles) {
    s <- p$profile_key$subject
    chk <- list(
      cmax = "cmax", tmax = "tmax", auclast = "auclast",
      aucinf = "aucinf.obs", lambda_z = "lambda.z", half_life = "half.life",
      cl = "cl.obs", vz = "vz.obs", pct_aucextrap = "aucpext.obs"
    )
    for (nm in names(chk)) {
      if (nm %in% skip) next
      exp <- p$parameters[[nm]]
      if (is.null(exp) || is.na(exp)) next
      got <- getp(df, s, chk[[nm]])
      if (is.na(got)) { cat(sprintf("  %s subj %s %s: PKNCA NA!\n", label, s, nm)); next }
      rel <- if (nm == "tmax") abs(got - exp) else abs(got - exp) / abs(exp)
      if (rel > maxerr) { maxerr <- rel; worst <- sprintf("%s/subj%s", nm, s) }
    }
  }
  cat(sprintf("[%s] max rel error vs committed: %.2e (worst: %s)\n",
              label, maxerr, worst))
}

collect_new <- function(df, subjects, route) {
  out <- list()
  for (s in subjects) {
    aumclast <- getp(df, s, "aumclast")
    aumcinf <- getp(df, s, "aumcinf.obs")
    mrt <- getp(df, s, if (route == "intravascular") "mrt.iv.obs" else "mrt.obs")
    vss <- if (route == "intravascular") getp(df, s, "vss.iv.obs") else NA_real_
    # Tlag is an absorption concept — N/A for IV routes (matches the core's
    # route gate). PKNCA's tlag on inserted-c0 IV data is a spurious artifact.
    tlag <- if (route == "intravascular") NA_real_ else getp(df, s, "tlag")
    pct_aumc <- if (is.na(aumcinf) || aumcinf == 0) NA_real_ else
      (aumcinf - aumclast) / aumcinf * 100
    out[[as.character(s)]] <- list(
      aumclast = aumclast, aumcinf_obs = aumcinf, mrt = mrt,
      vss = vss, tlag = tlag, pct_aumcextrap = pct_aumc
    )
  }
  out
}

# 01 Theoph — extravascular, dose = Dose(mg/kg) * Wt(kg).
th <- read_ds(here("datasets", "01_theoph.csv"))
th_dose <- aggregate(cbind(Dose, Wt) ~ Subject, data = th, FUN = function(x) x[1])
th_dose$dose <- th_dose$Dose * th_dose$Wt
th_dose$time <- 0
th_df <- run_nca(th[, c("Subject", "time", "conc")],
                 th_dose[, c("Subject", "time", "dose")], "extravascular")
validate_old("01_theoph", th_df, here("fixtures", "01_theoph.json"))
new_params[["01_theoph"]] <- collect_new(th_df, unique(th$Subject), "extravascular")

# 02 Indometh — IV bolus, dose = 25 mg. PKNCA does not back-extrapolate c0 for
# the AUC interval, so we insert the SAME log-linear c0 the core's insertC0
# produces (stored in the committed fixture provenance) at t=0 — making PKNCA's
# AUC/AUMC integrate over the identical augmented profile the core uses.
ind <- read_ds(here("datasets", "02_indometh.csv"))
ind_fx <- fromJSON(here("fixtures", "02_indometh.json"), simplifyVector = FALSE)
ind_c0 <- setNames(
  lapply(ind_fx$profiles, function(p) p$provenance$c0_extrapolated),
  vapply(ind_fx$profiles, function(p) p$profile_key$subject, character(1)))
ind_aug <- do.call(rbind, lapply(split(ind, ind$Subject), function(g) {
  s <- g$Subject[1]
  rbind(data.frame(Subject = s, time = 0, conc = as.numeric(ind_c0[[s]])),
        g[, c("Subject", "time", "conc")])
}))
ind_dose <- data.frame(Subject = unique(ind$Subject), time = 0, dose = 25)
# Already augmented with t=0 → pass as "extravascular" so PKNCA integrates from
# the inserted c0 without its own bolus extrapolation; cl/vz/vss below come
# from the IV-route call. (vss/mrt.iv use the IV-route run.)
ind_df <- run_nca(ind_aug, ind_dose, "intravascular", duration = 0,
                  extra = c("mrt.iv.obs", "vss.obs", "vss.iv.obs"))
# cmax/tmax differ here by construction: the inserted c0 is the array max, but
# the committed fixture reports the OBSERVED Cmax (core excludes inserted c0).
validate_old("02_indometh", ind_df, here("fixtures", "02_indometh.json"),
             skip = c("cmax", "tmax"))
new_params[["02_indometh"]] <- collect_new(ind_df, unique(ind$Subject), "intravascular")

# 03 Rat synthetic — extravascular, dose = Dose (per subject).
rat <- read_ds(here("datasets", "03_rat_simple.csv"))
rat_dose <- aggregate(Dose ~ Subject, data = rat, FUN = function(x) x[1])
rat_dose$dose <- rat_dose$Dose
rat_dose$time <- 0
rat_df <- run_nca(rat[, c("Subject", "time", "conc")],
                  rat_dose[, c("Subject", "time", "dose")], "extravascular")
validate_old("03_rat_simple", rat_df, here("fixtures", "03_rat_simple.json"))
new_params[["03_rat_simple"]] <- collect_new(rat_df, unique(rat$Subject), "extravascular")

write_json(new_params, here("fixtures", "_new_params.json"),
           pretty = TRUE, auto_unbox = TRUE, digits = 12, na = "null")
cat("wrote fixtures/_new_params.json\n")

# ---------------------------------------------------------------------------
# 04 IV-infusion — NEW fixture from a known 1-compartment model.
#   CL = 2 L/h, V = 10 L  ->  k = 0.2 /h ; dose = 100 mg, T_inf = 1 h.
#   Analytic truth: AUCinf = dose/CL = 50 ; MRT_iv = 1/k = 5 ; Vss = V = 10.
# ---------------------------------------------------------------------------
CL <- 2; V <- 10; k <- CL / V; dose <- 100; Tinf <- 1; R0 <- dose / Tinf
cfun <- function(t) ifelse(
  t <= Tinf,
  (R0 / CL) * (1 - exp(-k * t)),
  (R0 / CL) * (1 - exp(-k * Tinf)) * exp(-k * (t - Tinf))
)
inf_times <- c(0, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 12, 16, 24)
inf_conc <- cfun(inf_times)
inf_csv <- data.frame(Subject = "1", time = inf_times, conc = round(inf_conc, 10))
write.csv(inf_csv, here("datasets", "04_iv_infusion.csv"), row.names = FALSE)

inf_conc_df <- data.frame(Subject = "1", time = inf_times, conc = inf_conc)
inf_dose_df <- data.frame(Subject = "1", time = 0, dose = dose)
inf_df <- run_nca(inf_conc_df, inf_dose_df, "intravascular", duration = Tinf,
                  extra = c("mrt.iv.obs", "vss.obs", "vss.iv.obs"))

g <- function(code) getp(inf_df, "1", code)
cat(sprintf("[04 infusion] PKNCA: AUCinf=%.4f MRT.iv=%.4f Vss.iv=%.4f (truth 50/5/10)\n",
            g("aucinf.obs"), g("mrt.iv.obs"), g("vss.iv.obs")))

inf_fixture <- list(
  dataset = "04_iv_infusion",
  pknca_version = "0.12.1",
  config = list(
    auc_method = "lin up/log down", min_points = 3, min_r_squared = 0.85,
    exclude_cmax = TRUE, min_span_ratio = 2, extrap_warn = 20, extrap_error = 50
  ),
  dataset_meta = list(
    dose = list(value = dose, unit = "mg"),
    route = "iv-infusion", infusion_duration = list(value = Tinf, unit = "h"),
    model = "1-compartment: CL=2 L/h, V=10 L, k=0.2 /h (analytic AUCinf=50, MRT_iv=5, Vss=10)"
  ),
  profiles = list(list(
    profile_key = list(subject = "1", route = "iv-infusion"),
    parameters = list(
      cmax = g("cmax"), tmax = g("tmax"), auclast = g("auclast"),
      aucinf = g("aucinf.obs"), pct_aucextrap = g("aucpext.obs"),
      lambda_z = g("lambda.z"), half_life = g("half.life"),
      cl = g("cl.obs"), vz = g("vz.obs"),
      aumclast = g("aumclast"), aumcinf_obs = g("aumcinf.obs"),
      mrt = g("mrt.iv.obs"), vss = g("vss.iv.obs"), tlag = NA_real_,
      pct_aumcextrap = (g("aumcinf.obs") - g("aumclast")) / g("aumcinf.obs") * 100
    )
  ))
)
write_json(inf_fixture, here("fixtures", "04_iv_infusion.json"),
           pretty = TRUE, auto_unbox = TRUE, digits = 12, na = "null")
cat("wrote datasets/04_iv_infusion.csv + fixtures/04_iv_infusion.json\n")

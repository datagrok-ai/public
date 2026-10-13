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
    "aumclast", "aumcinf.obs", "mrt.obs", "tlag",
    # Terminal-phase span ratio — PKNCA's own `span.ratio`, the rule-18 oracle
    # for LambdaZResult.spanRatio. Requested here rather than derived from
    # lambda.z.time.first/last + half.life so the fixture carries the value
    # PKNCA itself reports, not our arithmetic on its inputs.
    "span.ratio"
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

# `options` is passed straight to PKNCAdata (e.g. `list(conc.blq = ...)`); the
# default `list()` keeps PKNCA's global defaults — the configuration every
# committed fixture was produced under. `prepend_t0` = FALSE runs an
# extravascular profile exactly as given (no (0, 0) insertion).
run_nca <- function(conc_df, dose_df, route, duration = 0, extra = character(0),
                    options = list(), prepend_t0 = TRUE) {
  if (route == "extravascular" && prepend_t0) conc_df <- ensure_t0_zero(conc_df)
  o_conc <- PKNCAconc(conc_df, conc ~ time | Subject)
  dose_df$duration <- duration
  o_dose <- PKNCAdose(dose_df, dose ~ time | Subject, route = route,
                      duration = "duration")
  d <- PKNCAdata(o_conc, o_dose, intervals = interval_template(extra),
                 options = options)
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
      vss = vss, tlag = tlag, pct_aumcextrap = pct_aumc,
      # Merged into `provenance` (not `parameters`) — it is a fit diagnostic,
      # sitting with the other lambda_z_* fields. See merge-fixtures.mjs.
      span_ratio = getp(df, s, "span.ratio")
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

# --- 02 Indometh: PKNCA's OWN c0 on the RAW profiles (GROK-20960, F8) -------
# Independent oracle for the core's back-extrapolation: `c0` is a real PKNCA
# PPTESTCD (pk.calc.c0, method chain c0 -> logslope -> c1 -> cmin -> set0), so
# the fixture carries PKNCA's number rather than the core's own c0 fed back in
# (the circularity the augmented run above has by construction). Requested on
# the raw profiles (no inserted row) with route = intravascular.
ind_raw <- ind[, c("Subject", "time", "conc")]
ind_raw_df <- run_nca(ind_raw, ind_dose, "intravascular", duration = 0, extra = "c0")
for (s in unique(ind$Subject)) {
  c0_pknca <- getp(ind_raw_df, s, "c0")
  # Back-extrapolated share of AUCinf (Phoenix AUC_%Back_Ext analogue; no PKNCA
  # equivalent, so the oracle is the stated formula on PKNCA inputs): the
  # dose-time -> first-observation segment, lin-up/log-down (log-down since
  # c0 > C1 — linear fallback otherwise), over PKNCA's aucinf.obs, x 100.
  g <- ind_raw[ind_raw$Subject == s, ]
  g <- g[order(g$time), ]
  t1 <- g$time[1]; c1 <- g$conc[1]
  seg <- if (c0_pknca > c1) (c0_pknca - c1) * t1 / log(c0_pknca / c1) else (c0_pknca + c1) / 2 * t1
  aucinf <- getp(ind_df, s, "aucinf.obs")
  new_params[["02_indometh"]][[s]]$c0_pknca <- c0_pknca
  new_params[["02_indometh"]][[s]]$pct_auc_back_extrap <- seg / aucinf * 100
  cat(sprintf("[02 indometh] subj %s: PKNCA c0 = %.10f (fixture c0_extrapolated %.10f), back-extrap %.2f%% of AUCinf\n",
              s, c0_pknca, as.numeric(ind_c0[[s]]), seg / aucinf * 100))
}

# --- 02 Indometh: STOCK PKNCA AUClast (printed, not asserted) ----------------
# Documents the IV-bolus AUC convention divergence (REGEN.md "IV-bolus AUC
# convention"): sci-comp integrates FROM the back-extrapolated c0 (Phoenix
# WinNonlin), stock PKNCA does NOT — it reports NA without a t=0 datum and
# integrates from the observed (0, 0) when a pre-dose row is present.
ind_zero <- do.call(rbind, lapply(split(ind_raw, ind_raw$Subject), function(g) {
  rbind(data.frame(Subject = g$Subject[1], time = 0, conc = 0), g)
}))
ind_zero_df <- run_nca(ind_zero, ind_dose, "intravascular", duration = 0)
ind_zero_drop_df <- run_nca(ind_zero, ind_dose, "intravascular", duration = 0,
                            options = list(conc.blq = list(first = "drop", middle = "drop", last = "drop")))
for (s in unique(ind$Subject)) {
  cat(sprintf("[02 indometh stock PKNCA] subj %s: auclast raw = %s | with (0,0) row = %s | (0,0) + first=drop = %s | sci-comp/fixture = %.6f\n",
              s, format(getp(ind_raw_df, s, "auclast")), format(getp(ind_zero_df, s, "auclast")),
              format(getp(ind_zero_drop_df, s, "auclast")), getp(ind_df, s, "auclast")))
}

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
    ),
    provenance = list(
      lambda_z_time_first = g("lambda.z.time.first"),
      lambda_z_time_last = g("lambda.z.time.last"),
      span_ratio = g("span.ratio")
    )
  ))
)
write_json(inf_fixture, here("fixtures", "04_iv_infusion.json"),
           pretty = TRUE, auto_unbox = TRUE, digits = 12, na = "null")
cat("wrote datasets/04_iv_infusion.csv + fixtures/04_iv_infusion.json\n")

# ---------------------------------------------------------------------------
# 06 BLQ rules — per-rule PKNCA oracle (GROK-20960 slice U2, Fx-3 / rule 18).
#
# Four hand-authored subjects (datasets/06_blq_rules.csv; BLQ encoded as
# conc = 0 with blq = 1, single LLOQ 0.05) x FIVE rule blocks. Each block maps a
# sci-comp BlqStrategy to the PKNCA `conc.blq` option that means the same
# thing, so every block is a genuine PKNCA run, never a transcribed literal:
#   R-A set-zero x4       -> first/middle/last = "keep"   (the 0 integrates)
#   R-B exclude x4        -> "drop" x3                     (the pre-fix numbers)
#   R-C set-half-lloq x4  -> numeric LLOQ/2 x3             (substitutes enter lambda_z, move tlast)
#   R-D missing x4        -> BLQ rows physically removed   (== PKNCA conc.na = "drop"; asserted == R-B)
#   R-E nca-studio default (set-zero, set-zero, set-zero, exclude) -> "keep" x3
#       PKNCA cannot express the afterLast/consecutiveAfterLast split; with
#       BLQ = 0 the trailing run is trimmed either way, so R-E is expected to
#       equal R-A on every parameter — the GAP-W3 assertion (the shipped default
#       finally has a reference assertion).
# I1 (indometh subject 1 + a FLAGGED (0, 0, blq = 1) pre-dose row): the c0
# oracle is PKNCA's own `c0` on the raw profile (flagged row encoded 0) and the
# AUC oracle is PKNCA on the profile with that row REPLACED by (0, c0_pknca) —
# the same WinNonlin-convention construction the 02 fixture uses. Labelling:
# independent VALUE (PKNCA's c0), same TECHNIQUE (integrate from c0 is
# sci-comp's/WinNonlin's convention imposed on PKNCA's integrator) — see
# REGEN.md "IV-bolus AUC convention". Under every block the flagged row is
# absent for c0 (a substitution rule never makes a pre-dose BLQ the C0), so the
# five I1 blocks are expected identical.
# ---------------------------------------------------------------------------
blq_ds <- read.csv(here("datasets", "06_blq_rules.csv"), stringsAsFactors = FALSE)
blq_ds$Subject <- as.character(blq_ds$Subject)
BLQ_LLOQ <- 0.05
stopifnot(all(blq_ds$lloq == BLQ_LLOQ))

blq_blocks <- list(
  list(id = "R-A",
       sci_comp_blq = list(preFirstMeasurable = "set-zero", embedded = "set-zero",
                           afterLast = "set-zero", consecutiveAfterLast = "set-zero"),
       conc_blq = list(first = "keep", middle = "keep", last = "keep"), remove_rows = FALSE,
       note = "D1/D2 semantics: the substituted 0 integrates; lambda_z excludes conc <= 0 itself"),
  # MEASURED 2026-09-22: PKNCA's first = "drop" leaves the [0, Inf) interval
  # without a t=0 datum and PKNCA does not extrapolate to the origin, so its
  # auclast / aucinf are NA on every PO subject (the same NA as a raw IV-bolus
  # profile). sci-comp's `exclude` at t=0 drops the row and then prepends
  # (0, 0) by the extravascular convention — which is exactly the R-D
  # construction below (rows removed + predose 0). R-D is therefore the AUC
  # oracle for `exclude` (`auc_oracle_block`); R-B itself pins Cmax / Tmax /
  # lambda_z / Tlag parity under PKNCA's own drop option.
  list(id = "R-B",
       sci_comp_blq = list(preFirstMeasurable = "exclude", embedded = "exclude",
                           afterLast = "exclude", consecutiveAfterLast = "exclude"),
       conc_blq = list(first = "drop", middle = "drop", last = "drop"), remove_rows = FALSE,
       auc_oracle_block = "R-D",
       note = "the pre-GROK-20960 numbers (every rule used to integrate as exclude); PKNCA drop leaves AUC NA without a t=0 datum -> AUC oracle = R-D"),
  list(id = "R-C",
       sci_comp_blq = list(preFirstMeasurable = "set-half-lloq", embedded = "set-half-lloq",
                           afterLast = "set-half-lloq", consecutiveAfterLast = "set-half-lloq"),
       conc_blq = list(first = BLQ_LLOQ / 2, middle = BLQ_LLOQ / 2, last = BLQ_LLOQ / 2), remove_rows = FALSE,
       note = "positive substitutes enter lambda_z and move tlast"),
  list(id = "R-D",
       sci_comp_blq = list(preFirstMeasurable = "missing", embedded = "missing",
                           afterLast = "missing", consecutiveAfterLast = "missing"),
       conc_blq = NULL, remove_rows = TRUE,
       note = "BLQ rows physically removed (== PKNCA conc.na = 'drop'); asserted identical to R-B"),
  list(id = "R-E",
       sci_comp_blq = list(preFirstMeasurable = "set-zero", embedded = "set-zero",
                           afterLast = "set-zero", consecutiveAfterLast = "exclude"),
       conc_blq = list(first = "keep", middle = "keep", last = "keep"), remove_rows = FALSE,
       note = "nca-studio shipped default; expected == R-A on every parameter (GAP-W3)")
)

# PKNCA's own c0 for I1 on the RAW profile (flagged pre-dose row encoded 0).
i1_raw <- blq_ds[blq_ds$Subject == "I1", c("Subject", "time", "conc")]
i1_dose <- data.frame(Subject = "I1", time = 0, dose = blq_ds$Dose[blq_ds$Subject == "I1"][1])
i1_raw_df <- run_nca(i1_raw, i1_dose, "intravascular", duration = 0, extra = "c0")
i1_c0 <- getp(i1_raw_df, "I1", "c0")
cat(sprintf("[06 I1] PKNCA c0 on the raw profile (flagged (0,0) row) = %.10f\n", i1_c0))
# The AUC oracle profile: flagged row REPLACED by (0, c0_pknca).
i1_aug <- i1_raw[i1_raw$time > 0, ]
i1_aug <- rbind(data.frame(Subject = "I1", time = 0, conc = i1_c0), i1_aug)

profile_from <- function(df, s, route) {
  iv <- route == "intravascular"
  aumclast <- getp(df, s, "aumclast"); aumcinf <- getp(df, s, "aumcinf.obs")
  list(
    profile_key = list(subject = s, route = if (iv) "iv-bolus" else "po"),
    parameters = list(
      cmax = getp(df, s, "cmax"), tmax = getp(df, s, "tmax"),
      auclast = getp(df, s, "auclast"), aucinf = getp(df, s, "aucinf.obs"),
      pct_aucextrap = getp(df, s, "aucpext.obs"),
      lambda_z = getp(df, s, "lambda.z"), half_life = getp(df, s, "half.life"),
      cl = getp(df, s, "cl.obs"), vz = getp(df, s, "vz.obs"),
      aumclast = aumclast, aumcinf_obs = aumcinf,
      mrt = getp(df, s, if (iv) "mrt.iv.obs" else "mrt.obs"),
      vss = if (iv) getp(df, s, "vss.iv.obs") else NA_real_,
      tlag = if (iv) NA_real_ else getp(df, s, "tlag"),
      pct_aumcextrap = if (is.na(aumcinf) || aumcinf == 0) NA_real_ else (aumcinf - aumclast) / aumcinf * 100
    ),
    provenance = list(
      lambda_z_n_points = getp(df, s, "lambda.z.n.points"),
      lambda_z_time_first = getp(df, s, "lambda.z.time.first"),
      lambda_z_time_last = getp(df, s, "lambda.z.time.last"),
      # PKNCA's fit statistics. pk.calc.half.life (0.12.1 source) selects by
      # lambda.z > 0 + the adj.r.squared tie-break ONLY — min.hl.r.squared is
      # not consulted there — so a PKNCA-reported fit can sit below sci-comp's
      # adj-R^2 floor; the suite reads these to assert that documented case.
      lambda_z_r_squared = getp(df, s, "r.squared"),
      lambda_z_adj_r_squared = getp(df, s, "adj.r.squared"),
      tlast = getp(df, s, "tlast"),
      clast_obs = getp(df, s, "clast.obs"),
      span_ratio = getp(df, s, "span.ratio")
    )
  )
}

po_subjects <- unique(blq_ds$Subject[blq_ds$route == "PO"])
po_dose <- data.frame(Subject = po_subjects, time = 0,
                      dose = vapply(po_subjects, function(s) blq_ds$Dose[blq_ds$Subject == s][1], numeric(1)))
blq_extra <- c("lambda.z.n.points", "lambda.z.time.first", "lambda.z.time.last", "tlast", "clast.obs",
               "r.squared", "adj.r.squared")
blq_out_blocks <- list()
for (b in blq_blocks) {
  po <- blq_ds[blq_ds$route == "PO", c("Subject", "time", "conc", "blq")]
  if (b$remove_rows) po <- po[po$blq == 0, ]
  opts <- if (is.null(b$conc_blq)) list() else list(conc.blq = b$conc_blq)
  po_df <- run_nca(po[, c("Subject", "time", "conc")], po_dose, "extravascular",
                   extra = blq_extra, options = opts)
  iv_df <- run_nca(i1_aug, i1_dose, "intravascular", duration = 0,
                   extra = c("mrt.iv.obs", "vss.obs", "vss.iv.obs", blq_extra), options = opts)
  profiles <- c(lapply(po_subjects, function(s) profile_from(po_df, s, "extravascular")),
                list(profile_from(iv_df, "I1", "intravascular")))
  # I1: cmax/tmax on the augmented profile are the inserted c0 (as in the 02
  # run, where they are skipped). The reported peak is the OBSERVED one, so the
  # oracle is PKNCA's cmax/tmax on the RAW profile (flagged row encoded 0 —
  # not the max).
  profiles[[length(profiles)]]$parameters$cmax <- getp(i1_raw_df, "I1", "cmax")
  profiles[[length(profiles)]]$parameters$tmax <- getp(i1_raw_df, "I1", "tmax")
  profiles[[length(profiles)]]$provenance$c0_pknca <- i1_c0
  for (p in profiles) {
    cat(sprintf("[06 %s] %s: cmax=%s tmax=%s auclast=%s aucinf=%s lambda_z=%s n=%s tlast=%s tlag=%s\n",
                b$id, p$profile_key$subject, format(p$parameters$cmax), format(p$parameters$tmax),
                format(p$parameters$auclast), format(p$parameters$aucinf), format(p$parameters$lambda_z),
                format(p$provenance$lambda_z_n_points), format(p$provenance$tlast), format(p$parameters$tlag)))
  }
  blq_out_blocks[[length(blq_out_blocks) + 1]] <- list(
    id = b$id,
    sci_comp_blq = b$sci_comp_blq,
    config = list(conc_blq = if (is.null(b$conc_blq)) "rows-removed" else b$conc_blq,
                  auc_oracle_block = if (is.null(b$auc_oracle_block)) b$id else b$auc_oracle_block),
    note = b$note,
    profiles = profiles
  )
}

blq_fixture <- list(
  dataset = "06_blq_rules",
  pknca_version = "0.12.1",
  config = list(
    auc_method = "lin up/log down", min_points = 3, min_r_squared = 0.85,
    exclude_cmax = TRUE, min_span_ratio = 2, extrap_warn = 20, extrap_error = 50,
    lloq = BLQ_LLOQ, blq_encoding = "conc = 0 with blq = 1 (PKNCA treats conc == 0 as BLQ)"
  ),
  dataset_meta = list(
    dose = list(po = list(value = po_dose$dose[1], unit = "mg"),
                iv_bolus = list(value = i1_dose$dose[1], unit = "mg")),
    subjects = list(
      P1 = "PO: leading BLQ at t=0 and 0.5, embedded BLQ at t=6, positive tail (preFirstMeasurable + embedded; EV t=0 on a substituted BLQ)",
      P2 = "PO: leading (0, BLQ), clean middle, two trailing BLQ at t=12, 24 (afterLast vs consecutiveAfterLast; trailing trim vs set-half-lloq)",
      P3 = "PO: positive to t=4, then three trailing BLQ, only 2 post-Cmax positives (status 'partial' under set-zero/exclude; under set-half-lloq PKNCA fits the flat LLOQ/2 tail below sci-comp's adj-R^2 floor)",
      P4 = "PO: clean log-linear decay with ONE trailing BLQ that lands on the line under LLOQ/2 (set-half-lloq: the substitute enters an ACCEPTED lambda_z fit and becomes cLast/tLast — the AUCinf-family divergence from PKNCA's clast.obs extrapolation)",
      I1 = "IV bolus: indometh subject 1 + a FLAGGED (0, 0, blq=1) pre-dose row (the dose-time gate under every rule)"
    ),
    iv_bolus_oracle = paste(
      "independent VALUE, same TECHNIQUE: c0_pknca is PKNCA's own c0 on the raw profile;",
      "the AUC oracle integrates PKNCA over the profile with the flagged row replaced by (0, c0_pknca)",
      "- sci-comp's/WinNonlin's C0-in-AUC convention imposed on PKNCA's integrator, so it validates",
      "the integration of the augmented profile, NOT the choice of convention (REGEN.md)."),
    tlag_note = paste(
      "sci-comp Tlag is mask-based (a BLQ sample before the first measurable is the lag boundary,",
      "whatever the rule); PKNCA computes tlag on the CLEANED profile, so the two diverge under",
      "exclude / missing / set-half-lloq for a pre-first BLQ. The fixture records PKNCA's value;",
      "reference-suite.test.ts asserts sci-comp's documented value and PKNCA parity only where they agree.")
  ),
  blocks = blq_out_blocks
)
write_json(blq_fixture, here("fixtures", "06_blq_rules.json"),
           pretty = TRUE, auto_unbox = TRUE, digits = 12, na = "null")
cat("wrote fixtures/06_blq_rules.json\n")

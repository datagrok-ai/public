# =============================================================================
# regen-aucint-fixture.R — PKNCA oracle for nca.aucInterval (interval AUC).
#
# Rule-18 oracle. Produces PKNCA 0.12.1 `aucint.last` (pk.calc.aucint with
# auc.type = "AUClast": interpolation inside the observed range, NO
# extrapolation beyond tlast) for every profile of the committed fixtures
# theoph (01), indometh (02) and rat_simple (03), over the intervals below, under
# BOTH integration methods PKNCA shares with sci-comp ("lin up/log down" and
# "linear"). Each value is a genuine pk.nca run, never a transcribed literal.
#
# Output: __tests__/fixtures/08_aucint.json (read by
# core/__tests__/auc-interval.test.ts). Run from libraries/sci-comp/:
#   Rscript src/nca/__tests__/regen-aucint-fixture.R
#
# The profiles PKNCA integrates are the SAME ones computeNca integrates (see
# REGEN.md "Route handling per dataset"):
#   - theoph, rat (extravascular): a predose (0, 0) is prepended when absent;
#   - indometh (IV bolus): the core's back-extrapolated c0 (committed fixture
#     provenance `c0_extrapolated`) is inserted at t = 0 — sci-comp's
#     WinNonlin-convention IV-bolus AUC, imposed on PKNCA's integrator exactly
#     as in regen-fixtures.R.
#
# Intervals:
#   [0, 12], [0, 24]  — the TK study-day exposure windows (AC-4.3);
#   [2, 10]           — a start bound strictly between samples on a descending
#                       segment (theoph), so start-side log interpolation is
#                       exercised, not only the end side.
# =============================================================================

suppressMessages({
  library(PKNCA)
  library(jsonlite)
})

stopifnot(as.character(packageVersion("PKNCA")) == "0.12.1")

here <- function(...) file.path("src", "nca", "__tests__", ...)

read_ds <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  names(d)[names(d) == "Time"] <- "time"
  d$Subject <- as.character(d$Subject)
  d
}

ensure_t0_zero <- function(conc_df) {
  parts <- lapply(split(conc_df, conc_df$Subject), function(g) {
    if (any(g$time == 0)) return(g)
    rbind(data.frame(Subject = g$Subject[1], time = 0, conc = 0), g)
  })
  do.call(rbind, parts)
}

INTERVALS <- list(c(0, 12), c(0, 24), c(2, 10))
# PKNCA auc.method -> sci-comp AucMethod. PKNCA's third method, "lin-log"
# (linear before tmax, log after), has no sci-comp counterpart, and sci-comp's
# "log-linear" (log on every interval) has no PKNCA counterpart — neither is
# in the oracle.
METHODS <- c("lin up/log down" = "linear-up-log-down", "linear" = "linear")

run_aucint <- function(conc_df, dose_df, route, method) {
  o_conc <- PKNCAconc(conc_df, conc ~ time | Subject)
  dose_df$duration <- 0
  o_dose <- PKNCAdose(dose_df, dose ~ time | Subject, route = route,
                      duration = "duration")
  # PKNCA restricts auclast to its interval, so the self-check below reads
  # auclast from a separate [0, Inf) row.
  iv <- rbind(
    do.call(rbind, lapply(INTERVALS, function(x)
      data.frame(start = x[1], end = x[2], aucint.last = TRUE, auclast = FALSE))),
    data.frame(start = 0, end = Inf, aucint.last = FALSE, auclast = TRUE))
  d <- PKNCAdata(o_conc, o_dose, intervals = iv,
                 options = list(auc.method = method))
  as.data.frame(pk.nca(d))
}

rows <- list()
collect <- function(dataset, df, subjects, method) {
  for (s in subjects) for (x in INTERVALS) {
    sel <- df$Subject == s & df$start == x[1] & df$end == x[2]
    v <- df$PPORRES[sel & df$PPTESTCD == "aucint.last"]
    stopifnot(length(v) == 1, !is.na(v))
    rows[[length(rows) + 1]] <<- list(
      dataset = dataset, subject = s, auc_method = METHODS[[method]],
      start = x[1], end = x[2], aucint_last = as.numeric(v))
  }
}

th <- read_ds(here("datasets", "01_theoph.csv"))
th_dose <- aggregate(cbind(Dose, Wt) ~ Subject, data = th, FUN = function(x) x[1])
th_dose$dose <- th_dose$Dose * th_dose$Wt
th_dose$time <- 0
th_conc <- ensure_t0_zero(th[, c("Subject", "time", "conc")])

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

rat <- read_ds(here("datasets", "03_rat_simple.csv"))
rat_dose <- aggregate(Dose ~ Subject, data = rat, FUN = function(x) x[1])
rat_dose$dose <- rat_dose$Dose
rat_dose$time <- 0
rat_conc <- ensure_t0_zero(rat[, c("Subject", "time", "conc")])

for (m in names(METHODS)) {
  th_df <- run_aucint(th_conc, th_dose[, c("Subject", "time", "dose")], "extravascular", m)
  collect("01_theoph", th_df, unique(th$Subject), m)
  ind_df <- run_aucint(ind_aug, ind_dose, "intravascular", m)
  collect("02_indometh", ind_df, unique(ind$Subject), m)
  rat_df <- run_aucint(rat_conc, rat_dose[, c("Subject", "time", "dose")], "extravascular", m)
  collect("03_rat_simple", rat_df, unique(rat$Subject), m)

  # Self-check: under the fixtures' own method, auclast from this run must
  # reproduce the committed fixtures' auclast — otherwise these profiles are
  # not the ones the reference suite validates.
  if (m == "lin up/log down") {
    for (spec in list(list("01_theoph", th_df), list("02_indometh", ind_df), list("03_rat_simple", rat_df))) {
      fx <- fromJSON(here("fixtures", paste0(spec[[1]], ".json")), simplifyVector = FALSE)
      maxerr <- 0
      for (p in fx$profiles) {
        s <- p$profile_key$subject
        sel <- spec[[2]]$Subject == s & spec[[2]]$PPTESTCD == "auclast" & is.infinite(spec[[2]]$end)
        stopifnot(sum(sel) == 1)
        got <- spec[[2]]$PPORRES[sel]
        maxerr <- max(maxerr, abs(got - p$parameters$auclast) / abs(p$parameters$auclast))
      }
      cat(sprintf("[%s] auclast vs committed fixture: max rel error %.2e\n", spec[[1]], maxerr))
      stopifnot(maxerr < 1e-9)
    }
  }
}

fixture <- list(
  dataset = "08_aucint",
  pknca_version = "0.12.1",
  pknca_parameter = "aucint.last",
  note = paste(
    "PKNCA pk.calc.aucint(auc.type = 'AUClast') per profile of 01/02/03:",
    "a bound strictly between two samples is interpolated (linear on a linear interval,",
    "log on a log interval of the chosen method); past tlast the concentration is 0",
    "(no extrapolation), so an end beyond tlast integrates to tlast only.",
    "02 indometh is the augmented profile (core c0 at t = 0), as in regen-fixtures.R."),
  intervals = lapply(INTERVALS, function(x) list(start = x[1], end = x[2])),
  rows = rows
)
write_json(fixture, here("fixtures", "08_aucint.json"),
           pretty = TRUE, auto_unbox = TRUE, digits = 12)
cat(sprintf("wrote fixtures/08_aucint.json (%d rows)\n", length(rows)))

# Generate PKNCA reference fixtures for NCA validation.
#
# Reads CSV datasets from src/nca/__tests__/datasets/, runs PKNCA with a fixed
# configuration, and writes JSON reference values to src/nca/__tests__/fixtures/.
# These JSON files are the source of truth for the Phase 1 validation gate
# (see Tasks 1.7.5 and 1.9.5 in docs/nca_development_plan_v2.md).
#
# Pinned versions (Phase 0):
#   R       4.6.0
#   PKNCA   0.12.1
#
# Run from libraries/sci-comp/ working directory:
#   Rscript scripts/generate-nca-fixtures.R
#
# Output JSON schema (per dataset file):
#   {
#     "dataset": "<name>",
#     "pknca_version": "<version>",
#     "config": { ... },
#     "dataset_meta": { "dose": {...}, "route": "..." },
#     "profiles": [
#       { "profile_key": { "subject": "...", "route": "..." },
#         "parameters": { "cmax", "tmax", "auclast", "aucinf",
#                         "pct_aucextrap", "lambda_z", "half_life",
#                         "cl", "vz" },
#         "provenance": { ... } },
#       ...
#     ]
#   }
#
# Per Variant-1 convention: `cl` and `vz` are PKNCA-style universal keys.
# Their physical interpretation (CL vs CL/F, Vz vs Vz/F) follows from the
# `route` field in profile_key / dataset_meta.
#
# PKNCA 0.12 specifics:
# - PKNCA requires AUC interval start >= first measurable time. Datasets that
#   don't include t=0 (Indometh, rat_simple) need a t=0 point inserted to
#   capture the [0, t_first] portion of AUC.
# - Insertion strategy by route:
#     * extravascular: insert (0, 0) — by convention pre-dose conc is 0.
#     * IV bolus    : insert (0, c0) where c0 is back-extrapolated from
#                     log-linear regression of the first two observable points.
#                     For IV bolus we then post-process cmax/tmax to be the
#                     OBSERVED peak (not the inserted c0).
# - Result column was renamed PPSTRES -> PPORRES in PKNCA 0.12.
# - Some interval `end = NA` and `end = Inf` semantics changed; we use
#   `end = Inf` together with a valid `start` in [min(time), ...].

suppressPackageStartupMessages({
  library(PKNCA)
  library(jsonlite)
})

# ------------------------------------------------------------------
# Pinned PKNCA configuration (mirrored in NcaRules defaults on TS side)
# ------------------------------------------------------------------
PKNCA_CONFIG <- list(
  auc_method     = "lin up/log down",   # PKNCA 0.12 naming
  min_points     = 3L,
  min_r_squared  = 0.85,
  exclude_cmax   = TRUE,
  min_span_ratio = 2,
  extrap_warn    = 20,
  extrap_error   = 50
)

PKNCA.options(
  auc.method              = PKNCA_CONFIG$auc_method,
  min.span.ratio          = PKNCA_CONFIG$min_span_ratio,
  max.aucinf.pext         = PKNCA_CONFIG$extrap_error,
  conc.na                 = "drop",
  conc.blq                = list(first = "keep", middle = "drop", last = "drop"),
  adj.r.squared.factor    = 1e-4,
  min.hl.points           = PKNCA_CONFIG$min_points,
  min.hl.r.squared        = PKNCA_CONFIG$min_r_squared,
  allow.tmax.in.half.life = !PKNCA_CONFIG$exclude_cmax,
  tau.choices             = NA
)

# Parameters we want extracted for each profile, after t=0 insertion has
# normalised the interval to [0, Inf].
TARGET_INTERVAL <- data.frame(
  start              = 0,
  end                = Inf,
  cmax               = TRUE,
  tmax               = TRUE,
  tlast              = TRUE,
  clast.obs          = TRUE,
  auclast            = TRUE,
  aucinf.obs         = TRUE,
  aucpext.obs        = TRUE,
  lambda.z           = TRUE,
  lambda.z.n.points  = TRUE,
  lambda.z.time.first = TRUE,
  half.life          = TRUE,
  r.squared          = TRUE,
  adj.r.squared      = TRUE,
  cl.obs             = TRUE,
  vz.obs             = TRUE
)

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------

# Estimate c0 (back-extrapolated initial concentration) for an IV bolus profile
# by running PKNCA with NA at t=0 and a c0-only interval.
estimate_c0_per_subject <- function(conc_df_no_zero, dose_df) {
  zero_rows <- data.frame(
    subject = unique(conc_df_no_zero$subject),
    time    = 0,
    conc    = NA_real_
  )
  conc_df_with_na <- rbind(zero_rows, conc_df_no_zero)
  conc_df_with_na <- conc_df_with_na[order(conc_df_with_na$subject,
                                           conc_df_with_na$time), ]
  conc_obj <- PKNCAconc(conc_df_with_na, conc ~ time | subject)
  dose_obj <- PKNCAdose(dose_df, dose ~ time | subject, route = "intravascular")
  intervals <- data.frame(start = 0, end = Inf, c0 = TRUE)
  data_obj <- PKNCAdata(conc_obj, dose_obj, intervals = intervals)
  res <- as.data.frame(pk.nca(data_obj))
  setNames(res$PPORRES[res$PPTESTCD == "c0"],
           res$subject[res$PPTESTCD == "c0"])
}

# Pivot one PKNCA result subset (long format) to a named numeric list.
profile_params_list <- function(prof_df) {
  ok <- is.na(prof_df$exclude) | prof_df$exclude == ""
  prof_df <- prof_df[ok, , drop = FALSE]
  out <- as.list(prof_df$PPORRES)
  names(out) <- prof_df$PPTESTCD
  out
}

# Convert one profile's parameter list into the JSON shape.
# jsonlite serializes NULL inside a list as {}, so we use unbox(NA) for
# missing values — combined with toJSON(na = "null") this yields JSON null.
to_json_profile <- function(subject_id, params, route, observed_cmax, observed_tmax) {
  get_v <- function(k) {
    v <- params[[k]]
    if (is.null(v) || length(v) == 0L || is.na(v)) unbox(NA) else unbox(v)
  }
  pknca_cmax <- get_v("cmax")  # for IV bolus this equals c0 (max of inserted+observed)
  cmax_value <- if (route == "iv-bolus") unbox(observed_cmax) else pknca_cmax
  tmax_value <- if (route == "iv-bolus") unbox(observed_tmax) else get_v("tmax")
  list(
    profile_key = list(
      subject = unbox(as.character(subject_id)),
      route   = unbox(route)
    ),
    parameters = list(
      cmax           = cmax_value,
      tmax           = tmax_value,
      auclast        = get_v("auclast"),
      aucinf         = get_v("aucinf.obs"),
      pct_aucextrap  = get_v("aucpext.obs"),
      lambda_z       = get_v("lambda.z"),
      half_life      = get_v("half.life"),
      cl             = get_v("cl.obs"),
      vz             = get_v("vz.obs")
    ),
    provenance = list(
      lambda_z_n_points       = get_v("lambda.z.n.points"),
      lambda_z_r_squared      = get_v("r.squared"),
      lambda_z_adj_r_squared  = get_v("adj.r.squared"),
      lambda_z_time_first     = get_v("lambda.z.time.first"),
      lambda_z_time_last      = get_v("tlast"),
      clast_obs               = get_v("clast.obs"),
      c0_extrapolated         = if (route == "iv-bolus") pknca_cmax else unbox(NA)
    )
  )
}

# Run full NCA pipeline for one dataset.
# - conc_df / dose_df contain only the original observations (no t=0 insertion).
# - For IV bolus we estimate c0, then insert (0, c0) and re-run.
# - For extravascular we insert (0, 0) directly (or skip if already has t=0).
#
# `route` uses our naming: "iv-bolus", "iv-infusion", "extravascular".
# PKNCA naming: "intravascular" / "extravascular".
run_pknca_dataset <- function(dataset_name, conc_df, dose_df, route, dose_meta) {
  route_pknca <- if (route == "iv-bolus" || route == "iv-infusion") "intravascular"
                 else if (route == "extravascular") "extravascular"
                 else stop("Unknown route: ", route)

  has_t0 <- 0 %in% conc_df$time

  dose_obj_pre <- PKNCAdose(dose_df, dose ~ time | subject, route = route_pknca)

  if (route == "iv-bolus" && !has_t0) {
    c0_map <- estimate_c0_per_subject(conc_df, dose_df)
    zero_rows <- data.frame(
      subject = names(c0_map),
      time    = 0,
      conc    = as.numeric(c0_map)
    )
    conc_df_proc <- rbind(zero_rows, conc_df)
  } else if (!has_t0) {
    # extravascular (or IV infusion) without t=0 → conc(0) = 0 by convention
    zero_rows <- data.frame(
      subject = unique(conc_df$subject),
      time    = 0,
      conc    = 0
    )
    conc_df_proc <- rbind(zero_rows, conc_df)
  } else {
    conc_df_proc <- conc_df
  }
  conc_df_proc <- conc_df_proc[order(conc_df_proc$subject, conc_df_proc$time), ]

  conc_obj <- PKNCAconc(conc_df_proc, conc ~ time | subject)
  dose_obj <- PKNCAdose(dose_df, dose ~ time | subject, route = route_pknca)

  data_obj <- PKNCAdata(conc_obj, dose_obj, intervals = TARGET_INTERVAL)
  res_df   <- as.data.frame(pk.nca(data_obj))

  subjects <- unique(res_df$subject)
  profiles <- lapply(subjects, function(sid) {
    prof <- res_df[res_df$subject == sid, ]
    params <- profile_params_list(prof)
    # Observed Cmax/Tmax come from the ORIGINAL data (before insertion).
    orig <- conc_df[conc_df$subject == sid, ]
    observed_cmax <- max(orig$conc, na.rm = TRUE)
    observed_tmax <- orig$time[which.max(orig$conc)]
    to_json_profile(sid, params, route, observed_cmax, observed_tmax)
  })

  list(
    dataset       = unbox(dataset_name),
    pknca_version = unbox(as.character(packageVersion("PKNCA"))),
    config        = lapply(PKNCA_CONFIG, unbox),
    dataset_meta  = list(
      dose  = list(value = unbox(dose_meta$value), unit = unbox(dose_meta$unit)),
      route = unbox(route)
    ),
    profiles = profiles
  )
}

write_json <- function(payload, out_path) {
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  json_text <- toJSON(payload, pretty = TRUE, digits = 10, na = "null")
  writeLines(json_text, out_path)
  cat("Wrote:", out_path, " (", length(payload$profiles), "profiles)\n")
}

# ------------------------------------------------------------------
# Dataset 01 — Theophylline (oral, has t=0)
# ------------------------------------------------------------------
process_theoph <- function() {
  csv <- read.csv("src/nca/__tests__/datasets/01_theoph.csv",
                  stringsAsFactors = FALSE)
  conc_df <- data.frame(
    subject = as.character(csv$Subject),
    time    = csv$Time,
    conc    = csv$conc
  )
  one_per <- !duplicated(csv$Subject)
  dose_df <- data.frame(
    subject = as.character(csv$Subject[one_per]),
    time    = 0,
    dose    = csv$Dose[one_per] * csv$Wt[one_per]
  )
  payload <- run_pknca_dataset(
    dataset_name = "01_theoph",
    conc_df      = conc_df,
    dose_df      = dose_df,
    route        = "extravascular",
    dose_meta    = list(value = "per-subject (Dose_mg_per_kg * Wt_kg)",
                        unit  = "mg")
  )
  write_json(payload, "src/nca/__tests__/fixtures/01_theoph.json")
}

# ------------------------------------------------------------------
# Dataset 02 — Indomethacin (IV bolus, 25 mg, no t=0)
# ------------------------------------------------------------------
process_indometh <- function() {
  csv <- read.csv("src/nca/__tests__/datasets/02_indometh.csv",
                  stringsAsFactors = FALSE)
  conc_df <- data.frame(
    subject = as.character(csv$Subject),
    time    = csv$time,
    conc    = csv$conc
  )
  subjects <- unique(as.character(csv$Subject))
  dose_df <- data.frame(
    subject = subjects,
    time    = 0,
    dose    = 25
  )
  payload <- run_pknca_dataset(
    dataset_name = "02_indometh",
    conc_df      = conc_df,
    dose_df      = dose_df,
    route        = "iv-bolus",
    dose_meta    = list(value = unbox(25), unit = "mg")
  )
  write_json(payload, "src/nca/__tests__/fixtures/02_indometh.json")
}

# ------------------------------------------------------------------
# Dataset 03 — Synthetic rat (oral, 2.5 mg, no t=0)
# ------------------------------------------------------------------
process_rat_simple <- function() {
  csv <- read.csv("src/nca/__tests__/datasets/03_rat_simple.csv",
                  stringsAsFactors = FALSE)
  conc_df <- data.frame(
    subject = as.character(csv$Subject),
    time    = csv$Time,
    conc    = csv$conc
  )
  one_per <- !duplicated(csv$Subject)
  dose_df <- data.frame(
    subject = as.character(csv$Subject[one_per]),
    time    = 0,
    dose    = csv$Dose[one_per]
  )
  payload <- run_pknca_dataset(
    dataset_name = "03_rat_simple",
    conc_df      = conc_df,
    dose_df      = dose_df,
    route        = "extravascular",
    dose_meta    = list(value = unbox(2.5), unit = "mg")
  )
  write_json(payload, "src/nca/__tests__/fixtures/03_rat_simple.json")
}

# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------
cat("PKNCA version:", as.character(packageVersion("PKNCA")), "\n")
cat("R version:    ", R.version.string, "\n\n")

process_theoph()
process_indometh()
process_rat_simple()

cat("\nDone. Fixtures written to src/nca/__tests__/fixtures/\n")

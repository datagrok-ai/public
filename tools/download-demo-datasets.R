# tools/download-demo-datasets.R
#
# Выгружает три демо-датасета для ANOVA UI:
#   - InsectSprays (R datasets package) → CSV
#   - ToothGrowth  (R datasets package) → CSV
#   - SiRstv       (NIST StRD)          → CSV
#
# Запуск из корня репозитория:
#   Rscript tools/download-demo-datasets.R
#
# Результат сохраняется в packages/EDA/demo-data/.

OUT_DIR <- "packages/EDA/demo-data"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ─── 1. InsectSprays ──────────────────────────────────────────────────
data(InsectSprays)
write.csv(InsectSprays,
          file.path(OUT_DIR, "insect-sprays.csv"),
          row.names = FALSE,
          fileEncoding = "UTF-8")
cat("Wrote insect-sprays.csv (", nrow(InsectSprays), "rows )\n")

# ─── 2. ToothGrowth ───────────────────────────────────────────────────
data(ToothGrowth)
write.csv(ToothGrowth,
          file.path(OUT_DIR, "tooth-growth.csv"),
          row.names = FALSE,
          fileEncoding = "UTF-8")
cat("Wrote tooth-growth.csv (", nrow(ToothGrowth), "rows )\n")

# ─── 3. SiRstv (NIST StRD) ────────────────────────────────────────────
# NIST .dat файл имеет большой текстовый header; данные идут после
# строки "Data:  Instrument   Resistance" (заголовок встроен в маркер).
nist_url <- "https://www.itl.nist.gov/div898/strd/anova/SiRstv.dat"
# На Windows libcurl падает на NIST SSL (CRYPT_E_NO_REVOCATION_CHECK).
# wininet использует системный HTTP-стек и работает; на других ОС — стандартный путь.
read_nist <- function(u) {
  if (.Platform$OS.type == "windows") {
    tmp <- tempfile()
    suppressWarnings(download.file(u, tmp, quiet = TRUE, method = "wininet"))
    readLines(tmp)
  } else {
    readLines(url(u))
  }
}
lines <- read_nist(nist_url)
data_start <- grep("^\\s*Data:\\s+Instrument", lines)
if (length(data_start) != 1)
  stop("Cannot locate 'Data: Instrument' marker in SiRstv.dat")

# Пропустить только строку "Data:..." — следующая уже содержит данные.
data_lines <- lines[(data_start + 1):length(lines)]
data_lines <- data_lines[nzchar(trimws(data_lines))]

# В каждой строке два числа через whitespace: Instrument Resistance
sirstv <- read.table(text = data_lines,
                     col.names = c("Instrument", "Resistance"))
sirstv$Instrument <- as.integer(sirstv$Instrument)

write.csv(sirstv,
          file.path(OUT_DIR, "silicon-resistivity.csv"),
          row.names = FALSE,
          fileEncoding = "UTF-8")
cat("Wrote silicon-resistivity.csv (", nrow(sirstv), "rows )\n")

cat("\nAll three demo datasets written to", OUT_DIR, "\n")

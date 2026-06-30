#!/usr/bin/env bash
# Collapse a long-form Spectronaut Report into the protein-group × sample shape
# the Proteomics package's parseSpectronautText parser consumes.
#
# This is the documented MANUAL FALLBACK for a Spectronaut report too large for
# the in-browser streaming importer: run it locally, then import the small
# aggregated .tsv it produces. It is also the D-04 equivalence oracle the Wave-0
# golden test pins to (see files/demo/README.md).
#
# Run: tools/spectronaut-aggregate.sh <input.tsv> [output.tsv]
#
# Requires the `duckdb` CLI on PATH (v1.3.0 used to derive the committed golden).

set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <input.tsv> [output.tsv]" >&2
  exit 1
fi

IN="$1"
OUT="${2:-${IN%.tsv}.aggregated.tsv}"
SQL="$(dirname "$0")/spectronaut-aggregate.sql"

if [[ ! -f "$IN" ]]; then
  echo "Input not found: $IN" >&2
  exit 1
fi

echo "Aggregating: $IN"
echo "Output:      $OUT"

# Substitute paths into the SQL template. DuckDB's COPY ... TO doesn't accept
# variables, so we materialize the SQL first. SQL single-quote escaping: double up
# any embedded ' in the paths.
IN_ESC="${IN//\'/\'\'}"
OUT_ESC="${OUT//\'/\'\'}"
TMP_SQL="$(mktemp -t spectronaut-aggregate.XXXXXX.sql)"
trap 'rm -f "$TMP_SQL"' EXIT
sed -e "s|__IN_PATH__|${IN_ESC}|g" -e "s|__OUT_PATH__|${OUT_ESC}|g" "$SQL" > "$TMP_SQL"

duckdb < "$TMP_SQL"

echo
echo "Input size:  $(du -h "$IN"  | cut -f1)"
echo "Output size: $(du -h "$OUT" | cut -f1)"

/**
 * Float64-faithful CSV column re-parse (BUG-01 mitigation).
 *
 * `DG.DataFrame.fromCsv` parses numeric text through a single-precision (Float32)
 * intermediate: `1.12` becomes `1.1200000047683716`, breaking PKNCA bit-identity
 * (NFR-04) on direct table reads like Cmax / Tmax. The `doublePrecision` import
 * option only changes the column's *backing array* (Float32Array → Float64Array);
 * the value is already corrupted before storage, so no `fromCsv` option recovers
 * it (verified empirically 2026-06-03). The only fix is to re-parse the original
 * CSV text with JS `Number()` (IEEE-754 double) and rebuild the affected columns.
 *
 * **Scope.** This works only where we hold the source CSV text — `loadFixture`
 * (reference suite) and `loadSampleByName` (bundled samples). The production
 * user-selected-table path arrives already parsed by the platform (text gone),
 * so the Float32 drift remains open there.
 *
 * **Disjoint from #15.** This targets columns DG already typed FLOAT (clean
 * numeric, Float32-parsed). `normalizeNumericNaMarkers` targets STRING columns
 * carrying `.`/empty markers — and already rebuilds those from text via `Number`
 * (full precision). The two operate on non-overlapping column sets and compose.
 *
 * **Parser scope.** A minimal RFC-4180-ish tokenizer (handles quoted fields and
 * `""` escapes) sufficient for the well-formed fixture/sample CSVs. It does NOT
 * support newlines embedded inside quoted cells; if the row count after a plain
 * line split disagrees with the parsed DataFrame, the function aborts (returns
 * `[]`) rather than risk misaligning columns.
 */

import type * as DG from 'datagrok-api/dg';
import * as DGNS from 'datagrok-api/dg';

/** Tokenize one CSV record, honoring double-quoted fields and `""` escapes. */
function splitCsvRow(line: string): string[] {
  const out: string[] = [];
  let cur = '';
  let inQuotes = false;
  for (let i = 0; i < line.length; i++) {
    const ch = line[i];
    if (inQuotes) {
      if (ch === '"') {
        if (line[i + 1] === '"') {
          cur += '"';
          i++;
        } else
          inQuotes = false;
      } else
        cur += ch;
    } else if (ch === '"')
      inQuotes = true;
    else if (ch === ',') {
      out.push(cur);
      cur = '';
    } else
      cur += ch;
  }
  out.push(cur);
  return out;
}

/**
 * Rebuild every FLOAT-typed column of `df` as a Float64 column re-parsed from
 * `csvText`, recovering the precision DG's Float32 CSV parse discards. Empty and
 * `.` cells map to `FLOAT_NULL`. Mutates `df` in place, preserving column order.
 *
 * @returns names of the columns that were re-parsed (for logging / tests).
 */
export function reparseFloatColumnsToFloat64(df: DG.DataFrame, csvText: string): string[] {
  const converted: string[] = [];
  const n = df.rowCount;

  const lines = csvText.split(/\r?\n/);
  if (lines.length < 2)
    return converted;

  const header = splitCsvRow(lines[0]);
  const dataRows = lines.slice(1).map(splitCsvRow);
  // Reconcile the trailing-newline ambiguity: DG may or may not count a final
  // empty line as a null row. Drop trailing fully-empty rows only while we have
  // MORE than the DataFrame's row count, so we land exactly on `n` either way.
  while (dataRows.length > n && dataRows[dataRows.length - 1].every((t) => t.trim() === ''))
    dataRows.pop();
  // Misalignment guard: a quoted newline (unsupported) would desync rows from
  // the parsed DataFrame. Abort rather than corrupt column values.
  if (dataRows.length !== n)
    return converted;

  const floatCols = Array.from(df.columns as Iterable<DG.Column>)
    .filter((c) => c.type === DGNS.COLUMN_TYPE.FLOAT);

  for (const col of floatCols) {
    const idx = header.indexOf(col.name);
    if (idx < 0)
      continue;
    const out = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      const tokens = dataRows[i];
      const cell = idx < tokens.length ? tokens[idx].trim() : '';
      if (cell === '' || cell === '.')
        out[i] = DGNS.FLOAT_NULL;
      else {
        const v = Number(cell);
        out[i] = Number.isFinite(v) ? v : DGNS.FLOAT_NULL;
      }
    }
    df.columns.replace(col, DGNS.Column.fromFloat64Array(col.name, out));
    converted.push(col.name);
  }
  return converted;
}

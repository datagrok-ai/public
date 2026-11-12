/**
 * Numeric NA-marker normalization for ingested DataFrames.
 *
 * NONMEM-canonical CSVs use `.` (literal period) as the missing-value marker:
 * dose rows have no DV, observation rows have no AMT/RATE, so those cells carry
 * `.`. Datagrok's CSV type inference sees the `.` mixed with numbers and resolves
 * the whole column to STRING — which then fails the `isNumeric` filter on the
 * mapping pane's Concentration / Dose / LLOQ / Time pickers, silently hiding the
 * column from the user (Ed observation #15, NCA_OBSERVATIONS.md §#15).
 *
 * The fix runs *after* parse, on whatever DataFrame reaches the analysis — both
 * the bundled-sample path (`loadSampleByName`) and, critically, the
 * user-selected-table path (`ui.input.table`), where the platform already parsed
 * the file before our code sees it. A parse-time `nullStrings` option on
 * `DataFrame.fromCsv` would miss that second path entirely.
 *
 * **Marker set is intentionally narrow.** Only `.` and the empty string (after
 * trim) count as NA, matching Ed's spec ("Treat `.` and empty string as NA in
 * otherwise-numeric columns. Standard pharmacometric convention."). We do NOT
 * treat `NA` / `NaN` / `NULL` as markers — those can be legitimate values of a
 * genuine string category column, and widening the set would silently retype it.
 *
 * **A column converts only when it is unambiguously numeric-with-gaps:** every
 * non-marker cell parses as a finite number AND at least one such value exists.
 * - A subject column like `R001` never converts (non-numeric content).
 * - An all-`.` column (e.g. file 03's empty LLOQ) stays STRING — converting it
 *   to an all-null FLOAT would surface a phantom selectable-but-empty column in
 *   the numeric pickers.
 *
 * Caveat: the heuristic will convert a string column whose only
 * non-numeric cells are `.`/empty even if it was conceptually an identifier
 * (e.g. integer IDs with a stray `.`). NONMEM identifier columns never carry
 * `.`, and pure-integer IDs are parsed as INT (not STRING) so they are never
 * scanned here — so the risk is residual, not zero. Load-bearing: the gain is
 * making the dominant pharmacometric format ingestible without pre-processing.
 */

import type * as DG from 'datagrok-api/dg';
import * as DGNS from 'datagrok-api/dg';

/** NA markers recognized in otherwise-numeric string columns (post-trim). */
const NA_MARKERS: ReadonlySet<string> = new Set(['', '.']);

/**
 * Classify a raw string cell: `null` for an NA marker, the parsed finite number
 * for a numeric value, or `NaN` to signal "not coercible — abort this column".
 */
function classifyCell(raw: string | null): number | null {
  if (raw === null)
    return null;
  const s = raw.trim();
  if (NA_MARKERS.has(s))
    return null;
  // Number('') is 0 and Number(' ') is 0, but both are filtered above as
  // markers; Number('.') is NaN. Reject any non-finite parse (incl. words).
  const n = Number(s);
  return Number.isFinite(n) ? n : NaN;
}

/**
 * In place, rebuild every STRING column that represents numeric data with `.`/
 * empty gaps as a FLOAT column (markers → `FLOAT_NULL`), preserving column order.
 *
 * Mutates `df` directly — the analysis operates on the user's selected table, so
 * retyping there is the intended outcome (decision: normalize in place).
 *
 * @returns the names of the columns that were converted (for logging / tests).
 */
export function normalizeNumericNaMarkers(df: DG.DataFrame): string[] {
  const converted: string[] = [];
  const n = df.rowCount;
  // Snapshot the column list — we mutate it inside the loop.
  const stringCols = Array.from(df.columns as Iterable<DG.Column>)
    .filter((c) => c.type === DGNS.COLUMN_TYPE.STRING);

  for (const col of stringCols) {
    const out = new Float64Array(n);
    let coercible = true;
    let sawValue = false;
    for (let i = 0; i < n; i++) {
      const v = classifyCell(col.get(i) as string | null);
      if (v === null)
        out[i] = DGNS.FLOAT_NULL;
      else if (Number.isFinite(v)) {
        out[i] = v;
        sawValue = true;
      } else {
        coercible = false;
        break;
      }
    }
    if (!coercible || !sawValue)
      continue;
    df.columns.replace(col, DGNS.Column.fromFloat64Array(col.name, out));
    converted.push(col.name);
  }
  return converted;
}

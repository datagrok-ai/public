import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {SEMTYPE} from '../utils/proteomics-types';
import {
  log2TransformColumns, copyAsLog2Columns, detectLog2Status,
  addPrimaryColumnIfNeeded,
} from './shared-utils';
import {setGroups} from '../analysis/experiment-setup';
import {resolveGeneLabels} from '../utils/gene-label-resolver';

/** Required Spectronaut columns for detection. */
const REQUIRED_COLUMNS = ['R.Condition', 'R.Replicate', 'PG.ProteinGroups'];

/** Protein-group quantity column candidates, in preference order. `PG.IBAQ` is the
 * historical Spectronaut export; modern long-form reports use `PG.Quantity` (MaxLFQ
 * or top-N peptide sum). Both are constant per (protein-group, run). */
const QUANTITY_COLUMNS = ['PG.IBAQ', 'PG.Quantity'];

/** Result of pivoting Spectronaut long-format data. */
interface PivotResult {
  proteinMap: Map<string, Map<string, number>>;
  sampleKeys: string[];
  sampleFileNames: Map<string, string>;
  organisms: Map<string, string>;
  /** Phase 16 SPC seed (D-01): earliest observed R.RunDate as ISO-8601, null if absent/unparseable. */
  earliestRunDate?: string | null;
  /** Phase 16 SPC seed (D-01): first non-empty R.InstrumentMethod cell, null if absent. */
  firstInstrumentMethod?: string | null;
}

/** Tries to parse a raw R.RunDate cell into a comparable epoch-ms; returns
 *  null when the value is empty or `Date` cannot parse it. */
function tryParseRunDate(raw: string): number | null {
  if (!raw || raw.trim().length === 0) return null;
  const t = Date.parse(raw);
  return Number.isNaN(t) ? null : t;
}

/** Pivots Spectronaut long-format rows into protein -> sample -> quantity map.
 * Filters CON__/REV__ proteins and rows exceeding q-value threshold.
 * Non-numeric q-values (e.g., 'Profiled', 'NaN') are treated as passing. */
function pivotSpectronaut(
  longDf: DG.DataFrame,
  quantityColName: string,
  qValueThreshold: number,
): PivotResult {
  const condCol = longDf.col('R.Condition')!;
  const replCol = longDf.col('R.Replicate')!;
  const protCol = longDf.col('PG.ProteinGroups')!;
  const ibaqCol = longDf.col(quantityColName)!;
  const qvalCol = longDf.col('EG.Qvalue');
  const fileCol = longDf.col('R.FileName');
  const orgCol = longDf.col('PG.Organisms');
  // Phase 16 SPC-seed columns (D-01): OPTIONAL; never added to REQUIRED_COLUMNS
  // so v1.3 inputs without them still parse.
  const runDateCol = longDf.col('R.RunDate');
  const instrMethodCol = longDf.col('R.InstrumentMethod');

  const proteinMap = new Map<string, Map<string, number>>();
  const sampleSet = new Set<string>();
  const sampleFileNames = new Map<string, string>();
  const organisms = new Map<string, string>();
  let earliestRunDateMs: number | null = null;
  let earliestRunDateIso: string | null = null;
  let firstInstrumentMethod: string | null = null;

  for (let i = 0; i < longDf.rowCount; i++) {
    // Q-value filter: numeric values > threshold are excluded; non-numeric pass
    if (qvalCol) {
      if (qvalCol.isNone(i)) {
        // null q-value: treat as passing
      } else {
        const raw = qvalCol.get(i);
        const qval = Number(raw);
        if (!isNaN(qval) && qval > qValueThreshold)
          continue;
      }
    }

    const protein = protCol.get(i) as string;
    if (!protein) continue;
    if (protein.startsWith('CON__') || protein.startsWith('REV__')) continue;

    const condition = String(condCol.get(i));
    const replicate = String(replCol.get(i));
    const sampleKey = `${condition}_${replicate}`;
    sampleSet.add(sampleKey);

    if (!proteinMap.has(protein))
      proteinMap.set(protein, new Map<string, number>());

    // First-encountered value wins (PG.IBAQ is constant per protein+sample)
    if (!proteinMap.get(protein)!.has(sampleKey) && !ibaqCol.isNone(i))
      proteinMap.get(protein)!.set(sampleKey, Number(ibaqCol.get(i)));

    if (fileCol && !sampleFileNames.has(sampleKey) && !fileCol.isNone(i))
      sampleFileNames.set(sampleKey, fileCol.get(i) as string);
    if (orgCol && !organisms.has(protein) && !orgCol.isNone(i))
      organisms.set(protein, orgCol.get(i) as string);

    // SPC seed: single-pass accumulators (Plan 16-04 Perf Decision Option A).
    if (runDateCol && !runDateCol.isNone(i)) {
      const raw = String(runDateCol.get(i));
      const ms = tryParseRunDate(raw);
      if (ms !== null && (earliestRunDateMs === null || ms < earliestRunDateMs)) {
        earliestRunDateMs = ms;
        earliestRunDateIso = new Date(ms).toISOString();
      }
    }
    if (instrMethodCol && firstInstrumentMethod === null && !instrMethodCol.isNone(i)) {
      const raw = String(instrMethodCol.get(i));
      if (raw.trim().length > 0) firstInstrumentMethod = raw;
    }
  }

  return {
    proteinMap,
    sampleKeys: Array.from(sampleSet).sort(),
    sampleFileNames,
    organisms,
    earliestRunDate: earliestRunDateIso,
    firstInstrumentMethod,
  };
}

/** Builds wide protein-by-sample DataFrame from pivot result. */
function buildWideDataFrame(result: PivotResult): DG.DataFrame {
  const proteins = Array.from(result.proteinMap.keys());
  const n = proteins.length;

  const proteinCol = DG.Column.fromStrings('PG.ProteinGroups', proteins);
  proteinCol.semType = SEMTYPE.PROTEIN_ID;

  const cols: DG.Column[] = [proteinCol];

  // Add organisms column if any were collected
  if (result.organisms.size > 0) {
    const orgValues = proteins.map((p) => result.organisms.get(p) ?? '');
    const orgCol = DG.Column.fromStrings('PG.Organisms', orgValues);
    cols.push(orgCol);
  }

  for (const sampleKey of result.sampleKeys) {
    const values = new Float32Array(n);
    for (let i = 0; i < n; i++) {
      const sampleMap = result.proteinMap.get(proteins[i])!;
      values[i] = sampleMap.has(sampleKey) ? sampleMap.get(sampleKey)! : DG.FLOAT_NULL;
    }
    const col = DG.Column.fromFloat32Array(sampleKey, values);
    col.semType = SEMTYPE.INTENSITY;
    if (result.sampleFileNames.has(sampleKey))
      col.setTag('spectronaut.fileName', result.sampleFileNames.get(sampleKey)!);
    cols.push(col);
  }

  return DG.DataFrame.fromColumns(cols);
}

/** Auto-populates experimental groups from R.Condition when exactly 2 conditions exist.
 * Groups use log2-prefixed column names for downstream analysis compatibility. */
function autoPopulateGroups(df: DG.DataFrame, sampleKeys: string[]): void {
  const conditionMap = new Map<string, string[]>();
  for (const key of sampleKeys) {
    const lastUnderscore = key.lastIndexOf('_');
    const condition = key.substring(0, lastUnderscore);
    if (!conditionMap.has(condition))
      conditionMap.set(condition, []);
    conditionMap.get(condition)!.push(`log2(${key})`);
  }

  const conditions = Array.from(conditionMap.keys());
  if (conditions.length === 2) {
    setGroups(df, {
      group1: {name: conditions[0], columns: conditionMap.get(conditions[0])!},
      group2: {name: conditions[1], columns: conditionMap.get(conditions[1])!},
    });
  }
}

/** Shared post-pivot tail for BOTH the text and streaming Spectronaut paths.
 *
 * Builds the wide DataFrame, derives the primary protein column, detects/applies
 * the log2 transform, sets the `proteomics.*` tags, and auto-populates the two
 * experimental groups. `parseSpectronautText` and `parseSpectronautStream` both
 * end by handing their `PivotResult` here, so their output is byte-identical in
 * shape, semTypes, log2 columns, tags, and groups (R3).
 *
 * NOTE: this is the verbatim former tail of `parseSpectronautText` extracted
 * unchanged — do not edit tag strings, semType assignments, the log2 branch, or
 * column construction here without re-locking the existing Spectronaut tests. */
async function finalizeSpectronaut(result: PivotResult): Promise<DG.DataFrame> {
  // Build wide DataFrame
  const df = buildWideDataFrame(result);

  // Parse semicolon-delimited protein groups into primary column
  addPrimaryColumnIfNeeded(df, 'PG.ProteinGroups', 'Primary Protein ID', SEMTYPE.PROTEIN_ID);

  // Get intensity column names (the sample columns)
  const intensityColNames = result.sampleKeys;

  // Detect log2 status
  const log2Status = detectLog2Status(df, intensityColNames);

  if (log2Status.isLog2) {
    // Pre-normalized: copy as log2 columns without transformation
    copyAsLog2Columns(df, intensityColNames);
  } else {
    // Raw: apply log2 transformation, keeping originals
    log2TransformColumns(df, intensityColNames);
  }

  // Spectronaut always performs its own normalization regardless of export format
  df.setTag('proteomics.preNormalized', 'true');

  // Tag source
  df.setTag('proteomics.source', 'spectronaut');

  // Auto-populate groups from R.Condition
  autoPopulateGroups(df, result.sampleKeys);

  await resolveGeneLabels(df);

  return df;
}

/** Mirrors duckdb `nullstr=['','NaN','NA']` + `TRY_CAST(... AS DOUBLE)` from
 * tools/spectronaut-aggregate.sql. Returns `null` for the duckdb null tokens
 * and for anything that is not a finite number; otherwise the parsed number.
 *
 * This is intentionally STRICTER than `pivotSpectronaut`'s loose
 * `Number(raw)`/`isNaN` (which treats `''` as `0` via `Number('')`): the
 * streaming path must match the duckdb D-04 golden oracle exactly. The two
 * paths only diverge on degenerate inputs (`''` q-value/quantity) that the
 * real Spectronaut exports and the synthetic fixture never contain — the text
 * path stays byte-identical because it does not route precursor reports. */
function tryCastDouble(s: string | null | undefined): number | null {
  if (s === null || s === undefined) return null;
  const t = s.trim();
  if (t === '' || t === 'NaN' || t === 'NA') return null;
  const n = Number(t);
  return Number.isFinite(n) ? n : null;
}

/** Outcome of feeding one already-split streamed data line to `handleFields`:
 *
 * - `'kept'`      — the row was aggregated into `agg` (counts toward `rowCount`).
 * - `'malformed'` — the row is structurally unparseable (too few tab-separated
 *                   fields). ONLY this category is surfaced to the user as
 *                   "malformed/unparseable line(s)".
 * - `'filtered'`  — the row was dropped by a correct, by-design duckdb-parity
 *                   filter (empty `PG.ProteinGroups`, CON__/REV__ decoy/
 *                   contaminant, numeric `EG.Qvalue` > threshold, or both
 *                   numeric casts null). This is NOT malformed and is SILENT —
 *                   exactly as the text path (`pivotSpectronaut`) is for the
 *                   identical drops.
 */
type LineOutcome = 'kept' | 'malformed' | 'filtered';

/** Per-(protein × condition × replicate) running aggregate accumulated by the
 * streaming path. Mirrors the duckdb GROUP BY: `max(quantity)`, `min(qvalue)`,
 * first-non-null `R.FileName` / `PG.Organisms`. */
interface AggRow {
  protein: string;
  condition: string;
  replicate: string;
  maxQuantity: number;
  minQvalue: number;
  fileName: string | null;
  organism: string | null;
}

/** Group-key field separator. ASCII Unit Separator (0x1F) is not a TSV-legal
 * character, so it cannot collide with any condition/replicate/protein value. */
const AGG_KEY_SEP = '\x1f';

/** Folds the streaming aggregate map into the same `PivotResult` shape that
 * `pivotSpectronaut` produces, so `finalizeSpectronaut` is path-agnostic.
 * `sampleKey` MUST be `${condition}_${replicate}` to match `pivotSpectronaut`
 * (L65) and `autoPopulateGroups`'s `lastIndexOf('_')` split. */
function aggToPivotResult(agg: Map<string, AggRow>): PivotResult {
  const proteinMap = new Map<string, Map<string, number>>();
  const sampleSet = new Set<string>();
  const sampleFileNames = new Map<string, string>();
  const organisms = new Map<string, string>();

  for (const row of agg.values()) {
    // Never-set groups (no surviving numeric quantity) carry -Infinity; skip
    // them exactly as duckdb would emit no quantity for an all-null group.
    if (!Number.isFinite(row.maxQuantity)) continue;

    const sampleKey = `${row.condition}_${row.replicate}`;
    sampleSet.add(sampleKey);

    if (!proteinMap.has(row.protein))
      proteinMap.set(row.protein, new Map<string, number>());
    proteinMap.get(row.protein)!.set(sampleKey, row.maxQuantity);

    if (row.fileName !== null && !sampleFileNames.has(sampleKey))
      sampleFileNames.set(sampleKey, row.fileName);
    if (row.organism !== null && !organisms.has(row.protein))
      organisms.set(row.protein, row.organism);
  }

  return {
    proteinMap,
    sampleKeys: Array.from(sampleSet).sort(),
    sampleFileNames,
    organisms,
  };
}

/** Streams a precursor/fragment-level Spectronaut long-format TSV `File` and
 * single-pass aggregates it down to the same wide protein-by-sample DataFrame
 * `parseSpectronautText` produces — without ever materializing the whole file
 * as a string (the V8 ~512 MB string ceiling is why a 2.6 GB report OOMs the
 * text path). Only state retained is the bounded carry-over line buffer and the
 * aggregate Map keyed by (protein × condition × replicate); raw input rows are
 * NEVER buffered (T-12-05).
 *
 * Filter / aggregate parity with tools/spectronaut-aggregate.sql:
 * - drop empty / `CON__` / `REV__` `PG.ProteinGroups`
 * - drop a row only if `tryCastDouble(EG.Qvalue) > threshold` (null/non-numeric
 *   passes — duckdb `IS NULL OR <= 0.01`)
 * - per group: `max(tryCastDouble(quantity))`, `min(tryCastDouble(EG.Qvalue))`,
 *   first-non-null `R.FileName` / `PG.Organisms`
 * - the DMD↔WT `R.Condition` flip in the SQL is REFERENCE-FILE-ONLY and is
 *   deliberately NOT ported here (RESEARCH Pitfall 1).
 *
 * Malformed lines (too few fields, or both numeric casts null) are skipped and
 * counted (duckdb `ignore_errors=true` parity), surfaced via the progress
 * message rather than aborting. A `TaskBarProgressIndicator` advances by
 * bytes-read and an explicit `setTimeout(0)` macrotask yield runs on a ~16 ms
 * cadence so the tab stays responsive (T-12-09; the stream-read await alone
 * does not yield on OS-buffered data). */
export async function parseSpectronautStream(
  file: File, qValueThreshold: number = 0.01,
): Promise<DG.DataFrame> {
  const pi = DG.TaskBarProgressIndicator.create('Streaming Spectronaut report...');
  try {
    const reader = file.stream()
      .pipeThrough(new TextDecoderStream('utf-8'))
      .getReader();

    const agg = new Map<string, AggRow>();
    let buffer = '';
    let headerParsed = false;
    let colIdx: Map<string, number> = new Map();
    let protI = -1; let condI = -1; let repI = -1;
    let qvalI = -1; let qtyI = -1; let fileI = -1; let orgI = -1;
    // Phase 16 SPC seed (D-01) — OPTIONAL columns, never required.
    let runDateI = -1; let instrMethodI = -1;
    let expectedFields = 0;
    let earliestRunDateMs: number | null = null;
    let earliestRunDateIso: string | null = null;
    let firstInstrumentMethod: string | null = null;

    let bytesRead = 0;
    let rowCount = 0;
    // Structurally unparseable lines (too few tab-separated fields). ONLY this
    // counter feeds the user-facing "malformed line(s)" message.
    let malformed = 0;
    // Rows dropped by a correct by-design duckdb-parity filter (CON__/REV__ or
    // numeric q > threshold). Tracked for completeness but SILENT to the user —
    // the text path (pivotSpectronaut) drops the identical rows with no message.
    let filtered = 0;
    let lastYield = performance.now();
    const fileSize = Math.max(1, file.size);

    /** Aggregate one already-split data line into `agg` (parity with the
     * duckdb WHERE + GROUP BY). Returns a discriminated outcome:
     * `'malformed'` for a structurally unparseable row (too few tab-separated
     * fields), `'filtered'` for a correct by-design duckdb-parity drop (empty
     * `PG.ProteinGroups`, CON__/REV__ decoy/contaminant, numeric `EG.Qvalue` >
     * threshold, or both numeric casts null), `'kept'` once the row is
     * aggregated. Empty-protein and both-casts-null are duckdb-silent drops
     * that the text path (`pivotSpectronaut`) also drops without messaging —
     * categorizing them `'filtered'` keeps streaming↔text user-message parity.
     * The predicates below are duckdb-parity-correct and UNCHANGED — only the
     * outcome value carries the malformed-vs-filtered distinction. */
    const handleFields = (f: string[]): LineOutcome => {
      if (f.length < expectedFields) return 'malformed';
      const protein = f[protI];
      // Drop empty / decoy / contaminant protein groups — all silent, matching
      // `pivotSpectronaut`'s `if (!protein) continue` + CON__/REV__ skip.
      if (!protein) return 'filtered';
      if (protein.startsWith('CON__') || protein.startsWith('REV__')) return 'filtered';

      const q = qvalI >= 0 ? tryCastDouble(f[qvalI]) : null;
      // duckdb: (qvalue IS NULL OR qvalue <= threshold) — null/non-numeric pass.
      if (q !== null && q > qValueThreshold) return 'filtered';

      const qty = tryCastDouble(f[qtyI]);
      // Both numeric casts null → unusable for aggregation; drop silently
      // (duckdb `ignore_errors` semantics — the text path likewise omits this row).
      if (qty === null && q === null) return 'filtered';

      const condition = f[condI] ?? '';
      const replicate = f[repI] ?? '';
      const key = protein + AGG_KEY_SEP + condition + AGG_KEY_SEP + replicate;
      let row = agg.get(key);
      if (row === undefined) {
        row = {
          protein, condition, replicate,
          maxQuantity: -Infinity, minQvalue: Infinity,
          fileName: null, organism: null,
        };
        agg.set(key, row);
      }
      if (qty !== null && qty > row.maxQuantity) row.maxQuantity = qty;
      if (q !== null && q < row.minQvalue) row.minQvalue = q;
      if (row.fileName === null && fileI >= 0) {
        const fn = f[fileI];
        if (fn) row.fileName = fn;
      }
      if (row.organism === null && orgI >= 0) {
        const og = f[orgI];
        if (og) row.organism = og;
      }
      // Phase 16 SPC seed (D-01) — accumulate inline so the streaming-path adds
      // no second pass over the file.
      if (runDateI >= 0) {
        const rd = f[runDateI];
        if (rd) {
          const ms = tryParseRunDate(rd);
          if (ms !== null && (earliestRunDateMs === null || ms < earliestRunDateMs)) {
            earliestRunDateMs = ms;
            earliestRunDateIso = new Date(ms).toISOString();
          }
        }
      }
      if (instrMethodI >= 0 && firstInstrumentMethod === null) {
        const im = f[instrMethodI];
        if (im && im.trim().length > 0) firstInstrumentMethod = im;
      }
      return 'kept';
    };

    /** Parse the first complete line as the header; resolve column indices and
     * the single quantity column once, matching the text-path validation. */
    const parseHeader = (line: string): void => {
      const cols = line.split('\t');
      colIdx = new Map();
      for (let i = 0; i < cols.length; i++)
        colIdx.set(cols[i], i);
      for (const colName of REQUIRED_COLUMNS) {
        if (!colIdx.has(colName))
          throw new Error(`Missing required Spectronaut column: ${colName}`);
      }
      const quantityColName = QUANTITY_COLUMNS.find((n) => colIdx.has(n));
      if (!quantityColName) {
        throw new Error(`Missing protein-group quantity column ` +
          `(expected one of ${QUANTITY_COLUMNS.join(', ')})`);
      }
      protI = colIdx.get('PG.ProteinGroups')!;
      condI = colIdx.get('R.Condition')!;
      repI = colIdx.get('R.Replicate')!;
      qtyI = colIdx.get(quantityColName)!;
      qvalI = colIdx.has('EG.Qvalue') ? colIdx.get('EG.Qvalue')! : -1;
      fileI = colIdx.has('R.FileName') ? colIdx.get('R.FileName')! : -1;
      orgI = colIdx.has('PG.Organisms') ? colIdx.get('PG.Organisms')! : -1;
      // Phase 16 SPC seed: optional indices, contribute to expectedFields ONLY
      // when present so v1.3 headers without them still parse.
      runDateI = colIdx.has('R.RunDate') ? colIdx.get('R.RunDate')! : -1;
      instrMethodI = colIdx.has('R.InstrumentMethod') ? colIdx.get('R.InstrumentMethod')! : -1;
      // Highest index we read; a shorter split means a truncated/malformed row.
      expectedFields = Math.max(
        protI, condI, repI, qtyI, qvalI, fileI, orgI, runDateI, instrMethodI,
      ) + 1;
      headerParsed = true;
    };

    while (true) {
      const {value, done} = await reader.read();
      if (done) break;
      bytesRead += value.length;
      buffer += value;

      let nl: number;
      while ((nl = buffer.indexOf('\n')) >= 0) {
        let line = buffer.slice(0, nl);
        buffer = buffer.slice(nl + 1);
        if (line.endsWith('\r')) line = line.slice(0, -1);
        if (!headerParsed) {
          if (line.length === 0) continue;
          parseHeader(line);
          continue;
        }
        if (line.length === 0) continue;
        switch (handleFields(line.split('\t'))) {
        case 'kept': rowCount++; break;
        case 'malformed': malformed++; break;
        case 'filtered': filtered++; break;
        }
      }

      // ~16 ms wall-clock yield cadence (RESEARCH Q3 resolution / D-02). The
      // explicit setTimeout(0) macrotask yield is REQUIRED — await reader.read()
      // does not reliably yield on OS-buffered data (T-12-09).
      const now = performance.now();
      if (now - lastYield >= 16) {
        const pct = Math.min(99, (bytesRead / fileSize) * 100);
        // Only the genuine-malformed count is surfaced; by-design-filtered rows
        // (CON__/REV__, q > threshold) stay silent — text-path-consistent.
        const msg = malformed > 0
          ? `${rowCount} rows (${malformed} malformed skipped)`
          : `${rowCount} rows`;
        pi.update(pct, msg);
        await new Promise<void>((r) => setTimeout(r, 0));
        lastYield = performance.now();
      }
    }

    // Flush a trailing partial line (file without a final newline).
    if (buffer.length > 0) {
      let line = buffer;
      if (line.endsWith('\r')) line = line.slice(0, -1);
      if (line.length > 0) {
        if (!headerParsed) {
          parseHeader(line);
        } else {
          switch (handleFields(line.split('\t'))) {
          case 'kept': rowCount++; break;
          case 'malformed': malformed++; break;
          case 'filtered': filtered++; break;
          }
        }
      }
    }

    if (!headerParsed)
      throw new Error('Missing required Spectronaut column: R.Condition');

    // Genuine-malformed only. By-design-filtered rows (CON__/REV__, q >
    // threshold) are intentionally SILENT — pivotSpectronaut drops the
    // identical rows with no message; matching that silence is the locked,
    // path-consistent UX decision. `filtered` is tracked but never surfaced.
    void filtered;
    if (malformed > 0)
      grok.shell.info(`Spectronaut import: skipped ${malformed} malformed line(s)`);

    const pivot = aggToPivotResult(agg);
    pivot.earliestRunDate = earliestRunDateIso;
    pivot.firstInstrumentMethod = firstInstrumentMethod;
    const df = await finalizeSpectronaut(pivot);
    // Phase 16 D-01: write the `proteomics.spc_run_meta_seed` tag when either
    // capture was made (streaming path).
    applySpcRunMetaSeed(df, pivot);
    return df;
  } finally {
    pi.close();
  }
}

/** Parses Spectronaut long-format TSV text into a wide protein-by-sample DataFrame.
 *
 * Applies:
 * - Q-value filtering (default threshold 0.01; non-numeric q-values pass)
 * - CON__/REV__ contaminant/decoy filtering
 * - Long-to-wide pivot using PG.IBAQ (preferred) or PG.Quantity (newer exports)
 *   as protein-level intensity
 * - Semantic type assignment (PROTEIN_ID, INTENSITY)
 * - Log2 transformation (or copy if pre-normalized)
 * - Auto-group population from R.Condition when exactly 2 conditions
 * - Pre-normalization detection and tagging */
export async function parseSpectronautText(text: string, qValueThreshold: number = 0.01): Promise<DG.DataFrame> {
  // Parse long-format TSV
  const longDf = DG.DataFrame.fromCsv(text, {delimiter: '\t'});

  // Validate required columns
  for (const colName of REQUIRED_COLUMNS) {
    if (!longDf.col(colName))
      throw new Error(`Missing required Spectronaut column: ${colName}`);
  }

  // Pick the protein-group quantity column (PG.IBAQ preferred, PG.Quantity accepted).
  const quantityColName = QUANTITY_COLUMNS.find((n) => longDf.col(n) !== null);
  if (!quantityColName) {
    throw new Error(`Missing protein-group quantity column ` +
      `(expected one of ${QUANTITY_COLUMNS.join(', ')})`);
  }

  // Pivot long-to-wide with filtering
  const result = pivotSpectronaut(longDf, quantityColName, qValueThreshold);

  // Shared post-pivot tail (identical for the streaming path).
  const df = await finalizeSpectronaut(result);
  // Phase 16 D-01: write the `proteomics.spc_run_meta_seed` tag when either
  // capture was made (text path).
  applySpcRunMetaSeed(df, result);
  return df;
}

/** Phase 16 D-01: Attach the SPC run-meta seed to the parsed df when either
 *  R.RunDate or R.InstrumentMethod was captured. Plan 16-05's Annotate
 *  Experiment dialog reads this seed to prefill its two new inputs. */
function applySpcRunMetaSeed(df: DG.DataFrame, result: PivotResult): void {
  const hasDate = result.earliestRunDate != null;
  const hasInstr = result.firstInstrumentMethod != null;
  if (!hasDate && !hasInstr) return;
  df.setTag('proteomics.spc_run_meta_seed', JSON.stringify({
    instrument_id: result.firstInstrumentMethod ?? '',
    acquisition_datetime: result.earliestRunDate ?? '',
  }));
}

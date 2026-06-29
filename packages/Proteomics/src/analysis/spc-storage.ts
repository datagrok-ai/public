/**
 * Phase 16 SPC AppData I/O layer.
 *
 * Sole owner of `System:AppData/Proteomics/spc/runs.csv`,
 * `runs-meta.json`, and `baseline-<slug>.json`. The math layer in
 * `./spc.ts` stays AppData-free; the dashboard in `viewers/spc-dashboard.ts`
 * stays I/O-free. This module is the only file in the package that touches
 * the SPC AppData paths.
 *
 * ## Phase 17 additivity contract (CONTEXT.md D-02)
 *
 * Phase 17 reads runs.csv with the SAME 13-column shape and appends a
 * `campaign_id` column when saving a run into a campaign. Existing rows
 * back-fill with NULL — no migration script. The Phase 16 dashboard ignores
 * `campaign_id` if present so a v1.4 reader stays compatible with v1.5
 * writers.
 *
 * ## Concurrency
 *
 * Single-writer assumption per CONTEXT.md D-02. No file lock. A second
 * concurrent appendRun would silently overwrite — Cytokinetics's usage
 * pattern (one Datagrok user per machine) makes this acceptable for v1.4.
 *
 * ## Test root override
 *
 * Every public function accepts an optional `{root}` parameter pointing at
 * an alternate `System:AppData/Proteomics/spc-test/<uuid>` sandbox for tests
 * — production code never passes `root`, so the default `RUNS_DIR_DEFAULT`
 * is used.
 */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {BaselineSnapshot, SpcRulesEnabled, SpcMetrics, MetricBaseline} from './spc';

export type {BaselineSnapshot, SpcRulesEnabled, SpcMetrics, MetricBaseline};

export const RUNS_DIR_DEFAULT = 'System:AppData/Proteomics/spc';
export const RUNS_CSV_PATH = 'System:AppData/Proteomics/spc/runs.csv';
export const RUNS_META_PATH = 'System:AppData/Proteomics/spc/runs-meta.json';
export const RUNS_CSV_SCHEMA_VERSION = '1';

export const RUNS_CSV_COLUMNS: readonly string[] = [
  'run_id',
  'instrument_id',
  'acquisition_datetime',
  'run_label',
  'median_intensity',
  'missing_pct',
  'control_corr',
  'protein_count',
  'status',
  'rules_tripped',
  'source_project_id',
  'source_df_name',
  'computed_at',
];

export interface SpcRunRow {
  run_id: string;
  instrument_id: string;
  acquisition_datetime: string;
  run_label: string;
  median_intensity: number;
  missing_pct: number;
  control_corr: number;
  protein_count: number;
  status: 'pass' | 'flagged' | 'out_of_control';
  rules_tripped: string[];
  source_project_id: string | null;
  source_df_name: string;
  computed_at: string;
}

export interface IterationTraceEntry {
  iter: number;
  dropped: string[];
}

export interface LockedBaseline {
  instrument_id: string;
  locked_at: string;
  included_run_ids: string[];
  excluded_run_ids: string[];
  iteration_trace: IterationTraceEntry[];
  metrics: BaselineSnapshot;
  rules_enabled: SpcRulesEnabled;
}

export interface StorageOpts {
  /** Override the AppData root (used by tests; production callers omit). */
  root?: string;
}

const SCHEMA_VERSION_KEY = 'spc.runs.schema_version';

function root(opts?: StorageOpts): string {
  return (opts && opts.root) ? opts.root : RUNS_DIR_DEFAULT;
}

function runsCsvPath(opts?: StorageOpts): string {
  return `${root(opts)}/runs.csv`;
}

function runsMetaPath(opts?: StorageOpts): string {
  return `${root(opts)}/runs-meta.json`;
}

function baselinePath(instrumentId: string, opts?: StorageOpts): string {
  return `${root(opts)}/baseline-${slugifyInstrumentId(instrumentId)}.json`;
}

// =====================================================================================
// Slugify — T-16-01 mitigation. Mirrors publishing/publish-state.ts:slugifyTarget but
// PRESERVES case so instrument names like 'QExactive-01' round-trip exactly.
// =====================================================================================

export function slugifyInstrumentId(raw: string): string {
  if (raw == null) return 'unnamed';
  let s = String(raw)
    .replace(/[^A-Za-z0-9._-]+/g, '-')
    .replace(/-{2,}/g, '-')
    .replace(/^[-.]+|[-.]+$/g, '')
    .slice(0, 64)
    .replace(/[-.]+$/g, '');
  if (s.length === 0) s = 'unnamed';
  return s;
}

// =====================================================================================
// CSV helpers — small and explicit so the upsert-by-key pattern stays simple.
// DG.DataFrame.fromObjects + toCsv is avoided because it muddies number↔string
// type inference on the round-trip.
// =====================================================================================

function csvEscape(value: string): string {
  if (value.indexOf(',') < 0 && value.indexOf('"') < 0 && value.indexOf('\n') < 0)
    return value;
  return '"' + value.replace(/"/g, '""') + '"';
}

function parseCsvRow(line: string): string[] {
  const out: string[] = [];
  let i = 0;
  while (i < line.length) {
    if (line[i] === '"') {
      i++;
      let cell = '';
      while (i < line.length) {
        if (line[i] === '"') {
          if (line[i + 1] === '"') { cell += '"'; i += 2; }
          else { i++; break; }
        } else { cell += line[i]; i++; }
      }
      out.push(cell);
      if (line[i] === ',') i++;
    } else {
      let cell = '';
      while (i < line.length && line[i] !== ',') { cell += line[i]; i++; }
      out.push(cell);
      if (line[i] === ',') i++;
    }
  }
  return out;
}

function parseCsv(text: string): string[][] {
  const rows: string[][] = [];
  // Strip trailing newline; rows separated by \n. CSV cells may include \n if
  // they're quoted, but the SPC schema never carries embedded newlines.
  const lines = text.split('\n').filter((l) => l.length > 0);
  for (const line of lines) rows.push(parseCsvRow(line));
  return rows;
}

function serializeRow(row: SpcRunRow): string {
  const ordered = [
    row.run_id,
    row.instrument_id,
    row.acquisition_datetime,
    row.run_label,
    String(row.median_intensity),
    String(row.missing_pct),
    String(row.control_corr),
    String(row.protein_count),
    row.status,
    JSON.stringify(row.rules_tripped),
    row.source_project_id == null ? '' : row.source_project_id,
    row.source_df_name,
    row.computed_at,
  ];
  return ordered.map(csvEscape).join(',');
}

function deserializeRow(cells: string[]): SpcRunRow {
  const idx = (n: string) => RUNS_CSV_COLUMNS.indexOf(n);
  const get = (n: string) => cells[idx(n)] ?? '';
  let rulesTripped: string[] = [];
  try {
    const raw = get('rules_tripped');
    if (raw) rulesTripped = JSON.parse(raw) as string[];
  } catch { /* malformed → empty */ }
  const sourceProjectId = get('source_project_id');
  return {
    run_id: get('run_id'),
    instrument_id: get('instrument_id'),
    acquisition_datetime: get('acquisition_datetime'),
    run_label: get('run_label'),
    median_intensity: Number(get('median_intensity')),
    missing_pct: Number(get('missing_pct')),
    control_corr: Number(get('control_corr')),
    protein_count: Number(get('protein_count')),
    status: get('status') as SpcRunRow['status'],
    rules_tripped: rulesTripped,
    source_project_id: sourceProjectId.length === 0 ? null : sourceProjectId,
    source_df_name: get('source_df_name'),
    computed_at: get('computed_at'),
  };
}

function emptyCsv(): string {
  return RUNS_CSV_COLUMNS.map(csvEscape).join(',') + '\n';
}

function rowsToCsv(rows: SpcRunRow[]): string {
  const lines = [RUNS_CSV_COLUMNS.map(csvEscape).join(',')];
  for (const r of rows) lines.push(serializeRow(r));
  return lines.join('\n') + '\n';
}

async function readRowsFromCsv(opts?: StorageOpts): Promise<SpcRunRow[]> {
  const path = runsCsvPath(opts);
  if (!await grok.dapi.files.exists(path)) return [];
  const text = await grok.dapi.files.readAsText(path);
  const parsed = parseCsv(text);
  if (parsed.length <= 1) return [];
  // First row is the header. Skip it.
  return parsed.slice(1).map(deserializeRow);
}

// =====================================================================================
// Schema versioning
// =====================================================================================

export async function readRunsMeta(opts?: StorageOpts): Promise<Record<string, string>> {
  const path = runsMetaPath(opts);
  if (!await grok.dapi.files.exists(path)) return {};
  try {
    return JSON.parse(await grok.dapi.files.readAsText(path)) as Record<string, string>;
  } catch (err) {
    // Storage layer is shell-free; user-facing warnings live in the viewer.
    // eslint-disable-next-line no-console
    console.warn(`[spc-storage] failed to parse ${path}: ${err}`);
    return {};
  }
}

export async function writeRunsMeta(
  meta: Record<string, string>, opts?: StorageOpts,
): Promise<void> {
  await grok.dapi.files.writeAsText(runsMetaPath(opts), JSON.stringify(meta, null, 2));
}

export async function readSchemaVersion(opts?: StorageOpts): Promise<string | null> {
  const meta = await readRunsMeta(opts);
  return meta[SCHEMA_VERSION_KEY] ?? null;
}

export async function writeSchemaVersion(version: string, opts?: StorageOpts): Promise<void> {
  const meta = await readRunsMeta(opts);
  meta[SCHEMA_VERSION_KEY] = version;
  await writeRunsMeta(meta, opts);
}

// =====================================================================================
// runs.csv upsert (primary write path)
// =====================================================================================

function genRunId(): string {
  if (typeof crypto !== 'undefined' && (crypto as any).randomUUID)
    return (crypto as any).randomUUID();
  // Fallback for environments without crypto.randomUUID.
  return 'r-' + Math.random().toString(16).slice(2) + '-' + Date.now().toString(16);
}

export interface UpsertInput {
  run_id?: string;
  instrument_id: string;
  acquisition_datetime: string;
  run_label: string;
  median_intensity: number;
  missing_pct: number;
  control_corr: number;
  protein_count: number;
  status: SpcRunRow['status'];
  rules_tripped: string[];
  source_project_id: string | null;
  source_df_name: string;
  computed_at: string;
}

/**
 * Upserts on the composite key (instrument_id, acquisition_datetime). Second
 * invocation with the same key OVERWRITES the prior row (preserving its
 * run_id). First-ever invocation creates runs.csv and pins the schema
 * version in runs-meta.json.
 */
export async function upsertRun(
  input: UpsertInput, opts?: StorageOpts,
): Promise<SpcRunRow> {
  const rows = await readRowsFromCsv(opts);
  const idx = rows.findIndex(
    (r) => r.instrument_id === input.instrument_id &&
      r.acquisition_datetime === input.acquisition_datetime,
  );
  let persisted: SpcRunRow;
  if (idx >= 0) {
    // Preserve the existing run_id on overwrite.
    persisted = {
      ...rows[idx],
      ...input,
      run_id: rows[idx].run_id,
      rules_tripped: input.rules_tripped,
      source_project_id: input.source_project_id,
    } as SpcRunRow;
    rows[idx] = persisted;
  } else {
    persisted = {
      run_id: input.run_id ?? genRunId(),
      instrument_id: input.instrument_id,
      acquisition_datetime: input.acquisition_datetime,
      run_label: input.run_label,
      median_intensity: input.median_intensity,
      missing_pct: input.missing_pct,
      control_corr: input.control_corr,
      protein_count: input.protein_count,
      status: input.status,
      rules_tripped: input.rules_tripped,
      source_project_id: input.source_project_id,
      source_df_name: input.source_df_name,
      computed_at: input.computed_at,
    };
    rows.push(persisted);
  }
  await grok.dapi.files.writeAsText(runsCsvPath(opts), rowsToCsv(rows));
  // First-ever write: pin the schema version. Idempotent on subsequent calls
  // because writeSchemaVersion merges into existing runs-meta.json.
  const version = await readSchemaVersion(opts);
  if (version !== RUNS_CSV_SCHEMA_VERSION)
    await writeSchemaVersion(RUNS_CSV_SCHEMA_VERSION, opts);
  return persisted;
}

/** Alias for Plan 16-05's wireup, which calls `appendRun` per the plan contract. */
export const appendRun = upsertRun;

// =====================================================================================
// runs.csv read (with optional instrument filter; always sorted by acquisition_datetime ASC)
// =====================================================================================

export async function loadRuns(
  instrumentId?: string, opts?: StorageOpts,
): Promise<SpcRunRow[]> {
  const rows = await readRowsFromCsv(opts);
  const filtered = instrumentId
    ? rows.filter((r) => r.instrument_id === instrumentId)
    : rows;
  // Pitfall 7 lever: ALWAYS sort by acquisition_datetime, never by computed_at
  // or insertion order. ISO-8601 lex-sorts correctly.
  filtered.sort((a, b) => a.acquisition_datetime.localeCompare(b.acquisition_datetime));
  return filtered;
}

// =====================================================================================
// Per-instrument baseline JSON
// =====================================================================================

export async function loadBaseline(
  instrumentId: string, opts?: StorageOpts,
): Promise<LockedBaseline | null> {
  const path = baselinePath(instrumentId, opts);
  if (!await grok.dapi.files.exists(path)) return null;
  try {
    return JSON.parse(await grok.dapi.files.readAsText(path)) as LockedBaseline;
  } catch (err) {
    // eslint-disable-next-line no-console
    console.warn(`[spc-storage] failed to parse baseline at ${path}: ${err}`);
    return null;
  }
}

export async function saveBaseline(
  instrumentId: string, baseline: LockedBaseline, opts?: StorageOpts,
): Promise<void> {
  await grok.dapi.files.writeAsText(
    baselinePath(instrumentId, opts),
    JSON.stringify(baseline, null, 2),
  );
}

// =====================================================================================
// Iterative outlier removal (Pitfall 6 lever, capped at 2 iterations).
// =====================================================================================

const MAX_BASELINE_ITERATIONS = 2;

function meanSdValues(values: number[]): {mean: number; sd: number} {
  const n = values.length;
  if (n === 0) return {mean: NaN, sd: NaN};
  let sum = 0;
  for (const v of values) sum += v;
  const mean = sum / n;
  let sse = 0;
  for (const v of values) {
    const d = v - mean;
    sse += d * d;
  }
  // Sample SD (n-1 denominator) for Shewhart 3σ consistency.
  const sd = n > 1 ? Math.sqrt(sse / (n - 1)) : 0;
  return {mean, sd};
}

export function iterativeOutlierRemoval(
  series: Array<{run_id: string; value: number}>,
  maxIter: number = MAX_BASELINE_ITERATIONS,
): {
  retained_run_ids: string[];
  excluded_run_ids: string[];
  iteration_trace: IterationTraceEntry[];
} {
  const retained = series.slice();
  const excluded: Array<{run_id: string; value: number}> = [];
  const trace: IterationTraceEntry[] = [];

  if (retained.length < 3) {
    return {
      retained_run_ids: retained.map((s) => s.run_id),
      excluded_run_ids: [],
      iteration_trace: [],
    };
  }

  for (let iter = 1; iter <= maxIter; iter++) {
    const vals = retained.map((s) => s.value).filter((v) => Number.isFinite(v));
    const {mean, sd} = meanSdValues(vals);
    if (!Number.isFinite(sd) || sd === 0) break;
    const dropped: typeof retained = [];
    for (let i = retained.length - 1; i >= 0; i--) {
      const v = retained[i].value;
      if (!Number.isFinite(v)) continue;
      const z = (v - mean) / sd;
      if (Math.abs(z) > 3) {
        dropped.push(retained[i]);
        excluded.push(retained[i]);
        retained.splice(i, 1);
      }
    }
    if (dropped.length === 0) break;
    trace.push({iter, dropped: dropped.map((d) => d.run_id)});
  }

  return {
    retained_run_ids: retained.map((s) => s.run_id),
    excluded_run_ids: excluded.map((s) => s.run_id),
    iteration_trace: trace,
  };
}

// =====================================================================================
// Baseline statistics (mean+sd per metric over retained rows)
// =====================================================================================

const METRIC_KEYS: Array<keyof BaselineSnapshot> = [
  'median_intensity', 'missing_pct', 'control_corr', 'protein_count',
];

function metricBaseline(values: number[]): MetricBaseline {
  const clean = values.filter((v) => Number.isFinite(v));
  if (clean.length < 2) return {mean: NaN, sd: NaN};
  const {mean, sd} = meanSdValues(clean);
  return {mean, sd};
}

export function computeBaselineStats(retainedRows: SpcRunRow[]): BaselineSnapshot {
  const snapshot = {} as BaselineSnapshot;
  for (const m of METRIC_KEYS)
    snapshot[m] = metricBaseline(retainedRows.map((r) => r[m] as number));
  return snapshot;
}

// =====================================================================================
// Reuse a DG.DataFrame helper for downstream viewers — small but keeps the dashboard
// from re-implementing the CSV walk. Plan 16-06 calls this in createSpcChartPanel.
// =====================================================================================

export function runsToDataFrame(rows: SpcRunRow[]): DG.DataFrame {
  const get = <T>(fn: (r: SpcRunRow) => T) => rows.map(fn);
  const cols = [
    DG.Column.fromStrings('run_id', get((r) => r.run_id)),
    DG.Column.fromStrings('instrument_id', get((r) => r.instrument_id)),
    DG.Column.fromList(DG.COLUMN_TYPE.DATE_TIME, 'acquisition_datetime',
      get((r) => Date.parse(r.acquisition_datetime))),
    DG.Column.fromStrings('run_label', get((r) => r.run_label)),
    DG.Column.fromFloat32Array('median_intensity', new Float32Array(get((r) => r.median_intensity))),
    DG.Column.fromFloat32Array('missing_pct', new Float32Array(get((r) => r.missing_pct))),
    DG.Column.fromFloat32Array('control_corr', new Float32Array(get((r) => r.control_corr))),
    DG.Column.fromFloat32Array('protein_count', new Float32Array(get((r) => r.protein_count))),
    DG.Column.fromStrings('status', get((r) => r.status)),
    DG.Column.fromStrings('rules_tripped', get((r) => JSON.stringify(r.rules_tripped))),
    DG.Column.fromStrings('source_project_id', get((r) => r.source_project_id ?? '')),
    DG.Column.fromStrings('source_df_name', get((r) => r.source_df_name)),
    DG.Column.fromStrings('computed_at', get((r) => r.computed_at)),
  ];
  return DG.DataFrame.fromColumns(cols);
}

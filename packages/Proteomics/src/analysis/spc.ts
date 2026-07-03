/**
 * Phase 16 SPC math + tag helpers.
 *
 * Pure TS layer — zero shell / UI side effects. AppData I/O lives in
 * `src/analysis/spc-storage.ts` (Plan 16-03); viewers live in
 * `src/viewers/spc-dashboard.ts` (Plan 16-06). Public-domain math: Nelson 1984
 * for the eight rules and Wheeler 1992 for the n=2 MR-chart constants.
 *
 * Belt-and-braces (Phase 15 D-05): {@link setSpcStatus} persists a single-row
 * string column `~spc_metrics_meta` carrying `JSON.stringify(metrics)`. One
 * column keeps the published-project allowlist small and survives the
 * serializer-strip pattern that loses df tags.
 */
import * as DG from 'datagrok-api/dg';

import {GroupAssignment} from './experiment-setup';

/** Run-identity payload persisted in `proteomics.spc_run_meta` JSON tag. */
export interface RunMeta {
  instrument_id: string;
  acquisition_datetime: string;
}

/** Per-run rollup written to AppData runs.csv and surfaced on every analyzed df. */
export interface SpcMetrics {
  median_intensity: number;
  missing_pct: number;
  control_corr: number;
  protein_count: number;
  sample_count: number;
  computed_at: string;
}

export type SpcStatus = 'pass' | 'flagged' | 'out_of_control';

/** One run's evaluated rules, with status classification per D-04. */
export interface SpcRuleResult {
  status: SpcStatus;
  rulesTripped: string[];
}

export type NelsonRuleId = 'nelson_1' | 'nelson_2' | 'nelson_3' | 'nelson_4'
                          | 'nelson_5' | 'nelson_6' | 'nelson_7' | 'nelson_8';

/** Per-metric enabled flags (used by the simple form of evaluateNelsonRules). */
export type FlatRulesEnabled = Record<NelsonRuleId, boolean>;

/** Four-metric rules-enabled grid persisted in baseline-<instrument_id>.json. */
export interface SpcRulesEnabled {
  median_intensity: FlatRulesEnabled;
  missing_pct: FlatRulesEnabled;
  control_corr: FlatRulesEnabled;
  protein_count: FlatRulesEnabled;
}

/** Mean + standard deviation for one metric over a locked baseline window. */
export interface MetricBaseline {
  mean: number;
  sd: number;
}

/** Four-metric baseline snapshot stored on baseline-<instrument_id>.json. */
export interface BaselineSnapshot {
  median_intensity: MetricBaseline;
  missing_pct: MetricBaseline;
  control_corr: MetricBaseline;
  protein_count: MetricBaseline;
}

/** Result of iteratively pruning 3σ outliers from a baseline window. */
export interface IterativeBaselineResult {
  mean: number;
  sd: number;
  includedCount: number;
  excludedIndices: number[];
  iterationTrace: Array<{iter: number; dropped: number[]}>;
}

/** Wheeler 1992 control-chart constants for n=2 (used by Plan 16-06 MR-chart). */
export const WHEELER_D2 = 1.128;
export const WHEELER_D4 = 3.267;

/** Phase 16 explicit ban (Pitfall 5) — Spectronaut Candidates carry pre-computed DE
 *  with no per-sample intensities, so SPC cannot be computed from them. */
const SPECTRONAUT_CANDIDATES_REFUSAL =
  'SPC requires per-sample intensities. Re-import this analysis from the Spectronaut PG ' +
  'report (not the Candidates report) to compute SPC.';

const BELT_AND_BRACES_COLUMN = '~spc_metrics_meta';

// =====================================================================================
// Run-meta tag helpers (mirror setGroups/getGroups in experiment-setup.ts).
// =====================================================================================

export function setRunMeta(df: DG.DataFrame, meta: RunMeta): void {
  df.setTag('proteomics.spc_run_meta', JSON.stringify(meta));
}

export function getRunMeta(df: DG.DataFrame): RunMeta | null {
  const raw = df.getTag('proteomics.spc_run_meta');
  if (!raw) return null;
  try {
    return JSON.parse(raw) as RunMeta;
  } catch {
    return null;
  }
}

// =====================================================================================
// Source-gate (Pitfall 5 — Spectronaut Candidates hard refusal).
// =====================================================================================

export function assertSpcEligible(df: DG.DataFrame): void {
  const source = df.getTag('proteomics.source');
  if (source === 'spectronaut-candidates')
    throw new Error(SPECTRONAUT_CANDIDATES_REFUSAL);
}

// =====================================================================================
// Metric math (NaN-tolerant; matches qc-computations.ts semantics).
// =====================================================================================

function isMissing(col: DG.Column, i: number): boolean {
  if (col.isNone(i)) return true;
  const v = col.get(i);
  if (typeof v === 'number' && Number.isNaN(v)) return true;
  return false;
}

function median(values: number[]): number {
  if (values.length === 0) return NaN;
  const sorted = values.slice().sort((a, b) => a - b);
  const n = sorted.length;
  const mid = Math.floor(n / 2);
  return (n % 2 === 0) ? (sorted[mid - 1] + sorted[mid]) / 2 : sorted[mid];
}

function pearson(xs: number[], ys: number[]): number {
  const n = xs.length;
  if (n < 3) return NaN;
  let sumX = 0, sumY = 0;
  for (let i = 0; i < n; i++) { sumX += xs[i]; sumY += ys[i]; }
  const mx = sumX / n;
  const my = sumY / n;
  let num = 0, dx2 = 0, dy2 = 0;
  for (let i = 0; i < n; i++) {
    const dx = xs[i] - mx;
    const dy = ys[i] - my;
    num += dx * dy;
    dx2 += dx * dx;
    dy2 += dy * dy;
  }
  const denom = Math.sqrt(dx2 * dy2);
  return denom === 0 ? NaN : num / denom;
}

function resolveCols(df: DG.DataFrame, names: string[]): DG.Column[] {
  const out: DG.Column[] = [];
  for (const name of names) {
    const c = df.col(name);
    if (c !== null) out.push(c);
  }
  return out;
}

export function computeSpcMetrics(
  df: DG.DataFrame, groups: GroupAssignment, runMeta: RunMeta,
): SpcMetrics {
  // Use runMeta so a future refactor can persist it alongside the metrics
  // without breaking the public signature; intentionally not stored here.
  void runMeta;

  const g1Cols = resolveCols(df, groups.group1.columns);
  const g2Cols = resolveCols(df, groups.group2.columns);
  const intensityCols = g1Cols.concat(g2Cols);
  const sampleCount = groups.group1.columns.length + groups.group2.columns.length;

  // Median: flat array of all non-null intensity cells.
  const allValues: number[] = [];
  const rowCount = df.rowCount;
  let missingCells = 0;
  let totalCells = 0;
  let rowsWithData = 0;
  for (let r = 0; r < rowCount; r++) {
    let rowHasData = false;
    for (const col of intensityCols) {
      totalCells++;
      if (isMissing(col, r)) {
        missingCells++;
      } else {
        allValues.push(col.get(r) as number);
        rowHasData = true;
      }
    }
    if (rowHasData) rowsWithData++;
  }

  const medianIntensity = median(allValues);
  const missingPct = totalCells === 0 ? NaN : (missingCells / totalCells) * 100;

  // Pairwise Pearson over Group 1 controls.
  let controlCorr = NaN;
  if (g1Cols.length >= 2) {
    let sum = 0;
    let pairCount = 0;
    for (let i = 0; i < g1Cols.length; i++) {
      for (let j = i + 1; j < g1Cols.length; j++) {
        const ci = g1Cols[i];
        const cj = g1Cols[j];
        const xs: number[] = [];
        const ys: number[] = [];
        for (let r = 0; r < rowCount; r++) {
          if (!isMissing(ci, r) && !isMissing(cj, r)) {
            xs.push(ci.get(r) as number);
            ys.push(cj.get(r) as number);
          }
        }
        const r = pearson(xs, ys);
        if (!Number.isNaN(r)) {
          sum += r;
          pairCount++;
        }
      }
    }
    controlCorr = pairCount === 0 ? NaN : sum / pairCount;
  }

  return {
    median_intensity: medianIntensity,
    missing_pct: missingPct,
    control_corr: controlCorr,
    protein_count: rowsWithData,
    sample_count: sampleCount,
    computed_at: new Date().toISOString(),
  };
}

// =====================================================================================
// Nelson rules engine (1984 definitions, locked in CONTEXT.md canonical_refs).
// =====================================================================================

export function defaultRulesEnabled(): FlatRulesEnabled {
  return {
    nelson_1: true, nelson_2: false, nelson_3: false, nelson_4: false,
    nelson_5: true, nelson_6: false, nelson_7: false, nelson_8: false,
  };
}

export function defaultRulesEnabledAllMetrics(): SpcRulesEnabled {
  const d = defaultRulesEnabled();
  return {
    median_intensity: {...d}, missing_pct: {...d},
    control_corr: {...d}, protein_count: {...d},
  };
}

interface RuleTripContext {
  series: number[];
  z: number[];
  mean: number;
  sd: number;
}

function rule1(ctx: RuleTripContext): boolean {
  // Single point > 3σ. Evaluated at the tail.
  const i = ctx.z.length - 1;
  return Math.abs(ctx.z[i]) > 3;
}

function rule2(ctx: RuleTripContext): boolean {
  // 9 consecutive points on the same side of CL.
  if (ctx.series.length < 9) return false;
  const tail = ctx.series.slice(-9);
  const above = tail.every((v) => v > ctx.mean);
  const below = tail.every((v) => v < ctx.mean);
  return above || below;
}

function rule3(ctx: RuleTripContext): boolean {
  // 6 consecutive points steadily increasing or decreasing.
  if (ctx.series.length < 6) return false;
  const t = ctx.series.slice(-6);
  let inc = true, dec = true;
  for (let i = 1; i < 6; i++) {
    if (!(t[i] > t[i - 1])) inc = false;
    if (!(t[i] < t[i - 1])) dec = false;
  }
  return inc || dec;
}

function rule4(ctx: RuleTripContext): boolean {
  // 14 consecutive points alternating above/below CL.
  if (ctx.series.length < 14) return false;
  const t = ctx.series.slice(-14);
  const sides = t.map((v) => v > ctx.mean ? 1 : (v < ctx.mean ? -1 : 0));
  for (let i = 1; i < sides.length; i++) {
    if (sides[i] === 0 || sides[i - 1] === 0) return false;
    if (sides[i] === sides[i - 1]) return false;
  }
  return true;
}

function rule5(ctx: RuleTripContext): boolean {
  // 2-of-3 consecutive points beyond 2σ on the same side. Evaluated at the tail.
  if (ctx.z.length < 3) {
    // With fewer than 3 points, count current + prior(s).
    const z = ctx.z;
    if (z.length < 2) return false;
    const above2 = z.filter((v) => v > 2).length;
    const below2 = z.filter((v) => v < -2).length;
    return above2 >= 2 || below2 >= 2;
  }
  const last3 = ctx.z.slice(-3);
  const above = last3.filter((v) => v > 2).length;
  const below = last3.filter((v) => v < -2).length;
  return above >= 2 || below >= 2;
}

function rule6(ctx: RuleTripContext): boolean {
  // 4-of-5 consecutive points beyond 1σ on the same side.
  if (ctx.z.length < 5) return false;
  const last5 = ctx.z.slice(-5);
  const above = last5.filter((v) => v > 1).length;
  const below = last5.filter((v) => v < -1).length;
  return above >= 4 || below >= 4;
}

function rule7(ctx: RuleTripContext): boolean {
  // 15 consecutive points within 1σ of CL.
  if (ctx.z.length < 15) return false;
  return ctx.z.slice(-15).every((v) => Math.abs(v) < 1);
}

function rule8(ctx: RuleTripContext): boolean {
  // 8 consecutive points outside 1σ on either side.
  if (ctx.z.length < 8) return false;
  return ctx.z.slice(-8).every((v) => Math.abs(v) > 1);
}

const RULE_FNS: Record<NelsonRuleId, (ctx: RuleTripContext) => boolean> = {
  nelson_1: rule1, nelson_2: rule2, nelson_3: rule3, nelson_4: rule4,
  nelson_5: rule5, nelson_6: rule6, nelson_7: rule7, nelson_8: rule8,
};

const RULE_IDS: NelsonRuleId[] = [
  'nelson_1', 'nelson_2', 'nelson_3', 'nelson_4',
  'nelson_5', 'nelson_6', 'nelson_7', 'nelson_8',
];

/**
 * Simple form — evaluates the enabled Nelson rules on a single metric's series.
 * `enabled` is a flat record (no per-metric grid). For the aggregate form used by
 * the dashboard layer, see {@link evaluateNelsonRulesAllMetrics}.
 */
export function evaluateNelsonRules(
  series: number[],
  baselineMean: number,
  baselineSd: number,
  enabled: FlatRulesEnabled,
  metricLabel?: string,
): SpcRuleResult {
  if (series.length === 0 || !Number.isFinite(baselineSd) || baselineSd === 0)
    return {status: 'pass', rulesTripped: []};
  const z = series.map((v) => (v - baselineMean) / baselineSd);
  const ctx: RuleTripContext = {series, z, mean: baselineMean, sd: baselineSd};
  const rulesTripped: string[] = [];
  for (const id of RULE_IDS) {
    if (!enabled[id]) continue;
    if (RULE_FNS[id](ctx))
      rulesTripped.push(metricLabel ? `${id}@${metricLabel}` : id);
  }
  return {status: classifyStatus(rulesTripped), rulesTripped};
}

/**
 * Aggregate form — evaluates rules across all four metrics. Used by the dashboard
 * and storage layers; the simple form above is the test-facing entry point.
 */
export function evaluateNelsonRulesAllMetrics(
  newMetrics: SpcMetrics,
  priorSeries: {
    median_intensity: number[]; missing_pct: number[];
    control_corr: number[]; protein_count: number[];
  },
  baseline: BaselineSnapshot,
  enabled: SpcRulesEnabled,
): SpcRuleResult {
  const tripped: string[] = [];
  const metricNames: Array<keyof BaselineSnapshot> = [
    'median_intensity', 'missing_pct', 'control_corr', 'protein_count',
  ];
  for (const m of metricNames) {
    const newValue = newMetrics[m];
    if (!Number.isFinite(newValue)) continue;     // D-03: NaN metric = unavailable, skip
    const series = priorSeries[m].concat([newValue]);
    const r = evaluateNelsonRules(series, baseline[m].mean, baseline[m].sd, enabled[m], m);
    for (const trip of r.rulesTripped) tripped.push(trip);
  }
  return {status: classifyStatus(tripped), rulesTripped: tripped};
}

export function classifyStatus(rulesTripped: string[]): SpcStatus {
  if (rulesTripped.length === 0) return 'pass';
  for (const t of rulesTripped) {
    if (t.startsWith('nelson_1@') || t === 'nelson_1') return 'out_of_control';
  }
  return 'flagged';
}

// =====================================================================================
// Iterative baseline (3σ outlier removal, capped at 2 iterations per Pitfall 6).
// =====================================================================================

function meanSd(values: number[]): {mean: number; sd: number} {
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
  // Sample standard deviation (n-1) when n>1; population (n) when n=1.
  const sd = n > 1 ? Math.sqrt(sse / (n - 1)) : 0;
  return {mean, sd};
}

const MAX_BASELINE_ITERATIONS = 2;

export function computeIterativeBaseline(
  values: number[], iterateOutliers: boolean,
): IterativeBaselineResult {
  const includedIdx = new Set<number>();
  for (let i = 0; i < values.length; i++) includedIdx.add(i);
  const iterationTrace: Array<{iter: number; dropped: number[]}> = [];
  let mean = NaN;
  let sd = NaN;

  if (!iterateOutliers) {
    ({mean, sd} = meanSd(values));
    return {
      mean, sd, includedCount: values.length, excludedIndices: [], iterationTrace,
    };
  }

  for (let iter = 1; iter <= MAX_BASELINE_ITERATIONS; iter++) {
    const current = Array.from(includedIdx).map((i) => values[i]);
    const stats = meanSd(current);
    mean = stats.mean;
    sd = stats.sd;
    if (!Number.isFinite(sd) || sd === 0) break;
    const dropped: number[] = [];
    for (const i of includedIdx) {
      const z = (values[i] - mean) / sd;
      if (Math.abs(z) > 3) dropped.push(i);
    }
    if (dropped.length === 0) break;
    for (const i of dropped) includedIdx.delete(i);
    iterationTrace.push({iter, dropped});
    // Recompute after dropping for the final returned mean/sd.
    const updated = Array.from(includedIdx).map((i) => values[i]);
    ({mean, sd} = meanSd(updated));
  }

  const excluded: number[] = [];
  for (let i = 0; i < values.length; i++)
    if (!includedIdx.has(i)) excluded.push(i);

  return {
    mean, sd,
    includedCount: includedIdx.size,
    excludedIndices: excluded,
    iterationTrace,
  };
}

// =====================================================================================
// Status persistence on the analyzed df (tags + belt-and-braces column).
// =====================================================================================

function ensureFreshStringColumn(df: DG.DataFrame, name: string): DG.Column {
  if (df.columns.contains(name))
    df.columns.remove(name);
  return df.columns.addNewString(name);
}

export function setSpcStatus(
  df: DG.DataFrame, metrics: SpcMetrics, ruleResult: SpcRuleResult,
): void {
  df.setTag('proteomics.spc_metrics', JSON.stringify(metrics));
  df.setTag('proteomics.spc_status', ruleResult.status);
  df.setTag('proteomics.spc_rules_tripped', JSON.stringify(ruleResult.rulesTripped));

  // Belt-and-braces (Phase 15 D-05) — single column carries the JSON so a
  // published Project's serializer-strip can be recovered from this column
  // even if df tags don't survive.
  const col = ensureFreshStringColumn(df, BELT_AND_BRACES_COLUMN);
  const payload = JSON.stringify({
    metrics, status: ruleResult.status, rulesTripped: ruleResult.rulesTripped,
  });
  if (df.rowCount === 0) {
    // Empty df has no row to write — the column still pins the schema.
    return;
  }
  col.set(0, payload);
}

export function getSpcStatus(df: DG.DataFrame): {
  metrics: SpcMetrics; status: SpcStatus; rulesTripped: string[];
} | null {
  const metricsRaw = df.getTag('proteomics.spc_metrics');
  const statusRaw = df.getTag('proteomics.spc_status');
  const rulesRaw = df.getTag('proteomics.spc_rules_tripped');
  if (!metricsRaw || !statusRaw) return null;
  try {
    return {
      metrics: JSON.parse(metricsRaw) as SpcMetrics,
      status: statusRaw as SpcStatus,
      rulesTripped: rulesRaw ? (JSON.parse(rulesRaw) as string[]) : [],
    };
  } catch {
    return null;
  }
}

/**
 * Recovers metrics from the belt-and-braces column when tags are missing
 * (Phase 15 D-05 fallback — protects against the serializer-strip pattern).
 */
export function recoverSpcMetricsFromColumn(df: DG.DataFrame): SpcMetrics | null {
  if (!df.columns.contains(BELT_AND_BRACES_COLUMN)) return null;
  const col = df.col(BELT_AND_BRACES_COLUMN);
  if (!col || df.rowCount === 0) return null;
  const raw = col.get(0);
  if (typeof raw !== 'string' || raw.length === 0) return null;
  try {
    const parsed = JSON.parse(raw);
    return parsed && parsed.metrics ? (parsed.metrics as SpcMetrics) : null;
  } catch {
    return null;
  }
}

import * as DG from 'datagrok-api/dg';

import {findColumn} from './column-detection';
import {getGroups} from '../analysis/experiment-setup';

/** Per-condition abundance, in log10(MS intensity) space, one value per protein
 * row (null when the protein was not quantified in that condition). */
export interface ConditionAbundance {
  name: string;
  log10Intensity: (number | null)[];
}

/** Abundance for the two experiment conditions — the input the rank-abundance
 * (dynamic-range) plot ranks and plots. */
export interface AbundanceByCondition {
  group1: ConditionAbundance;
  group2: ConditionAbundance;
}

/** Spectronaut Candidates carries pre-computed per-group mean abundance in these
 * two columns (kept verbatim by the Candidates parser). group1 = Numerator. */
const QTY_NUM_HINTS = ['avg group quantity numerator'];
const QTY_DEN_HINTS = ['avg group quantity denominator'];

/** Below this max value a column reads as log2-scale (log2 MS intensities rarely
 * exceed ~40; linear intensities are >>1000). Mirrors the log2-status heuristic
 * in `parsers/shared-utils`. */
const LOG2_MAX_THRESHOLD = 45;

/** log10 of a strictly-positive value, else null. */
function log10OrNull(v: number): number | null {
  return Number.isFinite(v) && v > 0 ? Math.log10(v) : null;
}

/** Single linear-abundance column → per-row log10. */
function columnToLog10(col: DG.Column): (number | null)[] {
  const out: (number | null)[] = new Array(col.length);
  for (let i = 0; i < col.length; i++)
    out[i] = col.isNone(i) ? null : log10OrNull(col.get(i) as number);
  return out;
}

/** How a group's per-sample intensity columns are collapsed to one per-protein
 * abundance on the Report path:
 * - `arithmetic` — delinearize (2^v) log2 columns, average in linear space, then
 *   log10. A true intensity mean; outlier-sensitive. Default (rank-abundance).
 * - `geometric` — average in log space (mean of log10 values). log10 of the
 *   geometric mean; the conventional summary for log-normal MS intensities and
 *   the one the abundance-correlation viewer uses. */
export type ReportMean = 'arithmetic' | 'geometric';

/** Per-row abundance (log10 space) across a group's intensity columns, collapsed
 * per `mode`. When the columns are log2-scale (Report analysis columns are
 * log2-transformed) each value is treated accordingly so the axis is consistent
 * with the Candidates path regardless of the source scale. */
function groupMeanToLog10(
  df: DG.DataFrame, colNames: string[], mode: ReportMean,
): (number | null)[] {
  const cols = colNames.map((n) => df.col(n)).filter((c): c is DG.Column => c != null);
  const n = df.rowCount;
  const out: (number | null)[] = new Array(n).fill(null);
  if (cols.length === 0) return out;

  // Detect log2 scale from the observed maximum across the group's columns.
  let maxVal = -Infinity;
  for (const c of cols) {
    for (let i = 0; i < n; i++) {
      if (!c.isNone(i)) {
        const v = c.get(i) as number;
        if (Number.isFinite(v) && v > maxVal) maxVal = v;
      }
    }
  }
  const isLog2 = Number.isFinite(maxVal) && maxVal < LOG2_MAX_THRESHOLD;
  const LOG10_2 = Math.log10(2);

  for (let i = 0; i < n; i++) {
    let sum = 0; let count = 0;
    for (const c of cols) {
      if (c.isNone(i)) continue;
      const v = c.get(i) as number;
      if (!Number.isFinite(v)) continue;
      if (mode === 'geometric') {
        // Mean of log10 values = log10(geometric mean). A log2 column is already
        // a log, so rebase (v * log10 2); a linear column is log10'd (skip <=0).
        const lv = isLog2 ? v * LOG10_2 : log10OrNull(v);
        if (lv == null) continue;
        sum += lv;
      } else {
        // Arithmetic mean of linear intensities, log10'd after the loop.
        sum += isLog2 ? Math.pow(2, v) : v;
      }
      count++;
    }
    if (count === 0) { out[i] = null; continue; }
    out[i] = mode === 'geometric' ? sum / count : log10OrNull(sum / count);
  }
  return out;
}

/**
 * Resolves per-condition abundance for the rank-abundance / dynamic-range plot,
 * working on both pipeline shapes:
 *
 * - **Candidates** — the pre-computed `AVG Group Quantity Numerator/Denominator`
 *   columns (linear); group1 = Numerator, group2 = Denominator. This is the
 *   shape CK-omics builds this plot from.
 * - **Report** — the per-group mean of each group's intensity columns (from the
 *   `proteomics.groups` assignment), delinearized when log2-scale.
 *
 * Returns null when neither shape yields abundance (e.g. a Report with no group
 * assignment yet, or an older Candidates export lacking group quantities) — the
 * caller warns and no-ops.
 *
 * `opts.reportMean` selects how the Report path collapses each group's per-sample
 * columns (default `arithmetic`, so rank-abundance is unchanged; the
 * abundance-correlation viewer passes `geometric`). It does not affect the
 * Candidates path, which reads the vendor's pre-computed per-group quantity.
 */
export function findAbundanceByCondition(
  df: DG.DataFrame, opts?: {reportMean?: ReportMean},
): AbundanceByCondition | null {
  const groups = getGroups(df);
  const mode: ReportMean = opts?.reportMean ?? 'arithmetic';

  // Candidates path — precomputed per-group means.
  const numCol = findColumn(df, '', QTY_NUM_HINTS);
  const denCol = findColumn(df, '', QTY_DEN_HINTS);
  if (numCol && denCol) {
    return {
      group1: {name: groups?.group1.name || 'Numerator', log10Intensity: columnToLog10(numCol)},
      group2: {name: groups?.group2.name || 'Denominator', log10Intensity: columnToLog10(denCol)},
    };
  }

  // Report path — per-group summary of each group's intensity columns.
  if (groups && groups.group1.columns.length > 0 && groups.group2.columns.length > 0) {
    return {
      group1: {name: groups.group1.name || 'Group 1', log10Intensity: groupMeanToLog10(df, groups.group1.columns, mode)},
      group2: {name: groups.group2.name || 'Group 2', log10Intensity: groupMeanToLog10(df, groups.group2.columns, mode)},
    };
  }

  return null;
}

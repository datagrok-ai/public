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

/** Per-row log10 of the linear mean across a group's intensity columns. When the
 * columns are log2-scale (Report analysis columns are log2-transformed), each
 * value is delinearized (2^v) before averaging so the mean is a true intensity
 * mean, then log10'd for a consistent axis with the Candidates path. */
function groupMeanToLog10(df: DG.DataFrame, colNames: string[]): (number | null)[] {
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

  for (let i = 0; i < n; i++) {
    let sum = 0; let count = 0;
    for (const c of cols) {
      if (c.isNone(i)) continue;
      const v = c.get(i) as number;
      if (!Number.isFinite(v)) continue;
      sum += isLog2 ? Math.pow(2, v) : v;
      count++;
    }
    out[i] = count > 0 ? log10OrNull(sum / count) : null;
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
 */
export function findAbundanceByCondition(df: DG.DataFrame): AbundanceByCondition | null {
  const groups = getGroups(df);

  // Candidates path — precomputed per-group means.
  const numCol = findColumn(df, '', QTY_NUM_HINTS);
  const denCol = findColumn(df, '', QTY_DEN_HINTS);
  if (numCol && denCol) {
    return {
      group1: {name: groups?.group1.name || 'Numerator', log10Intensity: columnToLog10(numCol)},
      group2: {name: groups?.group2.name || 'Denominator', log10Intensity: columnToLog10(denCol)},
    };
  }

  // Report path — mean of each group's intensity columns.
  if (groups && groups.group1.columns.length > 0 && groups.group2.columns.length > 0) {
    return {
      group1: {name: groups.group1.name || 'Group 1', log10Intensity: groupMeanToLog10(df, groups.group1.columns)},
      group2: {name: groups.group2.name || 'Group 2', log10Intensity: groupMeanToLog10(df, groups.group2.columns)},
    };
  }

  return null;
}

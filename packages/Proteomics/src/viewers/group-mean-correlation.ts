import * as DG from 'datagrok-api/dg';
import {findColumn} from '../utils/column-detection';
import {SEMTYPE} from '../utils/proteomics-types';
import {getGroups, GroupAssignment} from '../analysis/experiment-setup';

/**
 * R4 / D-12 — Group-Mean Correlation viewer.
 *
 * Adds two derived columns to the host DataFrame (Numerator Mean, Denominator
 * Mean — means of the group1/group2 intensity columns per row) and renders
 * them as a scatter colored by the existing `direction` column (magenta/cyan/
 * gray per D-04). The viewer title carries inline Pearson r and Spearman ρ
 * over the current filter; a y=x diagonal reference line is drawn in #888888.
 *
 * The derived columns are re-run-safe via `ensureFreshFloat` — a second
 * invocation replaces the columns instead of duplicating them.
 *
 * Statistics: `@datagrok-libraries/statistics` exports only Kendall's tau,
 * so Pearson and Spearman are implemented inline (Pearson over raw values,
 * Spearman = Pearson over fractional ranks).
 */

const NUMERATOR_MEAN_COL = 'Numerator Mean';
const DENOMINATOR_MEAN_COL = 'Denominator Mean';
const DIRECTION_COL = 'direction';
const DIAGONAL_FORMULA = `\${${NUMERATOR_MEAN_COL}} = \${${DENOMINATOR_MEAN_COL}}`;
const DIAGONAL_COLOR = '#888888';

/** Bulk-init pattern — removes an existing column first to avoid duplicates on
 * re-run. Mirrors src/viewers/qc-computations.ts:28-32. */
function ensureFreshFloat(df: DG.DataFrame, name: string): DG.Column {
  if (df.columns.contains(name))
    df.columns.remove(name);
  return df.columns.addNewFloat(name);
}

/** Mean of non-null values across `cols` at `rowIdx`. Returns NaN if every
 * value is null. Mirrors qc-computations.ts:15-25. */
function groupMean(cols: DG.Column[], rowIdx: number): number {
  let sum = 0;
  let n = 0;
  for (const c of cols) {
    if (c.isNone(rowIdx)) continue;
    sum += c.get(rowIdx) as number;
    n++;
  }
  return n > 0 ? sum / n : NaN;
}

/**
 * Derives Numerator Mean (group1) and Denominator Mean (group2) columns on
 * `df`. The join key is the protein row index — never a derived key like
 * gene name (project memory project_proteomics_deqms_result_misalignment).
 *
 * The columns carry the dedicated SEMTYPE.NUMERATOR_MEAN / SEMTYPE.DENOMINATOR_MEAN
 * tags so downstream tooling can find them via `findColumn`.
 */
export function computeGroupMeans(df: DG.DataFrame, groups: GroupAssignment): void {
  const g1Cols = groups.group1.columns
    .map((n) => df.col(n))
    .filter((c): c is DG.Column => c != null);
  const g2Cols = groups.group2.columns
    .map((n) => df.col(n))
    .filter((c): c is DG.Column => c != null);

  const nRows = df.rowCount;
  const numArr = new Float32Array(nRows);
  const denArr = new Float32Array(nRows);
  for (let i = 0; i < nRows; i++) {
    const m1 = groupMean(g1Cols, i);
    const m2 = groupMean(g2Cols, i);
    numArr[i] = isNaN(m1) ? DG.FLOAT_NULL : m1;
    denArr[i] = isNaN(m2) ? DG.FLOAT_NULL : m2;
  }

  const numCol = ensureFreshFloat(df, NUMERATOR_MEAN_COL);
  numCol.init((i) => numArr[i]);
  numCol.semType = SEMTYPE.NUMERATOR_MEAN;

  const denCol = ensureFreshFloat(df, DENOMINATOR_MEAN_COL);
  denCol.init((i) => denArr[i]);
  denCol.semType = SEMTYPE.DENOMINATOR_MEAN;
}

/** Pearson correlation. Skips paired entries containing NaN or FLOAT_NULL. */
export function pearson(xs: number[], ys: number[]): number {
  let n = 0;
  let sumX = 0;
  let sumY = 0;
  for (let i = 0; i < xs.length; i++) {
    const x = xs[i];
    const y = ys[i];
    if (!Number.isFinite(x) || !Number.isFinite(y) || x === DG.FLOAT_NULL || y === DG.FLOAT_NULL) continue;
    sumX += x;
    sumY += y;
    n++;
  }
  if (n === 0) return NaN;
  const meanX = sumX / n;
  const meanY = sumY / n;
  let num = 0;
  let denX = 0;
  let denY = 0;
  for (let i = 0; i < xs.length; i++) {
    const x = xs[i];
    const y = ys[i];
    if (!Number.isFinite(x) || !Number.isFinite(y) || x === DG.FLOAT_NULL || y === DG.FLOAT_NULL) continue;
    const dx = x - meanX;
    const dy = y - meanY;
    num += dx * dy;
    denX += dx * dx;
    denY += dy * dy;
  }
  const den = Math.sqrt(denX * denY);
  return den === 0 ? NaN : num / den;
}

/** Fractional ranks (1-based). Ties get the average of the positions they span. */
function rank(arr: number[]): number[] {
  const indexed = arr.map((v, i) => ({v, i}));
  indexed.sort((a, b) => a.v - b.v);
  const ranks = new Array<number>(arr.length);
  let i = 0;
  while (i < indexed.length) {
    let j = i;
    while (j + 1 < indexed.length && indexed[j + 1].v === indexed[i].v) j++;
    const avg = (i + j) / 2 + 1; // 1-based average of positions [i..j]
    for (let k = i; k <= j; k++) ranks[indexed[k].i] = avg;
    i = j + 1;
  }
  return ranks;
}

/** Spearman = Pearson over fractional ranks. */
export function spearman(xs: number[], ys: number[]): number {
  return pearson(rank(xs), rank(ys));
}

/** Computes Pearson + Spearman over the filtered subset of two columns. */
function computePearsonSpearman(
  df: DG.DataFrame, xName: string, yName: string,
): {r: number; rho: number} {
  const xCol = df.col(xName);
  const yCol = df.col(yName);
  if (!xCol || !yCol) return {r: NaN, rho: NaN};
  const xRaw = xCol.getRawData() as Float32Array | Float64Array;
  const yRaw = yCol.getRawData() as Float32Array | Float64Array;

  const xs: number[] = [];
  const ys: number[] = [];
  for (let i = 0; i < df.rowCount; i++) {
    if (!df.filter.get(i)) continue;
    const x = xRaw[i];
    const y = yRaw[i];
    if (x === DG.FLOAT_NULL || y === DG.FLOAT_NULL) continue;
    if (!Number.isFinite(x) || !Number.isFinite(y)) continue;
    xs.push(x);
    ys.push(y);
  }
  return {r: pearson(xs, ys), rho: spearman(xs, ys)};
}

/**
 * Main factory — derives Numerator/Denominator Mean, plots them as a scatter
 * colored by direction, draws the y=x diagonal, and writes an inline
 * Pearson/Spearman annotation into the title.
 */
export function createGroupMeanCorrelation(
  df: DG.DataFrame,
  options?: {title?: string},
): DG.ScatterPlotViewer {
  const groups = getGroups(df);
  if (!groups)
    throw new Error('Annotate experimental groups first (Proteomics | Annotate Experiment)');

  computeGroupMeans(df, groups);

  const sp = df.plot.scatter({
    x: NUMERATOR_MEAN_COL,
    y: DENOMINATOR_MEAN_COL,
    color: DIRECTION_COL,
  });

  // Display Name primary, Gene name fallback (Plan 02 pattern).
  const labelCol = findColumn(df, SEMTYPE.DISPLAY_NAME, ['display name']) ??
    findColumn(df, SEMTYPE.GENE_SYMBOL, ['gene name', 'gene symbol']);
  if (labelCol) {
    sp.props.labelColumnNames = [labelCol.name];
    sp.props.showLabelsFor = 'MouseOverRow';
  }

  // y = x diagonal reference — replace any prior diagonal so re-run doesn't stack.
  df.meta.formulaLines.items = df.meta.formulaLines.items.filter((line) => {
    const f = (line.formula ?? '') as string;
    return !(typeof f === 'string' &&
      f.includes(NUMERATOR_MEAN_COL) && f.includes(DENOMINATOR_MEAN_COL));
  });
  df.meta.formulaLines.addLine({
    formula: DIAGONAL_FORMULA,
    color: DIAGONAL_COLOR,
    width: 1,
    visible: true,
  });
  sp.props.showViewerFormulaLines = true;

  // Inline correlation annotation over the current filter.
  const {r, rho} = computePearsonSpearman(df, NUMERATOR_MEAN_COL, DENOMINATOR_MEAN_COL);
  const baseTitle = options?.title ?? 'Group-Mean Correlation';
  const rStr = Number.isFinite(r) ? r.toFixed(2) : 'NA';
  const rhoStr = Number.isFinite(rho) ? rho.toFixed(2) : 'NA';
  sp.setOptions({
    title: `${baseTitle} — r=${rStr} (Pearson), ρ=${rhoStr} (Spearman)`,
    xColumnLabel: `${groups.group1.name} mean (log2 intensity)`,
    yColumnLabel: `${groups.group2.name} mean (log2 intensity)`,
  } as any);

  return sp;
}

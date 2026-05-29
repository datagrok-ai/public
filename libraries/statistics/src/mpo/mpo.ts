/* eslint-disable guard-for-in */
/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import {TemplateFunction} from '../compute-functions/types';
import {AggregationCode, AGG_CODE, ColumnDesirability, HoistedColumn} from './mpo-types';

/// An array of [x, y] points representing the desirability line
/// [x, y] pairs are sorted by x in ascending order
export type DesirabilityLine = number[][];

export type DesirabilityMode = 'freeform' | 'gaussian' | 'sigmoid';

export type MissingValueConfig =
  | { strategy: 'exclude' }
  | { strategy: 'default'; score: number }
  | { strategy: 'skip' };

type BasePropertyDesirability = {
  weight: number; /// 0-1
  weightColumn?: string; /// optional column name to read per-row weights from (0-1)
  missingValues?: MissingValueConfig;
  computeFunction?: TemplateFunction;
}

export type NumericalDesirability = BasePropertyDesirability & {
  functionType: 'numerical';
  line: DesirabilityLine;
  min?: number; /// min value of the property (optional; used for editing the line)
  max?: number; /// max value of the property (optional; used for editing the line)

  /// When 'log', the x-domain (line x-coords, min/max, mean/x0) is stored in log10 units, so the curve is
  /// designed and evaluated against a log10 x-axis; non-positive raw values count as missing. Absent → linear.
  xScale?: 'linear' | 'log';

  mode?: DesirabilityMode;

  /// Gaussian mode parameters
  mean?: number;
  sigma?: number;

  /// Sigmoid mode parameters
  x0?: number;
  k?: number;

  freeformLine?: DesirabilityLine;
}

export type CategoricalDesirability = BasePropertyDesirability & {
  functionType: 'categorical';
  categories: { name: string; desirability: number }[];
}

export type PropertyDesirability = NumericalDesirability | CategoricalDesirability;

export function isNumerical(p: PropertyDesirability): p is NumericalDesirability {
  return p.functionType === 'numerical';
}

export function createDefaultNumerical(weight = 1, min = 0, max = 1): NumericalDesirability {
  return {functionType: 'numerical', weight, mode: 'freeform', min, max, line: []};
}

export function createDefaultCategorical(
  weight = 1,
  categories?: {name: string; desirability: number}[],
): CategoricalDesirability {
  return {functionType: 'categorical', weight, categories: categories ?? [{name: 'Category 1', desirability: 1}]};
}

export function migrateDesirability(raw: any): PropertyDesirability {
  if (raw.functionType)
    return raw;
  if (raw.categories)
    return {...raw, functionType: 'categorical'};
  return {...raw, functionType: 'numerical'};
}

export const DESIRABILITY_PROFILE_TYPE = 'MPO Desirability Profile';
export const CURRENT_MPO_VERSION = 1;

/// A map of desirability lines with their weights
export type DesirabilityProfile = {
  type: typeof DESIRABILITY_PROFILE_TYPE;
  version?: number;
  name: string;
  description: string;
  aggregation?: WeightedAggregation;
  properties: { [key: string]: PropertyDesirability };
}

export function isDesirabilityProfile(x: any): x is DesirabilityProfile {
  return x != null && typeof x === 'object' && x.type === DESIRABILITY_PROFILE_TYPE;
}

export function migrateProfile(raw: DesirabilityProfile): DesirabilityProfile {
  const version = raw.version ?? 0;
  if (version >= CURRENT_MPO_VERSION)
    return raw;

  // v0 → v1: backfill functionType on properties
  if (version < 1) {
    for (const key in raw.properties)
      raw.properties[key] = migrateDesirability(raw.properties[key]);
    raw.version = CURRENT_MPO_VERSION;
  }

  return raw;
}

export const WEIGHTED_AGGREGATIONS = ['Average', 'Sum', 'Product', 'Geomean', 'Min', 'Max'] as const;
export const WEIGHTED_AGGREGATIONS_LIST: WeightedAggregation[] = [...WEIGHTED_AGGREGATIONS];
export type WeightedAggregation = typeof WEIGHTED_AGGREGATIONS[number];
export const DEFAULT_AGGREGATION: WeightedAggregation = 'Average';

/// Calculates the desirability score for a given x value
/// Returns 0 if x is outside the range of the desirability line
/// Otherwise, returns the y value of the desirability line at x
export function desirabilityScore(x: number, desirabilityLine: DesirabilityLine): number {
  // If the line is empty or x is outside the range, return 0
  if (desirabilityLine.length === 0 || x < desirabilityLine[0][0] || x > desirabilityLine[desirabilityLine.length - 1][0])
    return 0;

  // Find the two points that x lies between
  for (let i = 0; i < desirabilityLine.length - 1; i++) {
    const [x1, y1] = desirabilityLine[i];
    const [x2, y2] = desirabilityLine[i + 1];

    if (x >= x1 && x <= x2) {
      // Linear interpolation between the two points
      if (x1 === x2) return y1;
      const slope = (y2 - y1) / (x2 - x1);
      return y1 + slope * (x - x1);
    }
  }

  // Should not happen if x is within bounds, but return 0 as fallback
  return 0;
}

export function categoricalDesirabilityScore(
  value: string,
  prop: CategoricalDesirability,
): number | null {
  const found = prop.categories.find((c) => c.name === value);
  return found?.desirability ?? null;
}

export type MpoResult = {
  scoreColumn: DG.Column;
  desirabilityColumns?: DG.Column[];
}

/// Fingerprint of the per-cell desirability mapping for a property — the line/categories and missing-value
/// handling, but NOT weight/weightColumn (those feed the reduction, not the mapping). Used as the
/// MpoCalculator cache key so editing one property's curve only recomputes that one column.
export function desirabilityKey(d: PropertyDesirability): string {
  return isNumerical(d) ?
    `n|${d.xScale ?? 'linear'}|${JSON.stringify(d.line)}|${JSON.stringify(d.missingValues ?? null)}` :
    `c|${JSON.stringify((d as CategoricalDesirability).categories)}|${JSON.stringify(d.missingValues ?? null)}`;
}

/// Hoists all column access out of the per-row loops: one getRawData() per column instead of per-cell
/// get()/isNone() (which cross the JS↔Dart boundary). The template is read from the 'desirabilityTemplate' tag.
export function hoistColumns(dataFrame: DG.DataFrame, columns: DG.Column[]): HoistedColumn[] {
  const hoisted: HoistedColumn[] = [];
  for (const col of columns) {
    const template = migrateDesirability(JSON.parse(col.getTag('desirabilityTemplate')));
    const isCat = !isNumerical(template);
    const wc = template.weightColumn ? dataFrame.col(template.weightColumn) : null;
    hoisted.push({
      col,
      template,
      raw: col.getRawData(),
      valNull: col.type === DG.COLUMN_TYPE.INT ? DG.INT_NULL : DG.FLOAT_NULL,
      isCat,
      cats: isCat ? col.categories : null,
      catScore: isCat ?
        new Map((template as CategoricalDesirability).categories.map((c) => [c.name, c.desirability])) : null,
      staticW: template.weight,
      rawW: wc ? wc.getRawData() : null,
      wNull: wc && wc.type === DG.COLUMN_TYPE.INT ? DG.INT_NULL : DG.FLOAT_NULL,
    });
  }
  return hoisted;
}

/// The expensive, cacheable phase: maps one column's values to desirability scores. Column-major sweep of
/// the contiguous raw array, with desirabilityScore inlined (no per-cell function call). Depends only on the
/// column data + line + missing-value handling — never on weight/aggregation — so its result is cacheable.
export function mapColumnDesirability(h: HoistedColumn, rowCount: number): ColumnDesirability {
  const {raw, valNull, isCat, cats, catScore, col} = h;
  const line = isCat ? null : (h.template as NumericalDesirability).line;
  const segN = line ? line.length : 0;
  const loX = segN ? line![0][0] : 0;
  const hiX = segN ? line![segN - 1][0] : 0;
  const logScale = !isCat && (h.template as NumericalDesirability).xScale === 'log';
  const mv = h.template.missingValues;
  const excludeMissing = !mv || mv.strategy === 'exclude';
  const skipMissing = mv?.strategy === 'skip';
  const defaultScore = mv && mv.strategy === 'default' ? mv.score : 0;

  const D = new Float32Array(rowCount);
  let state: Uint8Array | null = null; // lazily allocated on the first skip/bail row

  for (let i = 0; i < rowCount; i++) {
    const isMissing = isCat ? col.isNone(i) : (raw[i] === valNull || (logScale && raw[i] <= 0));
    if (isMissing) {
      if (excludeMissing) {
        (state ??= new Uint8Array(rowCount))[i] = 2;
        D[i] = NaN;
      }
      else if (skipMissing) {
        (state ??= new Uint8Array(rowCount))[i] = 1;
        D[i] = NaN;
      }
      else
        D[i] = defaultScore;
      continue;
    }
    if (isCat) {
      const s = catScore!.get(cats![raw[i]]);
      if (s == null) {
        (state ??= new Uint8Array(rowCount))[i] = 2;
        D[i] = NaN;
      }
      else
        D[i] = s;
      continue;
    }
    // Inlined desirabilityScore over the column's line.
    const x = logScale ? Math.log10(raw[i]) : raw[i];
    let score = 0;
    if (segN !== 0 && x >= loX && x <= hiX) {
      for (let k = 0; k < segN - 1; k++) {
        const x1 = line![k][0];
        const x2 = line![k + 1][0];
        if (x >= x1 && x <= x2) {
          const y1 = line![k][1];
          score = x1 === x2 ? y1 : y1 + (line![k + 1][1] - y1) * (x - x1) / (x2 - x1);
          break;
        }
      }
    }
    D[i] = score;
  }
  return {D, state};
}

/// The cheap phase (always re-run; depends on weights + aggregation): folds the per-column desirabilities
/// into the final score, column-major into per-row accumulators.
export function reduceMpo(
  maps: ColumnDesirability[], hoisted: HoistedColumn[], aggregation: WeightedAggregation, rowCount: number,
): Float64Array {
  const aggCode = AGG_CODE[aggregation];
  if (aggCode === undefined)
    throw new Error(`Unknown aggregation type: ${aggregation}`);

  const needWeightSum = aggCode === AggregationCode.Average || aggCode === AggregationCode.Geomean;
  const out = new Float64Array(rowCount); // final per-row score
  const bail = new Uint8Array(rowCount); // 1 → 'exclude' missing or unmatched category → null score
  const cnt = new Int32Array(rowCount); // contributing columns (0 → no data → null score)
  const acc1 = new Float64Array(rowCount); // primary accumulator
  const acc2 = needWeightSum ? new Float64Array(rowCount) : null; // Σweight (Average) / total weight (Geomean)

  if (aggCode === AggregationCode.Product || aggCode === AggregationCode.Geomean)
    acc1.fill(1);
  else if (aggCode === AggregationCode.Min)
    acc1.fill(Infinity);
  else if (aggCode === AggregationCode.Max)
    acc1.fill(-Infinity);

  for (let j = 0; j < maps.length; j++) {
    const D = maps[j].D;
    const st = maps[j].state;
    const rw = hoisted[j].rawW;
    const wN = hoisted[j].wNull;
    const sw = hoisted[j].staticW;

    for (let i = 0; i < rowCount; i++) {
      if (bail[i])
        continue;
      if (st !== null) {
        const s = st[i];
        if (s === 2) {
          bail[i] = 1;
          continue;
        }
        if (s === 1)
          continue;
      }
      const score = D[i];
      const w = rw && rw[i] !== wN ? (rw[i] < 0 ? 0 : rw[i] > 1 ? 1 : rw[i]) : sw;
      cnt[i]++;
      switch (aggCode) {
      case AggregationCode.Sum:
        acc1[i] += score * w;
        break;
      case AggregationCode.Average:
        acc1[i] += score * w;
        acc2![i] += w;
        break;
      case AggregationCode.Product:
        acc1[i] *= Math.pow(score, w);
        break;
      case AggregationCode.Geomean:
        acc1[i] *= Math.pow(score, w);
        acc2![i] += w;
        break;
      case AggregationCode.Min:
        acc1[i] = Math.min(acc1[i], Math.pow(score, w));
        break;
      case AggregationCode.Max:
        acc1[i] = Math.max(acc1[i], Math.pow(score, w));
        break;
      }
    }
  }

  for (let i = 0; i < rowCount; i++) {
    if (bail[i] || cnt[i] === 0)
      out[i] = DG.FLOAT_NULL;
    else if (aggCode === AggregationCode.Average)
      out[i] = acc1[i] / acc2![i];
    else if (aggCode === AggregationCode.Geomean)
      out[i] = Math.pow(acc1[i], 1 / acc2![i]); // (Π sᵏ^wᵏ)^(1/Σw) ≡ Π sᵏ^(wᵏ/Σw)
    else
      out[i] = acc1[i];
  }
  return out;
}

/// Builds the per-property desirability output columns from the mapping results (state 1/2 rows → null).
export function buildDesirabilityColumns(
  columns: DG.Column[], maps: ColumnDesirability[], profileName: string, rowCount: number,
): DG.Column[] {
  const result: DG.Column[] = [];
  for (let j = 0; j < columns.length; j++) {
    const {D, state} = maps[j];
    const dvals = new Float32Array(rowCount);
    for (let i = 0; i < rowCount; i++)
      dvals[i] = state && state[i] !== 0 ? DG.FLOAT_NULL : D[i];
    const c = DG.Column.float(`${columns[j].name} (${profileName} desirability)`, rowCount);
    c.setRawData(dvals, false);
    result.push(c);
  }
  return result;
}

interface CacheEntry {
  lineKey: string;
  D: Float32Array;
  state: Uint8Array | null;
}

/// Stateful MPO engine for the interactive preview path. Same compute() contract as the stateless mpo(),
/// but caches each column's desirability mapping — the expensive phase — keyed by df.id|column and validated
/// by the line/missing-value fingerprint. Editing one property's curve therefore recomputes only that column;
/// the others are reused and only the cheap reduction re-runs. Hold one per editing session and call release()
/// on teardown (dialog close / dataframe switch) to drop the cached arrays.
export class MpoCalculator {
  private cache = new Map<string, CacheEntry>();

  compute(
    dataFrame: DG.DataFrame,
    columns: DG.Column[],
    profileName: string,
    aggregation: WeightedAggregation,
    isDifferent: boolean = false,
    createDesirabilityColumns: boolean = false,
  ): MpoResult {
    if (columns.length === 0)
      throw new Error('No columns provided for MPO calculation.');

    const rowCount = columns[0].length;
    const resultColumn = isDifferent ?
      DG.Column.float(profileName, rowCount) :
      (dataFrame.col(profileName) ?? DG.Column.float(profileName, rowCount));

    const hoisted = hoistColumns(dataFrame, columns);
    const maps: ColumnDesirability[] = [];
    for (const h of hoisted) {
      const lineKey = desirabilityKey(h.template);
      const key = `${dataFrame.id}|${h.col.name}`;
      const hit = this.cache.get(key);
      if (hit && hit.lineKey === lineKey && hit.D.length === rowCount) {
        maps.push({D: hit.D, state: hit.state});
        continue;
      }
      const m = mapColumnDesirability(h, rowCount);
      this.cache.set(key, {lineKey, D: m.D, state: m.state});
      maps.push(m);
    }

    resultColumn.setRawData(reduceMpo(maps, hoisted, aggregation, rowCount), false);
    resultColumn.fireValuesChanged();

    const desirabilityColumns = createDesirabilityColumns ?
      buildDesirabilityColumns(columns, maps, profileName, rowCount) : undefined;
    if (desirabilityColumns) {
      for (const c of desirabilityColumns)
        c.fireValuesChanged();
    }

    return {scoreColumn: resultColumn, desirabilityColumns};
  }

  /// Drops all cached desirability arrays. Call on dialog close / dataframe switch.
  release(): void {
    this.cache = new Map();
  }
}

/** Calculates the multi parameter optimization score, 0-100, 100 is the maximum */
export function mpo(
  dataFrame: DG.DataFrame,
  columns: DG.Column[],
  profileName: string,
  aggregation: WeightedAggregation,
  isDifferent: boolean = false,
  createDesirabilityColumns: boolean = false,
): MpoResult {
  return new MpoCalculator()
    .compute(dataFrame, columns, profileName, aggregation, isDifferent, createDesirabilityColumns);
}

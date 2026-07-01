/* eslint-disable guard-for-in */
/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import {
  AggregationCode, AGG_CODE, CacheEntry, CategoricalDesirability, ColumnDesirability, CURRENT_MPO_VERSION,
  DESIRABILITY_PROFILE_TYPE, DesirabilityProfile, HoistedColumn, MpoResult,
  NumericalDesirability, PropertyDesirability, RowState, WeightedAggregation,
} from './mpo-types';

// mpo-types is the types/constants barrel for this module; re-export it so consumers keep importing from './mpo'.
export * from './mpo-types';

export function isNumerical(p: PropertyDesirability): p is NumericalDesirability {
  return p.functionType === 'numerical';
}

export function createDefaultNumerical(weight = 1, min = 0, max = 1): NumericalDesirability {
  return {functionType: 'numerical', weight, mode: 'freeform', min, max, line: []};
}

export const MPO_NUMERIC_TYPES = new Set<string>([DG.COLUMN_TYPE.INT, DG.COLUMN_TYPE.FLOAT]);

export function isMpoNumericColumn(col: DG.Column): boolean {
  return MPO_NUMERIC_TYPES.has(col.type);
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

/// Fingerprint of a property's value→desirability mapping (line/categories + missing-value handling, but NOT
/// weight). Validates the MpoCalculator cache so editing one curve only recomputes that column.
export function desirabilityKey(d: PropertyDesirability): string {
  return isNumerical(d) ?
    `n|${JSON.stringify(d.line)}|${JSON.stringify(d.missingValues ?? null)}|${d.inverted ? 1 : 0}` :
    `c|${JSON.stringify((d as CategoricalDesirability).categories)}|${JSON.stringify(d.missingValues ?? null)}`;
}

/// Hoists all column access out of the per-row loops: one getRawData() per column instead of per-cell
/// get()/isNone() (which cross the JS↔Dart boundary). The template is read from the 'desirabilityTemplate' tag.
export function hoistColumns(dataFrame: DG.DataFrame, columns: DG.Column[]): HoistedColumn[] {
  const hoisted: HoistedColumn[] = [];
  for (const col of columns) {
    const template = migrateDesirability(JSON.parse(col.getTag('desirabilityTemplate')));
    const isCat = !isNumerical(template);
    const cats = isCat ? col.categories : null;
    const wc = template.weightColumn ? dataFrame.col(template.weightColumn) : null;
    hoisted.push({
      col,
      template,
      raw: col.getRawData(),
      // For categorical columns the null sentinel is the raw index of the '' category — isNone(i) ⇔ raw[i] === valNull.
      valNull: isCat ? cats!.indexOf('') : (col.type === DG.COLUMN_TYPE.INT ? DG.INT_NULL : DG.FLOAT_NULL),
      isCat,
      cats,
      catScore: isCat ?
        new Map((template as CategoricalDesirability).categories.map((c) => [c.name, c.desirability])) : null,
      staticW: template.weight,
      rawW: wc ? wc.getRawData() : null,
      wNull: wc && wc.type === DG.COLUMN_TYPE.INT ? DG.INT_NULL : DG.FLOAT_NULL,
    });
  }
  return hoisted;
}

/// The expensive, cacheable phase: column-major sweep of the raw array mapping values to desirability scores
/// (desirabilityScore inlined). Depends only on column data + line + missing-value handling, never on weights.
export function mapColumnDesirability(h: HoistedColumn, rowCount: number): ColumnDesirability {
  const {raw, valNull, isCat, cats, catScore} = h;
  const line = isCat ? null : (h.template as NumericalDesirability).line;
  const inv = !isCat && (h.template as NumericalDesirability).inverted === true;
  const segN = line ? line.length : 0;
  const loX = segN ? line![0][0] : 0;
  const hiX = segN ? line![segN - 1][0] : 0;
  const mv = h.template.missingValues;
  // The RowState a missing value resolves to (Bail/Skip mark the row; Contribute means substitute defaultScore).
  const missState = !mv || mv.strategy === 'exclude' ? RowState.Bail :
    mv.strategy === 'skip' ? RowState.Skip : RowState.Contribute;
  const defaultScore = mv && mv.strategy === 'default' ? mv.score : 0;

  const D = new Float32Array(rowCount);
  let state: Uint8Array | null = null; // lazily allocated on the first skip/bail row

  for (let i = 0; i < rowCount; ++i) {
    if (raw[i] === valNull) {
      if (missState !== RowState.Contribute) {
        (state ??= new Uint8Array(rowCount))[i] = missState;
        D[i] = NaN;
      }
      else
        D[i] = defaultScore;
      continue;
    }
    if (isCat) {
      const s = catScore!.get(cats![raw[i]]);
      if (s == null) {
        (state ??= new Uint8Array(rowCount))[i] = RowState.Bail;
        D[i] = NaN;
      }
      else
        D[i] = s;
      continue;
    }
    // Inlined desirabilityScore over the column's line.
    const x = raw[i];
    let score = 0;
    let matched = false;
    if (segN !== 0 && x >= loX && x <= hiX) {
      for (let k = 0; k < segN - 1; ++k) {
        const x1 = line![k][0];
        const x2 = line![k + 1][0];
        if (x >= x1 && x <= x2) {
          const y1 = line![k][1];
          score = x1 === x2 ? y1 : y1 + (line![k + 1][1] - y1) * (x - x1) / (x2 - x1);
          matched = true;
          break;
        }
      }
    }
    // Out-of-range / no matching segment always scores 0, regardless of inversion; only the defined curve flips.
    D[i] = matched ? (inv ? 1 - score : score) : 0;
  }
  return {D, state};
}

/// The cheap phase (always re-run; depends on weights + aggregation): folds the per-column desirabilities
/// into the final score, column-major into per-row accumulators.
export function reduceMpo(
  maps: ColumnDesirability[], hoisted: HoistedColumn[], aggregation: WeightedAggregation, rowCount: number,
): {scores: Float64Array, minPositive: number} {
  const aggCode = AGG_CODE[aggregation];
  if (aggCode === undefined)
    throw new Error(`Unknown aggregation type: ${aggregation}`);

  const needWeightSum = aggCode === AggregationCode.Average || aggCode === AggregationCode.Geomean;
  const out = new Float64Array(rowCount); // final per-row score
  const bail = new Uint8Array(rowCount); // 1 → 'exclude' missing or unmatched category → null score
  const cnt = new Int32Array(rowCount); // contributing columns (0 → no data → null score)
  const acc1 = new Float64Array(rowCount); // primary accumulator
  const acc2 = needWeightSum ? new Float64Array(rowCount) : null; // Σweight (Average) / total weight (Geomean)

  switch (aggCode) {
  case AggregationCode.Product:
  case AggregationCode.Geomean:
    acc1.fill(1);
    break;
  case AggregationCode.Min:
    acc1.fill(Infinity);
    break;
  case AggregationCode.Max:
    acc1.fill(-Infinity);
    break;
  }

  for (let j = 0; j < maps.length; ++j) {
    const {D, state: st} = maps[j];
    const {rawW: rw, wNull: wN, staticW: sw} = hoisted[j];

    for (let i = 0; i < rowCount; ++i) {
      if (bail[i])
        continue;
      if (st !== null) {
        const s = st[i];
        if (s === RowState.Bail) {
          bail[i] = 1;
          continue;
        }
        if (s === RowState.Skip)
          continue;
      }
      const score = D[i];
      const w = rw && rw[i] !== wN ? (rw[i] < 0 ? 0 : rw[i] > 1 ? 1 : rw[i]) : sw;
      ++cnt[i];
      switch (aggCode) {
      case AggregationCode.Sum:
      case AggregationCode.Average:
        acc1[i] += score * w;
        break;
      case AggregationCode.Product:
      case AggregationCode.Geomean:
        acc1[i] *= Math.pow(score, w);
        break;
      case AggregationCode.Min:
        acc1[i] = Math.min(acc1[i], Math.pow(score, w));
        break;
      case AggregationCode.Max:
        acc1[i] = Math.max(acc1[i], Math.pow(score, w));
        break;
      }
      if (acc2 !== null) // Average / Geomean: accumulate Σweight
        acc2[i] += w;
    }
  }

  let minPositive = Infinity;
  for (let i = 0; i < rowCount; ++i) {
    if (bail[i] || cnt[i] === 0)
      out[i] = DG.FLOAT_NULL;
    else if (aggCode === AggregationCode.Average)
      out[i] = acc1[i] / acc2![i];
    else if (aggCode === AggregationCode.Geomean)
      out[i] = Math.pow(acc1[i], 1 / acc2![i]); // (Π sᵏ^wᵏ)^(1/Σw) ≡ Π sᵏ^(wᵏ/Σw)
    else
      out[i] = acc1[i];
    if (out[i] > 0 && out[i] !== DG.FLOAT_NULL && out[i] < minPositive)
      minPositive = out[i];
  }
  return {scores: out, minPositive};
}

/// Builds the per-property desirability output columns from the mapping results (state 1/2 rows → null).
export function buildDesirabilityColumns(
  columns: DG.Column[], maps: ColumnDesirability[], profileName: string, rowCount: number,
): DG.Column[] {
  const result: DG.Column[] = new Array(columns.length);
  for (let j = 0; j < columns.length; ++j) {
    const {D, state} = maps[j];
    const dvals = new Float32Array(rowCount);
    for (let i = 0; i < rowCount; ++i)
      dvals[i] = state && state[i] !== RowState.Contribute ? DG.FLOAT_NULL : D[i];
    const c = DG.Column.float(`${columns[j].name} (${profileName} desirability)`, rowCount);
    c.setRawData(dvals, false);
    result[j] = c;
  }
  return result;
}

/// Stateful MPO engine for the interactive preview: same contract as mpo(), but caches each column's
/// desirability mapping (the expensive phase) keyed by df.id|column, so editing one curve recomputes only that
/// column while the rest reuse the cache and only the cheap reduction re-runs. Call release() on teardown.
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

    const {scores, minPositive} = reduceMpo(maps, hoisted, aggregation, rowCount);
    resultColumn.setRawData(scores); // notify=true: refreshes viewers
    resultColumn.meta.format = minPositive < 1e-5 ? 'scientific' : '0.00000';

    const desirabilityColumns = createDesirabilityColumns ?
      buildDesirabilityColumns(columns, maps, profileName, rowCount) : undefined;

    return {scoreColumn: resultColumn, desirabilityColumns};
  }

  /// Drops all cached desirability arrays. Call on dialog close / dataframe switch.
  release(): void {
    this.cache = new Map();
  }
}

/// Calculates the multi parameter optimization score, 0-100, 100 is the maximum.
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

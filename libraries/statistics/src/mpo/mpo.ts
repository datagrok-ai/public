/* eslint-disable guard-for-in */
/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import {TemplateFunction} from '../compute-functions/types';

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
}

export type NumericalDesirability = BasePropertyDesirability & {
  functionType: 'numerical';
  line: DesirabilityLine;
  min?: number; /// min value of the property (optional; used for editing the line)
  max?: number; /// max value of the property (optional; used for editing the line)

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

/** Calculates the multi parameter optimization score, 0-100, 100 is the maximum */
export function mpo(
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

  const desirabilityTemplates: PropertyDesirability[] = [];
  const weightColumns: (DG.Column | null)[] = [];
  for (const column of columns) {
    const tag = column.getTag('desirabilityTemplate');
    const d = migrateDesirability(JSON.parse(tag));
    desirabilityTemplates.push(d);
    weightColumns.push(d.weightColumn ? dataFrame.col(d.weightColumn) ?? null : null);
  }

  const desirabilityColumns = createDesirabilityColumns ?
    columns.map((col) => DG.Column.float(`${col.name} (${profileName} desirability)`, rowCount)) : undefined;

  resultColumn.init((i) => {
    const scores: number[] = [];
    const weights: number[] = [];

    for (let j = 0; j < columns.length; j++) {
      const desirability = desirabilityTemplates[j];
      const value = columns[j].get(i);

      let score: number | null;

      if (columns[j].isNone(i)) {
        const mv = desirability.missingValues;
        if (!mv || mv.strategy === 'exclude')
          return NaN;
        if (mv.strategy === 'skip')
          continue;
        score = mv.score;
      }
      else {
        score = isNumerical(desirability) ?
          desirabilityScore(value, desirability.line) :
          categoricalDesirabilityScore(String(value), desirability);
      }

      if (score === null)
        return NaN;

      if (createDesirabilityColumns)
        desirabilityColumns![j].set(i, score, false);

      scores.push(score);
      const wCol = weightColumns[j];
      const w = wCol && !wCol.isNone(i) ? Math.max(0, Math.min(1, wCol.get(i))) : desirability.weight;
      weights.push(w);
    }

    if (scores.length === 0)
      return NaN;

    return aggregate(scores, weights, aggregation);
  });

  return {scoreColumn: resultColumn, desirabilityColumns};
}

export function aggregate(scores: number[], weights: number[], aggregation: WeightedAggregation): number {
  switch (aggregation) {
  case 'Sum':
    return scores.reduce((sum, s, idx) => sum + s * weights[idx], 0);

  case 'Average':
    const totalWeight = weights.reduce((sum, w) => sum + w, 0);
    return scores.reduce((sum, s, idx) => sum + s * weights[idx], 0) / totalWeight;

  case 'Product':
    return scores.reduce((prod, s, idx) => prod * Math.pow(s, weights[idx]), 1);

  case 'Geomean':
    const totalW = weights.reduce((sum, w) => sum + w, 0);
    return scores.reduce((prod, s, idx) => prod * Math.pow(s, weights[idx] / totalW), 1);

  case 'Min':
    return Math.min(...scores.map((s, idx) => Math.pow(s, weights[idx])));

  case 'Max':
    return Math.max(...scores.map((s, idx) => Math.pow(s, weights[idx])));

  default:
    throw new Error(`Unknown aggregation type: ${aggregation}`);
  }
}

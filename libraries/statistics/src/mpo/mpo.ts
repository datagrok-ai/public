/* eslint-disable guard-for-in */
/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';

/// An array of [x, y] points representing the desirability line
/// [x, y] pairs are sorted by x in ascending order
export type DesirabilityLine = number[][];

export type DesirabilityMode = 'freeform' | 'gaussian' | 'sigmoid';

type BasePropertyDesirability = {
  weight: number; /// 0-1
  defaultScore?: number;
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

export function migrateDesirability(raw: any): PropertyDesirability {
  if (raw.functionType)
    return raw;
  if (raw.categories)
    return {...raw, functionType: 'categorical'};
  return {...raw, functionType: 'numerical'};
}

/// A map of desirability lines with their weights
export type DesirabilityProfile = {
  name: string;
  description: string;
  properties: { [key: string]: PropertyDesirability };
}

export const WEIGHTED_AGGREGATIONS = ['Average', 'Sum', 'Product', 'Geomean', 'Min', 'Max'] as const;
export const WEIGHTED_AGGREGATIONS_LIST: WeightedAggregation[] = [...WEIGHTED_AGGREGATIONS];
export type WeightedAggregation = typeof WEIGHTED_AGGREGATIONS[number];

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
  return found?.desirability ?? prop.defaultScore ?? null;
}

/** Calculates the multi parameter optimization score, 0-100, 100 is the maximum */
export function mpo(
  dataFrame: DG.DataFrame,
  columns: DG.Column[],
  profileName: string,
  aggregation: WeightedAggregation,
): DG.Column {
  if (columns.length === 0)
    throw new Error('No columns provided for MPO calculation.');

  const resultColumn = dataFrame.col(profileName) ?? DG.Column.float(profileName, columns[0].length);

  const desirabilityTemplates = columns.map((column) => {
    const tag = column.getTag('desirabilityTemplate');
    return migrateDesirability(JSON.parse(tag));
  });

  resultColumn.init((i) => {
    const scores: number[] = [];
    const weights: number[] = [];

    for (let j = 0; j < columns.length; j++) {
      const desirability = desirabilityTemplates[j];
      const value = columns[j].get(i);

      if (columns[j].isNone(i))
        return desirability.defaultScore ?? NaN;

      const score = isNumerical(desirability) ?
        desirabilityScore(value, desirability.line) :
        categoricalDesirabilityScore(value, desirability);

      if (score === null)
        return NaN;

      scores.push(score);
      weights.push(desirability.weight);
    }

    return aggregate(scores, weights, aggregation);
  });

  return resultColumn;
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

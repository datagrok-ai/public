/* eslint-disable no-unused-vars */
/* eslint-disable max-len */
/* eslint-disable valid-jsdoc */

// General metrics tools

import * as DG from 'datagrok-api/dg';

import {DEFAULT_WEIGHT} from './constants';

/** General metric types. */
enum GENERAL_METRIC_TYPE {
  EQUALITY = 'coincidence',
};
const GEN_METRIC_TYPES_ARR = [GENERAL_METRIC_TYPE.EQUALITY];

/** String metric types. */
export enum STR_METRIC_TYPE {
  STEMMING_BASED = 'stemming',
  EQUALITY = 'coincidence',
};
const STR_METRIC_TYPES_ARR = [STR_METRIC_TYPE.EQUALITY, STR_METRIC_TYPE.STEMMING_BASED];

/** Numerical metric types. */
enum NUM_METRIC_TYPE {
  DIFFERENCE = 'difference',
  EQUALITY = 'coincidence',
};
const NUM_METRIC_TYPES_ARR = [NUM_METRIC_TYPE.EQUALITY, NUM_METRIC_TYPE.DIFFERENCE];

/** Default metric types. */
enum DEFAULT_METRIC {
  STR = STR_METRIC_TYPE.EQUALITY,
  NUM = NUM_METRIC_TYPE.DIFFERENCE,
  GENERAL = GENERAL_METRIC_TYPE.EQUALITY,
};

type MetricType = GENERAL_METRIC_TYPE | STR_METRIC_TYPE | NUM_METRIC_TYPE | DEFAULT_METRIC;

/** Metric specification. */
export type MetricInfo = {
  weight: number,
  type: MetricType,
};

/** Distance types. */
export enum DISTANCE_TYPE {
  EUCLIDEAN = 'Euclidian',
  MANHATTAN = 'Manhattan',
};

export const DIST_TYPES_ARR = [DISTANCE_TYPE.EUCLIDEAN, DISTANCE_TYPE.MANHATTAN];

/** Return choices list of metric types with respect to the column type. */
export function getMetricTypesChoicesList(col: DG.Column): string[] {
  switch (col.type) {
  case DG.COLUMN_TYPE.STRING:
    return STR_METRIC_TYPES_ARR;

  case DG.COLUMN_TYPE.FLOAT:
  case DG.COLUMN_TYPE.INT:
    return NUM_METRIC_TYPES_ARR;

  default:
    return GEN_METRIC_TYPES_ARR;
  };
}

/** Return default metric specification with respect to the column type. */
export function getDefaultMetric(col: DG.Column): MetricInfo {
  let type: MetricType = DEFAULT_METRIC.GENERAL;

  switch (col.type) {
  case DG.COLUMN_TYPE.STRING:
    type = DEFAULT_METRIC.STR;
    break;

  case DG.COLUMN_TYPE.FLOAT:
  case DG.COLUMN_TYPE.INT:
    type = DEFAULT_METRIC.NUM;
    break;

  default:
    break;
  };

  return {weight: DEFAULT_WEIGHT, type: type};
}

/** Return a map with a single metric corresponding to the input column. */
export function getSingleSourceMetricsMap(col: DG.Column): Map<string, MetricInfo> {
  const map = new Map<string, MetricInfo>();
  let type: MetricType = DEFAULT_METRIC.GENERAL;

  switch (col.type) {
  case DG.COLUMN_TYPE.STRING:
    type = STR_METRIC_TYPE.STEMMING_BASED;
    break;

  case DG.COLUMN_TYPE.FLOAT:
  case DG.COLUMN_TYPE.INT:
    type = DEFAULT_METRIC.NUM;
    break;

  default:
    break;
  };

  map.set(col.name, {weight: DEFAULT_WEIGHT, type: type});

  return map;
}

/** Return default type of distance between elements with the specified features. */
export function getDefaultDistnce(): DISTANCE_TYPE {
  return DIST_TYPES_ARR[0];
}

/** Return a function that computes the distance between elements with the specified features (actually, norm of the input vector). */
export function getDistanceFn(type: DISTANCE_TYPE, dimensionality: number): (vector: Float32Array) => number {
  if (dimensionality < 1)
    throw new Error('Incorrect dimensionality');

  switch (type) {
  case DISTANCE_TYPE.EUCLIDEAN:
    return (vector: Float32Array) => {
      let s = 0;

      for (let i = 0; i < dimensionality; ++i)
        s += vector[i] ** 2;

      return Math.sqrt(s);
    };

  case DISTANCE_TYPE.MANHATTAN:
    return (vector: Float32Array) => {
      let s = 0;

      for (let i = 0; i < dimensionality; ++i)
        s += Math.abs(vector[i]);

      return s;
    };

  default:
    throw new Error('Non-supported distance');
  }
}

/** Return metric function for computing distance between elements (within a single feature). */
export function getMetricFn(info: MetricInfo): (a: any, b: any) => number {
  switch (info.type) {
  case NUM_METRIC_TYPE.DIFFERENCE:
    return (a: any, b: any) => info.weight * (a - b);

  default:
    return (a: any, b: any) => {
      if (a === b)
        return 0;
      return info.weight;
    };
  }
}

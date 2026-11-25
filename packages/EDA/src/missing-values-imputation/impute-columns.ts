import * as DG from 'datagrok-api/dg';

import {ERROR_MSG} from './ui-constants';
import {METRIC_TYPE, DISTANCE_TYPE, MetricInfo, DEFAULT, impute, getMissingValsIndices} from './knn-imputer';

/** Setting of the feature metric inputs */
type FeatureInputSettings = {
  defaultWeight: number,
  defaultMetric: METRIC_TYPE,
  availableMetrics: METRIC_TYPE[],
};

/** Return default setting of the feature metric inputs */
export function getFeatureInputSettings(type: DG.COLUMN_TYPE): FeatureInputSettings {
  switch (type) {
  case DG.COLUMN_TYPE.STRING:
  case DG.COLUMN_TYPE.DATE_TIME:
    return {
      defaultWeight: DEFAULT.WEIGHT,
      defaultMetric: METRIC_TYPE.ONE_HOT,
      availableMetrics: [METRIC_TYPE.ONE_HOT],
    };

  case DG.COLUMN_TYPE.INT:
  case DG.COLUMN_TYPE.FLOAT:
  case DG.COLUMN_TYPE.QNUM:
    return {
      defaultWeight: DEFAULT.WEIGHT,
      defaultMetric: METRIC_TYPE.DIFFERENCE,
      availableMetrics: [METRIC_TYPE.DIFFERENCE, METRIC_TYPE.ONE_HOT],
    };

  default:
    throw new Error(ERROR_MSG.UNSUPPORTED_COLUMN_TYPE);
  }
}

/** Fill missing values of the columns in place with the default metrics; return the cells that stay empty */
export function imputeColumns(df: DG.DataFrame, columns: string[], features: string[], neighbors: number,
  distance: DISTANCE_TYPE): Map<string, number[]> {
  const misValsInds = getMissingValsIndices(df.columns.byNames(columns));
  const targets = columns.filter((name) => misValsInds.has(name));

  if (targets.length === 0)
    return new Map();

  const featuresMetrics = new Map<string, MetricInfo>(features.map((name) => [name, {
    weight: DEFAULT.WEIGHT,
    type: getFeatureInputSettings(df.getCol(name).type as DG.COLUMN_TYPE).defaultMetric,
  }]));

  return impute(df, targets, featuresMetrics, misValsInds, distance, neighbors, true);
}

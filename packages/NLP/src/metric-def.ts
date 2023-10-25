/** */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const DEFAULT_WEIGHT = 1;

enum GENERAL_METRIC_TYPE {
  EQUALITY = 'coincidence',
};
const GENERAL_METRIC_TYPES = [GENERAL_METRIC_TYPE.EQUALITY];

enum STR_METRIC_TYPE {
  STEMMING_BASED = 'stemming',
  EQUALITY = 'coincidence',
};
const STRING_METRIC_TYPES = [STR_METRIC_TYPE.EQUALITY, STR_METRIC_TYPE.STEMMING_BASED];

enum NUM_METRIC_TYPE {
  ABS_DIFFERENCE = 'difference',
  EQUALITY = 'coincidence',
};
const NUM_METRIC_TYPES = [NUM_METRIC_TYPE.EQUALITY, NUM_METRIC_TYPE.ABS_DIFFERENCE];

enum DEFAULT_METRIC {
  STR = STR_METRIC_TYPE.STEMMING_BASED,
  NUM = NUM_METRIC_TYPE.ABS_DIFFERENCE,
  GENERAL = GENERAL_METRIC_TYPE.EQUALITY,
};

export type MetricType = GENERAL_METRIC_TYPE | STR_METRIC_TYPE | NUM_METRIC_TYPE | DEFAULT_METRIC;

export type MetricInfo = {
  weight: number,
  type: MetricType,
};

export enum DISTANCE_TYPE {
  EUCLIDEAN = 'Euclidian',
  MANHATTAN = 'Manhattan',
};

export const DISTANCE_TYPES = [DISTANCE_TYPE.EUCLIDEAN, DISTANCE_TYPE.MANHATTAN];

export function getChoicesList(col: DG.Column): string[] {
  switch (col.type) {
    case DG.COLUMN_TYPE.STRING:
      return STRING_METRIC_TYPES;

    case DG.COLUMN_TYPE.FLOAT:
    case DG.COLUMN_TYPE.INT:
      return NUM_METRIC_TYPES;
      
    default:
      return GENERAL_METRIC_TYPES;
  };  
}

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

export function getOneHotMetricsMap(df: DG.DataFrame, target: DG.Column): Map<string, MetricInfo> {
  const map = new Map<string, MetricInfo>();
  let type: MetricType = DEFAULT_METRIC.GENERAL;
  
  switch (target.type) {
  
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

  map.set(target.name, {weight: DEFAULT_WEIGHT, type: type});
    
  return map;
}

export function getDefaultDistnce(): DISTANCE_TYPE { return DISTANCE_TYPES[0]; }

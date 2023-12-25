import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {ERROR_MSG, COPY_SUFFIX} from './ui-constants';

/** */
function getNullValue(col: DG.Column): number {
  switch (col.type) {
    case DG.COLUMN_TYPE.INT:        
      return DG.INT_NULL;

    case DG.COLUMN_TYPE.FLOAT:
      return DG.FLOAT_NULL;
    
    case DG.COLUMN_TYPE.QNUM:
      return DG.FLOAT_NULL;

    case DG.COLUMN_TYPE.DATE_TIME:
      return DG.FLOAT_NULL;

    case DG.COLUMN_TYPE.STRING:
      return col.max;
  
    default:
      throw new Error(ERROR_MSG.UNSUPPORTED_COLUMN_TYPE);
  }  
}

/** */
export function simpleImpute(col: DG.Column, strategy: string, inPlace: boolean): DG.Column {
  if (col.stats.missingValueCount === 0)
    throw new Error(ERROR_MSG.NO_MISSING_VALUES);

  const nullValue = getNullValue(col);
  const fillValue = 0;//getFillValue(col, strategy);

  const len = col.length;
  const source = col.getRawData();    

  if (inPlace) {
    for (let i = 0; i < len; ++i)
      if (source[i] === nullValue)
        source[i] = fillValue;

    const buf = col.get(0); 
    col.set(0, col.get(1));
    col.set(0, buf);

    return col;
  }

  //@ts-ignore
  const copy = col.clone();
  copy.name = `${col.name}${COPY_SUFFIX}`;

  const copyRaw = copy.getRawData();

  for (let i = 0; i < len; ++i)
    copyRaw[i] = (source[i] !== nullValue) ? source[i] : fillValue;

  /*console.log(source);
  console.log(copyRaw);*/

  return copy;
}

/** */
export enum METRIC_TYPE {
  ONE_HOT = 'One-hot',
  DIFFERENCE = 'Difference',
};

/** Distance types. */
export enum DISTANCE_TYPE {
  EUCLIDEAN = 'Euclidian',
  MANHATTAN = 'Manhattan',
};

/** Metric specification. */
export type MetricInfo = {
  weight: number,
  type: METRIC_TYPE,
};

/** */
export enum DEFAULT {
  WEIGHT = 1,
  NEIGHBORS = 2,
  IN_PLACE = 1,
};

/** */
export const MIN_NEIGHBORS = 1;

/** */
export function impute(df: DG.DataFrame, targetCols: DG.Column[], featuresMetrics: Map<string, MetricInfo>,
  distance: DISTANCE_TYPE, neighbors: number, inPlace: boolean) 
{
  console.log(df);
  console.log('Target');
  console.log(targetCols);
  console.log('Features');
  console.log(featuresMetrics);
  console.log(distance);
  console.log(neighbors);
  console.log(inPlace);
}

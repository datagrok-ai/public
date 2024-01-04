// Simple method for missing values imputation: mean, median, most frequent

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import { SUPPORTED_COLUMN_TYPES, getNullValue } from './knn-imputer';
import { ERROR_MSG } from './ui-constants';

/** Fill value type */
export enum FILL_VALUE {
  MEAN = 'Mean',
  MEDIAN = 'Median',
  MOST_FREQUENT = 'First non-null', // 'Most frequent' TODO: replace
};

/** Simple imputer task: indeces of elements with missing values & fill value */
export type SimpleImputTask = {
  indeces: number[],
  fillValueType: FILL_VALUE,
};

/** Returns first non-null value */
function getFirstNonNull<T>(col: DG.Column<T>): T {
  const nullValue = getNullValue(col);
  const raw = col.getRawData();
  const len = raw.length;

  for (let i = 0; i < len; ++i)
    if (raw[i] !== nullValue)
      return col.get(i)!;

  throw new Error(ERROR_MSG.EMPTY_COLUMN);  
}

/** */
function getFillValue<T>(col: DG.Column<T>, fillValueType: FILL_VALUE): T {
  switch (fillValueType) {
    case FILL_VALUE.MEAN:      
      return col.stats.avg as T;

    case FILL_VALUE.MEDIAN:
      return col.stats.med as T;

    case FILL_VALUE.MOST_FREQUENT:
      return getFirstNonNull(col); // TODO: to replace by most frequent, implemented by Davit
  
    default:
      throw new Error(ERROR_MSG.UNSUPPORTED_FILL_VALUE_TYPE);
  }
}

/** Perform missing values imputation using the simple approach */
export function simpleImput(df: DG.DataFrame, taskMap: Map<string, SimpleImputTask>): void {
  taskMap.forEach((task, colName) => {
    const col = df.col(colName);
    if (col !== null) {
      if (!SUPPORTED_COLUMN_TYPES.includes(col.type))
        throw new Error(ERROR_MSG.UNSUPPORTED_COLUMN_TYPE);        

      const fillVal = getFillValue(col, task.fillValueType);
      task.indeces.forEach((idx) => col.set(idx, fillVal));
    }
  });
}

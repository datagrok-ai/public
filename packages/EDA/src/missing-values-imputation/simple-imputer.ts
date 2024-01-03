// Simple method for missing values imputation: mean, median, most frequent

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

/** Fill value type */
export enum FILL_VALUE {
  MEAN = 'Mean',
  MEDIAN = 'Median',
  MOST_FREQUENT = 'Most frequent',
};

/** Simple imputer task: indeces of elements with missing values & fill value */
export type SimpleImputTask = {
  indeces: number[],
  fillValue: FILL_VALUE,
};

/** Perform missing values imputation using the simple approach */
export function simpleImput(df: DG.DataFrame, task: Map<string, SimpleImputTask>): void {
  console.log(task);
}

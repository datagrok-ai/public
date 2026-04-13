/**
 * Statistical utility functions (wrappers around jStat).
 */

import {jStat} from 'jstat';

/** Cumulative distribution function of Student's t-distribution. */
export function tdistCdf(t: number, df: number): number {
  return jStat.studentt.cdf(t, df);
}

/** Two-tailed p-value for a given t-statistic and degrees of freedom. */
export function twoTailPvalue(tStat: number, df: number): number {
  return 2 * (1 - tdistCdf(Math.abs(tStat), df));
}

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

//@ts-ignore: no types
import * as jStat from 'jstat';

import {DescriptorStatistics} from './pmpo-defs';

export function getDesiredTables(df: DG.DataFrame, desirability: DG.Column) {
  const groups = df.groupBy([desirability.name]).getGroups() as any;
  let desired: DG.DataFrame;
  let nonDesired: DG.DataFrame;

  for (const name in groups) {
    if (name.toLowerCase().includes('true'))
      desired = groups[name];
    else
      nonDesired = groups[name];
  }

  //@ts-ignore
  return {desired, nonDesired};
} // getDesiredTables

/* Welch two-sample t-test (two-sided) */
export function getDescriptorStatistics(des: DG.Column, nonDes: DG.Column): DescriptorStatistics {
  const desLen = des.length;
  const nonDesLen = nonDes.length;
  if (desLen < 2 || nonDesLen < 2) {
    throw new Error(`Failed to compute the "${des.name}" descriptor statistics: 
      both samples must have at least two observations.`);
  }

  const desAvg = des.stats.avg;
  const nonDesAvg = nonDes.stats.avg;
  const desVar = des.stats.variance;
  const nonDesVar = nonDes.stats.variance;
  const desStd = des.stats.stdev;
  const nonDesStd = nonDes.stats.stdev;

  const se = Math.sqrt(desVar / desLen + nonDesVar / nonDesLen);
  if (se === 0) {
    throw new Error(`Failed to compute the "${des.name}" descriptor statistics: 
      zero variance.`);
  }

  const t = (desAvg - nonDesAvg) / se;

  // Welch–Satterthwaite degrees of freedom
  const numerator = (desVar / desLen + nonDesVar / nonDesLen) ** 2;
  const denom = (desVar * desVar) / (desLen * desLen * (desLen - 1)) +
    (nonDesVar * nonDesVar) / (nonDesLen * nonDesLen * (nonDesLen - 1));
  const df = numerator / denom;

  // two-sided p-value
  const cdf = jStat.studentt.cdf(Math.abs(t), df);
  const pValue = 2 * (1 - cdf);

  return {
    desAvg: desAvg,
    desStd: desStd,
    desLen: desLen,
    nonDesAvg: nonDesAvg,
    nonDesStd: nonDesStd,
    nonSesLen: nonDesLen,
    tstat: t,
    pValue: pValue,
    zScore: Math.abs(desAvg - nonDesAvg) / (desStd + nonDesStd),
  };
} // getDescriptorStatistics

// Statistic tools

/* REFERENCES
   
     [1] One-way analysis of variance, https://en.wikipedia.org/wiki/One-way_analysis_of_variance

     [2] G.W. Heiman. Basic Statistics for the Behavioral Sciences, 6th ed. Wadsworth Publishing, 2010
     
     [3] F-test of equality of variances, https://en.wikipedia.org/wiki/F-test_of_equality_of_variances

     [4] S. McKillup. Statistics Explained, Cambridge University Press, 2005

*/

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

//@ts-ignore: no types
import * as jStat from 'jstat';

enum ERROR_MSG {
  NON_EQUAL_FACTORS_VALUES_SIZE = 'non-equal sizes of factor and values arrays. INPUT ERROR.',
  INCORRECT_SIGNIFICANCE_LEVEL = 'incorrect significance level. It must be from the interval (0, 1). INPUT ERROR.',
  INCORRECT_SAMPLE_SIZE = 'incorrect size of sample. DATA FACTORIZAING ERROR.',
};

export type SampleData = {
  sum: number,
  sumOfSquares: number,
  size: number,
  data: number[] | undefined | null,
};

/** One-way ANOVA computation results. The classic notations are used (see [2], p. 290). */
export type OneWayAnovaTable = {
  /** sum of squares between groups */
  SSbn: number,
  /** sum of squares within groups */ 
  SSwn: number,
  /** total sum of squares */ 
  SStot: number,
  /** degrees of freedom between groups */
  DFbn: number,
  /** degrees of freedom within groups */
  DFwn: number,
  /** total degrees of freedom */
  DFtot: number,
  /** mean square between groups */
  MSbn: number,
  /** mean square within groups */
  MSwn: number,
  /** Fobt, value of F-statistics */
  Fstat: number,
  /** p-value corresponding to F-statistics */
  pValue: number,
};

type ValuesType = number[] | Float32Array | Int32Array;

export type AnovaResults = {
  summaryTable: DG.DataFrame,
  conclusion: string;
};

/** Split values by factor levels and compute sums & sums of squares. */
export function getFactorizedData(factors: any[], values: ValuesType): Map<any, SampleData> {
  const size = factors.length;
  
  if (size !== values.length)
    throw new Error(ERROR_MSG.NON_EQUAL_FACTORS_VALUES_SIZE);
  
  const factorized = new Map<any, SampleData>();  
  
  for (let i = 0; i < size; ++i) {
    const cur = factorized.get(factors[i]) ?? {sum: 0, sumOfSquares: 0, size: 0, data: []};
  
    cur.sum += values[i];
    cur.sumOfSquares += values[i] ** 2;
    ++cur.size;
    cur.data?.push(values[i]);
  
    factorized.set(factors[i], cur);
  }
  
  return factorized;
} // getFactorizedData

/** Compute sum & sums of squares with respect to factor levels. */
export function getLevelsStat(factors: any[], values: ValuesType): Map<any, SampleData> {
  const size = factors.length;
    
  if (size !== values.length)
    throw new Error(ERROR_MSG.NON_EQUAL_FACTORS_VALUES_SIZE);
    
  const factorized = new Map<any, SampleData>();  
    
  for (let i = 0; i < size; ++i) {
    const cur = factorized.get(factors[i]) ?? {sum: 0, sumOfSquares: 0, size: 0, data: null};
    
    cur.sum += values[i];
    cur.sumOfSquares += values[i] ** 2;
    ++cur.size;
    
    factorized.set(factors[i], cur);
  }
    
  return factorized;
 } // getLevelsStat

/** Perform one-way ANOVA computations using factorized data. */
export function computeOneWayAnovaTable(factorizedData: Map<any, SampleData>): OneWayAnovaTable {
  // Further, notations and formulas from (see [2], p. 290) are used.

  const K = factorizedData.size;
  
  let sum = 0;
  let sumOfSquares = 0;
  let N = 0;
  let buf = 0; 
  
  factorizedData.forEach((levelData, key, map) => {
    sum += levelData.sum;
    sumOfSquares += levelData.sumOfSquares;
    N += levelData.size;
    buf += (levelData.sum) ** 2 / levelData.size;
  });
  
  const SStot = sumOfSquares - sum ** 2 / N;  
  const SSbn = buf - sum ** 2 / N;
  const SSwn = SStot - SSbn;
 
  const DFbn = K - 1;
  const DFwn = N - K;
  const DFtot = N - 1;
  
  const MSbn = SSbn / DFbn;
  const MSwn = SSwn / DFwn;
  
  const Fstat = MSbn / MSwn;
  
  return {
    SSbn: SSbn,
    SSwn: SSwn,
    SStot: SStot,
    DFbn: DFbn,
    DFwn: DFwn,
    DFtot: DFtot,
    MSbn: MSbn,
    MSwn: MSwn,
    Fstat: Fstat,
    pValue: 1 - jStat.centralF.cdf(Fstat, DFbn, DFwn)
  };
} // computeOneWayAnovaTable

/** Create dataframe with one-way ANOVA results. */
function getOneWayAnovaDataframe(anovaTable: OneWayAnovaTable, alpha: number, Fcritical: number): DG.DataFrame {
  return DG.DataFrame.fromColumns([
    DG.Column.fromStrings('Source of variance', ['Between groups', 'Within groups', 'Total']),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'Sum of squares', [anovaTable.SSbn, anovaTable.SSwn, anovaTable.SStot]),
    DG.Column.fromList(DG.COLUMN_TYPE.INT, 'Degrees of freedom', [anovaTable.DFbn, anovaTable.DFwn, anovaTable.DFtot]),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'Mean square', [anovaTable.MSbn, anovaTable.MSwn, null]),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'F-statistics', [anovaTable.Fstat, null, null]),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'p-value', [anovaTable.pValue, null, null]),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, `${alpha}-critical value`, [Fcritical, null, null]),
  ]);
} // getOneWayAnovaDataframe

/** Check correctness of significance level. */
function checkSignificanceLevel(alpha: number) {
  if ((alpha <= 0) || (alpha >= 1))
    throw new Error(ERROR_MSG.INCORRECT_SIGNIFICANCE_LEVEL);
}

/** Compute unbiased variance.*/
export function getVariance(sampleData: SampleData): number {
  // The applied formulas can be found in [4] (see p. 63)
  const size = sampleData.size;

  if (size <= 0)
    throw new Error(ERROR_MSG.INCORRECT_SAMPLE_SIZE);

  if (size === 1)
    return 0;
  
  return (sampleData.sumOfSquares - (sampleData.sum) ** 2 / size) / (size - 1);
} // getVariance

/** Checks equality of variances of 2 samples. F-test is performed.*/
export function areVariancesOfTwoSamplesEqual(xSampleData: SampleData, ySampleData: SampleData, alpha: number  = 0.05): boolean {
  // The applied approach can be found in [3]
  const xSampleVariance = getVariance(xSampleData);
  const ySampleVariance = getVariance(ySampleData);

  if (ySampleVariance === 0)
    return (xSampleVariance === ySampleVariance);

  const Fstat = xSampleVariance / ySampleVariance;
  const Fcritical = jStat.centralF.inv(1 - alpha, xSampleData.size - 1, ySampleData.size - 1);

  return (Fstat < Fcritical);
} // areVariancesOfTwoSamplesEqual

/** Checks equality of variances of factorized data. */
export function areVariancesOfFactorizedDataEqual(factorizedData: Map<any, SampleData>): boolean {
  const iter = factorizedData.keys();
  let key = iter.next();
  const xSampleData = factorizedData.get(key.value);

  while (true) {
    key = iter.next();

    if (key.done)
      break;

    if(!areVariancesOfTwoSamplesEqual(xSampleData!, factorizedData.get(key.value)!))
      return false;
  }

  return true;
} // areVariancesOfFactorizedDataEqual

/** Perform one-way analysis of variances. */
export function oneWayAnova(
  factorsCol: DG.Column, 
  valuesCol: DG.Column,
  alpha: number  = 0.05,
  toCheckNormality: boolean = false, 
  toCheckVariancesEqaulity: boolean = false): AnovaResults
{
  checkSignificanceLevel(alpha);

  const values = ((valuesCol.type === DG.COLUMN_TYPE.FLOAT) || (valuesCol.type === DG.COLUMN_TYPE.INT)) 
    ? valuesCol.getRawData() 
    : valuesCol.toList();

  const factors = factorsCol.toList();

  const factorizedData = getLevelsStat(factors, values as ValuesType);

  if (toCheckNormality) {
    // TODO: add this
  }

  if (toCheckVariancesEqaulity)
    if(!areVariancesOfFactorizedDataEqual(factorizedData))    
      throw new Error('Variances are not equal!');
  
  const anovaTable = computeOneWayAnovaTable(factorizedData);

  const Fcritical = jStat.centralF.inv(1 - alpha, anovaTable.DFbn, anovaTable.DFwn);

  const conclusion = (anovaTable.Fstat > Fcritical)
    ? `The null hypothesis is rejected. The "${factorsCol.name}"-factor produces a significant difference in "${valuesCol.name}"-values.`
    : `Failed to reject the null hypothesis. The "${factorsCol.name}"-factor does not produce a significant difference in "${valuesCol.name}"-values.`;

  return {
    summaryTable: getOneWayAnovaDataframe(anovaTable, alpha, Fcritical),
    conclusion: conclusion
  };
} // oneWayAnova

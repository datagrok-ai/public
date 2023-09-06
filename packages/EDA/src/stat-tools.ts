// Statistic tools

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

//@ts-ignore: no types
import * as jStat from 'jstat';

enum ERROR_MSGS {
  NON_EQUAL_FACTORS_VALUES_SIZE = 'non-equal sizes of factor and values arrays.',
  INCORRECT_SIGNIFICANCE_LEVEL = 'incorrect significance level. It must be from the interval (0, 1).',
};

export type FactorLevelData = {
  sum: number,
  sumOfSquares: number,
  size: number,
  data: number[] | undefined | null,
};
 
export type OneWayAnovaTable = {
  SSbn: number, // sum of squares between groups
  SSwn: number, // sum of squares within groups
  SStot: number, // total sum of squares
  DFbn: number, // degrees of freedom between groups
  DFwn: number, // degrees of freedom within groups
  DFtot: number, // total degrees of freedom
  MSbn: number, // mean square between groups
  MSwn: number, // mean square within groups
  Fstat: number, // Fobt, value of F-statistics
  pValue: number, // p-value corresponding to F-statistics
};

type ValuesType = number[] | Float32Array | Int32Array;

export type AnovaResults = {
  summaryTable: DG.DataFrame,
  conclusion: string;
};

export function getFactorizedData(factors: any[], values: ValuesType): Map<any, FactorLevelData> {
  const size = factors.length;
  
  if (size !== values.length)
    throw new Error(ERROR_MSGS.NON_EQUAL_FACTORS_VALUES_SIZE);
  
  const factorized = new Map<any, FactorLevelData>();  
  
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
  
export function getLevelsStat(factors: any[], values: ValuesType): Map<any, FactorLevelData> {
  const size = factors.length;
    
  if (size !== values.length)
    throw new Error(ERROR_MSGS.NON_EQUAL_FACTORS_VALUES_SIZE);
    
  const factorized = new Map<any, FactorLevelData>();  
    
  for (let i = 0; i < size; ++i) {
    const cur = factorized.get(factors[i]) ?? {sum: 0, sumOfSquares: 0, size: 0, data: null};
    
    cur.sum += values[i];
    cur.sumOfSquares += values[i] ** 2;
    ++cur.size;
    
    factorized.set(factors[i], cur);
  }
    
  return factorized;
 } // getLevelsStat
  
export function computeOneWayAnovaTable(factorizedData: Map<any, FactorLevelData>): OneWayAnovaTable {
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

function checkSignificanceLevel(alpha: number) {
  if ((alpha <= 0) || (alpha >= 1))
    throw new Error(ERROR_MSGS.INCORRECT_SIGNIFICANCE_LEVEL);
}

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
    // TODO: to add normality test    
  }

  if (toCheckVariancesEqaulity) {
    // TODO: to equality of variances test
  }

  const anovaTable = computeOneWayAnovaTable(factorizedData);

  const Fcritical = jStat.centralF.inv(1 - alpha, anovaTable.DFbn, anovaTable.DFwn);

  const conclusion = (anovaTable.Fstat > Fcritical)
    ? `The null hypothesis is rejected. The "${factorsCol.name}"-factor produces a significant difference in mean "${valuesCol.name}"-values.`
    : `The null hypothesis is not rejected. The "${factorsCol.name}"-factor does not produce a significant difference in mean "${valuesCol.name}"-values.`;

  return {
    summaryTable: getOneWayAnovaDataframe(anovaTable, alpha, Fcritical),
    conclusion: conclusion
  };
} // oneWayAnova

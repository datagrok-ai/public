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
  NON_EQUAL_VARIANCES = 'variances are not equal.',
  NON_NORMAL_DISTRIB = 'non-normal distribution.',
  UNSUPPORTED_COLUMN_TYPE = 'unsupported column type.',
};

export type SampleData = {
  sum: number,
  sumOfSquares: number,
  size: number,
  data: number[] | undefined,
};

/** One-way ANOVA computation results. The classic notations are used (see [2], p. 290). */
export type OneWayAnova = {
  /** sum of squares between groups, SSbn */
  ssBn: number,
  /** sum of squares within groups, SSnn */ 
  ssWn: number,
  /** total sum of squares, SStot */ 
  ssTot: number,
  /** degrees of freedom between groups, DFbn */
  dfBn: number,
  /** degrees of freedom within groups, DFwn */
  dfWn: number,
  /** total degrees of freedom, DFtot */
  dfTot: number,
  /** mean square between groups, MSbn */
  msBn: number,
  /** mean square within groups, MSwn */
  msWn: number,
  /** Fobt, value of F-statistics, Fstat */
  fStat: number,
  /** p-value corresponding to F-statistics, pValue */
  pValue: number,
};

type ValuesType = DG.Column<DG.COLUMN_TYPE.FLOAT> | DG.Column<DG.COLUMN_TYPE.INT>;

/** Create dataframe with one-way ANOVA results. */
export function getOneWayAnovaDF(anova: OneWayAnova, alpha: number, fCritical: number, hypothesis: string, testResult: string): DG.DataFrame {
  return DG.DataFrame.fromColumns([
    DG.Column.fromStrings('Source of variance', ['Between groups', 'Within groups', 'Total', '', hypothesis, '', testResult]),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'Sum of squares', [anova.ssBn, anova.ssWn, anova.ssTot, null, null, null, null]),
    DG.Column.fromList(DG.COLUMN_TYPE.INT, 'Degrees of freedom', [anova.dfBn, anova.dfWn, anova.dfTot, null, null, null, null]),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'Mean square', [anova.msBn, anova.msWn, null, null, null, null, null]),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'F-statistics', [anova.fStat, null, null, null, null, null, null]),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'p-value', [anova.pValue, null, null, null, null, null, null]),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, `${alpha}-critical value`, [fCritical, null, null, null, null, null, null]),
  ]);
} // getOneWayAnovaDF

/** Check correctness of significance level. */
export function checkSignificanceLevel(alpha: number) {
  if ((alpha <= 0) || (alpha >= 1))
    throw new Error(ERROR_MSG.INCORRECT_SIGNIFICANCE_LEVEL);
}

/** Compute unbiased variance.*/
export function getVariance(data: SampleData): number {
  // The applied formulas can be found in [4] (see p. 63)
  const size = data.size;

  if (size <= 0)
    throw new Error(ERROR_MSG.INCORRECT_SAMPLE_SIZE);

  if (size === 1)
    return 0;
  
  return (data.sumOfSquares - (data.sum) ** 2 / size) / (size - 1);
} // getVariance

/** Check equality of variances of 2 samples. F-test is performed.*/
export function areVarsEqual(xData: SampleData, yData: SampleData, alpha: number = 0.05): boolean {
  // The applied approach can be found in [3]
  checkSignificanceLevel(alpha);
  
  const xVar = getVariance(xData);
  const yVar = getVariance(yData);

  if (yVar === 0)
    return (xVar === yVar);

  const fStat = xVar / yVar;
  const fCrit = jStat.centralF.inv(1 - alpha, xData.size - 1, yData.size - 1);

  return (fStat < fCrit);
} // areVarsEqual

export class FactorizedData extends Map<any, SampleData> {
  private isNormDistrib: boolean | undefined = undefined;

  constructor(factors: DG.Column, values: ValuesType, extractData: boolean = false, checkNormality: boolean = false, alpha: number = 0.05) { 
    super();

    if (factors.length !== values.length)
      throw new Error(ERROR_MSG.NON_EQUAL_FACTORS_VALUES_SIZE);      

    if (extractData)
      this.setData(factors, values, checkNormality, alpha) // factorized data is stored & statistics are computed
    else
      this.setStats(factors, values, checkNormality, alpha); // just statistics are computed
  }

  public isNormal(): boolean | undefined {
    return true; // TODO: after check normality feature will be implemented, replace with the following:
    //return this.isNormDistrib;
  }

  /** Check equality of variances of factorized data. */
  public areVarsEqual(alpha: number = 0.05): boolean {
    const iter = this.keys();
    let key = iter.next();
    const xData = this.get(key.value);
    
    while (true) {
      key = iter.next();
      
      if (key.done)
        break;
      
      if(!areVarsEqual(xData!, this.get(key.value)!, alpha))
        return false;
    }

    return true;    
  } // areVarsEqual

  /** Perform one-way ANOVA computations. */
  public getOneWayAnova(): OneWayAnova {
    // Further, notations and formulas from (see [2], p. 290) are used.

    const K = this.size;
  
    let sum = 0;
    let sumOfSquares = 0;
    let N = 0;
    let buf = 0; 
  
    this.forEach((levelData, key, map) => {
      sum += levelData.sum;
      sumOfSquares += levelData.sumOfSquares;
      N += levelData.size;
      buf += (levelData.sum) ** 2 / levelData.size;
    });
  
    const ssTot = sumOfSquares - sum ** 2 / N;
    const ssBn = buf - sum ** 2 / N;
    const ssWn = ssTot - ssBn;

    const dfBn = K - 1;
    const dfWn = N - K;
    const dfTot = N - 1;
    
    const msBn = ssBn / dfBn;
    const msWn = ssWn / dfWn;

    const fStat = msBn / msWn;
  
    return {
      ssBn: ssBn,
      ssWn: ssWn,
      ssTot: ssTot,
      dfBn: dfBn,
      dfWn: dfWn,
      dfTot: dfTot,
      msBn: msBn,
      msWn: msWn,
     fStat: fStat,
      pValue: 1 - jStat.centralF.cdf(fStat, dfBn, dfWn)
    };
  } // getOneWayAnova

  /** Compute sum & sums of squares with respect to factor levels. */
  private setStats(factors: DG.Column, values: ValuesType, checkNormality: boolean = false, alpha: number = 0.05): void {
    // TODO: provide check normality feature
    const size = factors.length;

    const fact = factors.toList(); // making this copy improves performance 

    switch (values.type) {
      case DG.COLUMN_TYPE.INT:
      case DG.COLUMN_TYPE.FLOAT:
        const buf = values.getRawData();

        for (let i = 0; i < size; ++i) {
          const cur = this.get(fact[i]) ?? {sum: 0, sumOfSquares: 0, size: 0, data: undefined};
    
          cur.sum += buf[i];
          cur.sumOfSquares += buf[i] ** 2;
          ++cur.size;
    
          this.set(fact[i], cur);
        }        
        break;

      default:
        throw new Error(ERROR_MSG.UNSUPPORTED_COLUMN_TYPE);
    }    
  } // setStats

  /** Compute sum & sums of squares with respect to factor levels & store factorized data. */
  private setData(factors: DG.Column, values: ValuesType, checkNormality: boolean = false, alpha: number = 0.05): void {
    // TODO: provide check normality feature
    const size = factors.length;

    const fact = factors.toList(); // making this copy improves performance

    switch (values.type) {
      case DG.COLUMN_TYPE.INT:
      case DG.COLUMN_TYPE.FLOAT:
        const buf = values.getRawData();

        for (let i = 0; i < size; ++i) {
          const cur = this.get(fact[i]) ?? {sum: 0, sumOfSquares: 0, size: 0, data: []};
    
          cur.sum += buf[i];
          cur.sumOfSquares += buf[i] ** 2;
          ++cur.size;
          cur.data?.push(buf[i]);
    
          this.set(fact[i], cur);
        }        
        break;

      default:
        throw new Error(ERROR_MSG.UNSUPPORTED_COLUMN_TYPE);
    }    
  } // setData

} // FactorizedData

/** Perform one-way analysis of variances. */
export function oneWayAnova(factors: DG.Column, values: DG.Column, alpha: number = 0.05, validate: boolean = false): DG.DataFrame {
  checkSignificanceLevel(alpha);  

  const factorized = new FactorizedData(factors, values, false, validate, alpha);

  if (validate) {
    if(!factorized.areVarsEqual(alpha))
      throw new Error(ERROR_MSG.NON_EQUAL_VARIANCES);
  
    if (!factorized.isNormal())
      throw new Error(ERROR_MSG.NON_NORMAL_DISTRIB);
  }

  const anova = factorized.getOneWayAnova();
  const fCrit = jStat.centralF.inv(1 - alpha, anova.dfBn, anova.dfWn);

  const hypothesis = `THE NULL HYPOTHESIS: the "${factors.name}" factor does not produce a significant difference in the "${values.name}" feature.`;
  const testResult = `Test result: ${(anova.fStat > fCrit) ? 'REJECTED.' : 'FAILED TO REJECT.'}`;

  return getOneWayAnovaDF(anova, alpha, fCrit, hypothesis, testResult);
} // oneWayAnova

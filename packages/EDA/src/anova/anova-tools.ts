// Analysis of Variances (ANOVA): computations

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

import {getNullValue} from '../missing-values-imputation/knn-imputer';

enum ERROR_MSG {
  NON_EQUAL_FACTORS_VALUES_SIZE = 'non-equal sizes of factor and values arrays',
  INCORRECT_SIGNIFICANCE_LEVEL = 'incorrect significance level',
  NON_EQUAL_VARIANCES = 'non-equal variances',
  NON_NORMAL_DISTRIB = 'non-normal distribution',
  UNSUPPORTED_COLUMN_TYPE = 'unsupported feature column type',
  INCORRECT_CATEGORIES_COL_TYPE = 'incorrect categories column type',
  SINGLE_FACTOR = 'single category features',
  CATS_EQUAL_SIZE = 'single value in each category',
  NO_FEATURE_VARIATION = 'no feature variation',
  NO_FEATURE_VARIATION_WITHIN_GROUPS = 'no feature variation within groups',
};

type SampleData = {
  sum: number,
  sumOfSquares: number,
  size: number,
};

/** One-way ANOVA computation results. The classic notations are used (see [2], p. 290). */
type OneWayAnova = {
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

/** One-way ANOVA report */
export type OneWayAnovaReport = {
  anovaTable: OneWayAnova,
  fCritical: number,
  significance: number,
};

/** Categorical column */
type CatCol = DG.Column<DG.COLUMN_TYPE.STRING | DG.COLUMN_TYPE.BOOL>;

/** Numerical column */
type NumCol = DG.Column<DG.COLUMN_TYPE.FLOAT> | DG.Column<DG.COLUMN_TYPE.INT>;

/** Check correctness of significance level. */
export function checkSignificanceLevel(alpha: number) {
  if ((alpha <= 0) || (alpha >= 1))
    throw new Error(ERROR_MSG.INCORRECT_SIGNIFICANCE_LEVEL);
}

/** Compute unbiased variance.*/
export function getVariance(data: SampleData): number {
  // The applied formulas can be found in [4] (see p. 63)
  const size = data.size;

  if (size <= 1)
    return 0;

  return (data.sumOfSquares - (data.sum) ** 2 / size) / (size - 1);
} // getVariance

/** Check equality of variances of 2 samples. F-test is performed.*/
function areVarsEqual(xData: SampleData, yData: SampleData, alpha: number): boolean {
  // The applied approach can be found in [3]
  checkSignificanceLevel(alpha);

  const xVar = getVariance(xData);
  const yVar = getVariance(yData);

  if ((xVar === 0) || (yVar === 0))
    return (xVar === yVar);

  const fStat = xVar / yVar;
  const fCrit = jStat.centralF.inv(1 - alpha, xData.size - 1, yData.size - 1);

  return (fStat < fCrit);
} // areVarsEqual

export class FactorizedData {
  private sums!: Float64Array;
  private sumsOfSquares!: Float64Array;
  private subSampleSizes!: Int32Array;
  private size!: number;
  private catCount!: number;
  private nullsCount = 0;

  constructor(categories: CatCol, values: NumCol, uniqueCount: number) {
    if (categories.length !== values.length)
      throw new Error(ERROR_MSG.NON_EQUAL_FACTORS_VALUES_SIZE);

    if (values.stats.stdev > 0)
      this.setStats(categories, values, uniqueCount);
    else
      throw new Error(ERROR_MSG.NO_FEATURE_VARIATION);
  }

  /** Check equality of variances of factorized data. */
  public areVarsEqual(alpha: number): boolean {
    const K = this.catCount;

    if (K === 1)
      return true;

    const first: SampleData = {sum: this.sums[0], sumOfSquares: this.sumsOfSquares[0], size: this.subSampleSizes[0]};

    for (let i = 1; i < K; ++i) {
      if (!areVarsEqual(first, {sum: this.sums[i], sumOfSquares: this.sumsOfSquares[i],
        size: this.subSampleSizes[i]}, alpha))
        return false;
    }

    return true;
  } // areVarsEqual

  /** Perform one-way ANOVA computations. */
  public getOneWayAnova(): OneWayAnova {
    // Further, notations and formulas from (see [2], p. 290) are used.

    let sum = 0;
    let sumOfSquares = 0;
    let buf = 0;
    let K = this.catCount;
    let nonEmptyCategories = K;

    for (let i = 0; i < K; ++i) {
      if (this.subSampleSizes[i] !== 0) {
        sum += this.sums[i];
        sumOfSquares += this.sumsOfSquares[i];
        buf += this.sums[i] ** 2 / this.subSampleSizes[i];
      } else
        --nonEmptyCategories;
    }

    K = nonEmptyCategories;

    if (K === 1)
      throw new Error(ERROR_MSG.SINGLE_FACTOR);

    const N = this.size - this.nullsCount;
    if (N === K)
      throw new Error(ERROR_MSG.CATS_EQUAL_SIZE);

    const ssTot = sumOfSquares - sum ** 2 / N;
    const ssBn = buf - sum ** 2 / N;
    const ssWn = ssTot - ssBn;

    if (ssWn === 0)
      throw new Error(ERROR_MSG.NO_FEATURE_VARIATION_WITHIN_GROUPS);

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
      pValue: 1 - jStat.centralF.cdf(fStat, dfBn, dfWn),
    };
  } // getOneWayAnova

  /** Compute sum & sums of squares with respect to factor levels. */
  private setStats(categories: CatCol, features: NumCol, uniqueCount: number): void {
    const type = features.type;
    const size = features.length;
    const featuresNull = getNullValue(features);

    switch (type) {
    case DG.COLUMN_TYPE.INT:
    case DG.COLUMN_TYPE.FLOAT:
      const catCount = uniqueCount;
      this.catCount = catCount;
      this.size = size;

      const vals = features.getRawData();
      const cats = categories.getRawData();

      const sums = new Float64Array(catCount).fill(0);
      const sumsOfSquares = new Float64Array(catCount).fill(0);
      const subSampleSizes = new Int32Array(catCount).fill(0);

      let cat: number;

      if (categories.type == DG.COLUMN_TYPE.BOOL) {
        let catIdx = 0;
        let shift = 0;
        let packed = cats[0];
        const MAX_SHIFT = 8 * cats.BYTES_PER_ELEMENT - 1;

        for (let i = 0; i < size; ++i) {
          cat = 1 & (packed >> shift);

          if (vals[i] !== featuresNull) {
            sums[cat] += vals[i];
            sumsOfSquares[cat] += vals[i] ** 2;
            ++subSampleSizes[cat];
          } else
            ++this.nullsCount;


          ++shift;

          if (shift > MAX_SHIFT) {
            shift = 0;
            ++catIdx;
            packed = cats[catIdx];
          }
        }
      } else {
        const categoriesNull = categories.stats.missingValueCount > 0 ? getNullValue(categories) : -1;

        for (let i = 0; i < size; ++i) {
          cat = cats[i];

          if ((cat === categoriesNull) || (vals[i] === featuresNull)) {
            ++this.nullsCount;
            continue;
          }

          sums[cat] += vals[i];
          sumsOfSquares[cat] += vals[i] ** 2;
          ++subSampleSizes[cat];
        }
      }

      this.sums = sums;
      this.sumsOfSquares = sumsOfSquares;
      this.subSampleSizes = subSampleSizes;

      break;

    default:
      throw new Error(ERROR_MSG.UNSUPPORTED_COLUMN_TYPE);
    }
  } // setStats
} // FactorizedData

/** Perform one-way analysis of variances. */
export function oneWayAnova(categores: CatCol, values: NumCol, alpha: number,
  toValidate: boolean = true): OneWayAnovaReport {
  checkSignificanceLevel(alpha);

  const uniqueCount = categores.stats.uniqueCount;

  if (uniqueCount < 2)
    throw new Error(ERROR_MSG.SINGLE_FACTOR);

  const factorized = new FactorizedData(categores, values, uniqueCount);

  if (toValidate) {
    if (!factorized.areVarsEqual(alpha))
      throw new Error(ERROR_MSG.NON_EQUAL_VARIANCES);
  }

  const anova = factorized.getOneWayAnova();

  return {
    anovaTable: anova,
    fCritical: jStat.centralF.inv(1 - alpha, anova.dfBn, anova.dfWn),
    significance: alpha,
  };
} // oneWayAnova

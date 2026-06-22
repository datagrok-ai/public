// Analysis of Variances (ANOVA) - computations

/* REFERENCES

     [1] One-way analysis of variance, https://en.wikipedia.org/wiki/One-way_analysis_of_variance

     [2] G.W. Heiman. Basic Statistics for the Behavioral Sciences, 6th ed. Wadsworth Publishing, 2010

     [3] F-test of equality of variances, https://en.wikipedia.org/wiki/F-test_of_equality_of_variances

     [4] S. McKillup. Statistics Explained, Cambridge University Press, 2005

     [5] B.L. Welch. On the comparison of several mean values: an alternative approach.
         Biometrika 38 (3/4), 1951, pp. 330-336.

     [6] B.P. Welford. Note on a method for calculating corrected sums of squares and products.
         Technometrics 4 (3), 1962, pp. 419-420.

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

/** Per-group statistics produced by the Welford accumulator (see [6]). */
export type GroupStats = {
  /** Sample mean. */
  mean: number,
  /** Sum of squared deviations from the mean (M2). */
  m2: number,
  /** Sample size. */
  size: number,
};

/** One-way Fisher ANOVA computation results. The classic notations are used (see [2], p. 290). */
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

/** One-way Welch ANOVA (W-test) result (see [5]).
 *  No SS/MS decomposition — Welch is not built on a Fisher-style partition. */
export type WelchAnova = {
  /** Welch's W-statistic (the F-distributed test statistic). */
  fStat: number,
  /** Numerator degrees of freedom: k - 1. Integer. */
  dfBn: number,
  /** Welch-Satterthwaite denominator degrees of freedom. Fractional. */
  dfWn: number,
  /** p-value, P(F_{dfBn, dfWn} > fStat). */
  pValue: number,
  /** Per-group sample means (display only). */
  groupMeans: Float64Array,
  /** Per-group unbiased sample variances (display only). */
  groupVariances: Float64Array,
  /** Per-group sample sizes (display only). */
  groupSizes: Int32Array,
};

/** One-way ANOVA report. Discriminated by `method`. */
export type OneWayAnovaReport =
  | {method: 'Fisher', anovaTable: OneWayAnova, fCritical: number, significance: number}
  | {method: 'Welch', anovaTable: WelchAnova, fCritical: number, significance: number};

/** Options for `oneWayAnova`. */
export interface OneWayAnovaOptions {
  /** Method: 'Welch' (default; robust to unequal variances) or 'Fisher' (classical). */
  method?: 'Fisher' | 'Welch',
  /** When true (default), Fisher runs an F-test of equality of variances and throws on failure.
   *  Ignored for Welch — Welch does not require equal variances. */
  toValidate?: boolean,
}

/** Categorical column */
export type CatCol = DG.Column<DG.COLUMN_TYPE.STRING | DG.COLUMN_TYPE.BOOL>;

/** Numerical column */
export type NumCol = DG.Column<DG.COLUMN_TYPE.FLOAT> | DG.Column<DG.COLUMN_TYPE.INT>;

/** Check correctness of significance level. */
export function checkSignificanceLevel(alpha: number) {
  if ((alpha <= 0) || (alpha >= 1))
    throw new Error(ERROR_MSG.INCORRECT_SIGNIFICANCE_LEVEL);
}

/** Unbiased sample variance from a Welford accumulator (see [6]). */
export function getVariance(data: GroupStats): number {
  if (data.size <= 1)
    return 0;

  return data.m2 / (data.size - 1);
} // getVariance

/** Check equality of variances of 2 samples. F-test is performed.*/
export function areVarsEqual(xData: GroupStats, yData: GroupStats, alpha: number): boolean {
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

/**
 * Upper-tail probability P(F > f | df1, df2) for the central F-distribution.
 *
 * Numerically stable in the tail (unlike `1 - jStat.centralF.cdf(f, df1, df2)`,
 * which collapses to 0 once cdf rounds to 1).
 *
 * Identity: P(F > f) = I_x(df2/2, df1/2),  where x = df2 / (df2 + df1 * f),
 * and `jStat.ibeta` retains precision for small x (= large f).
 */
export function fSurvival(f: number, df1: number, df2: number): number {
  if (!Number.isFinite(f))
    return f > 0 ? 0 : 1;
  if (f <= 0)
    return 1;
  const x = df2 / (df2 + df1 * f);
  return jStat.ibeta(x, df2 / 2, df1 / 2);
}

/** Result of factorizing a feature column by a categorical column.
 *  Per-group Welford accumulators are kept as typed arrays (ANOVA hot path). */
export interface Factorization {
  /** Per-group sample means. */
  means: Float64Array;
  /** Per-group sum of squared deviations from the mean (M2). */
  m2: Float64Array;
  /** Per-group sample sizes. */
  sizes: Int32Array;
  /** Number of categories (groups). */
  catCount: number;
  /** Count of skipped rows with a null feature or null category. */
  nullsCount: number;
}

/** Split `values` by `categories` into per-group Welford stats (mean, M2, size).
 *
 *  Single source of truth for the factorization, shared by ANOVA (`FactorizedData`)
 *  and the t-test. Handles INT/FLOAT features and STRING/BOOL categories; rows with a
 *  null feature or null category are skipped (counted in `nullsCount`). */
export function factorize(categories: CatCol, values: NumCol, uniqueCount: number): Factorization {
  const type = values.type;
  const size = values.length;
  const featuresNull = getNullValue(values);

  switch (type) {
  case DG.COLUMN_TYPE.INT:
  case DG.COLUMN_TYPE.FLOAT:
    const catCount = uniqueCount;

    const vals = values.getRawData();
    const cats = categories.getRawData();

    const means = new Float64Array(catCount).fill(0);
    const m2 = new Float64Array(catCount).fill(0);
    const subSampleSizes = new Int32Array(catCount).fill(0);
    let nullsCount = 0;

    let cat: number;

    if (categories.type == DG.COLUMN_TYPE.BOOL) {
      let catIdx = 0;
      let shift = 0;
      let packed = cats[0];
      const MAX_SHIFT = 8 * cats.BYTES_PER_ELEMENT - 1;

      for (let i = 0; i < size; ++i) {
        cat = 1 & (packed >> shift);

        if (vals[i] !== featuresNull) {
          const x = vals[i];
          const n = ++subSampleSizes[cat];
          const delta = x - means[cat];
          means[cat] += delta / n;
          const delta2 = x - means[cat];
          m2[cat] += delta * delta2;
        } else
          ++nullsCount;


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
          ++nullsCount;
          continue;
        }

        const x = vals[i];
        const n = ++subSampleSizes[cat];
        const delta = x - means[cat];
        means[cat] += delta / n;
        const delta2 = x - means[cat];
        m2[cat] += delta * delta2;
      }
    }

    return {means, m2, sizes: subSampleSizes, catCount, nullsCount};

  default:
    throw new Error(ERROR_MSG.UNSUPPORTED_COLUMN_TYPE);
  }
} // factorize

/** Per-group stats for the i-th group of a factorization. */
export function groupStats(f: Factorization, i: number): GroupStats {
  return {mean: f.means[i], m2: f.m2[i], size: f.sizes[i]};
}

export class FactorizedData {
  private means!: Float64Array;
  private m2!: Float64Array;
  private subSampleSizes!: Int32Array;
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

    const first: GroupStats = {mean: this.means[0], m2: this.m2[0], size: this.subSampleSizes[0]};

    for (let i = 1; i < K; ++i) {
      if (!areVarsEqual(first, {mean: this.means[i], m2: this.m2[i],
        size: this.subSampleSizes[i]}, alpha))
        return false;
    }

    return true;
  } // areVarsEqual

  /** Perform one-way Fisher ANOVA computations.
   *  Numerically stable: uses Welford means/M2 instead of the naive Σx² − (Σx)²/N form. */
  public getOneWayAnova(): OneWayAnova {
    let N = 0;
    let nonEmpty = 0;

    for (let i = 0; i < this.catCount; ++i) {
      if (this.subSampleSizes[i] !== 0) {
        N += this.subSampleSizes[i];
        ++nonEmpty;
      }
    }

    const K = nonEmpty;

    if (K === 1)
      throw new Error(ERROR_MSG.SINGLE_FACTOR);
    if (N === K)
      throw new Error(ERROR_MSG.CATS_EQUAL_SIZE);

    // Grand mean weighted by group sizes.
    let grandMeanNumer = 0;
    for (let i = 0; i < this.catCount; ++i) {
      if (this.subSampleSizes[i] !== 0)
        grandMeanNumer += this.subSampleSizes[i] * this.means[i];
    }
    const grandMean = grandMeanNumer / N;

    // ssWn = Σ M2[k] — Welford identity for within-group SS.
    // ssBn = Σ n[k] * (mean[k] − grandMean)² — between-group SS.
    // ssTot = ssBn + ssWn — exact identity, no independent recomputation.
    let ssWn = 0;
    let ssBn = 0;
    for (let i = 0; i < this.catCount; ++i) {
      const n = this.subSampleSizes[i];
      if (n === 0) continue;
      ssWn += this.m2[i];
      const diff = this.means[i] - grandMean;
      ssBn += n * diff * diff;
    }
    const ssTot = ssBn + ssWn;

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
      pValue: fSurvival(fStat, dfBn, dfWn),
    };
  } // getOneWayAnova

  /** Perform one-way Welch ANOVA (W-test, see [5]).
   *  Robust to unequal variances across groups. */
  public getWelchAnova(): WelchAnova {
    const groupMeans: number[] = [];
    const groupVars: number[] = [];
    const groupSizes: number[] = [];

    for (let i = 0; i < this.catCount; ++i) {
      const n = this.subSampleSizes[i];
      if (n < 2) continue;
      const v = this.m2[i] / (n - 1);
      if (v <= 0)
        throw new Error(ERROR_MSG.NO_FEATURE_VARIATION_WITHIN_GROUPS);
      groupMeans.push(this.means[i]);
      groupVars.push(v);
      groupSizes.push(n);
    }

    const K = groupMeans.length;
    if (K < 2)
      throw new Error(ERROR_MSG.SINGLE_FACTOR);

    // Welch (1951): weights w_i = n_i / s_i², weighted grand mean, numerator and lambda.
    let W = 0;
    let wxSum = 0;
    for (let i = 0; i < K; ++i) {
      const w = groupSizes[i] / groupVars[i];
      W += w;
      wxSum += w * groupMeans[i];
    }
    const wGrandMean = wxSum / W;

    let numer = 0;
    let lambda = 0;
    for (let i = 0; i < K; ++i) {
      const w = groupSizes[i] / groupVars[i];
      const diff = groupMeans[i] - wGrandMean;
      numer += w * diff * diff;
      const r = 1 - w / W;
      lambda += (r * r) / (groupSizes[i] - 1);
    }
    numer /= (K - 1);
    const denom = 1 + 2 * (K - 2) / (K * K - 1) * lambda;

    const fStat = numer / denom;
    const dfBn = K - 1;
    const dfWn = (K * K - 1) / (3 * lambda);

    return {
      fStat: fStat,
      dfBn: dfBn,
      dfWn: dfWn,
      pValue: fSurvival(fStat, dfBn, dfWn),
      groupMeans: Float64Array.from(groupMeans),
      groupVariances: Float64Array.from(groupVars),
      groupSizes: Int32Array.from(groupSizes),
    };
  } // getWelchAnova

  /** Compute per-group mean and M2 (sum of squared deviations from the mean) via Welford (see [6]).
   *  Thin wrapper over the shared `factorize` — single source of truth. */
  private setStats(categories: CatCol, features: NumCol, uniqueCount: number): void {
    const f = factorize(categories, features, uniqueCount);
    this.means = f.means;
    this.m2 = f.m2;
    this.subSampleSizes = f.sizes;
    this.catCount = f.catCount;
    this.nullsCount = f.nullsCount;
  } // setStats
} // FactorizedData

/** Perform one-way analysis of variances.
 *  Default method is Welch (robust to unequal variances).
 *  Pass `{method: 'Fisher'}` for the classical equal-variance ANOVA. */
export function oneWayAnova(categores: CatCol, values: NumCol, alpha: number,
  opts: OneWayAnovaOptions = {}): OneWayAnovaReport {
  checkSignificanceLevel(alpha);

  const method = opts.method ?? 'Welch';
  const toValidate = opts.toValidate ?? true;

  const uniqueCount = categores.stats.uniqueCount;

  if (uniqueCount < 2)
    throw new Error(ERROR_MSG.SINGLE_FACTOR);

  const factorized = new FactorizedData(categores, values, uniqueCount);

  if (method === 'Fisher') {
    if (toValidate && !factorized.areVarsEqual(alpha))
      throw new Error(ERROR_MSG.NON_EQUAL_VARIANCES);

    const anova = factorized.getOneWayAnova();

    return {
      method: 'Fisher',
      anovaTable: anova,
      fCritical: jStat.centralF.inv(1 - alpha, anova.dfBn, anova.dfWn),
      significance: alpha,
    };
  } else {
    const anova = factorized.getWelchAnova();

    return {
      method: 'Welch',
      anovaTable: anova,
      fCritical: jStat.centralF.inv(1 - alpha, anova.dfBn, anova.dfWn),
      significance: alpha,
    };
  }
} // oneWayAnova

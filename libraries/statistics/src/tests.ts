//@ts-ignore: no types
import * as jStat from 'jstat';

type testStats = {
  'p-value': number,
  'Mean difference'?: number,
  'Median difference'?: number,
  'p-value more': number,
  'p-value less': number,
};

type Population = number[] | Float32Array | Int32Array;

export function tTest(sample1: Population, sample2: Population, devKnown=false, devEqual=false): testStats {
  if (sample1.length <= 1 || sample2.length <= 1)
    throw new Error(`StatisticsError: Wrong sample size; expected at least 2, got ${Math.min(sample1.length, sample2.length)})`);

  const mean1: number = jStat.mean(sample1);
  const mean2: number = jStat.mean(sample2);
  const variance1: number = jStat.variance(sample1);
  const variance2: number = jStat.variance(sample2);
  const length1 = sample1.length;
  const length2 = sample2.length;

  let pMore;
  let pLess;
  let pTot;

  if (!devKnown) {
    if (!devEqual) {
      const sampleVariance1 = variance1 / length1;
      const sampleVariance2 = variance2 / length2;
      const criticalValue = (mean1 - mean2) / Math.sqrt(sampleVariance1 + sampleVariance2);
      const dof = Math.pow((sampleVariance1 + sampleVariance2), 2) /
        (Math.pow(sampleVariance1, 2) / (length1 - 1) + Math.pow(sampleVariance2, 2) / (length2 - 1));

      pLess = jStat.studentt.cdf(criticalValue, dof);
      pMore = 1 - pLess;
      pTot = 2 * (pLess < pMore ? pLess : pMore);
    } else {
      const dof = length1 + length2 - 2;
      const totalVariance = (variance1 * (length1 - 1) + variance2 * (length2 - 1)) / dof;
      const criticalValue = Math.sqrt(length1 * length2 / (length1 + length2)) * (mean1 - mean2) / totalVariance;

      pMore = 1 - jStat.studentt.cdf(criticalValue, dof);
      pLess = jStat.studentt.cdf(criticalValue, dof);
      pTot = 2 * (pLess < pMore ? pLess : pMore);
    }
  } else {
    const sampleVariance1 = variance1 / length1;
    const sampleVariance2 = variance2 / length2;
    const criticalValue = (mean1 - mean2) / Math.sqrt(sampleVariance1 + sampleVariance2);

    pLess = jStat.normal.pdf(criticalValue, 0, 1);
    pMore = 1 - pLess;
    pTot = 2 * (pLess < pMore ? pLess : pMore);
  }
  return {'p-value': pTot, 'Mean difference': mean1 - mean2, 'p-value more': pMore, 'p-value less': pLess};
}

export function uTest(x: number[], y: number[], continuity=true): testStats {
  const xy = x.concat(y);
  const n1 = x.length;
  const n2 = y.length;
  const med1 = jStat.median(x);
  const med2 = jStat.median(y);

  const ranks = jStat.rank(xy);

  const R1 = jStat.sum(ranks.slice(0, n1));
  const U1 = R1 - n1 * (n1 + 1) / 2;
  const U2 = n1 * n2 - U1;
  const U = U1 > U2 ? U1 : U2;

  const mu = n1 * n2 / 2;
  const n = n1 + n2;

  const tieTerm = _tieTerm(ranks);
  const s = Math.sqrt(n1 * n2 / 12 * ((n + 1) - tieTerm / (n* (n - 1))));

  let numerator = U - mu;

  if (continuity)
    numerator -= 0.5;


  const z = numerator / s;

  const p = 2 * (1 - jStat.normal.cdf(z, 0, 1));

  return {'p-value': p, 'Median difference': med1 - med2, 'p-value more': p, 'p-value less': p};
}

function _tieTerm(ranks: number[]): number {
  const ties: {[key: number]: number} = {};

  ranks.forEach((num) => {
    ties[num] = (ties[num] || 0) + 1;
  });

  return jStat.sum(Object.values(ties));
}

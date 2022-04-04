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

export function tTest(arr1: Population, arr2: Population, devKnown=false, devEqual=false): testStats {
  const m1: number = jStat.mean(arr1);
  const m2: number = jStat.mean(arr2);
  const v1: number = jStat.variance(arr1);
  const v2: number = jStat.variance(arr2);
  const n1 = arr1.length;
  const n2 = arr2.length;

  let wv1;
  let wv2;
  let wv;
  let Z;
  let K;
  let pMore;
  let pLess;
  let pTot;

  if (!devKnown) {
    if (!devEqual) {
      wv1 = v1 / n1;
      wv2 = v2 / n2;
      Z = (m1 - m2) / Math.sqrt(wv1 + wv2);
      K = Math.pow((wv1 + wv2), 2) / (wv1 * wv1 / (n1 - 1) + wv2 * wv2 / (n2 - 1));

      pLess = jStat.studentt.cdf(Z, K);
      pMore = 1 - pLess;
      pTot = 2 * (pLess < pMore ? pLess : pMore);
    } else {
      K = n1 + n2 - 2;
      wv = (v1 * (n1 - 1) + v2 * (n2 - 1)) / K;
      Z = Math.sqrt(n1 * n2 / (n1 + n2)) * (m1 - m2) / wv;

      pMore = 1 - jStat.studentt.cdf(Z, K);
      pLess = jStat.studentt.cdf(Z, K);
      pTot = 2 * (pLess < pMore ? pLess : pMore);
    }
  } else {
    wv1 = v1 / n1;
    wv2 = v2 / n2;
    Z = (m1 - m2) / Math.sqrt(wv1 + wv2);

    pLess = jStat.normal.pdf(Z, 0, 1);
    pMore = 1 - pLess;
    pTot = 2 * (pLess < pMore ? pLess : pMore);
  }
  return {'p-value': pTot, 'Mean difference': m1 - m2, 'p-value more': pMore, 'p-value less': pLess};
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

  if (continuity) {
    numerator -= 0.5;
  }

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

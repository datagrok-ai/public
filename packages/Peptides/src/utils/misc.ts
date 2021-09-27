//@ts-ignore: no types
import * as jStat from 'jstat';
import binarySearch from 'binary-search';

type testStats = {'p-value': number, 't-score'?: number, 'u-score'?: number};

export function decimalAdjust(type: 'floor' | 'ceil' | 'round', value: number, exp: number): number {
  // If the exp is undefined or zero...
  if (typeof exp === 'undefined' || +exp === 0) {
    return Math[type](value);
  }
  value = +value;
  exp = +exp;
  // If the value is not a number or the exp is not an integer...
  if (isNaN(value) || !(typeof exp === 'number' && exp % 1 === 0)) {
    return NaN;
  }
  // Shift
  let valueArr = value.toString().split('e');
  value = Math[type](+(valueArr[0] + 'e' + (valueArr[1] ? (+valueArr[1] - exp) : -exp)));
  // Shift back
  valueArr = value.toString().split('e');
  return +(valueArr[0] + 'e' + (valueArr[1] ? (+valueArr[1] + exp) : exp));
}

export function tTest(arr1: number[], arr2: number[], alpha=0.05, devKnown=false, devEqual=false): testStats {
  const m1: number = jStat.mean(arr1);
  const m2: number = jStat.mean(arr2);
  const v1: number = jStat.variance(arr1);
  const v2: number = jStat.variance(arr2);
  const n1 = arr1.length;
  const n2 = arr2.length;

  let wv1, wv2, wv,
    Z, K,
    p_more, p_less, p_tot,
    Vk_more_ts, Vk_less_ts, Vk_more, Vk_less;

  if (!devKnown) {
    if (!devEqual) {
      wv1 = v1 / n1;
      wv2 = v2 / n2;
      Z = (m1 - m2) / Math.sqrt(wv1 + wv2);
      K = Math.pow((wv1 + wv2), 2) / (wv1 * wv1 / (n1 - 1) + wv2 * wv2 / (n2 - 1));

      p_less = jStat.studentt.cdf(Z, K);
      p_more = 1 - p_less;
      p_tot = 2 * (p_less < p_more ? p_less : p_more);
    } else {
      K = n1 + n2 - 2;
      wv = (v1 * (n1 - 1) + v2 * (n2 - 1)) / K;
      Z = Math.sqrt(n1 * n2 / (n1 + n2)) * (m1 - m2) / wv;

      p_more = 1 - jStat.studentt.cdf(Z, K);
      p_less = jStat.studentt.cdf(Z, K);
      p_tot =  2 * (p_less < p_more ? p_less : p_more);
    }
    Vk_more_ts = jStat.studentt.inv(1 - alpha / 2, K);
    Vk_less_ts = jStat.studentt.inv(alpha / 2, K);
    Vk_more = jStat.studentt.inv(1 - alpha, K);
    Vk_less = jStat.studentt.inv(alpha, K);
  } else {
    wv1 = v1 / n1;
    wv2 = v2 / n2;
    Z = (m1 - m2) / Math.sqrt(wv1 + wv2);

    p_less = jStat.normal.pdf(Z, 0, 1);
    p_more = 1 - p_less;
    p_tot = 2 * (p_less < p_more ? p_less : p_more);

    Vk_more_ts = jStat.normal.inv(1 - alpha / 2);
    Vk_less_ts = jStat.normal.inv(alpha / 2);
    Vk_more = jStat.normal.inv(1 - alpha);
    Vk_less = jStat.normal.inv(alpha);
  }
  return {'p-value': p_tot, 't-score': Z};
}

export function uTest(x1: number[], x2: number[]): testStats {
  const ranks = ranking([x1, x2]);
  const n1 = x1.length;
  const n2 = x2.length;
  const U = n1 * n2;

  // simple U-test
  const u1 = x1.reduce((pv, cv) => pv + cv) - (n1 * (n1 + 1)) / 2;
  const u2 = U - u1;

  // Rank-biserial correlation
  // const rankBiserial = 1 - (2 * u2) / U;

  // Stabdard Deviation and Absolute Value
  const T = tieCorrection(ranks);
  let sd = Math.sqrt((T * U * (n1 + n2 + 1)) / 12.0);
  const mRank = 0.5 + U / 2;
  const absolutValue = (Math.max(u1, u2) - mRank) / sd;

  // Effect strength
  // const effectStrength = absolutValue / Math.sqrt(n1 + n2); // effect strength

  //p value
  const p = jStat.normal.cdf(-Math.abs(((u1 > u2 ? u1 : u2) - mRank) / sd), 0, 1) * 2;

  // return {
  //   u1,
  //   u2,
  //   p,
  //   effectStrength,
  //   rankBiserial,
  //   absolutValue,
  // };
  return {'p-value': p, 'u-score': absolutValue};
}

function ranking(arrays: number[][]) {
  let concatArray: number[] = [];
  arrays.forEach((item) => {
    concatArray = concatArray.concat(item);
  });

  const sorted = concatArray.slice().sort((a, b) => b - a);
  const ranks = concatArray.map(
    (value) =>
      binarySearch(sorted, value, function (element: number, needle: number) {
        return needle - element;
      }) + 1,
  );

  return ranks;
}

function tieCorrection(rankValues: number[]) {
  if (rankValues.length === 0) {
    throw new Error('tieCorrection: array length should be greater than 0');
  }
  if (rankValues.length === 1) {
    return 1;
  }

  const sortedArr = rankValues.slice().sort((a, b) => b - a);
  const leftArr = sortedArr.slice(1, sortedArr.length);
  const rightArr = sortedArr.slice(0, sortedArr.length - 1);
  const nonNegative = [0, leftArr.join('').localeCompare(rightArr.join('')), 0];

  const nonNegativeIdxs = nonNegative
    .map((a, i) => (a === 0 ? i : -1))
    .filter((a) => a !== -1);

  let diffCounter = nonNegativeIdxs
    .slice(1, nonNegativeIdxs.length)
    .map(function (num, idx) {
      return num - nonNegativeIdxs.slice(0, nonNegativeIdxs.length - 1)[idx];
    })
    .length;

  return (
    1 -
    (Math.pow(diffCounter, 3) - diffCounter) /
      (Math.pow(sortedArr.length, 3) - sortedArr.length)
  );
}

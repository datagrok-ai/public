import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {BitArrayMetricsNames} from './typed-metrics/consts';
import {MmDistanceFunctionsNames} from './macromolecule-distance-functions';

export const similarityMetric: { [name: string]: (x: BitArray, y: BitArray) => number } = {
  [BitArrayMetricsNames.Tanimoto]: tanimotoSimilarity,
  [BitArrayMetricsNames.Dice]: diceSimilarity,
  [BitArrayMetricsNames.Asymmetric]: asymmetricSimilarity,
  [BitArrayMetricsNames.BraunBlanquet]: braunBlanquetSimilarity,
  [BitArrayMetricsNames.Cosine]: cosineSimilarity,
  [BitArrayMetricsNames.Kulczynski]: kulczynskiSimilarity,
  [BitArrayMetricsNames.McConnaughey]: mcConnaugheySimilarity,
  [BitArrayMetricsNames.RogotGoldberg]: rogotGoldbergSimilarity,
  [BitArrayMetricsNames.Russel]: russelSimilarity,
  [BitArrayMetricsNames.Sokal]: sokalSimilarity,
  [BitArrayMetricsNames.Hamming]: hammingSimilarity,
  [BitArrayMetricsNames.Euclidean]: euclideanSimilarity,
};

export const distanceMetrics: { [name: string]: (x: BitArray, y: BitArray) => number } = {
  [BitArrayMetricsNames.Tanimoto]: tanimotoDistance,
  [BitArrayMetricsNames.Dice]: diceDistance,
  [BitArrayMetricsNames.Asymmetric]: asymmetricDistance,
  [BitArrayMetricsNames.BraunBlanquet]: braunBlanquetDistance,
  [BitArrayMetricsNames.Cosine]: cosineDistance,
  [BitArrayMetricsNames.Kulczynski]: kulczynskiDistance,
  [BitArrayMetricsNames.McConnaughey]: mcConnaugheyDistance,
  [BitArrayMetricsNames.RogotGoldberg]: rogotGoldbergDistance,
  [BitArrayMetricsNames.Russel]: russelDistance,
  [BitArrayMetricsNames.Sokal]: sokalDistance,
  [BitArrayMetricsNames.Hamming]: hammingDistance,
  [BitArrayMetricsNames.Euclidean]: euclideanDistanceBitArray,
};

export const CHEM_SIMILARITY_METRICS = [
  BitArrayMetricsNames.Tanimoto,
  BitArrayMetricsNames.Dice,
  BitArrayMetricsNames.Cosine];
export const SEQ_SPACE_SIMILARITY_METRICS = [
  BitArrayMetricsNames.Tanimoto,
  BitArrayMetricsNames.Asymmetric,
  BitArrayMetricsNames.Cosine,
  BitArrayMetricsNames.Sokal];
export const MACROMOLECULE_SIMILARITY_METRICS = [
  MmDistanceFunctionsNames.HAMMING,
  MmDistanceFunctionsNames.LEVENSHTEIN,
  MmDistanceFunctionsNames.MONOMER_CHEMICAL_DISTANCE,
  MmDistanceFunctionsNames.NEEDLEMANN_WUNSCH
];


export function tanimotoSimilarity(x: BitArray, y: BitArray): number {
  const total = x.trueCount() + y.trueCount();
  if (total == 0) return 1.0;
  const common = x.andWithCountBits(y, true);
  return common / (total - common);
}

export function tanimotoDistance(x: BitArray, y: BitArray): number {
  return 1 - tanimotoSimilarity(x, y);
}

export function tanimotoDistanceIntArray(x: Uint32Array, y: Uint32Array): number {
  const xb = new BitArray(x, x.length * 32);
  const yb = new BitArray(y, y.length * 32);
  return getDistanceFromSimilarity(tanimotoSimilarity(xb, yb));
}

export function vectorEuclideanDistance(x: ArrayLike<number>, y: ArrayLike<number>): number {
  let sum = 0;
  for (let i = 0; i < x.length; i++)
    sum += Math.pow(x[i] - y[i], 2);
  return Math.sqrt(sum);
}

export function vectorManhattenDistance(x: ArrayLike<number>, y: ArrayLike<number>): number {
  let sum = 0;
  for (let i = 0; i < x.length; i++)
    sum += Math.abs(x[i] - y[i]);
  return sum;
}

export function vectorCosineDistance(x: ArrayLike<number>, y: ArrayLike<number>): number {
  let multSum = 0; let xSquareSum = 0; let ySquareSum = 0;
  for (let i = 0; i < x.length; i++) {
    multSum += x[i] * y[i];
    xSquareSum += x[i] * x[i];
    ySquareSum += y[i] * y[i];
  }
  const sim = multSum / (Math.sqrt(xSquareSum) * Math.sqrt(ySquareSum));
  // similarity is in range [-1, 1], but we need distance in range [0, 1]
  return (1 - sim) / 2; 

}

export function diceSimilarity(x: BitArray, y: BitArray): number {
  const total = x.trueCount() + y.trueCount();
  if (total == 0) return 0.0;
  const common = x.andWithCountBits(y, true);
  return 2 * common / total;
}

export function diceDistance(x: BitArray, y: BitArray): number {
  return 1 - diceSimilarity(x, y);
}

export function cosineSimilarity(x: BitArray, y: BitArray): number {
  const total = x.trueCount() * y.trueCount();
  if (total == 0) return 0.0;
  const common = x.andWithCountBits(y, true);
  return common / Math.sqrt(total);
}

export function cosineDistance(x: BitArray, y: BitArray): number {
  return 1 - cosineSimilarity(x, y);
}

export function euclideanSimilarity(x: BitArray, y: BitArray): number {
  return getSimilarityFromDistance(euclideanDistanceBitArray(x, y));
}

export function euclideanDistanceBitArray(x: BitArray, y: BitArray): number {
  return Math.sqrt(x.trueCount() + y.trueCount() - 2 * x.andWithCountBits(y, true));
}

export function hammingSimilarity(x: BitArray, y: BitArray): number {
  return getSimilarityFromDistance(hammingDistance(x, y));
}

export function hammingDistance(x: BitArray, y: BitArray): number {
  return x.trueCount() + y.trueCount() - 2 * x.andWithCountBits(y, true);
}

export function sokalSimilarity(x: BitArray, y: BitArray): number {
  const total = x.trueCount() + y.trueCount();
  const common = x.andWithCountBits(y, true);
  return common / (2 * total - 3 * common);
}

export function sokalDistance(x: BitArray, y: BitArray): number {
  return 1 - sokalSimilarity(x, y);
}

export function kulczynskiSimilarity(x: BitArray, y: BitArray): number {
  const total = x.trueCount() + y.trueCount();
  const totalProd = x.trueCount() * y.trueCount();
  if (totalProd == 0) return 0.0;
  const common = x.andWithCountBits(y, true);
  return (common * total) / (2 * totalProd);
}

export function kulczynskiDistance(x: BitArray, y: BitArray): number {
  return getDistanceFromSimilarity(kulczynskiSimilarity(x, y));
}

export function mcConnaugheySimilarity(x: BitArray, y: BitArray): number {
  const total = x.trueCount() + y.trueCount();
  const totalProd = x.trueCount() * y.trueCount();
  if (totalProd == 0) return 0.0;
  const common = x.andWithCountBits(y, true);
  return (common * total - totalProd) / totalProd;
}

export function mcConnaugheyDistance(x: BitArray, y: BitArray): number {
  return getDistanceFromSimilarity(mcConnaugheySimilarity(x, y));
}

export function asymmetricSimilarity(x: BitArray, y: BitArray): number {
  const min = Math.min(x.trueCount(), y.trueCount());
  if (min == 0) return 0.0;
  const common = x.andWithCountBits(y, true);
  return common / min;
}

export function asymmetricDistance(x: BitArray, y: BitArray): number {
  return 1 - asymmetricSimilarity(x, y);
}

export function braunBlanquetSimilarity(x: BitArray, y: BitArray): number {
  const max = Math.max(x.trueCount(), y.trueCount());
  if (max == 0) return 0.0;
  const common = x.andWithCountBits(y, true);
  return common / max;
}

export function braunBlanquetDistance(x: BitArray, y: BitArray): number {
  return getDistanceFromSimilarity(braunBlanquetSimilarity(x, y));
}

export function russelSimilarity(x: BitArray, y: BitArray): number {
  if (x.length == 0) return 0.0;
  const common = x.andWithCountBits(y, true);
  return common / x.length;
}

export function russelDistance(x: BitArray, y: BitArray): number {
  return getDistanceFromSimilarity(russelSimilarity(x, y));
}

export function rogotGoldbergSimilarity(x: BitArray, y: BitArray): number {
  const common = x.andWithCountBits(y, true);
  const total = x.countBits(true) + y.countBits(true);
  const len = x.length;
  const diff = len - total + common;
  if ((common == len) || (diff == len)) return 1.0;
  else return common / total + diff / (2 * len - total);
}

export function rogotGoldbergDistance(x: BitArray, y: BitArray): number {
  return getDistanceFromSimilarity(rogotGoldbergSimilarity(x, y));
}

export function getSimilarityFromDistance(distance: number) {
  return 1 / (1 + distance);
}

export function getDistanceFromSimilarity(similarity: number) { //in case similarity is 0, use max number for float32
  return similarity <= 0 ? 3.402823E+38 : (1 / similarity) - 1;
}

export function numericDistance(args?: {range?: number}) {
  if (args && args.range != undefined && args.range > 0) {
    const range = args.range;
    return (a: number, b: number) => Math.abs(a - b) / range;
  }

  return (a: number, b: number) => Math.abs(a - b);
}

export function commonItemsCount(args?: {mostCommon?: Set<number>}) {
  const mostCommon = args?.mostCommon ?? new Set<number>();
  return (arr1: ArrayLike<number>, arr2: ArrayLike<number>) => {
    const len1 = arr1.length;
    const len2 = arr2.length;
    let count = 0;
    let i1 = 0;
    let i2 = 0;

    while ((i1 < len1) && (i2 < len2)) {
      if (arr1[i1] === arr2[i2]) {
        if (!mostCommon?.has(arr1[i1]))
          ++count;
        ++i1;
        ++i2;
      } else if (arr1[i1] < arr2[i2]) { ++i1; } else { ++i2; }
    }

    return count;
  };
}

export function inverseCommonItemsCount(args?: {mostCommon?: Set<number>}) {
  const f = commonItemsCount(args);
  return (arr1: ArrayLike<number>, arr2: ArrayLike<number>) => {
    if (arr2.length === 0 || arr1.length === 0)
      return 10000;

    return Math.min(arr1.length, arr2.length) / (f(arr1, arr2) + 0.0001);
  };
}

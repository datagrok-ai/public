/* eslint-disable valid-jsdoc */

const MIN_DIM = 1;

/** Return clipped value */
export function clip(val: number, min: number, max: number): number {
  return (val < min) ? min : (val > max) ? max : val;
}

/** Return the Euclidean distance between vectors */
export function euclideanDistance(w1: Float32Array, w2: Float32Array) {
  let sum = 0;
  const dim = w1.length;

  for (let i = 0; i < dim; ++i)
    sum += (w1[i] - w2[i])**2;

  return Math.sqrt(sum);
}

/** Return two different items randomly from array */
export function pickTwo(arr: any[]) {
  const i = Math.floor(Math.random() * arr.length);
  let j;
  do
    j = Math.floor(Math.random() * arr.length);
  while (j === i);
  return [arr[i], arr[j]];
}

export function getFloatArrays(count: number, length: number): Float32Array[] {
  const arrs = new Array<Float32Array>(count);
  for (let i = 0; i < count; ++i)
    arrs[i] = new Float32Array(length);
  return arrs;
}

export function isDimValid(dim: number): boolean {
  if (!Number.isInteger(dim))
    return false;

  return dim >= MIN_DIM;
}

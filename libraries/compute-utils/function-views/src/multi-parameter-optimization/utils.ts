export function clip(val: number, min: number, max: number): number {
  return (val < min) ? min : (val > max) ? max : val;
}

export function euclideanDistance(w1: Float32Array, w2: Float32Array) {
  let sum = 0;
  const dim = w1.length;

  for (let i = 0; i < dim; ++i)
    sum += (w1[i] - w2[i])**2;

  return Math.sqrt(sum);
}

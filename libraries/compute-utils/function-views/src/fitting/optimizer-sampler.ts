function sampleUniform(samplesCount: number, top: number, bottom: number): number[] {
  const scale = top - bottom;
  const sample = new Array<number>(samplesCount);

  for (let i = 0; i < samplesCount; i ++) {
    const r = Math.random();
    sample[i] = bottom + r*scale;
  }

  return sample;
}

export function sampleParams(samplesCount: number, top: Float32Array, bottom: Float32Array): Float32Array[] {
  const dim = top.length;
  const params = new Array<Float32Array>(samplesCount);
  for (let i = 0; i < samplesCount; i ++)
    params[i] = new Float32Array(dim);

  for (let i = 0; i < dim; i ++) {
    if (top[i] === bottom[i]) {
      for (let j = 0; j < samplesCount; j ++)
        params[j][i] = top[i];
    } else {
      const paramVariations = sampleUniform(samplesCount, top[i], bottom[i]);
      for (let j = 0; j < samplesCount; j ++)
        params[j][i] = paramVariations[j];
    }
  }

  return params;
}

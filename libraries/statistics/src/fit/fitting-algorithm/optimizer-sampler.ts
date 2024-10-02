import seedRandom from 'seedrandom';

function sampleUniform(samplesCount: number, top: number, bottom: number, randomFn: seedRandom.PRNG): number[] {
  const scale = top - bottom;
  const sample = new Array<number>(samplesCount);

  for (let i = 0; i < samplesCount; i ++) {
    const r = randomFn();
    sample[i] = bottom + r * scale;
  }

  return sample;
}

export function sampleParams(samplesCount: number, top: Float32Array, bottom: Float32Array): Float32Array[] {
  const randomFn = seedRandom('12345');
  const dim = top.length;
  const params = new Array<Float32Array>(samplesCount);
  for (let i = 0; i < samplesCount; i ++)
    params[i] = new Float32Array(dim);

  for (let i = 0; i < dim; i ++) {
    if (top[i] === bottom[i]) {
      for (let j = 0; j < samplesCount; j ++)
        params[j][i] = top[i];
    }
    else {
      const paramVariations = sampleUniform(samplesCount, top[i], bottom[i], randomFn);
      for (let j = 0; j < samplesCount; j ++)
        params[j][i] = paramVariations[j];
    }
  }

  return params;
}

import {runFormula} from './formulas-resolver';
import {ChangingValue, ConstValue, ValueBoundsData} from './optimizer-misc';

function sampleUniform(samplesCount: number, top: number, bottom: number,
  rand: () => number = Math.random): number[] {
  const scale = top - bottom;
  const sample = new Array<number>(samplesCount);

  for (let i = 0; i < samplesCount; i ++) {
    const r =rand();
    sample[i] = bottom + r*scale;
  }

  return sample;
}

export function sampleParams(samplesCount: number, top: Float64Array, bottom: Float64Array,
  rand: () => number = Math.random): Float64Array[] {
  const dim = top.length;
  const params = new Array<Float64Array>(samplesCount);
  for (let i = 0; i < samplesCount; i ++)
    params[i] = new Float64Array(dim);

  for (let i = 0; i < dim; i ++) {
    if (top[i] === bottom[i]) {
      for (let j = 0; j < samplesCount; j ++)
        params[j][i] = top[i];
    } else {
      const paramVariations = sampleUniform(samplesCount, top[i], bottom[i], rand);
      for (let j = 0; j < samplesCount; j ++)
        params[j][i] = paramVariations[j];
    }
  }

  return params;
}


type AccData = {
  constValues: [string, ConstValue][],
  nonFormulaBounds: [string, number, ChangingValue][],
  formulaBounds: [string, number, ChangingValue][],
  boundsIdx: number,
};

export type SamplerInputsConfig = {
  inputs: Record<string, ValueBoundsData>,
}

export function sampleParamsWithFormulaBounds(
  samplesCount: number,
  inputs: Record<string, ValueBoundsData>,
  rand: () => number = Math.random,
): [Float64Array[], Float64Array[], Float64Array[]] {
  // partition inputs based on the type: fixed, numeric or formula bounds,
  // keeping original order index of varried inputs
  const {constValues, nonFormulaBounds, formulaBounds} = Object.entries(inputs).reduce((acc, [name, val]) => {
    if (val.type === 'const')
       acc.constValues.push([name, val]);
    else if (val.top.type === 'value' && val.bottom.type === 'value') {
      acc.nonFormulaBounds.push([name, acc.boundsIdx, val]);
      acc.boundsIdx++;
    } else {
      acc.formulaBounds.push([name, acc.boundsIdx, val]);
      acc.boundsIdx++;
    }
    return acc;
  }, {constValues: [], nonFormulaBounds: [], formulaBounds: [], boundsIdx: 0} as AccData);

  const variedInputsCount = nonFormulaBounds.length + formulaBounds.length;

  const params = new Array<Float64Array>(samplesCount).fill(null as any).map(() => new Float64Array(variedInputsCount));
  const bottoms = new Array<Float64Array>(samplesCount).fill(null as any).map(() => new Float64Array(variedInputsCount));
  const tops = new Array<Float64Array>(samplesCount).fill(null as any).map(() => new Float64Array(variedInputsCount));

  const contextFixed: Record<string, any> = {};
  // add fixes inputs (can be any type) to formula context
  for (const [name, data] of constValues) {
    contextFixed[name] = data.value;
  }

  // iteration over samples
  for (let nsample = 0; nsample < samplesCount; nsample++) {
    // reset context for each sample to fixed inputs only
    const context = {...contextFixed};
    // iteration over a single sample point varied dimensions
    // non-formula inputs first, so formulas can refer to already picked values
    for (const [name, boundsIdx, bound] of [...nonFormulaBounds, ...formulaBounds]) {
      const min = bound.bottom.type === 'value' ? bound.bottom.value : runFormula(bound.bottom.formula, context);
      const max = bound.top.type === 'value' ? bound.top.value : runFormula(bound.top.formula, context);
      if (min == null || max == null || min > max) {
        const msg = `Sampling staring point failed for input ${name}, min: ${min}, max: ${max}`;
        console.error(msg);
        throw new Error(msg);
      }
      // pick value and add to context
      const r = rand();
      const scale = max - min;
      const val = min + r*scale;
      params[nsample][boundsIdx] = val;
      context[name] = val;
      // add bounds
      bottoms[nsample][boundsIdx] = min;
      tops[nsample][boundsIdx] = max;
    }
  }

  return [params, bottoms, tops] as const;
}

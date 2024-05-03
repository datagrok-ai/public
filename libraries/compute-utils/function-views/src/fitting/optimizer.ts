/** Monte Carlo optimizer */

import {Extremum, OptimizationResult} from './optimizer-misc';
import {optimizeNM, NelderMeadSettings} from './optimizer-nelder-mead';
import {sampleParams} from './optimizer-sampler';

export async function performNelderMeadOptimization(
  objectiveFunc: (x: Float32Array) => Promise<number>,
  paramsBottom: Float32Array,
  paramsTop: Float32Array,
  settings: NelderMeadSettings,
  samplesCount: number = 1,
): Promise<OptimizationResult> {
  const params = sampleParams(samplesCount, paramsTop, paramsBottom);

  const extremums: Extremum[] = [];
  const warnings: string[] = [];

  for (let i = 0; i < samplesCount; ++i) {
    try {
      extremums.push(await optimizeNM(objectiveFunc, params[i], settings, paramsBottom, paramsTop));
    } catch (e) {
      warnings.push((e instanceof Error) ? e.message : 'Platform issue');
    }
  }

  return {
    extremums: extremums,
    warnings: warnings,
  };
}

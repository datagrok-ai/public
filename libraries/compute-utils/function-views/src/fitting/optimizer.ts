/** the Nelder-Mead optimizer */
import * as DG from 'datagrok-api/dg';

import {Extremum, OptimizationResult, InconsistentTables, sleep} from './optimizer-misc';
import {optimizeNM} from './optimizer-nelder-mead';
import {sampleParams} from './optimizer-sampler';
import {MS_TO_SLEEP} from './constants';

export async function performNelderMeadOptimization(
  objectiveFunc: (x: Float32Array) => Promise<number>,
  paramsBottom: Float32Array,
  paramsTop: Float32Array,
  settings: Map<string, number>,
  samplesCount: number = 1,
): Promise<OptimizationResult> {
  const params = sampleParams(samplesCount, paramsTop, paramsBottom);

  const extremums: Extremum[] = [];
  const warnings: string[] = [];
  const failedInitPoint: Float32Array[] = [];
  let failsCount = 0;
  let failsDF: DG.DataFrame | null = null;

  let i: number;
  let percentage = 0;
  const pi = DG.TaskBarProgressIndicator.create(`Fitting... (${percentage}%)`);

  for (i = 0; i < samplesCount; ++i) {
    try {
      extremums.push(await optimizeNM(objectiveFunc, params[i], settings, paramsBottom, paramsTop));

      pi.update(100 * (i + 1) / samplesCount, `Fitting...`);

      percentage = Math.floor(100 * (i + 1) / samplesCount);
      pi.update(percentage, `Fitting... (${percentage}%)`);

      await sleep(MS_TO_SLEEP);

      if ((pi as any).canceled)
        break;
    } catch (e) {
      pi.close();

      if (e instanceof InconsistentTables)
        throw new Error(`Inconsistent dataframes: ${e.message}`);

      ++failsCount;
      warnings.push((e instanceof Error) ? e.message : 'Platform issue');
      failedInitPoint.push(params[i]);
    }
  }

  pi.close();

  if (failsCount > 0) {
    const dim = paramsTop.length;
    const raw = new Array<Float32Array>(dim);

    for (let i = 0; i < dim; ++i)
      raw[i] = new Float32Array(failsCount);

    failedInitPoint.forEach((point, idx) => point.forEach((val, jdx) => raw[jdx][idx] = val));

    failsDF = DG.DataFrame.fromColumns(raw.map((arr, idx) => DG.Column.fromFloat32Array(`arg${idx}`, arr)));
    failsDF.columns.add(DG.Column.fromStrings('Issue', warnings));
  }

  return {
    extremums: extremums,
    fails: failsDF,
  };
}

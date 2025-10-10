/** the Nelder-Mead optimizer */
import * as DG from 'datagrok-api/dg';

import {Extremum, OptimizationResult, InconsistentTables, ValueBoundsData, throttle} from './optimizer-misc';
import {optimizeNM} from './optimizer-nelder-mead';
import {sampleParamsWithFormulaBounds} from './optimizer-sampler';
import {TIMEOUT, ReproSettings, EarlyStoppingSettings} from './constants';
import {seededRandom} from './fitting-utils';

export async function performNelderMeadOptimization(
  {
    objectiveFunc,
    inputsBounds,
    samplesCount = 1,
    settings,
    reproSettings,
    earlyStoppingSettings,
  }: {
    objectiveFunc: (x: Float64Array) => Promise<number>;
    inputsBounds: Record<string, ValueBoundsData>;
    samplesCount: number,
    settings: Map<string, number>;
    reproSettings: ReproSettings;
    earlyStoppingSettings: EarlyStoppingSettings;
  },
): Promise<OptimizationResult> {
  const rand = reproSettings.reproducible ? seededRandom(reproSettings.seed) : Math.random;
  const [params, paramsBottom, paramsTop] = sampleParamsWithFormulaBounds(samplesCount, inputsBounds, rand);

  let extremums: Extremum[] = [];
  const warnings: string[] = [];
  const failedInitPoint: Float64Array[] = [];
  let failsCount = 0;
  let failsDF: DG.DataFrame | null = null;

  let i: number;
  let percentage = 0;
  const pi = DG.TaskBarProgressIndicator.create(`Fitting... (${percentage}%)`);

  const useEarlyStopping = earlyStoppingSettings.useEarlyStopping;
  const threshold = useEarlyStopping ? earlyStoppingSettings.costFuncThreshold : undefined;
  const maxValidPoints = earlyStoppingSettings.stopAfter;

  let validPointsCount = 0;
  let lastWorkStartTs: number | undefined;

  for (i = 0; i < samplesCount; ++i) {
    try {
      lastWorkStartTs = await throttle(TIMEOUT.MS_TO_SLEEP * 4, TIMEOUT.MS_TO_SLEEP, lastWorkStartTs);

      if ((pi as any).canceled)
        break;

      const extremum = await optimizeNM(objectiveFunc, params[i], settings, paramsBottom[i], paramsTop[i], threshold);

      if (useEarlyStopping) {
        if (extremum.cost <= threshold!) {
          extremums.push(extremum);
          ++validPointsCount;
        }

        if (validPointsCount >= maxValidPoints)
          break;
      } else
        extremums.push(extremum);

      percentage = Math.floor(100 * (i + 1) / samplesCount);
      pi.update(percentage, `Fitting... (${percentage}%)`);
    } catch (e) {
      if (e instanceof InconsistentTables) {
        pi.close();
        throw new Error(`Inconsistent dataframes: ${e.message}`);
      }

      ++failsCount;
      warnings.push((e instanceof Error) ? e.message : 'Platform issue');
      failedInitPoint.push(params[i]);
    }
  }

  pi.close();

  if (failsCount > 0) {
    const dim = paramsTop.length;
    const raw = new Array<Float64Array>(dim);

    for (let i = 0; i < dim; ++i)
      raw[i] = new Float64Array(failsCount);

    failedInitPoint.forEach((point, idx) => point.forEach((val, jdx) => raw[jdx][idx] = val));

    failsDF = DG.DataFrame.fromColumns(raw.map((arr, idx) => DG.Column.fromFloat64Array(`arg${idx}`, arr)));
    failsDF.columns.add(DG.Column.fromStrings('Issue', warnings));
  }

  return {
    extremums: extremums,
    fails: failsDF,
  };
}

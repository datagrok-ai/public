import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {LOSS, TIMEOUT} from './constants';
import {getInputsData, makeGetCalledFuncCall} from './fitting-utils';
import {OutputTargetItem, throttle, ValueBoundsData} from './optimizer-misc';
import {makeBoundsChecker} from './optimizer-sampler';
import {accumulateLoss, ColLike, DfLike, FrameTarget, ScalarTarget} from './worker/cost-math';


export function makeConstFunction(
  type: LOSS,
  func: DG.Func,
  bounds: Record<string, ValueBoundsData>,
  outputTargets: OutputTargetItem[],
) {
  if (type !== LOSS.MAD && type !== LOSS.RMSE)
    throw new Error(`Unknown type ${type}`);
  const {variedInputNames, fixedInputs} = getInputsData(bounds);
  return makeCostFunc(type === LOSS.RMSE, func, bounds, fixedInputs, variedInputNames, outputTargets);
}

function makeCostFunc(useRmse: boolean, func: DG.Func, bounds: Record<string, ValueBoundsData>,
  inputs: Record<string, any>, variedInputNames: string[], outputTargets: OutputTargetItem[]) {
  const getCalledFuncCall = makeGetCalledFuncCall(func, inputs, variedInputNames, false);
  const boundsChecker = makeBoundsChecker(bounds, variedInputNames);
  let lastWorkStartTs: number | undefined = undefined;

  return async (x: Float64Array): Promise<number | undefined> => {
    if (!boundsChecker(x))
      return;

    const calledFuncCall = await getCalledFuncCall(x);
    lastWorkStartTs = await throttle(TIMEOUT.MS_TO_SLEEP * 19, TIMEOUT.MS_TO_SLEEP, lastWorkStartTs);

    const scalars: ScalarTarget[] = [];
    const frames: FrameTarget[] = [];
    for (const output of outputTargets) {
      if (output.type !== DG.TYPE.DATA_FRAME) {
        scalars.push({
          target: output.target as number,
          sim: calledFuncCall.getParamValue(output.propName) as number,
        });
      } else {
        const df = output.target as DG.DataFrame;
        frames.push({
          argCol: df.col(output.argName) as ColLike | null,
          funcCols: output.cols as unknown as ColLike[],
          simDf: calledFuncCall.getParamValue(output.propName) as DfLike,
        });
      }
    }
    return accumulateLoss(useRmse, scalars, frames);
  };
}

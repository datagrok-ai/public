import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {LOSS} from './constants';
import {getErrors, getInputsData, makeGetCalledFuncCall} from './fitting-utils';
import {OutputTargetItem, ValueBoundsData} from './optimizer-misc';
import {makeBoundsChecker} from './optimizer-sampler';


export function makeConstFunction(
  type: LOSS,
  func: DG.Func,
  bounds: Record<string, ValueBoundsData>,
  outputTargets: OutputTargetItem[],
) {
  const {variedInputNames, fixedInputs} = getInputsData(bounds);
  if (type === LOSS.MAD)
    return makeMadCostFunc(func, bounds, fixedInputs, variedInputNames, outputTargets);
  if (type === LOSS.RMSE)
    return makeRmseCostFunc(func, bounds, fixedInputs, variedInputNames, outputTargets);
  throw new Error(`Unknown type ${type}`);
}

function makeMadCostFunc(func: DG.Func, bounds: Record<string, ValueBoundsData>, inputs: Record<string, any>, variedInputNames: string[], outputTargets: OutputTargetItem[]) {
  /** Maximum absolute deviation (MAD) cost function */
  const getCalledFuncCall = makeGetCalledFuncCall(func, inputs, variedInputNames);
  const boundsChecker = makeBoundsChecker(bounds, variedInputNames)
  const madCostFunc = async (x: Float64Array): Promise<number | undefined> => {
    if(!boundsChecker(x))
      return;

    const calledFuncCall = await getCalledFuncCall(x);

    let mad = 0;

    outputTargets.forEach((output) => {
      if (output.type !== DG.TYPE.DATA_FRAME)
        mad = Math.max(mad, Math.abs(output.target as number - calledFuncCall.getParamValue(output.propName)));
      else {
        const df = output.target as DG.DataFrame;
        getErrors(df.col(output.argName), output.cols, calledFuncCall.getParamValue(output.propName), false)
          .forEach((err) => mad = Math.max(mad, Math.abs(err)));
      }
    });

    return mad;
  };
  return madCostFunc;
}

function makeRmseCostFunc(func: DG.Func,  bounds: Record<string, ValueBoundsData>, inputs: Record<string, any>, variedInputNames: string[], outputTargets: OutputTargetItem[]) {
  /** Root mean sqaure error (RMSE) cost function */
  const getCalledFuncCall = makeGetCalledFuncCall(func, inputs, variedInputNames);
  const boundsChecker = makeBoundsChecker(bounds, variedInputNames)
  const rmseCostFunc = async (x: Float64Array): Promise<number| undefined> => {
    if(!boundsChecker(x))
      return;

    const calledFuncCall = await getCalledFuncCall(x);

    let sumOfSquaredErrors = 0;
    let outputsCount = 0;
    let cur = 0;

    outputTargets.forEach((output) => {
      if (output.type !== DG.TYPE.DATA_FRAME) {
        cur = output.target as number;
        sumOfSquaredErrors += ((cur - calledFuncCall.getParamValue(output.propName)) / (cur !== 0 ? cur : 1)) ** 2;
        ++outputsCount;
      } else {
        const df = output.target;
        getErrors(df.col(output.argName), output.cols, calledFuncCall.getParamValue(output.propName), true)
          .forEach((err) => {
            sumOfSquaredErrors += err ** 2;
            ++outputsCount;
          });
      }
    });

    return Math.sqrt(sumOfSquaredErrors / outputsCount);
  };
  return rmseCostFunc;
}

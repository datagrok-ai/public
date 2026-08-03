/* eslint-disable valid-jsdoc */
import * as DG from 'datagrok-api/dg';

import {Extremum, OptimizationResult, OptimizerInputsConfig, OptimizerOutputsConfig,
  TargetTableOutput} from './optimizer-misc';
import {getInputsData, makeGetCalledFuncCall} from './fitting-utils';
import {getNonSimilar} from './similarity-utils';

export type FinalizedFitting = {
  allExtremums: Extremum[];
  selectedExtremums: Extremum[];
  calls: DG.FuncCall[];
  fails: DG.DataFrame | null;
};

export type FinalizeContext = {
  func: DG.Func;
  inputBounds: OptimizerInputsConfig;
  outputTargets: OptimizerOutputsConfig;
  similarity: number;
};

export async function finalizeOptimizationResult(
  optResult: OptimizationResult,
  ctx: FinalizeContext,
): Promise<FinalizedFitting> {
  const {variedInputNames, fixedInputs} = getInputsData(ctx.inputBounds);

  const allExtremums = optResult.extremums.slice();
  allExtremums.sort((a: Extremum, b: Extremum) => a.cost - b.cost);

  const getCalledFuncCall = makeGetCalledFuncCall(ctx.func, fixedInputs, variedInputNames, true);
  const targetDfs: TargetTableOutput[] = ctx.outputTargets
    .filter((output) => output.type === DG.TYPE.DATA_FRAME)
    .map((output) => ({
      name: output.propName,
      target: output.target as DG.DataFrame,
      argColName: output.argName,
    }));

  const selectedExtremums = await getNonSimilar(allExtremums, ctx.similarity, getCalledFuncCall, targetDfs);

  const calls: DG.FuncCall[] = [];
  for (const extr of selectedExtremums)
    calls.push(await getCalledFuncCall(extr.point));

  return {allExtremums, selectedExtremums, calls, fails: optResult.fails};
}

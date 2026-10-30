/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {
  HeadlessResult, ItemMetadata, PipelineInstanceConfigInput, RunState, normalizePipelineInstanceConfig,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/index';

type ConstructorTypeOf<T> = new (...args:any[]) => T;

declare global {
  interface Window {
    compute: {
      fileInput: typeof UiUtils.fileInput,
      historyInput: typeof UiUtils.historyInput,
      historyInputJSON: typeof UiUtils.historyInputJSON,
      historyPanel: typeof UiUtils.historyPanel,

      CompView: ConstructorTypeOf<ComputationView>,
      RFV: ConstructorTypeOf<RichFunctionView>,
      Pipeline: ConstructorTypeOf<PipelineView>,
      CFV: ConstructorTypeOf<CustomFunctionView>,

      makeValidationResult: typeof makeValidationResult,
      makeAdvice: typeof makeAdvice,
      makeRevalidation: typeof makeRevalidation,
      mergeValidationResults: typeof mergeValidationResults,

      deepCopy: typeof deepCopy,
    },
  }
}

import * as UiUtils from './src/ui-utils';
export {UiUtils};

import {deepCopy} from '@datagrok-libraries/compute-utils/shared-utils/utils';
export {deepCopy};

import {
  makeValidationResult, makeAdvice, makeRevalidation, mergeValidationResults,
  ValidationInfo,
} from './src/validation';
export {
  makeValidationResult, makeAdvice, makeRevalidation, mergeValidationResults,
  ValidationInfo,
};

export type * from '@datagrok-libraries/compute-utils/reactive-tree-driver/index';
export {isFuncCallState, isStaticPipelineState, isDynamicPipelineState} from '@datagrok-libraries/compute-utils/reactive-tree-driver/index';

import {
  ComputationView, RichFunctionView, PipelineView, CustomFunctionView,
  createCompView, createRFV, createPipeline, createCFV
} from './src/views';
export {
  ComputationView, RichFunctionView, PipelineView, CustomFunctionView,
  createCompView, createRFV, createPipeline, createCFV
};

export async function initComputeApi() {
  const initFunc = DG.Func.find({package: 'Compute', name: 'init'})[0];
  if (initFunc)
    await initFunc.prepare().call();
}

export async function startWorkflow<T=any>(nqName: string, version: string, instanceConfig: PipelineInstanceConfigInput): Promise<T> {
  const startFunc = DG.Func.find({package: 'Compute2', name: 'StartWorkflow'})[0];
  const call = startFunc.prepare({nqName, version, instanceConfig: normalizePipelineInstanceConfig(instanceConfig)});
  await call.call();
  return call.getOutputParamValue();
}

/** Runs a workflow without UI: creates it from the initial configuration and, unless a step ready to run
 * has validation errors, runs every step that can run. */
export async function runWorkflow(
  nqName: string, version: string, instanceConfig: PipelineInstanceConfigInput = {},
): Promise<HeadlessResult> {
  const runFunc = DG.Func.find({package: 'Compute2', name: 'RunWorkflow'})[0];
  const call = runFunc.prepare({nqName, version, instanceConfig: normalizePipelineInstanceConfig(instanceConfig)});
  await call.call();
  return call.getOutputParamValue();
}

/** Saves the state of a {@link runWorkflow} result as a new workflow run in history. */
export async function saveWorkflowRun(
  nqName: string, version: string, state: RunState, metadata: ItemMetadata = {},
): Promise<DG.FuncCall> {
  const saveFunc = DG.Func.find({package: 'Compute2', name: 'SaveWorkflowRun'})[0];
  const call = saveFunc.prepare({nqName, version, state, metadata});
  await call.call();
  return call.getOutputParamValue();
}

import {OptimizerParams} from '@datagrok-libraries/compute-utils/function-views/src/fitting/optimizer-api';

export async function runOptimizer(params: OptimizerParams): Promise<DG.FuncCall[]> {
  const optimFunc = DG.Func.find({package: 'Compute2', name: 'RunOptimizer'})[0];
  const call = optimFunc.prepare({params});
  await call.call();
  return call.getOutputParamValue();
}

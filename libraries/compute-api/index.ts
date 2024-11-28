/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

type ConstructorTypeOf<T> = new (...args:any[]) => T;

declare global {
  interface Window {
    compute: {
      fileInput: typeof UiUtils.fileInput,
      historyInput: typeof UiUtils.historyInput,
      historyInputJSON: typeof UiUtils.historyInputJSON,
      historyPanel: typeof UiUtils.historyPanel,
      testPipeline: typeof testPipeline,

      CompView: ConstructorTypeOf<ComputationView>,
      RFV: ConstructorTypeOf<RichFunctionView>,
      Pipeline: ConstructorTypeOf<PipelineView>,
      CompositionPipeline: {
        new (...args:any[]): CompositionPipeline,
        compose: typeof composeCompositionPipeline
      }

      makeValidationResult: typeof makeValidationResult,
      makeAdvice: typeof makeAdvice,
      makeRevalidation: typeof makeRevalidation,
      mergeValidationResults: typeof mergeValidationResults,
    }
  }
}

import {testPipeline} from './src/utils';
export {testPipeline};

import * as UiUtils from './src/ui-utils';
export {UiUtils};

import {
  makeValidationResult, makeAdvice, makeRevalidation, mergeValidationResults,
  ValidationInfo,
} from './src/validation';
export {
  makeValidationResult, makeAdvice, makeRevalidation, mergeValidationResults,
  ValidationInfo,
};

import {
  ComputationView, RichFunctionView, PipelineView, CompositionPipeline,
  createCompView, createRFV, createPipeline, createCompositionPipeline, composeCompositionPipeline,
  PipelineCompositionConfiguration, PipelineConfiguration,
} from './src/views';
export {
  ComputationView, RichFunctionView, PipelineView, CompositionPipeline,
  createCompView, createRFV, createPipeline, createCompositionPipeline, composeCompositionPipeline,
  PipelineCompositionConfiguration, PipelineConfiguration,
};

export async function initComputeApi() {
  const initFunc = DG.Func.find({package: 'Compute', name: 'init'})[0];
  if (initFunc)
    await initFunc.prepare().call();
}

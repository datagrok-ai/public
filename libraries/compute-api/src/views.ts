/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import type {
  ComputationView as ComputationViewType,
  RichFunctionView as RichFunctionViewType,
  PipelineView as PipelineViewType,
} from '@datagrok-libraries/compute-utils';

import type {
  CustomFunctionView as CustomFunctionViewType,
} from '@datagrok-libraries/compute-utils/function-views/src/custom-function-view';
import type {
  HistoryInput as HistoryInputType, historyInput as historyInputType,
} from '@datagrok-libraries/compute-utils/function-views/src/history-input';
export type {HistoryInputOptions} from '@datagrok-libraries/compute-utils/function-views/src/history-input';

// a custom view subclasses window.compute2.CustomFunctionView, so every package runs Compute2's copy
export type CustomFunctionView = CustomFunctionViewType;
export type HistoryInput = HistoryInputType;

export function historyInput(...args: Parameters<typeof historyInputType>): ReturnType<typeof historyInputType> {
  return window.compute2.historyInput(...args);
}

export type ComputationView = ComputationViewType;
export type RichFunctionView = RichFunctionViewType;
export type PipelineView = PipelineViewType;

export function createCompView(
  ...args: ConstructorParameters<typeof ComputationViewType>
): ComputationViewType {
  return new window.compute.CompView(...args);
};

export function createRFV(
  ...args: ConstructorParameters<typeof RichFunctionViewType>
): RichFunctionViewType {
  return new window.compute.RFV(...args);
};

export function createPipeline(
  ...args: ConstructorParameters<typeof PipelineViewType>
): PipelineViewType {
  return new window.compute.Pipeline(...args);
};

/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import type {
  ComputationView as ComputationViewType,
  RichFunctionView as RichFunctionViewType,
  PipelineView as PipelineViewType,
} from '@datagrok-libraries/compute-utils';

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

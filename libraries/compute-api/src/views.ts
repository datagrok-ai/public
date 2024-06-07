/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import type {
  PipelineView as PipelineViewType,
  RichFunctionView as RichFunctionViewType,
} from '@datagrok-libraries/compute-utils';

export type PipelineView = PipelineViewType;
export type RichFunctionView = RichFunctionViewType;

//@ts-ignore
export const Pipeline = (window.compute.Pipeline) as
  (new (...args: ConstructorParameters<typeof PipelineViewType>) => PipelineViewType);

//@ts-ignore
export const RFV = (window.compute.RFV) as
  (new (...args: ConstructorParameters<typeof RichFunctionViewType>) => RichFunctionViewType);

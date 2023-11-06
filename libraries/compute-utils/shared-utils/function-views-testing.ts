import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {FunctionView, PipelineView, RichFunctionView} from '../function-views';
import {ExpectDeepEqualOptions, expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {delay, last, takeUntil, takeWhile} from 'rxjs/operators';
import {of} from 'rxjs';

export interface FunctionViewTestOptions extends ExpectDeepEqualOptions {
  validatorsWaitTime?: number;
}
export type PipelineTestOptions = Record<string, FunctionViewTestOptions>

export async function testPipeline(spec: any, view: PipelineView, options: PipelineTestOptions = {}) {
  for (const [name, step] of Object.entries(view.steps)) {
    const stepSpec = spec[name];
    if (!stepSpec)
      continue;
    const prefix = options.prefix ? `${options.prefix}: ${name}` : name;
    await testFunctionView(stepSpec, step.view, {validatorsWaitTime: 1000, ...options, prefix});
  }
}

export async function testFunctionView(
  spec: any, view: FunctionView | RichFunctionView, options: FunctionViewTestOptions) {
  for (const [name, data] of Object.entries(spec.inputs))
    view.funcCall.inputs[name] = data;
  if (view instanceof RichFunctionView) {
    const validatorsWaitTime = options.validatorsWaitTime ?? 1000;
    const state = await view.pendingValidations.pipe(
      takeUntil(of(null).pipe(delay(validatorsWaitTime))),
      takeWhile(() => !view.isValid()),
      last(),
    ).toPromise();
    if (state)
      throw new Error(`{options.prefix}: validators failed to accept input in ${validatorsWaitTime}ms, state ${state}`);
  }
  await view.run();
  expectDeepEqual(view.funcCall.outputs, spec, options);
}

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {FunctionView, PipelineView, RichFunctionView} from '../function-views';
import {ExpectDeepEqualOptions, expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {delay, last, takeUntil, takeWhile} from 'rxjs/operators';
import {of} from 'rxjs';
import {fcInputFromSerializable} from './utils';

export interface FunctionViewTestOptions extends ExpectDeepEqualOptions {
  initWaitTimeout?: number;
  validatorsWaitTimeout?: number;
  nextStepTimeout?: number;
  updateMode?: boolean;
}
export interface PipelineTestOptions {
  initWaitTimeout?: number;
  updateMode?: boolean;
  stepOptions?: Record<string, FunctionViewTestOptions>
}

export async function testPipeline(
  spec: any, view: PipelineView, options: PipelineTestOptions = {}) {
  await waitForReady(view, options.initWaitTimeout ?? 5000);
  for (const [name, step] of Object.entries(view.steps)) {
    const stepSpec = spec[name];
    if (!stepSpec)
      continue;
    const stepOptions = options.stepOptions?.[name] ?? {};
    const prefix = stepOptions.prefix ? `${stepOptions.prefix}: ${name}` : name;
    const {updateMode} = options;
    const defaultSettings = {validatorsWaitTimeout: 1000, nextStepTimeout: 100, initWaitTimeout: 100, updateMode};
    await testFunctionView(stepSpec, step.view, {...defaultSettings, ...stepOptions, prefix});
  }
}

export async function testFunctionView(
  spec: any, view: FunctionView | RichFunctionView, options: FunctionViewTestOptions = {}) {
  await waitForReady(view, options.initWaitTimeout ?? 1000);
  for (const [name, data] of Object.entries(spec.inputs)) {
    const propertyType = view.funcCall.inputParams[name].property.propertyType;
    if (!options.updateMode)
      view.funcCall.inputs[name] = await fcInputFromSerializable(propertyType, data);
    else {
      const currentValue = view.funcCall.inputs[name];
      if (currentValue == null)
        view.funcCall.inputs[name] = await fcInputFromSerializable(propertyType, data);
    }
  }
  if (view instanceof RichFunctionView) {
    const validatorsWaitTimeout = options.validatorsWaitTimeout ?? 1000;
    const state = await view.pendingValidations.pipe(
      takeUntil(of(null).pipe(delay(validatorsWaitTimeout))),
      takeWhile(() => !view.isValid(), true),
      last(),
    ).toPromise();
    if (state) {
      const msg = `{options.prefix}: validators failed to accept input in ${validatorsWaitTimeout}ms, state ${state}`;
      throw new Error(msg);
    }
  }
  console.log(`running ${view.func.nqName}`);
  await view.run();
  if (!options.updateMode) {
    console.log(`checking ${view.func.nqName} results`);
    expectDeepEqual(view.funcCall.outputs, spec.outputs, options);
    console.log(`${view.func.nqName} ok`);
  }
}

async function waitForReady(view: FunctionView, timeout: number) {
  const isReady = await view.isReady.pipe(
    takeUntil(of(null).pipe(delay(timeout))),
    takeWhile((ready) => !ready, true),
    last(),
  ).toPromise();
  if (!isReady)
    throw new Error(`Waiting for view ${view.name} timeout`);
}

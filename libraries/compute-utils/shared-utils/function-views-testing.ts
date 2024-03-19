import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {FunctionView, PipelineView, RichFunctionView} from '../function-views';
import {ExpectDeepEqualOptions, expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {delay, last, takeUntil, takeWhile, map, debounceTime} from 'rxjs/operators';
import {of} from 'rxjs';
import {fcInputFromSerializable} from './utils';

export interface FunctionViewTestOptions extends ExpectDeepEqualOptions {
  initWaitTimeout?: number;
  validatorsWaitTimeout?: number;
  updateMode?: boolean;
  parent?: PipelineView;
  parentWaitTimeout?: number;
}

export interface PipelineTestOptions {
  initWaitTimeout?: number;
  updateMode?: boolean;
  stepOptions?: Record<string, FunctionViewTestOptions>
}

const defaultInitTimeout = 20000;
const defaultValidatorsTimeout = 5000;
const defaultParentTimeout = 5000;

export async function testPipeline(
  spec: any, view: PipelineView, options: PipelineTestOptions = {}) {
  await waitForViewReady(view, options.initWaitTimeout ?? defaultInitTimeout);
  for (const [name, step] of Object.entries(view.steps)) {
    const stepSpec = spec[name];
    if (!stepSpec)
      continue;
    const stepOptions = options.stepOptions?.[name] ?? {};
    const prefix = stepOptions.prefix ? `${stepOptions.prefix}: ${name}` : name;
    const {updateMode} = options;
    const defaultSettings = {
      validatorsWaitTimeout: defaultValidatorsTimeout,
      initWaitTimeout: defaultInitTimeout,
      parent: view,
      updateMode
    };
    await testFunctionView(stepSpec, step.view, {...defaultSettings, ...stepOptions, prefix});
  }
}

export async function testFunctionView(
  spec: any, view: FunctionView | RichFunctionView, options: FunctionViewTestOptions = {}
) {
  await waitForViewReady(view, options.initWaitTimeout ?? defaultInitTimeout);
  if (options.parent) {
    const parentWaitTimeout = options.parentWaitTimeout ?? defaultParentTimeout;
    const runningUpdates = await options.parent.isUpdating.pipe(
      map(() => options.parent!.getRunningUpdates()),
      debounceTime(100),
      takeUntil(of(null).pipe(delay(parentWaitTimeout))),
      takeWhile((updates) => updates.length !== 0, true),
      last(),
    ).toPromise();
    if (runningUpdates?.length) {
      const msg = `{options.prefix}: parent updates failed to complete in ${parentWaitTimeout}ms, updates: ${runningUpdates.join(',')}`;
      throw new Error(msg);
    }
  }
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
    const validatorsWaitTimeout = options.validatorsWaitTimeout ?? defaultValidatorsTimeout;
    const pendingValidators = await view.pendingValidations.pipe(
      map(vals => Object.keys(vals)),
      debounceTime(100),
      takeUntil(of(null).pipe(delay(validatorsWaitTimeout))),
      takeWhile(() => !view.isValid(), true),
      last(),
    ).toPromise();
    if (pendingValidators?.length) {
      const msg = `{options.prefix}: validators failed to accept input in ${validatorsWaitTimeout}ms, inputs: ${pendingValidators.join(',')}`;
      throw new Error(msg);
    }
  }
  console.log(`running ${view.func.nqName}`);
  await view.run();
  if (!options.updateMode) {
    console.log(`checking ${view.func.nqName} results`);
    try {
      expectDeepEqual(view.funcCall.outputs, spec.outputs, options);
    } catch (e) {
      grok.shell.error(String(e));
      throw e;
    }
    console.log(`${view.func.nqName} ok`);
  }
}

async function waitForViewReady(view: FunctionView, timeout: number) {
  const isReady = await view.isReady.pipe(
    takeUntil(of(null).pipe(delay(timeout))),
    takeWhile((ready) => !ready, true),
    last(),
  ).toPromise();
  if (!isReady)
    throw new Error(`Waiting for view ${view.name} timeout`);
}

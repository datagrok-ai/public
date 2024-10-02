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
  interactive?: boolean;
  parent?: PipelineView;
  parentWaitTimeout?: number;
}

export interface PipelineTestOptions {
  initWaitTimeout?: number;
  updateMode?: boolean;
  interactive?: boolean;
  stepOptions?: Record<string, FunctionViewTestOptions>
}

const defaultInitTimeout = 20000;
const defaultValidatorsTimeout = 5000;
const defaultParentTimeout = 5000;

export async function testPipeline(
  spec: any, view: PipelineView, options: PipelineTestOptions = {},
) {
  await waitForViewReady(view, options.initWaitTimeout ?? defaultInitTimeout);
  view.disableSetters(true);
  try {
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
        interactive: options.interactive,
        updateMode,
      };
      await testFunctionView(stepSpec, step.view, {...defaultSettings, ...stepOptions, prefix});
    }
  } finally {
    view.disableSetters(false);
  }
  if (options.interactive)
    grok.shell.info(`${view.name} test passed`);
}

export async function testFunctionView(
  spec: any, view: FunctionView | RichFunctionView, options: FunctionViewTestOptions = {},
) {
  // For RFV step inside of a pipeline:
  // 0. waiting for step view ready
  // 1. waiting for previous step releted parent logic to complete
  // 2. setting current step inputs
  // 3. waiting for current step inputs related parent logic to complete
  // 4. waiting for RFV validators to complete
  // 5. run the step

  await waitForViewReady(view, options.initWaitTimeout ?? defaultInitTimeout);
  const waitForParentReady = async () => {
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
        const err = new Error(msg);
        if (options.interactive)
          grok.shell.error(String(err));
        throw err;
      }
    }
  };
  await waitForParentReady();
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
  await waitForParentReady();
  const waitForValidators = async () => {
    if (view instanceof RichFunctionView) {
      const validatorsWaitTimeout = options.validatorsWaitTimeout ?? defaultValidatorsTimeout;
      const pendingValidators = await view.pendingInputValidations.pipe(
        map((vals) => Object.keys(vals)),
        debounceTime(100),
        takeUntil(of(null).pipe(delay(validatorsWaitTimeout))),
        takeWhile(() => !view.isValid(), true),
        last(),
      ).toPromise();
      if (pendingValidators?.length) {
        const msg = `{options.prefix}: validators failed to accept input in ${validatorsWaitTimeout}ms, inputs: ${pendingValidators.join(',')}`;
        const err = new Error(msg);
        if (options.interactive)
          grok.shell.error(String(err));
        throw err;
      }
    }
  };
  await waitForValidators();
  console.log(`running ${view.func.nqName}`);
  await view.run();
  if (!options.updateMode) {
    console.log(`checking ${view.func.nqName} results`);
    try {
      expectDeepEqual(view.funcCall.outputs, spec.outputs, options);
      if (!options.parent && options.interactive)
        grok.shell.info(`${view.name} test passed`);
    } catch (e) {
      if (options.interactive)
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

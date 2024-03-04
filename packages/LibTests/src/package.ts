/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
// eslint-disable-next-line
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {FuncCallInput} from '@datagrok-libraries/compute-utils/shared-utils/input-wrappers';
import {BehaviorSubject} from 'rxjs';
import {distinctUntilChanged} from 'rxjs/operators';
import equal from 'deep-equal';
import {ValidationInfo, makeAdvice, makeRevalidation, makeValidationResult}
  from '@datagrok-libraries/compute-utils/shared-utils/validation';
import { CompositionPipeline } from '@datagrok-libraries/compute-utils/composition-pipeline/src/composition-pipeline';
import { delay } from '@datagrok-libraries/utils/src/test';

export const _package = new DG.Package();

class InputMock implements FuncCallInput {
  _value = new BehaviorSubject<any>(null);
  notify = true;
  enabled = true;
  root = ui.div('', {style: {width: '100%'}});
  // fake input
  input = ui.switchInput('test', true);

  constructor() {
    this.root.append(this.input.root);
  }

  set value(v: any) {
    const nv = v ? JSON.parse(v) : '';
    this._value.next(nv);
  }

  get value() {
    if (!this._value.value)
      return '';

    return JSON.stringify(this._value.value);
  }

  onInput(fn: Function) {
    return this._value.pipe(distinctUntilChanged(equal)).subscribe(() => {
      if (this.notify)
        fn(this.value);
    });
  }
}

//name: CustomInputMock
//input: object params
//output: object input
export async function CustomInputMock(params: any): Promise<FuncCallInput> {
  console.log(params);
  return new InputMock();
}

//name: RangeValidatorFactory
//input: object params
//output: object validator
export function RangeValidatorFactory(params: any) {
  const {min, max} = params;
  return (val: number) => {
    if (val < min || val > max)
      return makeValidationResult({errors: [`Out of range [${min}, ${max}] value: ${val}`]});

  };
}

//name: AsyncValidatorDemoFactory
//input: object params
//output: object validator
export function AsyncValidatorDemoFactory(params: any) {
  return async (val: number) => {
    await new Promise((resolve) => setTimeout(resolve, 100));
    if (val === 0)
      return makeValidationResult({warnings: [`Try non-null value`]});

  };
}

//name: GlobalValidatorDemoFactory
//input: object params
//output: object validator
export function GlobalValidatorDemoFactory(params: any) {
  const {max} = params;
  return async (_val: number, info: ValidationInfo) => {
    if (info.isRevalidation) {
      if (info.context?.isOk)
        return makeValidationResult();
      return makeValidationResult({warnings: [`Try lowering a value as well`]});
    }
    await new Promise((resolve) => setTimeout(resolve, 100));
    const {a, b, c} = info.funcCall.inputs;
    const s = a + b + c;
    const isOk = s <= max;
    const valRes = isOk ? makeValidationResult() : makeValidationResult({warnings: [`Try lowering a value`]});
    const fields = ['a', 'b', 'c'].filter((p) => p !== info.param);
    return makeRevalidation(fields, {isOk}, valRes);
  };
}

//name: ValidatorActionsDemoFactory
//input: object params
//output: object validator
export function ValidatorActionsDemoFactory(params: any) {
  return async (_val: number, info: ValidationInfo) => {
    const notifications = [makeAdvice('Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.', [
      {actionName: 'First action', action: () => grok.shell.info('First action')},
      {actionName: 'Another action', action: () => {
        ui.dialog({title: 'Another action'}).show({center: true, fullScreen: true});
      }},
    ])];
    if (info.lastCall && info.funcCall.inputs.x !== info.lastCall.inputs.x) {
      const delta = info.funcCall.inputs.x - info.lastCall.inputs.x;
      const warnings = [makeAdvice(`Param delta change ${delta}`)];
      return makeValidationResult({warnings, notifications});
    }
    return makeValidationResult({notifications});
  };
}


//name: MockHook1
//input: object params
export function MockHook1(params: any) {
  return;
}

//name: MockHook2
//input: object params
export function MockHook2(params: any) {
  return;
}

//name: MockHandler1
//input: object params
export function MockHandler1(params: any) {
  return;
}

//name: MockHandler2
//input: object params
export function MockHandler2(params: any) {
  return;
}


//name: MockPopup1
//input: double a
//input: double b
//input: double c
//output: double res
export function MockPopup1(a: number, b: number, c: number) {
  return a * b * c;
}

//name: MockPopup2
//input: double a
//input: double b
//input: double c
//output: double res
export function MockPopup2(a: number, b: number, c: number) {
  return a + b + c;
}


//name: TestFn1
//input: double a
//input: double b
//input: double c
//output: double res
export function TestFn1(a: number, b: number, c: number) {
  return a + b + c;
}

//name: TestFn2
//input: double a
//input: double b
//input: double c
//output: double res
export function TestFn2(a: number, b: number, c: number) {
  return a * b * c;
}

//name: TestFn3
//input: double a
//input: double b
//input: double c
//input: double d
//output: double res
export function TestFn3(a: number, b: number, c: number, d: number) {
  return a + b + c + d;
}

//name: TestFn4
//input: double a
//input: double b
//input: double c
//input: double d
//output: double res
export function TestFn4(a: number, b: number, c: number, d: number) {
  return a * b * c * d;
}

//name: TestFn5
//input: double a
//input: double b
//input: double c
//input: double d
//input: double e
//output: double res
export function TestFn5(a: number, b: number, c: number, d: number, e: number) {
  return a + b + c + d + e;
}

//name: TestFn6
//input: double a
//input: double b
//input: double c
//input: double d
//input: double e
//output: double res
export function TestFn6(a: number, b: number, c: number, d: number, e: number) {
  return a * b * c * d * e;
}


//name: TestCompositionPipeline1
export async function TestCompositionPipeline1() {
  const pipeline = new CompositionPipeline({
    id: 'testPipeline',
    nqName: 'LibTests:MockWrapper',
    steps: [
      {
        id: 'step1',
        nqName: 'LibTests:AddMock',
      },
      {
        id: 'step2',
        nqName: 'LibTests:MulMock',
      },
    ],
    links: [{
      id: 'link1',
      from: ['step1', 'res'],
      to: ['step2', 'a'],
    }]
  });
  grok.shell.addView(pipeline.makePipelineView());
  await pipeline.init();
}

//name: TestCompositionPipeline2
export async function TestCompositionPipeline2() {
  const pipeline = new CompositionPipeline({
    id: 'testPipeline',
    nqName: 'LibTests:MockWrapper',
    steps: [
      {
        id: 'step1',
        nqName: 'LibTests:AddMock',
      },
      {
        id: 'step2',
        nqName: 'LibTests:MulMock',
      },
    ],
    links: [{
      id: 'link1',
      from: ['step1', 'a'],
      to: ['step2', 'a'],
    }]
  });
  grok.shell.addView(pipeline.makePipelineView());
  await pipeline.init();
}

//name: TestCompositionPipeline3
export async function TestCompositionPipeline3() {
  const pipeline = new CompositionPipeline({
    id: 'testPipeline',
    nqName: 'LibTests:MockWrapper',
    steps: [
      {
        id: 'step1',
        nqName: 'LibTests:AddMock',
      },
      {
        id: 'step2',
        nqName: 'LibTests:MulMock',
      },
    ],
    links: [{
      id: 'link1',
      from: ['step1', 'a'],
      to: ['step2', 'a'],
      handler: async ({controller}) => {
        const val = controller.getState(['step1', 'a']);
        controller.updateState(['step2', 'a'], val * 2);
      }
    }]
  });
  grok.shell.addView(pipeline.makePipelineView());
  await pipeline.init();
}

//name: TestCompositionPipeline4
export async function TestCompositionPipeline4() {
  const pipeline = new CompositionPipeline({
    id: 'testPipeline',
    nqName: 'LibTests:MockWrapper',
    steps: [
      {
        id: 'step1',
        nqName: 'LibTests:AddMock',
      },
      {
        id: 'step2',
        nqName: 'LibTests:MulMock',
      },
    ],
    links: [{
      id: 'link1',
      from: ['step1', 'a'],
      to: ['step2', 'a'],
      handler: async ({controller}) => {
        const val = controller.getState(['step1', 'a']);
        await delay(2000);
        controller.updateState(['step2', 'a'], val * 2);
      }
    }]
  });
  grok.shell.addView(pipeline.makePipelineView());
  await pipeline.init();
}

/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';

import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import type {ViewerT, InputFormT} from '@datagrok-libraries/webcomponents';

export const _package = new DG.Package();

//tags: test
export async function TestViewerComponent() {
  await DG.Func.byName('WebComponents:init').prepare().call();
  const view = new DG.ViewBase();
  const viewerComponent = document.createElement('dg-viewer') as ViewerT;

  const setSrcBtn1 = ui.button('Set source demog', () => {
    viewerComponent.dataFrame = grok.data.demo.demog();
  });
  const setSrcBtn2 = ui.button('Set source doseResponse', () => {
    viewerComponent.dataFrame = grok.data.demo.doseResponse();
  });

  const remSrcBtn = ui.button('Remove source', () => {
    viewerComponent.dataFrame = undefined;
  });

  const setViewerTypeBtn1 = ui.button('Line chart', () => {
    viewerComponent.type = 'Line chart';
  });
  const setViewerTypeBtn2 = ui.button('Grid', () => {
    viewerComponent.type = 'Grid';
  });

  const setViewerTypeBtn3 = ui.button('Remove type', () => {
    viewerComponent.type = undefined;
  });
  viewerComponent.style.height = '100%';
  view.root.insertAdjacentElement('beforeend', setSrcBtn1);
  view.root.insertAdjacentElement('beforeend', setSrcBtn2);
  view.root.insertAdjacentElement('beforeend', remSrcBtn);
  view.root.insertAdjacentElement('beforeend', setViewerTypeBtn1);
  view.root.insertAdjacentElement('beforeend', setViewerTypeBtn2);
  view.root.insertAdjacentElement('beforeend', setViewerTypeBtn3);
  view.root.insertAdjacentElement('beforeend', viewerComponent);

  grok.shell.addView(view);
}

//tags: test
export async function TestFromComponent() {
  await DG.Func.byName('WebComponents:init').prepare().call();
  const func: DG.Func = await grok.functions.eval('LibTests:simpleInputs');
  const fc1 = func.prepare({
    a: 1,
    b: 2,
    c: 3,
  });
  const formComponent = document.createElement('dg-input-form') as InputFormT;
  formComponent.funcCall = fc1;

  const view = new DG.ViewBase();

  const replaceFnBtn = ui.button('Replace funcall', () => {
    const fc2 = func.prepare({
      a: 1,
      b: 2,
      c: 3,
    });
    formComponent.funcCall = fc2;
  });

  const showFormFcInputsBtn = ui.button('Log funcall inputs', () => {
    console.log(Object.entries(formComponent.funcCall!.inputs));
  });

  view.root.insertAdjacentElement('beforeend', showFormFcInputsBtn);
  view.root.insertAdjacentElement('beforeend', replaceFnBtn);
  view.root.insertAdjacentElement('beforeend', formComponent);
  grok.shell.addView(view);
}

//tags: test
export async function TestElements() {
  await DG.Func.byName('WebComponents:init').prepare().call();
  const bnt = document.createElement('button', {is: 'dg-button'});
  bnt.textContent = 'Click me';
  const bigBtn = document.createElement('button', {is: 'dg-big-button'});
  bigBtn.textContent = 'Click me';
  const view = new DG.ViewBase();
  view.root.insertAdjacentElement('beforeend', bnt);
  view.root.insertAdjacentElement('beforeend', bigBtn);
  grok.shell.addView(view);
}

// pipeline driver testing

//input: double a
//input: double b
//output: double res
export async function TestAdd2(a: number, b: number) {
  return a + b;
}

//input: double a
//input: double b
//output: double res
export async function TestSub2(a: number, b: number) {
  return a - b;
}

//input: double a
//input: double b
//output: double res
export async function TestMul2(a: number, b: number) {
  return a * b;
}

//input: double a
//input: double b
//output: double res
export async function TestDiv2(a: number, b: number) {
  return a / b;
}

//input: dataframe df
//output: dataframe res
export async function TestDF1(df: DG.DataFrame) {
  return df;
}

//output: dataframe res
export async function TestPresets() {
  return DG.DataFrame.fromColumns([
    DG.Column.fromList('string', 'preset', ['fast', 'exact']),
    DG.Column.fromList('double', 'a', [1, 10]),
  ]);
}

//input: file inputFile
//output: string result
export async function TestFileInput(inputFile: DG.FileInfo) {
  const bytes = await inputFile.readAsBytes();
  return `${inputFile.fileName ?? inputFile.name}: ${bytes.length} bytes`;
}

//input: double a
//input: double b
//output: double res
export async function TestAdd2Error(a: number, b: number) {
  if (a < 0 || b < 0)
    throw new Error('Test error');
  return a + b;
}

//input: double a
//input: double b
//input: double c
//input: double d
//input: double e
//output: double res
export async function TestMultiarg5(a: number, b: number, c: number, d: number, e: number) {
  return a + b + c + d + e;
}

/* Test fixture for template name-pairing: TestIONamesA's outputs and
   TestIONamesB's inputs share names, so deferred outputs()/inputs() links
   can be exercised end-to-end. */
//input: double seed
//output: double x
//output: double y
export async function TestIONamesA(seed: number) {
  return {x: seed, y: seed + 1};
}

//input: double x
//input: double y
//output: double res
export async function TestIONamesB(x: number, y: number) {
  return x + y;
}

//input: double y
//input: double x
//output: double res
export async function TestIONamesBReversed(y: number, x: number) {
  return x + y;
}

//input: double seed
//output: double x
//output: double y
//output: double z
export async function TestIONamesAExtra(seed: number) {
  return {x: seed, y: seed + 1, z: seed + 2};
}

//input: object params
//output: object result
export async function MockWrapper1(params: any) {
  const c: PipelineConfiguration = {
    id: 'pipeline1',
    nqName: 'LibTests:MockWrapper1',
    version: '1.0',
    type: 'static',
    steps: [
      {
        id: 'step1',
        nqName: 'LibTests:TestAdd2',
      },
      {
        id: 'step2',
        nqName: 'LibTests:TestMul2',
      },
    ],
    links: [{
      id: 'link1',
      from: 'in1:step1/res',
      to: 'out1:step2/a',
    }],
  };
  return c;
}

//input: object params
//output: object result
export async function MockWrapper2(params: any) {
  const c: PipelineConfiguration = {
    id: 'pipelinePar',
    nqName: 'LibTests:MockWrapper2',
    version: '1.0',
    type: 'parallel',
    stepTypes: [{
      id: 'stepAdd',
      nqName: 'LibTests:TestAdd2',
      friendlyName: 'add',
    }, {
      id: 'stepMul',
      nqName: 'LibTests:TestMul2',
      friendlyName: 'mul',
    }, {
      type: 'ref',
      provider: 'LibTests:MockWrapper1',
      version: '1.0',
    }],
    initialSteps: [
      {
        id: 'stepAdd',
      }, {
        id: 'pipeline1',
      },
    ],
  };
  return c;
}

//input: object params
//output: object result
export async function MockWrapper3(params: any) {
  const c: PipelineConfiguration = {
    id: 'pipelinePar',
    nqName: 'LibTests:MockWrapper3',
    version: '1.0',
    type: 'parallel',
    stepTypes: [{
      type: 'ref',
      provider: 'LibTests:MockWrapper2',
      version: '1.0',
    }],
    initialSteps: [
      {
        id: 'pipelinePar',
      },
    ],
  };
  return c;
}

//input: object params
//output: object result
export async function MockWrapper4(params: any) {
  const config2: PipelineConfiguration = {
    id: 'pipeline1',
    type: 'static',
    nqName: 'LibTests:MockWrapper4',
    version: '1.0',
    steps: [
      {
        id: 'step1',
        nqName: 'LibTests:TestAdd2Error',
      },
      {
        id: 'step2',
        nqName: 'LibTests:TestMul2',
      },
    ],
    links: [{
      id: 'link1',
      from: 'in1:step1/a',
      to: 'out1:step2/a',
      handler({controller}) {
        controller.setAll('out1', 2, 'restricted');
        return;
      },
    }],
  };
  return config2;
}


//input: object params
//output: object result
export async function MockWrapper5(params: any) {
  const config2: PipelineConfiguration = {
    id: 'pipeline1',
    type: 'sequential',
    nqName: 'LibTests:MockWrapper5',
    version: '1.0',
    approversGroup: 'MockGroup',
    stepTypes: [
      {
        id: 'step1',
        nqName: 'LibTests:TestAdd2Error',
        disableUIAdding: true,
        disableUIDragging: true,
        disableUIRemoving: true,
      },
      {
        id: 'step2',
        nqName: 'LibTests:TestMul2',
      },
      {
        id: 'pipeline2',
        type: 'static',
        disableUIDragging: true,
        disableUIRemoving: true,
        disableUIAdding: true,
        steps: [
          {
            id: 'step3',
            nqName: 'LibTests:TestSub2',

          },
          {
            id: 'step4',
            nqName: 'LibTests:TestDiv2',
          },
        ],
      },
    ],
    initialSteps: [{
      id: 'step1',
    }, {
      id: 'step2',
    }, {
      id: 'pipeline2',
    }],
    links: [{
      id: 'link1',
      from: 'in1:step1/a',
      to: 'out1:step2/a',
      handler({controller}) {
        controller.setAll('out1', 2, 'restricted');
        return;
      },
    }],
  };
  return config2;
}

//input: object params
//output: object result
export async function MockWrapperAction(params: any) {
  const c: PipelineConfiguration = {
    id: 'pipelineAct',
    nqName: 'LibTests:MockWrapperAction',
    version: '1.0',
    type: 'static',
    steps: [
      {
        id: 'step1',
        nqName: 'LibTests:TestAdd2',
      },
      {
        id: 'act1',
        type: 'action',
        friendlyName: 'My action',
      },
      {
        id: 'step2',
        nqName: 'LibTests:TestMul2',
      },
    ],
  };
  return c;
}

//input: object params
//output: object result
export async function MockWrapperDF(params: any) {
  const c: PipelineConfiguration = {
    id: 'pipelineDF',
    nqName: 'LibTests:MockWrapperDF',
    version: '1.0',
    type: 'static',
    steps: [
      {
        id: 'step1',
        nqName: 'LibTests:TestDF1',
      },
      {
        id: 'step2',
        nqName: 'LibTests:TestDF1',
      },
    ],
  };
  return c;
}

// annotation checks

//input: double a {nullable: true}
//input: double b {optional: true}
//input: double c
//input: int v = 5 {min: 0; max: 10}
//input: string code = "1234" {validator: /^[0-9]{4}$/i}
//input: string mode = "fast" {choices: ["fast", "exact"]}
//input: dataframe df
//input: column col {type: numerical; allowNulls: false}
//input: column mol {semType: Molecule; table: df}
//output: double res
export async function TestAnnotatedInputs(a: number, b: number, c: number, v: number, code: string, mode: string,
  df: DG.DataFrame, col: DG.Column, mol: DG.Column) {
  return 1;
}

//input: int x
//output: string res
export function MockValidator(x: number): string | null {
  return x < 10 ? null : 'too big';
}

//input: int x = 1 {validators: ["LibTests:MockValidator"]}
//input: int y = 1
//output: int res
export function TestNamedValidators(x: number, y: number): number {
  return x + y;
}

// annotation values

//input: string region
//output: list<string> res
export function MockCities(region: string): string[] {
  return region ? [`${region}-1`, `${region}-2`] : [];
}

//output: dataframe res
export function MockCars(): DG.DataFrame {
  return DG.DataFrame.fromColumns([
    DG.Column.fromList('string', 'model', ['Mazda', 'Volvo']),
    DG.Column.fromList('int', 'mpg', [21, 30]),
    DG.Column.fromList('int', 'CYL', [6, 4]),
  ]);
}

//input: int calc = 2 + 2
//input: string bare = high
//input: string metric = minkowski {choices: ["euclidean", "minkowski"]}
//input: string speed {choices: ["slow", "fast"]}
//input: string region = "EU"
//input: string city {choices: LibTests:MockCities}
//input: string model {choices: LibTests:MockCars(); propagateChoice: all}
//input: int mpg
//input: int cyl
//output: string res
export function TestValueAnnotations(calc: number, bare: string, metric: string, speed: string, region: string,
  city: string, model: string, mpg: number, cyl: number): string {
  return `${model} ${mpg} ${cyl}`;
}

//output: dataframe res
export function MockCarsTyped(): DG.DataFrame {
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromList('string', 'model', ['Mazda', 'Volvo']),
    DG.Column.fromList('string', 'cyl', ['4', 'abc']),
    DG.Column.fromList('double', 'mpg', [21.5, 30]),
    DG.Column.fromList('int', 'name', [1, 2]),
    DG.Column.fromList('string', 'flag', ['true', 'no']),
    DG.Column.fromList('string', 'engine', ['E1', 'E2']),
    DG.Column.fromList('string', 'made', ['2020-05-01T00:00:00Z', 'not a date']),
  ]);
  const dates = ['2021-03-04T05:06:07Z', '2022-01-02T00:00:00Z'];
  df.columns.addNewDateTime('when').init((i: number) => dayjs(dates[i]));
  return df;
}

//output: dataframe res
export function MockEngines(): DG.DataFrame {
  return DG.DataFrame.fromColumns([
    DG.Column.fromList('string', 'engine', ['E1', 'E2']),
    DG.Column.fromList('int', 'cyl', [8, 10]),
  ]);
}

//input: string model {choices: LibTests:MockCarsTyped(); propagateChoice: all}
//input: string engine {choices: LibTests:MockEngines(); propagateChoice: all}
//input: int cyl
//input: int mpg
//input: string name
//input: bool flag
//input: datetime when
//input: datetime made
//output: string res
export function TestLookupAnnotations(model: string, engine: string, cyl: number, mpg: number, name: string,
  flag: boolean, when: dayjs.Dayjs, made: dayjs.Dayjs): string {
  return `${model} ${engine}`;
}

//input: int x
//output: bool res
export function MockValidatorBool(x: number): boolean {
  return x < 10;
}

//input: int x
//output: string res
export function MockValidatorThrow(x: number): string {
  throw new Error('boom');
}

//input: int k = 2
//input: int hv = 1 {visible: k > 1}
//input: int foo = 5 {validator: bar > 3}
//input: double bar = 2
//input: string code = "12ab" {validator: startsWith(value, "12")}
//output: int res
export function TestExpressionInputs(k: number, hv: number, foo: number, bar: number, code: string): number {
  return k;
}

/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {getEachOutputItemInfo} from './outputTools';
import {VariedNumericalInputInfo, FixedInputItem} from './inputTools';

import { VarianceBasedSenstivityAnalysis } from './sensitivityAnalysis';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//name: Ishigami func
//input: double x1 = 1 {caption: x1; category: Varibles}
//input: double x2 = 2 {caption: x2; category: Varibles}
//input: double x3 = 3 {caption: x3; category: Varibles}
//input: double a = 7.0 {caption: a; category: Parameters}
//input: double b = 0.1 {caption: b; category: Parameters}
//output: double result
//editor: Compute:RichFunctionViewEditor
export function ishigamiFunc(x1: number, x2: number, x3: number, a : number, b : number): number {  
  const sin = Math.sin;
  const pow = Math.pow;
  /*let res = sin(x1) + a * pow(sin(x2), 2) + b * pow(x3, 4) * sin(x1);
  console.log(res);
  return res;*/
  return sin(x1) + a * pow(sin(x2), 2) + b * pow(x3, 4) * sin(x1);
}

//name: Parabola
//input: double x = 1 {caption: x; category: Varibles}
//input: double a = 7.0 {caption: a; category: Parameters}
//input: double b = 0.1 {caption: b; category: Parameters}
//output: double result
//editor: Compute:RichFunctionViewEditor
export function parabola(x: number, a : number, b : number): number {      
  return a * x * x + b;
}

//name: Cubic parabola
//input: double x = 1 {caption: x; category: Varibles}
//input: double a = 7.0 {caption: a; category: Parameters}
//input: double b = 0.1 {caption: b; category: Parameters}
//output: double result
export function cubicParabola(x: number, a : number, b : number): number {      
  return a * x * x * x + b;
}

//name: Quadric parabola
//input: double x = 1 {caption: x; category: Varibles}
//input: double a = 7.0 {caption: a; category: Parameters}
//input: double b = 0.1 {caption: b; category: Parameters}
//output: double result
//editor: Compute:RichFunctionViewEditor
export function quadricParabola(x: number, a : number, b : number): number {      
  return a * x * x * x * x + b;
}

//name: oneDfOutput
//input: double x1 = 1 {caption: x1; category: Varibles}
//input: double x2 = 2 {caption: x2; category: Varibles}
//input: double x3 = 3 {caption: x3; category: Varibles}
//input: double a = 7.0 {caption: a; category: Parameters}
//input: double b = 0.1 {caption: b; category: Parameters}
//output: dataframe res {caption: Solution; viewer: Grid}
//editor: Compute:RichFunctionViewEditor
export function oneDfOutput(x1: number, x2: number, x3: number, a : number, b : number): DG.DataFrame {
  return DG.DataFrame.fromColumns([
    DG.Column.fromList('double', 'Ishigami', [ishigamiFunc(x1, x2, x3, a, b)]),
    DG.Column.fromList('double', 'square', [parabola(x1, a, b)]),
    DG.Column.fromList('double', 'cubic', [cubicParabola(x2, a, b)]),
    DG.Column.fromList('double', 'quadric', [quadricParabola(x3, a, b)])
  ]);
}

//name: twoDfOutput
//input: double x1 = 1 {caption: x1; category: Varibles}
//input: double x2 = 2 {caption: x2; category: Varibles}
//input: double x3 = 3 {caption: x3; category: Varibles}
//input: double a = 7.0 {caption: a; category: Parameters}
//input: double b = 0.1 {caption: b; category: Parameters}
//output: dataframe df1 {caption: Ishigami; viewer: Grid}
//output: dataframe df2 {caption: Polynomials; viewer: Grid | Line chart}
//editor: Compute:RichFunctionViewEditor
export function twoDfOutput(x1: number, x2: number, x3: number, a : number, b : number): any {
  return {df1: DG.DataFrame.fromColumns([
    DG.Column.fromList('double', 'Ishigami', [ishigamiFunc(x1, x2, x3, a, b)])
  ]), 
  df2: DG.DataFrame.fromColumns([
    DG.Column.fromList('double', 'x', [x1, x2, x3]),
    DG.Column.fromList('double', 'funcs', [parabola(x1, a, b),
    cubicParabola(x2, a, b), quadricParabola(x3, a, b)])
  ])};
}

//name: try
export async function foo() {
  const someFunc: DG.Func = await grok.functions.eval('SA:twoDfOutput');
  console.log(someFunc);

  console.log('INPUTS:');
  for (const prop of someFunc.inputs)
    console.log(prop);

  console.log('OUTPUTS:');
  for (const prop of someFunc.outputs)
    console.log(prop);

  const someFuncCall = someFunc.prepare();
  console.log(someFuncCall);
}

//name: Perform Sensitivity Analysis
//input: string funcName {choices: ['ishigamiFunc', 'oneDfOutput', 'twoDfOutput']}
//input: int samplesCount = 1000
export async function performSA(funcName: string, samplesCount: number): Promise<void> {
  const func: DG.Func = await grok.functions.eval('SA:' + funcName);
  
  /*console.log(someFunc);

  console.log('INPUTS:');
  for (const prop of someFunc.inputs)
    console.log(prop);

  console.log('OUTPUTS:');
  for (const prop of someFunc.outputs)
    console.log(prop);

  const someFuncCall = someFunc.prepare(
    {
      'x1': 1,
      'x2': 2,
      'x3': 3,
      'a': 7,
      'b': 0.1
    }
  );*/

  const fixed1: FixedInputItem = {name: 'a', value: 7};
  const fixed2: FixedInputItem = {name: 'b', value: 0.1};

  const fixedInputs: Array<FixedInputItem> = [fixed1, fixed2];

  /*const inp1: VariedNumericalInputInfo = {
    name: 'x1', caption: undefined, type: DG.COLUMN_TYPE.FLOAT, min: 10, max: 20, column: undefined};
  
  const inp2: VariedNumericalInputInfo = {
    name: 'x2', caption: undefined, type: DG.COLUMN_TYPE.FLOAT, min: 30, max: 40, column: undefined};

  const inp3: VariedNumericalInputInfo = {
    name: 'x3', caption: undefined, type: DG.COLUMN_TYPE.FLOAT, min: -10, max: 0, column: undefined};*/

  const inp1: VariedNumericalInputInfo = {
    prop: func.inputs[0], min: 10, max: 20};
    
  const inp2: VariedNumericalInputInfo = {
    prop: func.inputs[1], min: 30, max: 40};
  
  const inp3: VariedNumericalInputInfo = {
    prop: func.inputs[2], min: -10, max: 0};

  const variedInputs: Array<VariedNumericalInputInfo> = [inp1, inp2, inp3];
  
  const sa = new VarianceBasedSenstivityAnalysis(func, fixedInputs, variedInputs, samplesCount);

  //await sa.perform();

  grok.shell.addTableView(await sa.perform());

  /*const runsCount = samplesCount * (variedInputs.length + 2);

  createVariedNumericalInputColumns(samplesCount, variedInputs);

  const funcCalls: Array<DG.FuncCall> = [];

  for (let i = 0; i < runsCount; ++i) {
    let inputs: any = {};

    for (const input of fixedInputs)
      inputs[input.name] = input.value;

    for (const input of variedInputs)
      inputs[input.name] = input.column?.get(i);
      
    funcCalls.push(someFunc.prepare(inputs));
  }

  for (const funcCall of funcCalls)
    await funcCall.call();*/



  /*const columns: DG.Column[] = [];

  for (const item of variedInputs)
    columns.push(item.column!);

  const result = DG.DataFrame.fromColumns(columns);

  grok.shell.addTableView(result);*/

  /*console.log(someFuncCall);

  await someFuncCall.call();

  console.log(getEachOutputItemInfo(someFuncCall));*/

  /*console.log(someFuncCall);
  
  console.log('=================');

  //console.log(someFuncCall.outputs['result']);

  //console.log(someFuncCall.outputParams.values());

  for (const el of [...someFuncCall.outputs])
    console.log(`${el[0]}  <-->  ${el[1]}`);

  console.log('=================');

  const names = [];

  for (const el of [...someFuncCall.inputs]) {
    console.log(`${el[0]}  <-->  ${el[1]}`);
    names.push(el[0]);
  }

  console.log(names);



  /*someFuncCall.outputs['...'] // !!!
  [...someFuncCall.outputParams] // array DG.FUNCALLPARAM name 
  [...lastCall.outputParams.values()] as DG.FuncCallParam[];*/
}

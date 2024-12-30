/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {getEachOutputItemInfo} from './output-tools';
import {VariedNumericalInputInfo, FixedInputItem} from './input-tools';

import {SobolAnalysis} from './sobol-sensitivity-analysis';

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
    DG.Column.fromList('double', 'quadric', [quadricParabola(x3, a, b)]),
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
    DG.Column.fromList('double', 'Ishigami', [ishigamiFunc(x1, x2, x3, a, b)]),
  ]),
  df2: DG.DataFrame.fromColumns([
    DG.Column.fromList('double', 'x', [x1, x2, x3]),
    DG.Column.fromList('double', 'funcs', [parabola(x1, a, b),
      cubicParabola(x2, a, b), quadricParabola(x3, a, b)]),
  ])};
}

//name: try
export async function foo() {
  const someFunc = DG.Func.byName('SA:twoDfOutput');
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

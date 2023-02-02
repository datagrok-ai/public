/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

import { callWasm } from '../wasm/callWasm';

//tags: init
export async function init() {
  await initTestExample();
}

//name: solveTestExample
//input: double initial = 0.0 {caption: initial; category: time, minutes}
//input: double final = 5.0 {caption: final; category: time, minutes}
//input: double step = 0.1 {caption: step; category: time, minutes}
//input: double _xInitial = 2.0 {units: x y.o.; caption: x; category: initial values}
//input: double _yInitial = 0.0 {units: y y.o.; caption: y; category: initial values}
//input: double _param1Val = 1.0 {units: param1 y.o.; caption: param1; category: parameters}
//input: double _param2Val = -1.0 {units: param2 y.o.; caption: param2; category: parameters}
//output: dataframe dfSolution {caption: Solution; viewer: Line chart(x: "t, time (minutes)", sharex: "true", multiAxis: "true", multiAxisLegendPosition: "RightCenter") | Grid(block: 100) }
//editor: Compute:RichFunctionViewEditor
export async function solveTestExample(initial, final, step, _xInitial, _yInitial, _param1Val, _param2Val) {

  let _tCount = Math.trunc((final - initial) / step) + 1;
  let _varsCount = 3;

  return callWasm(TestExample, 'solveTestExample', [initial, final, step, _xInitial, _yInitial, _param1Val, _param2Val, _tCount, _varsCount]);
}

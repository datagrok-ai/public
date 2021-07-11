/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {MultiPlotViewer} from './multiplot.ts';

export const _package = new DG.Package();

//name: test
export function test() {
  // grok.shell.info(_package.webRoot);
}

//name: MultiPlot
//description: Creates MultiPlot viewer
//tags: viewer
//output: viewer result
export function MultiPlot() {
  console.log('package.js     export function MultiPlot()');
  return new MultiPlotViewer();
}


/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {FitGridCellHandler} from './fit/fit-grid-cell-handler';
import {FitChartCellRenderer} from './fit/fit-renderer';
import {MultiCurveViewer} from './fit/multi-curve-viewer';
import {createDemoDataFrame, curveDemo} from './fit/fit-demo';


export const _package = new DG.Package();


//name: Fit
//tags: cellRenderer
//meta.cellType: fit
//meta.virtual: true
//output: grid_cell_renderer result
export function fitCellRenderer() {
  return new FitChartCellRenderer();
}

//name: MultiCurveViewer
//description: A viewer that superimposes multiple in-cell curves on one chart
//tags: viewer
//output: viewer result
export function _FitViewer() {
  return new MultiCurveViewer();
}

//tags: app
//name: Curves Demo
export function curveFitDemoApp() {
  grok.shell.addTableView(createDemoDataFrame(30, 5, 2));
}

//name: Curve fitting
//description: Curve fitting is the process of constructing a curve, or mathematical function, that has the best fit to a series of data points
//meta.demoPath: Curves | Curve fitting
export async function curveFitDemo() {
  await curveDemo();
}

//tags: autostart
export function _autoCurves(): void {
  DG.ObjectHandler.register(new FitGridCellHandler());
}

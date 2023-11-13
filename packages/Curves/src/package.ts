/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';

import {FitGridCellHandler} from './fit/fit-grid-cell-handler';
import {FitChartCellRenderer} from './fit/fit-renderer';
import {MultiCurveViewer} from './fit/multi-curve-viewer';
import {curveDemo} from './fit/fit-demo';


export const _package = new DG.Package();


//name: Fit
//tags: cellRenderer
//meta.cellType: fit
//meta.virtual: true
//output: grid_cell_renderer result
export function fitCellRenderer(): FitChartCellRenderer {
  return new FitChartCellRenderer();
}

//name: MultiCurveViewer
//description: A viewer that superimposes multiple in-cell curves on one chart
//tags: viewer
//output: viewer result
export function _FitViewer(): MultiCurveViewer {
  return new MultiCurveViewer();
}

//name: Curve fitting
//description: Curve fitting is the process of constructing a curve, or mathematical function, that has the best fit to a series of data points
//meta.demoPath: Curves | Curve fitting
//test: curveFitDemo() //wait: 2000
export async function curveFitDemo(): Promise<void> {
  await curveDemo();
}

//tags: init
export function _initCurves(): void {
  DG.ObjectHandler.register(new FitGridCellHandler());
}

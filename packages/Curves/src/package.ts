/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {FitGridCellHandler} from './fit/fit-grid-cell-handler';
import {FitChartCellRenderer} from './fit/fit-renderer';
import {FitViewer} from './fit/fit-viewer';
import {curveFitDemo} from './fit/fit-demo';


export const _package = new DG.Package();


//name: Fit
//tags: cellRenderer
//meta.cellType: fit
//meta.virtual: true
//output: grid_cell_renderer result
export function fitCellRenderer() {
  return new FitChartCellRenderer();
}

//name: FitViewer
//description: Creates a fit viewer
//tags: viewer
//output: viewer result
export function _FitViewer() {
  return new FitViewer();
}

//tags: app
//name: Curve Fit Demo
export function curveFitDemoApp() {
  curveFitDemo();
}

//tags: autostart
export function _autoCurves(): void {
  DG.ObjectHandler.register(new FitGridCellHandler());
}

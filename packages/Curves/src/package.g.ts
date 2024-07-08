import {FitChartCellRenderer} from './fit/fit-renderer';
import {MultiCurveViewer} from './fit/multi-curve-viewer';

//name: Fit
//tags: cellRenderer
//output: grid_cell_renderer renderer
//meta.cellType: fit
//meta.virtual: true
export function _FitChartCellRenderer() {
  return new FitChartCellRenderer();
}

//name: MultiCurveViewer
//description: A viewer that superimposes multiple in-cell curves on one chart
//tags: viewer
//output: viewer result
//meta.icon: icons/multi-curve-viewer.png
export function _MultiCurveViewer() {
  return new MultiCurveViewer();
}


/* The hand-written stand-in for the generated viewer factories (V5 / VP-18): one line per typed
   js-api static, each under its `I*Settings` with every option a value or a signal. */
import * as DG from 'datagrok-api/dg';
import type {Control} from '../../core/component.js';
import {propertyForm, ObjectForm} from '../forms/object-form.js';
import {viewerOf, Bindable, TableRef} from './viewer-control.js';

const api = globalThis as unknown as {grok_Viewer_Get_Look: (dart: unknown) => object};

const ofType = (type: string) => (df: DG.DataFrame, look: object) => DG.Viewer.fromType(type, df, look);

/**
 * @example
 * const x = signal('age');
 * const plot = viewers.scatterPlot(orders, {xColumnName: x, yColumnName: 'weight', markerMinSize: 4});
 * x.value = 'height';      // the plot follows; dragging the axis in the plot writes x back
 * plot.dispose();          // Dart detach through the kill-walk; the effects with it
 */
export const viewers = {
  grid: (t: TableRef, o?: Bindable<DG.IGridSettings>) => viewerOf(DG.Viewer.grid, t, o),
  histogram: (t: TableRef, o?: Bindable<DG.IHistogramSettings>) => viewerOf(DG.Viewer.histogram, t, o),
  barChart: (t: TableRef, o?: Bindable<DG.IBarChartSettings>) => viewerOf(DG.Viewer.barChart, t, o),
  confusionMatrix: (t: TableRef, o?: Bindable<DG.IConfusionMatrixSettings>) =>
    viewerOf(DG.Viewer.confusionMatrix, t, o),
  rocCurve: (t: TableRef, o?: Bindable<DG.IRocCurveSettings>) => viewerOf(DG.Viewer.rocCurve, t, o),
  heatMap: (t: TableRef, o?: Bindable<DG.IGridSettings>) => viewerOf(DG.Viewer.heatMap, t, o),
  boxPlot: (t: TableRef, o?: Bindable<DG.IBoxPlotSettings>) => viewerOf(DG.Viewer.boxPlot, t, o),
  /** js-api types the static `Viewer<IFiltersSettings>` although it builds a `FilterGroup` (P14). */
  filters: (t: TableRef, o?: Bindable<DG.IFiltersSettings>) =>
    viewerOf(DG.Viewer.filters, t, o) as DG.FilterGroup & Control,
  scatterPlot: (t: TableRef, o?: Bindable<DG.IScatterPlotSettings>) => viewerOf(DG.Viewer.scatterPlot, t, o),
  lineChart: (t: TableRef, o?: Bindable<DG.ILineChartSettings>) => viewerOf(DG.Viewer.lineChart, t, o),
  networkDiagram: (t: TableRef, o?: Bindable<DG.INetworkDiagramSettings>) => viewerOf(DG.Viewer.network, t, o),
  calendar: (t: TableRef, o?: Bindable<DG.ICalendarSettings>) => viewerOf(DG.Viewer.calendar, t, o),
  correlationPlot: (t: TableRef, o?: Bindable<DG.ICorrelationPlotSettings>) =>
    viewerOf(DG.Viewer.correlationPlot, t, o),
  densityPlot: (t: TableRef, o?: Bindable<DG.IDensityPlotSettings>) => viewerOf(DG.Viewer.densityPlot, t, o),
  form: (t: TableRef, o?: Bindable<DG.IFormSettings>) => viewerOf(DG.Viewer.form, t, o),
  markup: (t: TableRef, o?: Bindable<DG.IMarkupViewerSettings>) => viewerOf(DG.Viewer.markup, t, o),
  matrixPlot: (t: TableRef, o?: Bindable<DG.IMatrixPlotSettings>) => viewerOf(DG.Viewer.matrixPlot, t, o),
  pcPlot: (t: TableRef, o?: Bindable<DG.IPcPlotSettings>) => viewerOf(DG.Viewer.pcPlot, t, o),
  pieChart: (t: TableRef, o?: Bindable<DG.IPieChartSettings>) => viewerOf(DG.Viewer.pieChart, t, o),
  scatterPlot3d: (t: TableRef, o?: Bindable<DG.IScatterPlot3dSettings>) => viewerOf(DG.Viewer.scatterPlot3d, t, o),
  statistics: (t: TableRef, o?: Bindable<DG.IStatsViewerSettings>) => viewerOf(DG.Viewer.statistics, t, o),
  tileViewer: (t: TableRef, o?: Bindable<DG.ITileViewerSettings>) => viewerOf(DG.Viewer.tile, t, o),
  treeMap: (t: TableRef, o?: Bindable<DG.ITreeMapSettings>) => viewerOf(DG.Viewer.treeMap, t, o),
  trellisPlot: (t: TableRef, o?: Bindable<DG.ITrellisPlotSettings>) => viewerOf(DG.Viewer.trellisPlot, t, o),
  pivotTable: (t: TableRef, o?: Bindable<DG.IPivotViewerSettings>) => viewerOf(ofType('Pivot table'), t, o),
  byType: (type: string, t: TableRef, o?: Bindable<Record<string, unknown>>) => viewerOf(ofType(type), t, o),
};

/** The settings form of a viewer: its user-editable properties over its look — the target a
 * viewer property reads and writes, never the viewer or its dart handle (P6, VP-19). */
export function viewerSettings(viewer: DG.Viewer): ObjectForm {
  return propertyForm(viewer.getProperties().filter((p) => p.userEditable !== false),
    api.grok_Viewer_Get_Look(viewer.dart), {editors: 'auto'});
}

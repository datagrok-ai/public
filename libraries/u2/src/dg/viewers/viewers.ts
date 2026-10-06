/* The hand-written stand-in for the generated viewer factories (V5 / VP-18): one line per viewer
   type, built through `fromType` (the cached wrapper), each under its `I*Settings` with every
   option a value or a signal. */
import * as DG from 'datagrok-api/dg';
import {propertyForm, ObjectForm} from '../forms/object-form.js';
import {viewerOf, Bindable, TableRef} from './viewer-control.js';

const api = globalThis as unknown as {grok_Viewer_Get_Look: (dart: unknown) => object};

const ofType = (type: string) => (df: DG.DataFrame, look: object) => DG.Viewer.fromType(type, df, look);

/** Every factory builds through `fromType` — the cached wrapper `DG.toJs(v.dart)` and `DG.Widget.find`
 * answer; the typed js-api statics mint an uncached second one (VP-32). */
const typed = <S, V extends DG.Viewer = DG.Viewer>(type: string) => (t: TableRef, o?: Bindable<S>) =>
  viewerOf(ofType(type) as (df: DG.DataFrame, look: Partial<S>) => V, t, o);

/**
 * @example
 * const x = signal('age');
 * const plot = viewers.scatterPlot(orders, {xColumnName: x, yColumnName: 'weight', markerMinSize: 4});
 * x.value = 'height';      // the plot follows; dragging the axis in the plot writes x back
 * plot.dispose();          // Dart detach through the kill-walk; the effects with it
 */
export const viewers = {
  grid: typed<DG.IGridSettings, DG.Grid>('Grid'),
  histogram: typed<DG.IHistogramSettings, DG.HistogramViewer>('Histogram'),
  barChart: typed<DG.IBarChartSettings>('Bar chart'),
  confusionMatrix: typed<DG.IConfusionMatrixSettings, DG.ConfusionMatrix>('Confusion matrix'),
  rocCurve: typed<DG.IRocCurveSettings, DG.RocCurve>('Roc curve'),
  heatMap: typed<DG.IGridSettings, DG.Grid>('Heat map'),
  boxPlot: typed<DG.IBoxPlotSettings, DG.BoxPlot>('Box plot'),
  /** js-api types the static `Viewer<IFiltersSettings>` although it builds a `FilterGroup` (P14). */
  filters: typed<DG.IFiltersSettings, DG.FilterGroup>('Filters'),
  scatterPlot: typed<DG.IScatterPlotSettings, DG.ScatterPlotViewer>('Scatter plot'),
  lineChart: typed<DG.ILineChartSettings, DG.LineChartViewer>('Line chart'),
  networkDiagram: typed<DG.INetworkDiagramSettings, DG.NetworkDiagramViewer>('Network diagram'),
  calendar: typed<DG.ICalendarSettings, DG.CalendarViewer>('Calendar'),
  correlationPlot: typed<DG.ICorrelationPlotSettings>('Correlation plot'),
  densityPlot: typed<DG.IDensityPlotSettings, DG.DensityPlotViewer>('Density plot'),
  form: typed<DG.IFormSettings, DG.FormViewer>('Form'),
  markup: typed<DG.IMarkupViewerSettings>('Markup'),
  matrixPlot: typed<DG.IMatrixPlotSettings, DG.MatrixPlot>('Matrix plot'),
  pcPlot: typed<DG.IPcPlotSettings>('PC Plot'),
  pieChart: typed<DG.IPieChartSettings>('Pie chart'),
  scatterPlot3d: typed<DG.IScatterPlot3dSettings>('3d scatter plot'),
  statistics: typed<DG.IStatsViewerSettings>('Statistics'),
  tileViewer: typed<DG.ITileViewerSettings, DG.TileViewer>('Tile Viewer'),
  treeMap: typed<DG.ITreeMapSettings, DG.TreeMap>('Tree map'),
  trellisPlot: typed<DG.ITrellisPlotSettings>('Trellis plot'),
  pivotTable: typed<DG.IPivotViewerSettings>('Pivot table'),
  byType: (type: string, t: TableRef, o?: Bindable<Record<string, unknown>>) => viewerOf(ofType(type), t, o),
};

/** The settings form of a viewer: its user-editable properties over its look — the target a
 * viewer property reads and writes, never the viewer or its dart handle (P6, VP-19). */
export function viewerSettings(viewer: DG.Viewer): ObjectForm {
  return propertyForm(viewer.getProperties().filter((p) => p.userEditable !== false),
    api.grok_Viewer_Get_Look(viewer.dart), {editors: 'auto'});
}

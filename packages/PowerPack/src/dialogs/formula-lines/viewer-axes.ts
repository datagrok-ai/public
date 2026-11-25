import * as DG from 'datagrok-api/dg';

import {AxisColumns, AxisNames} from './types';

// How each viewer type names its axis columns: the look properties differ, and single-axis
// viewers (box plot, histogram, bar chart) have only the value axis.

/** The look properties that hold the X and the Y column of the [viewer]. */
export function axisProps(viewer: DG.Viewer): {x?: string, y?: string} {
  switch (viewer.type) {
  case DG.VIEWER.LINE_CHART: return {x: 'xColumnName', y: 'yColumnNames'};
  case DG.VIEWER.BOX_PLOT: return {y: 'valueColumnName'};
  case DG.VIEWER.HISTOGRAM: return {x: 'valueColumnName'};
  case DG.VIEWER.BAR_CHART: return barChartIsVertical(viewer) ? {y: 'valueColumnName'} : {x: 'valueColumnName'};
  default: return {x: 'xColumnName', y: 'yColumnName'};
  }
}

/** True for the viewers that have only the value axis. */
export function isSingleAxis(viewer: DG.Viewer): boolean {
  return viewer.type === DG.VIEWER.BOX_PLOT || viewer.type === DG.VIEWER.HISTOGRAM ||
    viewer.type === DG.VIEWER.BAR_CHART;
}

/** Resolves bar-chart `verticalMode` from the attached host viewer, mirroring the Dart-side getter. */
function barChartIsVertical(viewer: DG.Viewer): boolean {
  const orientation = viewer.getOptions().look.orientation;
  if (orientation === 'vertical')
    return true;
  if (orientation === 'horizontal')
    return false;
  // auto — derive from the host's canvas aspect ratio.
  const w = viewer.root.clientWidth;
  const h = viewer.root.clientHeight;
  return w > 0 && h / w > 5;
}

/** Pick the first numeric column other than `exclude`. Returns undefined if no such column. */
function pickFallbackNumericColumn(df: DG.DataFrame, exclude: string | undefined | null): string | undefined {
  return df.columns.toList().find((c) => c.isNumerical && c.name !== exclude)?.name;
}

/** The table the preview of the [src] viewer shows: a line chart previews its aggregated table. */
export function previewDataFrame(src: DG.Viewer): DG.DataFrame {
  return src.type === DG.VIEWER.LINE_CHART ?
    new DG.DataFrame(new DG.LineChartViewer(src.dart).activeFrame!) :
    src.dataFrame!;
}

/** The axes of the [src] viewer, in the columns of its preview [dataFrame]. */
export function sourceAxes(src: DG.Viewer, dataFrame: DG.DataFrame): AxisNames {
  if (src.type === DG.VIEWER.LINE_CHART) {
    const yCols: string[] = src.props.yColumnNames;
    const yCol = dataFrame.columns.toList()
      .find((col) => col.isNumerical && col.name !== src.props.xColumnName && yCols.some((n) => col.name.includes(n)));
    return {x: src.props.xColumnName, xMap: src.props.xMap, y: yCol === undefined ? src.props.xColumnName : yCol.name};
  }

  if (src.type === DG.VIEWER.TRELLIS_PLOT || src.type === DG.VIEWER.DENSITY_PLOT) {
    const look = src.type === DG.VIEWER.TRELLIS_PLOT ? src.getOptions().look['innerViewerLook'] : src.getOptions().look;
    return {x: look.xColumnName, y: look.yColumnName, xMap: look.xMap, yMap: look.yMap};
  }

  if (isSingleAxis(src)) {
    const look = src.getOptions().look;
    // The other axis falls back to another numeric column; the categorical axis of a bar chart
    // (`splitColumnName`) is undefined for the dialog.
    const fallback = pickFallbackNumericColumn(dataFrame, look.valueColumnName);
    return axisProps(src).y ?
      {x: fallback, y: look.valueColumnName, yMap: look.yMap} :
      {x: look.valueColumnName, y: fallback, xMap: look.xMap};
  }

  return {y: src.props.yColumnName, x: src.props.xColumnName, yMap: src.props.yMap, xMap: src.props.xMap};
}

/** The current columns pair of the [viewer]. Single-axis viewers report only their value axis;
 *  the other side is `null`. */
export function readAxisColumns(viewer: DG.Viewer, dataFrame: DG.DataFrame): AxisColumns {
  const props = axisProps(viewer);
  const isLineChart = viewer.type === DG.VIEWER.LINE_CHART;
  // Some viewers (e.g. density plot) don't declare xMap/yMap; the props proxy
  // throws on unknown names, so check first.
  const safeProp = (name: string): any =>
    viewer.props.hasProperty(name) ? viewer.props[name] : undefined;

  const xColName: string | undefined = props.x ? viewer.props[props.x] : undefined;
  let yColName: string | undefined = props.y && !isLineChart ? viewer.props[props.y] : undefined;
  if (isLineChart) {
    const yCols: string[] = viewer.props.yColumnNames;
    yColName = dataFrame.columns.toList().find((col) => yCols.some((n) => col.name.includes(n)))?.name;
  }

  // Bar chart's value axis renders aggregated values (count/sum/avg of
  // valueColumnName grouped by splitColumnName), not the raw column. Carry
  // the raw column for the formula's column reference (`${AGE} = …`) but
  // override stats so defaults like q2 fall inside the displayed range.
  const aggStats = viewer instanceof DG.BarChartViewer ? viewer.getAggregatedValueColumn()?.stats : undefined;

  return {
    y: yColName ? dataFrame.col(yColName) : null,
    x: xColName ? dataFrame.col(xColName) : null,
    yMap: props.y && !isLineChart ? safeProp('yMap') : undefined,
    xMap: props.x ? safeProp('xMap') : undefined,
    yStats: props.y ? aggStats : undefined,
    xStats: props.x ? aggStats : undefined,
  };
}

/** The look options that put the [names] columns on the axes of the [viewer].
 *  Single-axis viewers only honor the matching side. */
export function axesToOptions(viewer: DG.Viewer, dataFrame: DG.DataFrame, names: AxisNames): {[x: string]: any} {
  const props = axisProps(viewer);
  const options: {[x: string]: any} = {};

  if (props.y && names?.y && dataFrame.col(names.y)) {
    if (viewer.type === DG.VIEWER.LINE_CHART)
      options[props.y] = [names.y];
    else {
      options[props.y] = names.y;
      options['yMap'] = names.yMap;
    }
  }
  if (props.x && names?.x && dataFrame.col(names.x)) {
    options[props.x] = names.x;
    options['xMap'] = names.xMap;
  }

  return options;
}

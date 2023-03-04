/* eslint-disable valid-jsdoc */
/* eslint-disable no-multi-spaces */

import * as DG from 'datagrok-api/dg';
import {Property} from 'datagrok-api/src/entities';
import {TYPE} from 'datagrok-api/src/const';
import {
  fit,
  FitErrorModel,
  FitResult, fitResultProperties, IFitOptions,
  sigmoid
} from '@datagrok-libraries/statistics/src/parameter-estimation/fit-curve';

/**
 *  Datagrok curve fitting
 *
 * - Fitting: computing parameters of the specified function to best fit the data
 *   - Uses BFGS optimization algorithm (multi-threading for performance).
 *     For dose-response curves, we are typically fitting the sigmoid function
 *   - Ability to dynamically register custom fitting functions
 *   - Ability to get fitting performance characteristics (r-squared, classification, etc)
 * - Deep integration with the Datagrok grid
 *   - Either fitting on the fly, or using the supplied function + parameters
 *   - Multiple series in one cell
 *   - Ability to define chart, marker, or fitting options (such as fit function or marker color)
 *     on the column level, with the ability to override it on a grid cell or point level
 *   - Clicking a point in a chart within a grid makes it an outlier -> curve is re-fitted on the fly
 *   - Ability to specify a chart as "reference" so that it is shown on every other chart for comparison
 * - Ability to overlay curves from multiple grid cells (special viewer)
 * - Work with series stored in multiple formats (binary for performance, json for flexibility, etc)
*/

export const FIT_SEM_TYPE = 'fit';
export const FIT_CELL_TYPE = 'fit';
export const TAG_FIT = '.fit';

export const CONFIDENCE_INTERVAL_STROKE_COLOR = 'rgba(255,191,63,0.9)';
export const CONFIDENCE_INTERVAL_FILL_COLOR = 'rgba(255,238,204,0.3)';

export type FitMarkerType = 'circle' | 'triangle up' | 'triangle down' | 'cross';

/** A point in the fit series. Only x and y are required. Can override some fields defined in IFitSeriesOptions. */
export interface IFitPoint {
  x: number;
  y: number;
  outlier?: boolean;       // if true, renders as 'x' and gets ignored for curve fitting
  minY?: number;           // when defined, the marker renders as a candlestick with whiskers [minY, maxY]
  maxY?: number;           // when defined, the marker renders as a candlestick with whiskers [minY, maxY]
  marker?: FitMarkerType;  // overrides the marker type defined in IFitSeriesOptions
  color?: string;          // overrides the marker color defined in IFitSeriesOptions
}


/** Series options can be either applied globally on a column level, or partially overridden in particular series */
export interface IFitSeriesOptions {
  name?: string;
  fitFunction?: string;
  parameters?: number[];         // auto-fitting when not defined
  pointColor?: string;
  fitLineColor?: string;
  confidenceIntervalColor?: string;
  showFitLine?: boolean;
  showPoints?: boolean;
  showCurveConfidenceInterval?: boolean;   // show ribbon
  showBoxPlot?: boolean;      // if true, multiple values with the same X are rendered as a candlestick
  showConfidenceForX?: number;

  clickToToggle?: boolean;    // If true, clicking on the point toggles its outlier status and causes curve refitting
}


/** A series consists of points, has a name, and options.
 * If defined, seriesOptions are merged with {@link IFitChartData.seriesOptions} */
export interface IFitSeries extends IFitSeriesOptions {
  points: IFitPoint[];
}


/** Chart options. For fitted curves, this object is stored in the grid column tags and is
 * used by the renderer. */
export interface IFitChartOptions {
  minX? : number;
  minY? : number;
  maxX? : number;
  maxY? : number;

  xAxisName?: string;
  yAxisName?: string;

  logX?: boolean;
  logY?: boolean;

  showStatistics?: string[];
}


/** Data for the fit chart. */
export interface IFitChartData {
  chartOptions?: IFitChartOptions;
  seriesOptions?: IFitSeriesOptions;  // Default series options. Individual series can override it.
  series?: IFitSeries[];
}


/** Properties that describe {@link IFitChartOptions}. Useful for editing, initialization, transformations, etc. */
export const fitChartDataProperties: Property[] = [
  // Style and zoom
  Property.js('minX', TYPE.FLOAT, {description: 'Minimum value of the X axis', nullable: true}),
  Property.js('minY', TYPE.FLOAT, {description: 'Minimum value of the Y axis', nullable: true}),
  Property.js('maxX', TYPE.FLOAT, {description: 'Maximum value of the X axis', nullable: true}),
  Property.js('maxY', TYPE.FLOAT, {description: 'Maximum value of the Y axis', nullable: true}),
  Property.js('xAxisName', TYPE.STRING, {description:
    'Label to show on the X axis. If not specified, corresponding data column name is used', nullable: true}),
  Property.js('yAxisName', TYPE.STRING, {description:
    'Label to show on the Y axis. If not specified, corresponding data column name is used', nullable: true}),
  Property.js('logX', TYPE.BOOL, {defaultValue: false}),
  Property.js('logY', TYPE.BOOL, {defaultValue: false}),
  Property.js('showStatistics', TYPE.STRING_LIST, { choices: fitResultProperties.map(frp => frp.name) }),
];


/** Properties that describe {@link IFitSeriesOptions}. Useful for editing, initialization, transformations, etc. */
export const fitSeriesProperties: Property[] = [
  Property.js('name', TYPE.STRING),
  Property.js('fitFunction', TYPE.STRING,
    {category: 'Fitting', choices: ['sigmoid', 'linear'], defaultValue: 'sigmoid'}),
  Property.js('pointColor', TYPE.STRING,
    {category: 'Rendering', defaultValue: DG.Color.toHtml(DG.Color.scatterPlotMarker), nullable: true}),
  Property.js('fitLineColor', TYPE.STRING,
    {category: 'Rendering', defaultValue: DG.Color.toHtml(DG.Color.scatterPlotMarker), nullable: true}),
  Property.js('clickToToggle', TYPE.BOOL, {category: 'Fitting', description:
    'If true, clicking on the point toggles its outlier status and causes curve refitting', nullable: true}),
  Property.js('autoFit', TYPE.BOOL,
    {category: 'Fitting', description: 'Perform fitting on-the-fly', defaultValue: true}),
  Property.js('showFitLine', TYPE.BOOL,
    {category: 'Fitting', description: 'Whether the fit line should be rendered', defaultValue: true}),
  Property.js('showPoints', TYPE.BOOL,
    {category: 'Fitting', description: 'Whether points should be rendered', defaultValue: true}),
];


/** Creates new object with the default values specified in {@link properties} */
function createFromProperties(properties: Property[]): any {
  const o: any = {};
  for (const p of properties) {
    if (p.defaultValue != null)
      o[p.name] = p.defaultValue;
  }
  return o;
}


/** Merges properties of the two objects by iterating over the specified {@link properties}
 * and assigning properties from {@link source} to {@link target} only when
 * the property is not defined in target and is defined in source. */
function mergeProperties(properties: Property[],
  source: IFitChartOptions | IFitSeriesOptions | undefined, target: IFitChartOptions | IFitSeries | undefined) {
  if (!source || !target)
    return;

  for (const p of properties) {
    if (p.name !in target && p.name in source)
      p.set(target, p.get(source));
  }
}


function createDefaultChartData(): IFitChartData {
  return {
    chartOptions: createFromProperties(fitChartDataProperties),
    seriesOptions: createFromProperties(fitSeriesProperties)
  };
}


/** Returns existing, or creates new column default chart options. */
export function getColumnChartOptions(gridColumn: DG.GridColumn): IFitChartData {
  return gridColumn.tags[TAG_FIT] ??= createDefaultChartData();
}


export function getChartBounds(chartData: IFitChartData): DG.Rect {
  const o = chartData.chartOptions;
  if (o?.minX && o.minY && o.maxX && o.maxY)
    return new DG.Rect(o.minX, o.minY, o.maxX - o.minX, o.maxY - o.minY);
  if (!chartData.series?.length || chartData.series.length == 0) {
    return new DG.Rect(0, 0, 1, 1);
  } else {
    let bounds = getDataBounds(chartData.series[0].points);
    for (let i = 1; i < chartData.series!.length; i++)
      bounds = union(bounds, getDataBounds(chartData.series[i].points));
    return bounds;
  }
}


//TODO: move to DG.Rect
export function getDataBounds(points: IFitPoint[]): DG.Rect {
  let minX = points[0].x;
  let minY = points[0].y;
  let maxX = points[0].x;
  let maxY = points[0].y;

  for (let i = 1; i < points.length; i++) {
    minX = Math.min(minX, points[i].x);
    minY = Math.min(minY, points[i].y);
    maxX = Math.max(maxX, points[i].x);
    maxY = Math.max(maxY, points[i].y);
  }

  return new DG.Rect(minX, minY, maxX - minX, maxY - minY);
}


//TODO: move to DG.Rect
function fromPoints(x1: number, y1: number, x2: number, y2: number): DG.Rect {
  const minX = Math.min(x1, x2);
  const minY = Math.min(y1, y2);
  return new DG.Rect(minX, minY, Math.max(x1, x2) - minX, Math.max(y1, y2) - minY);
}


//TODO: move to DG.Rect
function union(r1: DG.Rect, r2: DG.Rect): DG.Rect {
  return fromPoints(
    Math.min(r1.left, r2.left),
    Math.min(r1.top, r2.top),
    Math.max(r1.right, r2.right),
    Math.max(r1.bottom, r2.bottom));
}


/** Fits the series data according to the series fitting settings */
export function fitSeries(series: IFitSeries, statistics: boolean = false): FitResult {
  //TODO: optimize the calculation of the initial parameters
  const dataBounds = getDataBounds(series.points);
  const ys = series.points.map((p) => p.y);
  const medY = ys[Math.floor(series.points.length / 2)];
  let xAtMedY = -1;
  for (let i = 0; i < series.points.length; i++) {
    if (series.points[i].y === medY) {
      xAtMedY = series.points[i].x;
      break;
    }
  }
  const initialParams = [dataBounds.bottom, 1.2, xAtMedY, dataBounds.top];

  return fit(
    {x: series.points.map((p) => p.x), y: series.points.map((p) => p.x)},
    initialParams, sigmoid, FitErrorModel.Constant, 0.05, statistics);
}

/** Returns a curve function, either using the pre-computed parameters
 * or by fitting on-the-fly */
export function getFittedCurve(series: IFitSeries): (x: number) => number {
  if (series.parameters)
    return (x) => sigmoid(series.parameters!, x);
  else
    return fitSeries(series).fittedCurve;
}

/** Returns confidence interval functions */
export function getConfidenceIntrevals(series: IFitSeries): {top: (x: number) => number, bottom: (x: number) => number} {
  const confTop = fitSeries(series).confidenceTop;
  const confBottom = fitSeries(series).confidenceBottom;
  return {top: confTop, bottom: confBottom};
}

/** Constructs {@link IFitChartData} from the grid cell, taking into account
 * chart and fit settings potentially defined on the dataframe and column level. */
export function getChartData(gridCell: DG.GridCell): IFitChartData {
  const cellChartData: IFitChartData = gridCell.cell?.column?.type == DG.TYPE.STRING ?
    (JSON.parse(gridCell.cell.value ?? '{}') ?? {}) : createDefaultChartData();
  cellChartData.series ??= [];

  const columnChartOptions = getColumnChartOptions(gridCell.gridColumn);

  // merge cell options with column options
  mergeProperties(fitChartDataProperties, columnChartOptions.chartOptions, cellChartData.chartOptions);
  for (const series of cellChartData.series)
    mergeProperties(fitSeriesProperties, columnChartOptions.seriesOptions, series);

  return cellChartData;
}


const sample: IFitChartData = {
  // chartOptions could be retrieved either from the column, or from the cell
  'chartOptions': {
    'minX': 0, 'minY': 0, 'maxX': 5, 'maxY': 10,
    'xAxisName': 'concentration',
    'yAxisName': 'activity',
    'logX': false,
    'logY': false,
  },
  // These options are used as default options for the series. They could be overridden in series.
  'seriesOptions': {
    'fitFunction': 'sigmoid',
    // parameters not specified -> auto-fitting by default
    'pointColor': 'blue',
    'fitLineColor': 'red',
    'clickToToggle': true,
    'showFitLine': true,
    'showPoints': true
  },
  'series': [
    {
      'fitFunction': 'sigmoid',
      // parameters specified -> use them, no autofitting
      'parameters': [1.86011e-07, -0.900, 103.748, -0.001],
      'points': [
        {'x': 0, 'y': 0},
        {'x': 1, 'y': 0.5},
        {'x': 2, 'y': 1},
        {'x': 3, 'y': 10, 'outlier': true},
        {'x': 4, 'y': 0}
      ]
    }
  ]
};

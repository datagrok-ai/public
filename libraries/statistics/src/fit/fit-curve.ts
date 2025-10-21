/* eslint-disable max-len */
/* eslint-disable no-multi-spaces */
import * as DG from 'datagrok-api/dg';
//@ts-ignore: no types
import * as jStat from 'jstat';

export const FitErrorModel = {
  CONSTANT: 'constant',
  PROPORTIONAL: 'proportional',
  COMBINED: 'combined',
};

export type FitParamBounds = {
  min?: number;
  max?: number;
};

/** Fit function description. Applies to custom user fit functions.
 * Requires JS arrow functions for the fit functions and initial parameters. */
export interface IFitFunctionDescription {
  name: string;
  function: string;
  getInitialParameters: string;
  parameterNames: string[];
}

export type FitCurve = {
  fittedCurve: (x: number) => number;
  parameters: Float32Array;
};

export type FitConfidenceIntervals = {
  confidenceTop: (x: number) => number;
  confidenceBottom: (x: number) => number;
};

export type FitStatistics = {
  rSquared?: number,
  auc?: number,
  interceptX?: number, // parameters[2]
  interceptY?: number, // fittedCurve[parameters[2]]
  slope?: number, // parameters[1]
  top?: number, // parameters[0]
  bottom?: number, // parameters[3]
};

export type FitInvertedFunctions = {
  inverted: (y: number) => number,
  invertedTop: (y: number) => number,
  invertedBottom: (y: number) => number,
};

/**
 *  Datagrok curve fitting
 *
 * - Fitting: computing parameters of the specified function to best fit the data
 *   - Uses BFGS optimization algorithm (multi-threading for performance).
 *     For dose-response curves, we are typically fitting the sigmoid function
 *   - Ability to dynamically register custom fitting functions
 *     - Automatic fit function determination
 *     - Caching of custom fitting functions
 *   - Ability to get fitting performance characteristics (r-squared, classification, etc)
 * - Deep integration with the Datagrok grid
 *   - Either fitting on the fly, or using the supplied function + parameters
 *   - Multiple series in one cell
 *   - Candlesticks, confidence intervals, and droplines drawing
 *   - Ability to define chart, marker, or fitting options (such as fit function or marker color)
 *     on the column level, with the ability to override it on a grid cell or point level
 *   - Clicking a point in a chart within a grid makes it an outlier -> curve is re-fitted on the fly
 *   - Ability to specify a chart as "reference" so that it is shown on every other chart for comparison
 * - Ability to overlay curves from multiple grid cells (special viewer)
 * - Work with series stored in multiple formats (binary for performance, json for flexibility, etc)
*/

export const DROPLINES = ['IC50'];

export type FitMarkerType = 'asterisk' | 'circle' | 'cross border' | 'diamond' | 'square' | 'star' | 'triangle bottom' |
  'triangle left' | 'triangle right' | 'triangle top';

export type FitOutlierMarkerType = 'asterisk' | 'circle' | 'cross border' | 'diamond' | 'outlier' | 'square' | 'star' |
  'triangle bottom' | 'triangle left' | 'triangle right' | 'triangle top';

export type FitLineStyle = 'solid' | 'dotted' | 'dashed' | 'dashdotted';

export type FitErrorModelType = 'constant' | 'proportional' | 'combined';

/** A point in the fit series. Only x and y are required. Can override some fields defined in IFitSeriesOptions. */
export interface IFitPoint {
  x: number;
  y: number;
  outlier?: boolean;       // if true, renders as 'x' and gets ignored for curve fitting
  color?: string;          // overrides the marker color defined in IFitSeriesOptions
  outlierColor?: string;   // overrides the outlier color defined in IFitSeriesOptions
  marker?: FitMarkerType;  // overrides the marker type defined in IFitSeriesOptions
  outlierMarker?: FitOutlierMarkerType; // overrides the outlier marker type defined in IFitSeriesOptions
  size?: number;           // overrides the default marker size
  stdev?: number;          // when defined, renders an error bar candlestick
  // minY?: number;           // when defined, the marker renders as a candlestick with whiskers [minY, maxY]
  // maxY?: number;           // when defined, the marker renders as a candlestick with whiskers [minY, maxY]
  meta?: any;           // any additional data
}

/** A series consists of points, has a name, and options.
 * If defined, seriesOptions are merged with {@link IFitChartData.seriesOptions} */
export interface IFitSeries extends IFitSeriesOptions {
  points: IFitPoint[];
}

/** Chart labels options. Controls the chart labels. */
export interface IFitChartLabelOptions {
  visible: boolean;         // if true, renders the label on the plot
  color: string;            // defines the label color
  name: string;             // defines the label name
}

/** Chart options. For fitted curves, this object is stored in the grid column tags and is used by the renderer. */
export interface IFitChartOptions {
  minX?: number;
  minY?: number;
  maxX?: number;
  maxY?: number;

  title?: string;            // defines the plot title. If the plot size is enough, will render it.
  xAxisName?: string;        // defines the x axis name. If the plot size is enough, will render it.
  yAxisName?: string;        // defines the Y axis name. If the plot size is enough, will render it.

  logX?: boolean;            // defines whether the x data should be logarithmic or not
  logY?: boolean;            // defines whether the y data should be logarithmic or not

  allowXZeroes?: boolean;    // defines whether x zeroes allowed for logarithmic data or not
  mergeSeries?: boolean;     // defines whether to merge series or not

  showColumnLabel?: boolean; // defines whether to show the column label in the legend or not
  showStatistics?: string[]; // defines the statistics that would be shown on the plot
  labelOptions?: IFitChartLabelOptions[]; // controls the plot labels
  useAuxLegendNames?: boolean; // if true, uses aux legend names instead of series names
}

/** Data for the fit chart. */
export interface IFitChartData {
  chartOptions?: IFitChartOptions;
  seriesOptions?: IFitSeriesOptions;  // Default series options. Individual series can override it.
  series?: IFitSeries[];
}

/** Class that implements {@link IFitChartData} interface */
export class FitChartData implements IFitChartData {
  chartOptions: IFitChartOptions = {};
  seriesOptions: IFitSeriesOptions = {};  // Default series options. Individual series can override it.
  series: IFitSeries[] = [];
}

/** Series options can be either applied globally on a column level, or partially overridden in particular series */
export interface IFitSeriesOptions {
  [key: string]: any;                   // allows getting data by key
  name?: string;                        // controls the series name
  fitFunction?: string | IFitFunctionDescription; // controls the series fit function
  parameters?: number[];                // controls the series parameters, auto-fitting when not defined
  parameterBounds?: FitParamBounds[];   // defines the acceptable range of each parameter, which is taken into account during the fitting. See also `parameters`.
  markerType?: FitMarkerType;           // defines the series marker type
  outlierMarkerType?: FitOutlierMarkerType;  // defines the series outlier marker type
  lineStyle?: FitLineStyle;             // defines the series line style
  pointColor?: string;                  // overrides the standardized series point color
  fitLineColor?: string;                // overrides the standardized series fit line color
  confidenceIntervalColor?: string;     // overrides the standardized series confidence interval color
  outlierColor?: string;                // overrides the standardized series outlier color
  connectDots?: boolean;                // defines whether to connect the points with lines or not. If true and showFitLine is false - fitting is disabled - otherwise, it will be rendered accordingly to the parameter value.
  showFitLine?: boolean;                // defines whether to show the fit line or not
  showPoints?: string;                  // defines the data display mode
  showOutliers?: boolean;               // defines whether to show the outliers or not
  showCurveConfidenceInterval?: boolean;    // defines whether to show the confidence intervals or not
  errorModel?: FitErrorModelType;       // defines the series error model
  clickToToggle?: boolean;    // if true, clicking on the point toggles its outlier status and causes curve refitting
  labels?: {[key: string]: string | number | boolean}; // controlled by IFitChartData labelOptions, shows labels
  droplines?: string[];                 // defines the droplines that would be shown on the plot (IC50)
  columnName?: string;                  // defines the column name where the series is stored
  auxLegendName?: string;               // defines the auxiliary legend name for the series
}
// TODO: show labels in property panel if present, color by default from series


/** Properties that describe {@link FitStatistics}. Useful for editing, initialization, transformations, etc. */
export const statisticsProperties: DG.Property[] = [
  DG.Property.js('rSquared', DG.TYPE.FLOAT, {userEditable: false}),
  DG.Property.js('auc', DG.TYPE.FLOAT, {userEditable: false, friendlyName: 'AUC'}),
  DG.Property.js('interceptY', DG.TYPE.FLOAT, {userEditable: false}),
  DG.Property.js('interceptX', DG.TYPE.FLOAT, {userEditable: false}),
  DG.Property.js('slope', DG.TYPE.FLOAT, {userEditable: false}),
  DG.Property.js('top', DG.TYPE.FLOAT, {userEditable: false, friendlyName: 'Max Y'}),
  DG.Property.js('bottom', DG.TYPE.FLOAT, {userEditable: false, friendlyName: 'Min Y'}),
];

/** Properties that describe {@link IFitChartOptions}. Useful for editing, initialization, transformations, etc. */
export const fitChartDataProperties: DG.Property[] = [
  // Style and zoom
  // remove unnecessary tooltips
  DG.Property.js('minX', DG.TYPE.FLOAT, {nullable: true}),
  DG.Property.js('minY', DG.TYPE.FLOAT, {nullable: true}),
  DG.Property.js('maxX', DG.TYPE.FLOAT, {nullable: true}),
  DG.Property.js('maxY', DG.TYPE.FLOAT, {nullable: true}),
  DG.Property.js('title', DG.TYPE.STRING, {nullable: true}),
  DG.Property.js('xAxisName', DG.TYPE.STRING, {description:
    'Label to show on the X axis. If not specified, corresponding data column name is used', nullable: true}),
  DG.Property.js('yAxisName', DG.TYPE.STRING, {description:
    'Label to show on the Y axis. If not specified, corresponding data column name is used', nullable: true}),
  DG.Property.js('logX', DG.TYPE.BOOL, {description: 'Whether the X axis should be logarithmic', defaultValue: false}),
  DG.Property.js('logY', DG.TYPE.BOOL, {description: 'Whether the Y axis should be logarithmic', defaultValue: false}),
  DG.Property.js('allowXZeroes', DG.TYPE.BOOL, {description: 'Whether x zeroes allowed for logarithmic data or not', defaultValue: true}),
  DG.Property.js('mergeSeries', DG.TYPE.BOOL, {description: 'Merges all series within one cell into one series', defaultValue: false}),
  DG.Property.js('showColumnLabel', DG.TYPE.BOOL, {description: 'Whether to show the column label in the legend or not', defaultValue: false}),
  DG.Property.js('showStatistics', DG.TYPE.STRING_LIST, {description: 'Whether specific statistics should be rendered',
    choices: statisticsProperties.map((frp) => frp.name), inputType: 'MultiChoice',
    //@ts-ignore
    friendlyName: 'Statistics'}),
];

export const FIT_FUNCTION_SIGMOID = 'sigmoid';
export const FIT_FUNCTION_LINEAR = 'linear';
export const FIT_FUNCTION_LOG_LINEAR = 'log-linear';
export const FIT_FUNCTION_EXPONENTIAL = 'exponential';
export const FIT_FUNCTION_4PL_REGRESSION = '4pl-regression';
export const FIT_FUNCTION_4PL_DOSE_RESPONSE = '4pl-dose-response';
export const FIT_JS_FUNCTION = 'js-function';

export const FIT_STATS_RSQUARED = 'rSquared';
export const FIT_STATS_AUC = 'auc';

export function getFittedCurve(curveFunction: (params: Float32Array, x: number) => number, paramValues: Float32Array):
 (x: number) => number {
  return (x: number) => {
    return curveFunction(paramValues, x);
  };
}

// TODO: for linear - slope - A, interceptY - B
export function getStatistics(data: {x: number[], y: number[]}, paramValues: Float32Array,
  curveFunction: (params: Float32Array, x: number) => number, statistics: boolean = true): FitStatistics {
  const fittedCurve = getFittedCurve(curveFunction, paramValues);

  return {
    rSquared: statistics ? getDetCoeff(fittedCurve, data) : undefined,
    auc: statistics ? getAuc(fittedCurve, data) : undefined,
    interceptX: paramValues[2],
    interceptY: fittedCurve(paramValues[2]),
    slope: paramValues[1],
    top: paramValues[0],
    bottom: paramValues[3],
  };
}

export function getInvertedFunctions(data: {x: number[], y: number[]}, paramValues: number[],
  confidenceLevel: number = 0.05, statistics: boolean = true): FitInvertedFunctions | null {
  const studentQ = jStat.studentt.inv(1 - confidenceLevel / 2, data.x.length - paramValues.length);

  let inv: (y: number) => number = (y: number) => {
    return 0;
  };
  let invTop: (y: number) => number = (y: number) => {
    return 0;
  };
  let invBottom: (y: number) => number = (y: number) => {
    return 0;
  };

  if (statistics) {
    inv = (y: number) => {
      //should check if more than bottom and less than top
      return paramValues[2] / Math.pow((paramValues[0] - y) / (y - paramValues[3]), 1 / paramValues[1]);
    };

    const error = getInvError(inv, data);

    invTop = (y: number) => {
      const value = inv(y);
      return value + studentQ * error / Math.sqrt(data.y.length);
    };

    invBottom = (y: number) => {
      const value = inv(y);
      return value - studentQ * error / Math.sqrt(data.y.length);
    };

    return {
      inverted: inv,
      invertedTop: invTop,
      invertedBottom: invBottom,
    };
  }

  return null;
}

export function sigmoid(params: Float32Array, x: number): number {
  const A = params[0];
  const B = params[1];
  const C = params[2];
  const D = params[3];
  return (D + (A - D) / (1 + Math.pow(10, (x - C) * B)));
}

export function linear(params: Float32Array, x: number): number {
  const A = params[0];
  const B = params[1];
  return A * x + B;
}

export function logLinear(params: Float32Array, x: number): number {
  const A = params[0];
  const B = params[1];
  return A * Math.log(x + 1) + B;
}

export function exponential(params: Float32Array, x: number): number {
  const A = params[0];
  const B = params[1];
  return A * Math.exp(x * B);
}

export function fourPLRegression(params: Float32Array, x: number): number {
  const A = params[0]; // max
  const B = params[1]; // hill
  const C = params[2]; // IC50 / inflection point
  const D = params[3]; // min
  return D + (A - D) / (1 + Math.pow(x / C, B));
}

export function fourPLDoseResponse(params: Float32Array, x: number): number {
  const min = params[3];
  const max = params[0];
  const hill = params[1];
  const ic50 = params[2];
  // be ware, IC50 that comes here is already log10 transformed
  return min + (max - min) / (1 + Math.pow(10, hill * (ic50 - x)));
}

export function getAuc(fittedCurve: (x: number) => number, data: {x: number[], y: number[]}): number {
  let auc = 0;
  const min = Math.min(...data.x);
  const max = Math.max(...data.x);
  const integrationStep = (max - min) / 1000;
  for (let x = min; x < max; x+= integrationStep)
    auc += integrationStep * fittedCurve(x);

  return auc;
}

export function getDetCoeff(fittedCurve: (x: number) => number, data: {x: number[], y: number[]}): number {
  let ssRes = 0;
  let ssTot = 0;
  const yMean = jStat.mean(data.y);
  for (let i = 0; i < data.x.length; i++) {
    ssRes += Math.pow(data.y[i] - fittedCurve(data.x[i]), 2);
    ssTot += Math.pow(data.y[i] - yMean, 2);
  }

  return 1 - ssRes / ssTot;
}

function getInvError(targetFunc: (y: number) => number, data: {y: number[], x: number[]}): number {
  let sigma = 0;
  let sigmaSq = 0;
  const residuesSquares = new Float32Array(data.y.length);
  for (let i = 0; i < data.y.length; i++) {
    const obs = data.x[i];
    const pred = targetFunc(data.y[i]);
    residuesSquares[i] = Math.pow(obs - pred, 2);
  }

  for (let i = 0; i < residuesSquares.length; i++)
    sigmaSq += residuesSquares[i];

  sigmaSq /= residuesSquares.length;
  sigma = Math.sqrt(sigmaSq);

  return sigma;
}

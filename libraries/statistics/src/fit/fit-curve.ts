/* eslint-disable max-len */
/* eslint-disable no-multi-spaces */
import * as DG from 'datagrok-api/dg';

import {limitedMemoryBFGS} from '../../lbfgs/lbfgs';
//@ts-ignore: no types
import * as jStat from 'jstat';


type Optimizable = {
  getValue: (parameters: number[]) => number,
  getGradient: (parameters: number[], gradient: number[]) => number[],
}

type Likelihood = {
  value: number,
  const: number,
  mult: number
};

type ObjectiveFunction = (targetFunc: (params: number[], x: number) => number,
  data: {x: number[], y: number[]}, params: number[]) => Likelihood;

export const FitErrorModel = {
  CONSTANT: 'constant',
  PROPORTIONAL: 'proportional',
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
  parameters: number[];
};

export type FitConfidenceIntervals = {
  confidenceTop: (x: number) => number;
  confidenceBottom: (x: number) => number;
};

export type FitStatistics = {
  rSquared?: number,
  auc?: number,
  interceptX: number, // parameters[2]
  interceptY: number, // fittedCurve[parameters[2]]
  slope: number, // parameters[1]
  top: number, // parameters[0]
  bottom: number, // parameters[3]
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

export const FIT_SEM_TYPE = 'fit';
export const FIT_CELL_TYPE = 'fit';
export const TAG_FIT = '.fit';

export const CONFIDENCE_INTERVAL_STROKE_COLOR = 'rgba(255,191,63,0.4)';
export const CONFIDENCE_INTERVAL_FILL_COLOR = 'rgba(255,238,204,0.3)';

export const CURVE_CONFIDENCE_INTERVAL_BOUNDS = {
  TOP: 'top',
  BOTTOM: 'bottom',
};

export const DROPLINES = ['IC50'];

export type FitMarkerType = 'asterisk' | 'circle' | 'cross border' | 'diamond' | 'square' | 'star' | 'triangle bottom' |
  'triangle left' | 'triangle right' | 'triangle top';

export type FitLineStyle = 'solid' | 'dotted' | 'dashed' | 'dashdotted';

export type FitErrorModelType = 'constant' | 'proportional';

/** A point in the fit series. Only x and y are required. Can override some fields defined in IFitSeriesOptions. */
export interface IFitPoint {
  x: number;
  y: number;
  outlier?: boolean;       // if true, renders as 'x' and gets ignored for curve fitting
  color?: string;          // overrides the marker color defined in IFitSeriesOptions
  outlierColor?: string;   // overrides the outlier color defined in IFitSeriesOptions
  marker?: FitMarkerType;  // overrides the marker type defined in IFitSeriesOptions
  size?: number;           // overrides the default marker size
  stdev?: number;          // when defined, renders an error bar candlestick
  // minY?: number;           // when defined, the marker renders as a candlestick with whiskers [minY, maxY]
  // maxY?: number;           // when defined, the marker renders as a candlestick with whiskers [minY, maxY]
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

  showStatistics?: string[]; // defines the statistics that would be shown on the plot
  labelOptions?: IFitChartLabelOptions[]; // controls the plot labels
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
  name?: string;                        // controls the series name
  fitFunction?: string | IFitFunctionDescription; // controls the series fit function
  parameters?: number[];                // controls the series parameters, auto-fitting when not defined
  parameterBounds?: FitParamBounds[];   // defines the acceptable range of each parameter, which is taken into account during the fitting. See also `parameters`.
  markerType?: FitMarkerType;           // defines the series marker type
  lineStyle?: FitLineStyle;             // defines the series line style
  pointColor?: string;                  // overrides the standardized series point color
  fitLineColor?: string;                // overrides the standardized series fit line color
  confidenceIntervalColor?: string;     // overrides the standardized series confidence interval color
  outlierColor?: string;                // overrides the standardized series outlier color
  connectDots?: boolean;                // defines whether to connect the points with lines or not. If true and showFitLine is false - fitting is disabled - otherwise, it will be rendered accordingly to the parameter value.
  showFitLine?: boolean;                // defines whether to show the fit line or not
  showPoints?: string;                  // defines the data display mode
  showCurveConfidenceInterval?: boolean;    // defines whether to show the confidence intervals or not
  errorModel?: FitErrorModelType;       // defines the series error model
  clickToToggle?: boolean;    // if true, clicking on the point toggles its outlier status and causes curve refitting
  labels?: {[key: string]: string | number | boolean}; // controlled by IFitChartData labelOptions, shows labels
  droplines?: string[];                 // defines the droplines that would be shown on the plot (IC50)
}
// TODO: show labels in property panel if present, color by default from series


/** Properties that describe {@link FitStatistics}. Useful for editing, initialization, transformations, etc. */
export const statisticsProperties: DG.Property[] = [
  DG.Property.js('rSquared', DG.TYPE.FLOAT, {userEditable: false}),
  DG.Property.js('auc', DG.TYPE.FLOAT, {userEditable: false}),
  DG.Property.js('interceptY', DG.TYPE.FLOAT, {userEditable: false}),
  DG.Property.js('interceptX', DG.TYPE.FLOAT, {userEditable: false}),
  DG.Property.js('slope', DG.TYPE.FLOAT, {userEditable: false}),
  DG.Property.js('top', DG.TYPE.FLOAT, {userEditable: false}),
  DG.Property.js('bottom', DG.TYPE.FLOAT, {userEditable: false}),
];

/** Properties that describe {@link IFitChartOptions}. Useful for editing, initialization, transformations, etc. */
export const fitChartDataProperties: DG.Property[] = [
  // Style and zoom
  DG.Property.js('minX', DG.TYPE.FLOAT, {description: 'Minimum value of the X axis', nullable: true}),
  DG.Property.js('minY', DG.TYPE.FLOAT, {description: 'Minimum value of the Y axis', nullable: true}),
  DG.Property.js('maxX', DG.TYPE.FLOAT, {description: 'Maximum value of the X axis', nullable: true}),
  DG.Property.js('maxY', DG.TYPE.FLOAT, {description: 'Maximum value of the Y axis', nullable: true}),
  DG.Property.js('title', DG.TYPE.STRING, {nullable: true}),
  DG.Property.js('xAxisName', DG.TYPE.STRING, {description:
    'Label to show on the X axis. If not specified, corresponding data column name is used', nullable: true}),
  DG.Property.js('yAxisName', DG.TYPE.STRING, {description:
    'Label to show on the Y axis. If not specified, corresponding data column name is used', nullable: true}),
  DG.Property.js('logX', DG.TYPE.BOOL, {description: 'Whether the X axis should be logarithmic', defaultValue: false}),
  DG.Property.js('logY', DG.TYPE.BOOL, {description: 'Whether the Y axis should be logarithmic', defaultValue: false}),
  DG.Property.js('allowXZeroes', DG.TYPE.BOOL, {description: 'Whether x zeroes allowed for logarithmic data or not', defaultValue: true}),
  DG.Property.js('mergeSeries', DG.TYPE.BOOL, {description: 'Whether to merge series or not', defaultValue: false}),
  DG.Property.js('showStatistics', DG.TYPE.STRING_LIST, {description: 'Whether specific statistics should be rendered',
    choices: statisticsProperties.map((frp) => frp.name), inputType: 'MultiChoice',
    //@ts-ignore
    friendlyName: 'Statistics'}),
];

/** Properties that describe {@link IFitSeriesOptions}. Useful for editing, initialization, transformations, etc. */
export const fitSeriesProperties: DG.Property[] = [
  DG.Property.js('fitFunction', DG.TYPE.STRING,
    {category: 'Fitting', choices: ['sigmoid', 'linear'], defaultValue: 'sigmoid'}),
  DG.Property.js('pointColor', DG.TYPE.STRING,
    {category: 'Rendering', nullable: true, inputType: 'Color'}),
  DG.Property.js('fitLineColor', DG.TYPE.STRING,
    {category: 'Rendering', nullable: true, inputType: 'Color'}),
  DG.Property.js('outlierColor', DG.TYPE.STRING,
    {category: 'Rendering', nullable: true, inputType: 'Color'}),
  DG.Property.js('errorModel', DG.TYPE.STRING, {category: 'Fitting', defaultValue: 'constant',
    choices: ['constant', 'proportional'], nullable: false}),
  DG.Property.js('connectDots', DG.TYPE.BOOL, {category: 'Fitting', defaultValue: false}),
  DG.Property.js('clickToToggle', DG.TYPE.BOOL, {category: 'Fitting', description:
    'Click on a point to mark it as outlier and refit', nullable: true, defaultValue: false}),
  DG.Property.js('showFitLine', DG.TYPE.BOOL, {category: 'Fitting', defaultValue: true}),
  DG.Property.js('showPoints', DG.TYPE.STRING, // rewrite description
    {category: 'Fitting', description: 'Whether points/candlesticks/none should be rendered',
      defaultValue: 'points', choices: ['points', 'candlesticks', 'both']}),
  DG.Property.js('showCurveConfidenceInterval', DG.TYPE.BOOL,
    {category: 'Fitting', description: 'Whether confidence intervals should be rendered', defaultValue: false,
      //@ts-ignore
      friendlyName: 'Confidence Interval'}),
  DG.Property.js('markerType', DG.TYPE.STRING, {category: 'Rendering', defaultValue: 'circle',
    choices: ['asterisk', 'circle', 'cross border', 'diamond', 'square', 'star',
      'triangle bottom', 'triangle left', 'triangle right', 'triangle top'], nullable: false}),
  DG.Property.js('lineStyle', DG.TYPE.STRING, {category: 'Rendering', defaultValue: 'solid',
    choices: ['solid', 'dotted', 'dashed', 'dashdotted'], nullable: false}),
  DG.Property.js('droplines', DG.TYPE.STRING_LIST, {description: 'Whether specific droplines should be rendered',
    choices: DROPLINES, inputType: 'MultiChoice'}),
];

export const FIT_FUNCTION_SIGMOID = 'sigmoid';
export const FIT_FUNCTION_LINEAR = 'linear';

export const FIT_STATS_RSQUARED = 'rSquared';
export const FIT_STATS_AUC = 'auc';


// TODO?: add method to return parameters - get parameters from fit function
/** Class for the fit functions */
export abstract class FitFunction {
  abstract get name(): string;
  abstract get parameterNames(): string[];
  abstract y(params: number[], x: number): number;
  abstract getInitialParameters(x: number[], y: number[]): number[];
}

/** Class that implements the linear function */
export class LinearFunction extends FitFunction {
  get name(): string {
    return FIT_FUNCTION_LINEAR;
  }

  get parameterNames(): string[] {
    return ['Slope', 'Intercept'];
  }

  y(params: number[], x: number): number {
    return linear(params, x);
  }

  getInitialParameters(x: number[], y: number[]): number[] {
    let minIndex = 0;
    let maxIndex = 0;
    for (let i = 1; i < x.length; i++) {
      if (x[i] < x[minIndex])
        minIndex = i;
      if (x[i] > x[maxIndex])
        maxIndex = i;
    }

    const deltaX = x[maxIndex] - x[minIndex];
    const deltaY = y[maxIndex] - y[minIndex];
    const A = deltaY / deltaX;
    const B = y[maxIndex] - A * x[maxIndex];
    return [A, B];
  }
}

/** Class that implements the sigmoid function */
export class SigmoidFunction extends FitFunction {
  get name(): string {
    return FIT_FUNCTION_SIGMOID;
  }

  get parameterNames(): string[] {
    return ['Top', 'Bottom', 'Slope', 'IC50'];
  }

  y(params: number[], x: number): number {
    return sigmoid(params, x);
  }

  getInitialParameters(x: number[], y: number[]): number[] {
    const dataBounds = DG.Rect.fromXYArrays(x, y);
    const medY = (dataBounds.bottom - dataBounds.top) / 2 + dataBounds.top;
    let maxYInterval = dataBounds.bottom - dataBounds.top;
    let nearestXIndex = 0;
    for (let i = 0; i < x.length; i++) {
      const currentInterval = Math.abs(y[i] - medY);
      if (currentInterval < maxYInterval) {
        maxYInterval = currentInterval;
        nearestXIndex = i;
      }
    }
    const xAtMedY = x[nearestXIndex];
    const slope = y[0] > y[y.length - 1] ? 1 : -1;

    // params are: [max, tan, IC50, min]
    return [dataBounds.bottom, slope, xAtMedY, dataBounds.top];
  }
}

/** Class that implements user JS functions */
export class JsFunction extends FitFunction {
  private _name: string;
  private _parameterNames: string[];

  constructor(name: string, yFunc: (params: number[], x: number) => number,
    getInitParamsFunc: (x: number[], y: number[]) => number[], parameterNames: string[]) {
    super();

    this._name = name;
    this._parameterNames = parameterNames;

    this.y = yFunc;
    this.getInitialParameters = getInitParamsFunc;
  }

  get name(): string {
    return this._name;
  }

  get parameterNames(): string[] {
    return this._parameterNames;
  }

  y(params: number[], x: number): number {
    throw new Error('Not implemented');
  }

  getInitialParameters(x: number[], y: number[]): number[] {
    throw new Error('Not implemented');
  }
}

// Object with fit functions
export const fitFunctions: {[index: string]: FitFunction} = {
  'linear': new LinearFunction(),
  'sigmoid': new SigmoidFunction(),
};

export interface IFitOptions {
  errorModel: FitErrorModelType;
  confidenceLevel: number;
  statistics: boolean;
}


function createObjectiveFunction(errorModel: FitErrorModelType): ObjectiveFunction {
  let of: ObjectiveFunction;

  switch (errorModel) {
  case FitErrorModel.CONSTANT:
    of = objectiveNormalConstant;
    break;
  case FitErrorModel.PROPORTIONAL:
    of = objectiveNormalProportional;
    break;
  default:
    of = objectiveNormalConstant;
    break;
  }

  return of;
}

function createOptimizable(data: {x: number[], y: number[]}, curveFunction: (params: number[], x: number) => number,
  of: ObjectiveFunction, fixed: number[]): Optimizable {
  return {
    getValue: (parameters: number[]) => {
      return of(curveFunction, data, parameters).value;
    },
    getGradient: (parameters: number[], gradient: number[]) => {
      for (let i = 0; i < parameters.length; i++)
        gradient[i] = fixed.includes(i) ? 0 : getObjectiveDerivative(of, curveFunction, data, parameters, i);

      return gradient;
    },
  };
}

export function getOrCreateFitFunction(seriesFitFunc: string | IFitFunctionDescription): FitFunction {
  if (typeof seriesFitFunc === 'string')
    return fitFunctions[seriesFitFunc];
  else if (!fitFunctions[seriesFitFunc.name]) {
    const name = seriesFitFunc.name;
    const paramNames = seriesFitFunc.parameterNames;
    const fitFunctionParts = seriesFitFunc.function.split('=>').map((elem) => elem.trim());
    const getInitParamsParts = seriesFitFunc.getInitialParameters.split('=>').map((elem) => elem.trim());
    const fitFunction = new Function(fitFunctionParts[0].slice(1, fitFunctionParts[0].length - 1),
      `${fitFunctionParts[1].includes(';') ? '' : 'return '}${fitFunctionParts[1]}`);
    const getInitParamsFunc = new Function(getInitParamsParts[0].slice(1, getInitParamsParts[0].length - 1),
      `return ${getInitParamsParts[1]}`);
    const fitFunc = new JsFunction(name, (fitFunction as (params: number[], x: number) => number),
      (getInitParamsFunc as (x: number[], y: number[]) => number[]), paramNames);
    fitFunctions[name] = fitFunc;
  }

  return fitFunctions[seriesFitFunc.name];
}

export function fitData(data: {x: number[], y: number[]}, fitFunction: FitFunction, errorModel: FitErrorModelType,
  parameterBounds?: FitParamBounds[]): FitCurve {
  const curveFunction = fitFunction.y;
  const paramValues = fitFunction.getInitialParameters(data.x, data.y);

  const of = createObjectiveFunction(errorModel);
  const fixed: number[] = [];
  let overLimits = true;

  while (overLimits) {
    const optimizable = createOptimizable(data, curveFunction, of, fixed);
    limitedMemoryBFGS(optimizable, paramValues);
    limitedMemoryBFGS(optimizable, paramValues);

    overLimits = false;
    if (!parameterBounds)
      break;

    for (let i = 0; i < parameterBounds.length; i++) {
      if (parameterBounds[i]?.max !== undefined && paramValues[i] > parameterBounds[i].max!) {
        overLimits = true;
        fixed.push(i);
        paramValues[i] = parameterBounds[i].max!;
        break;
      }
      if (parameterBounds[i]?.min !== undefined && paramValues[i] < parameterBounds[i].min!) {
        overLimits = true;
        fixed.push(i);
        paramValues[i] = parameterBounds[i].min!;
        break;
      }
    }
  }

  const fittedCurve = getFittedCurve(curveFunction, paramValues);

  return {
    fittedCurve: fittedCurve,
    parameters: paramValues,
  };
}

export function getFittedCurve(curveFunction: (params: number[], x: number) => number, paramValues: number[]):
 (x: number) => number {
  return (x: number) => {
    return curveFunction(paramValues, x);
  };
}

export function getCurveConfidenceIntervals(data: {x: number[], y: number[]}, paramValues: number[],
  curveFunction: (params: number[], x: number) => number, confidenceLevel: number = 0.05, errorModel: FitErrorModelType):
  FitConfidenceIntervals {
  const of = createObjectiveFunction(errorModel);

  const error = errorModel === FitErrorModel.PROPORTIONAL ?
    of(curveFunction, data, paramValues).mult :
    of(curveFunction, data, paramValues).const;

  const quantile = jStat.normal.inv(1 - confidenceLevel/2, 0, 1);

  const top = (x: number) => {
    const value = curveFunction(paramValues, x);
    if (errorModel === FitErrorModel.CONSTANT)
      return value + quantile * error;
    else
      return value + quantile * Math.abs(value) * error;
  };

  const bottom = (x: number) => {
    const value = curveFunction(paramValues, x);
    if (errorModel === FitErrorModel.CONSTANT)
      return value - quantile * error;
    else
      return value - quantile * Math.abs(value) * error;
  };

  return {confidenceTop: top, confidenceBottom: bottom};
}

export function getStatistics(data: {x: number[], y: number[]}, paramValues: number[],
  curveFunction: (params: number[], x: number) => number, statistics: boolean = true): FitStatistics {
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

export function sigmoid(params: number[], x: number): number {
  const A = params[0];
  const B = params[1];
  const C = params[2];
  const D = params[3];
  return (D + (A - D) / (1 + Math.pow(10, (x - C) * B)));
}

export function linear(params: number[], x: number) {
  const A = params[0];
  const B = params[1];
  return A * x + B;
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

function getObjectiveDerivative(of: ObjectiveFunction, curveFunction: (params: number[], x: number) => number,
  data: {x: number[], y: number[]}, params: number[], selectedParam: number): number {
  const step = (params[selectedParam] * 0.0001) === 0 ? 0.001 : (params[selectedParam] * 0.0001);
  const paramsTop: number[] = [];
  const paramsBottom: number[] = [];
  for (let i = 0; i < params.length; i++) {
    if (i === selectedParam) {
      paramsTop.push(params[i] + step);
      paramsBottom.push(params[i] - step);
    } else {
      paramsTop.push(params[i]);
      paramsBottom.push(params[i]);
    }
  }
  const drvTop = of(curveFunction, data, paramsTop).value;
  const drvBottom = of(curveFunction, data, paramsBottom).value;

  return (drvTop - drvBottom) / (2 * step);
}

function objectiveNormalConstant(targetFunc: (params: number[], x: number) => number,
  data: {y: number[], x: number[]}, params: number[]): Likelihood {
  //assure observed and args same length
  const pi = Math.PI;
  let sigma = 0;
  let sigmaSq = 0;
  let likelihood = 0;

  const residuesSquares = new Float32Array(data.x.length);
  for (let i = 0; i < data.x.length; i++) {
    const obs = data.y[i];
    const pred = targetFunc(params, data.x[i]);
    residuesSquares[i] = Math.pow(obs - pred, 2);
  }

  for (let i = 0; i < residuesSquares.length; i++)
    sigmaSq += residuesSquares[i];

  sigmaSq /= residuesSquares.length;
  sigma = Math.sqrt(sigmaSq);

  for (let i = 0; i < residuesSquares.length; i++)
    likelihood += residuesSquares[i] / sigmaSq + Math.log(2 * pi * sigmaSq);

  return {value: -likelihood, const: sigma, mult: 0};
}

function objectiveNormalProportional(targetFunc: (params: number[], x: number) => number,
  data: {y: number[], x: number[]}, params: number[]): Likelihood {
  //assure observed and args same length
  const pi = Math.PI;
  let sigma = 0;
  let sigmaSq = 0;
  let likelihood = 0;

  const residuesSquares = new Float32Array(data.x.length);
  for (let i = 0; i < data.x.length; i++) {
    const obs = data.y[i];
    const pred = targetFunc(params, data.x[i]);
    residuesSquares[i] = Math.pow(obs - pred, 2);
  }

  for (let i = 0; i < residuesSquares.length; i++)
    sigmaSq += residuesSquares[i];

  sigmaSq /= residuesSquares.length;
  sigma = Math.sqrt(sigmaSq);

  for (let i = 0; i < residuesSquares.length; i++)
    likelihood += residuesSquares[i] / sigmaSq + Math.log(2 * pi * sigmaSq);

  return {value: -likelihood, const: 0, mult: sigma};
}

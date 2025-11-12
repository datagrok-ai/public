/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import {
  FIT_FUNCTION_LINEAR,
  linear,
  IFitSeries,
  IFitPoint,
  FIT_FUNCTION_SIGMOID,
  sigmoid,
  FIT_FUNCTION_LOG_LINEAR,
  logLinear,
  FIT_FUNCTION_EXPONENTIAL,
  exponential,
  FIT_FUNCTION_4PL_REGRESSION,
  fourPLRegression,
  IFitFunctionDescription,
  FitParamBounds,
  FitMarkerType,
  FitOutlierMarkerType,
  FitLineStyle,
  FitErrorModelType, FitCurve, getAuc, getDetCoeff,
  FIT_JS_FUNCTION,
  FIT_FUNCTION_4PL_DOSE_RESPONSE,
  fourPLDoseResponse,
  DROPLINES,
  FitErrorModel,
  FitConfidenceIntervals,
  getFittedCurve,
} from './fit-curve';
import {fitSeries, getDataPoints} from './fit-data';
//@ts-ignore: no types
import * as jStat from 'jstat';
import {Extremum, OptimizationResult} from './fitting-algorithm/optimizer-misc';
import {performNelderMeadOptimization} from './fitting-algorithm/optimizer';
import {NELDER_MEAD_DEFAULTS} from './fitting-algorithm/optimizer-nelder-mead';


/** Class for the fit functions */
export abstract class FitFunction<T = Fit> {
  abstract get name(): string;
  abstract get parameterNames(): string[];
  abstract fillParams(fitCurve: FitCurve, data: FitSeries): T;
  abstract y(params: Float32Array, x: number): number;
  abstract getInitialParameters(x: number[], y: number[]): Float32Array;
}

export const FitFunctionTypes = {
  SIGMOID: 'sigmoid',
  LINEAR: 'linear',
  LOG_LINEAR: 'log-linear',
  EXPONENTIAL: 'exponential',
  FOUR_PL_REGRESSION: '4pl-regression',
  FOUR_PL_DOSE_RESPONSE: '4pl-dose-response',
} as const;

export type FitFunctionType = typeof FitFunctionTypes[keyof typeof FitFunctionTypes];

interface IFit {
  auc: number;
  rSquared: number;
}

interface ILinearFit extends IFit {
  slope: number;
  intercept: number;
}

interface ILogLinearFit extends ILinearFit {}

interface ISigmoidFit extends IFit {
  top: number;
  slope: number;
  ic50: number;
  bottom: number;
}

interface IExponentialFit extends IFit {
  mantissa: number;
  power: number;
}

interface IFourPLRegressionFit extends IFit {
  top: number;
  bottom: number;
  slope: number;
  ec50: number;
}

export abstract class Fit implements IFit {
  auc: number;
  rSquared: number;
  series: FitSeries;

  abstract get name(): string;

  protected constructor(values: IFit, data: FitSeries) {
    this.auc = values.auc;
    this.rSquared = values.rSquared;
    this.series = data;
  }
}

export class JSFunctionFit extends Fit {
  constructor(values: IFit, data: FitSeries) {
    super(values, data);
  }
  get name(): string {
    return FIT_JS_FUNCTION;
  }
}

class LinearFit extends Fit implements ILinearFit {
  slope: number;
  intercept: number;

  get name(): string {
    return FIT_FUNCTION_LINEAR;
  }

  constructor(values: ILinearFit, data: FitSeries) {
    super(values, data);
    this.slope = values.slope;
    this.intercept = values.intercept;
  }
}

class LogLinearFit extends LinearFit implements ILogLinearFit {
  get name(): string {
    return FIT_FUNCTION_LOG_LINEAR;
  }
}

class SigmoidFit extends Fit implements ISigmoidFit {
  top: number;
  slope: number;
  ic50: number;
  bottom: number;

  get name(): string {
    return FIT_FUNCTION_SIGMOID;
  }

  constructor(values: ISigmoidFit, data: FitSeries) {
    super(values, data);
    this.top = values.top;
    this.slope = values.slope;
    this.ic50 = values.ic50;
    this.bottom = values.bottom;
  }
}

class ExponentialFit extends Fit implements IExponentialFit {
  mantissa: number;
  power: number;

  get name(): string {
    return FIT_FUNCTION_EXPONENTIAL;
  }

  constructor(values: IExponentialFit, data: FitSeries) {
    super(values, data);
    this.mantissa = values.mantissa;
    this.power = values.power;
  }
}

class FourPLRegressionFit extends Fit implements IFourPLRegressionFit {
  top: number;
  bottom: number;
  slope: number;
  ec50: number;

  get name(): string {
    return FIT_FUNCTION_4PL_REGRESSION;
  }

  constructor(values: IFourPLRegressionFit, data: FitSeries) {
    super(values, data);
    this.top = values.top;
    this.bottom = values.bottom;
    this.slope = values.slope;
    this.ec50 = values.ec50;
  }
}

function getAucAndRsquared(fitCurve: (x: number) => number, data: {x: number[], y: number[]}): IFit {
  return {
    auc: getAuc(fitCurve, data),
    rSquared: getDetCoeff(fitCurve, data),
  };
}

/** Class that implements the linear function */
export class LinearFunction extends FitFunction<LinearFit> {
  get name(): string {
    return FIT_FUNCTION_LINEAR;
  }

  get parameterNames(): string[] {
    return ['Slope', 'Intercept'];
  }

  fillParams(fitCurve: FitCurve, data: FitSeries): LinearFit {
    return new LinearFit({
      ...getAucAndRsquared(fitCurve.fittedCurve, getDataPoints(data)),
      slope: fitCurve.parameters[0],
      intercept: fitCurve.parameters[1],
    }, data);
  }

  y(params: Float32Array, x: number): number {
    return linear(params, x);
  }

  getInitialParameters(x: number[], y: number[]): Float32Array {
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
    const params = new Float32Array(2);
    params.set([A, B]);
    return params;
  }
}

/** Class that implements the sigmoid function */
export class SigmoidFunction extends FitFunction<SigmoidFit> {
  get name(): string {
    return FIT_FUNCTION_SIGMOID;
  }

  get parameterNames(): string[] {
    return ['Top', 'Slope', 'IC50', 'Bottom'];
  }

  fillParams(fitCurve: FitCurve, data: FitSeries): SigmoidFit {
    return new SigmoidFit({
      ...getAucAndRsquared(fitCurve.fittedCurve, getDataPoints(data)),
      top: fitCurve.parameters[0],
      slope: fitCurve.parameters[1],
      ic50: fitCurve.parameters[2],
      bottom: fitCurve.parameters[3],
    }, data);
  }

  y(params: Float32Array, x: number): number {
    return sigmoid(params, x);
  }

  getInitialParameters(x: number[], y: number[]): Float32Array {
    const dataBounds = DG.Rect.fromXYArrays(x, y);
    const medY = (dataBounds.maxY - dataBounds.minY) / 2 + dataBounds.minY;
    let maxYInterval = dataBounds.maxY - dataBounds.minY;
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
    const params = new Float32Array(4);
    params.set([dataBounds.maxY, slope, xAtMedY, dataBounds.minY]);
    return params;
  }
}

/** Class that implements the linear logarithmic function */
export class LogLinearFunction extends FitFunction<LogLinearFit> {
  get name(): string {
    return FIT_FUNCTION_LOG_LINEAR;
  }

  get parameterNames(): string[] {
    return ['Slope', 'Intercept'];
  }

  fillParams(fitCurve: FitCurve, data: FitSeries): LogLinearFit {
    return new LogLinearFit({
      ...getAucAndRsquared(fitCurve.fittedCurve, getDataPoints(data)),
      slope: fitCurve.parameters[0],
      intercept: fitCurve.parameters[1],
    }, data);
  }

  y(params: Float32Array, x: number): number {
    return logLinear(params, x);
  }

  getInitialParameters(x: number[], y: number[]): Float32Array {
    const params = new Float32Array(2);
    params.set([-5, 100]);
    return params;
  }
}

/** Class that implements the exponential function */
export class ExponentialFunction extends FitFunction<ExponentialFit> {
  get name(): string {
    return FIT_FUNCTION_EXPONENTIAL;
  }

  get parameterNames(): string[] {
    return ['Mantissa', 'Power'];
  }

  fillParams(fitCurve: FitCurve, data: FitSeries): ExponentialFit {
    return new ExponentialFit({
      ...getAucAndRsquared(fitCurve.fittedCurve, getDataPoints(data)),
      mantissa: fitCurve.parameters[0],
      power: fitCurve.parameters[1],
    }, data);
  }

  y(params: Float32Array, x: number): number {
    return exponential(params, x);
  }

  getInitialParameters(x: number[], y: number[]): Float32Array {
    const params = new Float32Array(2);
    params.set([100, -2]);
    return params;
  }
}

/** Class that implements the Four Parameter Logistic Regression function */
export class FourPLRegressionFunction extends FitFunction<FourPLRegressionFit> {
  get name(): string {
    return FIT_FUNCTION_4PL_REGRESSION;
  }

  get parameterNames(): string[] {
    return ['Top', 'Slope', 'EC50', 'Bottom'];
  }

  fillParams(fitCurve: FitCurve, data: FitSeries): FourPLRegressionFit {
    return new FourPLRegressionFit({
      ...getAucAndRsquared(fitCurve.fittedCurve, getDataPoints(data)),
      top: fitCurve.parameters[0],
      slope: fitCurve.parameters[1],
      ec50: fitCurve.parameters[2],
      bottom: fitCurve.parameters[3],
    }, data);
  }

  y(params: Float32Array, x: number): number {
    return fourPLRegression(params, x);
  }

  getInitialParameters(x: number[], y: number[]): Float32Array {
    const params = new Float32Array(4);
    const bottom = Math.min(...y);
    const top = Math.max(...y);
    const medY = (top - bottom) / 2 + bottom;
    let maxYInterval = top -bottom;
    let nearestXIndex = 0;
    for (let i = 0; i < x.length; i++) {
      const currentInterval = Math.abs(y[i] - medY);
      if (currentInterval < maxYInterval) {
        maxYInterval = currentInterval;
        nearestXIndex = i;
      }
    }
    const ec50 = x[nearestXIndex];
    const slope = y[0] > y[y.length - 1]? -10 : 10;
    params.set([top, slope, ec50, bottom]);
    return params;
  }
}

export class FourPLDoseResponseFunction extends FourPLRegressionFunction {
  override get name(): string {
    return FIT_FUNCTION_4PL_DOSE_RESPONSE;
  }

  override get parameterNames(): string[] {
    return ['Max', 'Hill', 'IC50', 'Min'];
  }

  override y(params: Float32Array, x: number): number {
    return fourPLDoseResponse(params, x);
  }
}

/** Class that implements user JS functions */
export class JsFunction extends FitFunction<Fit> {
  private _name: string;
  private _parameterNames: string[];

  constructor(name: string, yFunc: (params: Float32Array, x: number) => number,
    getInitParamsFunc: (x: number[], y: number[]) => Float32Array, parameterNames: string[]) {
    super();

    this._name = name;
    this._parameterNames = parameterNames;

    this.y = yFunc;
    this.getInitialParameters = getInitParamsFunc;
  }

  fillParams(fitCurve: FitCurve, data: FitSeries): Fit {
    const params = new Float32Array(this._parameterNames.length);
    params.set(fitCurve.parameters);
    return new JSFunctionFit({...getAucAndRsquared(fitCurve.fittedCurve, getDataPoints(data))}, data);
  }

  get name(): string {
    return this._name;
  }

  get parameterNames(): string[] {
    return this._parameterNames;
  }

  y(params: Float32Array, x: number): number {
    throw new Error('Not implemented');
  }

  getInitialParameters(x: number[], y: number[]): Float32Array {
    throw new Error('Not implemented');
  }
}

// Object with fit functions
export const fitFunctions: {[key: string]: FitFunction<any>} = {
  'linear': new LinearFunction(),
  'sigmoid': new SigmoidFunction(),
  'log-linear': new LogLinearFunction(),
  'exponential': new ExponentialFunction(),
  '4pl-regression': new FourPLRegressionFunction(),
  '4pl-dose-response': new FourPLDoseResponseFunction(),
};

class FitFunctions {
  series: FitSeries;

  constructor(series: FitSeries) {
    this.series = series;
  }

  linear(): LinearFit {
    return fitFunctions[FIT_FUNCTION_LINEAR]
      .fillParams(fitSeries(this.series, fitFunctions[FIT_FUNCTION_LINEAR]), this.series);
  }

  logLinear(): LogLinearFit {
    return fitFunctions[FIT_FUNCTION_LOG_LINEAR]
      .fillParams(fitSeries(this.series, fitFunctions[FIT_FUNCTION_LOG_LINEAR]), this.series);
  }

  sigmoid(): SigmoidFit {
    return fitFunctions[FIT_FUNCTION_SIGMOID]
      .fillParams(fitSeries(this.series, fitFunctions[FIT_FUNCTION_SIGMOID]), this.series);
  }

  exponential(): ExponentialFit {
    return fitFunctions[FIT_FUNCTION_EXPONENTIAL]
      .fillParams(fitSeries(this.series, fitFunctions[FIT_FUNCTION_EXPONENTIAL]), this.series);
  }

  fourPL(): FourPLRegressionFit {
    return fitFunctions[FIT_FUNCTION_4PL_REGRESSION]
      .fillParams(fitSeries(this.series, fitFunctions[FIT_FUNCTION_4PL_REGRESSION]), this.series);
  }

  fourPLDoseResponse(): FourPLRegressionFit {
    return fitFunctions[FIT_FUNCTION_4PL_DOSE_RESPONSE].fillParams(fitSeries(this.series, fitFunctions[FIT_FUNCTION_4PL_DOSE_RESPONSE]), this.series);
  }
}


export class FitSeries implements IFitSeries {
  [key: string]: any;

  fit: FitFunctions;
  points: IFitPoint[];

  constructor(points: IFitPoint[]) {
    this.points = points;
    this.fit = new FitFunctions(this);
  }

  name?: string; // controls the series name
  fitFunction?: string | IFitFunctionDescription; // controls the series fit function
  parameters?: number[]; // controls the series parameters, auto-fitting when not defined
  parameterBounds?: FitParamBounds[]; // defines the acceptable range of each parameter, which is taken into account during the fitting. See also `parameters`.
  markerType?: FitMarkerType; // defines the series marker type
  outlierMarkerType?: FitOutlierMarkerType; // defines the series outlier marker type
  lineStyle?: FitLineStyle; // defines the series line style
  pointColor?: string; // overrides the standardized series point color
  fitLineColor?: string; // overrides the standardized series fit line color
  confidenceIntervalColor?: string; // overrides the standardized series confidence interval color
  outlierColor?: string; // overrides the standardized series outlier color
  connectDots?: boolean; // defines whether to connect the points with lines or not. If true and showFitLine is false - fitting is disabled - otherwise, it will be rendered accordingly to the parameter value.
  showFitLine?: boolean; // defines whether to show the fit line or not
  showPoints?: string; // defines the data display mode
  showOutliers?: boolean; // defines whether to show the outliers or not
  showCurveConfidenceInterval?: boolean; // defines whether to show the confidence intervals or not
  errorModel?: FitErrorModelType; // defines the series error model
  clickToToggle?: boolean; // if true, clicking on the point toggles its outlier status and causes curve refitting
  labels?: {[key: string]: string | number | boolean}; // controlled by IFitChartData labelOptions, shows labels
  droplines?: string[]; // defines the droplines that would be shown on the plot (IC50)
  columnName?: string; // defines the column name where the series is stored
}


/** Properties that describe {@link IFitSeriesOptions}. Useful for editing, initialization, transformations, etc. */
export const fitSeriesProperties: DG.Property[] = [
  DG.Property.js('fitFunction', DG.TYPE.STRING,
    {category: 'Fitting', choices: Object.keys(fitFunctions), defaultValue: 'sigmoid'}),
  DG.Property.js('pointColor', DG.TYPE.STRING,
    {category: 'Rendering', nullable: true, inputType: 'Color'}),
  DG.Property.js('fitLineColor', DG.TYPE.STRING,
    {category: 'Rendering', nullable: true, inputType: 'Color'}),
  DG.Property.js('outlierColor', DG.TYPE.STRING,
    {category: 'Rendering', nullable: true, inputType: 'Color'}),
  DG.Property.js('errorModel', DG.TYPE.STRING, {category: 'Fitting', defaultValue: 'constant',
    choices: ['constant', 'proportional', 'combined'], nullable: false}),
  DG.Property.js('connectDots', DG.TYPE.BOOL, {category: 'Fitting', defaultValue: false}),
  DG.Property.js('clickToToggle', DG.TYPE.BOOL, {category: 'Fitting', description:
    'Click on a point to mark it as outlier and refit', nullable: true, defaultValue: false}),
  DG.Property.js('showFitLine', DG.TYPE.BOOL, {category: 'Fitting', defaultValue: true}),
  DG.Property.js('showPoints', DG.TYPE.STRING, // rewrite description
    {category: 'Fitting', description: 'Whether points/candlesticks/none should be rendered',
      defaultValue: 'points', choices: ['points', 'candlesticks', 'both']}),
  DG.Property.js('showOutliers', DG.TYPE.BOOL, {category: 'Fitting', defaultValue: true}),
  DG.Property.js('showCurveConfidenceInterval', DG.TYPE.BOOL,
    {category: 'Fitting', description: 'Whether confidence intervals should be rendered', defaultValue: false,
      //@ts-ignore
      friendlyName: 'Confidence Interval'}),
  DG.Property.js('markerType', DG.TYPE.STRING, {category: 'Rendering', defaultValue: 'circle',
    choices: ['asterisk', 'circle', 'cross border', 'diamond', 'square', 'star',
      'triangle bottom', 'triangle left', 'triangle right', 'triangle top'], nullable: false,
    //@ts-ignore
    friendlyName: 'Marker'}),
  DG.Property.js('outlierMarkerType', DG.TYPE.STRING, {category: 'Rendering', defaultValue: 'outlier',
    choices: ['asterisk', 'circle', 'cross border', 'diamond', 'outlier', 'square', 'star',
      'triangle bottom', 'triangle left', 'triangle right', 'triangle top'], nullable: false,
    //@ts-ignore
    friendlyName: 'Outlier Marker'}),
  DG.Property.js('lineStyle', DG.TYPE.STRING, {category: 'Rendering', defaultValue: 'solid',
    choices: ['solid', 'dotted', 'dashed', 'dashdotted'], nullable: false}),
  DG.Property.js('droplines', DG.TYPE.STRING_LIST, {description: 'Whether specific droplines should be rendered',
    choices: DROPLINES, inputType: 'MultiChoice'}),
  DG.Property.js('columnName', DG.TYPE.STRING, {description: 'Column name where the series is stored', defaultValue: ''}),
];


export function getOrCreateFitFunction(seriesFitFunc: string | IFitFunctionDescription): FitFunction<Fit> {
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
    const fitFunc = new JsFunction(name, (fitFunction as (params: Float32Array, x: number) => number),
      (getInitParamsFunc as (x: number[], y: number[]) => Float32Array), paramNames);
    fitFunctions[name] = fitFunc;
  }

  return fitFunctions[seriesFitFunc.name];
}

export function fitData(data: {x: number[], y: number[]}, fitFunction: FitFunction<Fit>, errorModel?: FitErrorModelType,
  parameterBounds?: FitParamBounds[]): FitCurve {
  errorModel ??= FitErrorModel.CONSTANT as FitErrorModelType;
  const curveFunction = fitFunction.y;
  let paramValues = fitFunction.getInitialParameters(data.x, data.y);

  const of = objectiveFactory(curveFunction, data, errorModel);
  const bottomParamBounds = new Float32Array(fitFunction.parameterNames.length);
  const topParamBounds = new Float32Array(fitFunction.parameterNames.length);
  for (let i = 0; i < fitFunction.parameterNames.length; i++) {
    bottomParamBounds[i] = paramValues[i] === 0 ? -1 : paramValues[i] - Math.abs(paramValues[i] * 0.5);
    topParamBounds[i] = paramValues[i] === 0 ? 1 : paramValues[i] + Math.abs(paramValues[i] * 0.5);
  }
  const parameterBoundsBitset: DG.BitSet = DG.BitSet.create(fitFunction.parameterNames.length * 2);
  if (parameterBounds && parameterBounds.length !== 0) {
    for (let i = 0; i < parameterBounds.length; i++) {
      if (parameterBounds[i].min !== undefined && parameterBounds[i].min !== null) {
        bottomParamBounds[i] = parameterBounds[i].min!;
        parameterBoundsBitset.set(i * 2, true);
      }
      if (parameterBounds[i].max !== undefined && parameterBounds[i].max !== null) {
        topParamBounds[i] = parameterBounds[i].max!;
        parameterBoundsBitset.set(i * 2 + 1, true);
      }
    }
  }
  const getConsistency = (extremum: Extremum) => {
    const residuals = of(extremum.point).residuals;
    let q1q4 = 0;
    let q2q3 = 0;
    // TODO: rewrite quartiles to a better condition
    for (let i = 0; i < residuals.length; i++) {
      if (residuals.length / 4 > i)
        q1q4 += Math.pow(residuals[i], 2);
      else if (3 * residuals.length / 4 > i)
        q2q3 += Math.pow(residuals[i], 2);
      else
        q1q4 += Math.pow(residuals[i], 2);
    }
    return q1q4 / q2q3;
  };

  let iter = -1;
  let continueOptimization = 0;
  let statistics: number | null = null;
  let optimization: OptimizationResult;
  let minIdx = 0;
  do {
    optimization = performNelderMeadOptimization(of, bottomParamBounds, topParamBounds, {
      tolerance: NELDER_MEAD_DEFAULTS.TOLERANCE,
      maxIter: NELDER_MEAD_DEFAULTS.MAX_ITER,
      nonZeroParam: NELDER_MEAD_DEFAULTS.NON_ZERO_PARAM,
      initialScale: NELDER_MEAD_DEFAULTS.INITIAL_SCALE,
      scaleReflection: NELDER_MEAD_DEFAULTS.SCALE_REFLECTION,
      scaleExpansion: NELDER_MEAD_DEFAULTS.SCALE_EXPANSION,
      scaleContraction: NELDER_MEAD_DEFAULTS.SCALE_CONTRACTION,
    }, 10, paramValues);

    minIdx = 0;
    for (let i = 1; i < optimization.extremums.length; i++) {
      if (optimization.extremums[i].cost < optimization.extremums[minIdx].cost)
        minIdx = i;
    }

    const newStatistics = getConsistency(optimization.extremums[minIdx]);
    if (statistics === null)
      statistics = newStatistics;
    else if (newStatistics < statistics) {
      statistics = newStatistics;
      continueOptimization = 0;
    }
    else
      continueOptimization++;

    paramValues = optimization.extremums[minIdx].point;

    if (iter > 40 || continueOptimization === 3)
      break;

    for (let i = 0; i < fitFunction.parameterNames.length; i++) {
      if (!parameterBoundsBitset.get(i * 2))
        bottomParamBounds[i] = paramValues[i] - Math.abs(paramValues[i] * (continueOptimization + 0.5));
      if (!parameterBoundsBitset.get(i * 2 + 1))
        topParamBounds[i] = paramValues[i] + Math.abs(paramValues[i] * (continueOptimization + 0.5));
    }

    iter++;
  } while (true);

  const fittedCurve = getFittedCurve(curveFunction, paramValues);

  return {
    fittedCurve: fittedCurve,
    parameters: paramValues,
  };
}


function objectiveFactory(targetFunc: (params: Float32Array, x: number) => number,
  data: {y: number[], x: number[]}, errorModel: FitErrorModelType):
    (params: Float32Array) => {likelihood: number, residuals: number[]} {
  return (params: Float32Array) => {
    let likelihood = 0;
    let sigmaA = 0;
    let sigmaB = 0;

    const obs = data.y;
    const pred = data.x.map((x) => targetFunc(params, x));
    const residuals = obs.map((obs, i) => obs - pred[i]);
    const residuesSquares = residuals.map((res) => Math.pow(res, 2));

    for (let i = 0; i < residuesSquares.length; i++) {
      if ([FitErrorModel.CONSTANT, FitErrorModel.COMBINED].includes(errorModel))
        sigmaA += residuesSquares[i];
      else if (errorModel === FitErrorModel.PROPORTIONAL)
        sigmaB += residuesSquares[i] / (pred[i] * pred[i]);
    }
    sigmaA = Math.sqrt(sigmaA / residuesSquares.length);
    sigmaB = Math.sqrt(sigmaB / residuesSquares.length);

    if (errorModel === FitErrorModel.COMBINED)
      ({sigmaA, sigmaB} = getSigmasForCombinedErrorModel(sigmaA, sigmaB, residuesSquares, pred));

    const sigmaSqI: number[] = [];
    for (let i = 0; i < residuesSquares.length; i++) {
      sigmaSqI[i] = Math.pow(sigmaA + sigmaB * pred[i], 2);
      likelihood += residuesSquares[i] / sigmaSqI[i] + Math.log(2 * Math.PI * sigmaSqI[i]);
    }

    return {likelihood, residuals};
  };
}


function getSigma(targetFunc: (params: Float32Array, x: number) => number,
  data: {y: number[], x: number[]}, errorModel: FitErrorModelType, params: Float32Array): {sigmaA: number, sigmaB: number} {
  let sigmaA = 0;
  let sigmaB = 0;

  const obs = data.y;
  const pred = data.x.map((x) => targetFunc(params, x));
  const residuals = obs.map((obs, i) => obs - pred[i]);
  const residuesSquares = residuals.map((res) => Math.pow(res, 2));

  for (let i = 0; i < residuesSquares.length; i++) {
    if ([FitErrorModel.CONSTANT, FitErrorModel.COMBINED].includes(errorModel))
      sigmaA += residuesSquares[i];
    else if (errorModel === FitErrorModel.PROPORTIONAL)
      sigmaB += Math.pow(residuals[i] / Math.abs(pred[i]), 2);
  }
  sigmaA = Math.sqrt(sigmaA / residuesSquares.length);
  sigmaB = Math.sqrt(sigmaB / residuesSquares.length);

  if (errorModel === FitErrorModel.COMBINED)
    ({sigmaA, sigmaB} = getSigmasForCombinedErrorModel(sigmaA, sigmaB, residuesSquares, pred));

  return {sigmaA, sigmaB};
}

function getSigmasForCombinedErrorModel(sigmaA: number, sigmaB: number, residuesSquares: number[], pred: number[]):
  {sigmaA: number, sigmaB: number} {
  let likelihoodCombined = Number.MAX_VALUE;
  let finalProportion = 0;
  for (let i = 0; i <= 10; i++) {
    let likelihoodCombinedLocal = 0;
    const proportion = i / 10;
    const sigmaACombined = sigmaA * proportion;
    const sigmaBCombined = sigmaA * (1 - proportion);
    const sigmaSqICombined: number[] = [];
    for (let j = 0; j < residuesSquares.length; j++) {
      sigmaSqICombined[j] = Math.pow(sigmaACombined + sigmaBCombined * pred[j], 2);
      likelihoodCombinedLocal += residuesSquares[j] / sigmaSqICombined[j] + Math.log(2 * Math.PI * sigmaSqICombined[j]);
    }
    if (likelihoodCombinedLocal < likelihoodCombined) {
      likelihoodCombined = likelihoodCombinedLocal;
      finalProportion = proportion;
    }
  }
  sigmaB = sigmaA * (1 - finalProportion);
  sigmaA = sigmaA * finalProportion;

  return {sigmaA, sigmaB};
}


export function getCurveConfidenceIntervals(data: {x: number[], y: number[]}, paramValues: Float32Array,
  curveFunction: (params: Float32Array, x: number) => number, confidenceLevel: number = 0.05, errorModel: FitErrorModelType):
  FitConfidenceIntervals {
  const error = getSigma(curveFunction, data, errorModel, paramValues);
  const quantile = jStat.normal.inv(1 - confidenceLevel/2, 0, 1);

  const top = (x: number) => {
    const value = curveFunction(paramValues, x);
    return value + quantile * (error.sigmaA + Math.abs(value) * error.sigmaB);
  };

  const bottom = (x: number) => {
    const value = curveFunction(paramValues, x);
    return value - quantile * (error.sigmaA + Math.abs(value) * error.sigmaB);
  };

  return {confidenceTop: top, confidenceBottom: bottom};
}


// const series: FitSeries = new FitSeries([
//   {'x': 0, 'y': 0},
//   {'x': 1, 'y': 0.5},
//   {'x': 2, 'y': 1},
//   {'x': 3, 'y': 10, 'outlier': true},
//   {'x': 4, 'y': 0},
// ]);

export interface FitCellOutlierToggleArgs {
  gridCell: DG.GridCell;
  series: IFitSeries;
  seriesIdx: number;
  pointIdx: number;
};

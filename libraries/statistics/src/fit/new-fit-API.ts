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
} from './fit-curve';
import {fitSeries, getDataPoints} from './fit-data';


/** Class for the fit functions */
export abstract class FitFunction<T> {
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


// Object with fit functions
export const fitFunctions = {
  'linear': new LinearFunction(),
  'sigmoid': new SigmoidFunction(),
  'log-linear': new LogLinearFunction(),
  'exponential': new ExponentialFunction(),
  '4pl-regression': new FourPLRegressionFunction(),
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
}


export class FitSeries implements IFitSeries {
  [key: string]: any;

  fit: FitFunctions;
  points: IFitPoint[];

  constructor(points: IFitPoint[]) {
    this.points = points;
    this.fit = new FitFunctions(this);
  }

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
  showCurveConfidenceInterval?: boolean;    // defines whether to show the confidence intervals or not
  errorModel?: FitErrorModelType;       // defines the series error model
  clickToToggle?: boolean;    // if true, clicking on the point toggles its outlier status and causes curve refitting
  labels?: {[key: string]: string | number | boolean}; // controlled by IFitChartData labelOptions, shows labels
  droplines?: string[];                 // defines the droplines that would be shown on the plot (IC50)
  columnName?: string;                  // defines the column name where the series is stored
}


const series: FitSeries = new FitSeries([
  {'x': 0, 'y': 0},
  {'x': 1, 'y': 0.5},
  {'x': 2, 'y': 1},
  {'x': 3, 'y': 10, 'outlier': true},
  {'x': 4, 'y': 0},
]);

import {limitedMemoryBFGS} from "../../lbfgs/lbfgs"
//@ts-ignore: no types
import * as jStat from 'jstat';
import {Property} from "datagrok-api/src/entities";
import {TYPE} from "datagrok-api/src/const";


type Likelihood = {
  value: number,
  const: number,
  mult: number
};


export type FitResult = {
  parameters: number[],
  fittedCurve: (x:number)=> number,
  confidenceTop: (x:number)=> number,
  confidenceBottom: (x:number)=> number,

  rSquared?: number,
  auc?: number,
  inverted?: (y:number)=> number,
  invertedTop?: (y:number)=> number,
  invertedBottom?: (y:number)=> number,
  interceptX: number, // parameters[2]
  interceptY: number, // fittedCurve[parameters[2]]
  slope: number, // parameters[1]
  top: number, // parameters[0]
  bottom: number, // parameters[3]
};


/** Properties that describe {@link FitResult}. Useful for editing, initialization, transformations, etc. */
export const fitResultProperties: Property[] = [
  Property.js('rSquared', TYPE.FLOAT, {userEditable: false}),
  Property.js('auc', TYPE.FLOAT, {userEditable: false}),
  Property.js('interceptY', TYPE.FLOAT, {userEditable: false}),
  Property.js('interceptX', TYPE.FLOAT, {userEditable: false}),
  Property.js('slope', TYPE.FLOAT, {userEditable: false}),
  Property.js('top', TYPE.FLOAT, {userEditable: false}),
  Property.js('bottom', TYPE.FLOAT, {userEditable: false}),
];


type ObjectiveFunction = (targetFunc: (params: number[], x: number) => number,
  data: {x: number[], y: number[]},
  params: number[]) => Likelihood;


export enum FitErrorModel {
  Constant,
  Proportional
}


export const FIT_FUNCTION_SIGMOID = 'Sigmoid';
export const FIT_FUNCTION_LINEAR = 'Linear';

// export const FIT_STATS_RSQUARED = 'rSquared';
// export const FIT_STATS_AUC = 'auc';


export type FitFunctionType = 'Sigmoid' | 'Linear';

abstract class FitFunction {
  abstract get name(): string;
  abstract get parameterNames(): string[];
  abstract y(params: number[], x: number): number;
  abstract getInitialParameters(x: number[], y: number[]): number[];
}


class LinearFunction extends FitFunction {
  get name(): string { return FIT_FUNCTION_LINEAR; }

  get parameterNames(): string[] {
    return ['Slope', 'Intercept'];
  }

  y(params: number[], x: number): number {
    throw 'Not implemented';
  }

  getInitialParameters(x: number[], y: number[]): number[] {
    throw 'Not implemented';
  }
}


class SigmoidFunction extends FitFunction {
  get name(): string { return FIT_FUNCTION_SIGMOID; }

  get parameterNames(): string[] {
    return ['Top', 'Bottom', 'Slope', 'IC50'];
  }

  y(params: number[], x: number): number {
    throw 'Not implemented';
  }

  getInitialParameters(x: number[], y: number[]): number[] {
    throw 'Not implemented';
  }
}


export const fitFunctions: {[index: string]: any} = {
  FIT_FUNCTION_LINEAR: new LinearFunction(),
  FIT_FUNCTION_SIGMOID: new SigmoidFunction(),
};


export interface IFitOptions {
  errorModel: FitErrorModel;
  confidenceLevel: number;
  statistics: boolean;
}


/**
 * statistics - whether or not to calculate fit statistics (potentially computationally intensive)
 * */
export function fit(data:{x: number[], y: number[]},
                    params: {paramvalue: number, min?:number, max?: number}[],
                    curveFunction: (paramValues: number[], x: number) => number,
                    errorModel: FitErrorModel,
                    confidenceLevel: number = 0.05,
                    statistics: boolean = true): FitResult {

  let paramValues = params.map(p => p.paramvalue);
  let of: ObjectiveFunction;
  switch(errorModel) {
    case FitErrorModel.Constant:
      of = objectiveNormalConstant;
      break;
    case FitErrorModel.Proportional:
      of = objectiveNormalProportional;
      break;
    default:
      of = objectiveNormalConstant;
      break;
  }

  let iterations = 0;

  let optimizable = {
    getValue: (parameters: number[]) => {
      return of(curveFunction, data, parameters).value;
    },
    getGradient: (parameters: number[], gradient: number[]) => {
      const length = Object.keys(parameters).length;
      iterations++;

      for (let i = 0; i < parameters.length; i++)
        gradient[i] = getObjectiveDerivative(of, curveFunction, data, parameters, i);

      return gradient;
    }
  };

  limitedMemoryBFGS(optimizable, paramValues);
  limitedMemoryBFGS(optimizable, paramValues);

  let fittedCurve = (x: number) => {
    return curveFunction(paramValues, x);
  }

  let error = errorModel == FitErrorModel.Proportional ?
  of(curveFunction, data, paramValues).mult :
  of(curveFunction, data, paramValues).const;

  let studentQ = jStat.studentt.inv(1 - confidenceLevel/2, data.x.length - paramValues.length);

  let top = (x: number) =>{
    let value = curveFunction(paramValues, x);
    if (errorModel == FitErrorModel.Constant)
      return  value + studentQ*error/Math.sqrt(data.x.length);
    else
      return  value + studentQ*(Math.abs(value)*error/Math.sqrt(data.x.length));
  }

  let bottom = (x: number) => {
    let value = curveFunction(paramValues, x);
    if (errorModel == FitErrorModel.Constant)
      return  value - studentQ*error/Math.sqrt(data.x.length);
    else
      return  value - studentQ*(Math.abs(value)*error/Math.sqrt(data.x.length));
  }

  let inv: (y: number) => number = (y: number) => {return 0;};
  let invTop: (y: number) => number = (y: number) => {return 0;};
  let invBottom: (y: number) => number = (y: number) => {return 0;};

  if (statistics) {
    inv = (y: number) => {
      //should check if more than bottom and less than top
      return paramValues[2]/Math.pow((paramValues[0] - y)/(y - paramValues[3]), 1/paramValues[1]);
    };

    let error = getInvError(inv, data);

    invTop = (y: number) =>{
      let value = inv(y);
      return  value + studentQ*error/Math.sqrt(data.y.length);
    }

    invBottom = (y: number) => {
      let value = inv(y);
      return  value - studentQ*error/Math.sqrt(data.y.length);
    }
  }

  let fitRes: FitResult = {
    parameters: paramValues,
    fittedCurve: fittedCurve,
    confidenceTop: top,
    confidenceBottom: bottom,
    rSquared: statistics ? getDetCoeff(fittedCurve, data) : undefined,
    auc: statistics ? getAuc(fittedCurve, data) : undefined,
    inverted: statistics ? inv : undefined,
    invertedTop: statistics ? invTop : undefined,
    invertedBottom: statistics ? invBottom : undefined,
    interceptX : paramValues[2],
    interceptY : fittedCurve(paramValues[2]),
    slope: paramValues[1],
    top: paramValues[0],
    bottom: paramValues[3]
  };

  return fitRes;
}

export function sigmoid(params: number[], x: number): number {
  const A = params[0];
  const B = params[1];
  const C = params[2];
  const D = params[3];
  const res = D + (A - D)/(1 + Math.pow(10, (x - C)*B));
  return res;
}

function getObjectiveDerivative(of: ObjectiveFunction, curveFunction: (params: number[], x: number) => number,
    data: {x: number[], y: number[]}, params: number[], selectedParam: number): number {
  let step = params[selectedParam]*0.0001;
  step = step == 0 ? 0.001 : step;
  let paramsTop: number[] = [];
  let paramsBottom: number[] = [];
  for (let i = 0; i < params.length; i++) {
    if(i == selectedParam) {
      paramsTop.push(params[i] + step);
      paramsBottom.push(params[i] - step);
    } else {
      paramsTop.push(params[i]);
      paramsBottom.push(params[i]);
    }
  }
  const drvTop = of(curveFunction, data, paramsTop).value;
  const drvBottom = of(curveFunction, data, paramsBottom).value;

  return (drvTop - drvBottom)/(2*step);
}

function getAuc(fittedCurve: (x: number) => number,
                data: {x: number[], y: number[]}): number {
  let auc = 0;

  let min = Math.min(...data.x);
  let max = Math.max(...data.x);
  const integrationStep = (max - min)/1000;

  for(let x = min; x < max; x+= integrationStep)
    auc += integrationStep*fittedCurve(x);

  return auc;
}

function getDetCoeff(fittedCurve: (x: number) => number,
                     data: {x: number[], y: number[]}): number {
  let ssRes = 0;
  let ssTot = 0;

  const yMean = jStat.mean(data.y);

  for(let i = 0; i < data.x.length; i++) {
    ssRes += Math.pow(data.y[i] - fittedCurve(data.x[i]), 2);
    ssTot += Math.pow(data.y[i] - yMean, 2);
  }

  return 1 - ssRes/ssTot;
}

function getInvError (
  targetFunc: (y: number) => number,
  data: {y: number[], x: number[]}
): number {
  const pi = Math.PI;
  let sigma = 0;
  let sigmaSq = 0;

  let residuesSquares = new Float32Array(data.y.length);
  for(let i = 0; i < data.y.length; i++) {
    const obs = data.x[i];
    const pred = targetFunc(data.y[i]);
    residuesSquares[i] = Math.pow(obs - pred, 2);
  }

  for(let i = 0; i < residuesSquares.length; i++)
    sigmaSq += residuesSquares[i];
  sigmaSq /= residuesSquares.length;
  sigma = Math.sqrt(sigmaSq);

  return sigma;
}

function objectiveNormalConstant (
  targetFunc: (params: number[], x: number) => number,
  data: {y: number[], x: number[]},
  params: number[]
): Likelihood {
  //assure observed and args same length
  const pi = Math.PI;
  let sigma = 0;
  let sigmaSq = 0;
  let likelihood = 0;

  let residuesSquares = new Float32Array(data.x.length);
  for(let i = 0; i < data.x.length; i++) {
    const obs = data.y[i];
    const pred = targetFunc(params, data.x[i]);
    residuesSquares[i] = Math.pow(obs - pred, 2);
  }

  for(let i = 0; i < residuesSquares.length; i++)
    sigmaSq += residuesSquares[i];
  sigmaSq /= residuesSquares.length;
  sigma = Math.sqrt(sigmaSq);

  for(let i = 0; i < residuesSquares.length; i++)
    likelihood += residuesSquares[i]/sigmaSq + Math.log(2 * pi * sigmaSq);

  return {value: -likelihood, const: sigma, mult: 0};
}

function objectiveNormalProportional (
targetFunc: (params: number[], x: number) => number,
data: {y: number[], x: number[]},
params: number[]
): Likelihood {
  //assure observed and args same length
  const pi = Math.PI;
  let sigma = 0;
  let sigmaSq = 0;
  let likelihood = 0;

  let residuesSquares = new Float32Array(data.x.length);
  for(let i = 0; i < data.x.length; i++) {
    const obs = data.y[i];
    const pred = targetFunc(params, data.x[i])
    residuesSquares[i] = Math.pow(obs - pred, 2);
  }

  for(let i = 0; i < residuesSquares.length; i++)
    sigmaSq += residuesSquares[i];
  sigmaSq /= residuesSquares.length;
  sigma = Math.sqrt(sigmaSq);

  for(let i = 0; i < residuesSquares.length; i++)
    likelihood += residuesSquares[i]/sigmaSq + Math.log(2*pi*sigmaSq);

  return {value: -likelihood, const: sigma, mult: 0};
}

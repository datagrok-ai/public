import {limitedMemoryBFGS} from "../../lbfgs/lbfgs"
//@ts-ignore: no types
import * as jStat from 'jstat';

type Likelihood = {
  value: number, 
  const: number, 
  mult: number
};

type FitResult = {
  parameters: number[],
  fittedCurve: (x:number)=> number,
  confidenceTop: (x:number)=> number,
  confidenceBottom: (x:number)=> number,
  rSquared: number,
  auc: number;
};

type ObjectiveFunction = (targetFunc: (params: number[], x: number) => number, 
data: {x: number[], y: number[]},
params: number[]) => Likelihood;

export enum FitErrorModel {
  Constant,
  Proportional
}

export function fit(data:{x: number[], y: number[]}, 
                    params: number[],
                    curveFunction: (params: number[], x: number) => number, 
                    errorModel: FitErrorModel,
                    confidenceLevel: number = 0.05): FitResult {

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

  limitedMemoryBFGS(optimizable, params);
  limitedMemoryBFGS(optimizable, params);

  let fittedCurve = (x: number) => {
    return curveFunction(params, x);
  }

  let error = errorModel == FitErrorModel.Proportional ?
  of(curveFunction, data, params).mult :
  of(curveFunction, data, params).const;

  let studentQ = jStat.studentt.inv(1 - confidenceLevel/2, data.x.length - params.length);

  let top = (x: number) =>{
    let value = curveFunction(params, x);
    if (errorModel == FitErrorModel.Constant)
      return  value + studentQ*error/Math.sqrt(data.x.length);
    else
      return  value + studentQ*(Math.abs(value)*error/Math.sqrt(data.x.length));
  }

  let bottom = (x: number) => {
    let value = curveFunction(params, x);
    if (errorModel == FitErrorModel.Constant)
      return  value - studentQ*error/Math.sqrt(data.x.length);
    else
      return  value - studentQ*(Math.abs(value)*error/Math.sqrt(data.x.length));
  }

  let fitRes: FitResult = {
    parameters: params,
    fittedCurve: fittedCurve,
    confidenceTop: top,
    confidenceBottom: bottom,
    rSquared: getDetCoeff(fittedCurve, data),
    auc: getAuc(fittedCurve, data)
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
  const integrationStep = 0.00001;
  let min = Math.min(...data.x);
  let max = Math.max(...data.x);


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

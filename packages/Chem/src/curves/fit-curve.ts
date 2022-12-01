import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

type Likelihood = {
  value: number, 
  const: number, 
  mult: number
};

interface ICurveFitter<TArg> {
  data: {observed: number[], args: TArg[]};
  params: number[];

  //curve function
  cf (params: number[], args: TArg): number
  //objective function
  of (targetFunc: (params: number[], args: TArg) => number, 
                                  data: {observed: number[], args: TArg[]},
                                  params: number[]): Likelihood
}
export class CurveFitter<TArg> implements ICurveFitter<TArg> {
  data: {observed: number[], args: TArg[]};
  params: number[];
  cf;
  of;
  optimizable;
  fitHistory: {iterations: number, parameters: {[iteration: number]: number[]}, likelihood: {[iteration: number]: number}} = {
    iterations: 0,
    parameters: {},
    likelihood: {}
  };

  constructor(curveFunction: (params: number[], args: TArg) => number, 
              data: {observed: number[], args: TArg[]},
              params: number[],
              errorModel: 'constant' | 'proportional' = 'constant'//| 'combinational1' | 'combinational2'
              ) {

    if(data.observed.length == 0 || data.observed.length != data.args.length)
      throw new Error('invalid data');
      
    this.data = data;
    this.params = params;

    this.cf = curveFunction;

    switch(errorModel) {
      case 'constant':
        this.of = objectiveNormalConstant;
      case 'proportional':
        this.of = objectiveNormalProportional;
      // case 'combinational1':
      //   this.of = objectiveNormalAdditive;
      // case 'combinational2':
      //   this.of = objectiveNormalAdditive;
      default:
        this.of = objectiveNormalConstant;
    }

    this.optimizable = {
      getValue: (parameters: number[]) => {
        return this.of(this.cf, this.data, parameters).value;
      },
      getGradient: (parameters: number[], gradient: number[]) => {
        const length = Object.keys(this.fitHistory.parameters).length;
        this.fitHistory.iterations++;
        this.fitHistory.parameters[length] = parameters.slice();
        this.fitHistory.likelihood[length] = this.of(this.cf, this.data, parameters).value;
        console.log(parameters.slice());
        
        for (let i = 0; i < parameters.length; i++)
          gradient[i] = this.getDerivative(parameters, i);

        return gradient;
      }
    };
  }

  getFittedCurve() {
    limitedMemoryBFGS(this.optimizable, this.params);
    limitedMemoryBFGS(this.optimizable, this.params);
    let res = (arg: TArg) =>{
      let params = this.fitHistory.parameters[this.fitHistory.iterations];

      return this.cf(params, arg);
    }

    return res;
  }

  getFittedParams() {
    return this.params;
  }

  //central difference
  getDerivative(params: number[], selectedParam: number): number {
    let step = params[selectedParam]*0.0001;
    step = step == 0 ? 0.001 : step; 
    let paramsTop = [];
    let paramsBottom = [];
    for (let i = 0; i < params.length; i++) {
      if(i == selectedParam) {
        paramsTop.push(params[i] + step);
        paramsBottom.push(params[i] - step);
      } else {
        paramsTop.push(params[i]);
        paramsBottom.push(params[i]);
      }
    }
    const drvTop = this.of(this.cf, this.data, paramsTop).value;
    const drvBottom = this.of(this.cf, this.data, paramsBottom).value;

    return (drvTop - drvBottom)/(2*step);
  }
}

export function sigmoid(params: number[], x: number): number {
  const A = params[0];
  const B = params[1];
  const C = params[2];
  const D = params[3];
  const res = D + (A - D)/(1 + Math.pow(10, (x - C)*B));
  return res;
}

function objectiveNormalConstant<TArg> (
                                        targetFunc: (params: number[], args: TArg) => number, 
                                        data: {observed: number[], args: TArg[]},
                                        params: number[]
): Likelihood {
  //assure observed and args same length
  const pi = Math.PI;
  let sigma = 0;
  let sigmaSq = 0;
  let likelihood = 0;

  let residuesSquares = new Float32Array(data.args.length);
  for(let i = 0; i < data.args.length; i++) {
    const obs = data.observed[i];
    const pred = targetFunc(params, data.args[i]);
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

function objectiveNormalProportional<TArg> (
  targetFunc: (params: number[], args: TArg) => number, 
  data: {observed: number[], args: TArg[]},
  params: number[]
): Likelihood {
//assure observed and args same length
const pi = Math.PI;
let sigma = 0;
let sigmaSq = 0;
let likelihood = 0;

let residuesSquares = new Float32Array(data.args.length);
for(let i = 0; i < data.args.length; i++) {
const obs = data.observed[i];
const pred = targetFunc(params, data.args[i])
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
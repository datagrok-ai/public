import {limitedMemoryBFGS} from "../../lbfgs/lbfgs"

type Likelihood = {
  value: number, 
  const: number, 
  mult: number
};

type FitResult = {
  parameters: number[],
  fittedCurve: (x:number)=> number,
  confidenceTop: (x:number)=> number,
  confidenceBottom: (x:number)=> number
};

export function fit (data:{x: number[], y: number[]}, params: number[],
                    curveFunction: (params: number[], x: number) => number, errorModel: string): FitResult {

  let of: (targetFunc: (params: number[], x: number) => number, 
  data: {y: number[], x: number[]},
  params: number[]) => Likelihood;
  switch(errorModel) {
    case 'constant':
      of = objectiveNormalConstant;
      break;
    case 'proportional':
      of = objectiveNormalProportional;
      break;
    default:
      of = objectiveNormalConstant;
      break;
  }

  function getDerivative(params: number[], selectedParam: number): number {
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


  let iterations = 0;

  let optimizable = {
    getValue: (parameters: number[]) => {
      return of(curveFunction, data, parameters).value;
    },
    getGradient: (parameters: number[], gradient: number[]) => {
      const length = Object.keys(parameters).length;
      iterations++;
      console.log(parameters.slice());
      
      for (let i = 0; i < parameters.length; i++)
        gradient[i] = getDerivative(parameters, i);

      return gradient;
    }
  };

  limitedMemoryBFGS(optimizable, params);
  limitedMemoryBFGS(optimizable, params);

  let fittedCurve = (x: number) =>{
    return curveFunction(params, x);
  }

  let top = (x: number) =>{
    let value = curveFunction(params, x);
    if (errorModel == "constant")
      return  value + 1.4*error;
    else
      return  value + 1.4*(Math.abs(value)*error);
  }

  let error = errorModel == "constant" ?
  of(curveFunction, data, params).const :
  of(curveFunction, data, params).mult;

  let bottom = (x: number) =>{
    let value = curveFunction(params, x);
    if (errorModel == "constant")
      return  value - 1.4*error;
    else
      return  value - 1.4*(Math.abs(value)*error);
  }

  let fitRes: FitResult = {
    parameters: params,
    fittedCurve: fittedCurve,
    confidenceTop: top,
    confidenceBottom: bottom
  };

  return fitRes;
}



// performFit(): FitResult {
//   limitedMemoryBFGS(this.optimizable, this.params);
//   limitedMemoryBFGS(this.optimizable, this.params);

//   let fit: FitResult = {
//     parameters = this.parameters;
//   };

//   return {parameters: , curveFunction: this.curveFunction, confidence: this.confidence};
// }


// interface ICurveFitter<TArg> {
//   data: {observed: number[], args: TArg[]};
//   params: number[];

//   //curve function
//   cf (params: number[], args: TArg): number
//   //objective function
//   of (targetFunc: (params: number[], args: TArg) => number, 
//                                   data: {observed: number[], args: TArg[]},
//                                   params: number[]): Likelihood
// }

// export class CurveFitter<TArg = number> implements ICurveFitter<TArg> {
//   errorModel: string;
//   data: {observed: number[], args: TArg[]};
//   params: number[];
//   cf;
//   of;
//   optimizable;
//   fitHistory: {iterations: number, parameters: {[iteration: number]: number[]}, 
//                 likelihood: {[iteration: number]: number}} = {
//     iterations: 0,
//     parameters: {},
//     likelihood: {}
//   };

//   constructor(curveFunction: (params: number[], args: TArg) => number, 
//               data: {observed: number[], args: TArg[]},
//               params: number[],
//               errorModel: 'constant' | 'proportional' = 'constant'//implement as enum
//               ) {

//     if(data.observed.length == 0 || data.observed.length != data.args.length)
//       throw new Error('invalid data');
      
//     this.data = data;
//     this.params = params;
//     this.errorModel = errorModel;
//     this.cf = curveFunction;

//     switch(errorModel) {
//       case 'constant':
//         this.of = objectiveNormalConstant;
//         break;
//       case 'proportional':
//         this.of = objectiveNormalProportional;
//         break;
//       default:
//         this.of = objectiveNormalConstant;
//         break;
//     }

//     this.optimizable = {
//       getValue: (parameters: number[]) => {
//         return this.of(this.cf, this.data, parameters).value;
//       },
//       getGradient: (parameters: number[], gradient: number[]) => {
//         const length = Object.keys(this.fitHistory.parameters).length;
//         this.fitHistory.iterations++;
//         this.fitHistory.parameters[length] = parameters.slice();
//         this.fitHistory.likelihood[length] = this.of(this.cf, this.data, parameters).value;
//         console.log(parameters.slice());
        
//         for (let i = 0; i < parameters.length; i++)
//           gradient[i] = this.getDerivative(parameters, i);

//         return gradient;
//       }
//     };
//   }

//   performFit(): void {
//     limitedMemoryBFGS(this.optimizable, this.params);
//     limitedMemoryBFGS(this.optimizable, this.params);
//   }

//   get parameters(): number[] {
//     return this.params;
//   }

//   get curveFunction(): Function {
//     let res = (arg: TArg) =>{
//       let params = this.fitHistory.parameters[this.fitHistory.iterations - 1];

//       return this.cf(params, arg);
//     }

//     return res;
//   }

//   get confidence(): {top:Function, bottom: Function} {
//     let params = this.fitHistory.parameters[this.fitHistory.iterations - 1];

//     let error = this.errorModel == "constant" ?
//                 this.of(this.cf, this.data, params).const :
//                 this.of(this.cf, this.data, params).mult;

//     let top = (arg: TArg) =>{
//       let value = this.cf(params, arg);
//       if (this.errorModel == "constant")
//         return  value + 1.4*error;
//       else
//         return  value + 1.4*(Math.abs(value)*error);
//     }

//     let bottom = (arg: TArg) =>{
//       let value = this.cf(params, arg);
//       if (this.errorModel == "constant")
//         return  value - 1.4*error;
//       else
//         return  value - 1.4*(Math.abs(value)*error);
//     }

//     return {top: top, bottom: bottom};
//   }

//   //central difference
//   private getDerivative(params: number[], selectedParam: number): number {
//     let step = params[selectedParam]*0.0001;
//     step = step == 0 ? 0.001 : step; 
//     let paramsTop: number[] = [];
//     let paramsBottom: number[] = [];
//     for (let i = 0; i < params.length; i++) {
//       if(i == selectedParam) {
//         paramsTop.push(params[i] + step);
//         paramsBottom.push(params[i] - step);
//       } else {
//         paramsTop.push(params[i]);
//         paramsBottom.push(params[i]);
//       }
//     }
//     const drvTop = this.of(this.cf, this.data, paramsTop).value;
//     const drvBottom = this.of(this.cf, this.data, paramsBottom).value;

//     return (drvTop - drvBottom)/(2*step);
//   }
// }

export function sigmoid(params: number[], x: number): number {
  const A = params[0];
  const B = params[1];
  const C = params[2];
  const D = params[3];
  const res = D + (A - D)/(1 + Math.pow(10, (x - C)*B));
  return res;
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



// function objectiveNormalConstant<TArg> (
//                                         targetFunc: (params: number[], args: TArg) => number, 
//                                         data: {observed: number[], args: TArg[]},
//                                         params: number[]
// ): Likelihood {
//   //assure observed and args same length
//   const pi = Math.PI;
//   let sigma = 0;
//   let sigmaSq = 0;
//   let likelihood = 0;

//   let residuesSquares = new Float32Array(data.args.length);
//   for(let i = 0; i < data.args.length; i++) {
//     const obs = data.observed[i];
//     const pred = targetFunc(params, data.args[i]);
//     residuesSquares[i] = Math.pow(obs - pred, 2);
//   }
  
//   for(let i = 0; i < residuesSquares.length; i++)
//     sigmaSq += residuesSquares[i];
//   sigmaSq /= residuesSquares.length;
//   sigma = Math.sqrt(sigmaSq);

//   for(let i = 0; i < residuesSquares.length; i++)
//     likelihood += residuesSquares[i]/sigmaSq + Math.log(2 * pi * sigmaSq);

//   return {value: -likelihood, const: sigma, mult: 0};                                              
// }

// function objectiveNormalProportional<TArg> (
//   targetFunc: (params: number[], args: TArg) => number, 
//   data: {observed: number[], args: TArg[]},
//   params: number[]
// ): Likelihood {
// //assure observed and args same length
// const pi = Math.PI;
// let sigma = 0;
// let sigmaSq = 0;
// let likelihood = 0;

// let residuesSquares = new Float32Array(data.args.length);
// for(let i = 0; i < data.args.length; i++) {
// const obs = data.observed[i];
// const pred = targetFunc(params, data.args[i])
// residuesSquares[i] = Math.pow(obs - pred, 2);
// }

// for(let i = 0; i < residuesSquares.length; i++)
// sigmaSq += residuesSquares[i];
// sigmaSq /= residuesSquares.length;
// sigma = Math.sqrt(sigmaSq);

// for(let i = 0; i < residuesSquares.length; i++)
// likelihood += residuesSquares[i]/sigmaSq + Math.log(2*pi*sigmaSq);

// return {value: -likelihood, const: sigma, mult: 0};                                              
// }

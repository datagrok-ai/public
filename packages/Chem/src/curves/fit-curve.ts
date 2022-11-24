import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as lbfgs from 'src/curves/lbfgs';



const dataX = [4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6, 6.25, 6.5, 6.75, 7, 7.25, 7.5, 7.75, 8];
const dataY = [1.005, 0.981, 0.921, 0.987, 0.948, 0.810, 0.730, 0.558, 0.501, 0.381, 0.271, 0.148, 0.139, 0.078, 0.095, -0.037];

type Likelihood = {
  value: number, 
  const: number, 
  mult: number
};

type ObjectiveFunction<TParam, TArg> = (targetFunc: (params: TParam, args: TArg) => number, 
                                        data: {observed: number[], args: TArg[]},
                                        params: TParam) => Likelihood

export function fitSigmoid(frame: DG.DataFrame) {
  //assure contains x and y


}

let objectiveNormalAdditive: ObjectiveFunction = ()=> {
  //assure observed and args same length
  const pi = Math.PI;
  let sigma = 0;
  let likelihood = 0;

  let residuesSquares = new Float32Array(data.args.length);
  for(let i = 0; i < data.args.length; i++)
    residuesSquares[i] = Math.pow(data.observed[i] - targetFunc(params, data.args[i]), 2);
  

  for(let i = 0; i < residuesSquares.length; i++)
    sigma += residuesSquares[i];
  sigma /= residuesSquares.length;
  let sigmaSq = sigma*sigma;

  for(let i = 0; i < residuesSquares.length; i++)
    likelihood += residuesSquares[i]/sigmaSq + Math.log(2*pi*sigmaSq);

  return {value: likelihood, const: sigma, mult: 0};                                              
}

export function fitCurve<TParam, TArg>(startingParameters: TParam, 
                         targetFunction: (params: TParam, args: TArg) => number, 
                         objectiveFunction: (targetFunc: (args: TArg) => number, data: TArg[]) => number): Function {
  
 
  
  return new Function;
}
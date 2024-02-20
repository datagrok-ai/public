
const INFINITY = Number.MAX_VALUE;
const PI = Math.PI;

export const FitErrorModel = {
  CONSTANT: 'constant',
  PROPORTIONAL: 'proportional',
};

export type Likelihood = {
  likelihood: number,
  likelihoodYs: number[],
  residualSquaresSums: number[]
};

export abstract class FitFunction {
  // abstract get name(): string;
  // abstract get parameterNames(): string[];
  abstract solve(params: Float32Array, xs: Float32Array[]): NumericalSolution;
  //abstract getInitialParameters(x: number[], y: number[]): number[];
}

export abstract class Optimizer {
  abstract optimize(
    params: Float32Array,
    objectiveFunc: (params: Float32Array) => Likelihood
  ): Likelihood
}

export type NumericalSolution = {
  ys: Float32Array [],
  devereges: boolean
}

function getResiduals(observed: Float32Array[], predicted: Float32Array[]): Float32Array[] {
  const dimYtypes = predicted.length;
  const residuals = new Array<Float32Array>(dimYtypes);

  for (let i = 0; i < dimYtypes; i++) {
    const dimY = predicted[i].length;
    residuals[i] = new Float32Array(dimY);
    for (let j = 0; j < dimY; j++)
      residuals[i][j] = observed[i][j] - predicted[i][j];
  }

  return residuals;
}

//sigma types number equal to data.ys.length, data.xs.length
export function getLikelihood(
  func: FitFunction, params: Float32Array, sigmaTypes: string[], sigmaValues: number[],
  data: {ys: Float32Array[], xs: Float32Array[]}): Likelihood {
  const solution = func.solve(params, data.xs);
  if (solution.devereges)
    return {likelihood: INFINITY, likelihoodYs: [], residualSquaresSums: []};

  const residuals = getResiduals(data.ys, solution.ys);
  const dimYtypes = data.ys.length;
  const residualSquaresSums = new Array<number>(dimYtypes);
  const likelihoodYs = new Array<number>(dimYtypes);
  let likelihood = 0;
  for (let i = 0; i < dimYtypes; i++) {
    const dimY = solution.ys[i].length;
    const sigmaSqare = sigmaValues[i]*sigmaValues[i];
    for (let j = 0; j < dimY; j++) {
      const square = residuals[i][j]*residuals[i][j];
      likelihoodYs[i] += square/sigmaSqare + Math.log(2*PI*sigmaSqare);
      if (sigmaTypes[i] == FitErrorModel.PROPORTIONAL)
        residualSquaresSums[i] += square/(solution.ys[i][j]*solution.ys[i][j]);
      else
        residualSquaresSums[i] += square;
    }
    likelihood += likelihoodYs[i];
  }

  return {likelihood, likelihoodYs, residualSquaresSums};
}

export function fit(func: FitFunction, params: Float32Array, sigmaTypes: string[],
  data: {ys: Float32Array[], xs: Float32Array[]}, optimizer: Optimizer) {
  let convergenceIndicator = INFINITY / 2;
  let convergenceIndicatorOld = INFINITY;
  const dimYtypes = data.ys.length;
  const sigmaValues = new Array<number>(dimYtypes).fill(1);

  const objectiveFunc = (params: Float32Array) => {
    return getLikelihood(func, params, sigmaTypes, sigmaValues, data);
  };

  while (Math.abs(convergenceIndicator - convergenceIndicatorOld) > 0.001) {
    const optimizedLikelihood = optimizer.optimize(params, objectiveFunc);

    for (let i = 0; i < dimYtypes; i++)
      sigmaValues[i] = Math.sqrt(optimizedLikelihood.residualSquaresSums[i]/data.ys[i].length);

    convergenceIndicatorOld = convergenceIndicator;
    convergenceIndicator = optimizedLikelihood.likelihood;
  }
}

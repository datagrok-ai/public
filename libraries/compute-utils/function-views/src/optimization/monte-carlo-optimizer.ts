/** Monte Carlo optimizer */

import {Extremum, OptimizationTask} from './definitions';

// export async function getMinimum(task: OptimizationTask): Promise<Extremum> {
//   // 1. Check dimensions
//   if (task.minVals.length !== task.maxVals.length)
//     throw new Error('Non-equal size of arrays with min & max values');

//   if (task.minVals.length === 0)
//     throw new Error('Empty arrays with min & max values');

//   // 2. Init optimization
//   const dim = task.minVals.length;
//   const argMin = new Float32Array(dim);
//   const arg = new Float32Array(dim);
//   const step = new Float32Array(dim);

//   task.minVals.forEach((val, idx) => {
//     argMin[idx] = val;
//     step[idx] = task.maxVals[idx] - val;
//   });

//   let minCost = await task.costFunc(argMin);
//   let cost: number;

//   // 3. Search for min
//   for (let i = 0; i < task.samplesCount; ++i) {
//     for (let j = 0; j < dim; ++j)
//       arg[j] = task.minVals[j] + step[j] * Math.random();

//     cost = await task.costFunc(arg);

//     if (cost < minCost) {
//       minCost = cost;
//       arg.forEach((val, idx) => argMin[idx] = val);
//     }
//   }

//   return {
//     arg: argMin,
//     cost: minCost,
//   };
// }

// costFunc: costFunc,
// minVals: minVals,
// maxVals: maxVals,
// samplesCount: this.samplesCount,


export async function optimizeNM (
  objectiveFunc: (x: Float32Array) => Promise<number>,
  paramsBottom: Float32Array,
  paramsTop: Float32Array,
  samplesCount: number): Promise<Extremum> {

  const params = new Float32Array(paramsBottom.length);
  for (let i = 0; i < paramsBottom.length; i++)
    params[i] =  (paramsBottom[i] + paramsTop[i])/2;


  const dim = params.length + 1;
  const dimParams = params.length;

  const optParams = new Array<Float32Array>(dim);
  const pointObjectives = new Array<number>(dim);

  for (let i = 0; i < dim; i++) {
    optParams[i] = new Float32Array(dimParams);
    for (let j = 0; j < dimParams; j++) {
      optParams[i][j] = params[j];
      if (i != 0) {
        if (params[i - 1] == 0)
          optParams[i][j] = 0.0001;
        else
          optParams[i][j] += 0.02 * params[i - 1];
      }
    }

    pointObjectives[i] = await objectiveFunc(optParams[i]);
  }

  const indexes = new Array<number>(dim);
  for (let i = 0; i < dim; i++)
    indexes[i] = i;

  const lastIndex = indexes.length - 1;

  let iteration = 0;
  const maxIter = 30;
  const infinitesemal = 5e-6;
  const tolerance = 1e-4;

  let best = 0;
  let previousBest = 0;
  let noImprovment = 0;

  const scaleReflection = 1;
  const scaleExpansion = 2;
  const scaleContraction = -0.5;

  const centroid = new Float32Array(dimParams);
  const reflectionPoint = new Float32Array(dimParams);
  const expansionPoint = new Float32Array(dimParams);
  const contractionPoint = new Float32Array(dimParams);


  if (dim > 1) {
    while (true) {
      indexes.sort((a:number, b:number) => {
        return pointObjectives[a] - pointObjectives[b];
      });
      if (iteration >= maxIter)
        break;

      if (iteration == 0) {
        best = pointObjectives[0];
        previousBest = 2*pointObjectives[indexes[0]];
      }

      iteration++;

      best = pointObjectives[indexes[0]];
      if ((best + infinitesemal)/(previousBest + infinitesemal) - 1 > tolerance)
        noImprovment = 0;
      else {
        ++noImprovment;
        if (noImprovment > 2 * dim)
          break;
      }

      previousBest = best;

      //centroid
      for (let i = 0; i < dimParams; i++)
        centroid[i] = params[i];
      for (let i = 0; i < dimParams; i++) {
        let val = 0;
        for (let j = 0; j < dim; j++) {
          if (j != indexes[lastIndex])
            val += optParams[j][i];
        }

        centroid[i] = val / (dim - 1);
      }

      // reflection
      for (let i = 0; i < dimParams; i++)
        reflectionPoint[i] = centroid[i];
      for (let i = 0; i < dimParams; i++)
        reflectionPoint[i] += scaleReflection * (centroid[i] - optParams[indexes[lastIndex]][i]);

      const reflectionScore = await objectiveFunc(reflectionPoint);

      // expansion
      if (reflectionScore < pointObjectives[indexes[lastIndex]]) {
        for (let i = 0; i < dimParams; i++)
          expansionPoint[i] = centroid[i];
        for (let i = 0; i < dimParams; i++)
          expansionPoint[i] += scaleExpansion * (centroid[i] - optParams[indexes[lastIndex]][i]);

        const expansionScore = await objectiveFunc(expansionPoint);


        if (expansionScore < reflectionScore) {
          pointObjectives[indexes[lastIndex]] = expansionScore;

          for (let i = 0; i < dimParams; i++)
            optParams[indexes[lastIndex]][i] = expansionPoint[i];

          continue;
        } else {
          pointObjectives[indexes[lastIndex]] = reflectionScore;

          for (let i = 0; i < dimParams; i++)
            optParams[indexes[lastIndex]][i] = reflectionPoint[i];

          continue;
        }
      }

      // Contraction
      for (let i = 0; i < dimParams; i++)
        contractionPoint[i] = centroid[i];
      for (let i = 0; i < dimParams; i++)
        contractionPoint[i] += scaleContraction * (centroid[i] - optParams[indexes[lastIndex]][i]);

      const contractionScore = await objectiveFunc(contractionPoint);

      if (contractionScore < pointObjectives[indexes[lastIndex]]) {
        pointObjectives[indexes[lastIndex]] = contractionScore;

        for (let i = 0; i < dimParams; i++)
          optParams[indexes[lastIndex]][i] = contractionPoint[i];

        continue;
      }

      break;
    }
  }

  for (let i = 0; i < dimParams; i++)
    params[i] = optParams[indexes[0]][i];

  return {arg: params, cost: pointObjectives[0]};
}

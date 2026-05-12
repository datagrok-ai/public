// GENERATED — do not edit by hand.
// Run `npm run update-codegen` to regenerate.
// Source: ./optimizer-nelder-mead.ts
import {Extremum} from './optimizer-misc';
import {fillCentroid, fillPoint} from './optimizer-nelder-mead';

function getInitialParams(
  objectiveFunc: (x: Float64Array) => number|undefined,
  settings: Map<string, number>,
  paramsInitial: Float64Array,
  costOutside: number,
): [Float64Array[], number[]] {
  const dim = paramsInitial.length + 1;
  const dimParams = paramsInitial.length;
  const nonZeroParam = settings.get('nonZeroParam')!;
  const initScale = settings.get('initialScale')!;

  const optParams = new Array<Float64Array>(dim);
  const pointObjectives = new Array<number>(dim);

  for (let i = 0; i < dim; i++) {
    optParams[i] = new Float64Array(dimParams);
    for (let j = 0; j < dimParams; j++) {
      optParams[i][j] = paramsInitial[j];
      if (i != 0 && i - 1 === j) {
        if (paramsInitial[j] == 0)
          optParams[i][j] = nonZeroParam;
        else
          optParams[i][j] += initScale * paramsInitial[i - 1];
      }
    }

    pointObjectives[i] = objectiveFunc(optParams[i]) ?? costOutside;
  }

  return [optParams, pointObjectives];
}

export const optimizeNMSync = function(
  objectiveFunc: (x: Float64Array) => number|undefined,
  paramsInitial: Float64Array,
  settings: Map<string, number>,
  threshold?: number,
) : Extremum {
  // Settings initialization
  const tolerance = settings.get('tolerance')!;
  const maxIter = settings.get('maxIter')!;
  const scaleReflection = settings.get('scaleReflaction')!;
  const scaleExpansion = settings.get('scaleExpansion')!;
  const scaleContraction = settings.get('scaleContraction')!;

  const dim = paramsInitial.length + 1;
  const dimParams = paramsInitial.length;

  const costOutside = 2*(objectiveFunc(paramsInitial) ?? Infinity);

  const [optParams, pointObjectives] =
    getInitialParams(objectiveFunc, settings, paramsInitial, costOutside);

  const indexes = new Array<number>(dim);
  for (let i = 0; i < dim; i++)
    indexes[i] = i;

  const lastIndex = indexes.length - 1;

  let iteration = 0;
  let best = 0;
  let previousBest = 0;
  let noImprovment = 0;

  const centroid = new Float64Array(dimParams);
  const reflectionPoint = new Float64Array(dimParams);
  const expansionPoint = new Float64Array(dimParams);
  const contractionPoint = new Float64Array(dimParams);
  const costs = new Array<number>(maxIter);

  if (dim > 1) {
    while (true) {
      indexes.sort((a:number, b:number) => {
        return pointObjectives[a] - pointObjectives[b];
      });
      if (iteration > maxIter)
        break;

      if (iteration == 0) {
        best = pointObjectives[0];
        previousBest = 2*pointObjectives[indexes[0]];
      }
      costs[iteration] = best;

      ++iteration;

      best = pointObjectives[indexes[0]];
      if (previousBest - best > tolerance)
        noImprovment = 0;
      else {
        ++noImprovment;
        if (noImprovment > 2 * dim)
          break;
      }

      if (threshold != null) {
        if (best <= threshold)
          break;
      }

      previousBest = best;

      //centroid
      fillCentroid(centroid, dimParams, indexes[lastIndex], optParams);

      // reflection
      fillPoint(centroid, reflectionPoint, indexes[lastIndex],
        optParams, scaleReflection, dimParams);
      const reflectionScore = objectiveFunc(reflectionPoint) ?? costOutside;

      // expansion
      if (reflectionScore < pointObjectives[indexes[lastIndex]]) {
        fillPoint(centroid, expansionPoint, indexes[lastIndex],
          optParams, scaleExpansion, dimParams);

        const expansionScore = objectiveFunc(expansionPoint) ?? costOutside;

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
      fillPoint(centroid, contractionPoint, indexes[lastIndex],
        optParams, scaleContraction, dimParams);

      const contractionScore = objectiveFunc(contractionPoint) ?? costOutside;

      if (contractionScore < pointObjectives[indexes[lastIndex]]) {
        pointObjectives[indexes[lastIndex]] = contractionScore;

        for (let i = 0; i < dimParams; i++)
          optParams[indexes[lastIndex]][i] = contractionPoint[i];

        continue;
      }

      break;
    }

    for (let i = iteration; i < maxIter; i++)
      costs[i] = pointObjectives[indexes[0]];
  }

  return {
    point: optParams[indexes[0]],
    cost: pointObjectives[indexes[0]],
    iterCosts: costs,
    iterCount: iteration - 1,
  };
};

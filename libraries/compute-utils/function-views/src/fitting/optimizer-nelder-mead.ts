import {Extremum, IOptimizer, Setting} from './optimizer-misc';

/** The Nelder Mead seetings vals */
export const nelderMeadSettingsOpts = new Map<string, Setting>([
  ['tolerance', {
    default: 0.000001,
    min: 1e-20,
    max: 1e-1,
    caption: 'tolerance',
    tooltipText: 'How precise the result should be. Lower value = more accurate result, but longer computation',
    inputType: 'Float',
  }],
  ['maxIter', {
    default: 50,
    min: 1,
    max: 10000,
    caption: 'max iterations',
    tooltipText: 'Maximum number of iterations. Higher = better fit',
    inputType: 'Int',
  }],
  ['nonZeroParam', {
    default: 0.0001,
    min: 1e-20,
    max: 1e-1,
    caption: 'non-zero param',
    tooltipText: 'Minimum parameter value (to avoid zero values)',
    inputType: 'Float',
  }],
  ['initialScale', {
    default: 0.02,
    min: 1e-20,
    max: 1e-1,
    caption: 'initial scale',
    tooltipText: 'Size of the initial search area. Higher value = wider search area',
    inputType: 'Float',
  }],
  ['scaleReflaction', {
    default: 1,
    min: 1,
    max: 2,
    caption: 'scale reflection',
    tooltipText: 'How far the algorithm \'bounces back\' from failed attempts',
    inputType: 'Float',
  }],
  ['scaleExpansion', {
    default: 2,
    min: 1,
    max: 2,
    caption: 'scale expansion',
    tooltipText: 'How much the algorithm expands search in a promising direction',
    inputType: 'Float',
  }],
  ['scaleContraction', {
    default: -0.5,
    min: -0.5,
    max: 0,
    caption: 'scale contraction',
    tooltipText: 'How much the algorithm narrows the search area when approaching the result',
    inputType: 'Float',
  }],
]);

async function getInitialParams(
  objectiveFunc: (x: Float64Array) => Promise<number|undefined>,
  settings: Map<string, number>,
  paramsInitial: Float64Array,
  costOutside: number
): Promise<[Float64Array[], number[]]> {
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

    pointObjectives[i] = await objectiveFunc(optParams[i]) ?? costOutside;
  }

  return [optParams, pointObjectives];
}

function fillCentroid(centroid: Float64Array, dimParams: number, lastIndex: number, optParams: Float64Array[]) {
  for (let i = 0; i < dimParams; i++) {
    let val = 0;
    for (let j = 0; j < dimParams + 1; j++) {
      if (j != lastIndex)
        val += optParams[j][i];
    }

    centroid[i] = val / dimParams;
  }
}

function fillPoint(
  centroid: Float64Array, point: Float64Array,
  lastIndex: number, optParams: Float64Array[],
  scale: number, dimParams: number,
) {
  for (let i = 0; i < dimParams; i++) {
    point[i] = centroid[i];
    point[i] += scale * (centroid[i] - optParams[lastIndex][i]);

  }
}

export const optimizeNM: IOptimizer = async function(
  objectiveFunc: (x: Float64Array) => Promise<number|undefined>,
  paramsInitial: Float64Array,
  settings: Map<string, number>,
  threshold?: number,
) : Promise<Extremum> {
  // Settings initialization
  const tolerance = settings.get('tolerance')!;
  const maxIter = settings.get('maxIter')!;
  const scaleReflection = settings.get('scaleReflaction')!;
  const scaleExpansion = settings.get('scaleExpansion')!;
  const scaleContraction = settings.get('scaleContraction')!;

  const dim = paramsInitial.length + 1;
  const dimParams = paramsInitial.length;

  const costOutside = 2*(await objectiveFunc(paramsInitial) ?? Infinity);

  const [optParams, pointObjectives] =
    await getInitialParams(objectiveFunc, settings, paramsInitial, costOutside);

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
      const reflectionScore = await objectiveFunc(reflectionPoint) ?? costOutside;

      // expansion
      if (reflectionScore < pointObjectives[indexes[lastIndex]]) {
        fillPoint(centroid, expansionPoint, indexes[lastIndex],
          optParams, scaleExpansion, dimParams);

        const expansionScore = await objectiveFunc(expansionPoint) ?? costOutside;

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

      const contractionScore = await objectiveFunc(contractionPoint) ?? costOutside;

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

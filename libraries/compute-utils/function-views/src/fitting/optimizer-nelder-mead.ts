import {Extremum, IOptimizer, Setting} from './optimizer-misc';

/** The Nelder Mead seetings vals */
export const nelderMeadSettingsVals = new Map<string, Setting>([
  ['tolerance', {default: 0.000001, min: 1e-20, max: 1e-1}],
  ['maxIter', {default: 50, min: 1, max: 10000}],
  ['nonZeroParam', {default: 0.0001, min: 1e-20, max: 1e-1}],
  ['initialScale', {default: 0.02, min: 1e-20, max: 1e-1}],
  ['scaleReflaction', {default: 1, min: 1, max: 2}],
  ['scaleExpansion', {default: 2, min: 1, max: 2}],
  ['scaleContraction', {default: -0.5, min: -0.5, max: 0}],
]);

/** Captions for the Nelder-Mead method settings */
export const nelderMeadCaptions = new Map([
  ['tolerance', 'tolerance'],
  ['maxIter', 'max iterations'],
  ['nonZeroParam', 'non-zero param'],
  ['initialScale', 'initial scale'],
  ['scaleReflaction', 'scale reflection'],
  ['scaleExpansion', 'scale expansion'],
  ['scaleContraction', 'scale contraction'],
]);

async function getInitialParams(
  objectiveFunc: (x: Float32Array) => Promise<number>,
  settings: Map<string, number>,
  paramsInitial: Float32Array,
  restrictionsBottom: Float32Array,
  restrictionsTop: Float32Array): Promise<[Float32Array[], number[]]> {
  const dim = paramsInitial.length + 1;
  const dimParams = paramsInitial.length;
  const nonZeroParam = settings.get('nonZeroParam')!;
  const initScale = settings.get('initialScale')!;

  const optParams = new Array<Float32Array>(dim);
  const pointObjectives = new Array<number>(dim);

  for (let i = 0; i < dim; i++) {
    optParams[i] = new Float32Array(dimParams);
    for (let j = 0; j < dimParams; j++) {
      optParams[i][j] = paramsInitial[j];
      if (i != 0 && i - 1 === j) {
        if (paramsInitial[j] == 0)
          optParams[i][j] = nonZeroParam;
        else
          optParams[i][j] += initScale * paramsInitial[i - 1];

        if (optParams[i][j] < restrictionsBottom[j])
          optParams[i][j] = restrictionsBottom[j];
        else if (optParams[i][j] > restrictionsTop[j])
          optParams[i][j] = restrictionsTop[j];
      }
    }

    pointObjectives[i] = await objectiveFunc(optParams[i]);
  }

  return [optParams, pointObjectives];
}

function fillCentroid(centroid: Float32Array, dimParams: number, lastIndex: number, optParams: Float32Array[]) {
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
  centroid: Float32Array, point: Float32Array,
  lastIndex: number, optParams: Float32Array[],
  scale: number, dimParams: number,
  restrictionsBottom: Float32Array,
  restrictionsTop: Float32Array) {
  for (let i = 0; i < dimParams; i++) {
    point[i] = centroid[i];
    point[i] += scale * (centroid[i] - optParams[lastIndex][i]);

    if (point[i] < restrictionsBottom[i])
      point[i] = restrictionsBottom[i];
    else if (point[i] > restrictionsTop[i])
      point[i] = restrictionsTop[i];
  }
}

export const optimizeNM:IOptimizer = async function(
  objectiveFunc: (x: Float32Array) => Promise<number>,
  paramsInitial: Float32Array,
  settings: Map<string, number>,
  restrictionsBottom: Float32Array,
  restrictionsTop: Float32Array) : Promise<Extremum> {
  // Settings initialization
  const tolerance = settings.get('tolerance')!;
  const maxIter = settings.get('maxIter')!;
  const scaleReflection = settings.get('scaleReflaction')!;
  const scaleExpansion = settings.get('scaleExpansion')!;
  const scaleContraction = settings.get('scaleContraction')!;

  const dim = paramsInitial.length + 1;
  const dimParams = paramsInitial.length;

  const [optParams, pointObjectives] =
    await getInitialParams(objectiveFunc, settings, paramsInitial, restrictionsBottom, restrictionsTop);

  const indexes = new Array<number>(dim);
  for (let i = 0; i < dim; i++)
    indexes[i] = i;

  const lastIndex = indexes.length - 1;

  let iteration = 0;
  let best = 0;
  let previousBest = 0;
  let noImprovment = 0;

  const centroid = new Float32Array(dimParams);
  const reflectionPoint = new Float32Array(dimParams);
  const expansionPoint = new Float32Array(dimParams);
  const contractionPoint = new Float32Array(dimParams);
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

      previousBest = best;

      //centroid
      fillCentroid(centroid, dimParams, indexes[lastIndex], optParams);

      // reflection
      fillPoint(centroid, reflectionPoint, indexes[lastIndex],
        optParams, scaleReflection, dimParams, restrictionsBottom, restrictionsTop);
      const reflectionScore = await objectiveFunc(reflectionPoint);

      // expansion
      if (reflectionScore < pointObjectives[indexes[lastIndex]]) {
        fillPoint(centroid, expansionPoint, indexes[lastIndex],
          optParams, scaleExpansion, dimParams, restrictionsBottom, restrictionsTop);

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
      fillPoint(centroid, contractionPoint, indexes[lastIndex],
        optParams, scaleContraction, dimParams, restrictionsBottom, restrictionsTop);

      const contractionScore = await objectiveFunc(contractionPoint);

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

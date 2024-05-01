import {Extremum, IOptimizer} from './optimizer-misc';

/** The Nelder-Mead method settings */
export type NelderMeadSettings = {
  tolerance: number,
  maxIter: number,
  nonZeroParam: number,
  initialScale: number,
  scaleReflaction: number,
  scaleExpansion: number,
  scaleContraction: number,
};

/** Default values for the Nelder-Mead method */
export enum NELDER_MEAD_DEFAULTS {
  TOLERANCE = 0.000001,
  MAX_ITER = 50,
  NON_ZERO_PARAM = 0.0001,
  INITIAL_SCALE = 0.02,
  SCALE_REFLECTION = 1,
  SCALE_EXPANSION = 2,
  SCALE_CONTRACTION = -0.5,
};

export const optimizeNM:IOptimizer = async function(
  objectiveFunc: (x: Float32Array) => Promise<number>,
  paramsInitial: Float32Array,
  settings: NelderMeadSettings) : Promise<Extremum> {
  // Settings initialization
  const tolerance = settings.tolerance;
  const maxIter = settings.maxIter;
  const nonZeroParam = settings.nonZeroParam;
  const initScale = settings.initialScale;
  const scaleReflection = settings.scaleReflaction;
  const scaleExpansion = settings.scaleExpansion;
  const scaleContraction = settings.scaleContraction;

  const dim = paramsInitial.length + 1;
  const dimParams = paramsInitial.length;

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
      }
    }

    pointObjectives[i] = await objectiveFunc(optParams[i]);
  }

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
      for (let i = 0; i < dimParams; i++)
        centroid[i] = paramsInitial[i];
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

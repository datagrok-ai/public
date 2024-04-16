import {Optimizer, FitFunction, getLikelihood, Likelihood} from './optimizer';

class OptimizerNelder implements Optimizer {
  optimize(
    params: Float32Array,
    objectiveFunc: (params: Float32Array) => Likelihood) : Likelihood {
    const dim = params.length + 1;
    const dimParams = params.length;

    const optParams = new Array<Float32Array>(dim);
    const pointLikelihoods = new Array<Likelihood>(dim);

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

      pointLikelihoods[i] = objectiveFunc(optParams[i]);
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

    // const reflectionScore: Likelihood;
    // const expansionScore: Likelihood;
    // const contractionScore: Likelihood;

    const centroid = new Float32Array(dimParams);
    const reflectionPoint = new Float32Array(dimParams);
    const expansionPoint = new Float32Array(dimParams);
    const contractionPoint = new Float32Array(dimParams);


    if (dim > 1) {
      while (true) {
        indexes.sort((a:number, b:number) => {
          return pointLikelihoods[a].likelihood - pointLikelihoods[b].likelihood;
        });
        if (iteration >= maxIter)
          break;

        if (iteration == 0) {
          best = pointLikelihoods[0].likelihood;
          previousBest = 2*pointLikelihoods[indexes[0]].likelihood;
        }

        iteration++;

        best = pointLikelihoods[indexes[0]].likelihood;
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

        const reflectionScore = objectiveFunc(reflectionPoint);

        // expansion
        if (reflectionScore.likelihood < pointLikelihoods[indexes[lastIndex]].likelihood) {
          for (let i = 0; i < dimParams; i++)
            expansionPoint[i] = centroid[i];
          for (let i = 0; i < dimParams; i++)
            expansionPoint[i] += scaleExpansion * (centroid[i] - optParams[indexes[lastIndex]][i]);

          const expansionScore = objectiveFunc(expansionPoint);


          if (expansionScore.likelihood < reflectionScore.likelihood) {
            pointLikelihoods[indexes[lastIndex]] = expansionScore;

            for (let i = 0; i < dimParams; i++)
              optParams[indexes[lastIndex]][i] = expansionPoint[i];

            continue;
          } else {
            pointLikelihoods[indexes[lastIndex]] = reflectionScore;

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

        const contractionScore = objectiveFunc(contractionPoint);

        if (contractionScore.likelihood < pointLikelihoods[indexes[lastIndex]].likelihood) {
          pointLikelihoods[indexes[lastIndex]] = contractionScore;

          for (let i = 0; i < dimParams; i++)
            optParams[indexes[lastIndex]][i] = contractionPoint[i];

          continue;
        }

        break;
      }
    }

    for (let i = 0; i < dimParams; i++)
      params[i] = optParams[indexes[0]][i];

    return pointLikelihoods[0];
  }
}

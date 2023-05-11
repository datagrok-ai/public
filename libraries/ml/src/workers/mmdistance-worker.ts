import {Matrix} from '@datagrok-libraries/utils/src/type-declarations';
import {calcDistanceMatrix} from '@datagrok-libraries/utils/src/vector-operations';
import {mmDistanceFunctions, MmDistanceFunctionsNames} from '../macromolecule-distance-functions';


onmessage = (event) => {
  const {values, fnName, options} = event.data;
  const data: { error?: any; distanceMatrixData?: Matrix } = {};
  try {
    const recalculatedDistances = calcDistanceMatrix(
      values, mmDistanceFunctions[fnName as MmDistanceFunctionsNames](options),
    );
    // normalize values
    let min = Number.MAX_VALUE;
    let max = Number.MIN_VALUE;
    for (let i = 0; i < values.length - 1; i++) {
      for (let j = i; j < values.length; j++) {
        const value = recalculatedDistances[i][j];
        if (value > max) max = value;
        if (value < min) min = value;
      }
    }

    for (let i = 0; i < values.length - 1; i++) {
      for (let j = i; j < values.length; j++)
        recalculatedDistances[i][j] = (recalculatedDistances[i][j] - min) / (max - min);
    }

    data.distanceMatrixData = recalculatedDistances;
  } catch (e) {
    data.error = e;
  }
  postMessage(data);
};

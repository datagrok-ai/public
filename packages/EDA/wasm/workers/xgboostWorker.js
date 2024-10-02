// Worker for XGBoost training

import {XGBoost} from '../../wasm/XGBoostAPIinWebWorker';

// Data constants
const INT_BYTES = 4;
const FLOAT_BYTES = 4;
const SIZE_IDX = 0;

onmessage = async function(evt) {
  // eslint-disable-next-line new-cap
  XGBoost().then((booster) => {
    const features = evt.data.features;
    const target = evt.data.target;
    const samplesCount = evt.data.samplesCount;
    const featuresCount = evt.data.featuresCount;
    const modelReserve = evt.data.modelReserve;
    const utilsLength = evt.data.utilsLength;
    const iterations = evt.data.iterations;
    const eta = evt.data.eta;
    const maxDepth = evt.data.maxDepth;
    const lambda = evt.data.lambda;
    const alpha = evt.data.alpha;
    const missingValue = evt.data.missingValue;

    // Allocate memory
    const featuresBuf = booster._malloc(samplesCount * featuresCount * FLOAT_BYTES);
    const targetBuf = booster._malloc(samplesCount * FLOAT_BYTES);
    const modelBuf = booster._malloc(modelReserve * INT_BYTES);
    const utilsBuf = booster._malloc(utilsLength * INT_BYTES);

    // Put features to wasm buffer
    const floatHeap = booster.HEAPF32;
    const featuresSize = samplesCount * featuresCount;
    for (let j = 0; j < featuresSize; ++j)
      floatHeap[featuresBuf / FLOAT_BYTES + j] = features[j];

    // Put targets to wasm buffer
    for (let i = 0; i < samplesCount; ++i)
      floatHeap[targetBuf / FLOAT_BYTES + i] = target[i];

    // Train model
    booster._train(
      featuresBuf, samplesCount, featuresCount, missingValue, // features data
      targetBuf, samplesCount, // target data
      iterations, eta, maxDepth, lambda, alpha, // hyperparameters
      utilsBuf, utilsLength, // utils
      modelBuf, modelReserve, // model params to be trained
    );

    // Extract model params from wasm buffer
    const intHeap = booster.HEAP32;
    const paramsCount = intHeap[utilsBuf / INT_BYTES + SIZE_IDX];
    const params = new Int32Array(paramsCount);

    for (let i = 0; i < paramsCount; ++i)
      params[i] = intHeap[modelBuf / INT_BYTES + i];

    // Free allocated memory
    booster._free(featuresBuf);
    booster._free(targetBuf);
    booster._free(utilsBuf);
    booster._free(modelBuf);

    postMessage({'params': params});
  });
};

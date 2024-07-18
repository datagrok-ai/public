// JavaScript API for call wasm-functions from the XGBoostAPI module

const INT_BYTES = 4;
const FLOAT_BYTES = 4;

export async function initXgboost() {
  await initXGBoostModule();
}

/** Train XGBoost model */
export function train(features, samplesCount, featuresCount, missingValue,
  labels,
  iterations, eta, maxDepth, lambda, alpha,
  utils,
  modelBytes, reserved) {
  return XGBoostModule._train(features.byteOffset, samplesCount, featuresCount, missingValue,
    labels.byteOffset, samplesCount,
    iterations, eta, maxDepth, lambda, alpha,
    utils.byteOffset, utils.length,
    modelBytes.byteOffset, reserved);
}

/** Predict by XGBoost model */
export function predict(features, samplesCount, featuresCount, missingValue,
  modelBytes, modelSize,
  predictions) {
  return XGBoostModule._predict(features.byteOffset, samplesCount, featuresCount, missingValue,
    modelBytes.byteOffset, modelSize,
    predictions.byteOffset, samplesCount);
}

/** Allocate memory for training */
export function allocTrainMemory(samplesCount, featuresCount, modelReserve, utilsLength) {
  return {
    int32Buffer: XGBoostModule.HEAP32.buffer,
    float32Buffer: XGBoostModule.HEAPF32.buffer,
    featuresOffset: XGBoostModule._malloc(samplesCount * featuresCount * FLOAT_BYTES),
    labelsOffset: XGBoostModule._malloc(samplesCount * FLOAT_BYTES),
    modelOffset: XGBoostModule._malloc(modelReserve * INT_BYTES),
    utilsOffset: XGBoostModule._malloc(utilsLength * INT_BYTES),
  };
}

/** Free memory allocated for training */
export function freeTrainMemory(trainMemory) {
  XGBoostModule._free(trainMemory.featuresOffset);
  XGBoostModule._free(trainMemory.labelsOffset);
  XGBoostModule._free(trainMemory.modelOffset);
  XGBoostModule._free(trainMemory.utilsOffset);
}

/** Allocate memory for predicting */
export function allocPredictMemory(samplesCount, featuresCount, modelReserve) {
  return {
    int32Buffer: XGBoostModule.HEAP32.buffer,
    float32Buffer: XGBoostModule.HEAPF32.buffer,
    featuresOffset: XGBoostModule._malloc(samplesCount * featuresCount * FLOAT_BYTES),
    predictOffset: XGBoostModule._malloc(samplesCount * FLOAT_BYTES),
    modelOffset: XGBoostModule._malloc(modelReserve * INT_BYTES),
  };
}

/** Free memory allocated for predicting */
export function freePredictMemory(predictMemory) {
  XGBoostModule._free(predictMemory.featuresOffset);
  XGBoostModule._free(predictMemory.predictOffset);
  XGBoostModule._free(predictMemory.modelOffset);
}

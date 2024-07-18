// JavaScript API for call wasm-functions from the XGBoostAPI module

const INT_BYTES = 4;
const FLOAT_BYTES = 4;

export async function initXgboost() {
  await initXGBoostModule();
}

/** Returns buffer & offset: the Int32Array case */
export function memAllocInt32Arr(size) {
  return {
    buf: XGBoostModule.HEAP32.buffer,
    off: XGBoostModule._malloc(size * INT_BYTES),
  };
}

/** Returns buffer & offset: the Float32Array case */
export function memAllocFloat32Arr(size) {
  return {
    buf: XGBoostModule.HEAPF32.buffer,
    off: XGBoostModule._malloc(size * FLOAT_BYTES),
  };
}

/** Free memory */
export function memFree(ptr) {
  XGBoostModule._free(ptr);
}

export function train(features, samplesCount, featuresCount, MISSING_VALUE,
  labels, labelsLength,
  iterations, eta, maxDepth, lambda, alpha,
  modelSizePtr,
  modelBytes, reserved) {
  return XGBoostModule._train(features.byteOffset, samplesCount, featuresCount, MISSING_VALUE,
    labels.byteOffset, labelsLength,
    iterations, eta, maxDepth, lambda, alpha,
    modelSizePtr.byteOffset, 1,
    modelBytes.byteOffset, reserved);
}

export function predict(features, samplesCount, featuresCount, MISSING_VALUE,
  modelBytes, modelSize,
  predictions, predictionsLength) {
  return XGBoostModule._predict(features.byteOffset, samplesCount, featuresCount, MISSING_VALUE,
    modelBytes.byteOffset, modelSize,
    predictions.byteOffset, predictionsLength);
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

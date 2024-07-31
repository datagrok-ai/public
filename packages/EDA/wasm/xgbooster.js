// JavaScript API for call wasm-functions from the XGBoostAPI module

// Data constants
const INT_BYTES = 4;
const FLOAT_BYTES = 4;
const SIZE_IDX = 0;

export async function initXgboost() {
  await initXGBoostModule();
}

/** Fit and return model params */
export function fit(features, target, missingValue, iterations, eta, maxDepth, lambda, alpha,
  modelReserve, utilsLength) {
  // Data size
  const samplesCount = target.length;
  const featuresCount = features.length;

  // Allocate memory
  const featuresBuf = XGBoostModule._malloc(samplesCount * featuresCount * FLOAT_BYTES);
  const targetBuf = XGBoostModule._malloc(samplesCount * FLOAT_BYTES);
  const modelBuf = XGBoostModule._malloc(modelReserve * INT_BYTES);
  const utilsBuf = XGBoostModule._malloc(utilsLength * INT_BYTES);

  // Wasm buffer
  const floatHeap = XGBoostModule.HEAPF32;
  const intHeap = XGBoostModule.HEAP32;

  let raw;

  // Put features to wasm buffer
  for (let j = 0; j < featuresCount; ++j) {
    raw = features.byIndex(j).getRawData();

    for (let i = 0; i < samplesCount; ++i)
      floatHeap[featuresBuf / FLOAT_BYTES + i + j * samplesCount] = raw[i];
  }

  // Put targets to wasm buffer
  raw = target.getRawData();
  for (let i = 0; i < samplesCount; ++i)
    floatHeap[targetBuf / FLOAT_BYTES + i] = raw[i];

  // Train model
  XGBoostModule._train(
    featuresBuf, samplesCount, featuresCount, missingValue, // features data
    targetBuf, samplesCount, // target data
    iterations, eta, maxDepth, lambda, alpha, // hyperparameters
    utilsBuf, utilsLength, // utils
    modelBuf, modelReserve, // model params to be trained
  );

  // Extract model params from wasm buffer
  const paramsCount = intHeap[utilsBuf / INT_BYTES + SIZE_IDX];
  const params = new Int32Array(paramsCount);

  for (let i = 0; i < paramsCount; ++i)
    params[i] = intHeap[modelBuf / INT_BYTES + i];

  // Free allocated memory
  XGBoostModule._free(featuresBuf);
  XGBoostModule._free(targetBuf);
  XGBoostModule._free(utilsBuf);
  XGBoostModule._free(modelBuf);

  return params;
} // fit

/** Return prediction by trained model */
export function predict(features, missingValue, params) {
  // Data & model sizes
  const samplesCount = features.byIndex(0).length;
  const featuresCount = features.length;
  const paramsCount = params.length;

  // Wasm buffer
  const floatHeap = XGBoostModule.HEAPF32;
  const intHeap = XGBoostModule.HEAP32;

  // Allocate memory
  const featuresBuf = XGBoostModule._malloc(samplesCount * featuresCount * FLOAT_BYTES);
  const targetBuf = XGBoostModule._malloc(samplesCount * FLOAT_BYTES);
  const modelBuf = XGBoostModule._malloc(paramsCount * INT_BYTES);

  // Put features to wasm buffer
  for (let j = 0; j < featuresCount; ++j) {
    const raw = features.byIndex(j).getRawData();

    for (let i = 0; i < samplesCount; ++i)
      floatHeap[featuresBuf / FLOAT_BYTES + i + j * samplesCount] = raw[i];
  }

  // Put model to wasm bufffer
  for (let i = 0; i < paramsCount; ++i)
    intHeap[modelBuf / INT_BYTES + i] = params[i];

  // Compute predictions
  XGBoostModule._predict(
    featuresBuf, samplesCount, featuresCount, missingValue, // features
    modelBuf, paramsCount, // model params
    targetBuf, samplesCount, // target to be predicted
  );

  // Extract predictions from wasm buffer
  const prediction = new Float32Array(samplesCount);

  for (let i = 0; i < samplesCount; ++i)
    prediction[i] = floatHeap[targetBuf / FLOAT_BYTES + i];

  // Free allocated memory
  XGBoostModule._free(featuresBuf);
  XGBoostModule._free(targetBuf);
  XGBoostModule._free(modelBuf);

  return prediction;
} // predict


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
export function predictOld(features, samplesCount, featuresCount, missingValue,
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

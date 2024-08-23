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

  // Wasm buffer routine
  const floatHeap = XGBoostModule.HEAPF32;
  let intHeap = XGBoostModule.HEAP32;
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
  intHeap = XGBoostModule.HEAP32;
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

/** Fit and return model params in webworker */
export async function fitInWebWorker(features, target, missingValue, iterations, eta, maxDepth, lambda, alpha,
  modelReserve, utilsLength) {
  return new Promise((resolve, reject) => {
    // Data size
    const samplesCount = target.length;
    const featuresCount = features.length;

    // Features raw data
    const featuresRaw = new Float32Array(samplesCount * featuresCount);
    let shift;
    let raw;
    for (let j = 0; j < featuresCount; ++j) {
      raw = features.byIndex(j).getRawData();
      shift = j * samplesCount;

      for (let i = 0; i < samplesCount; ++i)
        featuresRaw[i + shift] = raw[i];
    }

    const worker = new Worker(new URL('../wasm/workers/xgboostWorker.js', import.meta.url));

    worker.postMessage({
      features: featuresRaw,
      target: target.getRawData(),
      samplesCount: samplesCount,
      featuresCount: featuresCount,
      modelReserve: modelReserve,
      utilsLength: utilsLength,
      iterations: iterations,
      eta: eta,
      maxDepth: maxDepth,
      lambda: lambda,
      alpha: alpha,
      missingValue: missingValue,
    });

    worker.onmessage = function(e) {
      worker.terminate();
      resolve(e.data.params);
    };
  });
} // fitInWebWorker

/** Return prediction by trained model */
export function predict(features, missingValue, params) {
  // Data & model sizes
  const samplesCount = features.byIndex(0).length;
  const featuresCount = features.length;
  const paramsCount = params.length;

  // Wasm buffer routine
  let floatHeap = XGBoostModule.HEAPF32;
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
  floatHeap = XGBoostModule.HEAPF32;
  const prediction = new Float32Array(samplesCount);

  for (let i = 0; i < samplesCount; ++i)
    prediction[i] = floatHeap[targetBuf / FLOAT_BYTES + i];

  // Free allocated memory
  XGBoostModule._free(featuresBuf);
  XGBoostModule._free(targetBuf);
  XGBoostModule._free(modelBuf);

  return prediction;
} // predict

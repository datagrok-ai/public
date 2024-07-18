import {memAllocInt32Arr, memAllocFloat32Arr, memFree, train, predict,
  allocTrainMemory, freeTrainMemory, allocPredictMemory, freePredictMemory} from '../wasm/xgbooster';

export function basicTestXGBoost() {
  const iterations = 20;
  const eta = 0.3;
  const maxDepth = 6;
  const lambda = 1;
  const alpha = 0;

  const ROWS = 5;
  const COLS = 3;

  const featuresMemory = memAllocFloat32Arr(ROWS * COLS);
  const features = new Float32Array(featuresMemory.buf, featuresMemory.off, ROWS * COLS);

  const labelsMemory = memAllocFloat32Arr(ROWS);
  const labels = new Float32Array(labelsMemory.buf, labelsMemory.off, ROWS);

  const predictionsMemory = memAllocFloat32Arr(ROWS);
  const predictions = new Float32Array(predictionsMemory.buf, predictionsMemory.off, ROWS);

  const MISSING_VALUE = -1;

  const modelSizeMemory = memAllocInt32Arr(1);
  const modelSizePtr = new Int32Array(modelSizeMemory.buf, modelSizeMemory.off, 1);

  const reserved = 1000000;

  const modelBytesMemory = memAllocInt32Arr(reserved);
  const modelBytes = new Int32Array(modelBytesMemory.buf, modelBytesMemory.off, reserved);

  for (let i = 0; i < ROWS; ++i) {
    let sum = 0;

    for (let j = 0; j < COLS; ++j) {
      features[COLS * i + j] = Math.random() * 10;
      sum += features[COLS * i + j];
    }

    labels[i] = Math.round(sum) % 3;

    console.log(`${i}   <-->   ${labels[i]}`);
  }

  console.log(features);
  console.log(labels);
  console.log(modelSizePtr);
  console.log(modelBytes);

  console.log('Start training:');

  const trainRes = train(features, ROWS, COLS, MISSING_VALUE, labels, ROWS,
    iterations, eta, maxDepth, lambda, alpha,
    modelSizePtr,
    modelBytes, reserved);

  console.log(`Training: ${trainRes}`);

  const predictRes = predict(features, ROWS, COLS, MISSING_VALUE, modelBytes, modelSizePtr[0], predictions, ROWS);

  console.log(`Prediciting: ${predictRes}`);

  console.log('Predictions:');
  for (let i = 0; i < ROWS; i++)
    console.log(Math.round(predictions[i]));

  memFree(featuresMemory.off);
  memFree(labelsMemory.off);
  memFree(predictionsMemory.off);
  memFree(modelSizeMemory.off);
  memFree(modelBytesMemory.off);

  console.log('Done!');
}

export function testXGBoost() {
  const iterations = 20;
  const eta = 0.3;
  const maxDepth = 6;
  const lambda = 1;
  const alpha = 0;

  const samplesCount = 5;
  const featuresCount = 3;
  const reserved = 1000000;
  const MISSING_VALUE = -1;

  // 1. TRAINING

  const trainMemory = allocTrainMemory(samplesCount, featuresCount, reserved, 1);

  const trainFeatures = new Float32Array(
    trainMemory.float32Buffer,
    trainMemory.featuresOffset,
    samplesCount * featuresCount,
  );

  const trainLabels = new Float32Array(
    trainMemory.float32Buffer,
    trainMemory.labelsOffset,
    samplesCount,
  );

  const utils = new Int32Array(
    trainMemory.int32Buffer,
    trainMemory.utilsOffset,
    samplesCount * featuresCount,
  );

  const packedModelBytes = new Int32Array(
    trainMemory.int32Buffer,
    trainMemory.modelOffset,
    reserved,
  );

  for (let i = 0; i < samplesCount; ++i) {
    let sum = 0;

    for (let j = 0; j < featuresCount; ++j) {
      trainFeatures[featuresCount * i + j] = Math.random() * 10;
      sum += trainFeatures[featuresCount * i + j];
    }

    trainLabels[i] = Math.round(sum) % 3;

    console.log(`${i}   <-->   ${trainLabels[i]}`);
  }

  console.log('Start training:');

  const trainRes = train(trainFeatures, samplesCount, featuresCount, MISSING_VALUE, trainLabels, samplesCount,
    iterations, eta, maxDepth, lambda, alpha,
    utils,
    packedModelBytes, reserved);

  console.log(`Training result: ${trainRes}`);
  console.log(`Model size: ${utils[0]}`);

  // 2. PREDICTING
  const predictMemory = allocPredictMemory(samplesCount, featuresCount, utils[0]);

  const predictFeatures = new Float32Array(
    predictMemory.float32Buffer,
    predictMemory.featuresOffset,
    samplesCount * featuresCount,
  );

  for (let i = 0; i < featuresCount * samplesCount; ++i)
    predictFeatures[i] = trainFeatures[i];

  const predictLabels = new Float32Array(
    predictMemory.float32Buffer,
    predictMemory.predictOffset,
    samplesCount,
  );

  const unpackedModelBytes = new Int32Array(
    predictMemory.int32Buffer,
    predictMemory.modelOffset,
    utils[0],
  );

  for (let i = 0; i < utils[0]; ++i)
    unpackedModelBytes[i] = packedModelBytes[i];

  console.log('Start predicting:');

  const predictRes = predict(
    predictFeatures,
    samplesCount,
    featuresCount,
    MISSING_VALUE,
    unpackedModelBytes,
    utils[0],
    predictLabels,
    samplesCount,
  );

  console.log(`Predicting result: ${predictRes}`);

  console.log('Predictions:');
  for (let i = 0; i < samplesCount; i++)
    console.log(Math.round(predictLabels[i]));

  freeTrainMemory(trainMemory);
  freePredictMemory(predictMemory);

  console.log('Done!');
}

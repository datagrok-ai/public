// XGBooster modeling tools

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {train, predict, allocTrainMemory, freeTrainMemory, allocPredictMemory,
  freePredictMemory} from '../wasm/xgbooster';

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

  const trainRes = train(trainFeatures, samplesCount, featuresCount, MISSING_VALUE, trainLabels,
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
  );

  console.log(`Predicting result: ${predictRes}`);

  console.log('Predictions:');
  for (let i = 0; i < samplesCount; i++)
    console.log(Math.round(predictLabels[i]));

  freeTrainMemory(trainMemory);
  freePredictMemory(predictMemory);

  console.log('Done!');
}

/** Interactivity tresholds */
enum INTERACTIVITY {
  MAX_SAMLPES = 1000,
  MAX_FEATURES = 10,
};

/** Reserve sizes */
enum RESERVED {
  MODEL = 1000000,
  UTILS = 1,
};

/** XGBoost specific constants */
const MISSING_VALUE = DG.FLOAT_NULL;
const SIZE_IDX = 0;

/** Data structs for training */
type TrainStructs = {
  features: Float32Array,
  target: Float32Array,
  utils: Int32Array,
  model: Int32Array,
};

/** Data structs for predicting */
type PredictStructs = {
  features: Float32Array,
  predict: Float32Array,
  model: Int32Array,
};

/** XGBoost modeling */
export class XGBooster {
  /** Check applicability */
  static isApplicable(features: DG.ColumnList, predictColumn: DG.Column): boolean {
    for (const col of features) {
      if (!col.matches('numerical'))
        return false;
    }
    if (!predictColumn.matches('numerical'))
      return false;

    return true;
  }

  /** Check interactivity */
  static isInteractive(features: DG.ColumnList, predictColumn: DG.Column): boolean {
    return (features.length <= INTERACTIVITY.MAX_FEATURES) &&
      (predictColumn.length <= INTERACTIVITY.MAX_SAMLPES);
  }

  private modelParams: Int32Array | undefined = undefined;

  constructor(packedModel?: Uint8Array) {
    if (packedModel) {
      try {
        // unpacking the model

      } catch (error) {
        throw new Error(`Failed to load model: ${(error instanceof Error ? error.message : 'the platform issue')}`);
      }
    }
  }

  /** Fit model */
  public fit(features: DG.ColumnList, target: DG.Column, iterations: number, eta: number,
    maxDepth: number, lambda: number, alpha: number) {
    // Allocate memory in wasm-buffer & put training data there
    const wasmStructs = this.getWasmTrainStructs(features, target);

    // Train model
    train(
      wasmStructs.features, target.length, features.length, MISSING_VALUE, // features data
      wasmStructs.target, // target/labels data
      iterations, eta, maxDepth, lambda, alpha, // hyperparameters of training
      wasmStructs.utils, // model utils
      wasmStructs.model, RESERVED.MODEL, // trained model
    );

    // Extract model params from wasm-buffer
    this.modelParams = this.getModelParams(wasmStructs);
  }

  /** Predict using trained model */
  public predict(features: DG.ColumnList): DG.Column | undefined {
    // Allocate memory in wasm-buffer & put training data there
    const wasmStructs = this.getWasmPredictStructs(features);

    const samplesCount = features.byIndex(0).length;
    const featuresCount = features.length;

    // Train model
    predict(
      wasmStructs.features, samplesCount, featuresCount, MISSING_VALUE,
      wasmStructs.model, this.modelParams?.length,
      wasmStructs.predict,
    );

    // Extract prediction column from wasm-buffer
    return this.getPredictCol(wasmStructs);
  }

  /** Return packed model */
  //public toBytes(): Uint8Array {}

  /** Allocate structs for training at the wasm-side */
  private allocTrainStructs(samplesCount: number, featuresCount: number): TrainStructs {
    const memory = allocTrainMemory(samplesCount, featuresCount, RESERVED.MODEL, RESERVED.UTILS);

    return {
      features: new Float32Array(memory.float32Buffer, memory.featuresOffset, samplesCount * featuresCount),
      target: new Float32Array(memory.float32Buffer, memory.labelsOffset, samplesCount),
      model: new Int32Array(memory.int32Buffer, memory.modelOffset, RESERVED.MODEL),
      utils: new Int32Array(memory.int32Buffer, memory.utilsOffset, RESERVED.UTILS),
    };
  }

  /** Get wasm-structures for training model */
  private getWasmTrainStructs(features: DG.ColumnList, target: DG.Column): TrainStructs {
    const samplesCount = target.length;
    const featuresCount = features.length;

    // Allocate wasm-buffer
    const wasmStructs = this.allocTrainStructs(samplesCount, featuresCount);

    // Put targets to wasm-buffer
    const rawTarget = target.getRawData();
    const wasmTarget = wasmStructs.target;
    for (let i = 0; i < samplesCount; ++i)
      wasmTarget[i] = rawTarget[i];

    // Put features to wasm-buffer
    const wasmFeatures = wasmStructs.features;
    for (let j = 0; j < featuresCount; ++j) {
      const raw = features.byIndex(j).getRawData();

      for (let i = 0; i < samplesCount; ++i)
        wasmFeatures[j * samplesCount + i] = raw[i];
    }

    //console.log(wasmStructs);

    return wasmStructs;
  }

  /** Get model params from wasm-buffer */
  private getModelParams(wasmStructs: TrainStructs): Int32Array {
    const modelParamsCount = wasmStructs.utils[SIZE_IDX];
    const params = new Int32Array(modelParamsCount);
    const trained = wasmStructs.model;

    for (let i = 0; i < modelParamsCount; ++i)
      params[i] = trained[i];

    return params;
  }

  /** Allocate structs for predicting at the wasm-side */
  private allocPredictStructs(samplesCount: number, featuresCount: number): PredictStructs {
    if (this.modelParams === undefined)
      throw new Error('Non-trained model applying');

    const memory = allocPredictMemory(samplesCount, featuresCount, this.modelParams.length);

    return {
      features: new Float32Array(memory.float32Buffer, memory.featuresOffset, samplesCount * featuresCount),
      predict: new Float32Array(memory.float32Buffer, memory.predictOffset, samplesCount),
      model: new Int32Array(memory.int32Buffer, memory.modelOffset, this.modelParams.length),
    };
  }

  /** Get wasm-structures for predicting */
  private getWasmPredictStructs(features: DG.ColumnList): PredictStructs {
    if (this.modelParams === undefined)
      throw new Error('Non-trained model applying');

    const samplesCount = features.byIndex(0).length;
    const featuresCount = features.length;

    // Allocate wasm-buffer
    const wasmStructs = this.allocPredictStructs(samplesCount, featuresCount);

    // Put model params to wasm-buffer
    const wasmModel = wasmStructs.model;
    const model = this.modelParams;
    const paramsCount = model.length;

    for (let i = 0; i < paramsCount; ++i)
      wasmModel[i] = model[i];

    // Put features to wasm-buffer
    const wasmFeatures = wasmStructs.features;
    for (let j = 0; j < featuresCount; ++j) {
      const raw = features.byIndex(j).getRawData();

      for (let i = 0; i < samplesCount; ++i)
        wasmFeatures[j * samplesCount + i] = raw[i];
    }

    //console.log(wasmStructs);

    return wasmStructs;
  }

  /** Get perdiction colum */
  private getPredictCol(wasmStructs: PredictStructs): DG.Column {
    const wasmPredict = wasmStructs.predict;
    const size = wasmPredict.length;
    const raw = new Float32Array(size);

    for (let i = 0; i < size; ++i)
      raw[i] = wasmPredict[i];

    return DG.Column.fromFloat32Array('Predict', raw, size);
  }
}

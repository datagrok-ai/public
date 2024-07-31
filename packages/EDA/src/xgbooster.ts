// XGBooster modeling tools

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {train, predictOld, allocTrainMemory, freeTrainMemory, allocPredictMemory,
  freePredictMemory, fit, predict} from '../wasm/xgbooster';

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

  const predictRes = predictOld(
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
  SAMLPES_HIGH = 100000,
  SAMLPES_MID = 50000,
  SAMPLES_LOW = 10000,
  FEATURES_HIGH = 10,
  FEATURES_MID = 20,
  FEATURES_LOW = 100,
};

/** Reserve sizes */
enum RESERVED {
  MODEL = 10000000,
  UTILS = 1,
  PACK = 128,
  SIZE = 4,
};

/** XGBoost specific constants */
const MISSING_VALUE = DG.FLOAT_NULL;
const SIZE_IDX = 0;
const ALIGN_VAL = 4;
const BLOCK_SIZE = 64;

enum TITLES {
  PREDICT = 'Prediction',
  TYPE = 'Type',
  PARAMS = 'Params count',
  CATS = 'Categories',
  CATS_SIZE = 'Categories size',
}

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
    if (!predictColumn.matches('numerical') && !predictColumn.matches('string'))
      return false;

    return true;
  }

  /** Check interactivity */
  static isInteractive(features: DG.ColumnList, predictColumn: DG.Column): boolean {
    const featuresCount = features.length;
    const samplesCount = predictColumn.length;

    if (samplesCount <= INTERACTIVITY.SAMPLES_LOW)
      return featuresCount <= INTERACTIVITY.FEATURES_LOW;

    if (samplesCount <= INTERACTIVITY.SAMLPES_MID)
      return featuresCount <= INTERACTIVITY.FEATURES_MID;

    if (samplesCount <= INTERACTIVITY.SAMLPES_HIGH)
      return featuresCount <= INTERACTIVITY.FEATURES_HIGH;

    return false;
  }

  private modelParams: Int32Array | undefined = undefined;
  private targetType: string | undefined = undefined;
  private targetCategories: string[] | undefined = undefined;

  constructor(packedModel?: Uint8Array) {
    if (packedModel) {
      try {
        let offset = 0;

        // Unpack header size
        const headArr = new Uint32Array(packedModel.buffer, offset, 1);
        const headerBytesSize = headArr[0];
        offset += RESERVED.SIZE;

        // Unpack header
        const headerDf = DG.DataFrame.fromByteArray(new Uint8Array(packedModel.buffer, offset, headerBytesSize));
        offset += headerBytesSize;

        // Extract model specification
        this.targetType = headerDf.get(TITLES.TYPE, 0) as string;
        const modelParamsCount = headerDf.get(TITLES.PARAMS, 0) as number;
        const categoriesBytesSize = headerDf.get(TITLES.CATS_SIZE, 0) as number;

        // Unpack categories
        if (categoriesBytesSize > 0) {
          const categoriesDf = DG.DataFrame.fromByteArray(
            new Uint8Array(packedModel.buffer, offset, categoriesBytesSize),
          );

          this.targetCategories = categoriesDf.col(TITLES.CATS)?.toList();
        }
        offset += categoriesBytesSize;

        offset = Math.ceil(offset / ALIGN_VAL) * ALIGN_VAL;

        // Unpack model params
        this.modelParams = new Int32Array(packedModel.buffer, offset, modelParamsCount);
      } catch (error) {
        throw new Error(`Failed to load model: ${(error instanceof Error ? error.message : 'the platform issue')}`);
      }
    }
  }

  /** Fit model */
  public fit(features: DG.ColumnList, target: DG.Column, iterations: number, eta: number,
    maxDepth: number, lambda: number, alpha: number) {
    // Type of the target
    this.targetType = target.type;

    // Store categories of string target
    if (this.targetType === DG.COLUMN_TYPE.STRING)
      this.targetCategories = target.categories;

    /*
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
    console.log(this.modelParams);

    freeTrainMemory(wasmStructs);*/

    this.modelParams = fit(features, target, MISSING_VALUE,
      iterations, eta, maxDepth, lambda, alpha, RESERVED.MODEL, RESERVED.UTILS,
    );
  }

  /** Predict using trained model */
  public predict(features: DG.ColumnList): DG.Column {/*
    // Allocate memory in wasm-buffer & put training data there
    const wasmStructs = this.getWasmPredictStructs(features);

    const samplesCount = features.byIndex(0).length;
    const featuresCount = features.length;

    // Train model
    predictOld(
      wasmStructs.features, samplesCount, featuresCount, MISSING_VALUE,
      wasmStructs.model, this.modelParams?.length,
      wasmStructs.predict,
    );

    // Extract prediction column from wasm-buffer
    const prediction = this.getPredictColOld(wasmStructs);

    console.log(wasmStructs.predict);

    freePredictMemory(wasmStructs);
    console.log(predict(features, MISSING_VALUE, this.modelParams));*/
    return this.getPredictCol(predict(features, MISSING_VALUE, this.modelParams));
  }

  /** Return packed model */
  public toBytes(): Uint8Array {
    if ((this.modelParams === undefined) || (this.targetType === undefined))
      throw new Error('Failed to pack non-trained model');

    // Categories bytes
    const categoriesBytes = (this.targetCategories !== undefined) ?
      DG.DataFrame.fromColumns([DG.Column.fromStrings(TITLES.CATS, this.targetCategories)]).toByteArray():
      undefined;

    const categoriesBytesSize = (categoriesBytes !== undefined) ? categoriesBytes.length : 0;

    const modelParamsBytesSize = this.modelParams.length * this.modelParams.BYTES_PER_ELEMENT;

    // Header with model specification
    const headerDf = DG.DataFrame.fromColumns([
      DG.Column.fromStrings(TITLES.TYPE, [this.targetType]),
      DG.Column.fromInt32Array(TITLES.PARAMS, new Int32Array([this.modelParams.length])),
      DG.Column.fromInt32Array(TITLES.CATS_SIZE, new Int32Array([categoriesBytesSize])),
    ]);

    // Header bytes
    const headerBytes = headerDf.toByteArray();
    const headerBytesSize = headerBytes.length;

    // Packed model
    const reservedSize = Math.ceil((RESERVED.SIZE +
      headerBytesSize + categoriesBytesSize + modelParamsBytesSize + RESERVED.PACK) / BLOCK_SIZE) * BLOCK_SIZE;

    const packedModel = new Uint8Array(reservedSize);

    let offset = 0;

    // Pack header size
    const headArr = new Uint32Array(packedModel.buffer, offset, 1);
    headArr[0] = headerBytesSize;
    offset += RESERVED.SIZE;

    // Pack header
    packedModel.set(headerBytes, offset);
    offset += headerBytesSize;

    // Pack categories
    if (categoriesBytesSize > 0)
      packedModel.set(categoriesBytes!, offset);
    offset += categoriesBytesSize;

    offset = Math.ceil(offset / ALIGN_VAL) * ALIGN_VAL;

    // Pack model params
    packedModel.set(new Uint8Array(this.modelParams.buffer), offset);

    return packedModel;
  }

  /** Allocate structs for training at the wasm-side */
  private allocTrainStructs(samplesCount: number, featuresCount: number): TrainStructs {
    const memory = allocTrainMemory(samplesCount, featuresCount, RESERVED.MODEL, RESERVED.UTILS);
    console.log(memory);

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

    console.log(`Params count: ${modelParamsCount}`);

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

  /** Get predicted string column */
  private getPredictedStringColOld(wasmStructs: PredictStructs): DG.Column {
    const wasmPredict = wasmStructs.predict;
    const samplesCount = wasmPredict.length;

    if (this.targetCategories === undefined)
      throw new Error('Predicting fails: undefined categories');

    const predClass = new Array<string>(samplesCount);

    for (let i = 0; i < samplesCount; ++i)
      predClass[i] = this.targetCategories[Math.round(wasmPredict[i])];

    return DG.Column.fromStrings(TITLES.PREDICT, predClass);
  }

  /** Get predicted int column */
  private getPredictedIntColOld(wasmStructs: PredictStructs): DG.Column {
    const wasmPredict = wasmStructs.predict;
    const samplesCount = wasmPredict.length;

    const rawInts = new Int32Array(samplesCount);

    for (let i = 0; i < samplesCount; ++i)
      rawInts[i] = Math.round(wasmPredict[i]);

    return DG.Column.fromInt32Array(TITLES.PREDICT, rawInts, samplesCount);
  }

  /** Get predicted int column */
  private getPredictedFloatColOld(wasmStructs: PredictStructs): DG.Column {
    const wasmPredict = wasmStructs.predict;
    const samplesCount = wasmPredict.length;

    const rawFloats = new Float32Array(samplesCount);

    for (let i = 0; i < samplesCount; ++i)
      rawFloats[i] = wasmPredict[i];

    return DG.Column.fromFloat32Array(TITLES.PREDICT, rawFloats, samplesCount);
  }

  /** Get predicted bigint column */
  private getPredictedBigIntColOld(wasmStructs: PredictStructs): DG.Column {
    const wasmPredict = wasmStructs.predict;
    const samplesCount = wasmPredict.length;

    const rawInts = new BigInt64Array(samplesCount);

    for (let i = 0; i < samplesCount; ++i)
      rawInts[i] = BigInt(Math.round(wasmPredict[i]));

    return DG.Column.fromBigInt64Array(TITLES.PREDICT, rawInts);
  }

  /** Get perdiction colum */
  private getPredictColOld(wasmStructs: PredictStructs): DG.Column {
    switch (this.targetType) {
    case DG.COLUMN_TYPE.STRING:
      return this.getPredictedStringColOld(wasmStructs);

    case DG.COLUMN_TYPE.INT:
      return this.getPredictedIntColOld(wasmStructs);

    case DG.COLUMN_TYPE.BIG_INT:
      return this.getPredictedBigIntColOld(wasmStructs);

    default:
      return this.getPredictedFloatColOld(wasmStructs);
    }
  }

  /** Get predicted string column */
  private getPredictedStringCol(prediction: Float32Array): DG.Column {
    const samplesCount = prediction.length;

    if (this.targetCategories === undefined)
      throw new Error('Predicting fails: undefined categories');

    const predClass = new Array<string>(samplesCount);

    for (let i = 0; i < samplesCount; ++i)
      predClass[i] = this.targetCategories[Math.round(prediction[i])];

    return DG.Column.fromStrings(TITLES.PREDICT, predClass);
  }

  /** Get predicted int column */
  private getPredictedIntCol(prediction: Float32Array): DG.Column {
    const samplesCount = prediction.length;

    const rawInts = new Int32Array(samplesCount);

    for (let i = 0; i < samplesCount; ++i)
      rawInts[i] = Math.round(prediction[i]);

    return DG.Column.fromInt32Array(TITLES.PREDICT, rawInts, samplesCount);
  }

  /** Get predicted bigint column */
  private getPredictedBigIntCol(prediction: Float32Array): DG.Column {
    const samplesCount = prediction.length;

    const rawInts = new BigInt64Array(samplesCount);

    for (let i = 0; i < samplesCount; ++i)
      rawInts[i] = BigInt(Math.round(prediction[i]));

    return DG.Column.fromBigInt64Array(TITLES.PREDICT, rawInts);
  }

  /** Get perdiction colum */
  private getPredictCol(prediction: Float32Array): DG.Column {
    switch (this.targetType) {
    case DG.COLUMN_TYPE.STRING:
      return this.getPredictedStringCol(prediction);

    case DG.COLUMN_TYPE.INT:
      return this.getPredictedIntCol(prediction);

    case DG.COLUMN_TYPE.BIG_INT:
      return this.getPredictedBigIntCol(prediction);

    default:
      return DG.Column.fromFloat32Array(TITLES.PREDICT, prediction);
    }
  }
}

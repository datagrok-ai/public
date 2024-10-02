// XGBooster modeling tools

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {predict, fitInWebWorker} from '../wasm/xgbooster';

/** Default hyperparameters */
enum DEFAULT {
  ITERATIONS = 20,
  ETA = 0.3,
  MAX_DEPTH = 6,
  LAMBDA = 1,
  ALPHA = 0,
};

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
const ALIGN_VAL = 4;
const BLOCK_SIZE = 64;

enum TITLES {
  PREDICT = 'Prediction',
  TYPE = 'Type',
  PARAMS = 'Params count',
  CATS = 'Categories',
  CATS_SIZE = 'Categories size',
}

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
  public async fit(features: DG.ColumnList, target: DG.Column, iterations: number = DEFAULT.ITERATIONS,
    eta: number = DEFAULT.ETA, maxDepth: number = DEFAULT.MAX_DEPTH, lambda: number = DEFAULT.LAMBDA,
    alpha: number = DEFAULT.ALPHA) {
    // Type of the target
    this.targetType = target.type;

    // Store categories of string target
    if (this.targetType === DG.COLUMN_TYPE.STRING)
      this.targetCategories = target.categories;

    // Train model params
    this.modelParams = await fitInWebWorker(features, target, MISSING_VALUE,
      iterations, eta, maxDepth, lambda, alpha, RESERVED.MODEL, RESERVED.UTILS,
    );
  }

  /** Predict using trained model */
  public predict(features: DG.ColumnList): DG.Column {
    if (this.modelParams === undefined)
      throw new Error('Failed to apply non-trained model');

    // Get prediction
    const prediction = predict(features, MISSING_VALUE, this.modelParams);

    // Create an appropriate column
    switch (this.targetType) {
    case DG.COLUMN_TYPE.STRING:
      return this.stringColPrediction(prediction);

    case DG.COLUMN_TYPE.INT:
      return this.intColPrediction(prediction);

    case DG.COLUMN_TYPE.BIG_INT:
      return this.bigIntColPrediction(prediction);

    default:
      return DG.Column.fromFloat32Array(TITLES.PREDICT, prediction);
    }
  }

  /** Return packed model */
  public toBytes(): Uint8Array {
    if ((this.modelParams === undefined) || (this.targetType === undefined))
      throw new Error('Failed to pack non-trained model');

    // Categories bytes
    const categoriesBytes = (this.targetCategories !== undefined) ? DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, TITLES.CATS, this.targetCategories),
    ]).toByteArray(): undefined;

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
  } // toBytes

  /** Return predicted string column */
  private stringColPrediction(prediction: Float32Array): DG.Column {
    const samplesCount = prediction.length;

    if (this.targetCategories === undefined)
      throw new Error('Predicting fails: undefined categories');

    const predClass = new Array<string>(samplesCount);

    const maxCategory = this.targetCategories.length - 1;
    const categoryIdx = (val: number) => Math.max(0, Math.min(val, maxCategory));

    for (let i = 0; i < samplesCount; ++i)
      predClass[i] = this.targetCategories[categoryIdx(Math.round(prediction[i]))];

    return DG.Column.fromList(DG.COLUMN_TYPE.STRING, TITLES.PREDICT, predClass);
  }

  /** Return predicted int column */
  private intColPrediction(prediction: Float32Array): DG.Column {
    const samplesCount = prediction.length;

    const rawInts = new Int32Array(samplesCount);

    for (let i = 0; i < samplesCount; ++i)
      rawInts[i] = Math.round(prediction[i]);

    return DG.Column.fromInt32Array(TITLES.PREDICT, rawInts, samplesCount);
  }

  /** Return predicted bigint column */
  private bigIntColPrediction(prediction: Float32Array): DG.Column {
    const samplesCount = prediction.length;

    const rawInts = new BigInt64Array(samplesCount);

    for (let i = 0; i < samplesCount; ++i)
      rawInts[i] = BigInt(Math.round(prediction[i]));

    return DG.Column.fromBigInt64Array(TITLES.PREDICT, rawInts);
  }
}

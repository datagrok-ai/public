// XGBooster modeling tools

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {ColumnView, XgbHyperParams, XgbObjective,
  fitXgb, freeXgbModel, initXgboost, loadXgbModel, predictXgb} from '../wasm/xgbooster';

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

/** XGBoost specific constants */
const MISSING_VALUE = DG.FLOAT_NULL;
const ALIGN_VAL = 4;
const BLOCK_SIZE = 64;
const SIZE_BYTES = 4;
const PACK_RESERVE = 128;

/** Packed model container version. Models saved by pre-1.7.0 EDA (v1, no
 *  Version field) are NOT supported and must be retrained. */
const CONTAINER_VERSION = 2;

enum TITLES {
  PREDICT = 'Prediction',
  TYPE = 'Type',
  MODEL_SIZE = 'Model size',
  OBJECTIVE = 'Objective',
  VERSION = 'Version',
  CATS = 'Categories',
  CATS_SIZE = 'Categories size',
  WAS_BOOL = 'Was bool',
}

/** String labels of a boolean target after conversion to string */
const BOOL_TRUE = 'true';

/** Probability threshold for the binary objective */
const BINARY_THRESHOLD = 0.5;

/** Column views sliced to the column length: col.getRawData() may be LONGER
 *  than col.length (capacity mechanism) - never copy it unsliced. */
function columnView(col: DG.Column): ColumnView {
  return (col.getRawData() as ColumnView).subarray(0, col.length) as ColumnView;
}

function featureViews(features: DG.ColumnList): ColumnView[] {
  const views: ColumnView[] = [];
  for (const col of features)
    views.push(columnView(col));
  return views;
}

/** XGBoost modeling */
export class XGBooster {
  /** Check applicability */
  static isApplicable(features: DG.ColumnList, predictColumn: DG.Column): boolean {
    for (const col of features) {
      if (!col.matches('numerical'))
        return false;
    }
    if (!predictColumn.matches('numerical') && !predictColumn.matches('string') &&
      (predictColumn.type !== DG.COLUMN_TYPE.BOOL))
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

  private modelBytes: Uint8Array | undefined = undefined;
  private objective: XgbObjective = XgbObjective.Regression;
  private targetType: string | undefined = undefined;
  private targetCategories: string[] | undefined = undefined;
  /** True if the original target was a boolean column */
  private targetWasBool = false;
  /** Cached live wasm handle; recreated lazily from modelBytes. */
  private handle = 0;

  constructor(packedModel?: Uint8Array) {
    if (packedModel) {
      try {
        let offset = 0;

        // Unpack header size
        const headArr = new Uint32Array(packedModel.buffer, packedModel.byteOffset + offset, 1);
        const headerBytesSize = headArr[0];
        offset += SIZE_BYTES;

        // Unpack header
        const headerDf = DG.DataFrame.fromByteArray(
          new Uint8Array(packedModel.buffer, packedModel.byteOffset + offset, headerBytesSize));
        offset += headerBytesSize;

        // Extract model specification
        this.targetType = headerDf.get(TITLES.TYPE, 0) as string;
        const categoriesBytesSize = headerDf.get(TITLES.CATS_SIZE, 0) as number;

        const wasBoolCol = headerDf.col(TITLES.WAS_BOOL);
        this.targetWasBool = (wasBoolCol !== null) && (wasBoolCol.get(0) === 1);

        // Models saved before the Version field (pre-1.7.0) are unsupported.
        if (headerDf.col(TITLES.VERSION) === null) {
          throw new Error('this model was saved by an older EDA version ' +
            'and is no longer supported - retrain the model');
        }
        const modelBytesSize = headerDf.get(TITLES.MODEL_SIZE, 0) as number;
        this.objective = headerDf.get(TITLES.OBJECTIVE, 0) as number;

        // Unpack categories
        if (categoriesBytesSize > 0) {
          const categoriesDf = DG.DataFrame.fromByteArray(
            new Uint8Array(packedModel.buffer, packedModel.byteOffset + offset, categoriesBytesSize),
          );
          this.targetCategories = categoriesDf.col(TITLES.CATS)?.toList();
        }
        offset += categoriesBytesSize;

        offset = Math.ceil(offset / ALIGN_VAL) * ALIGN_VAL;

        // Unpack model bytes (copy: the container buffer may be reused by the caller)
        this.modelBytes = new Uint8Array(packedModel.buffer,
          packedModel.byteOffset + offset, modelBytesSize).slice();
      } catch (error) {
        throw new Error(`Failed to load model: ${(error instanceof Error ? error.message : 'the platform issue')}`);
      }
    }
  }

  /** Fit model */
  public async fit(features: DG.ColumnList, target: DG.Column, iterations: number = DEFAULT.ITERATIONS,
    eta: number = DEFAULT.ETA, maxDepth: number = DEFAULT.MAX_DEPTH, lambda: number = DEFAULT.LAMBDA,
    alpha: number = DEFAULT.ALPHA) {
    // The wasm build has C++ exceptions disabled: a bad call aborts the whole
    // instance. Validate everything here, before any wasm call.
    this.validateFitInputs(features, target, iterations, eta, maxDepth, lambda, alpha);

    // Ensure the main-thread module (used by the synchronous predict) is up.
    // Idempotent; needed where package init does not run (package-test bundle).
    await initXgboost();

    // A boolean target is converted to a 2-class string target ('true'/'false')
    this.targetWasBool = (target.type === DG.COLUMN_TYPE.BOOL);
    const tgt = this.targetWasBool ? target.convertTo(DG.COLUMN_TYPE.STRING) : target;

    this.targetType = tgt.type;

    // Objective by target type: string/bool targets are classified honestly
    // (logistic / softmax), numeric targets are regressed.
    let numClass = 0;
    if (this.targetType === DG.COLUMN_TYPE.STRING) {
      this.targetCategories = tgt.categories;
      numClass = this.targetCategories.length;
      if (numClass < 2)
        throw new Error('XGBoost: classification target must have at least 2 categories');
      this.objective = (numClass === 2) ? XgbObjective.Binary : XgbObjective.Multiclass;
    } else {
      this.targetCategories = undefined;
      this.objective = XgbObjective.Regression;
    }

    const hyper: XgbHyperParams = {iterations, eta, maxDepth, lambda, alpha};
    // String labels are category codes (0..K-1 by construction).
    this.modelBytes = await fitXgb(featureViews(features), columnView(tgt),
      target.length, MISSING_VALUE, this.objective,
      this.objective === XgbObjective.Multiclass ? numClass : 0, hyper);

    this.invalidateHandle();
  }

  /** Predict using trained model */
  public predict(features: DG.ColumnList): DG.Column {
    if (this.modelBytes === undefined)
      throw new Error('Failed to apply non-trained model');

    const samplesCount = features.byIndex(0).length;
    let prediction: Float32Array;
    try {
      if (this.handle === 0)
        this.handle = loadXgbModel(this.modelBytes);
      prediction = predictXgb(this.handle, featureViews(features), samplesCount, MISSING_VALUE);
    } catch (err) {
      // On a module crash the cached handle is dead; drop it so that the
      // next call reloads the model into the re-initialized module.
      this.invalidateHandle();
      throw err;
    }

    switch (this.objective) {
    case XgbObjective.Binary:
      return this.binaryPrediction(prediction);

    case XgbObjective.Multiclass:
      return this.multiclassPrediction(prediction);

    default:
      return this.regressionPrediction(prediction);
    }
  }

  /** Release the cached wasm model handle (model bytes are kept). */
  public dispose(): void {
    this.invalidateHandle();
  }

  /** Return packed model (container v2) */
  public toBytes(): Uint8Array {
    if ((this.modelBytes === undefined) || (this.targetType === undefined))
      throw new Error('Failed to pack non-trained model');

    // Categories bytes
    const categoriesBytes = (this.targetCategories !== undefined) ? DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, TITLES.CATS, this.targetCategories),
    ]).toByteArray() : undefined;

    const categoriesBytesSize = (categoriesBytes !== undefined) ? categoriesBytes.length : 0;

    // Header with model specification
    const headerDf = DG.DataFrame.fromColumns([
      DG.Column.fromStrings(TITLES.TYPE, [this.targetType]),
      DG.Column.fromInt32Array(TITLES.MODEL_SIZE, new Int32Array([this.modelBytes.length])),
      DG.Column.fromInt32Array(TITLES.CATS_SIZE, new Int32Array([categoriesBytesSize])),
      DG.Column.fromInt32Array(TITLES.WAS_BOOL, new Int32Array([this.targetWasBool ? 1 : 0])),
      DG.Column.fromInt32Array(TITLES.OBJECTIVE, new Int32Array([this.objective])),
      DG.Column.fromInt32Array(TITLES.VERSION, new Int32Array([CONTAINER_VERSION])),
    ]);

    const headerBytes = headerDf.toByteArray();
    const headerBytesSize = headerBytes.length;

    // Packed model
    const reservedSize = Math.ceil((SIZE_BYTES +
      headerBytesSize + categoriesBytesSize + this.modelBytes.length + PACK_RESERVE) / BLOCK_SIZE) * BLOCK_SIZE;

    const packedModel = new Uint8Array(reservedSize);

    let offset = 0;

    // Pack header size
    const headArr = new Uint32Array(packedModel.buffer, offset, 1);
    headArr[0] = headerBytesSize;
    offset += SIZE_BYTES;

    // Pack header
    packedModel.set(headerBytes, offset);
    offset += headerBytesSize;

    // Pack categories
    if (categoriesBytesSize > 0)
      packedModel.set(categoriesBytes!, offset);
    offset += categoriesBytesSize;

    offset = Math.ceil(offset / ALIGN_VAL) * ALIGN_VAL;

    // Pack model bytes
    packedModel.set(this.modelBytes, offset);

    return packedModel;
  } // toBytes

  private invalidateHandle(): void {
    if (this.handle !== 0) {
      freeXgbModel(this.handle);
      this.handle = 0;
    }
  }

  private validateFitInputs(features: DG.ColumnList, target: DG.Column,
    iterations: number, eta: number, maxDepth: number, lambda: number, alpha: number): void {
    if (features.length < 1)
      throw new Error('XGBoost: no feature columns');
    if (target.length < 1)
      throw new Error('XGBoost: empty target column');
    for (const col of features) {
      if (col.length !== target.length)
        throw new Error(`XGBoost: column "${col.name}" length differs from the target length`);
    }
    if (target.stats.missingValueCount > 0)
      throw new Error('XGBoost: target column contains missing values');
    if (!Number.isInteger(iterations) || iterations < 1)
      throw new Error('XGBoost: iterations must be a positive integer');
    if (!(eta > 0 && eta <= 1))
      throw new Error('XGBoost: eta must be in (0, 1]');
    if (!Number.isInteger(maxDepth) || maxDepth < 1)
      throw new Error('XGBoost: max depth must be a positive integer');
    if (!(lambda >= 0))
      throw new Error('XGBoost: lambda must be non-negative');
    if (!(alpha >= 0))
      throw new Error('XGBoost: alpha must be non-negative');
  }

  // ------------------------------------------------- prediction decoding --

  private categoryByIndex(idx: number): string {
    if (this.targetCategories === undefined)
      throw new Error('Predicting fails: undefined categories');
    const maxCategory = this.targetCategories.length - 1;
    return this.targetCategories[Math.max(0, Math.min(idx, maxCategory))];
  }

  /** binary:logistic output: probability of class 1 */
  private binaryPrediction(prediction: Float32Array): DG.Column {
    const samplesCount = prediction.length;

    if (this.targetWasBool) {
      const predClass = new Array<boolean>(samplesCount);
      for (let i = 0; i < samplesCount; ++i)
        predClass[i] = this.categoryByIndex(prediction[i] >= BINARY_THRESHOLD ? 1 : 0) === BOOL_TRUE;
      return DG.Column.fromList(DG.COLUMN_TYPE.BOOL, TITLES.PREDICT, predClass);
    }

    const predClass = new Array<string>(samplesCount);
    for (let i = 0; i < samplesCount; ++i)
      predClass[i] = this.categoryByIndex(prediction[i] >= BINARY_THRESHOLD ? 1 : 0);
    return DG.Column.fromList(DG.COLUMN_TYPE.STRING, TITLES.PREDICT, predClass);
  }

  /** multi:softmax output: class index */
  private multiclassPrediction(prediction: Float32Array): DG.Column {
    const samplesCount = prediction.length;
    const predClass = new Array<string>(samplesCount);
    for (let i = 0; i < samplesCount; ++i)
      predClass[i] = this.categoryByIndex(Math.round(prediction[i]));
    return DG.Column.fromList(DG.COLUMN_TYPE.STRING, TITLES.PREDICT, predClass);
  }

  private regressionPrediction(prediction: Float32Array): DG.Column {
    switch (this.targetType) {
    case DG.COLUMN_TYPE.INT:
      return this.intColPrediction(prediction);
    case DG.COLUMN_TYPE.BIG_INT:
      return this.bigIntColPrediction(prediction);
    default:
      return DG.Column.fromFloat32Array(TITLES.PREDICT, prediction);
    }
  }

  private intColPrediction(prediction: Float32Array): DG.Column {
    const samplesCount = prediction.length;
    const rawInts = new Int32Array(samplesCount);
    for (let i = 0; i < samplesCount; ++i)
      rawInts[i] = Math.round(prediction[i]);
    return DG.Column.fromInt32Array(TITLES.PREDICT, rawInts, samplesCount);
  }

  private bigIntColPrediction(prediction: Float32Array): DG.Column {
    const samplesCount = prediction.length;
    const rawInts = new BigInt64Array(samplesCount);
    for (let i = 0; i < samplesCount; ++i)
      rawInts[i] = BigInt(Math.round(prediction[i]));
    return DG.Column.fromBigInt64Array(TITLES.PREDICT, rawInts);
  }
}

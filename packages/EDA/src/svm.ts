// SVM (LIBSVM) modeling tools.
//
// Structure mirrors src/xgbooster.ts (packed container, cached handle,
// prediction decoding); feature standardization mirrors src/softmax-classifier.ts
// (the kernel never normalizes - "variant B"). Two tasks are chosen from the
// target column - string/bool gives classification (C_SVC, one-vs-one for 3+
// classes), numeric gives regression (EPSILON_SVR).

import * as DG from 'datagrok-api/dg';

import {ColumnView, KernelType, SvmHyperParams, SvmType,
  fitSvm, freeSvmModel, initSvm, loadSvmModel, predictSvm} from '../wasm/svm';

/** Fixed training hyperparameters not exposed in the UI. */
enum DEFAULT {
  EPS = 0.001, // stopping tolerance
  CACHE_MB = 100, // kernel cache size
  NU = 0.5, // for nu-variants (unused in v1)
}

/** Interactivity thresholds. SMO is superlinear in the sample count
 *  (RBF ~ O(n^2)-O(n^3)), so these are far below the XGBoost thresholds. */
enum INTERACTIVITY {
  MAX_SAMPLES = 5000,
  MAX_FEATURES = 50,
}

const PRED_NAME = 'Prediction';
const BOOL_TRUE = 'true';

/** Packed container layout constants (same scheme as xgbooster.ts). */
const ALIGN_VAL = 4;
const BLOCK_SIZE = 64;
const SIZE_BYTES = 4;
const PACK_RESERVE = 128;
const CONTAINER_VERSION = 1;

enum TITLES {
  TYPE = 'Type',
  MODEL_SIZE = 'Model size',
  CATS = 'Categories',
  CATS_SIZE = 'Categories size',
  STATS_SIZE = 'Stats size',
  WAS_BOOL = 'Was bool',
  VERSION = 'Version',
  AVGS = 'Avg-s',
  STDEVS = 'Stddev-s',
}

/** User-facing subset of hyperparameters; the rest are filled by fit(). */
export interface SvmFitOptions {
  kernelType: KernelType;
  C: number;
  gamma: number; // 0 = 1/features
  degree: number;
  coef0: number;
  epsilon: number; // SVR loss epsilon
}

/** Column view sliced to the column length: col.getRawData() may be LONGER
 *  than col.length (capacity mechanism) - never copy it unsliced. */
function columnView(col: DG.Column): ColumnView {
  return (col.getRawData() as ColumnView).subarray(0, col.length) as ColumnView;
}

/** SVM modeling */
export class SVM {
  /** Check applicability */
  static isApplicable(features: DG.ColumnList, predictColumn: DG.Column): boolean {
    for (const col of features) {
      if (!col.matches('numerical'))
        return false;
    }
    return predictColumn.matches('numerical') || predictColumn.matches('string') ||
      (predictColumn.type === DG.COLUMN_TYPE.BOOL);
  }

  /** Check interactivity */
  static isInteractive(features: DG.ColumnList, predictColumn: DG.Column): boolean {
    return (predictColumn.length <= INTERACTIVITY.MAX_SAMPLES) &&
      (features.length <= INTERACTIVITY.MAX_FEATURES);
  }

  private modelBytes: Uint8Array | undefined = undefined;
  private targetType: string | undefined = undefined;
  private targetCategories: string[] | undefined = undefined;
  private targetWasBool = false;
  private avgs: Float32Array = new Float32Array(0);
  private stdevs: Float32Array = new Float32Array(0);
  private featuresCount = 0;
  /** Cached live wasm handle; recreated lazily from modelBytes. */
  private handle = 0;

  constructor(packedModel?: Uint8Array) {
    if (packedModel === undefined)
      return;

    try {
      let offset = 0;

      // Header size
      const headerBytesSize = new Uint32Array(packedModel.buffer, packedModel.byteOffset + offset, 1)[0];
      offset += SIZE_BYTES;

      // Header
      const headerDf = DG.DataFrame.fromByteArray(
        new Uint8Array(packedModel.buffer, packedModel.byteOffset + offset, headerBytesSize));
      offset += headerBytesSize;

      if (headerDf.col(TITLES.VERSION) === null)
        throw new Error('unsupported model container - retrain the model');

      this.targetType = headerDf.get(TITLES.TYPE, 0) as string;
      const categoriesBytesSize = headerDf.get(TITLES.CATS_SIZE, 0) as number;
      const statsBytesSize = headerDf.get(TITLES.STATS_SIZE, 0) as number;
      const modelBytesSize = headerDf.get(TITLES.MODEL_SIZE, 0) as number;
      const wasBoolCol = headerDf.col(TITLES.WAS_BOOL);
      this.targetWasBool = (wasBoolCol !== null) && (wasBoolCol.get(0) === 1);

      // Categories
      if (categoriesBytesSize > 0) {
        const categoriesDf = DG.DataFrame.fromByteArray(
          new Uint8Array(packedModel.buffer, packedModel.byteOffset + offset, categoriesBytesSize));
        this.targetCategories = categoriesDf.col(TITLES.CATS)?.toList();
      }
      offset += categoriesBytesSize;

      // Standardization stats
      const statsDf = DG.DataFrame.fromByteArray(
        new Uint8Array(packedModel.buffer, packedModel.byteOffset + offset, statsBytesSize));
      this.avgs = statsDf.col(TITLES.AVGS)!.getRawData() as Float32Array;
      this.stdevs = statsDf.col(TITLES.STDEVS)!.getRawData() as Float32Array;
      this.featuresCount = this.avgs.length;
      offset += statsBytesSize;

      offset = Math.ceil(offset / ALIGN_VAL) * ALIGN_VAL;

      // Model bytes (copy: the container buffer may be reused by the caller)
      this.modelBytes = new Uint8Array(packedModel.buffer,
        packedModel.byteOffset + offset, modelBytesSize).slice();
    } catch (error) {
      throw new Error(`Failed to load model: ${(error instanceof Error ? error.message : 'the platform issue')}`);
    }
  }

  /** Fit model */
  public async fit(features: DG.ColumnList, target: DG.Column, options: SvmFitOptions): Promise<void> {
    this.validateFitInputs(features, target, options);

    // Ensure the main-thread module (used by the synchronous predict) is up.
    // Idempotent; needed where package init does not run (package-test bundle).
    await initSvm();

    // A boolean target is converted to a 2-class string target ('true'/'false').
    this.targetWasBool = (target.type === DG.COLUMN_TYPE.BOOL);
    const tgt = this.targetWasBool ? target.convertTo(DG.COLUMN_TYPE.STRING) : target;
    this.targetType = tgt.type;

    // Task from the target column - string/bool is classification, numeric is regression.
    let svmType: SvmType;
    let labels: ColumnView;
    if (tgt.type === DG.COLUMN_TYPE.STRING) {
      this.targetCategories = tgt.categories;
      if (this.targetCategories.length < 2)
        throw new Error('SVM: classification target must have at least 2 categories');
      svmType = SvmType.CSvc;
      labels = columnView(tgt); // Int32Array of category codes (0..K-1)
    } else {
      this.targetCategories = undefined;
      svmType = SvmType.EpsilonSvr;
      labels = this.regressionLabels(tgt);
    }

    this.featuresCount = features.length;
    this.extractStats(features);
    const cols = this.standardizedColumns(features);

    const hyper: SvmHyperParams = {
      svmType,
      kernelType: options.kernelType,
      C: options.C,
      gamma: options.gamma, // 0 -> 1/features on the wasm side
      degree: options.degree,
      coef0: options.coef0,
      epsilon: options.epsilon,
      nu: DEFAULT.NU,
      eps: DEFAULT.EPS,
      cacheSize: DEFAULT.CACHE_MB,
      shrinking: 1,
      probability: 0,
    };

    this.modelBytes = await fitSvm(cols, labels, target.length, hyper);
    this.invalidateHandle();
  }

  /** Predict using the trained model */
  public predict(features: DG.ColumnList): DG.Column {
    if (this.modelBytes === undefined)
      throw new Error('Failed to apply non-trained model');
    if (features.length !== this.featuresCount)
      throw new Error('SVM: prediction feature count differs from the trained model');

    const nRows = features.byIndex(0).length;
    const cols = this.standardizedColumns(features);

    let prediction: Float32Array;
    try {
      if (this.handle === 0)
        this.handle = loadSvmModel(this.modelBytes);
      prediction = predictSvm(this.handle, cols, nRows);
    } catch (err) {
      // On a module crash the cached handle is dead; drop it so the next call
      // reloads the model into the re-initialized module.
      this.invalidateHandle();
      throw err;
    }

    return (this.targetType === DG.COLUMN_TYPE.STRING) ?
      this.classificationPrediction(prediction) :
      this.regressionPrediction(prediction);
  }

  /** Release the cached wasm model handle (model bytes are kept). */
  public dispose(): void {
    this.invalidateHandle();
  }

  /** Return packed model (container v1) */
  public toBytes(): Uint8Array {
    if ((this.modelBytes === undefined) || (this.targetType === undefined))
      throw new Error('Failed to pack non-trained model');

    // Categories
    const categoriesBytes = (this.targetCategories !== undefined) ? DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, TITLES.CATS, this.targetCategories),
    ]).toByteArray() : undefined;
    const categoriesBytesSize = (categoriesBytes !== undefined) ? categoriesBytes.length : 0;

    // Standardization stats
    const statsBytes = DG.DataFrame.fromColumns([
      DG.Column.fromFloat32Array(TITLES.AVGS, this.avgs),
      DG.Column.fromFloat32Array(TITLES.STDEVS, this.stdevs),
    ]).toByteArray();
    const statsBytesSize = statsBytes.length;

    // Header
    const headerBytes = DG.DataFrame.fromColumns([
      DG.Column.fromStrings(TITLES.TYPE, [this.targetType]),
      DG.Column.fromInt32Array(TITLES.MODEL_SIZE, new Int32Array([this.modelBytes.length])),
      DG.Column.fromInt32Array(TITLES.CATS_SIZE, new Int32Array([categoriesBytesSize])),
      DG.Column.fromInt32Array(TITLES.STATS_SIZE, new Int32Array([statsBytesSize])),
      DG.Column.fromInt32Array(TITLES.WAS_BOOL, new Int32Array([this.targetWasBool ? 1 : 0])),
      DG.Column.fromInt32Array(TITLES.VERSION, new Int32Array([CONTAINER_VERSION])),
    ]).toByteArray();
    const headerBytesSize = headerBytes.length;

    const reservedSize = Math.ceil((SIZE_BYTES + headerBytesSize + categoriesBytesSize +
      statsBytesSize + this.modelBytes.length + PACK_RESERVE) / BLOCK_SIZE) * BLOCK_SIZE;
    const packedModel = new Uint8Array(reservedSize);

    let offset = 0;

    new Uint32Array(packedModel.buffer, offset, 1)[0] = headerBytesSize;
    offset += SIZE_BYTES;

    packedModel.set(headerBytes, offset);
    offset += headerBytesSize;

    if (categoriesBytesSize > 0)
      packedModel.set(categoriesBytes!, offset);
    offset += categoriesBytesSize;

    packedModel.set(statsBytes, offset);
    offset += statsBytesSize;

    offset = Math.ceil(offset / ALIGN_VAL) * ALIGN_VAL;

    packedModel.set(this.modelBytes, offset);

    return packedModel;
  }

  // ---------------------------------------------------------------- internals --

  private invalidateHandle(): void {
    if (this.handle !== 0) {
      freeSvmModel(this.handle);
      this.handle = 0;
    }
  }

  private validateFitInputs(features: DG.ColumnList, target: DG.Column, options: SvmFitOptions): void {
    if (features.length < 1)
      throw new Error('SVM: no feature columns');
    if (target.length < 1)
      throw new Error('SVM: empty target column');
    for (const col of features) {
      if (col.length !== target.length)
        throw new Error(`SVM: column "${col.name}" length differs from the target length`);
      if (!col.matches('numerical'))
        throw new Error(`SVM: feature column "${col.name}" is not numerical`);
      if (col.stats.missingValueCount > 0)
        throw new Error(`SVM: feature column "${col.name}" contains missing values`);
    }
    if (target.stats.missingValueCount > 0)
      throw new Error('SVM: target column contains missing values');
    if (!(options.C > 0))
      throw new Error('SVM: cost C must be positive');
    if (!(options.gamma >= 0))
      throw new Error('SVM: gamma must be non-negative');
    if (!Number.isInteger(options.degree) || options.degree < 1)
      throw new Error('SVM: degree must be a positive integer');
    if (!(options.epsilon >= 0))
      throw new Error('SVM: epsilon must be non-negative');
  }

  private extractStats(features: DG.ColumnList): void {
    this.avgs = new Float32Array(this.featuresCount);
    this.stdevs = new Float32Array(this.featuresCount);
    let j = 0;
    for (const col of features) {
      this.avgs[j] = col.stats.avg;
      this.stdevs[j] = col.stats.stdev;
      ++j;
    }
  }

  /** One standardized float32 view per feature column (zero when stdev == 0). */
  private standardizedColumns(features: DG.ColumnList): Float32Array[] {
    const nRows = features.byIndex(0).length;
    const out = new Array<Float32Array>(features.length);
    let j = 0;
    for (const col of features) {
      const raw = columnView(col);
      const avg = this.avgs[j];
      const stdev = this.stdevs[j];
      const c = new Float32Array(nRows);
      if (stdev > 0) {
        for (let i = 0; i < nRows; ++i)
          c[i] = (raw[i] - avg) / stdev;
      }
      out[j] = c;
      ++j;
    }
    return out;
  }

  /** Float32 regression labels; handles int/float/bigint targets. */
  private regressionLabels(target: DG.Column): Float32Array {
    const n = target.length;
    const labels = new Float32Array(n);
    for (let i = 0; i < n; ++i)
      labels[i] = Number(target.get(i));
    return labels;
  }

  private categoryByIndex(idx: number): string {
    if (this.targetCategories === undefined)
      throw new Error('Predicting fails: undefined categories');
    const maxCategory = this.targetCategories.length - 1;
    return this.targetCategories[Math.max(0, Math.min(idx, maxCategory))];
  }

  /** LIBSVM returns the class label as the category code we trained with. */
  private classificationPrediction(prediction: Float32Array): DG.Column {
    const n = prediction.length;
    if (this.targetWasBool) {
      const predClass = new Array<boolean>(n);
      for (let i = 0; i < n; ++i)
        predClass[i] = this.categoryByIndex(Math.round(prediction[i])) === BOOL_TRUE;
      return DG.Column.fromList(DG.COLUMN_TYPE.BOOL, PRED_NAME, predClass);
    }
    const predClass = new Array<string>(n);
    for (let i = 0; i < n; ++i)
      predClass[i] = this.categoryByIndex(Math.round(prediction[i]));
    return DG.Column.fromList(DG.COLUMN_TYPE.STRING, PRED_NAME, predClass);
  }

  private regressionPrediction(prediction: Float32Array): DG.Column {
    switch (this.targetType) {
    case DG.COLUMN_TYPE.INT:
      return this.intColPrediction(prediction);
    case DG.COLUMN_TYPE.BIG_INT:
      return this.bigIntColPrediction(prediction);
    default:
      return DG.Column.fromFloat32Array(PRED_NAME, prediction);
    }
  }

  private intColPrediction(prediction: Float32Array): DG.Column {
    const n = prediction.length;
    const rawInts = new Int32Array(n);
    for (let i = 0; i < n; ++i)
      rawInts[i] = Math.round(prediction[i]);
    return DG.Column.fromInt32Array(PRED_NAME, rawInts, n);
  }

  private bigIntColPrediction(prediction: Float32Array): DG.Column {
    const n = prediction.length;
    const rawInts = new BigInt64Array(n);
    for (let i = 0; i < n; ++i)
      rawInts[i] = BigInt(Math.round(prediction[i]));
    return DG.Column.fromBigInt64Array(PRED_NAME, rawInts);
  }
}

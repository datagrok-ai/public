import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

const ROWS_EXTRA = 1;
const COLS_EXTRA = 2;
const MIN_COLS_COUNT = 1 + COLS_EXTRA;
const AVGS_NAME = 'Avg-s';
const STDEVS_NAME = 'Stddev-s';
const PRED_NAME = 'predicted';
const DEFAULT_LEARNING_RATE = 0.01;
const DEFAULT_ITER_COUNT = 100;
const DEFAULT_REG_RATE = 0.1;

type DataSpecification = {
  classesCount: number,
  featuresCount: number,
};

type TargetLabelsData = {
  oneHot: Array<Uint8Array>,
  weights: Float32Array,
};

/** Softmax classifier */
export class SoftmaxClassifier {
  private isTrained = false;

  private avgs: Float32Array;
  private stdevs: Float32Array;
  private categories: string[];
  private params: Float32Array[];
  private classesCount = 1;
  private featuresCount = 1;

  constructor(specification?: DataSpecification, packedModel?: Uint8Array) {
    if (specification !== undefined) {
      /** features count */
      const n = specification.featuresCount;

      /** classes count */
      const c = specification.classesCount;

      if (n < 1)
        throw new Error('Incorrect features count');

      if (c < 1)
        throw new Error('Incorrect classes count');

      /** length of arrays */
      const len = n + ROWS_EXTRA;

      // Init model routine
      this.avgs = new Float32Array(len);
      this.stdevs = new Float32Array(len);
      this.categories = new Array<string>(len);
      this.featuresCount = n;
      this.classesCount = c;

      // Xavier initialization scale value
      const xavierScale = 2 * Math.sqrt(6.0 / (c + n));

      // Init params
      this.params = new Array<Float32Array>(c);
      for (let i = 0; i < c; ++i) {
        const current = new Float32Array(len);

        // initialize bias, b
        current[n] = 0;

        //Xavier initialization of weights, w
        for (let j = 0; j < n; ++j)
          current[j] = (Math.random() - 0.5) * xavierScale;

        this.params[i] = current;
      }

      this.isTrained = true; // DELETE THIS
    } else if (packedModel !== undefined) {
      try {
        const modelDf = DG.DataFrame.fromByteArray(packedModel);
        const columns = modelDf.columns;
        const colsCount = columns.length;

        if (colsCount < MIN_COLS_COUNT)
          throw new Error('incorrect columns count');

        this.classesCount = colsCount - COLS_EXTRA;
        this.featuresCount = modelDf.rowCount - ROWS_EXTRA;

        const c = this.classesCount;

        // extract params & categories
        this.params = new Array<Float32Array>(c);
        this.categories = new Array<string>(modelDf.rowCount);

        for (let i = 0; i < c; ++i) {
          const c = columns.byIndex(i);
          this.categories[i] = c.name;

          if (c.type !== DG.COLUMN_TYPE.FLOAT)
            throw new Error('incorrect params columns type');

          this.params[i] = c.getRawData() as Float32Array;
        }

        // extract averages
        const avgsCol = columns.byName(AVGS_NAME);
        if (avgsCol.type !== DG.COLUMN_TYPE.FLOAT)
          throw new Error('incorrect average values column type');
        this.avgs = avgsCol.getRawData() as Float32Array;

        // extract stdevs
        const stdevsCol = columns.byName(STDEVS_NAME);
        if (stdevsCol.type !== DG.COLUMN_TYPE.FLOAT)
          throw new Error('incorrect standard deviations column type');
        this.stdevs = stdevsCol.getRawData() as Float32Array;

        this.isTrained = true;
      } catch (e) {
        throw new Error(`Failed to load model: ${(e instanceof Error ? e.message : 'the platform issue')}`);
      }
    } else
      throw new Error('Softmax classifier not initialized');
  }; // constructor

  /** Return packed softmax classifier */
  public toBytes(): Uint8Array {
    if (!this.isTrained)
      throw new Error('Empty model');

    const c = this.classesCount;
    const columns = new Array<DG.Column>(c + COLS_EXTRA);

    // params columns
    for (let i = 0; i < c; ++i)
      columns[i] = DG.Column.fromFloat32Array(this.categories[i], this.params[i]);

    // averages
    columns[c] = DG.Column.fromFloat32Array(AVGS_NAME, this.avgs);

    // stdevs
    columns[c + 1] = DG.Column.fromFloat32Array(STDEVS_NAME, this.stdevs);

    const modelDf = DG.DataFrame.fromColumns(columns);
    grok.shell.addTableView(modelDf);

    return modelDf.toByteArray();
  } // toBytes

  public fit(features: DG.ColumnList, target: DG.Column, learningRate: number = DEFAULT_LEARNING_RATE,
    iterCount: number = DEFAULT_ITER_COUNT, regRate: number = DEFAULT_REG_RATE): void {
    const n = this.featuresCount;
    const m = target.length; // samples count
    const c = this.classesCount;

    if (features.length !== n)
      throw new Error('Training failes - incorrect features count');

    if ((learningRate <= 0) || (iterCount < 1) || (regRate <= 0))
      throw new Error('Training failes - incorrect fitting hyperparameters');

    // 1. Pre-process data

    // 1.1) Normalize data

    this.extractStats(features);

    /** Normalized data */
    const X = this.normalized(features);
    const transposedX = this.transposed(features);

    //console.log(X);

    //console.log('==============================');

    //console.log(transposedX);

    // 1.2) Classes
    const targetData = this.preprocessedTargets(target);
    const Y = targetData.oneHot;
    const classesWeights = targetData.weights;

    //console.log(Y);
    //console.log(classesWeights);

    // 2. Fitting

    // Routine
    let xBuf: Float32Array;
    let wBuf: Float32Array;
    let zBuf: Float32Array;
    let sum: number;
    let sumExp: number;
    let yTrue: Uint8Array;
    let yPred: Float32Array;
    let dWbuf: Float32Array;
    const Z = new Array<Float32Array>(m);
    for (let i = 0; i < m; ++i)
      Z[i] = new Float32Array(c);
    const dZ = new Array<Float32Array>(c);
    for (let i = 0; i < c; ++i)
      dZ[i] = new Float32Array(m);
    const dW = new Array<Float32Array>(c);
    for (let i = 0; i < c; ++i)
      dW[i] = new Float32Array(n + 1);

    //console.log(this.params);

    // Fitting
    for (let iter = 0; iter < iterCount; ++iter) {
      //console.log(`iter: ${iter}`);

      // 2.1) Forward propagation
      for (let j = 0; j < m; ++j) {
        xBuf = X[j];
        zBuf = Z[j];
        sum = 0;
        sumExp = 0;

        for (let i = 0; i < c; ++i) {
          wBuf = this.params[i];
          sum = wBuf[n];

          for (let k = 0; k < n; ++k)
            sum += wBuf[k] * xBuf[k];

          zBuf[i] = Math.exp(sum);
          sumExp += zBuf[i];
        }

        for (let i = 0; i < c; ++i)
          zBuf[i] /= sumExp;
      }

      // 2.2) Backward propagation

      // 2.2.1) dZ
      for (let j = 0; j < m; ++j) {
        yPred = Z[j];
        yTrue = Y[j];

        for (let i = 0; i < c; ++i)
          dZ[i][j] = yPred[i] - yTrue[i];
      }

      // 2.2.2) dB
      for (let i = 0; i < c; ++i) {
        sum = 0;
        zBuf = dZ[i];

        for (let j = 0; j < m; ++j)
          sum += zBuf[j];

        dW[i][n] = sum / m;
      }

      // 2.2.3) dW
      for (let i = 0; i < c; ++i) {
        zBuf = dZ[i];
        wBuf = dW[i];

        for (let j = 0; j < n; ++j) {
          xBuf = transposedX[j];

          sum = 0;
          for (let k = 0; k < m; ++k)
            sum += zBuf[k] * xBuf[k];

          wBuf[j] = sum / m;
        }
      }

      //console.log(dW);

      // 2.3) Update weights
      for (let i = 0; i < c; ++i) {
        wBuf = this.params[i];
        dWbuf = dW[i];

        for (let j = 0; j < n; ++j)
          wBuf[j] = (1 - learningRate * regRate / m) * wBuf[j] - learningRate * dWbuf[j];

        wBuf[n] -= learningRate * dWbuf[n];
      }
    } // for iter

    //console.log(this.params);

    this.isTrained = true;
  }; // fit

  private extractStats(features: DG.ColumnList): void {
    let j = 0;

    for (const col of features) {
      if ((col.type !== DG.COLUMN_TYPE.INT) && (col.type !== DG.COLUMN_TYPE.FLOAT))
        throw new Error('Training failes - incorrect features type');

      this.avgs[j] = col.stats.avg;
      this.stdevs[j] = col.stats.stdev;

      ++j;
    }
  } // extractStats

  /** Retrun normalized features */
  private normalized(features: DG.ColumnList): Array<Float32Array> {
    const m = features.byIndex(0).length;

    const X = new Array<Float32Array>(m);

    for (let i = 0; i < m; ++i)
      X[i] = new Float32Array(this.featuresCount);

    let j = 0;
    for (const col of features) {
      if ((col.type !== DG.COLUMN_TYPE.INT) && (col.type !== DG.COLUMN_TYPE.FLOAT))
        throw new Error('Training failes - incorrect features type');

      const raw = col.getRawData();
      const avg = this.avgs[j];
      const stdev = this.stdevs[j];

      if (stdev > 0) {
        for (let i = 0; i < m; ++i)
          X[i][j] = (raw[i] - avg) / stdev;
      } else {
        for (let i = 0; i < m; ++i)
          X[i][j] = 0;
      }

      ++j;
    }

    return X;
  } // normalized

  /** Retrun normalized & transposed features */
  private transposed(features: DG.ColumnList): Array<Float32Array> {
    const m = features.byIndex(0).length;
    const n = this.featuresCount;

    const X = new Array<Float32Array>(n);

    for (let i = 0; i < n; ++i)
      X[i] = new Float32Array(m);

    let j = 0;
    for (const col of features) {
      if ((col.type !== DG.COLUMN_TYPE.INT) && (col.type !== DG.COLUMN_TYPE.FLOAT))
        throw new Error('Training failes - incorrect features type');

      const raw = col.getRawData();
      const avg = this.avgs[j];
      const stdev = this.stdevs[j];

      if (stdev > 0) {
        for (let i = 0; i < m; ++i)
          X[j][i] = (raw[i] - avg) / stdev;
      } else {
        for (let i = 0; i < m; ++i)
          X[j][i] = 0;
      }

      ++j;
    }

    return X;
  } // transposed

  /** Return one-hot vectors and classes weights */
  private preprocessedTargets(target: DG.Column): TargetLabelsData {
    if (target.type !== DG.COLUMN_TYPE.STRING)
      throw new Error('Training failes - incorrect target type');

    const cats = target.categories;
    const c = this.classesCount;
    const m = target.length;
    const raw = target.getRawData();

    for (let i = 0; i < c; ++i)
      this.categories[i] = cats[i];

    const Y = new Array<Uint8Array>(m);
    const weights = new Float32Array(c).fill(0);

    for (let i = 0; i < m; ++i)
      Y[i] = new Uint8Array(c).fill(0);

    for (let i = 0; i < m; ++i) {
      Y[i][raw[i]] = 1;
      ++weights[raw[i]];
    }

    for (let i = 0; i < c; ++i)
      weights[i] = m / (c * weights[i]);

    return {
      oneHot: Y,
      weights: weights,
    };
  } // getOneHot

  /** Return prediction column */
  public predict(features: DG.ColumnList): DG.Column {
    if (!this.isTrained)
      throw new Error('Predcition fails: no fitted parameters');

    if (features.length !== this.featuresCount)
      throw new Error('Predcition fails: incorrect features count');

    // Normalize features
    const X = this.normalized(features);

    // Routine items
    const m = X.length;
    const n = this.featuresCount;
    const c = this.classesCount;
    let xBuf: Float32Array;
    let wBuf: Float32Array;
    const Z = new Float32Array(c);
    let sum: number;
    let max: number;
    let argMax: number;
    const predClass = new Array<string>(m);

    // get prediction for each sample
    for (let j = 0; j < m; ++j) {
      xBuf = X[j];
      sum = 0;

      for (let i = 0; i < c; ++i) {
        wBuf = this.params[i];
        sum = wBuf[n];

        for (let k = 0; k < n; ++k)
          sum += wBuf[k] * xBuf[k];

        Z[i] = Math.exp(sum);
      }

      max = Z[0];
      argMax = 0;

      for (let k = 1; k < c; ++k) {
        if (max < Z[k]) {
          max = Z[k];
          argMax = k;
        }
      }

      predClass[j] = this.categories[argMax];
    }

    return DG.Column.fromStrings(PRED_NAME, predClass);
  }
}; // SoftmaxClassifier

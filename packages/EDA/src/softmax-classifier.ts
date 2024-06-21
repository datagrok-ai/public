import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

const ROWS_EXTRA = 1;
const COLS_EXTRA = 2;
const MIN_COLS_COUNT = 1 + COLS_EXTRA;
const AVGS_NAME = 'Avg-s';
const STDEVS_NAME = 'Stddev-s';

type DataSpecification = {
  classesCount: number,
  featuresCount: number,
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

  public fit(features: DG.ColumnList, target: DG.Column): void {
    const n = this.featuresCount;
    const m = target.length; // samples count

    if (features.length !== n)
      throw new Error('Training failes - incorrect features count');

    // 1. Pre-process data

    // 1.1) Normalize data

    this.extractStats(features);

    /** Normalized data */
    const X = this.normalized(features);
    const transposedX = this.transposed(features);

    console.log(X);

    console.log('==============================');

    console.log(transposedX);

    // 1.2) Classes
    const Y = this.getOneHot(target);

    console.log(Y);

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

  private getOneHot(target: DG.Column): Array<Uint8Array> {
    if (target.type !== DG.COLUMN_TYPE.STRING)
      throw new Error('Training failes - incorrect target type');

    const cats = target.categories;
    const c = this.classesCount;
    const m = target.length;
    const raw = target.getRawData();

    for (let i = 0; i < c; ++i)
      this.categories[i] = cats[i];

    const Y = new Array<Uint8Array>(c);

    for (let i = 0; i < c; ++i)
      Y[i] = new Uint8Array(m).fill(0);

    for (let i = 0; i < m; ++i)
      Y[raw[i]][i] = 1;

    return Y;
  } // getOneHot
}; // SoftmaxClassifier

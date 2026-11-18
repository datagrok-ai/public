import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {accuracy, classificationDataset, copyCol, mulberry32, numericFeatures, r2, seededRegression} from './utils';
import {nonNullRows} from '../utils';
import {trainCatBoost, trainDecisionTree, trainExtraTrees, trainHistGradientBoosting, trainKnn, trainLda, trainMlp,
  trainNaiveBayes, trainRandomForest, trainRealMlp, trainSoftmax, trainSvm, trainTabm,
  trainXgboost} from '../crux/methods';
import {trainLinearRegression} from '../regression';
import {trainPlsRegression} from '../pls/pls-ml';

export const BASELINE_FILE = 'System:AppData/Eda/baselines/baseline.json';
export const FOLDS = 5;

const CARS = 'System:DemoFiles/cars.csv';
const IRIS = 'System:DemoFiles/iris.csv';
const WINE = 'System:DemoFiles/winequality.csv';
const DEMOG = 'System:DemoFiles/demog.csv';

export type Task = 'classification' | 'regression';
export type ScoreName = 'accuracy' | 'r2';
export type Scorer = (target: DG.Column, pred: DG.Column) => number;
export type Fit = (df: DG.DataFrame, target: DG.Column) => Promise<Uint8Array>;

export const SCORERS: Record<Task, Scorer> = {classification: accuracy, regression: r2};

export interface BaselineCase {
  id: string;
  file?: string;
  builder?: () => Promise<DG.DataFrame>;
  target: string;
  exclude: string[];
  task: Task;
}

export interface Baseline {
  created: string;
  eda: string;
  cases: Record<string, {score: number; fitMs: number; scoreName: ScoreName; trainScore?: number}>;
  pca: Record<string, number[][]>;
  pls: {tScores: number[][]; uScores: number[][]; xLoadings: number[][]; yLoadings: number[]; coef: number[];
    vip: number[]};
}

type Part = {x: DG.DataFrame, y: DG.Column};

/** The features and the target of the training rows and of the test rows `i % FOLDS === k` */
function split(df: DG.DataFrame, k: number, features: string[], target: string): {train: Part, test: Part} {
  const n = df.rowCount;
  const part = (mask: DG.BitSet) => ({x: df.clone(mask, features), y: df.col(target)!.clone(mask)});
  return {train: part(DG.BitSet.create(n, (i) => i % FOLDS !== k)),
    test: part(DG.BitSet.create(n, (i) => i % FOLDS === k))};
}

/** Mean test-fold score and mean fit time over the `FOLDS` folds `i % FOLDS`; with `trainScore`, also the mean
 * training-set score of the fold models */
export async function crossScore(df: DG.DataFrame, target: string, featureNames: string[],
  fit: (features: DG.ColumnList, target: DG.Column) => Promise<Uint8Array>,
  predict: (features: DG.ColumnList, bytes: Uint8Array) => Promise<DG.Column>,
  score: Scorer, trainScore?: Scorer): Promise<{score: number, fitMs: number, trainScore?: number}> {
  let total = 0;
  let trainTotal = 0;
  let fitMs = 0;
  for (let k = 0; k < FOLDS; ++k) {
    const {train, test} = split(df, k, featureNames, target);
    const t0 = performance.now();
    const bytes = await fit(train.x.columns, train.y);
    fitMs += performance.now() - t0;
    total += score(test.y, await predict(test.x.columns, bytes));
    if (trainScore !== undefined)
      trainTotal += trainScore(train.y, await predict(train.x.columns, bytes));
  }
  return {score: total / FOLDS, fitMs: fitMs / FOLDS,
    ...(trainScore !== undefined ? {trainScore: trainTotal / FOLDS} : {})};
}

async function demog(columns: string[]): Promise<DG.DataFrame> {
  const df = await grok.dapi.files.readCsv(DEMOG);
  return df.clone(nonNullRows(df.columns.byNames(columns)), columns);
}

async function oneHotDemog(): Promise<DG.DataFrame> {
  const df = await demog(['AGE', 'HEIGHT', 'WEIGHT', 'SEX', 'RACE']);
  const n = df.rowCount;
  const cols = ['AGE', 'HEIGHT', 'WEIGHT'].map((name) => copyCol(df.col(name)!));
  for (const name of ['SEX', 'RACE']) {
    const col = df.col(name)!;
    for (const cat of col.categories) {
      const raw = Int32Array.from({length: n}, (_, i) => col.get(i) === cat ? 1 : 0);
      if (raw.some((v) => v === 1))
        cols.push(DG.Column.fromInt32Array(`${name}=${cat}`, raw));
    }
  }
  return DG.DataFrame.fromColumns(cols);
}

function illConditioned(): DG.DataFrame {
  const n = 200;
  const rng = mulberry32(7);
  const t = new Float64Array(n);
  const tEps = new Float64Array(n);
  const z = new Float64Array(n);
  const y = new Float64Array(n);
  for (let i = 0; i < n; ++i) {
    const noise = 2 * rng() - 1;
    z[i] = 2 * rng() - 1;
    const noise2 = 2 * rng() - 1;
    t[i] = i + 1;
    tEps[i] = t[i] + 1e-6 * noise;
    y[i] = 1 + 2 * t[i] - z[i] + 0.05 * noise2;
  }
  return DG.DataFrame.fromColumns([DG.Column.fromFloat64Array('t', t), DG.Column.fromFloat64Array('t_eps', tEps),
    DG.Column.fromFloat64Array('z', z), DG.Column.fromFloat64Array('y', y)]);
}

export const CASES: BaselineCase[] = [
  {id: 'cars-price', file: CARS, target: 'price', exclude: [], task: 'regression'},
  {id: 'iris-species', file: IRIS, target: 'Species', exclude: ['col 1'], task: 'classification'},
  {id: 'iris-petal', file: IRIS, target: 'Petal.Length', exclude: ['col 1'], task: 'regression'},
  {id: 'wine-type', file: WINE, target: 'type', exclude: ['quality'], task: 'classification'},
  {id: 'wine-quality', file: WINE, target: 'quality', exclude: [], task: 'regression'},
  {id: 'demog-sex', builder: () => demog(['AGE', 'HEIGHT', 'WEIGHT', 'SEX']), target: 'SEX', exclude: [],
    task: 'classification'},
  {id: 'demog-weight', builder: () => demog(['AGE', 'HEIGHT', 'WEIGHT']), target: 'WEIGHT', exclude: [],
    task: 'regression'},
  {id: 'synth-cls-100k', builder: async () => classificationDataset(100000, 10, true, 42), target: 'Labels',
    exclude: [], task: 'classification'},
  {id: 'synth-reg-100k', builder: async () => seededRegression(100000, 10, 42), target: 'y', exclude: [],
    task: 'regression'},
  {id: 'synth-cls-10k', builder: async () => classificationDataset(10000, 10, true, 42), target: 'Labels',
    exclude: [], task: 'classification'},
  {id: 'ill-linreg', builder: async () => illConditioned(), target: 'y', exclude: [], task: 'regression'},
  {id: 'onehot-linreg', builder: oneHotDemog, target: 'AGE', exclude: [], task: 'regression'},
];

/** Baseline keys `<method>/<case id>` */
export const METHOD_CASES: Record<string, string[]> = {
  xgboost: ['cars-price', 'iris-species', 'wine-type', 'wine-quality', 'synth-cls-100k', 'synth-reg-100k'],
  svm: ['iris-species', 'iris-petal', 'wine-type', 'synth-cls-10k'],
  softmax: ['iris-species', 'wine-type', 'demog-sex'],
  linreg: ['cars-price', 'wine-quality', 'demog-weight', 'ill-linreg', 'onehot-linreg'],
  pls: ['cars-price'],
};

export const PLS_COMPONENTS = 3;
export const PCA_COMPONENTS = 4;

/** Every method with its default hyperparameters */
export const FITS: Record<string, Fit> = {
  xgboost: (df, t) => trainXgboost(df, t, 20, 0.3, 6, 1, 0),
  svm: (df, t) => trainSvm(df, t, 'RBF', 1, 0, 3, 0, 0.1),
  softmax: (df, t) => trainSoftmax(df, t, 100, 0.1, 0.001),
  linreg: (df, t) => trainLinearRegression(df, t, 0, 0, 1000),
  pls: (df, t) => trainPlsRegression(df, t, PLS_COMPONENTS),
  randomForest: (df, t) => trainRandomForest(df, t, 100, 10, 1, 'auto', 0),
  extraTrees: (df, t) => trainExtraTrees(df, t, 100, 10, 1, 'auto', 0),
  histGradientBoosting: (df, t) => trainHistGradientBoosting(df, t, 100, 0.1, 31, 0, 20, 0, 0),
  catBoost: (df, t) => trainCatBoost(df, t, 100, 0.1, 6, 3, 0),
  decisionTree: (df, t) => trainDecisionTree(df, t, 10, 1, 0),
  knn: (df, t) => trainKnn(df, t, 5, 'Uniform'),
  naiveBayes: trainNaiveBayes,
  lda: (df, t) => trainLda(df, t, 0),
  mlp: (df, t) => trainMlp(df, t, '100', 'ReLU', 200, 0.001, 0.0001, 0),
  realMlp: (df, t) => trainRealMlp(df, t, 'TD', 128, '256,256,256', 0, 256, 0.2, 0),
  tabm: (df, t) => trainTabm(df, t, 8, 128, 'Piecewise linear', 32, 16, 0.002, 0.1, 0.2, 0),
};
export const PCA_FLAGS: Record<string, [boolean, boolean]> =
  {FF: [false, false], TF: [true, false], FT: [false, true], TT: [true, true]};

/** The case's frame and its feature names: numeric columns except the target and `exclude` */
export async function loadCase(id: string): Promise<{c: BaselineCase, df: DG.DataFrame, features: string[]}> {
  const c = CASES.find((item) => item.id === id)!;
  const df = c.file !== undefined ? await grok.dapi.files.readCsv(c.file) : await c.builder!();
  return {c, df, features: numericFeatures(df, c.target, c.exclude).names()};
}

export function columnValues(col: DG.Column): number[] {
  return Array.from(col.getRawData().subarray(0, col.length) as Float32Array);
}

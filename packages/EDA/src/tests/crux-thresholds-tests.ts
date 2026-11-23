// The Crux methods against the pre-Crux implementation (files/baselines/baseline.json, PLAN-ru section 4)

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {before, category, expect, test} from '@datagrok-libraries/test/src/test';

import {Baseline, BASELINE_FILE, FITS, METHOD_CASES, PCA_COMPONENTS, PCA_FLAGS, PLS_COMPONENTS, SCORERS,
  columnValues, crossScore, loadCase} from './baseline';
import {getPlsAnalysis} from '../pls/pls-tools';
import {computePCA} from '../eda-tools';
import {predictColumn, r2} from './utils';

type Result = {key: string, old: number, now: number, threshold: string, pass: boolean};
type Run = {score: number, fitMs: number, trainScore?: number};
type Rule = (now: number, old: number) => {threshold: string, pass: boolean};

/** Cases where the old model's out-of-fold score is not comparable: n ~ p, the training-set score is compared */
const TRAIN_SCORE_CASES = ['linreg/cars-price'];

const TIMEOUT = 600000;

/** Baseline fit times below this are dominated by jitter: they are only logged */
const MIN_TIMED_MS = 100;

let baseline: Baseline;

/** Old XGBoost / SVM rounded regression predictions of an integer target; so does the comparison */
function rounds(method: string, df: DG.DataFrame, target: string): boolean {
  return (method === 'xgboost' || method === 'svm') && df.col(target)!.type === DG.COLUMN_TYPE.INT;
}

async function predict(features: DG.ColumnList, bytes: Uint8Array, round: boolean): Promise<DG.Column> {
  const pred = await predictColumn(features.byIndex(0).dataFrame, bytes);
  return round ? DG.Column.fromFloat32Array(pred.name, Float32Array.from(pred.getRawData(), Math.round)) : pred;
}

/** New out-of-fold score and fit time of every case of `method` (and the training-set R² of the fold models for
 * {@link TRAIN_SCORE_CASES}); the first fit warms the module up */
async function methodRuns(method: string): Promise<Map<string, Run>> {
  const fit = FITS[method];
  const out = new Map<string, Run>();
  for (const [idx, id] of METHOD_CASES[method].entries()) {
    const {c, df, features} = await loadCase(id);
    if (idx === 0)
      await fit(df.clone(null, features), df.col(c.target)!);
    const key = `${method}/${id}`;
    const round = rounds(method, df, c.target);
    out.set(key, await crossScore(df, c.target, features, (f, t) => fit(f.byIndex(0).dataFrame, t),
      (f, bytes) => predict(f, bytes, round), SCORERS[c.task], TRAIN_SCORE_CASES.includes(key) ? r2 : undefined));
  }
  return out;
}

function report(title: string, results: Result[]): void {
  for (const r of results)
    console.log(`THRESHOLD ${title} ${r.key} old=${r.old} new=${r.now} threshold=${r.threshold} pass=${r.pass}`);
  const failed = results.filter((r) => !r.pass);
  expect(failed.length, 0, failed.map((r) => `${r.key}: old ${r.old}, new ${r.now}, threshold ${r.threshold}`)
    .join('; '));
}

/** Every case of `method` against its baseline score (the training-set one for {@link TRAIN_SCORE_CASES}) */
async function compare(method: string, title: string, rule: Rule, extra: Result[] = []): Promise<void> {
  const results: Result[] = [];
  for (const [key, run] of await methodRuns(method)) {
    const train = run.trainScore !== undefined;
    const old = train ? baseline.cases[key].trainScore! : baseline.cases[key].score;
    const now = train ? run.trainScore! : run.score;
    results.push({key: train ? `${key} (training set)` : key, old, now, ...rule(now, old)});
  }
  report(title, [...results, ...extra]);
}

const withinPercent: Rule = (now, old) => {
  const tol = 0.01 * Math.abs(old);
  return {threshold: `|new - old| <= ${tol}`, pass: Math.abs(now - old) <= tol};
};

const notLowerThan = (margin: number): Rule => (now, old) => ({threshold: `>= ${old - margin}`,
  pass: now >= old - margin});

/** Wall-clock checks are machine-dependent: they assert only in a benchmark run, otherwise just log */
async function fitTimes(method: string, enforce: boolean): Promise<void> {
  const assert = enforce && DG.Test.isInBenchmark;
  const results: Result[] = [];
  for (const [key, {fitMs}] of await methodRuns(method)) {
    const old = baseline.cases[key].fitMs;
    const timed = assert && old >= MIN_TIMED_MS;
    const pass = !timed || fitMs <= 2 * old;
    const threshold = timed ? `<= ${2 * old}` : assert ? `informational (baseline < ${MIN_TIMED_MS} ms)` : 'none';
    results.push({key, old, now: fitMs, threshold, pass});
  }
  report(`${method} fitMs`, results);
}

function maxAbs(a: ArrayLike<number>): number {
  let m = 0;
  for (let i = 0; i < a.length; ++i)
    m = Math.max(m, Math.abs(a[i]));
  return m;
}

function maxDiff(a: ArrayLike<number>, b: ArrayLike<number>, sign: number = 1): number {
  let m = 0;
  for (let i = 0; i < a.length; ++i)
    m = Math.max(m, Math.abs(a[i] - sign * b[i]));
  return m;
}

function dot(a: ArrayLike<number>, b: ArrayLike<number>): number {
  let s = 0;
  for (let i = 0; i < a.length; ++i)
    s += a[i] * b[i];
  return s;
}

function closeResult(key: string, old: number[], now: ArrayLike<number>, sign: number = 1): Result {
  const tol = 1e-6 * Math.max(1, maxAbs(old));
  const err = maxDiff(now, old, sign);
  return {key, old: maxAbs(old), now: err, threshold: `max error <= ${tol}`, pass: err <= tol};
}

category('Crux thresholds', () => {
  before(async () => {
    baseline = JSON.parse(await grok.dapi.files.readAsText(BASELINE_FILE));
  });

  test('XGBoost metrics', async () => compare('xgboost', 'xgboost metrics', withinPercent), {timeout: TIMEOUT});

  test('SVM metrics', async () => compare('svm', 'svm metrics', withinPercent), {timeout: TIMEOUT});

  test('Softmax accuracy', async () => compare('softmax', 'softmax accuracy', notLowerThan(0.02)),
    {timeout: TIMEOUT});

  test('Linear regression R2', async () => {
    const {c, df, features} = await loadCase('onehot-linreg');
    const x = df.clone(null, features);
    const pred = await predictColumn(x, await FITS.linreg(x, df.col(c.target)!));
    const finite = pred.toList().every((v) => v !== null && Number.isFinite(v));
    await compare('linreg', 'linreg R2', notLowerThan(1e-5), [{key: 'linreg/onehot-linreg predictions', old: NaN,
      now: finite ? 1 : 0, threshold: 'all finite', pass: finite}]);
  }, {timeout: TIMEOUT});

  test('PLS regression R2', async () => compare('pls', 'pls R2', notLowerThan(1e-5)), {timeout: TIMEOUT});

  test('PLS analysis matches', async () => {
    const {c, df, features} = await loadCase('cars-price');
    const out = await getPlsAnalysis({table: df, features: df.clone(null, features).columns,
      predict: df.col(c.target)!, components: PLS_COMPONENTS, isQuadratic: false, names: undefined});
    const old = baseline.pls;
    const signs = old.tScores.map((t, k) => Math.sign(dot(t, columnValues(out.tScores[k]))));

    const results = [
      closeResult('pls/coef', old.coef, columnValues(out.regressionCoefficients)),
      closeResult('pls/vip', old.vip, columnValues(out.vip)),
      closeResult('pls/yLoadings', old.yLoadings, columnValues(out.yLoadings).map((v, k) => v * signs[k])),
    ];
    for (let k = 0; k < PLS_COMPONENTS; ++k) {
      results.push(closeResult(`pls/tScores/${k + 1}`, old.tScores[k], columnValues(out.tScores[k]), signs[k]));
      results.push(closeResult(`pls/uScores/${k + 1}`, old.uScores[k], columnValues(out.uScores[k]), signs[k]));
      results.push(closeResult(`pls/xLoadings/${k + 1}`, old.xLoadings[k], columnValues(out.xLoadings[k]), signs[k]));
    }
    report('pls analysis', results);
  }, {timeout: TIMEOUT});

  test('PCA matches up to sign', async () => {
    const iris = await grok.dapi.files.readCsv('System:DemoFiles/iris.csv');
    const cols = iris.clone(null, ['Sepal.Length', 'Sepal.Width', 'Petal.Length', 'Petal.Width']).columns;
    const results: Result[] = [];
    for (const [flag, [center, scale]] of Object.entries(PCA_FLAGS)) {
      const out = await computePCA(iris, cols, PCA_COMPONENTS, center, scale);
      for (let k = 0; k < PCA_COMPONENTS; ++k) {
        const old = baseline.pca[flag][k];
        const now = columnValues(out.columns.byIndex(k));
        const sign = maxDiff(now, old, 1) <= maxDiff(now, old, -1) ? 1 : -1;
        results.push(closeResult(`pca/${flag}/${k + 1}`, old, now, sign));
      }
    }
    report('pca', results);
  }, {timeout: TIMEOUT});

  test('XGBoost time', async () => fitTimes('xgboost', true), {timeout: TIMEOUT, benchmark: true});

  test('SVM time', async () => fitTimes('svm', true), {timeout: TIMEOUT, benchmark: true});

  test('Linear regression / Softmax time', async () => {
    await fitTimes('linreg', false);
    await fitTimes('softmax', false);
  }, {timeout: TIMEOUT, benchmark: true});
});

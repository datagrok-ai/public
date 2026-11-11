import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, expectArray, expectFloat, test} from '@datagrok-libraries/test/src/test';
import {defaultHyperparameters, enabledRules, Engine, hasProbabilities, Hyperparameters, hyperparametersOf, isEnabled}
  from '../engines/engine';
import {LoopProgress} from '../engines/engine-calls';
import {EngineRegistry} from '../engines/engine-registry';
import {ForgeError} from '../forge-error';
import {MetricId, metricsOf} from '../metrics/metrics';
import {replayPostprocessing} from '../preparation/pipeline';
import {browserSearch} from '../training/browser-search';
import {DEFAULT_VALIDATION, mulberry32, ValidationSettings} from '../training/k-fold';
import {fnv1a, FNV_OFFSET} from '../storage/dataset-fingerprint';
import {defaultFeatures} from '../training/default-features';
import {bestCandidates, bestSearch, bestSearchPlan, canContinue, chooseLeader, COMBINED_LIMITS, continueSearch,
  COST_RANK, costRankOf, drawnPoints,
  FLOAT_DIGITS, FROM_BELOW_TO, holdOutSplit, isInsideSpace, Leader, leadersOf, MAX_FAILED_WITHOUT_COMPLETE,
  NO_PROGRESS, PARALLEL_FITS, rangeProblem, requestFingerprint, runSearch, SearchCandidate, SearchChoice,
  searchableEngines, SearchPlan, searchPlanOf, SearchRecord, searchRecordFrom, searchRecordOf, SearchSpace,
  searchSpaceOf, SearchState, SearchStrategy, Trial, trialPoints, tuneCandidate, tuneSearch, tuneSearchPlan,
  userSpace}
  from '../training/search';
import {TrainingRequest, trainModel} from '../training/train-model';
import {columnsOf, edaEngine, IRIS, MEASUREMENTS, NO_PREPARATION_STEPS, PREDICT_PROBABILITY, requestOf, svm,
  twoSpeciesIris, xgboost} from './test-data';

const TIMEOUT = 120000;
const XGBOOST_SPACE = ['iterations int 10 100', 'eta float 0.03 0.5', 'maxDepth int 1 10', 'lambda float 0.1 10',
  'alpha float 0.01 10', 'subsample float 0.5 1', 'colsample float 0.5 1'];
/** The rules EDA declares for SVM (the work order's format). */
const SVM_RULES = '{"gamma": {"kernel": ["RBF", "Polynomial", "Sigmoid"]}, "degree": {"kernel": ["Polynomial"]}, ' +
  '"coef0": {"kernel": ["Polynomial", "Sigmoid"]}, "epsilon": {"task": ["regression"]}}';

/** The one limit: [maxTrials] trials, or without them 600 s. */
const limitOf = (maxTrials: number | undefined) => maxTrials === undefined ? {budgetSeconds: 600} : {maxTrials};

/** A Tune plan of [engine] on iris Species from its defaults, limited by [maxTrials] (600 s without them). */
async function irisPlan(metric: MetricId, seed: number, maxTrials: number | undefined, engine: Engine,
  validation: ValidationSettings = DEFAULT_VALIDATION): Promise<SearchPlan> {
  const iris = await grok.data.files.openTable(IRIS);
  const request = {...await requestOf(columnsOf(iris, MEASUREMENTS), iris.getCol('Species'), NO_PREPARATION_STEPS,
    engine), seed, validation};
  const {search} = holdOutSplit(request, seed);
  return searchPlanOf('tune', search, [tuneCandidate(engine, request.hyperparameters)],
    {metric, ...limitOf(maxTrials)}, seed);
}

/** [rows] rows of four seeded features and a two-class target (about two thirds a), or with [isRegression] a numeric
 * one. */
async function twoClassRequest(rows: number = 600, isRegression: boolean = false): Promise<TrainingRequest> {
  const random = mulberry32(3);
  const features = ['x1', 'x2', 'x3', 'x4'].map((name) =>
    DG.Column.fromFloat32Array(name, Float32Array.from({length: rows}, () => random())));
  const x1 = features[0];
  const target = isRegression ?
    DG.Column.fromFloat32Array('y', Float32Array.from({length: rows}, (_, i) => x1.getNumber(i) + 0.3 * random())) :
    DG.Column.fromStrings('class', Array.from({length: rows}, (_, i) =>
      x1.getNumber(i) + 0.3 * random() < 0.8 ? 'a' : 'b'));
  return requestOf(features, target);
}

/** A task-bar indicator stand-in whose cancel is clicked [ms] after it is created. */
function cancelledAfter(ms: number): LoopProgress {
  const indicator = {canceled: false, update: () => {}};
  setTimeout(() => indicator.canceled = true, ms);
  return indicator;
}

function share(col: DG.Column, cls: string): number {
  let n = 0;
  for (let i = 0; i < col.length; i++) {
    if (col.get(i) === cls)
      n++;
  }
  return n / col.length;
}

/** The values of [point] outside [space]: `<name>=<value>` for each. */
function outside(space: SearchSpace, point: Hyperparameters): string[] {
  return space.filter((d) => {
    const v = point[d.name];
    if (d.type === 'categorical')
      return !d.choices!.includes(v);
    return typeof v !== 'number' || v < d.low! || v > d.high! || (d.type === 'int' && !Number.isInteger(v));
  }).map((d) => `${d.name}=${point[d.name]}`);
}

const states = (trials: Trial[], state: Trial['state']) => trials.filter((t) => t.state === state);

/** The [n]-th of [points], counted from 0. */
function nth(points: Iterator<Hyperparameters>, n: number): Hyperparameters {
  for (let i = 0; i < n; i++)
    points.next();
  return points.next().value as Hyperparameters;
}

/** Point [n] of a candidate's stream: 0 is [start], then the drawn points. */
const pointOf = (space: SearchSpace, start: Hyperparameters, seed: number, name: string, n: number) =>
  n === 0 ? {...start} : nth(drawnPoints(space, start, seed, name), n - 1);

const trialPoint = (candidate: SearchCandidate, seed: number, task: 'classification', n: number) =>
  nth(trialPoints(candidate, seed, task), n);

/** A stand-in of [engine] whose train declares [mlspace] and, when given, [mlenabled]. */
function withSpace(engine: Engine, mlspace: string, mlenabled?: string): Engine {
  const train = engine.functions.train!;
  const stub: DG.Func = Object.create(train);
  Object.defineProperty(stub, 'options', {value: {...train.options, mlspace,
    ...(mlenabled === undefined ? {} : {mlenabled})}});
  return {...engine, functions: {...engine.functions, train: stub}};
}

/** A stand-in of [engine] named [name] whose train runs on the server (`meta.mlserver`). */
function onServer(engine: Engine, name: string): Engine {
  const train = engine.functions.train!;
  const stub: DG.Func = Object.create(train);
  Object.defineProperty(stub, 'options', {value: {...train.options, mlserver: 'true'}});
  return {...engine, name, functions: {...engine.functions, train: stub}};
}

/** A Find best model plan of the methods [names] on cars price by its default features, or on iris Species, limited
 * by [maxTrials] (600 s without them), by successive halving unless [strategy] says otherwise. */
async function bestPlan(table: 'cars' | 'iris', names: string[], metric: MetricId, seed: number,
  maxTrials: number | undefined, topK?: number, strategy: SearchStrategy = 'halving'): Promise<SearchPlan> {
  const isCars = table === 'cars';
  const frame = await grok.data.files.openTable(isCars ? 'System:DemoFiles/cars.csv' : IRIS);
  const target = frame.getCol(isCars ? 'price' : 'Species');
  const features = isCars ? defaultFeatures(frame, target) : columnsOf(frame, MEASUREMENTS);
  const candidates = bestCandidates(names.map(edaEngine), EngineRegistry.discover(), false);
  const request = {...await requestOf(features, target, NO_PREPARATION_STEPS, candidates[0].engine), seed};
  const {search} = holdOutSplit(request, seed);
  return searchPlanOf('best', search, candidates, {metric, ...limitOf(maxTrials), topK, strategy}, seed);
}

const CARS_METHODS = ['XGBoost', 'Linear Regression', 'PLS Regression'];
const HALVING_METHODS = ['Softmax', 'LDA', 'Decision Tree', 'k-Nearest Neighbors', 'XGBoost', 'Hist Gradient Boosting',
  'Extra Trees', 'Random Forest'];

/** The candidates successive halving with [topK] gives [trials] in turn (every candidate completing its first trial
 * and none exhausted): each round's candidates, 2^round trials each in turn, the best (by its best complete score so
 * far, ties by index) half kept but never fewer than topK; then the finalists in turn. Also the round each candidate
 * was eliminated in. */
function halvingOrder(plan: SearchPlan, trials: Trial[],
  topK: number): {order: number[]; eliminated: (number | undefined)[]} {
  const bestOf = (k: number, upTo: number) => Math.max(...trials.slice(0, upTo)
    .filter((t) => t.candidate === k && t.state === 'complete').map((t) => plan.direction * t.score!));
  const eliminated: (number | undefined)[] = plan.candidates.map(() => undefined);
  const order: number[] = [];
  let pool = plan.candidates.map((_, k) => k);
  for (let round = 0; order.length < trials.length; round++) {
    for (let i = 0; i < 2 ** round; i++)
      order.push(...pool);
    const ranked = [...pool].sort((a, b) => bestOf(b, order.length) - bestOf(a, order.length) || a - b);
    const keep = Math.min(ranked.length, Math.max(topK, Math.ceil(ranked.length / 2)));
    for (const k of ranked.slice(keep))
      eliminated[k] = round;
    pool = ranked.slice(0, keep);
    if (pool.length <= topK) {
      while (order.length < trials.length)
        order.push(...pool);
    }
  }
  return {order: order.slice(0, trials.length), eliminated};
}

const RECORD: SearchRecord = {mode: 'tune', metric: 'accuracy', budgetSeconds: 60, seed: 1, trials: 3, completed: 2,
  pruned: 0, failed: 1, elapsedMs: 1234, stopped: 'trials', holdOut: {rows: 120, fraction: 0.2},
  leaders: [{method: 'XGBoost', hyperparameters: {iterations: 20, eta: 0.3, booster: 'gbtree', on: true}, score: 0.95,
    se: 0.01, foldScores: [0.9, 1, 0.95, 0.95, 0.95], trials: 3, chosen: true}]};

category('Search', () => {
  test('holdOutSplit holds out a stratified part above 500 rows and nothing below', async () => {
    const iris = await grok.data.files.openTable(IRIS);
    const small = await requestOf(columnsOf(iris, MEASUREMENTS), iris.getCol('Species'));
    const none = holdOutSplit(small, 1);
    expect(none.holdOut === undefined && none.search === small, true, 'iris got a hold out');

    const request = await twoClassRequest();
    const {search, holdOut} = holdOutSplit(request, 1);
    if (holdOut === undefined)
      throw new Error('No hold out on 600 rows');
    expect(holdOut.rows, 120);
    expect(holdOut.target.length, 120);
    expect(search.target.length, 480);
    expect(search.prepared.columns[0].length, 480);
    expect(search.keptRows?.trueCount, 480);
    const expected = share(request.target, 'a') * 120;
    const held = share(holdOut.target, 'a') * 120;
    expect(Math.abs(held - expected) < 1, true, `${held} of class a held out, ${expected} expected`);
    const again = holdOutSplit(request, 1);
    expect(again.search.keptRows?.toBinaryString(), search.keptRows?.toBinaryString());
    expect(holdOutSplit(request, 2).search.keptRows?.toBinaryString() !== search.keptRows?.toBinaryString(), true,
      'Another seed holds out the same rows');

    for (const isRegression of [false, true]) {
      const large = await twoClassRequest(10999, isRegression);
      const split = holdOutSplit(large, 1);
      expect(split.holdOut?.rows, 2000, `${isRegression ? 'regression' : 'classification'} on 10,999 rows`);
      expect(split.search.target.length, 8999);
    }
  }, {timeout: TIMEOUT});

  test('trainModel measures metrics.test on the held out rows and keeps search', async () => {
    const {search, holdOut} = holdOutSplit(await twoClassRequest(), 1);
    const result = await trainModel({...search, holdOut, search: RECORD});
    const accuracy = result.metrics.test?.accuracy ?? NaN;
    expect(accuracy >= 0 && accuracy <= 1, true, `test accuracy ${accuracy}`);
    expect(JSON.stringify(result.search), JSON.stringify(RECORD));
    expect(result.rowCount, search.target.length);
    expect(result.splitting.isStratified, true);

    const iris = await grok.data.files.openTable(IRIS);
    const regression = await trainModel(await requestOf(columnsOf(iris, MEASUREMENTS.slice(0, 3)),
      iris.getCol('Petal.Width')));
    expect(regression.metrics.test === undefined, true, 'A test without a hold out');
    expect(regression.splitting.isStratified, false);
  }, {timeout: TIMEOUT});

  test('searchRecordOf keeps a record and drops a malformed one', async () => {
    expect(JSON.stringify(searchRecordOf(JSON.parse(JSON.stringify(RECORD)))), JSON.stringify(RECORD));
    expect(searchRecordOf({...RECORD, mode: 'x'}) === undefined, true, 'mode x kept');
    expect(searchRecordOf({mode: 'x'}) === undefined, true, '{mode: x} kept');
    expect(searchRecordOf({...RECORD, leaders: [{...RECORD.leaders[0], score: 'high'}]}) === undefined, true,
      'A text score kept');
    expect(searchRecordOf({...RECORD, holdOut: undefined})?.holdOut === undefined, true, 'A hold out made up');
  });

  test('searchSpaceOf reads EDA spaces and refuses a bad one', async () => {
    const space = searchSpaceOf(xgboost()) ?? [];
    expectArray(space.map((d) => `${d.name} ${d.type} ${d.low} ${d.high}`), XGBOOST_SPACE);
    expect(searchSpaceOf(edaEngine('Naive Bayes')) === undefined, true, 'Naive Bayes has a space');
    const errorOf = (engine: Engine) => {
      try {
        searchSpaceOf(engine);
        return '';
      } catch (e) {
        return e instanceof ForgeError ? e.message : String(e);
      }
    };
    const unknown = errorOf(withSpace(xgboost(), '{"speed": {"type": "int", "low": 1, "high": 2}}'));
    expect(unknown.includes('names no setting'), true, unknown);
    const aboveMax = errorOf(withSpace(xgboost(), '{"iterations": {"type": "int", "low": 200, "high": 300}}'));
    expect(aboveMax.includes('at most 100'), true, aboveMax);
    const shorthand = searchSpaceOf(withSpace(xgboost(), '{"maxDepth": [2, 5], "eta": [0.1, 0.5]}')) ?? [];
    expectArray(shorthand.map((d) => `${d.name} ${d.type}`), ['maxDepth int', 'eta float']);
  });

  test('pointOf is seeded per method and the first point is the start', async () => {
    const space = searchSpaceOf(xgboost())!;
    const start: Hyperparameters = {iterations: 20, eta: 0.3, maxDepth: 6, lambda: 1, alpha: 0, seed: 0};
    expect(JSON.stringify(pointOf(space, start, 7, 'XGBoost', 0)), JSON.stringify(start));
    for (let n = 1; n <= 5; n++) {
      const point = pointOf(space, start, 7, 'XGBoost', n);
      expect(JSON.stringify(pointOf(space, start, 7, 'XGBoost', n)), JSON.stringify(point));
      expectArray(outside(space, point), []);
      expect(point.seed, 0);
    }
    const first = (name: string) => JSON.stringify(pointOf(space, start, 7, name, 1));
    expect(first('SVM') !== first('XGBoost'), true, 'SVM and XGBoost share their points');
  });

  test('Drawn floats keep 3 significant digits within their bounds; ints and the start are not rounded', async () => {
    const engine = withSpace(xgboost(), '{"eta": [0.03, 0.5], "lambda": {"type": "float", "low": 0.1, "high": 10, ' +
      '"log": true}, "iterations": [10, 100]}');
    const space = searchSpaceOf(engine)!;
    expectArray(space.map((d) => `${d.name} ${d.type} ${d.log}`),
      ['eta float undefined', 'lambda float true', 'iterations int undefined']);
    const start = {...defaultHyperparameters(engine), eta: 0.2453, lambda: 2.1234, iterations: 20};
    // The stream drawnPoints draws from: one value per dimension per point.
    const random = mulberry32(7 ^ fnv1a(FNV_OFFSET, new TextEncoder().encode('XGBoost')));
    const points = drawnPoints(space, start, 7, 'XGBoost');
    const wrong: string[] = [];
    for (let n = 0; n < 200; n++) {
      const point = points.next().value;
      const [eta, lambda, iterations] = [point.eta, point.lambda, point.iterations] as number[];
      const raw = [0.03 + random() * 0.47, Math.exp(Math.log(0.1) + random() * (Math.log(10) - Math.log(0.1)))];
      const whole = 10 + Math.floor(random() * 91);
      for (const [i, [v, low, high]] of [[eta, 0.03, 0.5], [lambda, 0.1, 10]].entries()) {
        if (Number(v.toPrecision(FLOAT_DIGITS)) !== v || v < low || v > high || Math.abs(v - raw[i]) > 0.005 * raw[i])
          wrong.push(`point ${n}: ${v} for ${raw[i]} in ${low}-${high}`);
      }
      if (iterations !== whole)
        wrong.push(`point ${n}: iterations ${iterations}, ${whole} drawn`);
    }
    expectArray(wrong, []);
    const first = trialPoint(tuneCandidate(engine, start), 7, 'classification', 0);
    expect(`${first.eta} ${first.lambda} ${first.iterations}`, '0.2453 2.1234 20');
  });

  test('A float range narrower than its rounding exhausts the search', async () => {
    const plan = await irisPlan('accuracy', 1, 10, xgboost());
    const engine = plan.candidates[0].engine;
    const pinned = plan.candidates[0].space!.map((d) => d.name).filter((name) => name !== 'eta');
    const candidate = tuneCandidate(engine, plan.request.hyperparameters,
      {space: {pinned, ranges: {eta: {low: 0.1, high: 0.1001}}}, task: 'classification'});
    const {outcome} = await runSearch({...plan, candidates: [candidate]}, NO_PROGRESS, browserSearch);
    expect(outcome.stopped, 'exhausted');
    expectArray(outcome.trials.map((t) => t.hyperparameters.eta), [0.1]);
  }, {timeout: TIMEOUT});

  test('A tune search on iris is reproducible with a fixed seed', async () => {
    const plan = await irisPlan('balancedAccuracy', 7, 12, xgboost());
    const space = plan.candidates[0].space!;
    const runs = [await runSearch(plan, NO_PROGRESS, browserSearch), await runSearch(plan, NO_PROGRESS, browserSearch)];
    for (const {outcome, chosen} of runs) {
      expect(outcome.trials.length, 12);
      expect(outcome.stopped, 'trials');
      expect(chosen !== undefined, true, 'Nothing chosen');
      // The defaults lie outside the space (XGBoost alpha 0), so trial 0 is drawn too.
      expectArray(outcome.trials.flatMap((t) => outside(space, t.hyperparameters)), []);
    }
    const [a, b] = runs.map((r) => r.leaders[0]);
    expect(JSON.stringify(a.hyperparameters), JSON.stringify(b.hyperparameters));
    expectFloat(a.score, b.score, 1e-12);
    expectArray(a.foldScores.map((s, i) => Math.abs(s - b.foldScores[i]) < 1e-12), a.foldScores.map(() => true));
  }, {timeout: TIMEOUT});

  test('Tune leaves the other methods alone', async () => {
    const {outcome, leaders} = await runSearch(await irisPlan('accuracy', 3, 3, xgboost()), NO_PROGRESS,
      browserSearch);
    expectArray(outcome.trials.map((t) => t.candidate), outcome.trials.map(() => 0));
    expectArray(leaders.map((l) => l.method), ['XGBoost']);
  }, {timeout: TIMEOUT});

  test('The budget stops the search within one trial', async () => {
    const plan = {...await irisPlan('accuracy', 1, undefined, xgboost()), budgetMs: 1500};
    const {outcome} = await runSearch(plan, NO_PROGRESS, browserSearch);
    expect(outcome.stopped, 'budget');
    expect(outcome.elapsedMs < 1500 + 2000, true, `${outcome.elapsedMs} ms`);
    expect(outcome.trials.length >= 1, true, 'No trial');
    expectArray(states(outcome.trials, 'abandoned').map((t) => t.foldScores.length),
      states(outcome.trials, 'abandoned').map(() => 2));
  }, {timeout: TIMEOUT});

  test('Cancel between fits keeps the best so far', async () => {
    const plan = await irisPlan('accuracy', 1, undefined, xgboost());
    const {outcome, chosen} = await runSearch(plan, cancelledAfter(1000), browserSearch);
    expect(outcome.stopped, 'cancelled');
    expect(states(outcome.trials, 'complete').length >= 1, true, 'No complete trial');
    expect(chosen !== undefined, true, 'Nothing chosen');
    expect(states(outcome.trials, 'failed').length, 0);
  }, {timeout: TIMEOUT});

  test('Pruned trials stop after two folds and count against the trial cap', async () => {
    const {outcome, leaders} = await runSearch(await irisPlan('accuracy', 2, 30, xgboost()), NO_PROGRESS,
      browserSearch);
    const pruned = states(outcome.trials, 'pruned');
    console.log(`Search: ${pruned.length} of ${outcome.trials.length} trials pruned`);
    expect(pruned.length > 0, true, 'No trial pruned');
    expect(`${outcome.stopped} ${outcome.trials.length} ${leaders[0].trials}`, 'trials 30 30');
    expectArray(pruned.map((t) => t.foldScores.length), pruned.map(() => 2));
    const complete = states(outcome.trials, 'complete');
    expectArray(complete.map((t) => t.foldScores.length), complete.map(() => 5));
  }, {timeout: TIMEOUT});

  test('The winner equals a regular training with its settings', async () => {
    const plan = await irisPlan('balancedAccuracy', 7, 6, xgboost());
    const {chosen} = await runSearch(plan, NO_PROGRESS, browserSearch);
    if (chosen === undefined)
      throw new Error('Nothing chosen');
    const request = {...plan.request, hyperparameters: chosen.hyperparameters};
    const first = await trainModel(request);
    expectFloat(first.metrics.validation.balancedAccuracy ?? NaN, chosen.score, 1e-9);
    const second = await trainModel(request);
    expect(JSON.stringify(second.metrics.validation), JSON.stringify(first.metrics.validation));
  }, {timeout: TIMEOUT});

  test('An AUC search fails the trials of a method without probabilities, each counted by the cap', async () => {
    const plan = await irisPlan('auc', 1, 2, svm());
    const {outcome, chosen} = await runSearch(plan, NO_PROGRESS, browserSearch);
    expect(`${outcome.stopped} ${outcome.trials.length}`, 'trials 2');
    const wrong = outcome.trials.filter((t) => t.state !== 'failed' || !t.error?.includes('needs class probabilities'));
    expectArray(wrong.map((t) => `${t.state} ${t.error}`), []);
    expect(chosen === undefined, true, 'A leader was chosen');
  }, {timeout: TIMEOUT});

  test('A log loss search minimizes', async () => {
    const plan = await irisPlan('logLoss', 1, 6, xgboost());
    expect(plan.direction, -1);
    const {outcome, leaders, chosen} = await runSearch(plan, NO_PROGRESS, browserSearch);
    const complete = states(outcome.trials, 'complete');
    expect(complete.length > 0, true, 'No complete trial');
    expectArray(complete.filter((t) => !(t.score! > 0)).map((t) => t.score), []);
    expectArray(leaders.map((l) => l.score), [...leaders.map((l) => l.score)].sort((a, b) => a - b));
    const best = complete.reduce((a, b) => b.score! < a.score! ? b : a);
    const n = best.foldScores.length;
    const mean = best.foldScores.reduce((sum, s) => sum + s, 0) / n;
    const se = Math.sqrt(best.foldScores.reduce((sum, s) => sum + (s - mean) ** 2, 0) / (n - 1) / n);
    if (chosen === undefined)
      throw new Error('Nothing chosen');
    expect(chosen.score >= best.score! && chosen.score <= best.score! + se + 1e-12, true,
      `chosen ${chosen.score}, the lowest ${best.score} with SE ${se}`);
  }, {timeout: TIMEOUT});

  test('A brier search of SVM fails its trials', async () => {
    const plan = await irisPlan('brier', 1, 2, svm());
    const {outcome, chosen} = await runSearch(plan, NO_PROGRESS, browserSearch);
    expect(outcome.trials.length >= 1, true, 'No trial');
    const wrong = outcome.trials.filter((t) => t.state !== 'failed' || !t.error?.includes('needs class probabilities'));
    expectArray(wrong.map((t) => `${t.state} ${t.error}`), []);
    expect(chosen === undefined, true, 'A leader was chosen');
  }, {timeout: TIMEOUT});

  test('The trials score the labels the step gives', async () => {
    const iris = await twoSpeciesIris();
    const species = iris.getCol('Species');
    const seed = 3;
    const request = {...await requestOf(columnsOf(iris, MEASUREMENTS), species,
      {...PREDICT_PROBABILITY, threshold: 0.9}, xgboost()), seed};
    const {options} = request.prepared;
    expect(options.binaryClassificationThreshold, 0.9);
    const {search} = holdOutSplit(request, seed);
    const plan = searchPlanOf('tune', search, [tuneCandidate(xgboost(), request.hyperparameters)],
      {metric: 'accuracy', maxTrials: 1}, seed);
    const {outcome} = await runSearch(plan, NO_PROGRESS, browserSearch);
    const [trial] = states(outcome.trials, 'complete');
    if (trial === undefined)
      throw new Error(`No complete trial: ${outcome.trials.map((t) => t.error).join('; ')}`);
    const result = await trainModel({...request, hyperparameters: trial.hyperparameters});
    const {prediction, probabilities} = result.outOfFold;
    const labels = replayPostprocessing(prediction, options, probabilities);
    expectFloat(trial.score!, metricsOf('classification', species, labels, 'versicolor').accuracy ?? NaN, 1e-9);
  }, {timeout: TIMEOUT});

  test('A search where every trial fails stops after 10 failures', async () => {
    const plan = {...await irisPlan('auc', 1, undefined, svm()), budgetMs: 600000};
    const {outcome, leaders, chosen} = await runSearch(plan, NO_PROGRESS, browserSearch);
    expect(outcome.stopped, 'failed');
    expect(outcome.trials.length, MAX_FAILED_WITHOUT_COMPLETE);
    expect(states(outcome.trials, 'failed').length, MAX_FAILED_WITHOUT_COMPLETE);
    expect(leaders.length === 0 && chosen === undefined, true, 'A leader was chosen');
    expect(outcome.elapsedMs < 60000, true, `${outcome.elapsedMs} ms`);
  }, {timeout: TIMEOUT});

  test('The one-standard-error rule picks the simplest trial within one SE', async () => {
    const plan = await irisPlan('accuracy', 1, undefined, xgboost());
    const trial = (number: number, maxDepth: number, score: number, foldScores: number[],
      state: Trial['state'] = 'complete'): Trial => ({candidate: 0, number, state, score, foldScores, fitMs: 0,
      hyperparameters: {iterations: 20, eta: 0.3, maxDepth, lambda: 1, alpha: 1}});
    // The best: 0.95, SE about 0.018, so 0.932 and above are within one SE.
    const best = trial(0, 8, 0.95, [0.9, 1, 0.95, 0.92, 0.98]);
    const simpler = trial(1, 3, 0.94, [0.93, 0.95, 0.94, 0.94, 0.94]);
    const far = trial(2, 1, 0.9, [0.9, 0.9, 0.9, 0.9, 0.9]);
    const later = trial(3, 3, 0.945, [0.94, 0.95, 0.945, 0.945, 0.945]);
    const failed = trial(4, 1, NaN, [], 'failed');
    const leaders = leadersOf(plan, {trials: [best, later, simpler, far, failed], elapsedMs: 0,
      stopped: 'budget'});
    expect(leaders.length, 1);
    const [leader] = leaders;
    expect(leader.hyperparameters.maxDepth, 3);
    expect(leader.score, 0.94);
    expectArray(leader.foldScores, simpler.foldScores);
    expect(leader.trials, 5);
    expect(chooseLeader(plan, leaders) === leader && leader.chosen, true, 'The leader is not chosen');

    const alone = leadersOf(plan, {trials: [best, far], elapsedMs: 0, stopped: 'budget'});
    expect(alone[0].hyperparameters.maxDepth, 8);
  });

  test('enabledRules reads meta.mlenabled; isEnabled holds when any condition does', async () => {
    const rulesOf = (mlenabled?: string) => enabledRules(withSpace(xgboost(), '{}', mlenabled));
    expect(JSON.stringify(rulesOf('')), '{}');
    const malformed = ['{"gamma": ', '["gamma"]', '{"gamma": 1}', '{"gamma": {"kernel": {"ge": 1}}}', '{"x": []}'];
    for (const text of malformed)
      expect(rulesOf(text) === undefined, true, `${text} was read`);
    const svmRules = rulesOf(SVM_RULES)!;
    const cases: [string, Hyperparameters, 'classification' | 'regression' | undefined, boolean][] = [
      ['gamma', {kernel: 'Linear'}, 'classification', false], ['gamma', {kernel: 'RBF'}, 'classification', true],
      ['degree', {kernel: 'Sigmoid'}, 'regression', false], ['coef0', {kernel: 'Sigmoid'}, 'regression', true],
      ['epsilon', {kernel: 'RBF'}, 'classification', false], ['epsilon', {kernel: 'RBF'}, 'regression', true],
      ['epsilon', {}, undefined, true], ['cost', {kernel: 'Linear'}, 'classification', true]];
    expectArray(cases.map(([name, values, task]) => isEnabled(svmRules, name, values, task)),
      cases.map(([, , , expected]) => expected));
    const anyOf = rulesOf('{"seed": [{"subsample": {"lt": 1}}, {"colsample": {"lt": 1}}], ' +
      '"iterations": {"alpha": {"gt": 0}}}')!;
    expectArray([isEnabled(anyOf, 'seed', {subsample: 1, colsample: 1}, 'regression'),
      isEnabled(anyOf, 'seed', {subsample: 1, colsample: 0.8}, 'regression'),
      isEnabled(anyOf, 'iterations', {alpha: 0}, 'regression'),
      isEnabled(anyOf, 'iterations', {alpha: 0.5}, undefined)],
    [false, true, false, true]);
  });

  test('A point\'s dimensions its method ignores take the defaults, so such points are one', async () => {
    const engine = withSpace(xgboost(), '{"maxDepth": [1, 10], "lambda": [0.1, 10.0]}',
      '{"lambda": {"maxDepth": {"gt": 5}}}');
    const defaults = defaultHyperparameters(engine);
    const candidate = tuneCandidate(engine, defaults);
    const points = Array.from({length: 40}, (_, n) => trialPoint(candidate, 7, 'classification', n + 1));
    const shallow = points.filter((p) => Number(p.maxDepth) <= 5);
    expect(shallow.length > 0 && shallow.length < points.length, true, `${shallow.length} of 40 shallow points`);
    expectArray(shallow.map((p) => p.lambda), shallow.map(() => defaults.lambda));
    expect(points.filter((p) => Number(p.maxDepth) > 5).every((p) => p.lambda !== defaults.lambda), true,
      'A deep point lost its lambda');
    expect(new Set(shallow.map((p) => JSON.stringify(p))).size <= 5, true, 'Shallow points differ in lambda');
  });

  test('Tune fixes the settings SVM ignores to their defaults (EDA\'s meta.mlenabled)', async () => {
    const engine = svm();
    expect((enabledRules(engine) ?? {}).gamma !== undefined, true, 'EDA\'s SVM declares no meta.mlenabled');
    const defaults = defaultHyperparameters(engine);
    const candidate = tuneCandidate(engine, defaults);
    for (let n = 1; n <= 40; n++) {
      const p = trialPoint(candidate, 3, 'classification', n);
      const wrong = [...(p.kernel === 'Linear' && p.gamma !== defaults.gamma ? ['gamma'] : []),
        ...(p.kernel !== 'Polynomial' && p.degree !== defaults.degree ? ['degree'] : []),
        ...(p.epsilon !== defaults.epsilon ? ['epsilon'] : [])];
      expect(wrong.length, 0, `Point ${n} keeps ${wrong.join(', ')}: ${JSON.stringify(p)}`);
    }
  });

  test('userSpace leaves out pinned and rule-disabled settings and keeps the user\'s bounds and choices', async () => {
    const engine = svm();
    const space = searchSpaceOf(engine)!;
    const rules = enabledRules(engine)!;
    const linear = {...defaultHyperparameters(engine), kernel: 'Linear'};
    const pinned = userSpace(space, {pinned: ['kernel'], ranges: {cost: {low: 0.1, high: 10}}}, rules, linear,
      'classification');
    expectArray(pinned.map((d) => `${d.name} ${d.type} ${d.low} ${d.high} ${d.log} ${d.simpler}`),
      ['cost float 0.1 10 true low']);
    const tuned = userSpace(space, {pinned: [], ranges: {kernel: {choices: ['Polynomial', 'RBF']}}}, rules, linear,
      'regression');
    expectArray(tuned.map((d) => d.name), ['kernel', 'cost', 'gamma', 'degree', 'epsilon']);
    expectArray(tuned[0].choices ?? [], ['RBF', 'Polynomial']);
    expect(tuned[0].simpler, 'first');
    const candidate = tuneCandidate(engine, linear, {space: {pinned: ['kernel', 'cost'], ranges: {}},
      task: 'classification'});
    expect(`${candidate.space?.length} ${candidate.pinned}`, '0 kernel,cost');
  });

  test('rangeProblem checks the order, the min and max, whole numbers and the log scale', async () => {
    const engine = xgboost();
    const space = searchSpaceOf(engine)!;
    const property = (name: string) => hyperparametersOf(engine).find((p) => p.name === name)!;
    const problem = (name: string, low: number | null, high: number | null) =>
      rangeProblem(space.find((d) => d.name === name)!, property(name), low, high);
    expectArray([problem('iterations', 20, 50), problem('iterations', 50, 20), problem('iterations', 20, 20),
      problem('iterations', 2.5, 20), problem('iterations', null, 20), problem('eta', 0.1, 0.2)],
    [null, FROM_BELOW_TO, FROM_BELOW_TO, 'The bounds must be whole numbers', 'Enter both bounds', null]);
    expect(problem('iterations', 20, 500), '“To” must be at most 100');
    const log = withSpace(engine, '{"lambda": {"type": "float", "low": 0.1, "high": 10, "log": true}}');
    const lambda = searchSpaceOf(log)![0];
    expect(rangeProblem(lambda, property('lambda'), 0, 1), '“From” must be above 0');
  });

  test('The first trial is the start only when every value it uses lies inside the space', async () => {
    const engine = xgboost();
    const space = searchSpaceOf(engine)!;
    const rules = enabledRules(engine) ?? {};
    const defaults = defaultHyperparameters(engine);
    expect(isInsideSpace(space, defaults, rules, 'classification'), false, 'XGBoost alpha 0 lies inside');
    const drawn = trialPoint(tuneCandidate(engine, defaults), 7, 'classification', 0);
    expectArray(outside(space, drawn), []);
    const start = {...defaults, alpha: 1};
    expect(JSON.stringify(trialPoint(tuneCandidate(engine, start), 7, 'classification', 0)), JSON.stringify(start));
    const narrowed = tuneCandidate(engine, start, {space: {pinned: [], ranges: {alpha: {low: 2, high: 5}}},
      task: 'classification'});
    expect(JSON.stringify(trialPoint(narrowed, 7, 'classification', 0)) !== JSON.stringify(start), true,
      'A start outside the user\'s range is tried');

    const svmEngine = svm();
    const svmDefaults = {...defaultHyperparameters(svmEngine), kernel: 'Linear', gamma: 0};
    expect(isInsideSpace(searchSpaceOf(svmEngine)!, svmDefaults, enabledRules(svmEngine)!, 'classification'), true,
      'A gamma the linear kernel ignores counts');
  });

  test('A Tune record keeps the space searched and the pinned settings; a malformed space is dropped', async () => {
    const plan = await irisPlan('accuracy', 1, 1, xgboost());
    const engine = plan.candidates[0].engine;
    const tuned = {...plan, candidates: [tuneCandidate(engine, plan.request.hyperparameters,
      {space: {pinned: ['maxDepth'], ranges: {iterations: {low: 20, high: 40}}}, task: 'classification'})]};
    const record = searchRecordFrom(tuned, {trials: [], elapsedMs: 0, stopped: 'budget'}, []);
    expectArray(record.space?.pinned ?? [], ['maxDepth']);
    expectArray(record.space?.searched.map((d) => `${d.name} ${d.low} ${d.high}`) ?? [],
      ['iterations 20 40', 'eta 0.03 0.5', 'lambda 0.1 10', 'alpha 0.01 10', 'subsample 0.5 1', 'colsample 0.5 1']);
    const stored = JSON.parse(JSON.stringify({...RECORD, space: record.space}));
    expect(JSON.stringify(searchRecordOf(stored)?.space), JSON.stringify(record.space));
    expect(JSON.stringify({...searchRecordOf(stored), space: undefined}), JSON.stringify(RECORD));
    const malformed = searchRecordOf({...stored, space: {searched: [{name: 'x', type: 'tree'}], pinned: []}});
    expect(malformed !== undefined && malformed.space === undefined, true, 'A malformed space kept or the record lost');
  });

  test('Tune with a hold out scores each trial on the one split and prunes nothing', async () => {
    const plan = await irisPlan('accuracy', 2, 12, xgboost(), {...DEFAULT_VALIDATION, scheme: 'holdout'});
    const {outcome, chosen} = await runSearch(plan, NO_PROGRESS, browserSearch);
    expect(states(outcome.trials, 'complete').length, 12);
    expect(states(outcome.trials, 'pruned').length, 0);
    expectArray(outcome.trials.map((t) => t.foldScores.length), outcome.trials.map(() => 1));
    if (chosen === undefined)
      throw new Error('Nothing chosen');
    const result = await trainModel({...plan.request, hyperparameters: chosen.hyperparameters});
    expect(result.splitting.scheme, 'holdout');
    expectFloat(result.metrics.validation.accuracy ?? NaN, chosen.score, 1e-9);
  }, {timeout: TIMEOUT});

  test('Parallel and sequential searches give the same leaderboard', async () => {
    const plan = await irisPlan('accuracy', 2, 30, xgboost());
    const sequential = await runSearch({...plan, parallelism: 1}, NO_PROGRESS, browserSearch);
    const parallel = await runSearch({...plan, parallelism: Math.max(2, PARALLEL_FITS)}, NO_PROGRESS, browserSearch);
    const summary = (leaders: Leader[]) =>
      leaders.map((l) => `${l.method} ${JSON.stringify(l.hyperparameters)} ${l.score.toFixed(12)}`);
    expectArray(summary(parallel.leaders), summary(sequential.leaders));
    expect(states(sequential.outcome.trials, 'pruned').length > 0, true, 'No trial pruned');
    const stateList = (trials: Trial[]) => trials.map((t) => `${t.number} ${t.state}`);
    expectArray(stateList(parallel.outcome.trials), stateList(sequential.outcome.trials));
    console.log(`Search: parallelism ${PARALLEL_FITS}, ${Math.round(sequential.outcome.elapsedMs)} ms sequential, ` +
      `${Math.round(parallel.outcome.elapsedMs)} ms parallel`);
  }, {timeout: TIMEOUT});
});

// Find best model: the candidates, the executor over several methods and the record.
category('Search', () => {
  test('bestCandidates sorts by cost rank, keeps a malformed space as a note, carries the probability input and ' +
    'refuses an empty list', async () => {
    const discovery = EngineRegistry.discover();
    const candidates = bestCandidates(CARS_METHODS.map(edaEngine), discovery, false);
    expectArray(candidates.map((c) => `${c.engine.name} ${c.costRank}`),
      ['Linear Regression 0', 'PLS Regression 1', 'XGBoost 7']);
    for (const c of candidates)
      expect(JSON.stringify(c.start), JSON.stringify(defaultHyperparameters(c.engine)), c.engine.name);
    const [malformed] = bestCandidates([withSpace(xgboost(), '{"speed": [1, 2]}')], discovery, false);
    expect(malformed.space === undefined, true, 'A malformed space kept');
    expect(malformed.note?.includes('names no setting'), true, malformed.note);
    for (const on of [true, false]) {
      const [candidate] = bestCandidates([svm()], discovery, on);
      expect(candidate.start.probability, on);
      expect(trialPoint(candidate, 1, 'classification', 1).probability, on);
    }
    let error = '';
    try {
      bestCandidates([], discovery, false);
    } catch (e) {
      error = e instanceof ForgeError ? e.message : String(e);
    }
    expect(error.startsWith('No method can take part in the search'), true, error);
    const zeta = {...xgboost(), name: 'Zeta'};
    expect(costRankOf(zeta, [...discovery, zeta]), COST_RANK.length + discovery.length);
    expectArray(searchableEngines([xgboost(), onServer(xgboost(), 'Remote'), edaEngine('Linear Regression')],
      discovery).map((e) => e.name), ['Linear Regression', 'XGBoost']);
    expectArray([hasProbabilities(svm(), false), hasProbabilities(svm(), true), hasProbabilities(xgboost(), false)],
      [false, true, true]);
  });

  test('A best search tries every candidate in cost order, then tunes the top ones round-robin', async () => {
    const plan = await bestPlan('cars', CARS_METHODS, 'r2', 3, 9, 2);
    const runs = [await runSearch(plan, NO_PROGRESS, browserSearch), await runSearch(plan, NO_PROGRESS, browserSearch)];
    for (const {outcome, leaders, chosen} of runs) {
      const {trials} = outcome;
      console.log(`Search: ${trials.map((t) => `${t.candidate}:${t.number}`).join(' ')}`);
      expect(`${outcome.stopped} ${trials.length}`, 'trials 9');
      expectArray(trials.slice(0, 3).map((t) => `${t.candidate}:${t.number}`), ['0:0', '1:0', '2:0']);
      const screened = trials.slice(0, 3).filter((t) => t.state === 'complete' &&
        plan.candidates[t.candidate].space !== undefined);
      const cycle = screened.sort((a, b) => b.score! - a.score! || a.candidate - b.candidate).slice(0, 2)
        .map((t) => t.candidate);
      const later = trials.slice(3).map((t) => t.candidate);
      expect(trials.slice(3).every((t) => t.number >= 1), true, 'A later trial numbered 0');
      expectArray(later.filter((k) => !cycle.includes(k)), []);
      expect(later[0], cycle[0]);
      // A candidate repeats only once the other of the cycle has no point left.
      const wrong = later.filter((k, i) => i > 0 && k === later[i - 1] && later.slice(i).some((x) => x !== k));
      expectArray(wrong, []);
      expect(leaders.length, 3);
      expect(chosen !== undefined, true, 'Nothing chosen');
      const record = searchRecordFrom(plan, outcome, leaders);
      expectArray(record.methods?.map((m) => `${m.method} ${m.outcome}`) ?? [],
        ['Linear Regression complete', 'PLS Regression complete', 'XGBoost complete']);
      // Three above the top two: the worst is eliminated after the first round.
      expectArray(record.methods?.map((m) => m.eliminatedRound) ?? [],
        plan.candidates.map((_, k) => cycle.includes(k) ? undefined : 0));
    }
    const [a, b] = runs.map((r) => r.leaders);
    expectArray(a.map((l) => JSON.stringify(l.hyperparameters)), b.map((l) => JSON.stringify(l.hyperparameters)));
    expectArray(a.map((l, i) => Math.abs(l.score - b[i].score) < 1e-12), a.map(() => true));
  }, {timeout: TIMEOUT});

  test('A method without a space is screened once and not tuned', async () => {
    const plan = await bestPlan('iris', ['Naive Bayes', 'Decision Tree', 'Softmax'], 'accuracy', 1, 8);
    const {outcome, leaders} = await runSearch(plan, NO_PROGRESS, browserSearch);
    const trialsOf = (name: string) => outcome.trials.filter((t) => plan.candidates[t.candidate].engine.name === name);
    expect(trialsOf('Naive Bayes').length, 1);
    expect(trialsOf('Decision Tree').length >= 1 && trialsOf('Softmax').length >= 1, true, 'A method not screened');
    expect(leaders.some((l) => l.method === 'Naive Bayes'), true, 'Naive Bayes has no leader');
  }, {timeout: TIMEOUT});

  test('Parallel and sequential best searches give the same leaderboard', async () => {
    const plan = await bestPlan('cars', CARS_METHODS, 'r2', 3, 9, 2);
    const sequential = await runSearch({...plan, parallelism: 1}, NO_PROGRESS, browserSearch);
    const parallel = await runSearch({...plan, parallelism: Math.max(2, PARALLEL_FITS)}, NO_PROGRESS, browserSearch);
    const summary = (leaders: Leader[]) =>
      leaders.map((l) => `${l.method} ${JSON.stringify(l.hyperparameters)} ${l.score.toFixed(12)}`);
    expectArray(summary(parallel.leaders), summary(sequential.leaders));
  }, {timeout: TIMEOUT});

  test('A best record lists every candidate\'s outcome and survives a round trip', async () => {
    const plan = await bestPlan('iris', ['Naive Bayes', 'Decision Tree', 'Softmax'], 'accuracy', 1, 8);
    expectArray(plan.candidates.map((c) => c.engine.name), ['Softmax', 'Naive Bayes', 'Decision Tree']);
    const trial = (candidate: number, state: Trial['state'], error?: string): Trial => ({candidate, number: 0, state,
      hyperparameters: {}, foldScores: [0.9], fitMs: 1, ...(state === 'complete' ? {score: 0.9} : {}),
      ...(error === undefined ? {} : {error})});
    const outcome = {trials: [trial(0, 'complete'), trial(1, 'failed', 'Boom')], elapsedMs: 10,
      stopped: 'budget' as const};
    const record = searchRecordFrom(plan, outcome, leadersOf(plan, outcome));
    expect(JSON.stringify(record.methods), JSON.stringify([{method: 'Softmax', trials: 1, outcome: 'complete'},
      {method: 'Naive Bayes', trials: 1, outcome: 'failed', note: 'Boom'},
      {method: 'Decision Tree', trials: 0, outcome: 'unstarted'}]));
    expect(JSON.stringify(searchRecordOf(JSON.parse(JSON.stringify(record)))), JSON.stringify(record));
    const halving = {round: 1, pool: [0], given: [0, 0, 0], eliminated: [null, 0, null], next: 0, finalists: [0]};
    const halved = searchRecordFrom(plan, {...outcome, resume: {drawn: [], exhausted: [], passed: [], halving}}, []);
    expectArray(halved.methods?.map((m) => m.eliminatedRound) ?? [], [undefined, 0, undefined]);
    expect(JSON.stringify(searchRecordOf(JSON.parse(JSON.stringify(halved)))), JSON.stringify(halved));
    expect(searchRecordOf({...halved, methods: [{...halved.methods![1], eliminatedRound: 'first'}]})?.methods ===
      undefined, true, 'A text round kept');
    const malformed = searchRecordOf({...JSON.parse(JSON.stringify(record)), methods: [{method: 1}]});
    expect(malformed !== undefined && malformed.methods === undefined, true,
      'A malformed list kept or the record lost');
  });

  test('chooseLeader picks the best score less SE, a tie by the name as the chart sorts', async () => {
    const plan = await bestPlan('iris', ['LDA', 'XGBoost'], 'accuracy', 1, 8);
    const leader = (method: string, score: number, se: number): Leader => ({method, hyperparameters: {}, score, se,
      foldScores: [], trials: 1, chosen: false});
    // LDA is within one SE of XGBoost, which the one-standard-error rule would choose, but 0.93 beats 0.925.
    expect(chooseLeader(plan, [leader('XGBoost', 0.95, 0.02), leader('LDA', 0.94, 0.015)])?.method, 'XGBoost');
    expect(chooseLeader(plan, [leader('XGBoost', 0.95, 0.03), leader('LDA', 0.94, 0.01)])?.method, 'LDA');
    expect(chooseLeader(plan, [leader('XGBoost', 0.95, 0.02), leader('LDA', 0.94, 0.01)])?.method, 'LDA', 'A tie');
    // Capitals first, as the bar chart's code-unit order.
    expect(chooseLeader(plan, [leader('Softmax', 0.95, 0.02), leader('SVM', 0.95, 0.02),
      leader('Neural Network', 0.95, 0.02)])?.method, 'Neural Network');
    expect(chooseLeader(plan, [leader('Softmax', 0.95, 0.02), leader('SVM', 0.95, 0.02)])?.method, 'SVM');
    const chosen = chooseLeader(plan, [leader('XGBoost', 0.95, 0.02), leader('LDA', 0.9, 0.01)]);
    expect(`${chosen?.method} ${chosen?.chosen}`, 'XGBoost true');
    expect(chooseLeader(plan, []) === undefined, true, 'A leader chosen from none');
    const rmse = await bestPlan('cars', ['Linear Regression', 'XGBoost'], 'rmse', 1, 8);
    expect(rmse.direction, -1);
    expect(chooseLeader(rmse, [leader('XGBoost', 1, 0.1), leader('Linear Regression', 1.05, 0.02)])?.method,
      'Linear Regression');
    expect(chooseLeader(rmse, [leader('XGBoost', 1, 0.01), leader('Linear Regression', 1.05, 0.02)])?.method,
      'XGBoost');
  });

  test('The progress shows the search and what is left of the limit', async () => {
    const plan = await bestPlan('cars', CARS_METHODS, 'rmse', 3, 9, 2);
    const texts: string[] = [];
    await runSearch(plan, {canceled: false, update: (_, text) => texts.push(text)}, browserSearch);
    expect(texts[0], 'Finding best model (0/9)');
    expect(texts[texts.length - 1], 'Finding best model (9/9)');
    const counts = texts.map((t) => Number(/\((\d+)\/9\)/.exec(t)?.[1]));
    expect(counts.every((n, i) => i === 0 || n >= counts[i - 1]), true, texts.join('; '));

    const iris = await grok.data.files.openTable(IRIS);
    const request = await requestOf(columnsOf(iris, MEASUREMENTS), iris.getCol('Species'), NO_PREPARATION_STEPS,
      xgboost());
    const choice: SearchChoice = {metric: 'accuracy', budgetSeconds: 1};
    const timed: string[] = [];
    const first = await tuneSearch(request, xgboost(), request.hyperparameters, choice, 2,
      {canceled: false, update: (_, text) => timed.push(text)}, browserSearch);
    expect(timed[0], 'Tuning XGBoost (1 s left)');
    const timedText = /^Tuning XGBoost \(([01] s left|finishing ([2-9]|\d\d+) trials|finishing 1 trial)\)$/;
    expect(timed.every((t) => timedText.test(t)), true, timed.join('; '));
    const continued: string[] = [];
    await continueSearch(first.state, tuneSearchPlan(request, xgboost(), request.hyperparameters, choice, 2),
      {canceled: false, update: (_, text) => continued.push(text)}, browserSearch);
    expect(continued.length > 0 && continued.every((t) => timedText.test(t)), true, continued.join('; '));

    const long: string[] = [];
    const cancelling = {canceled: false, update: (_: number, text: string) => {
      long.push(text);
      cancelling.canceled = true;
    }};
    await runSearch(tuneSearchPlan(request, xgboost(), request.hyperparameters, {metric: 'accuracy',
      budgetSeconds: 90}, 2).plan, cancelling, browserSearch);
    expect(long[0], 'Tuning XGBoost (2 min left)');
  }, {timeout: TIMEOUT});

  test('Successive halving doubles the trials of the better half each round, keeps topK, then tunes them',
    async () => {
      const plan = {...await bestPlan('iris', HALVING_METHODS, 'accuracy', 5, 26, 1), parallelism: 1};
      expectArray(plan.candidates.map((c) => c.engine.name), [...HALVING_METHODS].sort((a, b) =>
        COST_RANK.indexOf(a) - COST_RANK.indexOf(b)));
      const {outcome, leaders} = await runSearch(plan, NO_PROGRESS, browserSearch);
      const {trials} = outcome;
      console.log(`Search: ${trials.map((t) => `${t.candidate}:${t.number}:${t.state}`).join(' ')}`);
      expect(`${outcome.stopped} ${trials.length}`, 'trials 26');
      expectArray(trials.slice(0, 8).map((t) => `${t.candidate}:${t.number}:${t.state}`),
        plan.candidates.map((_, k) => `${k}:0:complete`));
      // 8 candidates with topK 1: rounds of 8 x 1, 4 x 2 and 2 x 4 trials, then the one finalist.
      const expected = halvingOrder(plan, trials, 1);
      expectArray(trials.map((t) => t.candidate), expected.order);
      expectArray(expected.eliminated.map((r) => r ?? -1).sort(), [-1, 0, 0, 0, 0, 1, 1, 2]);
      const record = searchRecordFrom(plan, outcome, leaders);
      expect(`${record.strategy} ${record.maxTrials}`, 'halving 26');
      expectArray(record.methods?.map((m) => m.eliminatedRound) ?? [], expected.eliminated);
      const finalist = expected.eliminated.indexOf(undefined);
      expect(trials.slice(24).every((t) => t.candidate === finalist), true, 'A non-finalist tuned');
      expect(record.methods?.reduce((sum, m) => sum + m.trials, 0), 26);
    }, {timeout: TIMEOUT});

  test('Halving never goes below topK and leaves out a method without a space', async () => {
    const plan = await bestPlan('iris', ['Naive Bayes', ...HALVING_METHODS.slice(0, 4)], 'accuracy', 5, 9, 3);
    const {outcome, leaders} = await runSearch(plan, NO_PROGRESS, browserSearch);
    const record = searchRecordFrom(plan, outcome, leaders);
    const bayes = record.methods?.find((m) => m.method === 'Naive Bayes');
    expect(`${bayes?.trials} ${bayes?.eliminatedRound}`, '1 undefined');
    // Four with a space above three: ceil(4 / 2) = 2 would go below topK, so one is eliminated.
    expect(record.methods?.filter((m) => m.eliminatedRound === 0).length, 1);
    expect(outcome.trials.length, 9);
  }, {timeout: TIMEOUT});

  test('A continued best search equals one run of the summed trial cap', async () => {
    const plan = {...await bestPlan('iris', HALVING_METHODS, 'accuracy', 5, 10, 1), parallelism: 1};
    const first = await runSearch(plan, NO_PROGRESS, browserSearch);
    expect(first.outcome.trials.length, 10);
    const continued = await runSearch({...plan, maxTrials: 16}, NO_PROGRESS, browserSearch, first.outcome);
    const whole = await runSearch({...plan, maxTrials: 26}, NO_PROGRESS, browserSearch);
    const summary = (trials: Trial[]) => trials.map((t) => `${t.candidate}:${t.number}:${t.state}:` +
      JSON.stringify(t.hyperparameters));
    expectArray(summary(continued.outcome.trials.slice(0, 10)), summary(first.outcome.trials));
    expectArray(summary(continued.outcome.trials), summary(whole.outcome.trials));
    expect(continued.outcome.elapsedMs > first.outcome.elapsedMs, true, 'The time is not cumulative');
    const records = [continued, whole].map((r) => searchRecordFrom(plan, r.outcome, r.leaders).methods);
    expect(JSON.stringify(records[0]), JSON.stringify(records[1]));
    expectArray(continued.leaders.map((l) => `${l.method} ${l.trials}`), whole.leaders.map((l) =>
      `${l.method} ${l.trials}`));
  }, {timeout: TIMEOUT});

  test('A continued Tune keeps its trials, tries no point twice and records the sums', async () => {
    const iris = await grok.data.files.openTable(IRIS);
    const request = await requestOf(columnsOf(iris, MEASUREMENTS), iris.getCol('Species'), NO_PREPARATION_STEPS,
      xgboost());
    const choice: SearchChoice = {metric: 'accuracy', maxTrials: 6};
    const first = await tuneSearch(request, xgboost(), request.hyperparameters, choice, 2, NO_PROGRESS,
      browserSearch);
    const planned = tuneSearchPlan(request, xgboost(), request.hyperparameters, choice, 2);
    const continued = await continueSearch(first.state, planned, NO_PROGRESS, browserSearch);
    const whole = await runSearch({...first.plan, maxTrials: 12}, NO_PROGRESS, browserSearch);
    const points = (trials: Trial[]) => trials.map((t) => `${t.state} ${JSON.stringify(t.hyperparameters)}`);
    expectArray(points(continued.outcome.trials.slice(0, 6)), points(first.outcome.trials));
    expectArray(points(continued.outcome.trials), points(whole.outcome.trials));
    expect(new Set(continued.outcome.trials.map((t) => JSON.stringify(t.hyperparameters))).size, 12);
    const {record} = continued;
    expect(`${record.trials} ${record.maxTrials} ${record.budgetSeconds} ${record.leaders[0].trials}`,
      '12 12 undefined 12');
    expect(record.completed + record.pruned + record.failed, 12);
    expectFloat(continued.chosen!.score, whole.chosen!.score, 1e-12);
    expect(continued.state.plan.maxTrials, 12);
    // A time limit after a trial limit: the record holds both sums.
    const timed = await continueSearch(continued.state, tuneSearchPlan(request, xgboost(), request.hyperparameters,
      {metric: 'accuracy', budgetSeconds: 1}, 2), NO_PROGRESS, browserSearch);
    expect(`${timed.record.maxTrials} ${timed.record.budgetSeconds}`, '12 1');
  }, {timeout: TIMEOUT});

  test('canContinue refuses each change of the data, the validation and the search settings', async () => {
    const iris = await grok.data.files.openTable(IRIS);
    const species = iris.getCol('Species');
    const request = await requestOf(columnsOf(iris, MEASUREMENTS), species, NO_PREPARATION_STEPS, xgboost());
    const start = request.hyperparameters;
    const pinnedChoice: SearchChoice = {metric: 'accuracy', maxTrials: 5, space: {pinned: ['maxDepth'], ranges: {}}};
    const stateOf = (plan: SearchPlan): SearchState => ({plan, outcome: {trials: [], elapsedMs: 0, stopped: 'budget'},
      fingerprint: requestFingerprint(plan.request)});
    const tune = stateOf(tuneSearchPlan(request, xgboost(), start, pinnedChoice, 1).plan);
    const tunePlan = (r: TrainingRequest, s: Hyperparameters, choice: SearchChoice, seed: number = 1) =>
      tuneSearchPlan(r, xgboost(), s, choice, seed).plan;
    const discovery = EngineRegistry.discover();
    const candidates = (names: string[]) => bestCandidates(names.map(edaEngine), discovery, false);
    const best = stateOf(bestSearchPlan(request, candidates(['XGBoost', 'Softmax']), pinnedChoice, 1).plan);
    const twoSpecies = await twoSpeciesIris();
    const cases: [string, SearchState, SearchPlan, boolean][] = [
      ['same, another limit and searched start', tune, tunePlan(request, {...start, eta: 0.11},
        {...pinnedChoice, maxTrials: undefined, budgetSeconds: 300}), true],
      ['metric', tune, tunePlan(request, start, {...pinnedChoice, metric: 'balancedAccuracy'}), false],
      ['seed', tune, tunePlan(request, start, pinnedChoice, 2), false],
      ['validation', tune, tunePlan({...request, validation: {...DEFAULT_VALIDATION, folds: 3}}, start,
        pinnedChoice), false],
      ['features', tune, tunePlan(await requestOf(columnsOf(iris, MEASUREMENTS.slice(0, 3)), species,
        NO_PREPARATION_STEPS, xgboost()), start, pinnedChoice), false],
      ['rows', tune, tunePlan(await requestOf(columnsOf(twoSpecies, MEASUREMENTS), twoSpecies.getCol('Species'),
        NO_PREPARATION_STEPS, xgboost()), start, pinnedChoice), false],
      ['preparation', tune, tunePlan(await requestOf(columnsOf(iris, MEASUREMENTS), species, PREDICT_PROBABILITY,
        xgboost()), start, pinnedChoice), false],
      ['pinned value', tune, tunePlan(request, {...start, maxDepth: Number(start.maxDepth) + 1}, pinnedChoice), false],
      ['space', tune, tunePlan(request, start, {...pinnedChoice, space: {pinned: [], ranges: {}}}), false],
      ['range', tune, tunePlan(request, start, {...pinnedChoice, space: {pinned: ['maxDepth'],
        ranges: {iterations: {low: 20, high: 40}}}}), false],
      ['mode', tune, bestSearchPlan(request, candidates(['XGBoost']), pinnedChoice, 1).plan, false],
      ['same candidates', best, bestSearchPlan(request, candidates(['Softmax', 'XGBoost']),
        {...pinnedChoice, maxTrials: 9}, 1).plan, true],
      ['candidates', best, bestSearchPlan(request, candidates(['XGBoost', 'Softmax', 'LDA']), pinnedChoice, 1).plan,
        false],
      ['strategy', best, bestSearchPlan(request, candidates(['XGBoost', 'Softmax']),
        {...pinnedChoice, strategy: 'halving'}, 1).plan, false],
    ];
    expectArray(cases.map(([name, state, plan]) => `${name} ${canContinue(state, plan)}`),
      cases.map(([name, , , expected]) => `${name} ${expected}`));
    let error = '';
    try {
      await continueSearch(tune, {plan: cases[1][2]}, NO_PROGRESS, browserSearch);
    } catch (e) {
      error = e instanceof ForgeError ? e.message : String(e);
    }
    expect(error.startsWith('The search cannot be continued'), true, error);
  }, {timeout: TIMEOUT});
});

// Find best model's equal split and the one limit.
category('Search', () => {
  const EQUAL_METHODS = ['Softmax', 'LDA', 'Decision Tree', 'k-Nearest Neighbors'];
  const countsOf = (plan: SearchPlan, trials: Trial[]) =>
    plan.candidates.map((_, k) => trials.filter((t) => t.candidate === k).length);

  test('An equal split by trials gives each method its share in cost order, the remainder to the first, the ' +
    'defaults first', async () => {
    const plan = {...await bestPlan('iris', HALVING_METHODS, 'accuracy', 5, 19, undefined, 'equal'), parallelism: 1};
    expect(plan.strategy, 'equal');
    const {outcome, leaders} = await runSearch(plan, NO_PROGRESS, browserSearch);
    const {trials} = outcome;
    // 19 over 8: 2 each, the first 3 one more.
    expectArray(countsOf(plan, trials), [3, 3, 3, 2, 2, 2, 2, 2]);
    expect(`${outcome.stopped} ${trials.length}`, 'trials 19');
    const order = trials.map((t) => t.candidate);
    expectArray(order, [...order].sort((a, b) => a - b));
    for (const [k, candidate] of plan.candidates.entries()) {
      const first = trials.find((t) => t.candidate === k)!;
      expect(`${first.number} ${JSON.stringify(first.hyperparameters)}`,
        `0 ${JSON.stringify(trialPoint(candidate, plan.seed, 'classification', 0))}`, candidate.engine.name);
    }
    const record = searchRecordFrom(plan, outcome, leaders);
    expect(`${record.strategy} ${record.maxTrials} ${record.budgetSeconds}`, 'equal 19 undefined');
    expectArray(record.methods?.map((m) => m.eliminatedRound) ?? [], plan.candidates.map(() => undefined));
    expect(JSON.stringify(searchRecordOf(JSON.parse(JSON.stringify(record)))), JSON.stringify(record));
  }, {timeout: TIMEOUT});

  test('An equal split by trials passes the trials a method leaves unused on', async () => {
    // Naive Bayes, without a space, is tried once: LDA after it gets the other 9 trials, not 5.
    const plan = await bestPlan('iris', ['LDA', 'Naive Bayes'], 'accuracy', 1, 10, undefined, 'equal');
    const {outcome} = await runSearch(plan, NO_PROGRESS, browserSearch);
    expectArray(plan.candidates.map((c) => c.engine.name), ['Naive Bayes', 'LDA']);
    expect(`${outcome.stopped} ${countsOf(plan, outcome.trials).join(' ')}`, 'trials 1 9');
  }, {timeout: TIMEOUT});

  test('An equal split by time passes the time a method leaves unused on', async () => {
    // Naive Bayes, without a space, is tried once: LDA after it gets the whole time, not half.
    const plan = await bestPlan('iris', ['LDA', 'Naive Bayes'], 'accuracy', 1, undefined, undefined, 'equal');
    const timed = {...plan, budgetMs: 3000};
    const {outcome} = await runSearch(timed, NO_PROGRESS, browserSearch);
    expectArray(plan.candidates.map((c) => c.engine.name), ['Naive Bayes', 'LDA']);
    expect(`${outcome.stopped} ${countsOf(plan, outcome.trials)[0]}`, 'budget 1');
    expect(outcome.elapsedMs >= 3000, true, `${outcome.elapsedMs} ms`);

    const each = await bestPlan('iris', EQUAL_METHODS, 'accuracy', 1, undefined, undefined, 'equal');
    const shared = await runSearch({...each, budgetMs: 4000}, NO_PROGRESS, browserSearch);
    const counts = countsOf(each, shared.outcome.trials);
    expect(counts.every((n) => n >= 1), true, `${counts}`);
    const order = shared.outcome.trials.map((t) => t.candidate);
    expectArray(order, [...order].sort((a, b) => a - b));
  }, {timeout: TIMEOUT});

  test('A continued equal split splits the new limit equally again', async () => {
    const iris = await grok.data.files.openTable(IRIS);
    const request = await requestOf(columnsOf(iris, MEASUREMENTS), iris.getCol('Species'), NO_PREPARATION_STEPS,
      xgboost());
    const candidates = () => bestCandidates(EQUAL_METHODS.map(edaEngine), EngineRegistry.discover(), false);
    const choice: SearchChoice = {metric: 'accuracy', maxTrials: 8, strategy: 'equal'};
    const first = await bestSearch(request, candidates(), choice, 3, NO_PROGRESS, browserSearch);
    expectArray(countsOf(first.plan, first.outcome.trials), [2, 2, 2, 2]);
    const planned = bestSearchPlan(request, candidates(), {...choice, maxTrials: 4}, 3);
    const continued = await continueSearch(first.state, planned, NO_PROGRESS, browserSearch);
    expectArray(countsOf(planned.plan, continued.outcome.trials), [3, 3, 3, 3]);
    const summary = (trials: Trial[]) => trials.map((t) => `${t.candidate}:${t.number}`);
    expectArray(summary(continued.outcome.trials.slice(0, 8)), summary(first.outcome.trials));
    expect(`${continued.record.strategy} ${continued.record.maxTrials}`, 'equal 12');
  }, {timeout: TIMEOUT});

  test('A plan has one limit; an equal split by trials needs one per method', async () => {
    const iris = await grok.data.files.openTable(IRIS);
    const request = await requestOf(columnsOf(iris, MEASUREMENTS), iris.getCol('Species'), NO_PREPARATION_STEPS,
      xgboost());
    const candidates = bestCandidates(EQUAL_METHODS.map(edaEngine), EngineRegistry.discover(), false);
    const errorOf = (choice: SearchChoice) => {
      try {
        bestSearchPlan(request, candidates, choice, 1);
        return '';
      } catch (e) {
        return e instanceof ForgeError ? e.message : String(e);
      }
    };
    expect(errorOf({budgetSeconds: 60, maxTrials: 10}), COMBINED_LIMITS);
    expect(errorOf({maxTrials: 3}), 'At least one trial per method: 4.');
    expect(errorOf({maxTrials: 3, strategy: 'halving'}), '');
    const {plan} = bestSearchPlan(request, candidates, {}, 1);
    expect(`${plan.budgetMs} ${plan.maxTrials} ${plan.strategy}`, '60000 undefined equal');
    const tune = tuneSearchPlan(request, xgboost(), request.hyperparameters, {maxTrials: 2, strategy: 'equal'}, 1);
    expect(`${tune.plan.budgetMs} ${tune.plan.maxTrials} ${tune.plan.strategy}`, 'undefined 2 undefined');
    expect(searchRecordOf({...RECORD, budgetSeconds: undefined}) === undefined, true, 'A record without a limit read');
    expect(searchRecordOf({...RECORD, strategy: 'fast'}) === undefined, true, 'An unknown strategy read');
  });
});

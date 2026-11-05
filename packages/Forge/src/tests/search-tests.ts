import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, expectArray, expectFloat, test} from '@datagrok-libraries/test/src/test';
import {defaultHyperparameters, enabledRules, Engine, Hyperparameters, hyperparametersOf, isEnabled}
  from '../engines/engine';
import {LoopProgress} from '../engines/engine-calls';
import {ForgeError} from '../forge-error';
import {MetricId} from '../metrics/metrics';
import {browserSearch} from '../training/browser-search';
import {DEFAULT_VALIDATION, mulberry32, ValidationSettings} from '../training/k-fold';
import {fnv1a, FNV_OFFSET} from '../storage/dataset-fingerprint';
import {chooseLeader, drawnPoints, FLOAT_DIGITS, FROM_BELOW_TO, holdOutSplit, isInsideSpace, Leader, leadersOf,
  MAX_FAILED_WITHOUT_COMPLETE, NO_PROGRESS, PARALLEL_FITS, rangeProblem, runSearch, SearchCandidate, SearchPlan,
  searchPlanOf, SearchRecord, searchRecordFrom, searchRecordOf, SearchSpace, searchSpaceOf, Trial, trialPoints,
  tuneCandidate, userSpace} from '../training/search';
import {TrainingRequest, trainModel} from '../training/train-model';
import {columnsOf, edaEngine, IRIS, MEASUREMENTS, NO_PREPARATION_STEPS, requestOf, svm, xgboost} from './test-data';

const TIMEOUT = 120000;
const XGBOOST_SPACE = ['iterations int 10 100', 'eta float 0.03 0.5', 'maxDepth int 1 10', 'lambda float 0.1 10',
  'alpha float 0.01 10', 'subsample float 0.5 1', 'colsample float 0.5 1'];
/** The rules EDA declares for SVM (the work order's format). */
const SVM_RULES = '{"gamma": {"kernel": ["RBF", "Polynomial", "Sigmoid"]}, "degree": {"kernel": ["Polynomial"]}, ' +
  '"coef0": {"kernel": ["Polynomial", "Sigmoid"]}, "epsilon": {"task": ["regression"]}}';

/** A Tune plan of [engine] on iris Species from its defaults: 600 s unless [maxTrials] stops it first. */
async function irisPlan(metric: MetricId, seed: number, maxTrials: number | undefined, engine: Engine,
  validation: ValidationSettings = DEFAULT_VALIDATION): Promise<SearchPlan> {
  const iris = await grok.data.files.openTable(IRIS);
  const request = {...await requestOf(columnsOf(iris, MEASUREMENTS), iris.getCol('Species'), NO_PREPARATION_STEPS,
    engine), seed, validation};
  const {search} = holdOutSplit(request, seed);
  return searchPlanOf('tune', search, [tuneCandidate(engine, request.hyperparameters)],
    {metric, budgetSeconds: 600, maxTrials}, seed);
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
      expect(states(outcome.trials, 'complete').length, 12);
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

  test('Pruned trials stop after two folds', async () => {
    const {outcome} = await runSearch(await irisPlan('accuracy', 2, 30, xgboost()), NO_PROGRESS, browserSearch);
    const pruned = states(outcome.trials, 'pruned');
    console.log(`Search: ${pruned.length} of ${outcome.trials.length} trials pruned`);
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

  test('An AUC search fails the trials of a method without probabilities', async () => {
    const plan = await irisPlan('auc', 1, 2, svm());
    const {outcome, chosen} = await runSearch(plan, NO_PROGRESS, browserSearch);
    expect(outcome.trials.length >= 1, true, 'No trial');
    const wrong = outcome.trials.filter((t) => t.state !== 'failed' || !t.error?.includes('needs class probabilities'));
    expectArray(wrong.map((t) => `${t.state} ${t.error}`), []);
    expect(chosen === undefined, true, 'A leader was chosen');
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
    // A parallel search may launch a few more trials after the last complete one; the ones both ran agree.
    const stateList = (trials: Trial[]) => trials.slice(0, sequential.outcome.trials.length)
      .map((t) => `${t.number} ${t.state}`);
    expectArray(stateList(parallel.outcome.trials), stateList(sequential.outcome.trials));
    console.log(`Search: parallelism ${PARALLEL_FITS}, ${Math.round(sequential.outcome.elapsedMs)} ms sequential, ` +
      `${Math.round(parallel.outcome.elapsedMs)} ms parallel`);
  }, {timeout: TIMEOUT});
});

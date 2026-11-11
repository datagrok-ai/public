import * as DG from 'datagrok-api/dg';
import {defaultHyperparameters, enabledRules, EnabledRules, Engine, givesProbabilities, Hyperparameters, isChoice,
  isEnabled, isParallel, isServerEngine, probabilityInput, searchSpaceText, visibleHyperparametersOf,
  withProbabilityInput} from '../engines/engine';
import {LoopProgress} from '../engines/engine-calls';
import {errorMessage, ForgeError} from '../forge-error';
import {ModelTask} from '../generated/db';
import {METRIC_IDS, METRIC_LABELS, MetricId, taskOf} from '../metrics/metrics';
import {isScoreMetric} from '../metrics/report';
import {isRecord} from '../preparation/preparation-options';
import {datasetFingerprint, fnv1a, FNV_OFFSET} from '../storage/dataset-fingerprint';
import {classGroups, groupQuotas, mulberry32, shuffledRows} from './k-fold';
import {HoldOut, retrainsLive, TrainingRequest, validationSplit} from './train-model';

export type SearchMode = 'tune' | 'best';
/** Find best model's split of the limit: an equal share per method, or successive halving. */
export type SearchStrategy = 'equal' | 'halving';
export type SearchStop = 'budget' | 'cancelled' | 'exhausted' | 'trials' | 'failed';
export interface SpaceDimension {
  name: string;
  type: 'int' | 'float' | 'categorical';
  low?: number;
  high?: number;
  log?: boolean;
  step?: number;
  choices?: (string | number | boolean)[];
  simpler?: 'low' | 'high' | 'first';
}
export type SearchSpace = SpaceDimension[];
/** The settings left at their start values although the space has them ([pinned]) are kept with the candidate. */
export interface SearchCandidate { engine: Engine; space?: SearchSpace; start: Hyperparameters; costRank: number;
  note?: string; pinned?: string[] }
/** A tuned setting's part of the search: its bounds, or its choices. */
export interface TunedRange { low?: number; high?: number; choices?: (string | number | boolean)[] }
/** The user's part of Tune: the settings pinned to their values and the ranges of the tuned ones. */
export interface TuneSpace { pinned: string[]; ranges: {[name: string]: TunedRange} }
export interface SearchPlan {
  mode: SearchMode;
  request: TrainingRequest;
  fold: Int32Array;
  candidates: SearchCandidate[];
  metric: MetricId;
  direction: 1 | -1;
  /** The one limit: the time, or the trials launched. */
  budgetMs?: number;
  maxTrials?: number;
  /** Find best model only. */
  strategy?: SearchStrategy;
  topK: number;
  seed: number;
  parallelism: number;
}
export type TrialState = 'complete' | 'pruned' | 'failed' | 'abandoned';
export interface Trial {
  candidate: number;
  number: number;
  hyperparameters: Hyperparameters;
  state: TrialState;
  score?: number;
  foldScores: number[];
  fitMs: number;
  error?: string;
}
/** Successive halving of Find best model: the round under way (round r gives each candidate of [pool] 2^r trials),
 * each pool member's trials of the round (abandoned ones not counted), the round each candidate was eliminated in
 * (null while it was not); once [finalists] is set, those are tuned in turn from [next]. */
export interface Halving { round: number; pool: number[]; given: number[]; eliminated: (number | null)[];
  finalists?: number[]; next: number }
/** Where an executor left a search, per candidate: the points drawn from its stream, whether it is exhausted, its
 * trials' first fold scores when they passed them (else null); for Find best model, the halving. */
export interface ExecutorState { drawn: number[]; exhausted: boolean[]; passed: (number[] | null)[][];
  halving?: Halving }
/** The trials and the time of a search, a continued one's included; [resume] lets an executor continue it. */
export interface SearchOutcome { trials: Trial[]; elapsedMs: number; stopped: SearchStop; resume?: ExecutorState }
export interface Leader {
  method: string;
  hyperparameters: Hyperparameters;
  score: number;
  se: number;
  foldScores: number[];
  trials: number;
  chosen: boolean;
  note?: string;
}
/** A continued search adds up the limits of its runs: a record of runs with both kinds of limit holds both. */
export interface SearchRecord {
  mode: SearchMode;
  metric: MetricId;
  budgetSeconds?: number;
  maxTrials?: number;
  strategy?: SearchStrategy;
  seed: number;
  trials: number;
  completed: number;
  pruned: number;
  failed: number;
  elapsedMs: number;
  stopped: SearchStop;
  holdOut?: {rows: number; fraction: number};
  /** Tune: the dimensions searched and the settings pinned. */
  space?: {searched: SearchSpace; pinned: string[]};
  /** Find best model: every method tried, in cost order. */
  methods?: CandidateRecord[];
  leaders: Leader[];
}
export type CandidateOutcome = 'complete' | 'failed' | 'abandoned' | 'unstarted';
/** A method of Find best model; [eliminatedRound] is the round of successive halving it was dropped after. */
export interface CandidateRecord { method: string; trials: number; outcome: CandidateOutcome; note?: string;
  eliminatedRound?: number }
/** What the user asks of a search; without a metric, the task's default; one limit, [budgetSeconds] or [maxTrials]
 * (neither: {@link DEFAULT_BUDGET_SECONDS}); Find best model's [strategy], `equal` by default; Tune searches [space]
 * when given, else the method's whole space. */
export interface SearchChoice { metric?: MetricId; budgetSeconds?: number; methods?: string[]; topK?: number;
  maxTrials?: number; strategy?: SearchStrategy; space?: TuneSpace }
export interface SearchExecutor {
  /** Runs the trials of [plan] until its limit (the time, or the trials launched) or exhaustion,
   * {@link MAX_FAILED_WITHOUT_COMPLETE} failed trials before any completes, or until [progress] is cancelled (between
   * fits); reports `update(percent, text)` after every trial. With [previous] (an outcome of the same plan but its
   * limit), continues it: its trials kept first, the limit counting the new ones only. */
  run(plan: SearchPlan, progress: LoopProgress, previous?: SearchOutcome): Promise<SearchOutcome>;
}

export const SEARCH_CAPTIONS: Record<SearchMode, string> = {tune: 'Tune', best: 'Find best model'};
export const SEARCH_STRATEGIES: SearchStrategy[] = ['equal', 'halving'];
const HALVING_CAPTION = 'Focus on leaders';
const FIND_TITLE = 'Finding best model';
export const SEARCH_METRICS: Record<ModelTask, MetricId[]> = {regression: ['r2', 'rmse', 'mae'],
  classification: ['balancedAccuracy', 'accuracy', 'f1', 'mcc', 'auc', 'logLoss', 'brier']};
const METRIC_DIRECTIONS: Partial<Record<MetricId, -1>> = {mse: -1, rmse: -1, mae: -1, logLoss: -1, brier: -1};
export const DEFAULT_BUDGET_SECONDS = 60;
export const DEFAULT_STRATEGY: SearchStrategy = 'equal';
export const TOP_K = 3;
const HOLD_OUT_MIN_ROWS = 500;
const HOLD_OUT_FRACTION = 0.2;
const HOLD_OUT_ROWS = [100, 2000];
export const PARALLEL_FITS = Math.max(1, Math.min(4, Math.floor((navigator.hardwareConcurrency ?? 2) / 2)));
/** The methods from the cheapest to fit to the most expensive: the first round's order of Find best model. */
export const COST_RANK = ['Linear Regression', 'PLS Regression', 'Softmax', 'Naive Bayes', 'LDA', 'Decision Tree',
  'k-Nearest Neighbors', 'XGBoost', 'Hist Gradient Boosting', 'CatBoost', 'Extra Trees', 'Random Forest', 'SVM',
  'Neural Network', 'TabM', 'RealMLP'];
/** Above this time limit Find best model tries every method by default; at it, or with a trial limit, only the
 * methods that retrain live. */
export const QUICK_BUDGET_SECONDS = 60;
/** Points already tried in a row that exhaust a candidate (every space is finite: drawn floats are rounded). */
export const MAX_DUPLICATES = 20;
/** The significant digits of a drawn float. */
export const FLOAT_DIGITS = 3;
/** Failed trials that stop a search while none has completed (with an equal split, a method's own trials end its
 * share); no more are launched meanwhile. */
export const MAX_FAILED_WITHOUT_COMPLETE = 10;
/** Complete trials before a trial can be pruned, and the folds it is judged on. */
export const PRUNE_MIN_TRIALS = 5;
export const PRUNE_AFTER_FOLDS = 2;
export const NO_PROGRESS: LoopProgress = {canceled: false, update: () => {}};
const SIMPLER = ['low', 'high', 'first'] as const;

export function directionOf(metric: MetricId): 1 | -1 {
  return METRIC_DIRECTIONS[metric] ?? 1;
}

export function defaultMetric(task: ModelTask): MetricId {
  return SEARCH_METRICS[task][0];
}

/** The task-bar title of a search: `Tuning <method>` for Tune, `Finding best model` for Find best model. */
export function searchTitle(mode: SearchMode, method: string): string {
  return mode === 'tune' ? `Tuning ${method}` : FIND_TITLE;
}

/** The caption of the cautious score of [metric] ({@link cautiousScore}): `<Metric> − SE`, or `<Metric> + SE` for a
 * minimized one. */
export function cautiousCaption(metric: MetricId): string {
  return `${METRIC_LABELS[metric]} ${directionOf(metric) === 1 ? '−' : '+'} SE`;
}

/** Whether a method gives class probabilities: always, only through its probability input, or never. */
export type MethodProbabilities = 'always' | 'input' | 'none';

export function probabilitiesOf(engine: Engine): MethodProbabilities {
  return probabilityInput(engine) !== undefined ? 'input' : givesProbabilities(engine) ? 'always' : 'none';
}

/** Why a method of [probabilities] cannot take part in a search: `predictProbability` when it gives none while
 * Predict probability is on, `metric` when [metric] needs class probabilities it does not give then; null when it
 * can. The Find block and `Forge:trainModel` word each reason their own way. */
export function exclusionOf(probabilities: MethodProbabilities, predictProbability: boolean, metric: MetricId):
  'predictProbability' | 'metric' | null {
  if (predictProbability && probabilities === 'none')
    return 'predictProbability';
  const has = probabilities === 'input' ? predictProbability : probabilities === 'always';
  return isScoreMetric(metric) && !has ? 'metric' : null;
}

/** Whether a limit of [budgetSeconds] (undefined: a trial limit) is quick: no time limit above
 * {@link QUICK_BUDGET_SECONDS}, so Find best model takes only the quick methods by default. */
export function isQuickLimit(budgetSeconds: number | undefined): boolean {
  return (budgetSeconds ?? 0) <= QUICK_BUDGET_SECONDS;
}

/** Whether a method is quick: it retrains live ([isLive]) and its probabilities are not a train input while
 * [predictProbability] is on (it fits several times slower with them). */
export function isQuickMethod(isLive: boolean, probabilities: MethodProbabilities, predictProbability: boolean):
  boolean {
  return isLive && !(probabilities === 'input' && predictProbability);
}

/** The search space [engine] declares, checked: every name a hyperparameter of its train function, int / float with
 * numbers low < high (log needs low > 0, an int log low >= 1), bounds inside the input's min / max when declared, a
 * categorical with choices inside the input's when declared, `simpler` as allowed; the ml.tune shorthands ([low, high]
 * -> int when both whole else float; a list -> categorical). Undefined without a space; a malformed one throws.
 * [properties] are the engine's visible hyperparameters, read from it unless given. */
export function searchSpaceOf(engine: Engine, properties: DG.Property[] = visibleHyperparametersOf(engine)):
  SearchSpace | undefined {
  const text = searchSpaceText(engine);
  if (text === undefined)
    return undefined;
  const fail = (detail: string) => new ForgeError(`The search space of ${engine.name} ${detail}.`);
  let parsed: unknown;
  try {
    parsed = JSON.parse(text);
  } catch {
    throw fail('is not JSON');
  }
  if (!isRecord(parsed))
    throw fail('is not an object');
  return Object.entries(parsed).map(([name, value]) =>
    dimensionOf(name, value, properties.find((p) => p.name === name), fail));
}

function dimensionOf(name: string, value: unknown, property: DG.Property | undefined,
  fail: (detail: string) => ForgeError): SpaceDimension {
  if (property === undefined)
    throw fail(`names no setting ${name}`);
  const spec = Array.isArray(value) ? shorthandOf(value) : value;
  if (!isRecord(spec))
    throw fail(`gives ${name} no dimension`);
  const {type, low, high, log, step, choices} = spec;
  const simpler = SIMPLER.find((s) => s === spec.simpler);
  if (spec.simpler !== undefined && (simpler === undefined || (simpler === 'first') !== (type === 'categorical')))
    throw fail(`gives ${name} the simpler end ${String(spec.simpler)}`);
  const extra = simpler === undefined ? {} : {simpler};
  if (type === 'categorical') {
    if (!Array.isArray(choices) || choices.length === 0 || !choices.every(isChoice))
      throw fail(`gives ${name} no choices`);
    const allowed = property.choices ?? [];
    const outside = choices.find((c) => allowed.length > 0 && !allowed.includes(String(c)));
    if (outside !== undefined)
      throw fail(`gives ${name} the choice ${outside}, which is not one of: ${allowed.join(', ')}`);
    return {name, type, choices, ...extra};
  }
  if (type !== 'int' && type !== 'float')
    throw fail(`gives ${name} the type ${String(type)}; the types are int, float, categorical`);
  const isLog = log === true;
  const problem = rangeProblem({name, type, log: isLog}, property, typeof low === 'number' ? low : null,
    typeof high === 'number' ? high : null);
  if (problem !== null || typeof low !== 'number' || typeof high !== 'number')
    throw fail(`gives ${name} an invalid range: ${problem}`);
  if (step !== undefined && (typeof step !== 'number' || !(step > 0)))
    throw fail(`gives ${name} a step that is not a positive number`);
  return {name, type, low, high, ...(isLog ? {log: true} : {}), ...(step === undefined ? {} : {step}), ...extra};
}

function shorthandOf(value: unknown[]): {[key: string]: unknown} {
  const [low, high] = value;
  if (value.length === 2 && typeof low === 'number' && typeof high === 'number')
    return {type: Number.isInteger(low) && Number.isInteger(high) ? 'int' : 'float', low, high};
  return {type: 'categorical', choices: value};
}

/** The points of [space] drawn for a candidate, one after another from one stream (`mulberry32(seed ^ fnv1a(name))`),
 * one draw per dimension in space order (a float rounded, see `sampled`); the other settings keep their values of
 * [start]. */
export function* drawnPoints(space: SearchSpace, start: Hyperparameters, seed: number,
  name: string): Generator<Hyperparameters, never> {
  const random = mulberry32(seed ^ fnv1a(FNV_OFFSET, new TextEncoder().encode(name)));
  while (true) {
    const point: Hyperparameters = {...start};
    for (const dim of space)
      point[dim.name] = sampled(dim, random());
    yield point;
  }
}

/** [point] with each dimension of [space] whose `meta.mlenabled` rule ([rules]) fails for it at its value of
 * [defaults] (a dimension without a default keeps its value), so that points differing only in settings the method
 * ignores are one. */
function activePoint(space: SearchSpace, point: Hyperparameters, rules: EnabledRules, defaults: Hyperparameters,
  task: ModelTask): Hyperparameters {
  const active = {...point};
  for (const {name} of space) {
    if (defaults[name] !== undefined && !isEnabled(rules, name, point, task))
      active[name] = defaults[name];
  }
  return active;
}

/** The points of [candidate] for a [task] in turn, their inactive dimensions at the method's defaults
 * ({@link activePoint}): the start first when it lies inside the space ({@link isInsideSpace}), then the
 * {@link drawnPoints}; without a space, the start alone. */
export function* trialPoints(candidate: SearchCandidate, seed: number, task: ModelTask): Generator<Hyperparameters> {
  const {engine, space, start} = candidate;
  if (space === undefined) {
    yield {...start};
    return;
  }
  const rules = enabledRules(engine) ?? {};
  const defaults = defaultHyperparameters(engine);
  if (isInsideSpace(space, start, rules, task))
    yield activePoint(space, start, rules, defaults, task);
  for (const point of drawnPoints(space, start, seed, engine.name))
    yield activePoint(space, point, rules, defaults, task);
}

/** Whether each value of [point] the method uses ([rules], [task]) lies inside its dimension of [space]: a choice of a
 * categorical, a number within the bounds (whole for an int) otherwise. */
export function isInsideSpace(space: SearchSpace, point: Hyperparameters, rules: EnabledRules,
  task: ModelTask): boolean {
  return space.every((dim) => {
    const v = point[dim.name];
    if (!isEnabled(rules, dim.name, point, task))
      return true;
    if (dim.type === 'categorical')
      return dim.choices!.some((c) => String(c) === String(v));
    return typeof v === 'number' && v >= dim.low! && v <= dim.high! && (dim.type === 'float' || Number.isInteger(v));
  });
}

/** [space] narrowed by the user's [tune]: the pinned dimensions left out, and those whose `meta.mlenabled` rule
 * ([rules]) fails with the fixed [values] and the [task] whatever the tuned ones are; the others within the user's
 * bounds or choices, their type, scale, step and simpler end kept. */
export function userSpace(space: SearchSpace, tune: TuneSpace, rules: EnabledRules, values: Hyperparameters,
  task: ModelTask): SearchSpace {
  const tuned = space.filter((dim) => !tune.pinned.includes(dim.name));
  const free = new Set(tuned.map((dim) => dim.name));
  return tuned.filter((dim) => isEnabled(rules, dim.name, values, task, free)).map((dim) => {
    const range = tune.ranges[dim.name];
    if (range === undefined)
      return dim;
    if (dim.type === 'categorical')
      return range.choices === undefined ? dim : {...dim, choices: chosenOf(dim, range.choices)};
    return {...dim, low: range.low ?? dim.low, high: range.high ?? dim.high};
  });
}

/** The choices of the categorical [dim] among [values] (compared as text), in the dimension's order. */
export function chosenOf(dim: SpaceDimension, values: unknown[]): (string | number | boolean)[] {
  return dim.choices!.filter((c) => values.some((x) => String(x) === String(c)));
}

export const FROM_BELOW_TO = '“From” must be less than “to”';

/** Why [low] - [high] cannot bound the tuned [dim] of the setting [property]; null when it can: two numbers, whole for
 * an int, low below high, inside the setting's min and max, low above 0 (at least 1 for an int) on a log scale. */
export function rangeProblem(dim: SpaceDimension, property: DG.Property, low: number | null,
  high: number | null): string | null {
  if (low === null || high === null || Number.isNaN(low) || Number.isNaN(high))
    return 'Enter both bounds';
  if (dim.type === 'int' && !(Number.isInteger(low) && Number.isInteger(high)))
    return 'The bounds must be whole numbers';
  if (!(low < high))
    return FROM_BELOW_TO;
  const min: unknown = property.min;
  const max: unknown = property.max;
  if (typeof min === 'number' && low < min)
    return `“From” must be at least ${min}`;
  if (typeof max === 'number' && high > max)
    return `“To” must be at most ${max}`;
  if (dim.log && !(dim.type === 'int' ? low >= 1 : low > 0))
    return dim.type === 'int' ? '“From” must be at least 1' : '“From” must be above 0';
  return null;
}

/** The value of [dim] at [u] in 0..1: a float rounded to {@link FLOAT_DIGITS} significant digits within its bounds, so
 * the method trains on the value the form shows. */
function sampled(dim: SpaceDimension, u: number): number | string | boolean {
  if (dim.type === 'categorical') {
    const choices = dim.choices!;
    return choices[Math.min(choices.length - 1, Math.floor(u * choices.length))];
  }
  const low = dim.low!;
  const high = dim.high!;
  const logUniform = () => Math.exp(Math.log(low) + u * (Math.log(high) - Math.log(low)));
  if (dim.type === 'float') {
    const v = Number((dim.log ? logUniform() : low + u * (high - low)).toPrecision(FLOAT_DIGITS));
    return Math.min(high, Math.max(low, v));
  }
  const v = dim.log ? Math.round(logUniform()) : low + Math.floor(u * (high - low + 1));
  return dim.step === undefined ? v : Math.min(high, Math.round((v - low) / dim.step) * dim.step + low);
}

/** With at least {@link HOLD_OUT_MIN_ROWS} rows, `clamp(round(0.2 n), 100, 2000)` of [request]'s n rows held out before
 * a search: the first ones of the rows shuffled by `seed + 1`, for a classification each class's share of them (the
 * largest remainders round up, so the total is exact); the search request holds the other rows, with the whole
 * table's preparation. Below that, the request as it is and no hold out. */
export function holdOutSplit(request: TrainingRequest, seed: number): {search: TrainingRequest; holdOut?: HoldOut} {
  const rows = request.target.length;
  if (rows < HOLD_OUT_MIN_ROWS)
    return {search: request};
  const [fewest, most] = HOLD_OUT_ROWS;
  const testRows = Math.min(most, Math.max(fewest, Math.round(HOLD_OUT_FRACTION * rows)));
  const {prepared} = request;
  const order = shuffledRows(rows, seed + 1);
  const isTest = DG.BitSet.create(rows);
  const groups = taskOf(prepared.target) === 'classification' ? classGroups(order, prepared.target) : [[...order]];
  const quotas = groupQuotas(groups.map((g) => g.length), testRows);
  for (const [i, group] of groups.entries()) {
    for (const r of group.slice(0, quotas[i]))
      isTest.set(r, true, false);
  }
  const isSearch = isTest.clone().invert();
  const search: TrainingRequest = {...request, features: request.features.map((c) => c.clone(isSearch)),
    target: request.target.clone(isSearch), prepared: {...prepared,
      columns: prepared.columns.map((c) => c.clone(isSearch)), target: prepared.target.clone(isSearch)},
    keptRows: maskedKeptRows(request.keptRows, isSearch)};
  const holdOut: HoldOut = {columns: prepared.columns.map((c) => c.clone(isTest)),
    target: prepared.target.clone(isTest), rows: isTest.trueCount, fraction: isTest.trueCount / rows};
  return {search, holdOut};
}

/** The table's rows of [mask] over the request's rows: the i-th kept row of the table stays when [mask] holds i. */
function maskedKeptRows(kept: DG.BitSet | null, mask: DG.BitSet): DG.BitSet {
  if (kept === null)
    return mask.clone();
  const result = DG.BitSet.create(kept.length);
  for (const [i, r] of kept.getSelectedIndexes().entries()) {
    if (mask.get(i))
      result.set(r, true, false);
  }
  return result;
}

/** The one candidate of Tune: [engine] with [start] and its [space] (read from the engine unless given), narrowed by
 * the user's [tune] for a [task] when given ({@link userSpace}); without a space throws. */
export function tuneCandidate(engine: Engine, start: Hyperparameters, tune?: {space: TuneSpace; task: ModelTask},
  space: SearchSpace | undefined = searchSpaceOf(engine)): SearchCandidate {
  if (space === undefined)
    throw new ForgeError(`${engine.name} declares no search space, so it cannot be tuned.`);
  if (tune === undefined)
    return {engine, space, start: {...start}, costRank: 0};
  return {engine, space: userSpace(space, tune.space, enabledRules(engine) ?? {}, start, tune.task),
    start: {...start}, costRank: 0, pinned: [...tune.space.pinned]};
}

/** The cost rank of [engine]: its place in {@link COST_RANK}, or after it in [discovery] order. */
export function costRankOf(engine: Engine, discovery: Engine[]): number {
  const rank = COST_RANK.indexOf(engine.name);
  return rank >= 0 ? rank : COST_RANK.length + discovery.indexOf(engine);
}

/** The methods of [applicable] that run in the browser, cheapest first ({@link costRankOf}). */
export function searchableEngines(applicable: Engine[], discovery: Engine[]): Engine[] {
  return applicable.filter((e) => !isServerEngine(e))
    .sort((a, b) => costRankOf(a, discovery) - costRankOf(b, discovery));
}

/** The quick ones of [engines] ({@link isQuickMethod}), whether each retrains live on [columns] / [target] asked at
 * once. */
export async function quickEngines(engines: Engine[], columns: DG.Column[], target: DG.Column,
  predictProbability: boolean): Promise<Engine[]> {
  const live = await Promise.all(engines.map((e) => retrainsLive(e, columns, target)));
  return engines.filter((e, i) => isQuickMethod(live[i], probabilitiesOf(e), predictProbability));
}

/** The candidates of Find best model: each of [engines] with its space (a malformed one becomes a note and the
 * defaults only), its defaults with the probability input of [predictProbability] as the start and its cost rank,
 * cheapest first; none throws. */
export function bestCandidates(engines: Engine[], discovery: Engine[], predictProbability: boolean):
  SearchCandidate[] {
  if (engines.length === 0)
    throw new ForgeError('No method can take part in the search. Check the features and the target.');
  return engines.map((engine): SearchCandidate => {
    const start = withProbabilityInput(engine, defaultHyperparameters(engine), predictProbability);
    const costRank = costRankOf(engine, discovery);
    try {
      return {engine, space: searchSpaceOf(engine), start, costRank};
    } catch (e) {
      return {engine, start, costRank, note: errorMessage(e)};
    }
  }).sort((a, b) => a.costRank - b.costRank);
}

export const COMBINED_LIMITS = 'The search options \'budgetSeconds\' and \'maxTrials\' cannot be combined: ' +
  'choose one limit.';

/** The plan of a search of [mode] over [request] (the search rows) and [candidates] with the choice's metric (checked
 * against the task), its one limit (both refused; neither: {@link DEFAULT_BUDGET_SECONDS}), topK and, for Find best
 * model, strategy (an equal split by trials needs one per candidate); the split is the request's validation
 * ({@link validationSplit}): its folds, or its one hold out; the parallelism is {@link PARALLEL_FITS} when every
 * candidate is parallel, else 1. */
export function searchPlanOf(mode: SearchMode, request: TrainingRequest, candidates: SearchCandidate[],
  choice: SearchChoice, seed: number): SearchPlan {
  const {target} = request.prepared;
  const task = taskOf(target);
  const metric = choice.metric ?? defaultMetric(task);
  if (!SEARCH_METRICS[task].includes(metric)) {
    throw new ForgeError(`The search metric '${metric}' does not apply to a ${task}. ` +
      `Choose one of: ${SEARCH_METRICS[task].join(', ')}.`);
  }
  const {budgetSeconds, maxTrials} = choice;
  if (budgetSeconds !== undefined && maxTrials !== undefined)
    throw new ForgeError(COMBINED_LIMITS);
  const strategy = mode === 'best' ? choice.strategy ?? DEFAULT_STRATEGY : undefined;
  if (strategy === 'equal' && maxTrials !== undefined && maxTrials < candidates.length)
    throw new ForgeError(`At least one trial per method: ${candidates.length}.`);
  const limit = maxTrials !== undefined ? {maxTrials} : {budgetMs: (budgetSeconds ?? DEFAULT_BUDGET_SECONDS) * 1000};
  return {mode, request, fold: validationSplit({...request, seed}), candidates, metric, direction: directionOf(metric),
    ...limit, ...(strategy === undefined ? {} : {strategy}), topK: choice.topK ?? TOP_K, seed,
    parallelism: candidates.every((c) => isParallel(c.engine)) ? PARALLEL_FITS : 1};
}

/** Each candidate's leader, best first: by the one-standard-error rule among its complete trials, the one of the
 * smallest simplicity within one standard error of its best trial, then the earlier trial. */
export function leadersOf(plan: SearchPlan, outcome: SearchOutcome): Leader[] {
  const leaders: Leader[] = [];
  for (const [k, candidate] of plan.candidates.entries()) {
    const own = outcome.trials.filter((t) => t.candidate === k);
    const complete = own.filter((t) => t.state === 'complete').sort((a, b) => a.number - b.number);
    let best: Trial | undefined;
    for (const t of complete) {
      if (best === undefined || plan.direction * (t.score! - best.score!) > 0)
        best = t;
    }
    if (best === undefined)
      continue;
    const bestScore = best.score!;
    const se = standardError(best.foldScores);
    const simplicity = (t: Trial) => simplicityOf(candidate.space ?? [], t.hyperparameters);
    const within = complete.filter((t) => plan.direction * (t.score! - bestScore) >= -se);
    let chosen = within[0];
    for (const t of within) {
      if (simplicity(t) < simplicity(chosen))
        chosen = t;
    }
    leaders.push({method: candidate.engine.name, hyperparameters: {...chosen.hyperparameters}, score: chosen.score!,
      se: standardError(chosen.foldScores), foldScores: [...chosen.foldScores], trials: own.length, chosen: false,
      ...(candidate.note === undefined ? {} : {note: candidate.note})});
  }
  return leaders.sort((a, b) => plan.direction * (b.score - a.score));
}

function standardError(scores: number[]): number {
  const n = scores.length;
  if (n < 2)
    return 0;
  const mean = scores.reduce((sum, s) => sum + s, 0) / n;
  const variance = scores.reduce((sum, s) => sum + (s - mean) ** 2, 0) / (n - 1);
  return Math.sqrt(variance / n);
}

/** The cautious score of [leader]: its score less its standard error, plus it for a metric minimized ([direction]
 * -1). */
export function cautiousScore(leader: Pick<Leader, 'score' | 'se'>, direction: 1 | -1): number {
  return leader.score - direction * leader.se;
}

/** Whether [a] comes before [b] in the Search chart's order: the better cautious score as a frame stores it (32-bit),
 * a tie by the method's name in code-unit order (capitals first), as the bar chart sorts. */
export function isRankedBefore(a: Leader, b: Leader, direction: 1 | -1): boolean {
  const [x, y] = [a, b].map((l) => direction * Math.fround(cautiousScore(l, direction)));
  return x !== y ? x > y : a.method < b.method;
}

/** The choice across the methods' leaders (Tune's one leader is the one-standard-error rule's within its method): the
 * first in the Search chart's order ({@link isRankedBefore}); marks it `chosen`. */
export function chooseLeader(plan: SearchPlan, leaders: Leader[]): Leader | undefined {
  let chosen: Leader | undefined;
  for (const leader of leaders) {
    if (chosen === undefined || isRankedBefore(leader, chosen, plan.direction))
      chosen = leader;
  }
  if (chosen !== undefined)
    chosen.chosen = true;
  return chosen;
}

/** The sum over the dimensions with a simpler end of the point's distance from it, each in 0..1 (in log space for a
 * log dimension, the choice's index for a categorical); dimensions without one count 0. */
function simplicityOf(space: SearchSpace, values: Hyperparameters): number {
  let sum = 0;
  for (const dim of space) {
    const v = values[dim.name];
    if (dim.simpler === undefined || v === undefined)
      continue;
    if (dim.type === 'categorical') {
      const choices = dim.choices!;
      sum += choices.length < 2 ? 0 : Math.max(0, choices.indexOf(v)) / (choices.length - 1);
      continue;
    }
    const scale = (x: number) => dim.log ? Math.log(Math.max(dim.low!, x)) : x;
    const t = Math.min(1, Math.max(0, (scale(Number(v)) - scale(dim.low!)) / (scale(dim.high!) - scale(dim.low!))));
    sum += dim.simpler === 'high' ? 1 - t : t;
  }
  return sum;
}

export function searchRecordFrom(plan: SearchPlan, outcome: SearchOutcome, leaders: Leader[],
  holdOut?: HoldOut): SearchRecord {
  const count = (state: TrialState) => outcome.trials.filter((t) => t.state === state).length;
  const [candidate] = plan.candidates;
  const space = plan.mode === 'tune' && candidate?.space !== undefined ?
    {space: {searched: candidate.space.map((dim) => ({...dim})), pinned: [...candidate.pinned ?? []]}} : {};
  return {mode: plan.mode, metric: plan.metric,
    ...(plan.budgetMs === undefined ? {} : {budgetSeconds: plan.budgetMs / 1000}),
    ...(plan.maxTrials === undefined ? {} : {maxTrials: plan.maxTrials}),
    ...(plan.strategy === undefined ? {} : {strategy: plan.strategy}), seed: plan.seed,
    trials: outcome.trials.length, completed: count('complete'), pruned: count('pruned'), failed: count('failed'),
    elapsedMs: Math.round(outcome.elapsedMs), stopped: outcome.stopped,
    ...(holdOut === undefined ? {} : {holdOut: {rows: holdOut.rows, fraction: holdOut.fraction}}), ...space,
    ...(plan.mode === 'best' ? {methods: candidateOutcomes(plan, outcome)} : {}), leaders};
}

/** Each candidate of [plan] with its trials of [outcome]: complete when one completed, failed (the first error as the
 * note) when none did and one failed, abandoned when its trials were abandoned only, unstarted without a trial; a
 * malformed-space note is kept when there is no error; the round of successive halving it was eliminated in, if it
 * was. */
export function candidateOutcomes(plan: SearchPlan, outcome: SearchOutcome): CandidateRecord[] {
  const eliminated = outcome.resume?.halving?.eliminated ?? [];
  return plan.candidates.map((candidate, k) => {
    const own = outcome.trials.filter((t) => t.candidate === k);
    const failed = own.find((t) => t.state === 'failed');
    const state: CandidateOutcome = own.some((t) => t.state === 'complete') ? 'complete' :
      failed !== undefined ? 'failed' : own.length === 0 ? 'unstarted' : 'abandoned';
    const note = state === 'failed' ? failed?.error ?? candidate.note : candidate.note;
    const round = eliminated[k] ?? undefined;
    return {method: candidate.engine.name, trials: own.length, outcome: state, ...(note === undefined ? {} : {note}),
      ...(round === undefined ? {} : {eliminatedRound: round})};
  });
}

const CANDIDATE_OUTCOMES: CandidateOutcome[] = ['complete', 'failed', 'abandoned', 'unstarted'];

const SEARCH_MODES: SearchMode[] = ['tune', 'best'];
const SEARCH_STOPS: SearchStop[] = ['budget', 'cancelled', 'exhausted', 'trials', 'failed'];

/** A stored `search` value narrowed to {@link SearchRecord}; undefined for anything else. */
export function searchRecordOf(value: unknown): SearchRecord | undefined {
  if (!isRecord(value) || !Array.isArray(value.leaders))
    return undefined;
  const mode = SEARCH_MODES.find((m) => m === value.mode);
  const metric = METRIC_IDS.find((id) => id === value.metric);
  const stopped = SEARCH_STOPS.find((s) => s === value.stopped);
  const {budgetSeconds, maxTrials, seed, trials, completed, pruned, failed, elapsedMs} = value;
  const strategy = SEARCH_STRATEGIES.find((s) => s === value.strategy);
  const leaders = value.leaders.map(leaderOf).filter((l): l is Leader => l !== undefined);
  const isLimit = (x: unknown) => x === undefined || typeof x === 'number';
  if (mode === undefined || metric === undefined || stopped === undefined || !isLimit(budgetSeconds) ||
    !isLimit(maxTrials) || (budgetSeconds === undefined && maxTrials === undefined) ||
    (value.strategy !== undefined && strategy === undefined) || typeof seed !== 'number' ||
    typeof trials !== 'number' || typeof completed !== 'number' || typeof pruned !== 'number' ||
    typeof failed !== 'number' || typeof elapsedMs !== 'number' || leaders.length !== value.leaders.length)
    return undefined;
  const {holdOut} = value;
  const kept = isRecord(holdOut) && typeof holdOut.rows === 'number' && typeof holdOut.fraction === 'number' ?
    {holdOut: {rows: holdOut.rows, fraction: holdOut.fraction}} : {};
  const space = recordSpaceOf(value.space);
  const methods = Array.isArray(value.methods) ? value.methods.map(candidateRecordOf) : undefined;
  const listed = methods !== undefined && methods.every((m): m is CandidateRecord => m !== undefined) ?
    {methods} : {};
  return {mode, metric, ...(typeof budgetSeconds === 'number' ? {budgetSeconds} : {}),
    ...(typeof maxTrials === 'number' ? {maxTrials} : {}), ...(strategy === undefined ? {} : {strategy}), seed, trials,
    completed, pruned, failed, elapsedMs, stopped, ...kept,
    ...(space === undefined ? {} : {space}), ...listed, leaders};
}

function candidateRecordOf(value: unknown): CandidateRecord | undefined {
  if (!isRecord(value))
    return undefined;
  const {method, trials, note, eliminatedRound} = value;
  const outcome = CANDIDATE_OUTCOMES.find((o) => o === value.outcome);
  if (typeof method !== 'string' || typeof trials !== 'number' || outcome === undefined ||
    (note !== undefined && typeof note !== 'string') ||
    (eliminatedRound !== undefined && typeof eliminatedRound !== 'number'))
    return undefined;
  return {method, trials, outcome, ...(note === undefined ? {} : {note}),
    ...(eliminatedRound === undefined ? {} : {eliminatedRound})};
}

/** A stored record's `space`; undefined when it is missing or malformed. */
function recordSpaceOf(value: unknown): SearchRecord['space'] {
  if (!isRecord(value) || !Array.isArray(value.searched) || !Array.isArray(value.pinned) ||
    !value.pinned.every((name) => typeof name === 'string'))
    return undefined;
  const searched = value.searched.map(storedDimensionOf);
  return searched.every((dim): dim is SpaceDimension => dim !== undefined) ?
    {searched, pinned: [...value.pinned]} : undefined;
}

function storedDimensionOf(value: unknown): SpaceDimension | undefined {
  if (!isRecord(value) || typeof value.name !== 'string')
    return undefined;
  const {name, type, low, high, log, step, choices} = value;
  const simpler = SIMPLER.find((s) => s === value.simpler);
  const optional = (x: unknown, kind: string) => x === undefined || typeof x === kind;
  if ((type !== 'int' && type !== 'float' && type !== 'categorical') || !optional(low, 'number') ||
    !optional(high, 'number') || !optional(log, 'boolean') || !optional(step, 'number') ||
    (value.simpler !== undefined && simpler === undefined) ||
    (choices !== undefined && !(Array.isArray(choices) && choices.every(isChoice))))
    return undefined;
  const dim: SpaceDimension = {name, type};
  if (typeof low === 'number')
    dim.low = low;
  if (typeof high === 'number')
    dim.high = high;
  if (log === true)
    dim.log = true;
  if (typeof step === 'number')
    dim.step = step;
  if (Array.isArray(choices))
    dim.choices = [...choices];
  if (simpler !== undefined)
    dim.simpler = simpler;
  return dim;
}

function leaderOf(value: unknown): Leader | undefined {
  if (!isRecord(value) || !isRecord(value.hyperparameters) || !Array.isArray(value.foldScores))
    return undefined;
  const {method, score, se, trials, chosen, note} = value;
  const foldScores = value.foldScores.filter((s): s is number => typeof s === 'number');
  if (typeof method !== 'string' || typeof score !== 'number' || typeof se !== 'number' ||
    typeof trials !== 'number' || typeof chosen !== 'boolean' || foldScores.length !== value.foldScores.length ||
    (note !== undefined && typeof note !== 'string'))
    return undefined;
  const hyperparameters: Hyperparameters = {};
  for (const [name, v] of Object.entries(value.hyperparameters)) {
    if (typeof v === 'number' || typeof v === 'string' || typeof v === 'boolean')
      hyperparameters[name] = v;
  }
  return {method, hyperparameters, score, se, foldScores, trials, chosen, ...(note === undefined ? {} : {note})};
}

/** The limit of [search]: `limit 5 min` (`<n> s` for seconds that are not whole minutes) or `limit 100 trials` (both
 * of a continued search with either: `limit 1 min + 50 trials`), then `, Focus on leaders` for successive halving. */
export function limitCaption(search: Pick<SearchRecord, 'budgetSeconds' | 'maxTrials' | 'strategy'>): string {
  const {budgetSeconds: s, maxTrials} = search;
  const parts = [...(s === undefined ? [] : [s % 60 === 0 ? `${s / 60} min` : `${s} s`]),
    ...(maxTrials === undefined ? [] : [`${maxTrials} trials`])];
  return `limit ${parts.join(' + ')}${search.strategy === 'halving' ? `, ${HALVING_CAPTION}` : ''}`;
}

/** Runs [plan] on [executor], continuing [previous] when given, then the leaders and the choice over all trials. */
export async function runSearch(plan: SearchPlan, progress: LoopProgress, executor: SearchExecutor,
  previous?: SearchOutcome): Promise<{outcome: SearchOutcome; leaders: Leader[]; chosen?: Leader}> {
  const outcome = await executor.run(plan, progress, previous);
  const leaders = leadersOf(plan, outcome);
  return {outcome, leaders, chosen: chooseLeader(plan, leaders)};
}

/** A search before it runs: its plan over the search rows and the rows held out from them. */
export interface PlannedSearch { plan: SearchPlan; holdOut?: HoldOut }
/** What a search leaves to continue it with: its plan (the budget the sum of its runs'), held out rows, outcome and
 * the fingerprint of its search rows ({@link requestFingerprint}). In memory only: the plan holds the columns. */
export interface SearchState { plan: SearchPlan; holdOut?: HoldOut; outcome: SearchOutcome; fingerprint: string }
/** A finished search: its plan, the held out rows, the trials, the leaders, the chosen one, the record and the state
 * to continue it from. */
export interface SearchResult { plan: SearchPlan; holdOut?: HoldOut; outcome: SearchOutcome; leaders: Leader[];
  chosen?: Leader; record: SearchRecord; state: SearchState }

/** Tune's plan: [engine] searched from [start] on [request] without the held out rows; within the choice's space
 * when it has one. [space] is the engine's, read from it unless given. */
export function tuneSearchPlan(request: TrainingRequest, engine: Engine, start: Hyperparameters,
  choice: SearchChoice, seed: number, space?: SearchSpace): PlannedSearch {
  const {search, holdOut} = holdOutSplit(request, seed);
  const tune = choice.space === undefined ? undefined :
    {space: choice.space, task: taskOf(search.prepared.target)};
  return {plan: searchPlanOf('tune', search, [tuneCandidate(engine, start, tune, space)], choice, seed), holdOut};
}

/** Find best model's plan: [candidates] searched on [request] without the held out rows. */
export function bestSearchPlan(request: TrainingRequest, candidates: SearchCandidate[], choice: SearchChoice,
  seed: number): PlannedSearch {
  const {search, holdOut} = holdOutSplit(request, seed);
  return {plan: searchPlanOf('best', search, candidates, choice, seed), holdOut};
}

/** Runs [planned] on [executor], on top of [previous] when given (its limits added to the recorded ones); the state
 * keeps the request's [fingerprint] ({@link requestFingerprint}), computed unless given. */
export async function runPlanned({plan, holdOut}: PlannedSearch, progress: LoopProgress, executor: SearchExecutor,
  previous?: SearchState, fingerprint?: string): Promise<SearchResult> {
  const {outcome, leaders, chosen} = await runSearch(plan, progress, executor, previous?.outcome);
  const sum = (a?: number, b?: number) => a === undefined ? b : b === undefined ? a : a + b;
  const recorded = previous === undefined ? plan : {...plan, budgetMs: sum(previous.plan.budgetMs, plan.budgetMs),
    maxTrials: sum(previous.plan.maxTrials, plan.maxTrials)};
  return {plan, holdOut, outcome, leaders, chosen, record: searchRecordFrom(recorded, outcome, leaders, holdOut),
    state: {plan: recorded, holdOut, outcome,
      fingerprint: previous?.fingerprint ?? fingerprint ?? requestFingerprint(plan.request)}};
}

/** Tune: {@link tuneSearchPlan} run on [executor]. */
export async function tuneSearch(request: TrainingRequest, engine: Engine, start: Hyperparameters,
  choice: SearchChoice, seed: number, progress: LoopProgress, executor: SearchExecutor,
  space?: SearchSpace): Promise<SearchResult> {
  return runPlanned(tuneSearchPlan(request, engine, start, choice, seed, space), progress, executor);
}

/** Find best model: {@link bestSearchPlan} run on [executor]. */
export async function bestSearch(request: TrainingRequest, candidates: SearchCandidate[], choice: SearchChoice,
  seed: number, progress: LoopProgress, executor: SearchExecutor): Promise<SearchResult> {
  return runPlanned(bestSearchPlan(request, candidates, choice, seed), progress, executor);
}

/** The search rows of [request] as text: the selection's feature and target names, the prepared columns' and
 * target's data ({@link datasetFingerprint}) and the preparation options. */
export function requestFingerprint(request: TrainingRequest): string {
  const {prepared} = request;
  return JSON.stringify({features: request.features.map((c) => c.name), target: request.target.name,
    data: datasetFingerprint(prepared.columns, prepared.target).hash, rows: prepared.target.length,
    options: prepared.options});
}

/** Whether [plan] can continue the search of [previous]: the same mode, strategy, search rows
 * ({@link requestFingerprint}), validation, seed and metric, and the same candidates in order, each with the same
 * space, pinned settings and start values outside its space; [fingerprint] is the plan's request's, computed unless
 * given. */
export function canContinue(previous: SearchState, plan: SearchPlan, fingerprint?: string): boolean {
  const before = previous.plan;
  const candidateKey = ({engine, space, pinned, start}: SearchCandidate) => {
    const searched = new Set((space ?? []).map((d) => d.name));
    const fixed = Object.keys(start).filter((name) => !searched.has(name)).sort().map((name) => [name, start[name]]);
    return JSON.stringify([engine.name, space ?? null, pinned ?? [], fixed]);
  };
  const same = (a: unknown, b: unknown) => JSON.stringify(a) === JSON.stringify(b);
  return before.mode === plan.mode && before.strategy === plan.strategy && before.seed === plan.seed &&
    before.metric === plan.metric &&
    same(before.request.validation, plan.request.validation) &&
    same(before.candidates.map(candidateKey), plan.candidates.map(candidateKey)) &&
    previous.fingerprint === (fingerprint ?? requestFingerprint(plan.request));
}

/** Continues [previous] with [planned] (its budget and trial cap for the new trials), its trials kept and its points
 * not tried again; the leaders, the counts and the choice over all trials. Refuses a plan {@link canContinue}
 * rejects. */
export async function continueSearch(previous: SearchState, planned: PlannedSearch, progress: LoopProgress,
  executor: SearchExecutor): Promise<SearchResult> {
  if (!canContinue(previous, planned.plan))
    throw new ForgeError('The search cannot be continued: the data, the validation or the search settings changed.');
  return runPlanned(planned, progress, executor, previous);
}

/** `No trial completed: <error>` with the first failed trial's error; undefined when no trial failed. */
export function failedSearchMessage(outcome: SearchOutcome): string | undefined {
  const failed = outcome.trials.find((t) => t.error !== undefined);
  return failed === undefined ? undefined : `No trial completed: ${failed.error}`;
}

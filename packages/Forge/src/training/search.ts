import * as DG from 'datagrok-api/dg';
import {defaultHyperparameters, enabledRules, EnabledRules, Engine, Hyperparameters, hyperparametersOf, isChoice,
  isEnabled, isParallel, searchSpaceText} from '../engines/engine';
import {LoopProgress} from '../engines/engine-calls';
import {ForgeError} from '../forge-error';
import {ModelTask} from '../generated/db';
import {METRIC_IDS, MetricId, taskOf} from '../metrics/metrics';
import {isRecord} from '../preparation/preparation-options';
import {fnv1a, FNV_OFFSET} from '../storage/dataset-fingerprint';
import {classGroups, groupQuotas, mulberry32, shuffledRows} from './k-fold';
import {HoldOut, TrainingRequest, validationSplit} from './train-model';

export type SearchMode = 'tune' | 'best';
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
  budgetMs: number;
  maxTrials?: number;
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
export interface SearchOutcome { trials: Trial[]; elapsedMs: number; stopped: SearchStop }
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
export interface SearchRecord {
  mode: SearchMode;
  metric: MetricId;
  budgetSeconds: number;
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
  leaders: Leader[];
}
/** What the user asks of a search; without a metric, the task's default; Tune searches [space] when given, else the
 * method's whole space. */
export interface SearchChoice { metric?: MetricId; budgetSeconds: number; methods?: string[]; topK?: number;
  maxTrials?: number; space?: TuneSpace }
export interface SearchExecutor {
  /** Runs the trials of [plan] until its budget, trial cap or exhaustion, {@link MAX_FAILED_WITHOUT_COMPLETE} failed
   * trials before any completes, or until [progress] is cancelled (between fits); reports `update(percent, text)`
   * after every trial. */
  run(plan: SearchPlan, progress: LoopProgress): Promise<SearchOutcome>;
}

export const SEARCH_CAPTIONS: Record<SearchMode, string> = {tune: 'Tune', best: 'Find best model'};
export const SEARCH_METRICS: Record<ModelTask, MetricId[]> = {regression: ['r2', 'rmse', 'mae'],
  classification: ['balancedAccuracy', 'accuracy', 'f1', 'mcc', 'auc']};
const METRIC_DIRECTIONS: Partial<Record<MetricId, -1>> = {mse: -1, rmse: -1, mae: -1};
export const BUDGETS_SECONDS = [60, 300, 900];
export const DEFAULT_BUDGET_SECONDS = 60;
const TOP_K = 3;
const HOLD_OUT_MIN_ROWS = 500;
const HOLD_OUT_FRACTION = 0.2;
const HOLD_OUT_ROWS = [100, 2000];
export const PARALLEL_FITS = Math.max(1, Math.min(4, Math.floor((navigator.hardwareConcurrency ?? 2) / 2)));
/** Points already tried in a row that exhaust a candidate (every space is finite: drawn floats are rounded). */
export const MAX_DUPLICATES = 20;
/** The significant digits of a drawn float. */
export const FLOAT_DIGITS = 3;
/** Failed trials that stop a search while none has completed; no more are launched meanwhile. */
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

/** The search space [engine] declares, checked: every name a hyperparameter of its train function, int / float with
 * numbers low < high (log needs low > 0, an int log low >= 1), bounds inside the input's min / max when declared, a
 * categorical with choices inside the input's when declared, `simpler` as allowed; the ml.tune shorthands ([low, high]
 * -> int when both whole else float; a list -> categorical). Undefined without a space; a malformed one throws.
 * [properties] are the engine's hyperparameters, read from it unless given. */
export function searchSpaceOf(engine: Engine, properties: DG.Property[] = hyperparametersOf(engine)):
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

/** The plan of a search of [mode] over [request] (the search rows) and [candidates] with the choice's metric (checked
 * against the task), budget, topK and maxTrials; the split is the request's validation ({@link validationSplit}): its
 * folds, or its one hold out; the parallelism is {@link PARALLEL_FITS} when every candidate is parallel, else 1. */
export function searchPlanOf(mode: SearchMode, request: TrainingRequest, candidates: SearchCandidate[],
  choice: SearchChoice, seed: number): SearchPlan {
  const {target} = request.prepared;
  const task = taskOf(target);
  const metric = choice.metric ?? defaultMetric(task);
  if (!SEARCH_METRICS[task].includes(metric)) {
    throw new ForgeError(`The search metric '${metric}' does not apply to a ${task}. ` +
      `Choose one of: ${SEARCH_METRICS[task].join(', ')}.`);
  }
  return {mode, request, fold: validationSplit({...request, seed}), candidates, metric, direction: directionOf(metric),
    budgetMs: choice.budgetSeconds * 1000, maxTrials: choice.maxTrials, topK: choice.topK ?? TOP_K, seed,
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

/** The one-standard-error rule across the methods' leaders: the best leader, or the simpler one within one standard
 * error of it (the lower cost rank, then the smaller simplicity of its point, then the earlier place); marks it
 * `chosen`. */
export function chooseLeader(plan: SearchPlan, leaders: Leader[]): Leader | undefined {
  const best = leaders[0];
  if (best === undefined)
    return undefined;
  const rankOf = (leader: Leader) => {
    const candidate = plan.candidates.find((c) => c.engine.name === leader.method);
    return [candidate?.costRank ?? 0, simplicityOf(candidate?.space ?? [], leader.hyperparameters)];
  };
  let chosen = best;
  let chosenRank = rankOf(best);
  for (const leader of leaders.slice(1)) {
    if (plan.direction * (leader.score - best.score) < -best.se)
      continue;
    const rank = rankOf(leader);
    if (rank[0] < chosenRank[0] || (rank[0] === chosenRank[0] && rank[1] < chosenRank[1])) {
      chosen = leader;
      chosenRank = rank;
    }
  }
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
  return {mode: plan.mode, metric: plan.metric, budgetSeconds: plan.budgetMs / 1000, seed: plan.seed,
    trials: outcome.trials.length, completed: count('complete'), pruned: count('pruned'), failed: count('failed'),
    elapsedMs: Math.round(outcome.elapsedMs), stopped: outcome.stopped,
    ...(holdOut === undefined ? {} : {holdOut: {rows: holdOut.rows, fraction: holdOut.fraction}}), ...space, leaders};
}

const SEARCH_MODES: SearchMode[] = ['tune', 'best'];
const SEARCH_STOPS: SearchStop[] = ['budget', 'cancelled', 'exhausted', 'trials', 'failed'];

/** A stored `search` value narrowed to {@link SearchRecord}; undefined for anything else. */
export function searchRecordOf(value: unknown): SearchRecord | undefined {
  if (!isRecord(value) || !Array.isArray(value.leaders))
    return undefined;
  const mode = SEARCH_MODES.find((m) => m === value.mode);
  const metric = METRIC_IDS.find((id) => id === value.metric);
  const stopped = SEARCH_STOPS.find((s) => s === value.stopped);
  const {budgetSeconds, seed, trials, completed, pruned, failed, elapsedMs} = value;
  const leaders = value.leaders.map(leaderOf).filter((l): l is Leader => l !== undefined);
  if (mode === undefined || metric === undefined || stopped === undefined || typeof budgetSeconds !== 'number' ||
    typeof seed !== 'number' || typeof trials !== 'number' || typeof completed !== 'number' ||
    typeof pruned !== 'number' || typeof failed !== 'number' || typeof elapsedMs !== 'number' ||
    leaders.length !== value.leaders.length)
    return undefined;
  const {holdOut} = value;
  const kept = isRecord(holdOut) && typeof holdOut.rows === 'number' && typeof holdOut.fraction === 'number' ?
    {holdOut: {rows: holdOut.rows, fraction: holdOut.fraction}} : {};
  const space = recordSpaceOf(value.space);
  return {mode, metric, budgetSeconds, seed, trials, completed, pruned, failed, elapsedMs, stopped, ...kept,
    ...(space === undefined ? {} : {space}), leaders};
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

/** `1 minute`, `5 minutes`, else `<n> s`. */
export function budgetCaption(seconds: number): string {
  return seconds === 60 ? '1 minute' : seconds % 60 === 0 ? `${seconds / 60} minutes` : `${seconds} s`;
}

/** Runs [plan] on [executor], then the leaders and the choice. */
export async function runSearch(plan: SearchPlan, progress: LoopProgress, executor: SearchExecutor):
  Promise<{outcome: SearchOutcome; leaders: Leader[]; chosen?: Leader}> {
  const outcome = await executor.run(plan, progress);
  const leaders = leadersOf(plan, outcome);
  return {outcome, leaders, chosen: chooseLeader(plan, leaders)};
}

interface TuneOutcome { plan: SearchPlan; holdOut?: HoldOut; outcome: SearchOutcome; leaders: Leader[];
  chosen?: Leader; record: SearchRecord }

/** Tune: [engine] searched from [start] on [request] without the held out rows, on [executor]; within the choice's
 * space when it has one. [space] is the engine's, read from it unless given. */
export async function tuneSearch(request: TrainingRequest, engine: Engine, start: Hyperparameters,
  choice: SearchChoice, seed: number, progress: LoopProgress, executor: SearchExecutor,
  space?: SearchSpace): Promise<TuneOutcome> {
  const {search, holdOut} = holdOutSplit(request, seed);
  const tune = choice.space === undefined ? undefined :
    {space: choice.space, task: taskOf(search.prepared.target)};
  const plan = searchPlanOf('tune', search, [tuneCandidate(engine, start, tune, space)], choice, seed);
  const {outcome, leaders, chosen} = await runSearch(plan, progress, executor);
  return {plan, holdOut, outcome, leaders, chosen, record: searchRecordFrom(plan, outcome, leaders, holdOut)};
}

/** `No trial completed: <error>` with the first failed trial's error; undefined when no trial failed. */
export function failedSearchMessage(outcome: SearchOutcome): string | undefined {
  const failed = outcome.trials.find((t) => t.error !== undefined);
  return failed === undefined ? undefined : `No trial completed: ${failed.error}`;
}

import * as DG from 'datagrok-api/dg';
import {boolColumn} from '../columns';
import {defaultHyperparameters, defaultValuesOf, Engine, givesProbabilities, Hyperparameters, NUMBER_TYPES,
  visibleHyperparametersOf, withProbabilityInput} from '../engines/engine';
import {methodList} from '../engines/method-list';
import {errorMessage, ForgeError} from '../forge-error';
import {ModelStorageMode, ModelTask, TrainingRunStatus} from '../generated/db';
import {BEST_THRESHOLD_METRICS} from '../metrics/metrics';
import {imputeFunction, imputeSettings, imputeSettingsOf, MISSING_VALUES_MODES, MissingValuesSettings}
  from '../preparation/missing-values';
import {DEFAULT_STEPS, DEFAULT_THRESHOLD, DEFAULT_THRESHOLD_METRIC, defaultOneHot, PreparationSteps, thresholdProblem,
  twoClasses} from '../preparation/pipeline';
import {isRecord} from '../preparation/preparation-options';
import {datasetFingerprint} from '../storage/dataset-fingerprint';
import {datasetRefOf} from '../storage/dataset-ref';
import {modelNameOf, StorageChoice, tagsOf, trainingRunOf} from '../storage/model-fields';
import {saveTrainedModel, TrainedModel} from '../storage/model-store';
import {recordTrainingRun} from '../storage/training-run-store';
import {browserSearch} from './browser-search';
import {defaultFeatures} from './default-features';
import {DEFAULT_SEED, DEFAULT_VALIDATION, TEST_FRACTION_RANGE, ValidationScheme, ValidationSettings} from './k-fold';
import {DEFAULT_BUDGET_SECONDS, failedSearchMessage, NO_PROGRESS, SEARCH_METRICS, SearchChoice, SearchRecord,
  tuneSearch} from './search';
import {checkSelection, hasDataProblems, MAX_SEED, MetricsRecord, prepareTraining, problemMessages,
  SelectionCheck, trainingEngines, TrainingRequest, TrainingResult, TrainingSelection, trainModel}
  from './train-model';

const OPTION_NAMES = ['missingValues', 'neighbors', 'distance', 'oneHot', 'skipUniqueCategories',
  'predictProbability', 'threshold', 'thresholdMetric', 'seed', 'scheme', 'folds', 'testFraction', 'stratified',
  'shuffle', 'search'];
const SCHEMES: ValidationScheme[] = ['kfold', 'holdout'];
const SEARCH_OPTION_NAMES = ['metric', 'budgetSeconds', 'maxTrials'];
const SEARCH_METRIC_IDS = [...new Set([...SEARCH_METRICS.regression, ...SEARCH_METRICS.classification])];
const STORAGE_MODES: ModelStorageMode[] = ['none', 'reference', 'copy'];
const NO_SOURCE = 'The table has no known source, so the model cannot keep a reference to it. ' +
  'Use storage none or copy.';

/** The inputs of `Forge:trainModel`; a missing optional input is null or undefined. */
export interface TrainModelInput {
  table: DG.DataFrame;
  target: string;
  features?: string[] | null;
  method?: string | null;
  hyperparameters?: {[name: string]: unknown} | null;
  options?: {[name: string]: unknown} | null;
  name?: string | null;
  description?: string | null;
  tags?: string | null;
  storage?: string | null;
  showProgress?: boolean | null;
}

export interface TrainModelResult {
  /** The saved model's id; null when no name was given. */
  modelId: string | null;
  runId: string;
  name: string | null;
  method: string;
  task: ModelTask;
  target: string;
  features: string[];
  rowCount: number;
  skippedRows: number;
  seed: number;
  metrics: MetricsRecord;
  /** The search that chose the settings, after `options.search`. */
  search?: SearchRecord;
}

/** `missingValues` is undefined when the option was not given: the method decides it. */
interface TrainingOptions {
  missingValues?: MissingValuesSettings;
  steps: PreparationSteps;
  seed: number;
  validation: ValidationSettings;
  search?: SearchChoice;
}
type ModelSave = Pick<TrainedModel, 'name' | 'description' | 'tags' | 'storage'>;

/** The body of `Forge:trainModel`: the Train view's checks, defaults and validation without the view.
 * Every finished attempt is a `training_run` row with `source: api` (a failed or cancelled one is recorded before its
 * error is rethrown); a refused input or a failed preparation records nothing. With a name, the model is saved;
 * nothing is uploaded unless the storage is `copy`. Nothing is added to the table. */
export async function runTrainModel(input: TrainModelInput): Promise<TrainModelResult> {
  const {table} = input;
  const target = table.col(input.target);
  if (target === null)
    throw new ForgeError(`The table has no column '${input.target}'.`);
  const features = featuresOf(table, target, input.features);
  const {engines, engine: preferred} = trainingEngines();
  const method = input.method || null;
  const chosen = method !== null ? engineNamed(engines, method) : preferred;
  const options = optionsOf(input.options ?? {}, features, target);
  const save = modelSaveOf(input, table);
  const selection: TrainingSelection = {engine: chosen, features, target, hyperparameters: {}, seed: options.seed,
    validation: options.validation, missingValues: options.missingValues, steps: options.steps};
  const engine = engineOf(await checkSelection(selection, engines), chosen, method !== null);
  // A suggested method has its own rules (an explicit keep, for one): they are checked against it.
  if (engine !== chosen)
    engineOf(await checkSelection({...selection, engine}, [engine]), engine, true);
  const {predictProbability} = options.steps;
  if (predictProbability && !givesProbabilities(engine))
    throw new ForgeError(`${engine.name} gives no class probabilities, so Predict probability cannot be used with it.`);
  const hyperparameters = withProbabilityInput(engine, hyperparametersFor(engine, input.hyperparameters ?? {}),
    predictProbability);
  const progressOf = (text: string) => input.showProgress === false ? undefined :
    DG.TaskBarProgressIndicator.create(text, {cancelable: true});
  return await trainAndSave({...selection, engine, hyperparameters}, table.name, save, progressOf, options.search);
}

/** The body of `Forge:listMethods`: {@link methodList}; with a table, **Applicable** and **Suggested** as `trainModel`
 * decides them for the [target] column (by default the table's last column), its default features and options. */
export async function runListMethods(table: DG.DataFrame | null, target: string | null): Promise<DG.DataFrame> {
  if (!table)
    return methodList();
  const targetColumn = target ? table.col(target) : table.columns.toList().pop() ?? null;
  if (targetColumn === null)
    throw new ForgeError(target ? `The table has no column '${target}'.` : 'The table has no columns.');
  const features = defaultFeatures(table, targetColumn);
  const {engines, engine: preferred} = trainingEngines();
  const options = optionsOf({}, features, targetColumn);
  const check = await checkSelection({engine: preferred, features, target: targetColumn, hyperparameters: {},
    ...options}, engines);
  const suggested = engineOf(check, preferred, false);
  const frame = methodList(engines);
  frame.columns.add(boolColumn('Applicable', engines.map((e) => check.engines.includes(e))));
  frame.columns.add(boolColumn('Suggested', engines.map((e) => e === suggested)));
  return frame;
}

/** The method [check] trains with: the [chosen] one, or without one the suggested method; the check's problems are
 * thrown as one error. */
function engineOf(check: SelectionCheck, chosen: Engine, isChosen: boolean): Engine {
  const {problems} = check;
  // As the view: without a chosen method, the suggested one replaces the one the check started with.
  const suggests = !isChosen && !hasDataProblems(problems) && check.engines.length > 0;
  const messages = problemMessages(suggests ? {...problems, method: []} : problems);
  if (messages.length > 0)
    throw new ForgeError(messages.join(' '));
  return suggests ? check.best ?? check.engines[0] : chosen;
}

async function trainAndSave(selection: TrainingSelection, datasetName: string, save: ModelSave | null,
  progressOf: (text: string) => DG.TaskBarProgressIndicator | undefined, search?: SearchChoice):
  Promise<TrainModelResult> {
  const startedOn = Date.now();
  const training = `Training ${selection.engine.name} model`;
  let progress = progressOf(search === undefined ? training : `Tuning ${selection.engine.name}`);
  try {
    const prepared = await prepareTraining(selection);
    const fingerprint = datasetFingerprint(prepared.features, prepared.target);
    let request = prepared;
    if (search !== undefined) {
      request = await tuned(selection, prepared, search, progress);
      progress?.close();
      progress = progressOf(training);
    }
    const run = trainingRunOf({request, datasetName, fingerprint, status: 'completed', source: 'api',
      startedOn: new Date(startedOn).toISOString(), durationMs: 0});
    const record = (status: TrainingRunStatus, metrics?: MetricsRecord, error?: string) =>
      recordTrainingRun({...run, status, metrics, error, duration_ms: Date.now() - startedOn});
    let result: TrainingResult;
    try {
      result = await trainModel(request, progress);
    } catch (e) {
      // The caller gets the training's own error, even when its record cannot be written.
      await (progress?.canceled ? record('cancelled') : record('failed', undefined, errorMessage(e)))
        .catch(() => undefined);
      throw e;
    }
    const runId = await recordTrainingRun({...run, options: result.options, status: 'completed',
      metrics: result.metrics, duration_ms: Date.now() - startedOn});
    const {engine, features, target} = selection;
    const modelId = save === null ? null : await saveTrainedModel({...save, engine, datasetName, result, fingerprint,
      columns: [...features, target], runId});
    return {modelId, runId, name: save?.name ?? null, method: engine.name, task: result.task, target: target.name,
      features: features.map((c) => c.name), rowCount: result.rowCount,
      skippedRows: result.options.missingValues?.skippedRows ?? 0, seed: result.seed, metrics: result.metrics,
      ...(result.search === undefined ? {} : {search: result.search})};
  } finally {
    progress?.close();
  }
}

/** Tune: the selection's method searched from its settings on [request] without the held out rows; the request of
 * the winner's training, with the hold out and the search record. A search without a complete trial records
 * nothing. */
async function tuned(selection: TrainingSelection, request: TrainingRequest, choice: SearchChoice,
  progress?: DG.TaskBarProgressIndicator): Promise<TrainingRequest> {
  const {engine, hyperparameters, seed} = selection;
  const {plan, holdOut, outcome, chosen, record} = await tuneSearch(request, engine, hyperparameters, choice, seed,
    progress ?? NO_PROGRESS, browserSearch);
  if (chosen === undefined) {
    throw new ForgeError(progress?.canceled ? 'The search was cancelled before a trial completed.' :
      failedSearchMessage(outcome) ?? 'No trial completed within the budget. Give the search more time.');
  }
  return {...plan.request, hyperparameters: {...chosen.hyperparameters}, holdOut, search: record};
}

/** The [names] columns of [table]; without names, the view's default features. */
function featuresOf(table: DG.DataFrame, target: DG.Column, names: string[] | null | undefined): DG.Column[] {
  if (names === null || names === undefined)
    return defaultFeatures(table, target);
  const missing = names.filter((name) => table.col(name) === null).map((name) => `'${name}'`);
  if (missing.length === 1)
    throw new ForgeError(`The table has no column ${missing[0]}.`);
  if (missing.length > 1)
    throw new ForgeError(`The table has no columns ${missing.join(', ')}.`);
  return names.map((name) => table.getCol(name));
}

function engineNamed(engines: Engine[], name: string): Engine {
  const engine = engines.find((e) => e.name === name);
  if (engine === undefined) {
    throw new ForgeError(`No method '${name}' is installed. ` +
      `The methods are: ${engines.map((e) => e.name).join(', ')}.`);
  }
  return engine;
}

/** The method's defaults overlaid with [given], each value of the kind and in the range its input declares. */
function hyperparametersFor(engine: Engine, given: {[name: string]: unknown}): Hyperparameters {
  const properties = visibleHyperparametersOf(engine);
  const names = properties.map((p) => p.name);
  const unknown = firstUnknown(given, names);
  if (unknown !== undefined) {
    throw new ForgeError(`The method ${engine.name} has no setting '${unknown}'. ` +
      `Its settings are: ${names.join(', ')}.`);
  }
  return {...defaultHyperparameters(engine),
    ...settingValues(properties, given, (name) => `The setting '${name}' of ${engine.name}`)};
}

/** The preparation, the split seed and the validation of [options], with the Train view's defaults. */
function optionsOf(input: {[name: string]: unknown}, features: DG.Column[], target: DG.Column): TrainingOptions {
  const unknown = firstUnknown(input, [...OPTION_NAMES, 'cutoff']);
  if (unknown !== undefined)
    throw new ForgeError(`Unknown option '${unknown}'. The options are: ${OPTION_NAMES.join(', ')}.`);
  const {cutoff, ...options} = input;
  if (options.threshold === undefined)
    options.threshold = cutoff;
  const given = options.missingValues;
  const modeName = given === undefined ? null : MISSING_VALUES_MODES.find((m) => m === given);
  if (modeName === undefined) {
    throw new ForgeError(`The option 'missingValues' needs ${MISSING_VALUES_MODES.join(', ')} ` +
      '(keep: only for a method that handles missing values).');
  }
  const skipUniqueCategories = boolOption(options, 'skipUniqueCategories', DEFAULT_STEPS.skipUniqueCategories);
  const search = options.search === undefined ? undefined : searchChoiceOf(options.search);
  const predictProbability = boolOption(options, 'predictProbability', DEFAULT_STEPS.predictProbability);
  if (predictProbability && target.matches(DG.COLUMN_TYPE_FILTER.NUMERICAL_NO_DATE_TIME))
    throw new ForgeError('Predict probability needs a text or boolean target.');
  const threshold = numberOption(options, 'threshold', DEFAULT_THRESHOLD);
  const problem = predictProbability && twoClasses(target) !== null ? thresholdProblem(threshold) : null;
  if (problem !== null)
    throw new ForgeError(problem);
  const thresholdMetric = BEST_THRESHOLD_METRICS.find((m) => m === (options.thresholdMetric ??
    DEFAULT_THRESHOLD_METRIC));
  if (thresholdMetric === undefined)
    throw new ForgeError(`The option 'thresholdMetric' needs one of: ${BEST_THRESHOLD_METRICS.join(', ')}.`);
  const steps: PreparationSteps = {oneHot: boolOption(options, 'oneHot', defaultOneHot(features, skipUniqueCategories)),
    skipUniqueCategories, predictProbability, threshold,
    ...(options.threshold === undefined ? {thresholdMetric} : {})};
  const seed = wholeOption(options, 'seed', DEFAULT_SEED, 0, MAX_SEED);
  const scheme = SCHEMES.find((s) => s === (options.scheme ?? DEFAULT_VALIDATION.scheme));
  if (scheme === undefined)
    throw new ForgeError(`The option 'scheme' needs ${SCHEMES.join(' or ')}.`);
  const [fewest, most] = TEST_FRACTION_RANGE;
  const testFraction = numberOption(options, 'testFraction', DEFAULT_VALIDATION.testFraction);
  if (testFraction < fewest || testFraction > most)
    throw new ForgeError(`The option 'testFraction' needs a number from ${fewest} to ${most}.`);
  const validation: ValidationSettings = {scheme, folds: wholeOption(options, 'folds', DEFAULT_VALIDATION.folds, 2),
    testFraction, stratified: boolOption(options, 'stratified', DEFAULT_VALIDATION.stratified),
    shuffle: boolOption(options, 'shuffle', DEFAULT_VALIDATION.shuffle)};
  return {missingValues: modeName === null ? undefined : modeName === 'impute' ? imputeOf(options) : {mode: modeName},
    steps, seed, validation, search};
}

/** The `search` option: a metric of either task (the task's check is the plan's), a whole budget in seconds of at
 * least 1 and an optional whole trial cap. */
function searchChoiceOf(value: unknown): SearchChoice {
  if (!isRecord(value) || firstUnknown(value, SEARCH_OPTION_NAMES) !== undefined)
    throw new ForgeError(`The option 'search' needs an object with ${SEARCH_OPTION_NAMES.join(', ')}.`);
  const metric = value.metric === undefined ? undefined : SEARCH_METRIC_IDS.find((id) => id === value.metric);
  if (value.metric !== undefined && metric === undefined) {
    throw new ForgeError(`The search metric '${String(value.metric)}' is not one of: ` +
      `${SEARCH_METRIC_IDS.join(', ')}.`);
  }
  const budgetSeconds = wholeOption(value, 'budgetSeconds', DEFAULT_BUDGET_SECONDS, 1);
  const maxTrials = value.maxTrials === undefined ? undefined : wholeOption(value, 'maxTrials', 0, 1);
  return {metric, budgetSeconds, maxTrials};
}

/** Impute with the imputer's defaults overlaid with the `neighbors` and `distance` of [options]; without the
 * imputer, the check reports it. */
function imputeOf(options: {[name: string]: unknown}): MissingValuesSettings {
  const func = imputeFunction();
  if (func === undefined)
    return {mode: 'impute'};
  const properties = imputeSettingsOf(func);
  const given = Object.fromEntries(Object.entries(options)
    .filter(([key]) => key === 'neighbors' || key === 'distance'));
  const {neighbors, distance} = {...defaultValuesOf(properties),
    ...settingValues(properties, given, (name) => `The option '${name}'`)};
  return imputeSettings(neighbors, distance);
}

/** What saving needs, checked before any work: null without a name. Storage `reference` needs the table's origin. */
function modelSaveOf(input: TrainModelInput, table: DG.DataFrame): ModelSave | null {
  const mode = STORAGE_MODES.find((m) => m === (input.storage ?? 'none'));
  if (mode === undefined)
    throw new ForgeError(`Unknown storage '${input.storage}'. The storage modes are: ${STORAGE_MODES.join(', ')}.`);
  if (input.name === null || input.name === undefined)
    return null;
  const name = modelNameOf(input.name);
  let storage: StorageChoice = {mode: mode === 'copy' ? 'copy' : 'none'};
  if (mode === 'reference') {
    const ref = datasetRefOf(table);
    if (ref === null)
      throw new ForgeError(NO_SOURCE);
    storage = {mode, ref};
  }
  return {name, description: input.description ?? '', tags: tagsOf(input.tags), storage};
}

/** [given] checked against [properties]: each value of the kind and in the range its property declares;
 * [captionOf] names a setting in the messages. */
function settingValues(properties: DG.Property[], given: {[name: string]: unknown},
  captionOf: (name: string) => string): Hyperparameters {
  const values: Hyperparameters = {};
  for (const [name, value] of Object.entries(given)) {
    const property = properties.find((p) => p.name === name);
    if (property !== undefined)
      values[name] = settingValue(property, value, captionOf(name));
  }
  return values;
}

function settingValue(p: DG.Property, value: unknown, caption: string): number | string | boolean {
  if (NUMBER_TYPES.includes(p.propertyType)) {
    const n = numberValue(value, caption);
    if (p.propertyType === DG.TYPE.INT && !Number.isInteger(n))
      throw new ForgeError(`${caption} needs a whole number.`);
    const min: unknown = p.min;
    const max: unknown = p.max;
    if (typeof min === 'number' && n < min)
      throw new ForgeError(`${caption} must be at least ${min}.`);
    if (typeof max === 'number' && n > max)
      throw new ForgeError(`${caption} must be at most ${max}.`);
    return n;
  }
  if (p.propertyType === DG.TYPE.BOOL)
    return boolValue(value, caption);
  const choices = p.choices ?? [];
  if (typeof value === 'string' && (choices.length === 0 || choices.includes(value)))
    return value;
  throw new ForgeError(choices.length > 0 ? `${caption} needs one of: ${choices.join(', ')}.` :
    `${caption} needs a text.`);
}

/** The first key of [given] that is not one of [names]. */
function firstUnknown(given: {[name: string]: unknown}, names: readonly string[]): string | undefined {
  return Object.keys(given).find((key) => !names.includes(key));
}

function numberValue(value: unknown, caption: string): number {
  if (typeof value !== 'number' || !Number.isFinite(value))
    throw new ForgeError(`${caption} needs a number.`);
  return value;
}

function boolValue(value: unknown, caption: string): boolean {
  if (typeof value !== 'boolean')
    throw new ForgeError(`${caption} needs true or false.`);
  return value;
}

function boolOption(options: {[name: string]: unknown}, name: string, fallback: boolean): boolean {
  return boolValue(options[name] ?? fallback, `The option '${name}'`);
}

function numberOption(options: {[name: string]: unknown}, name: string, fallback: number): number {
  return numberValue(options[name] ?? fallback, `The option '${name}'`);
}

function wholeOption(options: {[name: string]: unknown}, name: string, fallback: number, min: number,
  max?: number): number {
  const value = numberOption(options, name, fallback);
  if (!Number.isInteger(value) || value < min || (max !== undefined && value > max)) {
    throw new ForgeError(`The option '${name}' needs a whole number ` +
      `${max === undefined ? `of at least ${min}` : `from ${min} to ${max}`}.`);
  }
  return value;
}

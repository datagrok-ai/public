import * as DG from 'datagrok-api/dg';
import {applicableEngines, EngineFailure} from '../engines/applicable-engines';
import {selectBestEngine} from '../engines/best-engine';
import {Engine, Hyperparameters} from '../engines/engine';
import {apply, isInteractive, LoopProgress, train, yieldToEventLoop} from '../engines/engine-calls';
import {ForgeError} from '../forge-error';
import {ModelTask} from '../generated/db';
import {aucOf, classificationMetrics, METRIC_IDS, MetricValues, metricsOf} from '../metrics/metrics';
import {missingColumnsOf, MissingValuesSettings, missingValuesProblems, prepareMissingValues}
  from '../preparation/missing-values';
import {PreparationSteps, PreparedFeatures, prepareFeatures, replayPostprocessing} from '../preparation/pipeline';
import {BINARY_CLASSIFICATION, IGNORE_MISSING, IMPUTE_MISSING, isRecord, MissingValuesRecord, PreparationOptions}
  from '../preparation/preparation-options';
import {bigIntProblem} from './default-features';
import {kFold} from './k-fold';

export interface ColumnSchema { name: string; type: string; semType?: string }
export interface TargetSchema extends ColumnSchema { categories?: string[] }
export interface FeaturesSchema { columns: ColumnSchema[] }
export interface MetricsRecord { train: MetricValues; validation: MetricValues; positiveClass?: string }
export type SplittingScheme = 'none' | 'kfold' | 'holdout';
export interface Splitting { scheme: SplittingScheme; folds?: number; trainFraction?: number; isStratified?: boolean }

/** What the user chose: the table's own feature columns, missing values not handled yet. */
export interface TrainingSelection {
  engine: Engine;
  features: DG.Column[];
  target: DG.Column;
  hyperparameters: Hyperparameters;
  seed: number;
  folds: number;
  missingValues: MissingValuesSettings;
  steps: PreparationSteps;
}

/** The data a model is trained on. `features` and `target` are the selection's columns with the missing values
 * handled (the user's own columns where nothing changed): the model's feature list, target and data summary describe
 * them. `prepared` is what the method gets, the preparation steps applied, every step recorded in its `options`. */
export interface TrainingRequest {
  engine: Engine;
  features: DG.Column[];
  target: DG.Column;
  prepared: PreparedFeatures;
  hyperparameters: Hyperparameters;
  seed: number;
  folds: number;
}

/** Predict probability: the method's scores of the train and the out-of-fold predictions, and the target they are
 * measured against, kept to cut them again at another cutoff. */
export interface ProbabilityScores {
  actual: DG.Column;
  positiveClass: string;
  train: DG.Column;
  validation: DG.Column;
  /** AUC-ROC of each split's scores, which no cutoff changes. */
  auc: {train?: number; validation?: number};
}

export interface TrainingSetup {
  task: ModelTask;
  target: TargetSchema;
  features: FeaturesSchema;
  options: PreparationOptions;
  splitting: Splitting;
}

export interface TrainingResult extends TrainingSetup {
  blob: Uint8Array;
  metrics: MetricsRecord;
  seed: number;
  hyperparameters: Hyperparameters;
  rowCount: number;
  scores?: ProbabilityScores;
}

export interface TrainingProblems { target: string[]; features: string[]; missingValues: string[]; method: string[] }

export interface SelectionCheck {
  problems: TrainingProblems;
  /** The methods that can learn from the selection, in discovery order; none while a target or features rule fails. */
  engines: Engine[];
  best: Engine | undefined;
  /** The selection's method is listed, live retraining is on for it, and it says it is interactive on this data. */
  isInteractive: boolean;
  /** Methods whose check threw; they are left out of `engines`. */
  failed: EngineFailure[];
  /** The columns and target the methods were asked about; none while a target or features rule fails. */
  prepared?: PreparedFeatures;
}

const NUMERICAL = DG.COLUMN_TYPE_FILTER.NUMERICAL_NO_DATE_TIME;
const NO_METHOD = 'No method can learn from this selection. Check the features and the target.';
const NO_FEATURE_LEFT = 'Skip unique categories leaves no feature. Check more features.';

/** The task follows the target the method learns (Predict probability makes it a regression); the target schema
 * describes the selection's target. */
export function trainingSetupOf(request: TrainingRequest): TrainingSetup {
  const {target, features, folds, prepared} = request;
  const task: ModelTask = prepared.target.matches(NUMERICAL) ? 'regression' : 'classification';
  const targetSchema: TargetSchema = columnSchemaOf(target);
  if (!target.matches(NUMERICAL))
    targetSchema.categories = target.categories;
  return {
    task,
    target: targetSchema,
    features: {columns: features.map(columnSchemaOf)},
    options: prepared.options,
    splitting: {scheme: 'kfold', folds, isStratified: false},
  };
}

/** The rules, then the methods of [engines] that can learn from the selection, the suggested one and whether the
 * selection's method retrains live; the methods are asked about the prepared columns, as they will get them. */
export async function checkSelection(selection: TrainingSelection, engines: Engine[]): Promise<SelectionCheck> {
  const {engine} = selection;
  // The missing values are handled only when it trains.
  const prepared = prepareFeatures(selection.features, selection.target, selection.steps,
    {preprocessingInfo: [], postprocessingInfo: []});
  const {columns, target} = prepared;
  const problems = ruleProblems(selection, columns);
  if (hasDataProblems(problems))
    return failedCheck(problems);
  const {applicable, failed} = await applicableEngines(engines, columns, target);
  const isListed = applicable.some((e) => e.name === engine.name);
  if (!isListed)
    problems.method.push(applicable.length === 0 ? NO_METHOD : cannotLearn(engine));
  return {problems, engines: applicable, best: selectBestEngine(applicable, columns, target),
    isInteractive: isListed && await retrainsLive(engine, columns, target), failed, prepared};
}

/** A check that lists no method: the selection cannot be trained for [problems]. */
export function failedCheck(problems: Partial<TrainingProblems>): SelectionCheck {
  return {problems: {...noProblems(), ...problems}, engines: [], best: undefined, isInteractive: false, failed: []};
}

function noProblems(): TrainingProblems {
  return {target: [], features: [], missingValues: [], method: []};
}

/** A target or features rule fails: no method is asked. */
export function hasDataProblems(problems: TrainingProblems): boolean {
  return problems.target.length > 0 || problems.features.length > 0;
}

/** Live retraining is on for [engine] and it says it is interactive on this data. */
export async function retrainsLive(engine: Engine, features: DG.Column[], target: DG.Column): Promise<boolean> {
  return engine.isLiveUpdate && await isInteractive(engine, features, target);
}

/** Rules 1-7 and 9: checks of the column list that call no method; rules 1b and 7 look at the [prepared] columns,
 * which only Skip unique categories can leave empty. */
function ruleProblems(selection: TrainingSelection, prepared: DG.Column[]): TrainingProblems {
  const {engine, features, target, folds} = selection;
  const problems = noProblems();
  const targetName = target.name;

  const names = features.map((c) => c.name);
  if (names.length === 0)
    problems.features.push('Choose at least one feature.');
  else if (prepared.length === 0)
    problems.features.push(NO_FEATURE_LEFT);
  if (names.includes(targetName))
    problems.features.push(`The target '${targetName}' is also a feature. Uncheck it in Features.`);
  const rows = target.length - target.stats.missingValueCount;
  if (rows < 2 * folds)
    problems.target.push(`Training needs at least ${2 * folds} rows; the table has ${rows}.`);
  const isNumerical = target.matches(NUMERICAL);
  if (target.type === DG.COLUMN_TYPE.BIG_INT)
    problems.target.push(bigIntProblem(targetName, engine.name, 'choose another target'));
  if (!isNumerical && !hasTwoValues(target))
    problems.target.push(`The target '${targetName}' has only one value; a classifier needs at least two.`);
  if (!isNumerical && !target.isCategorical) {
    problems.target.push(`${engine.name} cannot predict '${targetName}' (type ${target.type}). ` +
      'Choose a numerical, text or boolean target.');
  }
  // Text and yes/no columns are left only while One-hot encoding is off.
  const nonNumerical = prepared.filter((c) => !c.matches(NUMERICAL));
  const text = nonNumerical.filter((c) => c.isCategorical).map((c) => c.name);
  const others = nonNumerical.filter((c) => !c.isCategorical).map((c) => c.name);
  if (text.length > 0) {
    problems.features.push(`${engine.name} needs numerical features. Check One-hot encoding in Preparation, ` +
      `or uncheck: ${text.join(', ')}.`);
  }
  if (others.length > 0)
    problems.features.push(`${engine.name} needs numerical features. Uncheck: ${others.join(', ')}.`);
  for (const col of features.filter((c) => c.type === DG.COLUMN_TYPE.BIG_INT))
    problems.features.push(bigIntProblem(col.name, engine.name, 'uncheck it'));
  problems.missingValues.push(...missingValuesProblems(features, selection.missingValues));
  return problems;
}

function cannotLearn(engine: Engine): string {
  return `${engine.name} cannot learn from this selection. ` +
    'It needs numerical features and a numerical, text or boolean target.';
}

/** Skips the rows with a missing target, handles the features' missing values as chosen, then applies the
 * preparation steps. */
export async function prepareTraining(selection: TrainingSelection): Promise<TrainingRequest> {
  const {engine, hyperparameters, seed, folds, missingValues} = selection;
  const isSkip = missingValues.mode === 'skip';
  const hasFeatureGaps = isSkip && missingColumnsOf(selection.features).length > 0;
  const handled = await prepareMissingValues(selection.features, selection.target, missingValues);
  const {features, target, skippedRows} = handled;
  if (target === undefined || target.length < 2 * folds) {
    throw new ForgeError(`After skipping rows with missing values, ${target?.length ?? 0} rows remain; ` +
      `training needs at least ${2 * folds}.`);
  }
  const preprocessingInfo: string[] = [];
  if (handled.imputedColumns.length > 0)
    preprocessingInfo.push(IMPUTE_MISSING);
  if (isSkip ? hasFeatureGaps : handled.failedRows > 0)
    preprocessingInfo.push(IGNORE_MISSING);
  const impute = missingValues.mode === 'impute' ? missingValues.impute : undefined;
  const record: MissingValuesRecord = impute === undefined ? {mode: missingValues.mode, skippedRows} :
    {mode: 'impute', neighbors: impute.neighbors, distance: impute.distance, skippedRows};
  const prepared = prepareFeatures(features, target, selection.steps,
    {preprocessingInfo, postprocessingInfo: [], missingValues: record});
  return {engine, features, target, prepared, hyperparameters, seed, folds};
}

export async function trainModel(request: TrainingRequest, progress?: LoopProgress): Promise<TrainingResult> {
  const {engine, hyperparameters, seed, folds} = request;
  const {columns, target} = request.prepared;
  const setup = trainingSetupOf(request);
  const rowCount = target.length;
  const steps = folds + 1;

  const fold = kFold(rowCount, folds, seed);
  const foldPredictions: DG.Column[] = [];
  for (let i = 0; i < folds; i++) {
    // Fold copies keep the target's full category list, so every fold model has the final model's classes.
    const isFit = DG.BitSet.create(rowCount, (r) => fold[r] !== i);
    const isHeldOut = isFit.clone().invert();
    await checkCancelled(progress);
    const foldBlob = await train(engine, columns.map((c) => c.clone(isFit)), target.clone(isFit), hyperparameters);
    foldPredictions.push(await apply(engine, columns.map((c) => c.clone(isHeldOut)), foldBlob));
    progress?.update(100 * (i + 1) / steps, `Fold ${i + 1} of ${folds}`);
  }
  const outOfFold = DG.Column.fromType(foldPredictions[0].type, foldPredictions[0].name, rowCount);
  const next = new Int32Array(folds);
  for (let r = 0; r < rowCount; r++)
    outOfFold.set(r, foldPredictions[fold[r]].get(next[fold[r]]++), false);

  await checkCancelled(progress);
  const blob = await train(engine, columns, target, hyperparameters);
  const trainPrediction = await apply(engine, columns, blob);
  progress?.update(100, 'Final model');

  const result = {...setup, blob, seed, hyperparameters: {...hyperparameters}, rowCount};
  const {options} = setup;
  if (options.postprocessingInfo.includes(BINARY_CLASSIFICATION) && options.positiveClass !== undefined) {
    const {positiveClass} = options;
    const actual = request.target;
    const scores: ProbabilityScores = {actual, positiveClass, train: trainPrediction, validation: outOfFold,
      auc: {train: aucOf(actual, trainPrediction, positiveClass), validation: aucOf(actual, outOfFold, positiveClass)}};
    return {...result, metrics: probabilityMetrics(scores, options), scores};
  }
  const categories = setup.target.categories;
  const positiveClass = categories !== undefined && categories.length === 2 ? categories[0] : undefined;
  const metrics: MetricsRecord = {
    train: metricsOf(setup.task, target, trainPrediction, positiveClass),
    validation: metricsOf(setup.task, target, outOfFold, positiveClass),
  };
  if (positiveClass !== undefined)
    metrics.positiveClass = positiveClass;
  return {...result, metrics};
}

/** A Predict probability result cut again at [cutoff], without retraining: the metrics of the labels change, AUC-ROC
 * does not. Any other result is returned as it is. */
export function recutResult(result: TrainingResult, cutoff: number): TrainingResult {
  if (result.scores === undefined)
    return result;
  const options = {...result.options, binaryClassificationThreshold: cutoff};
  return {...result, options, metrics: probabilityMetrics(result.scores, options)};
}

/** The metrics of the labels the scores give at the cutoff of [options], and AUC-ROC of the scores themselves. */
function probabilityMetrics(scores: ProbabilityScores, options: PreparationOptions): MetricsRecord {
  const {actual, positiveClass} = scores;
  const valuesOf = (split: 'train' | 'validation'): MetricValues => {
    const values = classificationMetrics(actual, replayPostprocessing(scores[split], options), positiveClass);
    const auc = scores.auc[split];
    return auc === undefined ? values : {...values, auc};
  };
  return {train: valuesOf('train'), validation: valuesOf('validation'), positiveClass};
}

/** A stored `metrics` value: the known metric ids with numbers; null without a train and a validation object. */
export function metricsRecordOf(value: unknown): MetricsRecord | null {
  if (!isRecord(value))
    return null;
  const train = metricValuesOf(value.train);
  const validation = metricValuesOf(value.validation);
  if (train === null || validation === null)
    return null;
  const record: MetricsRecord = {train, validation};
  if (typeof value.positiveClass === 'string')
    record.positiveClass = value.positiveClass;
  return record;
}

function metricValuesOf(value: unknown): MetricValues | null {
  if (!isRecord(value))
    return null;
  const values: MetricValues = {};
  for (const id of METRIC_IDS) {
    const v = value[id];
    if (typeof v === 'number')
      values[id] = v;
  }
  return values;
}

/** Gives the event loop a turn first, so a click on the progress's cancel is seen before the next fit. */
async function checkCancelled(progress?: LoopProgress): Promise<void> {
  await yieldToEventLoop();
  if (progress?.canceled)
    throw new ForgeError('Training was cancelled.');
}

function columnSchemaOf(col: DG.Column): ColumnSchema {
  const semType = col.semType;
  return semType ? {name: col.name, type: col.type, semType} : {name: col.name, type: col.type};
}

function hasTwoValues(col: DG.Column): boolean {
  let first: string | undefined;
  for (let i = 0; i < col.length; i++) {
    if (col.isNone(i))
      continue;
    const value = String(col.get(i));
    if (first === undefined)
      first = value;
    else if (value !== first)
      return true;
  }
  return false;
}

import * as DG from 'datagrok-api/dg';
import {randomInt} from '@datagrok-libraries/utils/src/random';
import {applicableEngines, EngineFailure} from '../engines/applicable-engines';
import {selectBestEngine} from '../engines/best-engine';
import {defaultMissingValues, Engine, handlesMissingValues, Hyperparameters, isComplete} from '../engines/engine';
import {EngineRegistry} from '../engines/engine-registry';
import {apply, EnginePrediction, isInteractive, LoopProgress, train, yieldToEventLoop} from '../engines/engine-calls';
import {ForgeError} from '../forge-error';
import {ModelTask} from '../generated/db';
import {classesOf, classLabelMetrics, ClassMetrics, ClassProbabilities, Curves, METRIC_IDS, MetricValues, metricsOf,
  probabilityOf, Split, taskOf} from '../metrics/metrics';
import {classMetricsOf, curvesOf, measure, Measured, SCORE_METRIC_IDS} from '../metrics/report';
import {missingColumnsOf, MissingValuesSettings, missingValuesProblems, prepareMissingValues}
  from '../preparation/missing-values';
import {PreparationSteps, PreparedFeatures, prepareFeatures, replayPostprocessing} from '../preparation/pipeline';
import {BINARY_CLASSIFICATION, IGNORE_MISSING, IMPUTE_MISSING, isRecord, MissingValuesRecord, PreparationOptions}
  from '../preparation/preparation-options';
import {bigIntProblem} from './default-features';
import {dataInsights, Insight} from './insights';
import {kFold} from './k-fold';

export interface ColumnSchema { name: string; type: string; semType?: string }
export interface TargetSchema extends ColumnSchema { categories?: string[] }
export interface FeaturesSchema { columns: ColumnSchema[] }
export interface MetricsRecord {
  train: MetricValues;
  validation: MetricValues;
  positiveClass?: string;
  /** A classifier's per-class values, when EDA's report gave them. */
  classes?: ClassMetrics[];
}
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
  /** Undefined when none was chosen: {@link selectedMissingValues} gives the method's default. */
  missingValues?: MissingValuesSettings;
  steps: PreparationSteps;
}

/** The selection's missing-values handling, or the method's default when none was chosen. */
export function selectedMissingValues(selection: TrainingSelection): MissingValuesSettings {
  return selection.missingValues ?? defaultMissingValues(selection.engine);
}

/** The data a model is trained on. `features` and `target` are the selection's columns with the missing values
 * handled (the user's own columns where nothing changed): the model's feature list, target and data summary describe
 * them. `prepared` is what the method gets, the preparation steps applied, every step recorded in its `options`. */
export interface TrainingRequest {
  engine: Engine;
  features: DG.Column[];
  target: DG.Column;
  /** The selection's rows the request holds, or null when no row was skipped. */
  keptRows: DG.BitSet | null;
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
  /** The final model's prediction of the training rows, as the method returned it (with Predict probability, the
   * scores: `scores.train`); {@link trainPredictionOf} gives the labels. */
  prediction: DG.Column;
  metrics: MetricsRecord;
  seed: number;
  hyperparameters: Hyperparameters;
  rowCount: number;
  scores?: ProbabilityScores;
  /** The final model's class probabilities of the training rows, when the method returns them (`proba_<class>`
   * columns after its prediction). */
  probabilities?: ClassProbabilities;
  /** The ROC and PR curves of the final model's positive-class probability (or scores) on the training rows, for two
   * classes with a positive class. */
  curves?: Curves;
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
  /** The data insights of the prepared columns and the selection's target; none while a target or features rule
   * fails. */
  insights: Insight[];
}

/** The method a training starts with when it is installed; the check then suggests one for the data. */
const PREFERRED_ENGINE = 'XGBoost';
export const FOLDS = 5;
export const MAX_SEED = 2 ** 31 - 1;
const NUMERICAL = DG.COLUMN_TYPE_FILTER.NUMERICAL_NO_DATE_TIME;
const NO_METHOD = 'No method can learn from this selection. Check the features and the target.';
const NO_FEATURE_LEFT = 'Skip unique categories leaves no feature. Check more features.';

/** The task follows the target the method learns (Predict probability makes it a regression); the target schema
 * describes the selection's target. */
export function trainingSetupOf(request: TrainingRequest): TrainingSetup {
  const {target, features, folds, prepared} = request;
  const targetSchema: TargetSchema = columnSchemaOf(target);
  if (!target.matches(NUMERICAL))
    targetSchema.categories = target.categories;
  return {
    task: taskOf(prepared.target),
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
  return {problems, engines: applicable, best: selectBestEngine(applicable, selection.features, selection.target),
    isInteractive: isListed && await retrainsLive(engine, columns, target), failed, prepared,
    insights: dataInsights(columns, selection.target)};
}

/** A check that lists no method: the selection cannot be trained for [problems]. */
export function failedCheck(problems: Partial<TrainingProblems>): SelectionCheck {
  return {problems: {...noProblems(), ...problems}, engines: [], best: undefined, isInteractive: false, failed: [],
    insights: []};
}

function noProblems(): TrainingProblems {
  return {target: [], features: [], missingValues: [], method: []};
}

/** A target or features rule fails: no method is asked. */
export function hasDataProblems(problems: TrainingProblems): boolean {
  return problems.target.length > 0 || problems.features.length > 0;
}

/** Every message of [problems]: the target, the features, the missing values, the method. */
export function problemMessages({target, features, missingValues, method}: TrainingProblems): string[] {
  return [...target, ...features, ...missingValues, ...method];
}

/** The complete methods in discovery order, and the one a training starts with. */
export function trainingEngines(): {engines: Engine[]; engine: Engine} {
  const engines = EngineRegistry.discover().filter(isComplete);
  const engine = engines.find((e) => e.name === PREFERRED_ENGINE) ?? engines[0];
  if (engine === undefined)
    throw new ForgeError('No machine learning method is installed. Install the EDA package.');
  return {engines, engine};
}

/** A seed for a training that is given none. */
export function randomSeed(): number {
  return randomInt(MAX_SEED + 1);
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
  const missingValues = selectedMissingValues(selection);
  problems.missingValues.push(...missingValuesProblems(features, missingValues));
  if (missingValues.mode === 'keep' && !handlesMissingValues(engine))
    problems.missingValues.push(needsRowsWithoutGaps(engine));
  return problems;
}

/** Why `keep` cannot train [engine], a method without `meta.mlnativenans`. */
export function needsRowsWithoutGaps(engine: Engine): string {
  return `${engine.name} needs rows without missing values. Choose Skip rows or Impute.`;
}

function cannotLearn(engine: Engine): string {
  return `${engine.name} cannot learn from this selection. ` +
    'It needs numerical features and a numerical, text or boolean target.';
}

/** Skips the rows with a missing target, handles the features' missing values as chosen, then applies the
 * preparation steps. */
export async function prepareTraining(selection: TrainingSelection): Promise<TrainingRequest> {
  const {engine, hyperparameters, seed, folds} = selection;
  const missingValues = selectedMissingValues(selection);
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
  return {engine, features, target, keptRows: handled.keptRows, prepared, hyperparameters, seed, folds};
}

export async function trainModel(request: TrainingRequest, progress?: LoopProgress): Promise<TrainingResult> {
  const {engine, hyperparameters, seed, folds} = request;
  const {columns, target} = request.prepared;
  const setup = trainingSetupOf(request);
  const rowCount = target.length;
  const steps = folds + 1;

  const fold = kFold(rowCount, folds, seed);
  const foldResults: EnginePrediction[] = [];
  for (let i = 0; i < folds; i++) {
    // Fold copies keep the target's full category list, so every fold model has the final model's classes.
    const isFit = DG.BitSet.create(rowCount, (r) => fold[r] !== i);
    const isHeldOut = isFit.clone().invert();
    await checkCancelled(progress);
    const foldBlob = await train(engine, columns.map((c) => c.clone(isFit)), target.clone(isFit), hyperparameters);
    foldResults.push(await apply(engine, columns.map((c) => c.clone(isHeldOut)), foldBlob));
    progress?.update(100 * (i + 1) / steps, `Fold ${i + 1} of ${folds}`);
  }
  const next = new Int32Array(folds);
  const position = new Int32Array(rowCount);
  const first = foldResults[0].prediction;
  const outOfFold = DG.Column.fromType(first.type, first.name, rowCount);
  for (let r = 0; r < rowCount; r++) {
    position[r] = next[fold[r]]++;
    outOfFold.set(r, foldResults[fold[r]].prediction.get(position[r]), false);
  }

  await checkCancelled(progress);
  const blob = await train(engine, columns, target, hyperparameters);
  const {prediction: trainPrediction, probabilities} = await apply(engine, columns, blob);
  progress?.update(100, 'Final model');

  const result = {...setup, blob, prediction: trainPrediction, seed, hyperparameters: {...hyperparameters}, rowCount};
  const {options} = setup;
  if (options.postprocessingInfo.includes(BINARY_CLASSIFICATION) && options.positiveClass !== undefined) {
    const {positiveClass} = options;
    const scores: ProbabilityScores = {actual: request.target, positiveClass, train: trainPrediction,
      validation: outOfFold};
    const measured = (split: Split) => measure('classification', scores.actual,
      replayPostprocessing(scores[split], options), positiveClass, {classes: [positiveClass], columns: [scores[split]]},
      split === 'train');
    const [onTrain, onValidation] = await Promise.all([measured('train'), measured('validation')]);
    const {record, curves} = recordOf(onTrain, onValidation, positiveClass);
    return {...result, metrics: record, scores, curves};
  }
  const classes = setup.target.categories === undefined ? [] : classesOf(request.target);
  const positiveClass = classes.length === 2 ? classes[0] : undefined;
  const validationProbabilities = probabilities === undefined ? undefined :
    outOfFoldProbabilities(foldResults, probabilities, fold, position);
  const [onTrain, onValidation] = await Promise.all([
    measure(setup.task, target, trainPrediction, positiveClass, probabilities, true),
    measure(setup.task, target, outOfFold, positiveClass, validationProbabilities)]);
  const {record, curves} = recordOf(onTrain, onValidation, positiveClass);
  return {...result, metrics: record, probabilities, curves};
}

/** The record of a training's two measurements, with the per-class values and the train split's curves of their
 * reports. */
function recordOf(train: Measured, validation: Measured, positiveClass: string | undefined):
  {record: MetricsRecord; curves?: Curves} {
  const record: MetricsRecord = {train: train.values, validation: validation.values};
  if (positiveClass !== undefined)
    record.positiveClass = positiveClass;
  const classes = classMetricsOf(train.report, validation.report);
  if (classes !== undefined)
    record.classes = classes;
  return {record, curves: curvesOf(train.report)};
}

/** The out-of-fold probabilities of the classes of [final] (the final model's, whose column names they take),
 * gathered by class name from the fold models' results ([fold] and [position] map a row to its fold and its row
 * there); a class a fold model never saw is 0 in its rows; undefined when a fold returned none. */
export function outOfFoldProbabilities(results: EnginePrediction[], final: ClassProbabilities,
  fold: ArrayLike<number>, position: ArrayLike<number>): ClassProbabilities | undefined {
  if (results.some((r) => r.probabilities === undefined))
    return undefined;
  const columns = final.classes.map((cls, k) => {
    const sources = results.map((r) => probabilityOf(r.probabilities!, cls)?.getRawData() ?? null);
    const values = new Float32Array(position.length);
    for (let r = 0; r < position.length; r++) {
      const source = sources[fold[r]];
      values[r] = source === null ? 0 : source[position[r]];
    }
    return DG.Column.fromFloat32Array(final.columns[k].name, values);
  });
  return {classes: final.classes, columns};
}

/** A Predict probability result cut again at [cutoff], without retraining: the metrics of the labels change (overall
 * and per class), AUC-ROC, average precision and the curves do not. Any other result is returned as it is. */
export function recutResult(result: TrainingResult, cutoff: number): TrainingResult {
  const {scores, metrics} = result;
  if (scores === undefined)
    return result;
  const options = {...result.options, binaryClassificationThreshold: cutoff};
  const {actual, positiveClass} = scores;
  const labels = {train: replayPostprocessing(scores.train, options),
    validation: replayPostprocessing(scores.validation, options)};
  const recut: MetricsRecord = {...metrics, positiveClass,
    train: {...scoreValuesOf(metrics.train), ...metricsOf('classification', actual, labels.train, positiveClass)},
    validation: {...scoreValuesOf(metrics.validation),
      ...metricsOf('classification', actual, labels.validation, positiveClass)}};
  if (metrics.classes !== undefined) {
    const train = classLabelMetrics(actual, labels.train);
    const validation = classLabelMetrics(actual, labels.validation);
    recut.classes = metrics.classes.map((c) => ({...c,
      train: {...scoreValuesOf(c.train), ...train.get(c.label)},
      validation: {...scoreValuesOf(c.validation), ...validation.get(c.label)}}));
  }
  return {...result, options, metrics: recut};
}

/** The cutoff-independent values of [values]: AUC-ROC and average precision. */
function scoreValuesOf(values: MetricValues): MetricValues {
  const kept: MetricValues = {};
  for (const id of SCORE_METRIC_IDS) {
    if (values[id] !== undefined)
      kept[id] = values[id];
  }
  return kept;
}

/** The final model's labels of the training rows: with Predict probability, its scores cut at the result's cutoff. */
export function trainPredictionOf(result: TrainingResult): DG.Column {
  return result.scores === undefined ? result.prediction : replayPostprocessing(result.prediction, result.options);
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
  if (Array.isArray(value.classes)) {
    const classes = value.classes.map(classMetricsRecordOf).filter((c): c is ClassMetrics => c !== null);
    if (classes.length > 0)
      record.classes = classes;
  }
  return record;
}

function classMetricsRecordOf(value: unknown): ClassMetrics | null {
  if (!isRecord(value) || typeof value.label !== 'string' || typeof value.rows !== 'number')
    return null;
  const train = metricValuesOf(value.train);
  const validation = metricValuesOf(value.validation);
  return train === null || validation === null ? null : {label: value.label, rows: value.rows, train, validation};
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

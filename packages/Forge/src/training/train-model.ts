import * as DG from 'datagrok-api/dg';
import {Engine, Hyperparameters} from '../engines/engine';
import {apply, isApplicable, train} from '../engines/engine-calls';
import {ForgeError} from '../forge-error';
import {ModelTask} from '../generated/db';
import {MetricValues, metricsOf} from '../metrics/metrics';
import {kFold} from './k-fold';

export interface ColumnSchema { name: string; type: string; semType?: string }
export interface TargetSchema extends ColumnSchema { categories?: string[] }
export interface FeaturesSchema { columns: ColumnSchema[] }
export interface PreparationOptions { preprocessingInfo: string[]; postprocessingInfo: string[] }
export interface MetricsRecord { train: MetricValues; validation: MetricValues; positiveClass?: string }
export type SplittingScheme = 'none' | 'kfold' | 'holdout';
export interface Splitting { scheme: SplittingScheme; folds?: number; trainFraction?: number; isStratified?: boolean }

export interface TrainingRequest {
  engine: Engine;
  features: DG.DataFrame;
  target: DG.Column;
  hyperparameters: Hyperparameters;
  seed: number;
  folds: number;
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
}

export interface TrainingProblems { target: string[]; features: string[] }

export type TrainingProgress = Pick<DG.ProgressIndicator, 'canceled' | 'update'>;

export function trainingSetupOf(request: TrainingRequest): TrainingSetup {
  const {target, features, folds} = request;
  const task: ModelTask = target.matches('numerical') ? 'regression' : 'classification';
  const targetSchema: TargetSchema = columnSchemaOf(target);
  if (task === 'classification')
    targetSchema.categories = target.categories;
  return {
    task,
    target: targetSchema,
    features: {columns: features.columns.toList().map(columnSchemaOf)},
    options: {preprocessingInfo: [], postprocessingInfo: []},
    splitting: {scheme: 'kfold', folds, isStratified: false},
  };
}

export async function trainingProblems(request: TrainingRequest): Promise<TrainingProblems> {
  const {engine, features, target, folds} = request;
  const problems: TrainingProblems = {target: [], features: []};
  const targetName = target.name;

  const names = features.columns.names();
  if (names.length === 0)
    problems.features.push('Choose at least one feature.');
  if (names.includes(targetName))
    problems.features.push(`The target '${targetName}' is also a feature. Uncheck it in Features.`);
  if (target.length < 2 * folds)
    problems.target.push(`Training needs at least ${2 * folds} rows; the table has ${target.length}.`);
  const missing = target.stats.missingValueCount;
  if (missing > 0)
    problems.target.push(`The target '${targetName}' has ${missing} empty values. Remove those rows first.`);
  const isNumerical = target.matches('numerical');
  if (!isNumerical && !hasTwoValues(target))
    problems.target.push(`The target '${targetName}' has only one value; a classifier needs at least two.`);
  if (!isNumerical && !target.isCategorical) {
    problems.target.push(`${engine.name} cannot predict '${targetName}' (type ${target.type}). ` +
      'Choose a numerical, text or boolean target.');
  }
  const nonNumerical = features.columns.toList().filter((c) => !c.matches('numerical')).map((c) => c.name);
  if (nonNumerical.length > 0)
    problems.features.push(`${engine.name} needs numerical features. Uncheck: ${nonNumerical.join(', ')}.`);

  const hasProblems = problems.target.length > 0 || problems.features.length > 0;
  if (!hasProblems && !(await isApplicable(engine, features, target))) {
    problems.features.push(`${engine.name} cannot learn from this selection. ` +
      'It needs numerical features and a numerical, text or boolean target.');
  }
  return problems;
}

export async function checkTrainable(request: TrainingRequest): Promise<void> {
  const problems = await trainingProblems(request);
  const messages = [...problems.target, ...problems.features];
  if (messages.length > 0)
    throw new ForgeError(messages.join(' '));
}

export async function trainModel(request: TrainingRequest, progress?: TrainingProgress): Promise<TrainingResult> {
  const {engine, features, target, hyperparameters, seed, folds} = request;
  const setup = trainingSetupOf(request);
  const rowCount = target.length;
  const steps = folds + 1;

  const fold = kFold(rowCount, folds, seed);
  const foldPredictions: DG.Column[] = [];
  for (let i = 0; i < folds; i++) {
    const isFit = DG.BitSet.create(rowCount, (r) => fold[r] !== i);
    const isHeldOut = isFit.clone().invert();
    checkCancelled(progress);
    const foldBlob = await train(engine, features.clone(isFit), target.clone(isFit), hyperparameters);
    foldPredictions.push(await apply(engine, features.clone(isHeldOut), foldBlob));
    progress?.update(100 * (i + 1) / steps, `Fold ${i + 1} of ${folds}`);
  }
  const outOfFold = DG.Column.fromType(foldPredictions[0].type, foldPredictions[0].name, rowCount);
  const next = new Int32Array(folds);
  for (let r = 0; r < rowCount; r++)
    outOfFold.set(r, foldPredictions[fold[r]].get(next[fold[r]]++), false);

  checkCancelled(progress);
  const blob = await train(engine, features, target, hyperparameters);
  const trainPrediction = await apply(engine, features, blob);
  progress?.update(100, 'Final model');

  const categories = setup.target.categories;
  const positiveClass = categories !== undefined && categories.length === 2 ? categories[0] : undefined;
  const metrics: MetricsRecord = {
    train: metricsOf(setup.task, target, trainPrediction, positiveClass),
    validation: metricsOf(setup.task, target, outOfFold, positiveClass),
  };
  if (positiveClass !== undefined)
    metrics.positiveClass = positiveClass;

  return {...setup, blob, metrics, seed, hyperparameters: {...hyperparameters}, rowCount};
}

function checkCancelled(progress?: TrainingProgress): void {
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

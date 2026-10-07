import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {expect, expectArray, expectFloat} from '@datagrok-libraries/test/src/test';
import {defaultHyperparameters, Engine} from '../engines/engine';
import {EngineRegistry} from '../engines/engine-registry';
import {forgeDb, ModelInsert} from '../generated/db';
import {METRIC_IDS, MetricValues} from '../metrics/metrics';
import {MissingValuesSettings} from '../preparation/missing-values';
import {releaseFrame} from '../preparation/shared-frame';
import {datasetFingerprint} from '../storage/dataset-fingerprint';
import {ModelFields, modelFieldsOf} from '../storage/model-fields';
import {saveModel} from '../storage/model-store';
import {prepareTraining, TrainingRequest, TrainingResult, TrainingSelection, trainModel}
  from '../training/train-model';

export const IRIS = 'System:DemoFiles/iris.csv';
export const MEASUREMENTS = ['Sepal.Length', 'Sepal.Width', 'Petal.Length', 'Petal.Width'];
export const XGBOOST_FIELDS: Pick<ModelInsert, 'engine_name' | 'engine_namespace' | 'engine_kind'> =
  {engine_name: 'XGBoost', engine_namespace: 'Eda', engine_kind: 'function'};
export const IMPUTE: MissingValuesSettings = {mode: 'impute', impute: {neighbors: 4, distance: 'Euclidean'}};

export function engineByName(engines: Engine[], name: string): Engine {
  const engine = engines.find((e) => e.name === name);
  if (engine === undefined)
    throw new Error(`Engine '${name}' is not discovered`);
  return engine;
}

/** A model row without a model file: an XGBoost classifier of Species, no data stored, unless [fields] differ. */
export async function insertModelRow(fields: Pick<ModelInsert, 'name'> & Partial<ModelInsert>): Promise<string> {
  return (await forgeDb.models.insert({...XGBOOST_FIELDS, task: 'classification', target_name: 'Species',
    storage_mode: 'none', ...fields}))[0].id;
}

/** A category's shared fixture, saved by its `before`. */
export function savedFixture<T>(fixture: T | undefined): T {
  if (fixture === undefined)
    throw new Error('The test model was not saved');
  return fixture;
}

let xgboostEngine: Engine | undefined;

export function xgboost(): Engine {
  xgboostEngine ??= engineByName(EngineRegistry.discover(), 'XGBoost');
  return xgboostEngine;
}

export async function openIris(): Promise<DG.DataFrame> {
  const iris = await grok.data.files.openTable(IRIS);
  iris.name = `forge-test-iris-${Date.now()}`;
  return iris;
}

export function columnsOf(table: DG.DataFrame, names: string[]): DG.Column[] {
  return names.map((name) => table.getCol(name));
}

/** A selection of the columns themselves, as the Train view makes it. */
export function selectionOf(features: DG.Column[] | DG.DataFrame, target: DG.Column,
  missingValues: MissingValuesSettings = {mode: 'skip'}): TrainingSelection {
  const engine = xgboost();
  const columns = features instanceof DG.DataFrame ? features.columns.toList() : features;
  return {engine, features: columns, target, hyperparameters: defaultHyperparameters(engine), seed: 42, folds: 5,
    missingValues};
}

/** Runs [action]; returns its result and every frame `DG.DataFrame.fromColumns` built meanwhile over any of
 * [columns]. A frame given back with `releaseFrame` has no columns left. */
export async function framesSharing<T>(columns: DG.Column[], action: () => Promise<T>): Promise<[T, DG.DataFrame[]]> {
  const darts = new Set(columns.map((c) => c.dart));
  const fromColumns = DG.DataFrame.fromColumns;
  const frames: DG.DataFrame[] = [];
  DG.DataFrame.fromColumns = (list: DG.Column[]) => {
    const frame = fromColumns(list);
    if (list.some((c) => darts.has(c.dart)))
      frames.push(frame);
    return frame;
  };
  try {
    return [await action(), frames];
  } finally {
    DG.DataFrame.fromColumns = fromColumns;
  }
}

export function expectReleased(frames: DG.DataFrame[]): void {
  expect(frames.length > 0, true, 'No frame of the user\'s columns was built');
  expectArray(frames.map((f) => f.columns.length), frames.map(() => 0));
}

export async function requestOf(features: DG.Column[] | DG.DataFrame, target: DG.Column): Promise<TrainingRequest> {
  return prepareTraining(selectionOf(features, target));
}

/** Trains a model of [target] by [features] and saves it as a `forge-test-` model; [fields] override the row. */
export async function saveTestModel(features: DG.Column[], target: DG.Column, datasetName: string,
  fields: Partial<ModelFields> = {}): Promise<{id: string; name: string; result: TrainingResult}> {
  const request = await requestOf(features, target);
  const result = await trainModel(request);
  const fingerprint = datasetFingerprint(request.features, request.target);
  releaseFrame(request.features);
  const name = `forge-test-model-${Date.now()}`;
  const id = await saveModel({...modelFieldsOf({name, description: '', tags: [], engine: request.engine, datasetName,
    result, fingerprint}), ...fields}, result.blob);
  return {id, name, result};
}

/** {@link saveTestModel} of an iris classifier of Species by the measurements. */
export async function saveIrisModel(iris: DG.DataFrame): Promise<{id: string; name: string; result: TrainingResult}> {
  return saveTestModel(columnsOf(iris, MEASUREMENTS), iris.getCol('Species'), iris.name);
}

/** Every metric of [expected] within [tolerance] in [actual], and no metric [expected] lacks. */
export function expectMetrics(actual: MetricValues, expected: MetricValues, tolerance: number): void {
  for (const id of METRIC_IDS) {
    const value = expected[id];
    if (value === undefined)
      expect(actual[id] === undefined, true, id);
    else
      expectFloat(actual[id] ?? NaN, value, tolerance, id);
  }
}

/** Values 0..length-1 mapped by [value], empty at [nullRows]. */
export function valuesOf(length: number, value: (i: number) => number, nullRows: number[] = []): (number | null)[] {
  return Array.from({length}, (_, i) => nullRows.includes(i) ? null : value(i));
}

/** The column's raw values, to compare a column before and after an operation. */
export function rawValues(col: DG.Column): number[] {
  return Array.from(col.getRawData().subarray(0, col.length));
}

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {awaitCheck, delay, expect, expectArray, expectFloat} from '@datagrok-libraries/test/src/test';
import {defaultHyperparameters, Engine} from '../engines/engine';
import {EngineRegistry} from '../engines/engine-registry';
import {ForgeError} from '../forge-error';
import {forgeDb, ModelInsert, TrainingRunRow} from '../generated/db';
import {METRIC_IDS, MetricValues} from '../metrics/metrics';
import {MissingValuesSettings} from '../preparation/missing-values';
import {DEFAULT_CUTOFF, PreparationSteps} from '../preparation/pipeline';
import {PreparationOptions} from '../preparation/preparation-options';
import {datasetFingerprint} from '../storage/dataset-fingerprint';
import {ModelFields, modelFieldsOf} from '../storage/model-fields';
import {deleteModel, modelsChanged, saveModel} from '../storage/model-store';
import {prepareTraining, TrainingRequest, TrainingResult, TrainingSelection, trainModel}
  from '../training/train-model';

export const IRIS = 'System:DemoFiles/iris.csv';
export const WINE = 'System:DemoFiles/winequality.csv';
export const MEASUREMENTS = ['Sepal.Length', 'Sepal.Width', 'Petal.Length', 'Petal.Width'];
export const EDA_METHODS = ['Linear Regression', 'Softmax', 'PLS Regression', 'XGBoost', 'SVM', 'Random Forest',
  'Extra Trees', 'Hist Gradient Boosting', 'CatBoost', 'Decision Tree', 'k-Nearest Neighbors', 'Naive Bayes', 'LDA',
  'Neural Network'];
/** The EDA methods that learn a text target (iris Species), and those that learn a numerical one (Petal.Length). */
export const CLASSIFIERS = EDA_METHODS.filter((name) => !name.endsWith('Regression'));
export const REGRESSORS = EDA_METHODS.filter((name) => !['Softmax', 'Naive Bayes', 'LDA'].includes(name));
/** The EDA methods with `meta.mlnativenans`. */
export const NATIVE_NANS = ['XGBoost', 'Random Forest', 'Extra Trees', 'Hist Gradient Boosting', 'CatBoost',
  'Decision Tree'];
export const XGBOOST_FIELDS: Pick<ModelInsert, 'engine_name' | 'engine_namespace' | 'engine_kind'> =
  {engine_name: 'XGBoost', engine_namespace: 'Eda', engine_kind: 'function'};
export const IMPUTE: MissingValuesSettings = {mode: 'impute', impute: {neighbors: 4, distance: 'Euclidean'}};
export const NO_PREPARATION_STEPS: PreparationSteps =
  {oneHot: false, skipUniqueCategories: false, predictProbability: false, cutoff: DEFAULT_CUTOFF};
export const PREDICT_PROBABILITY: PreparationSteps = {...NO_PREPARATION_STEPS, predictProbability: true};
export const ONE_HOT_STEPS: PreparationSteps = {...NO_PREPARATION_STEPS, oneHot: true};
export const SKIP_UNIQUE_STEPS: PreparationSteps = {...NO_PREPARATION_STEPS, skipUniqueCategories: true};
export const NO_OPTIONS: PreparationOptions = {preprocessingInfo: [], postprocessingInfo: []};

/** An `expectExceptionAsync` check: a ForgeError with exactly [text]. */
export const failsWith = (text: string) => (e: unknown) => e instanceof ForgeError && e.message === text;

export function engineByName(engines: Engine[], name: string): Engine {
  const engine = engines.find((e) => e.name === name);
  if (engine === undefined)
    throw new Error(`Engine '${name}' is not discovered`);
  return engine;
}

/** [names] among the [discovered] engines, in discovery order. */
export function inDiscoveryOrder(discovered: Engine[], names: string[]): string[] {
  return discovered.map((e) => e.name).filter((name) => names.includes(name));
}

/** The open balloon that contains [text], or null. */
export function balloonWith(text: string): HTMLElement | null {
  return Array.from(document.querySelectorAll<HTMLElement>('.d4-balloon'))
    .find((b) => (b.textContent ?? '').includes(text)) ?? null;
}

/** Reads [read] every 100 ms until [done] holds of what it read or 5 s pass; the last reading. */
export async function readUntil<T>(read: () => Promise<T> | T, done: (value: T) => boolean): Promise<T> {
  let value = await read();
  for (let i = 0; i < 50 && !done(value); i++) {
    await delay(100);
    value = await read();
  }
  return value;
}

/** The task-bar entry whose text includes [text], or undefined. */
export const taskBarEntry = (text: string) => Array.from(document.querySelectorAll('.d4-task-bar-entry'))
  .find((e) => (e.textContent ?? '').includes(text));

/** Clicks the cancel icon of the task-bar [entry], as the user does. */
export function cancelTaskBarEntry(entry: Element | undefined): void {
  const cancel = entry?.querySelector('.fa-times');
  if (!(cancel instanceof HTMLElement) || !cancel.isConnected)
    throw new Error('No task-bar entry to cancel');
  cancel.click();
}

/** Waits up to [ms] for the task-bar entry that shows [text] and cancels it a moment later. The entry is held from
 * when it shows: a training's text becomes the fold progress after the first fit. */
export async function cancelOnceShown(text: string, ms: number): Promise<void> {
  await awaitCheck(() => taskBarEntry(text) !== undefined, 'No task-bar progress', ms);
  const entry = taskBarEntry(text);
  await delay(100);
  cancelTaskBarEntry(entry);
}

/** The demo table demog of [rows] generated rows, named `forge-test-demog-<time>`. */
export function demoDemog(rows: number): DG.DataFrame {
  const demog = grok.data.demo.demog(rows);
  demog.name = `forge-test-demog-${Date.now()}`;
  return demog;
}

/** Every grid of [grids] (by its pane's name) spans its accordion pane's content within 2 px, and still does once
 * [host], the accordion's container, is made wider. */
export async function expectGridsFillPanes(host: HTMLElement, grids: {[pane: string]: DG.Grid}): Promise<void> {
  const narrower = () => Object.entries(grids).filter(([, grid]) => {
    const pane = grid.root.closest('.d4-accordion-pane-content');
    return pane === null || Math.abs(grid.root.getBoundingClientRect().width - pane.getBoundingClientRect().width) >= 2;
  }).map(([name]) => name).join(', ');
  expect(narrower(), '', 'The grids narrower than their pane');
  host.style.width = `${host.getBoundingClientRect().width + 160}px`;
  await awaitCheck(() => narrower() === '', `Grids narrower than their pane once wider: ${narrower()}`, 5000);
}

/** The catalog's readiness attributes in the Forge view [root]: `<data-forge-catalog-load>/<aria-busy>`. */
export function catalogState(root: Element): string {
  const pane = root.querySelector('.forge-catalog');
  return `${pane?.getAttribute('data-forge-catalog-load')}/${pane?.getAttribute('aria-busy')}`;
}

/** Clicks the catalog's Refresh in the Forge view [root] and reports a change of the models at once. The catalog must
 * be busy then; its readiness attributes ({@link catalogState}) when it next publishes a load and 1.5 s later, and
 * whether the frame [frame] reads was replaced in between. */
export async function refreshThenChange(root: HTMLElement, frame: () => unknown):
  Promise<{first: string; later: string; isRebound: boolean}> {
  const refresh = root.querySelector('.forge-pane-header .fa-sync');
  if (!(refresh instanceof HTMLElement))
    throw new Error('No Refresh icon');
  const state = () => catalogState(root);
  refresh.click();
  modelsChanged.next();
  if (!state().endsWith('/true'))
    throw new Error(`The catalog is not busy during a Refresh and a change: ${state()}`);
  await awaitCheck(() => state().endsWith('/false'), 'The catalog does not publish a load', 5000);
  const first = state();
  const bound = frame();
  await delay(1500);
  return {first, later: state(), isRebound: frame() !== bound};
}

/** The training runs recorded for [tables] (frames or their names). */
export async function runsOf(tables: (DG.DataFrame | string)[]): Promise<TrainingRunRow[]> {
  return await forgeDb.trainingRuns.query().where('dataset_name', '=',
    tables.map((t) => typeof t === 'string' ? t : t.name));
}

export async function deleteRunsOf(tables: (DG.DataFrame | string)[]): Promise<void> {
  await Promise.all((await runsOf(tables)).map((run) => forgeDb.trainingRuns.delete(run.id)));
}

/** Deletes every model named [name] (a test's `finally` after a save that should not have happened). */
export async function deleteModelsNamed(name: string): Promise<void> {
  await Promise.all((await forgeDb.models.query().where('name', '=', name)).map((model) => deleteModel(model.id)));
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

const edaEngines = new Map<string, Engine>();

/** The discovered method [name], discovered once per bundle. */
export function edaEngine(name: string): Engine {
  let engine = edaEngines.get(name);
  if (engine === undefined) {
    engine = engineByName(EngineRegistry.discover(), name);
    edaEngines.set(name, engine);
  }
  return engine;
}

export const xgboost = () => edaEngine('XGBoost');
export const svm = () => edaEngine('SVM');

/** The demo file [path], named `forge-test-<prefix>-<time>`. */
export async function openDemoTable(path: string, prefix: string): Promise<DG.DataFrame> {
  const table = await grok.data.files.openTable(path);
  table.name = `forge-test-${prefix}-${Date.now()}`;
  return table;
}

export async function openIris(): Promise<DG.DataFrame> {
  return openDemoTable(IRIS, 'iris');
}

/** iris opened as Browse > Files opens a file (an unprocessed `OpenServerFile` call), so the platform records its
 * creation script; `grok.data.files.openTable` records none. The table joins the workspace: close it after use. */
export async function openIrisFromFile(): Promise<DG.DataFrame> {
  const call = DG.Func.byName('OpenServerFile').prepare({fullPath: IRIS});
  await call.call(false, undefined, {processed: false});
  const tables: unknown = call.getOutputParamValue();
  const iris: unknown = Array.isArray(tables) ? tables[0] : null;
  if (!(iris instanceof DG.DataFrame))
    throw new Error(`OpenServerFile gave no table for ${IRIS}`);
  return iris;
}

export function columnsOf(table: DG.DataFrame, names: string[]): DG.Column[] {
  return names.map((name) => table.getCol(name));
}

export function columnNamed(columns: DG.Column[], name: string): DG.Column {
  const col = columns.find((c) => c.name === name);
  if (col === undefined)
    throw new Error(`No column '${name}' among ${columns.map((c) => c.name).join(', ')}`);
  return col;
}

export function names(columns: DG.Column[]): string[] {
  return columns.map((c) => c.name);
}

/** A selection of the columns themselves, as the Train view makes it. */
export function selectionOf(features: DG.Column[] | DG.DataFrame, target: DG.Column,
  missingValues: MissingValuesSettings = {mode: 'skip'}, engine: Engine = xgboost()): TrainingSelection {
  const columns = features instanceof DG.DataFrame ? features.columns.toList() : features;
  return {engine, features: columns, target, hyperparameters: defaultHyperparameters(engine), seed: 42, folds: 5,
    missingValues, steps: NO_PREPARATION_STEPS};
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

/** No frame of the user's columns was built: every column the action passed on was a copy. */
export function expectNoFrames(frames: DG.DataFrame[]): void {
  expect(frames.length, 0, 'A frame of the user\'s columns was built');
}

export async function requestOf(features: DG.Column[] | DG.DataFrame, target: DG.Column,
  steps: PreparationSteps = NO_PREPARATION_STEPS, engine?: Engine,
  missingValues: MissingValuesSettings = {mode: 'skip'}): Promise<TrainingRequest> {
  return prepareTraining({...selectionOf(features, target, missingValues, engine), steps});
}

/** iris without setosa: two classes of Species that overlap, so a method's scores are not all 0 or 1. */
export async function twoSpeciesIris(): Promise<DG.DataFrame> {
  const iris = await openIris();
  const species = iris.getCol('Species');
  const twoSpecies = iris.clone(DG.BitSet.create(iris.rowCount, (i) => species.get(i) !== 'setosa'));
  twoSpecies.getCol('Species').compact();
  twoSpecies.name = iris.name;
  return twoSpecies;
}

/** The iris rows 0-119 under iris's name: 50 setosa, 50 versicolor, 20 virginica. */
export async function unbalancedIris(): Promise<DG.DataFrame> {
  const iris = await openIris();
  const unbalanced = iris.clone(DG.BitSet.create(iris.rowCount, (i) => i < 120));
  unbalanced.name = iris.name;
  return unbalanced;
}

/** Trains a model of [target] by [features] and saves it as a `forge-test-` model; [fields] override the row. */
export async function saveTestModel(features: DG.Column[], target: DG.Column, datasetName: string,
  fields: Partial<ModelFields> = {}, steps: PreparationSteps = NO_PREPARATION_STEPS, engine?: Engine,
  missingValues?: MissingValuesSettings): Promise<{id: string; name: string; result: TrainingResult}> {
  const request = await requestOf(features, target, steps, engine, missingValues);
  const result = await trainModel(request);
  const fingerprint = datasetFingerprint(request.features, request.target);
  const name = `forge-test-model-${Date.now()}`;
  const id = await saveModel({...modelFieldsOf({name, description: '', tags: [], engine: request.engine, datasetName,
    result, fingerprint, storage: {mode: 'none'}}), ...fields}, result.blob);
  return {id, name, result};
}

/** {@link saveTestModel} of y by x and 'subject' (all different), which Skip unique categories leaves out. */
export async function saveSkipUniqueModel(): Promise<{id: string; name: string; result: TrainingResult}> {
  const rows = 40;
  const x = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'x', valuesOf(rows, (i) => i / 4));
  const subject = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'subject', Array.from({length: rows}, (_, i) => `s${i}`));
  const y = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'y', valuesOf(rows, (i) => i / 2 + 1));
  return saveTestModel([x, subject], y, `forge-test-skip-unique-${Date.now()}`, {}, SKIP_UNIQUE_STEPS);
}

/** {@link saveTestModel} of an iris classifier of Species by the measurements. */
export async function saveIrisModel(iris: DG.DataFrame, fields: Partial<ModelFields> = {}, engine?: Engine,
  missingValues?: MissingValuesSettings): Promise<{id: string; name: string; result: TrainingResult}> {
  return saveTestModel(columnsOf(iris, MEASUREMENTS), iris.getCol('Species'), iris.name, fields,
    NO_PREPARATION_STEPS, engine, missingValues);
}

/** {@link saveIrisModel} of XGBoost trained under `keep` on a copy of [iris] with gaps at rows 3 and 70. */
export async function saveKeptGapsModel(iris: DG.DataFrame):
  Promise<{id: string; name: string; result: TrainingResult}> {
  const gapped = iris.clone();
  gapped.name = iris.name;
  withGaps(columnsOf(gapped, ['Sepal.Width', 'Petal.Width']));
  return saveIrisModel(gapped, {}, undefined, {mode: 'keep'});
}

/** The rows the gap tests clear: one in the first column, one in the second. */
export const GAP_ROWS = [3, 70];

/** [columns] with the value of the first at row 3 and of the second at row 70 cleared, in place. */
export function withGaps(columns: DG.Column[]): DG.Column[] {
  for (const [i, row] of GAP_ROWS.entries())
    columns[i].set(row, null);
  return columns;
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

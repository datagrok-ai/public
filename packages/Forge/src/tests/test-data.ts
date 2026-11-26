import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {expect, expectArray, expectFloat, waitUntil} from '@datagrok-libraries/test/src/test';
import {defaultHyperparameters, DOCKER_FUNC, Engine, EngineKind} from '../engines/engine';
import {EngineRegistry} from '../engines/engine-registry';
import {ForgeError} from '../forge-error';
import {forgeDb, ModelInsert, TrainingRunRow} from '../generated/db';
import {METRIC_IDS, MetricValues} from '../metrics/metrics';
import {MissingValuesSettings} from '../preparation/missing-values';
import {DEFAULT_THRESHOLD, PreparationSteps} from '../preparation/pipeline';
import {PreparationOptions} from '../preparation/preparation-options';
import {datasetFingerprint} from '../storage/dataset-fingerprint';
import {ModelFields, modelFieldsOf} from '../storage/model-fields';
import {deleteModel, saveModel} from '../storage/model-store';
import {DEFAULT_VALIDATION} from '../training/k-fold';
import {prepareTraining, TrainingRequest, TrainingResult, TrainingSelection, trainModel}
  from '../training/train-model';
import {Gallery} from '../ui/gallery';

export const IRIS = 'System:DemoFiles/iris.csv';
export const WINE = 'System:DemoFiles/winequality.csv';
export const SOLUBILITY = 'System:DemoFiles/chem/solubility_train.csv';
export const MEASUREMENTS = ['Sepal.Length', 'Sepal.Width', 'Petal.Length', 'Petal.Width'];
export const EDA_METHODS = ['Linear Regression', 'Softmax', 'PLS Regression', 'XGBoost', 'SVM', 'Random Forest',
  'Extra Trees', 'Hist Gradient Boosting', 'CatBoost', 'Decision Tree', 'k-Nearest Neighbors', 'Naive Bayes', 'LDA',
  'Neural Network', 'RealMLP', 'TabM'];
/** The EDA methods that learn a text target (iris Species), and those that learn a numerical one (Petal.Length). */
export const CLASSIFIERS = EDA_METHODS.filter((name) => !name.endsWith('Regression'));
export const REGRESSORS = EDA_METHODS.filter((name) => !['Softmax', 'Naive Bayes', 'LDA'].includes(name));
/** The EDA methods with `meta.mlnativenans`. */
export const NATIVE_NANS = ['XGBoost', 'Random Forest', 'Extra Trees', 'Hist Gradient Boosting', 'CatBoost',
  'Decision Tree', 'TabM'];
export const XGBOOST_FIELDS: {engine_name: string; engine_namespace: string; engine_kind: EngineKind} =
  {engine_name: 'XGBoost', engine_namespace: 'Eda', engine_kind: 'function'};
export const IMPUTE: MissingValuesSettings = {mode: 'impute', impute: {neighbors: 4, distance: 'Euclidean'}};
export const NO_PREPARATION_STEPS: PreparationSteps =
  {oneHot: false, skipUniqueCategories: false, predictProbability: false, threshold: DEFAULT_THRESHOLD};
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

/** Waits up to [ms] for the task-bar entry that shows [text] and cancels it. */
export async function cancelOnceShown(text: string, ms: number): Promise<void> {
  await waitUntil(() => taskBarEntry(text) !== undefined, `No task-bar progress ${text}`, ms);
  cancelTaskBarEntry(taskBarEntry(text));
}

/** The demo table demog of [rows] generated rows, named `forge-test-demog-<time>`. */
export function demoDemog(rows: number): DG.DataFrame {
  const demog = grok.data.demo.demog(rows);
  demog.name = `forge-test-demog-${Date.now()}`;
  return demog;
}

/** The tests' own storage of the galleries' modes and sorts: the user's stay as they are, and each category's views
 * start from the defaults. */
const TEST_GALLERY_STORAGE = 'Forge view tests';

const USER_GALLERY_STORAGE = Gallery.storage;

export function useTestGalleryStorage(): void {
  Gallery.storage = TEST_GALLERY_STORAGE;
  grok.userSettings.put(TEST_GALLERY_STORAGE, {});
}

/** Empties the tests' storage and gives the galleries the user's storage back. */
export function clearTestGalleryStorage(): void {
  grok.userSettings.put(TEST_GALLERY_STORAGE, {});
  Gallery.storage = USER_GALLERY_STORAGE;
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

/** A server method whose check answers [answer] (never: `new Promise(() => {})`), counting the calls; an engine stub,
 * not a platform mock. */
export function stuckEngine(name: string, answer: Promise<unknown>, calls = {count: 0}): Engine {
  const {isApplicable, train, apply} = xgboost().functions;
  if (isApplicable === undefined || train === undefined)
    throw new Error('XGBoost has no isApplicable or train');
  const containerFunc = (func: DG.Func): DG.Func => {
    const stub: DG.Func = Object.create(func);
    Object.defineProperty(stub, 'type', {value: DOCKER_FUNC});
    return stub;
  };
  const check = containerFunc(isApplicable);
  check.apply = () => {
    calls.count++;
    return answer;
  };
  return {name, namespace: 'Forge', kind: 'script', functions: {isApplicable: check, train: containerFunc(train),
    apply}, isLiveUpdate: true};
}

/** [base] as a method that never retrains live, so a small table reaches the slow-method paths; with [gate] its
 * [role] calls (train by default, or apply) wait for it, so a test cancels a training or an application at a known
 * point. An engine stub, not a platform mock. */
export function slowEngine(base: Engine, gate?: EngineGate, role: 'train' | 'apply' = 'train'): Engine {
  const func = base.functions[role];
  if (gate === undefined || func === undefined)
    return {...base, isLiveUpdate: false};
  const gated: DG.Func = Object.create(func);
  gated.apply = async (params) => {
    await gate.pass();
    return func.apply(params);
  };
  return {...base, isLiveUpdate: false, functions: {...base.functions, [role]: gated}};
}

/** A promise the test settles when it chooses: a stub's check answer. */
export function deferred<T = void>(): {promise: Promise<T>; resolve: (value: T) => void} {
  const holder: {resolve: (value: T) => void} = {resolve: () => {}};
  const promise = new Promise<T>((resolve) => holder.resolve = resolve);
  return {promise, resolve: (value) => holder.resolve(value)};
}

/** The gate of a {@link slowEngine}: its calls wait until [open]; [started] resolves when the first one reaches it. */
export class EngineGate {
  private readonly opened = deferred();
  private readonly reached = deferred();
  readonly started = this.reached.promise;

  open(): void {
    this.opened.resolve();
  }

  async pass(): Promise<void> {
    this.reached.resolve();
    await this.opened.promise;
  }
}

/** Waits until [gate] holds a call of [work] (or [work] ends first), runs [cancel], opens the gate and returns [work]'s
 * outcome. */
export async function cancelAtGate<T>(gate: EngineGate, work: Promise<T>, cancel: () => unknown): Promise<T> {
  const settled = Promise.allSettled([work]);
  try {
    await Promise.race([gate.started, settled]);
    await cancel();
  } finally {
    gate.open();
  }
  return work;
}

/** [base] with [stub] in place of the method of its name, or added after them. */
export function withEngine(stub: Engine, base: Engine[] = EngineRegistry.discover()): Engine[] {
  return [...base.filter((e) => e.name !== stub.name), stub];
}

/** [rows] rows of two alternating classes and two float features; a test's own table. */
export function twoClassTable(rows: number = 1000): DG.DataFrame {
  const a = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'a', valuesOf(rows, (i) => i % 2 + (i % 7) / 100));
  const b = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'b', valuesOf(rows, (i) => (i % 11) / 10));
  const label = DG.Column.fromStrings('label', Array.from({length: rows}, (_, i) => i % 2 === 0 ? 'even' : 'odd'));
  const table = DG.DataFrame.fromColumns([a, b, label]);
  table.name = `forge-test-two-class-${rows}-${Date.now()}`;
  return table;
}

/** A browser method that declares [semTypes] features (`mlsemtypes`), whose check answers [answer] and counts the
 * [calls]; an engine stub, not a platform mock. */
export function moleculeEngine(name: string, {calls = {count: 0}, answer = true, semTypes = 'Molecule'}:
  {calls?: {count: number}; answer?: boolean; semTypes?: string} = {}): Engine {
  const {isApplicable, train, apply} = xgboost().functions;
  if (isApplicable === undefined || train === undefined)
    throw new Error('XGBoost has no isApplicable or train');
  const trainStub: DG.Func = Object.create(train);
  Object.defineProperty(trainStub, 'options', {value: {mlname: name, mlrole: 'train', mlsemtypes: semTypes}});
  const check: DG.Func = Object.create(isApplicable);
  check.apply = () => {
    calls.count++;
    return Promise.resolve(answer);
  };
  return {name, namespace: 'Forge', kind: 'function', functions: {isApplicable: check, train: trainStub, apply},
    isLiveUpdate: true};
}

/** XGBoost under [name] whose every fit fails with [error]: an engine stub, not a platform mock. */
export function failingEngine(name: string, error: string): Engine {
  const engine = xgboost();
  const {train} = engine.functions;
  if (train === undefined)
    throw new Error('XGBoost has no train');
  const trainStub: DG.Func = Object.create(train);
  Object.defineProperty(trainStub, 'options', {value: {...train.options, mlname: name}});
  trainStub.apply = () => Promise.reject(new Error(error));
  return {...engine, name, functions: {...engine.functions, train: trainStub}};
}

/** A text column of SMILES [values] with the semantic type Molecule. */
export function smilesColumn(name: string, values: string[]): DG.Column {
  const col = DG.Column.fromStrings(name, values);
  col.semType = DG.SEMTYPE.MOLECULE;
  return col;
}

/** [rows] different SMILES: `CO`, `CCO`, `CCCO`, ... */
export function smilesSeries(rows: number): string[] {
  return Array.from({length: rows}, (_, i) => `C${'C'.repeat(i)}O`);
}

/** [rows] rows: `smiles` all different (Molecule), `color` of three values, floats `x` and `y`. */
export function moleculeFrame(rows: number): {smiles: DG.Column; color: DG.Column; x: DG.Column; y: DG.Column} {
  return {
    smiles: smilesColumn('smiles', smilesSeries(rows)),
    color: DG.Column.fromStrings('color', Array.from({length: rows}, (_, i) => ['red', 'green', 'blue'][i % 3])),
    x: DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'x', valuesOf(rows, (i) => i / 2)),
    y: DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'y', valuesOf(rows, (i) => i + 1)),
  };
}

/** The demo file [path], named `forge-test-<prefix>-<time>`. */
export async function openDemoTable(path: string, prefix: string): Promise<DG.DataFrame> {
  const table = await grok.data.files.openTable(path);
  table.name = `forge-test-${prefix}-${Date.now()}`;
  return table;
}

export async function openIris(): Promise<DG.DataFrame> {
  return openDemoTable(IRIS, 'iris');
}

/** solubility_train.csv with `smiles` marked Molecule (not relying on Chem's detectors). */
export async function openSolubility(): Promise<DG.DataFrame> {
  const table = await openDemoTable(SOLUBILITY, 'solubility');
  table.getCol('smiles').semType = DG.SEMTYPE.MOLECULE;
  return table;
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
  return {engine, features: columns, target, hyperparameters: defaultHyperparameters(engine), seed: 42,
    validation: DEFAULT_VALIDATION, missingValues, steps: NO_PREPARATION_STEPS};
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

/* eslint-disable max-lines */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {after, awaitCheck, before, category, delay, expect, expectArray, expectFloat, test, waitUntil,
  waitUntilAsync} from '@datagrok-libraries/test/src/test';
import {recordApplication} from '../apply/application-store';
import {applyAndRecord, DEFAULT_BATCH_SIZE, loadModel} from '../apply/apply-model';
import {exactMapping} from '../apply/column-matching';
import {updateModelInfo} from '../catalog/model-edit';
import {MODEL_TYPE, PANEL_TABLE_MODELS} from '../constants';
import {CHECK_METHOD} from '../engines/check-method';
import {defaultHyperparameters, enabledRules, Engine, hyperparametersOf, isEnabled, isServerEngine,
  visibleHyperparametersOf} from '../engines/engine';
import {EngineRegistry} from '../engines/engine-registry';
import {errorMessage, ForgeError} from '../forge-error';
import {forgeDb, ModelStorageMode} from '../generated/db';
import {METRIC_DESCRIPTIONS, METRIC_LABELS, MetricFinding, MetricId, overfittingOf, roundTo3}
  from '../metrics/metrics';
import {datasetRefOf, storedDatasetRef} from '../storage/dataset-ref';
import {deleteModel, modelsChanged} from '../storage/model-store';
import {CandidateRecord, COST_RANK, directionOf, Leader, SearchRecord} from '../training/search';
import {bestThresholdOf, trainingEngines, trainPredictionOf} from '../training/train-model';
import {applyModelDialog, ApplyDialogOptions, isOpen, modelLabels} from '../ui/apply-model-dialog';
import {ButtonGate} from '../ui/button-gate';
import {gridTooltip} from '../ui/data-grid';
import {editModelDialog} from '../ui/edit-model-dialog';
import {FIND_KEY, FindSettings} from '../ui/find-settings';
import {LimitInputs, TIME, TRIALS} from '../ui/limit-inputs';
import {EvaluationView, runEvaluation} from '../ui/evaluation-view';
import {ForgeApp} from '../ui/forge-app';
import {Gallery, MAX_ITEMS} from '../ui/gallery';
import {CATALOG_ACTIONS, MODEL_ACTIONS, openShareDialog} from '../ui/model-actions';
import {ModelCards} from '../ui/model-cards';
import {hasFormsViewer, ModelComparison} from '../ui/model-comparison';
import {ForgeModelHandler} from '../ui/model-handler';
import {metricIdsOf, MetricsTable, metricTooltip, storedSummary, tinyValueText} from '../ui/metrics-table';
import {modelAccordion, refreshSharing} from '../ui/model-panes';
import {PanelKey, resultsLayout} from '../ui/results-layout';
import {SessionRun} from '../ui/runs-panel';
import {saveModelDialog} from '../ui/save-model-dialog';
import {tableModelsPane} from '../ui/table-models-panel';
import {tagsOfInput} from '../ui/tags-input';
import {launchTraining} from '../ui/train-launcher';
import {CHECK_DELAY_MS, SHOWN_RUN, TrainView} from '../ui/train-view';
import {AUTO_TUNE, MANUAL, TuneRow, TuneSettings} from '../ui/tune-settings';
import {balloonWith, cancelOnceShown, cancelTaskBarEntry, catalogState, clearTestGalleryStorage, CLASSIFIERS,
  columnsOf, deleteRunsOf, demoDemog, expectReleased, expectGridsFillPanes, framesSharing,
  insertModelRow, IRIS, MEASUREMENTS, openDemoTable, openIris, openIrisFromFile, openSolubility, PREDICT_PROBABILITY,
  readUntil, refreshThenChange, REGRESSORS, runsOf, savedFixture, saveIrisModel, saveKeptGapsModel,
  saveSkipUniqueModel, saveTestModel, slowEngine, spreadRows, stuckEngine, svm, taskBarEntry, twoClassTable,
  twoSpeciesIris, unbalancedIris, useTestGalleryStorage, valuesOf, WINE, xgboost, XGBOOST_FIELDS} from './test-data';

const WAIT_MS = 5000;
const TIMEOUT = 60000;
/** The methods Method offers on iris, sorted: for Species, and for Petal.Length. */
const IRIS_CLASSIFIERS = [...CLASSIFIERS].sort();
const IRIS_REGRESSORS = [...REGRESSORS].sort();
const STUCK = 'forge-test-stuck';

/** The installed methods and {@link stuckEngine} [STUCK] answering [answer]. */
const withStuck = (answer: Promise<unknown>, calls?: {count: number}): Engine[] =>
  [...trainingEngines().engines.filter((e) => !isServerEngine(e)), stuckEngine(STUCK, answer, calls)];

/** A check answer the test gives when it chooses. */
function deferred(): {promise: Promise<boolean>; resolve: (value: boolean) => void} {
  const holder: {resolve: (value: boolean) => void} = {resolve: () => {}};
  const promise = new Promise<boolean>((resolve) => holder.resolve = resolve);
  return {promise, resolve: (value) => holder.resolve(value)};
}

/** The installed methods with SVM as {@link slowEngine} (its train calls waiting for [gate] when given). */
const withSlowSvm = (gate?: Promise<void>) =>
  [...EngineRegistry.discover().filter((e) => e.name !== 'SVM'), slowEngine(svm(), gate)];

/** Opens the Train view on [table] (with the [discovered] methods instead of the installed ones when given) and waits
 * for the training it starts by itself. */
async function openTrainView(table: DG.DataFrame, discovered?: Engine[]): Promise<TrainView> {
  const view = await TrainView.create(table, null, discovered);
  grok.shell.addView(view);
  await awaitTrained(view);
  return view;
}

/** `demog.csv` of the demo files: text SEX, RACE and USUBJID (all different), gaps in HEIGHT and WEIGHT. */
async function openDemog(): Promise<DG.DataFrame> {
  return openDemoTable('System:DemoFiles/demog.csv', 'demog');
}

async function awaitTrained(view: TrainView): Promise<void> {
  await awaitCheck(() => view.lastTraining !== undefined && !view.isTraining, 'The view does not train', TIMEOUT);
}

/** Closes [view], waits until no training of it is left, then deletes the runs of [tables]. */
async function closeTrainView(view: TrainView, tables: DG.DataFrame[]): Promise<void> {
  view.close();
  await awaitCheck(() => !view.isTraining, 'A training outlives the view', TIMEOUT);
  await deleteRunsOf(tables);
}

const isDisabled = (button: HTMLElement) => button.classList.contains('d4-disabled');
const isShown = (element: HTMLElement) => element.style.display !== 'none';
/** Every input of [view]'s form: the data, the method, the preparation, the output, the validation and the method's
 * settings. */
const formInputsOf = (view: TrainView) => [view.tableInput, view.targetInput, view.featuresInput, view.methodInput,
  ...view.preparationInputs.inputs, ...view.outputInputs.inputs, ...view.validationInputs.inputs,
  ...view.hyperparameterInputs.values()];
const formStates = (view: TrainView) => formInputsOf(view).map((input) => `${input.caption} ${input.enabled}`);
/** The enabled inputs of [view]'s form and Settings, then the enabled controls of Search options. */
const enabledInputs = (view: TrainView) => [...[...formInputsOf(view), ...view.tuneSettings.inputs]
  .filter((input) => input.enabled).map((input) => input.caption),
...view.tuneSettings.controls.filter((c) => !c.disabled).map((c) => c.getAttribute('aria-label') ?? c.className)];

/** The title of [view]'s ROC Curve panel (`ROC Curve (AUC = ...)`); '' when there is none. */
const rocTitleOf = (view: TrainView | EvaluationView) =>
  rowTitles(view).flat().find((t) => t.startsWith('ROC Curve')) ?? '';
/** The metrics grid under [root] has the row [label]. */
const hasMetric = (root: HTMLElement, label: string) =>
  gridWith(root, 'Metric').dataFrame.getCol('Metric').toList().includes(label);
/** **Train** is laid out: its row under the inputs is shown. */
const isTrainShown = (view: TrainView) => view.trainButton.offsetParent !== null;

/** The tooltip overlays the platform keeps over the disabled [button]. */
function overlaysOf(button: HTMLElement): HTMLElement[] {
  const rect = button.getBoundingClientRect();
  return Array.from(document.querySelectorAll<HTMLElement>('.d4-tooltip-overlays > span')).filter((o) => {
    const r = o.getBoundingClientRect();
    return Math.abs(r.left - rect.left) < 1 && Math.abs(r.top - rect.top) < 1;
  });
}

/** The tooltip a hover on the disabled [button] shows, through its overlay (an enabled button's tooltip answers only
 * a real mouse); '' if none. */
async function hoverText(button: HTMLElement): Promise<string> {
  ui.tooltip.hide();
  for (let i = 0; i < 30; i++) {
    const rect = button.getBoundingClientRect();
    const target = overlaysOf(button)[0] ?? button;
    target.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, clientX: rect.left + 2, clientY: rect.top + 2}));
    const text = ui.tooltip.isVisible ? ui.tooltip.root.textContent ?? '' : '';
    if (text !== '') {
      ui.tooltip.hide();
      return text;
    }
    await delay(100);
  }
  return '';
}

/** The tooltip the mouse entering the enabled [element] shows, its HTML with [isHtml]; '' if none. */
async function enterTooltip(element: HTMLElement, isHtml = false): Promise<string> {
  ui.tooltip.hide();
  const rect = element.getBoundingClientRect();
  element.dispatchEvent(new MouseEvent('mouseenter', {clientX: rect.left + 2, clientY: rect.top + 2}));
  // The platform shows a tooltip 250 ms after the mouse enters.
  await delay(600);
  const text = !ui.tooltip.isVisible ? '' : (isHtml ? ui.tooltip.root.innerHTML : ui.tooltip.root.textContent) ?? '';
  ui.tooltip.hide();
  return text;
}

/** XGBoost still retrains live on this many demog rows with three features, but each fit takes a while. */
const SLOW_ROWS = 100000;

/** Sets [view] to XGBoost's slowest valid settings on sex by age, height and weight of [demog]. */
function trainSlowly(view: TrainView, demog: DG.DataFrame): void {
  view.targetInput.value = demog.getCol('sex');
  view.featuresInput.value = columnsOf(demog, ['age', 'height', 'weight']);
  setHyperparameter(view, 'iterations', 100);
  setHyperparameter(view, 'maxDepth', 20);
}

/** The dialog **Save** of [view] opens, if any. */
function saveDialogOf(view: TrainView): DG.Dialog | undefined {
  const dialogs = new Set(DG.Dialog.getOpenDialogs().map((d) => d.root));
  view.save();
  return DG.Dialog.getOpenDialogs().find((d) => !dialogs.has(d.root));
}

const resultsText = (view: TrainView) => view.root.querySelector('.forge-train-results')?.textContent ?? '';

/** The readiness attributes of the results of [view]: `<aria-busy>/<aria-invalid>`. */
const resultsState = (view: DG.ViewBase, results = view.root.querySelector('.forge-train-results')) =>
  `${results?.getAttribute('aria-busy')}/${results?.getAttribute('aria-invalid')}`;
/** The run whose result the results of [view] show (`SHOWN_RUN`), or null. */
const shownRun = (view: DG.ViewBase) =>
  view.root.querySelector('.forge-train-results')?.getAttribute(SHOWN_RUN) ?? null;

/** The value of the hyperparameter input [name] of [view] becomes [value]. */
function setHyperparameter(view: TrainView, name: string, value: number): void {
  const input = view.hyperparameterInputs.get(name);
  if (input === undefined)
    throw new Error(`No ${name} input`);
  input.value = value;
}

/** The body and the summary of the collapsible group captioned [caption] inside [root]. */
function groupOf(root: HTMLElement, caption: string): {body: HTMLElement; summary: HTMLElement} {
  for (const group of Array.from(root.getElementsByClassName('forge-group'))) {
    const label = group.querySelector('.forge-group-header label');
    const body = group.querySelector('.forge-group-body');
    const summary = group.querySelector('.forge-group-summary');
    if (label?.textContent === caption && body instanceof HTMLElement && summary instanceof HTMLElement)
      return {body, summary};
  }
  throw new Error(`No group ${caption}`);
}

/** Checks that the captions of the groups [captions] inside [root] start at the same x. */
function expectAlignedCaptions(root: HTMLElement, captions: string[]): void {
  const lefts = captions.map((caption) => Array.from(root.querySelectorAll('.forge-group-header label'))
    .find((label) => label.textContent === caption)?.getBoundingClientRect().left ?? NaN);
  expect(lefts.every((left) => Math.abs(left - lefts[0]) < 0.5), true, `Caption positions ${lefts.join(', ')}`);
}

/** The option captions of the radio [input], in order. */
const optionsOfRadio = (input: DG.InputBase) =>
  Array.from(input.root.querySelectorAll('.ui-radio-button label'), (label) => label.textContent ?? '');

/** Checks that the two options of the shown radio [input] sit on one line. */
function expectOptionsInRow(input: DG.InputBase): void {
  const options = Array.from(input.root.querySelectorAll<HTMLElement>('.ui-radio-button'));
  expect(options.length === 2 && options.every((o) => o.offsetHeight > 0 && o.offsetTop === options[0].offsetTop),
    true, `Option positions ${options.map((o) => o.offsetTop).join(', ')}`);
}

/** Waits until [dialog] is closed, [table] has the Species prediction and the application of [modelId] is recorded. */
async function expectApplied(dialog: DG.Dialog, table: DG.DataFrame, modelId: string): Promise<void> {
  await awaitCheck(() => !dialog.root.isConnected, 'The dialog stays open', WAIT_MS);
  await awaitCheck(() => table.col('Species (predicted)') !== null, 'No prediction column', WAIT_MS);
  expect(await recordedApplications(modelId), 1);
}

/** The count of [modelId]'s applications once one is recorded, or after 5 s. */
const recordedApplications = (modelId: string) =>
  readUntil(() => forgeDb.applications.query().where('model_id', '=', modelId).count(), (count) => count > 0);

/** Waits (errors swallowed) for the application a dialog started, so that deleting [modelId] leaves no row behind. */
const applicationSettled = (modelId: string) => recordedApplications(modelId).catch(() => {});

/** The texts of the options of the choice input captioned [caption] in [dialog], in list order. */
const optionsOf = (dialog: DG.Dialog, caption: string) =>
  Array.from(dialog.input(caption).root.querySelectorAll('option')).map((o) => o.textContent ?? '');

const pressEnter = (dialog: DG.Dialog, caption: string) => dialog.input(caption).input.dispatchEvent(
  new KeyboardEvent('keydown', {key: 'Enter', keyCode: 13, bubbles: true}));

/** The platform grid inside [root] whose table has the column [column]. */
function gridWith(root: HTMLElement, column: string): DG.Grid {
  for (const element of Array.from(root.querySelectorAll('.d4-grid'))) {
    const grid: unknown = DG.toJs(DG.Widget.find(element));
    if (grid instanceof DG.Grid && grid.dataFrame.col(column) !== null)
      return grid;
  }
  throw new Error(`No grid with the column ${column}`);
}

function headerTooltip(grid: DG.Grid, column: string): string {
  const gridColumn = grid.columns.byName(column);
  if (gridColumn === null)
    throw new Error(`No grid column ${column}`);
  return gridTooltip(grid, DG.GridCell.createColHeader(gridColumn));
}

/** Checks that the attached [grid] shows its header and all its rows without scrolling; returns what was measured. */
function expectRowsFit(name: string, grid: DG.Grid): string {
  const needed = grid.colHeaderHeight + grid.dataFrame.rowCount * grid.props.rowHeight;
  const scroll = grid.vertScroll.root;
  const measured = `${name}: ${grid.root.clientHeight}px for ${needed}px (header ${grid.colHeaderHeight}px, ` +
    `${grid.dataFrame.rowCount} rows of ${grid.props.rowHeight}px; scroll bar ` +
    `${scroll.offsetWidth}x${scroll.offsetHeight} ${getComputedStyle(scroll).visibility})`;
  expect(grid.root.clientHeight >= needed, true, measured);
  return measured;
}

/** The browser autofill setting of the Tags text box inside [root]. */
const tagsAutofill = (root: Element) =>
  root.querySelector('input.d4-tags-selector-input')?.getAttribute('autocomplete') ?? null;

/** The ids of the models the context panel compares, sorted and joined; '' when it shows no comparison. */
function comparedIds(): string {
  const shown: unknown = grok.shell.o;
  return shown instanceof ModelComparison ? shown.rows.map((r) => r.id).sort().join() : '';
}

/** Whether the context panel shows the model [modelId]. */
function shows(modelId: string): boolean {
  const shown: unknown = grok.shell.o;
  return shown instanceof DG.DomainRow && shown.id === modelId;
}

/** The tests' own storage of the tune settings: the user's stay as they are, and each category's Train views start
 * from the defaults. */
const TEST_TUNE_STORAGE = 'Forge Tune tests';

function useTestTuneStorage(): void {
  TuneSettings.storage = TEST_TUNE_STORAGE;
  clearTestTuneStorage();
}

function clearTestTuneStorage(): void {
  grok.userSettings.put(TEST_TUNE_STORAGE, {});
}

category('UI', () => {
  test('the catalog refreshes itself after a save and a delete', async () => {
    const view = await ForgeApp.create();
    grok.shell.addView(view);
    const name = `forge-test-model-${Date.now()}`;
    const rowOf = () => view.models.getCol('name').toList().indexOf(name);
    let modelId: string | undefined;
    try {
      [{id: modelId}] = await forgeDb.models.insert({...XGBOOST_FIELDS, name, task: 'regression', target_name: 'y',
        storage_mode: 'none'});
      modelsChanged.next();
      await awaitCheck(() => rowOf() >= 0, 'The saved model does not appear', WAIT_MS);
      expect(catalogState(view.root), '1/false');

      await deleteModel(modelId);
      modelId = undefined;
      await awaitCheck(() => rowOf() < 0, 'The deleted model stays', WAIT_MS);
    } finally {
      if (modelId !== undefined)
        await deleteModel(modelId);
      view.close();
    }
  }, {timeout: 30000});

  test('catalog menu item', async () => {
    const topMenu = DG.Func.find({package: 'Forge', name: 'forgeModels'})[0]?.topMenu;
    expect(topMenu === null || topMenu === '', true, `Unexpected top menu: ${topMenu}`);
  });

  test('Train view opens with Save in the ribbon and one form with a hidden Train', async () => {
    const iris = await openIris();
    const view = await TrainView.create(iris);
    grok.shell.addView(view);
    try {
      expect(view.name, 'Predictive model');
      expect(view.getIcon().classList.contains('svg-model'), true, 'The tab icon is not the model icon');
      const text = view.root.textContent ?? '';
      for (const caption of ['Table', 'Target', 'Features', 'Method', 'Performance'])
        expect(text.includes(caption), true, `${caption} is missing`);
      expect(text.includes('Choose the target and the features; the model trains as you change them.'), true,
        'No hint before the first training');
      expect(text.includes('Engine'), false, 'The word Engine is shown');
      const ribbon = view.getRibbonPanels().flat();
      expect(ribbon.some((e) => e.contains(view.saveButton)), true, 'Save is not in the ribbon');
      expect(ribbon.some((e) => e.contains(view.trainButton)), false, 'Train is in the ribbon');
      expect(view.root.querySelector('.fa-sync') === null, true, 'A sync icon');
      const forms = view.root.querySelectorAll('.ui-form');
      expect(forms.length, 1, 'Data and Method are not one form');
      const trainRow = view.trainButton.closest('.forge-train-row');
      expect(trainRow?.parentElement === forms[0] && trainRow.previousElementSibling === view.groups[4].root, true,
        'Train is not right after Validation in the form');
      expect(view.root.querySelectorAll('.ui-split-h-divider').length, 1, 'No splitter between the form and Results');
      expect(isDisabled(view.saveButton), true, 'Save is enabled before training');
      expect(isTrainShown(view), false, 'Train is shown before the check');
      expect(view.targetInput.value?.name, 'Species');
      expectArray(view.featuresInput.value.map((c) => c.name), MEASUREMENTS);
      expect(CHECK_DELAY_MS, 200);

      await awaitTrained(view);
      expect(isTrainShown(view), false, 'Train is shown for a method that retrains by itself');
      const iterations = view.hyperparameterInputs.get('iterations');
      const labelOf = (input: DG.InputBase | undefined) => input?.root.querySelector('label')?.getBoundingClientRect();
      const [table, method] = [labelOf(view.tableInput), labelOf(iterations)];
      expect(table !== undefined && method !== undefined && Math.abs(table.right - method.right) < 0.5, true,
        `Label right edges ${table?.right}, ${method?.right}`);
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: TIMEOUT});

  test('Train view trains iris by itself and saves the model once', async () => {
    const iris = await openIris();
    const runs = () => runsOf([iris]);
    const models = () => forgeDb.models.query().where('dataset_name', '=', iris.name);
    const view = await TrainView.create(iris);
    grok.shell.addView(view);
    let measured = '';
    try {
      await awaitCheck(() => view.isTraining, 'The view does not start training by itself', WAIT_MS);
      expect(isDisabled(view.saveButton), true, 'Save is enabled while training');
      await awaitTrained(view);
      const recorded = await runs();
      expectArray(recorded.map((r) => `${r.status} ${r.engine_name}`), ['completed XGBoost']);
      expect(isDisabled(view.saveButton), false, 'Save stays disabled after training');
      const text = view.root.textContent ?? '';
      const results = gridWith(view.root, 'Metric');
      const metrics = results.dataFrame.getCol('Metric').toList();
      for (const caption of ['Accuracy', 'F1 score'])
        expect(metrics.includes(caption), true, `${caption} is missing: ${metrics}`);
      expect(metrics.includes('Recall (Sensitivity)'), false, 'Recall (Sensitivity) on three classes');
      expect(text.includes('Split seed: 42'), true, 'No split seed line');
      expect(text.includes('Save...'), false, 'Save... is still in Results');
      measured = expectRowsFit('Results', results);

      const name = `forge-test-model-${Date.now()}`;
      await Promise.all([view.saveModelAs({name, description: '', tags: ['x', ' y ', 'x']}),
        view.saveModelAs({name: `${name}-again`, description: '', tags: []})]);
      expect(isDisabled(view.saveButton), true, 'Save stays enabled after saving');
      await view.saveModelAs({name: `${name}-later`, description: '', tags: []});
      const saved = await models();
      expect(saved.length, 1, 'The same training was saved more than once');
      expect(saved[0].name, name);
      expect(saved[0].tags, 'x, y');
      expect((await runs())[0].model_id, saved[0].id);

      await view.train();
      expectArray((await runs()).map((r) => r.status), ['completed', 'completed']);
      expect(isDisabled(view.saveButton), false, 'Train gave no new training to save');
    } finally {
      for (const model of await models())
        await deleteModel(model.id);
      await closeTrainView(view, [iris]);
    }
    return measured;
  }, {timeout: TIMEOUT});

  test('a failing training records a failed run', async () => {
    const iris = await openIris();
    const view = await openTrainView(iris);
    try {
      // In the range EDA declares, but refused by its trainer.
      setHyperparameter(view, 'maxDepth', 0);
      const recorded = await readUntil(() => runsOf([iris]), (runs) => runs.some((r) => r.status === 'failed'));
      const failed = recorded.filter((r) => r.status === 'failed');
      expect(failed.length, 1);
      expect((failed[0].error ?? '') !== '', true, 'The failed run has no error');
      await awaitCheck(() => !view.isTraining, 'The training does not end', WAIT_MS);
      expect(view.lastTraining === undefined, true, 'A failed training left a result');
      expect(isDisabled(view.saveButton), true, 'Save is enabled after a failed training');
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: TIMEOUT});

  test('a bad selection disables Train, marks the input and trains nothing', async () => {
    const iris = await openIris();
    const view = await openTrainView(iris);
    const columns = (names: string[]) => names.map((name) => iris.getCol(name));
    try {
      view.targetInput.value = iris.getCol('Petal.Length');
      view.featuresInput.value = columns(['Species', 'Sepal.Length', 'Sepal.Width', 'Petal.Width']);
      await awaitCheck(() => isShown(view.preparationInputs.oneHot.root), 'One-hot encoding is not offered', WAIT_MS);
      view.preparationInputs.oneHot.value = false;
      await awaitCheck(() => view.featuresInput.validity !== null, 'Features is not marked', WAIT_MS);
      expect(view.featuresInput.validity?.includes('Species'), true, `${view.featuresInput.validity}`);
      expect(!isTrainShown(view) && isDisabled(view.trainButton), true, 'Train is shown for XGBoost or enabled');
      expect(view.targetInput.validity === null, true, `Target is marked: ${view.targetInput.validity}`);
      expect(view.isTraining || view.lastTraining !== undefined, false, 'A bad selection trains');
      expect((await runsOf([iris])).length, 1, 'A bad selection recorded a run');

      view.featuresInput.value = columns(['Sepal.Length', 'Sepal.Width', 'Petal.Width']);
      await awaitCheck(() => view.featuresInput.validity === null, 'Features stays marked', WAIT_MS);
      await awaitTrained(view);
      expect(view.methodInput.value, 'Linear Regression');
      expectArray(view.methodInput.items.slice(1).sort(), IRIS_REGRESSORS);
      expect(view.lastTraining?.engine.name, 'Linear Regression');
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: TIMEOUT});

  test('a change during a training supersedes it without a run row', async () => {
    const wine = await grok.data.files.openTable(WINE);
    wine.name = `forge-test-wine-${Date.now()}`;
    const view = await TrainView.create(wine);
    grok.shell.addView(view);
    try {
      await awaitCheck(() => view.isTraining, 'The view does not train by itself', TIMEOUT);
      const first = view.methodInput.value;
      const quality = wine.getCol('quality');
      view.targetInput.value = quality;
      await awaitCheck(() => view.featuresInput.validity !== null, 'Features is not marked', WAIT_MS);
      view.featuresInput.value = view.featuresInput.value.filter((c) => c.name !== 'quality');
      await awaitTrained(view);
      expect(view.methodInput.value, 'PLS Regression');
      const recorded = await runsOf([wine]);
      expect(recorded.map((r) => `${r.status} ${r.engine_name} ${r.target_name}`).join(),
        'completed PLS Regression quality', `The first training (${first}) was recorded`);
    } finally {
      await closeTrainView(view, [wine]);
    }
  }, {timeout: TIMEOUT});

  test('a cancelled training is recorded; a slow method waits for Train', async () => {
    const demog = demoDemog(SLOW_ROWS);
    const view = await TrainView.create(demog);
    grok.shell.addView(view);
    try {
      trainSlowly(view, demog);
      for (let i = 0; i < TIMEOUT / 10 && !view.isTraining; i++)
        await delay(10);
      expect(`${view.isTraining} ${view.methodInput.value} ${resultsState(view)}`, 'true XGBoost true/false');
      await cancelOnceShown('Training XGBoost model', TIMEOUT);
      await awaitCheck(() => !view.isTraining, 'The cancelled training does not end', TIMEOUT);
      expect(resultsState(view), 'false/true');
      const recorded = await runsOf([demog]);
      expectArray(recorded.map((r) => `${r.status} ${r.engine_name}`), ['cancelled XGBoost']);
      expect(view.lastTraining === undefined && isDisabled(view.saveButton), true, 'A cancelled training is saveable');

      view.methodInput.value = 'SVM';
      const hint = `SVM takes a while on ${SLOW_ROWS} rows, so it does not retrain on every change. Click Train.`;
      await awaitCheck(() => (view.root.textContent ?? '').includes(hint), 'No hint for a slow method', WAIT_MS);
      await delay(500);
      expect(`${view.isTraining} ${resultsState(view)}`, 'false false/false', 'SVM trains by itself');
      expect(isTrainShown(view) && !isDisabled(view.trainButton), true, 'Train is not offered for a slow method');
      expect((await runsOf([demog])).length, 1, 'A slow method recorded a run');

      setHyperparameter(view, 'cost', 5000);
      await awaitCheck(() => view.hyperparameterInputs.get('cost')?.validity !== null, 'No mark on Penalty', WAIT_MS);
      await awaitCheck(() => isDisabled(view.trainButton), 'Train stays enabled with Penalty 5000', WAIT_MS);
      await delay(CHECK_DELAY_MS * 2);
      expect(isTrainShown(view), true, 'Train is hidden for a slow method with a problem');
      const reason = `Penalty: ${view.hyperparameterInputs.get('cost')?.validity}`;
      expect(await hoverText(view.trainButton), reason);
      expect(resultsText(view).includes('Fix the settings.'), true, 'The Results hint is not Fix the settings.');
    } finally {
      await closeTrainView(view, [demog]);
    }
  }, {timeout: TIMEOUT});

  test('Method lists the applicable methods, keeps the user\'s choice and the values per method', async () => {
    const iris = await openIris();
    const other = await openIris();
    const tables = [grok.shell.addTable(iris), grok.shell.addTable(other)];
    const view = await openTrainView(iris);
    try {
      expectArray(view.methodInput.items.slice(1).sort(), IRIS_CLASSIFIERS);
      expect(view.methodInput.value, 'XGBoost');
      const iterations = () => view.hyperparameterInputs.get('iterations');
      setHyperparameter(view, 'iterations', 25);
      await awaitTrained(view);

      view.methodInput.value = 'SVM';
      await awaitTrained(view);
      expect(view.lastTraining?.engine.name, 'SVM');
      expect(iterations() === undefined, true, 'SVM shows the XGBoost hyperparameters');
      const save = saveDialogOf(view);
      expect((save?.root.textContent ?? '').includes('contains the training rows'), false, 'A warning for SVM');
      save?.close();
      view.featuresInput.value = columnsOf(iris, ['Sepal.Length', 'Petal.Length', 'Petal.Width']);
      await awaitTrained(view);
      expect(view.methodInput.value, 'SVM', 'A Features change dropped the user\'s choice');

      view.methodInput.value = 'XGBoost';
      await awaitCheck(() => iterations()?.value === 25, 'The XGBoost values are lost', WAIT_MS);
      view.tableInput.value = other;
      await awaitTrained(view);
      expect(iterations()?.value, 25, 'A Table change lost the values');
      expect(view.methodInput.value, 'XGBoost');

      view.methodInput.value = 'Softmax';
      await awaitTrained(view);
      view.featuresInput.value = columnsOf(other, ['Sepal.Length', 'Sepal.Width', 'Petal.Length']);
      await awaitTrained(view);
      view.targetInput.value = other.getCol('Petal.Width');
      const balloon = 'Softmax cannot learn from this selection; Linear Regression is chosen for now.';
      await awaitCheck(() => balloonWith(balloon) !== null, 'No balloon for a dropped choice', WAIT_MS);
      await awaitTrained(view);
      expect(view.methodInput.value, 'Linear Regression');
    } finally {
      await closeTrainView(view, tables);
      for (const table of tables)
        grok.shell.closeTable(table);
    }
  }, {timeout: 120000});

  test('Train view groups its inputs and opens a group with a problem', async () => {
    const iris = await openIris();
    const view = await openTrainView(iris);
    try {
      const data = groupOf(view.root, 'Data');
      const preparation = groupOf(view.root, 'Preparation');
      const method = groupOf(view.root, 'Method');
      const output = groupOf(view.root, 'Output');
      expect(isShown(data.body) && isShown(preparation.body) && isShown(method.body) && isShown(output.body), true,
        'A group starts collapsed');
      expect(data.body.contains(view.featuresInput.root), true, 'Features is not under Data');
      expect(preparation.body.contains(view.missingValuesInputs.choice.root), true,
        'Missing values is not under Preparation');
      expect(isShown(view.preparationInputs.group.root), false, 'Preparation is shown for iris');
      const shown = view.outputInputs.inputs.filter((i) => isShown(i.root)).map((i) => i.caption).join();
      expect(isShown(view.outputInputs.group.root) && shown === 'Predict probability', true, `Shown: ${shown}`);
      expect(output.body.contains(view.outputInputs.predictProbability.root), true,
        'Predict probability is not under Output');
      expectArray(view.groups.map((g) => Array.from(g.root.parentElement?.children ?? []).indexOf(g.root)),
        [0, 1, 2, 3, 4]);
      const iterations = view.hyperparameterInputs.get('iterations');
      expect(iterations !== undefined && method.body.contains(iterations.root), true, 'Iterations is not under Method');
      expect(method.body.contains(view.methodInput.root), true, 'Method is not under Method');
      expect(method.body.querySelector('.forge-note') === null, true, 'A server note for XGBoost');

      view.groups[0].setExpanded(false);
      expect(isShown(data.body), false, 'Data does not collapse');
      expectAlignedCaptions(view.root, ['Data', 'Method']);
      const target = view.targetInput.value;
      if (target === null)
        throw new Error('No target');
      view.featuresInput.value = [...view.featuresInput.value, target];
      await awaitCheck(() => view.featuresInput.validity !== null && isShown(data.body),
        'Data stays collapsed with a problem on Features', WAIT_MS);
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: TIMEOUT});

  test('Train view drops the old form\'s subscriptions on a Table change', async () => {
    const tables = [grok.shell.addTable(await openIris()), grok.shell.addTable(await openIris())];
    const view = await openTrainView(tables[0]);
    try {
      const oldData = view.groups[0];
      const oldFeatures = view.featuresInput;
      const subCount = view.subs.length;
      for (const table of [tables[1], tables[0]]) {
        const features = view.featuresInput;
        view.tableInput.value = table;
        await awaitCheck(() => view.featuresInput !== features, 'The form was not rebuilt', WAIT_MS);
      }
      expect(view.subs.length, subCount);

      const failing = () => 'forge-test-problem';
      for (const [group, input] of [[oldData, oldFeatures], [view.groups[0], view.featuresInput]] as const) {
        group.setExpanded(false);
        input.addValidator(failing);
        input.validate();
      }
      await awaitCheck(() => view.groups[0].isExpanded, 'The current Data group does not open on an error', WAIT_MS);
      expect(oldData.isExpanded, false, 'A replaced form still opens its group on an error');
    } finally {
      await closeTrainView(view, tables);
      for (const table of tables)
        grok.shell.closeTable(table);
    }
  }, {timeout: TIMEOUT});

  test('Save model needs a name', async () => {
    const dialog = saveModelDialog('Iris model', async () => {});
    dialog.show();
    const ok = dialog.getButton('OK');
    const name = dialog.input('Name');
    try {
      expect(isDisabled(ok), false, 'OK is disabled with a name');
      name.value = '';
      await awaitCheck(() => isDisabled(ok), 'OK is enabled without a name', WAIT_MS);
      name.value = '   ';
      await awaitCheck(() => isDisabled(ok) && name.validity !== null, 'A blank name is accepted', WAIT_MS);
      name.value = 'Iris model';
      await awaitCheck(() => !isDisabled(ok), 'OK stays disabled with a name', WAIT_MS);
    } finally {
      dialog.close();
    }
  });

  test('a disabled ButtonGate keeps one overlay that shows the latest reason', async () => {
    let reason = 'forge-test-reason-0';
    const button = ui.button('forge-test-gate', () => {});
    const gate = new ButtonGate(button, () => reason);
    document.body.append(button);
    try {
      for (let i = 1; i <= 5; i++) {
        gate.update();
        reason = `forge-test-reason-${i}`;
      }
      gate.update();
      // The overlay takes the button's place on the platform's next 100 ms tick.
      await delay(300);
      expect(overlaysOf(button).length, 1);
      expect(await hoverText(button), 'forge-test-reason-5');
    } finally {
      button.remove();
    }
  });

  test('an enabled ButtonGate drops the disabled reason from the button\'s tooltip', async () => {
    let reason: string | null = 'forge-test-reason';
    const plain = ui.button('forge-test-plain', () => {});
    const described = ui.button('forge-test-described', () => {});
    const gates = [new ButtonGate(plain, () => reason), new ButtonGate(described, () => reason, 'forge-test-enabled')];
    document.body.append(plain, described);
    try {
      for (const gate of gates)
        gate.update();
      reason = null;
      for (const gate of gates)
        gate.update();
      expectArray([await enterTooltip(plain), await enterTooltip(described)], ['', 'forge-test-enabled']);
    } finally {
      plain.remove();
      described.remove();
    }
  });

  test('the core\'s entry points exist without a menu of their own', async () => {
    for (const name of ['forgeTrain', 'forgeApply']) {
      const func = DG.Func.find({package: 'Forge', name})[0];
      expect(func !== undefined, true, `${name} is not registered`);
      expect(!func.topMenu, true, `Unexpected top menu of ${name}: ${func.topMenu}`);
    }
    const app = DG.Func.find({package: 'Forge', name: 'forgeApp'})[0];
    expectArray(app?.outputs.map((p) => p.propertyType) ?? [], [DG.TYPE.VIEW]);
    const inputs = DG.Func.find({package: 'Forge', name: 'applyModel'})[0]?.inputs.map((p) => p.name) ?? [];
    expectArray(inputs, ['model', 'table', 'columnNamesMap', 'showProgress']);
  });
});

category('UI', () => {
  test('Forge has no app and no view function', async () => {
    expect(DG.Func.find({package: 'Forge', tags: [DG.FUNC_TYPES.APP]}).length, 0);
    expect(DG.Func.find({package: 'Forge', tags: ['view']}).length, 0);
    expect(DG.Func.find({package: 'Forge', name: 'forgeApp'}).length, 1);
  });

  test('the view is at /models only', async () => {
    const view = await ForgeApp.create();
    grok.shell.addView(view);
    try {
      expect(`${view.name} ${view.type} ${view.path.endsWith('/models')}`, 'Models models true', view.path);
      const accepted = ['/models', '/MODELS/', '/Models', '/ml', '/ml/models', '/ml/methods', '/modelsx', '/models/x',
        '/'].map((path) => view.acceptsPath(path));
      expectArray(accepted, [true, true, true, false, false, false, false, false, false]);
      view.handlePath('/models');
      expect(`${grok.shell.v === view} ${view.path.endsWith('/models')}`, 'true true', view.path);
    } finally {
      view.close();
    }

    // forgeApp runs in the package's bundle, whose views this bundle's ForgeApp does not recognize: the view is checked
    // through its name, type and path.
    const opened: unknown = await grok.functions.call('Forge:forgeApp', {});
    if (!(opened instanceof DG.ViewBase))
      throw new Error(`forgeApp returned ${opened}`);
    grok.shell.addView(opened);
    try {
      expect(`${opened.name} ${opened.type} ${opened.path.endsWith('/models')}`, 'Models models true', opened.path);
    } finally {
      opened.close();
    }
  });
});

category('UI: Train', () => {
  test('Train view hides Missing values for a method that handles them', async () => {
    const iris = await openIris();
    iris.getCol('Sepal.Width').set(3, null);
    const view = await openTrainView(iris);
    const listItems = () => Array.from(view.root.querySelectorAll('ul > li'), (li) => li.textContent);
    try {
      const choice = view.missingValuesInputs.choice;
      expect(view.methodInput.value, 'XGBoost');
      expect(isShown(choice.root), false, 'Missing values is shown for XGBoost');
      expect(isShown(view.outputInputs.predictProbability.root), true, 'No Predict probability for Species');
      await awaitCheck(() => listItems().includes('Missing values: kept (the method handles them)'),
        `No kept line: ${listItems().join(' | ')}`, WAIT_MS);
      expect(view.lastTraining?.result.rowCount, 150);

      view.methodInput.value = 'SVM';
      await awaitTrained(view);
      expect(isShown(choice.root), true, 'Missing values is hidden for SVM');
      expect(choice.value, 'Skip rows');
      await awaitCheck(() => listItems().includes('Rows: 149 used, 1 skipped (missing values)') &&
        !listItems().includes('Missing values: kept (the method handles them)'), listItems().join(' | '), WAIT_MS);
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: TIMEOUT});

  test('Forge Check: checking a molecule feature switches Method to the method that takes it', async () => {
    const table = await openSolubility();
    const view = await openTrainView(table);
    try {
      view.featuresInput.value = [table.getCol('smiles')];
      view.targetInput.value = table.getCol('SOL');
      await waitUntil(() => view.methodInput.value === CHECK_METHOD, `Method is ${view.methodInput.value}`, WAIT_MS);
      await waitUntilAsync(async () => (await runsOf([table])).some((r) => r.engine_name === CHECK_METHOD &&
        r.status === 'completed'), `No completed run of ${CHECK_METHOD}`, TIMEOUT / 2);
      expect(view.methodInput.items.includes('XGBoost'), false, view.methodInput.items.join(', '));
      expect(isShown(view.preparationInputs.group.root), false, 'Preparation is shown');
      expect(isShown(view.outputInputs.group.root), false, 'Output is shown');
      expect(balloonWith('cannot learn') === null, true, 'A balloon named a refused method');
    } finally {
      await closeTrainView(view, [table]);
    }
  }, {timeout: TIMEOUT});

  test('Forge Check: a chosen method that does not take a molecule feature says so', async () => {
    const table = await openSolubility();
    const view = await TrainView.create(table, 'XGBoost');
    grok.shell.addView(view);
    const refusal = `XGBoost cannot learn from Molecule features. Uncheck: smiles. ${CHECK_METHOD} is chosen for now.`;
    try {
      await awaitTrained(view);
      view.featuresInput.value = [table.getCol('smiles'), table.getCol('MolLogP')];
      await waitUntil(() => balloonWith(refusal) !== null, 'No balloon for the refused method', WAIT_MS);
      await waitUntil(() => view.methodInput.value === CHECK_METHOD, `Method is ${view.methodInput.value}`, WAIT_MS);
    } finally {
      await closeTrainView(view, [table]);
    }
  }, {timeout: TIMEOUT});

  test('Target and Features follow a column removed from the table', async () => {
    const iris = await openIris();
    const view = await openTrainView(iris);
    try {
      iris.columns.remove('Petal.Width');
      expectArray(view.featuresInput.value.map((c) => c.name), ['Sepal.Length', 'Sepal.Width', 'Petal.Length']);
      await awaitTrained(view);
      expect(view.lastTraining?.result.features.columns.length, 3);
      iris.columns.remove('Species');
      expect(view.targetInput.value === null, true, `Target stays ${view.targetInput.value?.name}`);
      await awaitCheck(() => view.targetInput.validity !== null, 'An empty Target is not marked', WAIT_MS);
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: TIMEOUT});

  test('a training superseded while its task-bar progress is shown writes no run', async () => {
    const demog = demoDemog(SLOW_ROWS);
    const view = await TrainView.create(demog);
    grok.shell.addView(view);
    const progressShown = () => taskBarEntry('Training XGBoost model') !== undefined;
    try {
      trainSlowly(view, demog);
      for (let i = 0; i < TIMEOUT / 10 && !(view.isTraining && progressShown()); i++)
        await delay(10);
      expect(view.isTraining && progressShown(), true, 'XGBoost does not train with a progress');
      expect((view.root.textContent ?? '').includes('Training XGBoost'), false, 'Results names the training');
      setHyperparameter(view, 'maxDepth', 2);
      await awaitTrained(view);
      const recorded = await runsOf([demog]);
      expectArray(recorded.map((r) => `${r.status} ${r.engine_name} ${r.hyperparameters?.maxDepth}`),
        ['completed XGBoost 2']);
      expect(shownRun(view), recorded[0].id, 'The results show another run than the one recorded');
    } finally {
      await closeTrainView(view, [demog]);
    }
  }, {timeout: TIMEOUT});

  test('an invalid hyperparameter blocks training and saving', async () => {
    const iris = await openIris();
    const view = await openTrainView(iris);
    try {
      await awaitCheck(() => resultsState(view) === 'false/false', 'The results are not ready', WAIT_MS);
      setHyperparameter(view, 'iterations', 500);
      await awaitCheck(() => view.hyperparameterInputs.get('iterations')?.validity !== null,
        'Iterations is not marked', WAIT_MS);
      await delay(CHECK_DELAY_MS * 2);
      expect(resultsState(view), 'false/true');
      expect(!isTrainShown(view) && isDisabled(view.trainButton), true, 'Train is shown for XGBoost, which retrains');
      expect(resultsText(view).includes('Fix the settings.'), true, 'The Results hint is not Fix the settings.');
      await view.train();
      await delay(500);
      expect(view.isTraining || view.lastTraining !== undefined, false, 'An invalid value trains');
      expect(isDisabled(view.saveButton), true, 'Save is enabled with an invalid value');
      expect(saveDialogOf(view) === undefined, true, 'Save opened a dialog');
      expect((await runsOf([iris])).length, 1, 'An invalid value recorded a run');

      setHyperparameter(view, 'iterations', 50);
      await awaitTrained(view);
      expect(isTrainShown(view), false, 'Train stays shown with a valid value');
      expect((await runsOf([iris])).length, 2);
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: TIMEOUT});

  test('Train is shown for a slow method, disabled after a training and enabled by a change', async () => {
    const table = twoClassTable();
    const view = await openTrainView(table, withSlowSvm());
    const first = view.lastTraining?.engine.name;
    try {
      expect(isTrainShown(view), false, `Train is shown for ${first}`);
      view.methodInput.value = 'SVM';
      await awaitCheck(() => isTrainShown(view) && !isDisabled(view.trainButton), 'Train is not offered for SVM',
        WAIT_MS);
      view.trainButton.click();
      await awaitCheck(() => view.isTraining, 'Train does not train', WAIT_MS);
      expect(isDisabled(view.trainButton), true, 'Train is enabled while training');
      await awaitTrained(view);
      expect(view.lastTraining?.engine.name, 'SVM');
      expect(isDisabled(view.trainButton), true, 'Train is enabled after the training');
      setHyperparameter(view, 'cost', 2);
      await awaitCheck(() => !isDisabled(view.trainButton), 'Train stays disabled after a change', WAIT_MS);
      expect(view.isTraining, false, 'SVM retrains by itself');

      setHyperparameter(view, 'cost', 5000);
      await awaitCheck(() => view.hyperparameterInputs.get('cost')?.validity !== null, 'No mark on Penalty', WAIT_MS);
      await delay(CHECK_DELAY_MS * 2);
      expect(isTrainShown(view) && isDisabled(view.trainButton), true, 'Train is not disabled with Penalty 5000');
      await view.train();
      expect(view.isTraining || view.lastTraining !== undefined, false, 'Penalty 5000 trains');
      setHyperparameter(view, 'cost', 2);
      await awaitCheck(() => !isDisabled(view.trainButton), 'Train stays disabled with Penalty 2', WAIT_MS);
      expectArray((await runsOf([table])).map((r) => `${r.status} ${r.engine_name}`).sort(),
        [`completed ${first}`, 'completed SVM'].sort());
    } finally {
      await closeTrainView(view, [table]);
    }
  }, {timeout: 120000});

  test('Train disables the form and ignores Runs rows; each input\'s state comes back after the training', async () => {
    const table = twoClassTable();
    const view = await openTrainView(table, withSlowSvm());
    try {
      view.methodInput.value = 'SVM';
      await awaitCheck(() => isTrainShown(view) && !isDisabled(view.trainButton), 'Train is not offered for SVM',
        WAIT_MS);
      const degree = view.hyperparameterInputs.get('degree');
      expect(degree?.enabled, false, 'Degree is enabled for the RBF kernel');
      const before = formStates(view);
      view.trainButton.click();
      expect(view.isTraining, true, 'Train does not train');
      expectArray(enabledInputs(view), []);
      // Tune is gated only in Auto-tune: in Manual it is hidden.
      expect(isDisabled(view.trainButton) && (isDisabled(view.tuneButton) || view.tuneButton.offsetParent === null),
        true, 'Train or a shown Tune is enabled');
      view.restoreRun(view.runs[0]);
      expectArray(enabledInputs(view), []);
      await awaitTrained(view);
      expect(view.lastTraining?.engine.name, 'SVM', 'The training was superseded');
      expectArray(formStates(view), before);
      expect(degree?.enabled, false, 'Degree is enabled after the training');
    } finally {
      await closeTrainView(view, [table]);
    }
  }, {timeout: 120000});

  test('a cancel unlocks the form; a training completed after a cancelled one clears the invalid mark', async () => {
    const demog = demoDemog(1000);
    let open: () => void = () => {};
    const gate = new Promise<void>((resolve) => open = resolve);
    const view = await TrainView.create(demog, null, withSlowSvm(gate));
    view.targetInput.value = demog.getCol('sex');
    view.methodInput.value = 'SVM';
    grok.shell.addView(view);
    try {
      await awaitCheck(() => isTrainShown(view) && !isDisabled(view.trainButton), 'Train is not offered for SVM',
        WAIT_MS);
      const before = formStates(view);
      const cancelled = view.train();
      expectArray(enabledInputs(view), []);
      await cancelOnceShown('Training SVM model', WAIT_MS);
      open();
      await cancelled;
      expect(`${view.runs[0]?.status} ${resultsState(view)}`, 'cancelled false/true',
        'The cancel came after the training ended');
      expectArray(formStates(view), before);
      await view.train();
      await awaitCheck(() => resultsState(view).startsWith('false/'), 'The results stay busy', WAIT_MS);
      expect(`${view.runs[0]?.status} ${resultsState(view)}`, 'completed false/false');
      expect(shownRun(view), view.lastTraining?.runId ?? null, 'The results show another run than the last training');
    } finally {
      open();
      await closeTrainView(view, [demog]);
    }
  }, {timeout: TIMEOUT});

  test('Validation: collapsed, the schemes of the task, Folds or Test fraction, Split seed with Shuffle and no icons',
    async () => {
      const iris = await openIris();
      const view = await openTrainView(iris);
      const v = view.validationInputs;
      try {
        expectArray(view.groups.map((g) => g.root.querySelector('.forge-group-header label')?.textContent),
          ['Data', 'Preparation', 'Method', 'Output', 'Validation']);
        expect(view.groups[2].root.nextElementSibling === view.outputInputs.group.root &&
          view.groups[3].root.nextElementSibling === v.group.root &&
          view.groups[2].root.contains(view.tuneSettings.root), true,
        'Output and Validation are not right after Method and its Settings');
        expect(v.group.isExpanded, false, 'Validation starts expanded');
        expectAlignedCaptions(view.root, ['Data', 'Method', 'Output', 'Validation']);
        v.group.setExpanded(true);
        expect(v.seed.root.querySelector('.ui-input-options i, .ui-input-plus, .ui-input-minus') === null, true,
          'Split seed has icons');
        expectArray(v.scheme.items, ['K-fold, stratified', 'K-fold', 'Hold-out']);
        expectArray([v.scheme.value, v.folds.value, v.testFraction.value, v.shuffle.value, v.seed.value],
          ['K-fold, stratified', 5, 0.2, true, 42]);
        expect(isShown(v.folds.root) && !isShown(v.testFraction.root), true, 'Folds is hidden or Test fraction shown');
        expect(resultsText(view).includes('Validation: 5-fold cross-validation on 150 rows, stratified'), true,
          resultsText(view));

        v.shuffle.value = false;
        expect(v.seed.enabled, false, 'Split seed is enabled without Shuffle');
        await awaitTrained(view);
        expect(view.lastTraining?.result.splitting.isShuffled, false);
        expect(resultsText(view).includes('Split seed'), false, 'A seed line without Shuffle');
        expect(runsGrid(view).dataFrame.getCol('Split seed').isNone(0), true, 'A split seed without Shuffle');
        expect(headerTooltip(runsGrid(view), 'Split seed').startsWith('The seed that shuffled'), true,
          headerTooltip(runsGrid(view), 'Split seed'));
        v.shuffle.value = true;
        expect(v.seed.enabled, true, 'Split seed stays disabled with Shuffle');

        v.scheme.value = 'Hold-out';
        expect(!isShown(v.folds.root) && isShown(v.testFraction.root), true, 'Folds is shown for Hold-out');
        await awaitTrained(view);
        expect(view.lastTraining?.result.splitting.scheme, 'holdout');
        expect(resultsText(view).includes('Validation: hold out of 20% (30 rows), stratified'), true,
          resultsText(view));
        expect(resultsText(view).includes('Split seed: 42'), true, resultsText(view));
        const heldOut = view.runs[0];
        expect(heldOut.validation.scheme, 'holdout');

        v.scheme.value = 'K-fold, stratified';
        view.featuresInput.value = columnsOf(iris, ['Sepal.Length', 'Sepal.Width', 'Petal.Width']);
        view.targetInput.value = iris.getCol('Petal.Length');
        expectArray(v.scheme.items, ['K-fold', 'Hold-out']);
        expect(v.scheme.value, 'K-fold');
        await awaitTrained(view);
        expect(view.lastTraining?.result.splitting.isStratified, false);
        view.targetInput.value = iris.getCol('Species');
        view.featuresInput.value = columnsOf(iris, MEASUREMENTS);
        expect(v.scheme.value, 'K-fold, stratified');
        await awaitTrained(view);

        v.folds.value = 25;
        await awaitCheck(() => resultsState(view) === 'false/true', 'Folds 25 is not a problem', WAIT_MS);
        v.folds.value = 5;
        await awaitTrained(view);
        view.restoreRun(heldOut);
        expectArray([v.scheme.value, v.seed.value], ['Hold-out', 42]);
        await awaitTrained(view);
        expect(view.lastTraining?.result.splitting.scheme, 'holdout');
      } finally {
        await closeTrainView(view, [iris]);
      }
    }, {timeout: 120000});

  test('the Validation values stay across a Table change', async () => {
    const tables = [grok.shell.addTable(await openIris()), grok.shell.addTable(await openIris())];
    const view = await openTrainView(tables[0]);
    try {
      const v = view.validationInputs;
      v.scheme.value = 'Hold-out';
      v.testFraction.value = 0.3;
      v.shuffle.value = false;
      view.tableInput.value = tables[1];
      const kept = view.validationInputs;
      expect(kept !== v, true, 'The form was not rebuilt');
      expectArray([kept.scheme.value, kept.testFraction.value, kept.shuffle.value, kept.seed.enabled],
        ['Hold-out', 0.3, false, false]);
    } finally {
      await closeTrainView(view, tables);
      for (const table of tables)
        grok.shell.closeTable(table);
    }
  }, {timeout: TIMEOUT});

  test('SVM settings follow the kernel and the task: disabled, kept and sent (EDA\'s meta.mlenabled)', async () => {
    const iris = await openIris();
    const view = await openTrainView(iris);
    const input = (name: string) => {
      const found = view.hyperparameterInputs.get(name);
      if (found === undefined)
        throw new Error(`No ${name} input`);
      return found;
    };
    const enabled = (names: string[]) => names.map((name) => `${name} ${input(name).enabled}`);
    try {
      view.methodInput.value = 'SVM';
      await awaitTrained(view);
      expect(view.hyperparameterInputs.has('seed'), false, 'SVM has a seed');
      expectArray(enabled(['gamma', 'degree', 'coef0', 'epsilon']),
        ['gamma true', 'degree false', 'coef0 false', 'epsilon false']);
      input('kernel').value = 'Polynomial';
      expectArray(enabled(['degree', 'coef0']), ['degree true', 'coef0 true']);
      input('degree').value = 5;
      await awaitTrained(view);
      input('kernel').value = 'Linear';
      expectArray(enabled(['gamma', 'degree', 'coef0']), ['gamma false', 'degree false', 'coef0 false']);
      await awaitTrained(view);
      expect(input('degree').value, 5);
      expect(view.lastTraining?.result.hyperparameters.degree, 5);

      view.featuresInput.value = columnsOf(iris, ['Sepal.Length', 'Sepal.Width', 'Petal.Width']);
      view.targetInput.value = iris.getCol('Petal.Length');
      expectArray(enabled(['epsilon']), ['epsilon true']);
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: 120000});

  test('a server method that does not answer does not hold back the training', async () => {
    const iris = await openIris();
    const view = await TrainView.create(iris, null, withStuck(new Promise(() => {})));
    view.checkTimeoutMs = 3000;
    grok.shell.addView(view);
    try {
      await awaitTrained(view);
      expect(view.methodInput.items.includes(STUCK), false, 'Method lists the method still checking');
      expectArray((await runsOf([iris])).map((r) => r.status), ['completed']);
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: TIMEOUT});

  test('a server method\'s late answer joins Method without retraining', async () => {
    const iris = await openIris();
    const answer = deferred();
    const view = await TrainView.create(iris, null, withStuck(answer.promise));
    grok.shell.addView(view);
    try {
      await awaitTrained(view);
      const shown = shownRun(view);
      answer.resolve(true);
      await waitUntil(() => view.methodInput.items.includes(STUCK), 'The late method does not join Method', WAIT_MS);
      expect(view.isTraining, false, 'The late answer trains');
      expect(shownRun(view), shown, 'The results show another run');
      expect((await runsOf([iris])).length, 1, 'The late answer recorded a run');
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: TIMEOUT});

  test('a method that does not answer is not asked again', async () => {
    const iris = await openIris();
    const calls = {count: 0};
    const view = await TrainView.create(iris, null, withStuck(new Promise(() => {}), calls));
    view.checkTimeoutMs = 300;
    grok.shell.addView(view);
    try {
      await awaitTrained(view);
      await waitUntil(() => view.silentMethods.has(STUCK), `${STUCK} is not marked silent`, WAIT_MS);
      view.featuresInput.value = columnsOf(iris, ['Sepal.Length', 'Petal.Length', 'Petal.Width']);
      await awaitTrained(view);
      expect(calls.count, 1, 'The silent method was asked again');
      expect(balloonWith('did not answer') === null, true, 'A balloon names the silent method');
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: TIMEOUT});

  test('a chosen server method waits for its answer', async () => {
    const iris = await openIris();
    const answer = deferred();
    const view = await TrainView.create(iris, STUCK, withStuck(answer.promise));
    grok.shell.addView(view);
    try {
      await waitUntil(() => resultsState(view) === 'true/false' &&
        resultsText(view).includes(`Waiting for ${STUCK} to check the data...`),
      () => `Results: ${resultsState(view)} ${resultsText(view)}`, WAIT_MS);
      expect(view.methodInput.value, STUCK);
      expect(view.lastTraining === undefined, true, 'The view trained before the method answered');
      expect((await runsOf([iris])).length, 0, 'A run was recorded before the method answered');
      answer.resolve(false);
      await waitUntil(() => balloonWith(`${STUCK} cannot learn from this selection; XGBoost is chosen for now.`) !==
        null, 'No balloon for the refusing method', WAIT_MS);
      await awaitTrained(view);
      expect(view.methodInput.value, 'XGBoost');
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: TIMEOUT});
});

category('UI: Preparation', () => {
  test('Preparation follows demog: Missing values, One-hot encoding, Skip unique categories', async () => {
    const full = await openDemog();
    const demog = spreadRows(full, 1000, (i) => full.getCol('HEIGHT').isNone(i) || full.getCol('WEIGHT').isNone(i));
    expect(demog.rowCount <= 1000 && demog.getCol('HEIGHT').stats.missingValueCount +
      demog.getCol('WEIGHT').stats.missingValueCount > 0, true, 'The subset lost the gaps');
    const view = await TrainView.create(demog);
    grok.shell.addView(view);
    const preparation = view.preparationInputs;
    const {oneHot, skipUniqueCategories} = preparation;
    const {predictProbability, threshold} = view.outputInputs;
    try {
      view.methodInput.value = 'SVM';
      await awaitTrained(view);
      expect(isShown(preparation.group.root) && isShown(view.missingValuesInputs.choice.root), true,
        'Preparation does not show Missing values for HEIGHT and WEIGHT');
      expect([oneHot, skipUniqueCategories, threshold].some((input) => isShown(input.root)), false,
        'A step is shown for numerical features and a many-class target');
      expect(isShown(view.outputInputs.group.root) && isShown(predictProbability.root), true,
        'Predict probability is not offered for a many-class target');
      const forms = view.root.querySelectorAll('.ui-form');
      expect(forms.length === 1 && forms[0].contains(preparation.group.root), true, 'Preparation is not in the form');

      view.targetInput.value = demog.getCol('AGE');
      view.featuresInput.value = columnsOf(demog, ['HEIGHT', 'WEIGHT', 'SEX', 'RACE']);
      await awaitCheck(() => isShown(oneHot.root), 'One-hot encoding is not offered', WAIT_MS);
      expect(oneHot.value, true, 'One-hot encoding is not checked for SEX and RACE');
      await awaitTrained(view);
      const encoded = view.lastTraining?.result.options;
      expectArray(encoded?.preprocessingInfo ?? [], ['ignore-missing', 'one-hot']);
      expectArray(Object.keys(encoded?.oneHotCategories ?? {}), ['SEX', 'RACE']);
      expect(isShown(skipUniqueCategories.root) || isShown(predictProbability.root), false,
        'A step that does not apply is shown');
      expect(isShown(view.outputInputs.group.root), false, 'Output is shown for a numerical target');

      oneHot.value = false;
      await awaitCheck(() => view.featuresInput.validity !== null, 'Text features are not marked', TIMEOUT);
      expect(view.featuresInput.validity, `${view.methodInput.value} needs numerical features. Check One-hot ` +
        'encoding in Preparation, or uncheck: SEX, RACE.');

      oneHot.value = true;
      await awaitCheck(() => view.featuresInput.validity === null, 'One-hot encoding leaves Features marked', WAIT_MS);
      await awaitTrained(view);

      view.featuresInput.value = [...view.featuresInput.value, demog.getCol('USUBJID')];
      await awaitCheck(() => isShown(skipUniqueCategories.root), 'Skip unique categories is not offered', WAIT_MS);
      expect(skipUniqueCategories.value, true, 'Skip unique categories is not checked when it appears');
      await awaitTrained(view);
      expectArray(view.lastTraining?.result.options.preprocessingInfo ?? [],
        ['ignore-missing', 'skip-unique-categories', 'one-hot']);
      expectArray(view.lastTraining?.result.features.columns.map((c) => c.name).sort() ?? [],
        ['HEIGHT', 'RACE', 'SEX', 'USUBJID', 'WEIGHT']);

      view.featuresInput.value = columnsOf(demog, ['HEIGHT', 'WEIGHT', 'SEX', 'RACE']);
      await awaitCheck(() => !isShown(skipUniqueCategories.root), 'Skip unique categories stays', WAIT_MS);
      view.featuresInput.value = columnsOf(demog, ['HEIGHT', 'WEIGHT']);
      await awaitCheck(() => !isShown(oneHot.root), 'One-hot encoding stays without text features', WAIT_MS);
      expect(oneHot.value, false, 'A hidden One-hot encoding is not unchecked');
    } finally {
      await closeTrainView(view, [demog]);
    }
  }, {timeout: 120000});

  test('One-hot encoding is checked by default only for few categories, and keeps the user\'s choice', async () => {
    const rows = 60;
    const text = (name: string, value: (i: number) => string) =>
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, name, Array.from({length: rows}, (_, i) => value(i)));
    const table = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'x', valuesOf(rows, (i) => i)),
      text('few', (i) => i % 2 === 0 ? 'a' : 'b'), text('many', (i) => `k${i % 30}`), text('id', (i) => `id${i}`),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'y', valuesOf(rows, (i) => i * 1.5 + i % 7)),
    ]);
    table.name = `forge-test-categories-${Date.now()}`;
    const view = await TrainView.create(table);
    grok.shell.addView(view);
    const {oneHot, skipUniqueCategories} = view.preparationInputs;
    const use = async (names: string[], isChecked: boolean) => {
      view.featuresInput.value = columnsOf(table, names);
      await awaitCheck(() => isShown(oneHot.root) && oneHot.value === isChecked,
        `One-hot encoding is not ${isChecked ? '' : 'un'}checked for ${names.join(', ')}`, WAIT_MS);
    };
    try {
      await use(['x', 'few'], true);
      await use(['x', 'few', 'many'], false);
      await use(['x', 'few', 'id'], true);
      skipUniqueCategories.value = false;
      await awaitCheck(() => !oneHot.value, 'An unskipped id leaves One-hot encoding checked', WAIT_MS);
      skipUniqueCategories.value = true;
      await awaitCheck(() => oneHot.value, 'Skipping the id does not check One-hot encoding', WAIT_MS);

      await use(['x', 'many'], false);
      oneHot.value = true;
      await use(['x', 'few', 'many'], true);
      view.featuresInput.value = columnsOf(table, ['x']);
      await awaitCheck(() => !isShown(oneHot.root), 'One-hot encoding stays without text features', WAIT_MS);
      expect(oneHot.value, false, 'A hidden One-hot encoding is not unchecked');
      await use(['x', 'few', 'many'], false);
    } finally {
      await closeTrainView(view, [table]);
    }
  }, {timeout: 120000});

  test('Train keeps Predict probability while the form is locked', async () => {
    const iris = await twoSpeciesIris();
    const view = await openTrainView(iris);
    const {predictProbability} = view.outputInputs;
    try {
      view.methodInput.value = 'SVM';
      await awaitCheck(() => view.lastTraining?.engine.name === 'SVM' && !view.isTraining, 'No SVM training', TIMEOUT);
      expect(view.lastTraining?.request.hyperparameters.probability, true, 'SVM trained without probabilities');
      const count = view.runs.length;
      predictProbability.value = true;
      await waitUntil(() => view.lastTraining?.result.options.postprocessingInfo.includes('class-probabilities') ===
        true, 'Predict probability does not re-cut SVM', WAIT_MS);
      await delay(CHECK_DELAY_MS * 3);
      expect(view.runs.length, count, 'Predict probability retrained SVM');
      const live = view.lastTraining;

      await view.train();
      const training = view.lastTraining!;
      expect(training !== live && training.engine.name === 'SVM', true, 'Train gave no new SVM training');
      expect(predictProbability.value && predictProbability.enabled, true, 'Predict probability is off after Train');
      expect(training.request.hyperparameters.probability, true, 'SVM trained without probabilities');
      expect(training.result.hyperparameters.probability, true);
      expect(training.result.options.postprocessingInfo.includes('class-probabilities'), true,
        'The result has no class-probabilities step');
      expect(view.runs.length, count + 1, 'No new run');
      expect(view.runs[0].steps.predictProbability, true, 'The run row lost Predict probability');
      expect(view.runs[0].hyperparameters.probability, true);
      const record = await forgeDb.trainingRuns.get(training.runId);
      expect((record.options?.postprocessingInfo as string[] | undefined)?.includes('class-probabilities'), true,
        'The training run was recorded without the step');
      expect(record.hyperparameters?.probability, true);
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: 120000});
});

/** The Results rows of [view], in order. */
const rowsOf = (view: TrainView | EvaluationView) =>
  Array.from(view.resultsPane.root.querySelectorAll<HTMLElement>('.forge-results-row'));

const titleOf = (panel: Element) => (panel.firstElementChild?.textContent ?? '').trim();

/** The panel boxes of a Results row. */
const slotsOf = (row: HTMLElement) => Array.from(row.children)
  .filter((c): c is HTMLElement => c instanceof HTMLElement && !c.classList.contains('ui-split-h-divider'));

/** The panel titles of each Results row of [view], in order. */
const rowTitles = (view: TrainView | EvaluationView) =>
  rowsOf(view).map((row) => Array.from(row.querySelectorAll('.forge-results-panel'), titleOf));

const expectRows = (actual: unknown, expected: unknown) => expect(JSON.stringify(actual), JSON.stringify(expected));

/** The viewer of the Results panel titled [title] of [view]. */
function panelViewer(view: TrainView, title: string): DG.Viewer {
  const element = panelOf(view, title).querySelector('.d4-viewer');
  const viewer: unknown = element instanceof HTMLElement ? DG.toJs(DG.Widget.find(element)) : null;
  if (!(viewer instanceof DG.Viewer))
    throw new Error(`No viewer in the panel ${title}`);
  return viewer;
}

/** The text of each line of Insights & Tips of [view], its pointer line included. */
const insightLines = (view: TrainView | EvaluationView) =>
  Array.from(view.resultsPane.root.querySelectorAll('.forge-insights > ul > li'), (li) => li.textContent ?? '');

/** The line of Insights & Tips of [view] whose text starts with [text]. */
function insightItem(view: TrainView, text: string): HTMLElement {
  const item = Array.from(view.resultsPane.root.querySelectorAll<HTMLElement>('.forge-insights > ul > li'))
    .find((li) => (li.textContent ?? '').startsWith(text));
  if (item === undefined)
    throw new Error(`No insight ${text}: ${insightLines(view).join(' | ')}`);
  return item;
}

/** The panel titled [title] of [view]. */
function panelOf(view: TrainView, title: string): HTMLElement {
  const panel = Array.from(view.resultsPane.root.querySelectorAll<HTMLElement>('.forge-results-panel'))
    .find((p) => titleOf(p) === title);
  if (panel === undefined)
    throw new Error(`No panel ${title}`);
  return panel;
}

/** Checks that in every panel of [view]'s rows the title is shown, the content starts under it and ends inside the
 * panel, which does not scroll (its content does). */
function expectUnderTitles(view: TrainView): void {
  for (const panel of view.resultsPane.root.querySelectorAll<HTMLElement>('.forge-results-row .forge-results-panel')) {
    const title = panel.firstElementChild?.getBoundingClientRect();
    const content = panel.children[1]?.getBoundingClientRect();
    const box = panel.getBoundingClientRect();
    expect(title !== undefined && content !== undefined && title.height > 30 && content.top >= title.bottom - 0.5 &&
      content.bottom <= box.bottom + 0.5 && panel.scrollHeight <= panel.clientHeight + 1, true,
    `${titleOf(panel)}: title ${title?.top}-${title?.bottom}, content ${content?.top}-${content?.bottom}, panel ` +
      `${box.top}-${box.bottom}, scroll ${panel.scrollHeight} of ${panel.clientHeight}`);
  }
}

/** The row of the metric [label] in a metrics grid. */
const metricRow = (grid: DG.Grid, label: string) => grid.dataFrame.getCol('Metric').toList().indexOf(label);

/** The title of [view]'s PR Curve panel (`PR Curve (AP = ...)`); '' when there is none. */
const prTitleOf = (view: TrainView | EvaluationView) =>
  rowTitles(view).flat().find((t) => t.startsWith('PR Curve')) ?? '';

/** The row of a curve frame's `Threshold` that [threshold] labels: the last at or above it ([fromEnd], the ROC's
 * decreasing thresholds), else the first (the PR's increasing ones); -1 for none. */
function markRow(frame: DG.DataFrame, threshold: number, fromEnd = true): number {
  const thresholds = frame.getCol('Threshold').toList();
  const rows = thresholds.map((_, i) => i).filter((i) => thresholds[i] >= threshold);
  return rows.length === 0 ? -1 : fromEnd ? rows[rows.length - 1] : rows[0];
}

/** [viewer]'s current point is [row] of its frame, the only `cutoff` mark (circle) of its `Mark` column. The canvas
 * itself is not checked. */
function expectMarked(viewer: DG.Viewer, row: number, name: string): void {
  expect(viewer.dataFrame.currentRowIdx, row, `${name}: the threshold is not marked`);
  expect(viewer.props.showCurrentPoint, true, `${name}: no current point`);
  expect(viewer.props.markersColumnName, 'Mark', `${name}: no markers column`);
  const marks = viewer.dataFrame.getCol('Mark').toList();
  expect(marks.filter((m) => m === 'cutoff').length, row < 0 ? 0 : 1, `${name}: not one circle`);
  if (row >= 0)
    expect(marks[row], 'cutoff', `${name}: the circle is not at the threshold`);
}

/** [viewer] draws a 3 px line on automatic axes. */
function expectCurveStyle(viewer: DG.Viewer, name: string): void {
  expect(viewer.props.linesWidth, 3, `${name}: line width`);
  const {xMin, xMax, yMin, yMax} = viewer.props;
  expect([xMin, xMax, yMin, yMax].every((v) => v == null), true, `${name}: a fixed axis range`);
}

/** [frame]'s columns as lists, to tell a changed curve. */
const curveValues = (frame: DG.DataFrame) => JSON.stringify(frame.columns.toList().map((c) => c.toList()));

/** Predict probability on two iris species: ROC Curve and PR Curve over their own frames of Crux's curves, kept by a
 * threshold change (the mark moves), given new frames by a retrain. */
async function expectRocCurve(): Promise<void> {
  const iris = await twoSpeciesIris();
  const view = await openTrainView(iris);
  const {predictProbability, threshold} = view.outputInputs;
  try {
    predictProbability.value = true;
    const rocTitle = () => rocTitleOf(view);
    const prTitle = () => prTitleOf(view);
    const hasStep = () => view.lastTraining?.result.options.postprocessingInfo.includes('class-probabilities') === true;
    await awaitCheck(() => rocTitle() !== '' && hasStep() && !view.isTraining, 'No ROC Curve', TIMEOUT);
    await delay(CHECK_DELAY_MS * 3);
    const aucTitle = () => `ROC Curve (AUC = ${view.lastTraining?.result.metrics.train.auc?.toFixed(3)})`;
    const apTitle = () =>
      `PR Curve (AP = ${view.lastTraining?.result.metrics.train.averagePrecision?.toFixed(3)})`;
    expect(rocTitle(), aucTitle());
    expect(prTitle(), apTitle());
    expectRows(rowTitles(view)[0], ['Predicted Species vs Actual', rocTitle(), prTitle()]);
    expect(rowTitles(view)[1].includes('Performance'), true, 'Performance is not in row 2');
    const roc = panelViewer(view, rocTitle());
    const frame = roc.dataFrame;
    expect(roc.type, DG.VIEWER.SCATTER_PLOT);
    const [fprName, recallName] = ['FPR (1 − Specificity)', 'Recall (Sensitivity)'];
    expectArray([roc.props.xColumnName, roc.props.yColumnName, roc.props.linesOrderColumnName],
      [fprName, recallName, 'Threshold']);
    expectArray(frame.columns.names(), ['Threshold', fprName, recallName, 'Mark']);
    expectCurveStyle(roc, 'ROC Curve');
    expect(frame.rowCount, view.lastTraining?.result.curves?.roc.fpr.length);
    expect(view.resultsPane.frames.some((f) => f.dart === frame.dart), false, 'ROC Curve is on the preview frame');
    expect(frame.meta.formulaLines.items.some((line) => line.title === 'Non-informative baseline'), true,
      'No baseline');
    const optimized = view.lastTraining!.result.options.binaryClassificationThreshold!;
    expect(threshold.value, optimized, 'The threshold is not the optimized one');
    expectMarked(roc, markRow(frame, optimized), 'ROC Curve');
    const pr = panelViewer(view, prTitle());
    const prFrame = pr.dataFrame;
    expectArray(prFrame.columns.names(), ['Threshold', 'Recall (Sensitivity)', 'Precision', 'Mark']);
    expectCurveStyle(pr, 'PR Curve');
    const base = prFrame.meta.formulaLines.items.find((l) => l.title === 'Non-informative baseline');
    expect(base?.formula, `\${Precision} = ${view.lastTraining?.result.curves?.pr.precision[0].toFixed(4)}`);
    expectMarked(pr, markRow(prFrame, optimized, false), 'PR Curve');
    const content = (title: string) => panelOf(view, title).lastElementChild;
    const confusions = content('Confusions');
    expect(panelViewer(view, 'Confusions').props.descriptionPosition, 'Bottom');
    const points = () => JSON.stringify(frame.columns.toList().filter((c) => c.name !== 'Mark').map((c) => c.toList()));
    const pointsBefore = points();
    threshold.value = 0.9;
    await awaitCheck(() => view.lastTraining?.result.options.binaryClassificationThreshold === 0.9 &&
      !view.isTraining, 'The threshold does not re-cut the result', WAIT_MS);
    await delay(CHECK_DELAY_MS * 3);
    const recut = panelViewer(view, rocTitle());
    expect(recut.root === roc.root && roc.root.isConnected && recut.dataFrame.dart === frame.dart, true,
      'The threshold rebuilt ROC Curve');
    expect(points() === pointsBefore, true, 'The threshold changed the curve');
    expect(frame.currentRowIdx, markRow(frame, 0.9), 'The mark did not move to the threshold 0.9');
    expect(prFrame.currentRowIdx, markRow(prFrame, 0.9, false), 'The PR mark did not move to the threshold 0.9');
    expect(content('Confusions') === confusions, true, 'The threshold rebuilt Confusions');

    const values = curveValues(frame);
    const before = view.lastTraining;
    if (before?.engine.name === 'XGBoost')
      setHyperparameter(view, 'iterations', 1);
    else
      view.methodInput.value = 'XGBoost';
    await awaitCheck(() => view.lastTraining !== before && view.lastTraining?.engine.name === 'XGBoost' &&
      !view.isTraining, 'No XGBoost retrain', TIMEOUT);
    await delay(CHECK_DELAY_MS * 3);
    const retrained = panelViewer(view, rocTitle());
    const newFrame = retrained.dataFrame;
    expect(retrained.root === roc.root && newFrame.dart !== frame.dart, true,
      'The retrain rebuilt ROC Curve or kept its frame');
    expectArray([retrained.props.xColumnName, retrained.props.yColumnName, retrained.props.linesOrderColumnName],
      [fprName, recallName, 'Threshold']);
    expectCurveStyle(retrained, 'ROC Curve after a retrain');
    expect(newFrame.rowCount, view.lastTraining?.result.curves?.roc.fpr.length);
    expect(rocTitle(), aucTitle());
    expect(prTitle(), apTitle());
    expect(curveValues(newFrame) !== values, true, 'The retrain left the curve as it was');
    expectMarked(retrained, markRow(newFrame, 0.9), 'ROC Curve after a retrain');
    const newPr = panelViewer(view, prTitle());
    expectCurveStyle(newPr, 'PR Curve after a retrain');
    expectMarked(newPr, markRow(newPr.dataFrame, 0.9, false), 'PR Curve after a retrain');
  } finally {
    await closeTrainView(view, [iris]);
  }
}

/** XGBoost on two iris species: ROC Curve and PR Curve from the class probabilities, no mark, new frames on a
 * retrain; SVM has them too (always trained with class probabilities); none for a regression or three classes, which
 * get the AUC-ROC and Average precision rows and the Per-class parallel coordinates plot, refilled by a retrain. */
async function expectProbabilityRocCurve(): Promise<void> {
  const iris = await twoSpeciesIris();
  const view = await openTrainView(iris);
  const rocTitle = () => rocTitleOf(view);
  const prTitle = () => prTitleOf(view);
  let other: TrainView | undefined;
  let threeSpecies: DG.DataFrame | undefined;
  try {
    await awaitCheck(() => rocTitle() !== '' && !view.isTraining, 'No ROC Curve', TIMEOUT);
    expectRows(rowTitles(view)[0], ['Predicted Species vs Actual', rocTitle(), prTitle()]);
    expect(rowTitles(view)[1].includes('Performance'), true, 'Performance is not in row 2');
    expect(view.lastTraining?.result.outOfFold.probabilities !== undefined, true, 'No out-of-fold probabilities');
    expect(rocTitle(), `ROC Curve (AUC = ${view.lastTraining?.result.metrics.train.auc?.toFixed(3)})`);
    const roc = panelViewer(view, rocTitle());
    const frame = roc.dataFrame;
    expectArray(frame.columns.names(), ['Threshold', 'FPR (1 − Specificity)', 'Recall (Sensitivity)', 'Mark']);
    expect(frame.rowCount, view.lastTraining?.result.curves?.roc.fpr.length);
    expect(frame.currentRowIdx, -1, 'A curve without a threshold has a marked point');
    expect(frame.getCol('Mark').toList().includes('cutoff'), false, 'A curve without a threshold has a circle');
    expect(panelViewer(view, prTitle()).dataFrame.currentRowIdx, -1, 'A PR curve without a threshold has a mark');
    expect(view.resultsPane.frames.some((f) => f.dart === frame.dart), false, 'ROC Curve is on the preview frame');
    const values = curveValues(frame);

    const before = view.lastTraining;
    setHyperparameter(view, 'iterations', 1);
    await awaitCheck(() => view.lastTraining !== before && !view.isTraining, 'No retrain', TIMEOUT);
    await delay(CHECK_DELAY_MS * 3);
    const retrained = panelViewer(view, rocTitle());
    expect(retrained.root === roc.root && retrained.dataFrame.dart !== frame.dart, true,
      'The retrain rebuilt ROC Curve or kept its frame');
    expect(curveValues(retrained.dataFrame) !== values, true, 'The retrain left the curve as it was');

    view.methodInput.value = 'SVM';
    await awaitCheck(() => view.lastTraining?.engine.name === 'SVM' && !view.isTraining, 'No SVM training', TIMEOUT);
    const svmAuc = `ROC Curve (AUC = ${view.lastTraining?.result.metrics.train.auc?.toFixed(3)})`;
    await waitUntil(() => rocTitle() === svmAuc && !view.isTraining, () => `SVM's ROC Curve: '${rocTitle()}', ` +
      `expected '${svmAuc}'`, WAIT_MS);
    expect(prTitle().startsWith('PR Curve (AP = '), true, prTitle());
    expect(hasMetric(view.root, 'AUC-ROC'), true, 'SVM has no AUC-ROC row');
    expect(hasMetric(view.root, 'Average precision'), true, 'SVM has no Average precision row');
    expect(hasMetric(view.root, 'MCC'), true, 'SVM has no MCC row');

    view.methodInput.value = 'XGBoost';
    view.featuresInput.value = columnsOf(iris, ['Sepal.Length', 'Sepal.Width']);
    view.targetInput.value = iris.getCol('Petal.Length');
    await awaitCheck(() => view.lastTraining?.result.task === 'regression' && !view.isTraining,
      'No regression training', TIMEOUT);
    await delay(CHECK_DELAY_MS * 3);
    expect(rocTitle(), '', 'A regression has a ROC Curve');

    threeSpecies = await openIris();
    other = await openTrainView(threeSpecies);
    expect(rocTitleOf(other), '', 'Three classes have ROC Curve');
    expect(hasMetric(other.root, 'AUC-ROC'), true, 'Three classes have no AUC-ROC row');
    expect(hasMetric(other.root, 'Average precision'), true, 'Three classes have no Average precision row');
    const trained = other;
    await awaitCheck(() => rowTitles(trained)[0]?.includes('Per-class') === true, 'No Per-class panel', WAIT_MS);
    expectRows(rowTitles(other)[0], ['Predicted Species vs Actual', 'Confusions', 'Per-class']);
    const pcPlot = panelViewer(other, 'Per-class');
    expect(pcPlot.type, DG.VIEWER.PC_PLOT);
    const perClass = pcPlot.dataFrame;
    expect(perClass.rowCount, 3);
    expectArray(perClass.columns.names(), ['Class', 'Precision', 'Recall', 'F1 score', 'AUC-ROC']);
    expectArray(pcPlot.props.columnNames, ['Precision', 'Recall', 'F1 score', 'AUC-ROC']);
    expect(pcPlot.props.colorColumnName, 'Class');
    expect(pcPlot.props.legendVisibility, 'Always', 'The three species have no legend');
    expect(pcPlot.props.legendPosition, 'Top', 'The three species\' legend is not on top');
    expect(pcPlot.props.showFilters, false, 'Per-class shows its in-chart filters');
    expect(pcPlot.props.showColorSelector, false, 'Per-class shows its color selector');
    const target = other.lastTraining!.request.target;
    const texts = target.toList().map(String);
    for (const [k, label] of perClass.getCol('Class').toList().entries()) {
      expect(DG.Color.getRowColor(perClass.getCol('Class'), k), DG.Color.getRowColor(target, texts.indexOf(label)),
        `${label}: another color than the target's`);
    }
    expect(other.resultsPane.frames.some((f) => f.dart === perClass.dart), false, 'Per-class is on the preview frame');
    const perClassValues = curveValues(perClass);
    const beforeRetrain = other.lastTraining;
    setHyperparameter(other, 'iterations', 1);
    await awaitCheck(() => other!.lastTraining !== beforeRetrain && !other!.isTraining, 'No retrain', TIMEOUT);
    await delay(CHECK_DELAY_MS * 3);
    const refilled = panelViewer(other, 'Per-class');
    expect(refilled.root === pcPlot.root && refilled.dataFrame.dart === perClass.dart, true,
      'The retrain rebuilt Per-class');
    expect(curveValues(perClass) !== perClassValues, true, 'The retrain left the per-class values as they were');
    other.methodInput.value = 'SVM';
    await awaitCheck(() => other!.lastTraining?.engine.name === 'SVM' && !other!.isTraining, 'No SVM training',
      TIMEOUT);
    await delay(CHECK_DELAY_MS * 3);
    expectArray(panelViewer(other, 'Per-class').dataFrame.columns.names(),
      ['Class', 'Precision', 'Recall', 'F1 score', 'AUC-ROC']);
  } finally {
    if (other !== undefined && threeSpecies !== undefined)
      await closeTrainView(other, [threeSpecies]);
    await closeTrainView(view, [iris]);
  }
}

/** XGBoost on demog SEX: ROC Curve from the class probabilities without a mark; Predict probability on marks the
 * threshold, off drops the mark again; the metrics grid is kept throughout. */
async function expectProbabilityToggle(): Promise<void> {
  const demog = await openDemog();
  const view = await TrainView.create(demog);
  grok.shell.addView(view);
  const {predictProbability} = view.outputInputs;
  const rocTitle = () => rocTitleOf(view);
  const rocFrame = () => panelViewer(view, rocTitle()).dataFrame;
  try {
    view.targetInput.value = demog.getCol('SEX');
    view.featuresInput.value = columnsOf(demog, ['AGE', 'HEIGHT', 'WEIGHT']);
    await awaitCheck(() => view.lastTraining?.result.probabilities !== undefined && !view.isTraining &&
      rocTitle() !== '', 'No ROC Curve from the class probabilities', TIMEOUT);
    expect(rocFrame().currentRowIdx, -1, 'A curve without a threshold has a marked point');
    const grid = gridWith(view.root, 'Metric').root;
    const metricRows = () => JSON.stringify(gridWith(view.root, 'Metric').dataFrame.getCol('Metric').toList());
    const rows = metricRows();

    predictProbability.value = true;
    const hasStep = () => view.lastTraining?.result.options.postprocessingInfo.includes('class-probabilities') === true;
    await awaitCheck(() => hasStep() && !view.isTraining && rocTitle() !== '',
      'No ROC Curve with Predict probability', TIMEOUT);
    await delay(CHECK_DELAY_MS * 3);
    const shown = view.lastTraining!.result.options.binaryClassificationThreshold!;
    expectMarked(panelViewer(view, rocTitle()), markRow(rocFrame(), shown), 'ROC Curve');
    const pr = panelViewer(view, prTitleOf(view));
    expectMarked(pr, markRow(pr.dataFrame, shown, false), 'PR Curve');
    expect(gridWith(view.root, 'Metric').root === grid, true, 'Predict probability built a new metrics grid');
    expect(metricRows(), rows, 'Predict probability changed the metric rows');

    predictProbability.value = false;
    await awaitCheck(() => view.lastTraining?.result.probabilities !== undefined && !view.isTraining &&
      rocTitle() !== '', 'No ROC Curve after Predict probability is off', TIMEOUT);
    await delay(CHECK_DELAY_MS * 3);
    expect(rocFrame().currentRowIdx, -1, 'The mark stayed without Predict probability');
    expect(rocFrame().getCol('Mark').toList().includes('cutoff'), false,
      'The circle stayed without Predict probability');
    expect(gridWith(view.root, 'Metric').root === grid, true, 'Turning Predict probability off built a new grid');
    expect(metricRows(), rows, 'Turning Predict probability off changed the metric rows');
  } finally {
    await closeTrainView(view, [demog]);
  }
}

/** Run Evaluation of a two-class XGBoost model without Predict probability: ROC Curve from the class probabilities,
 * no mark, the AUC-ROC row with an Evaluation value. */
async function expectEvaluationRocCurve(): Promise<void> {
  const iris = await twoSpeciesIris();
  const {id} = await saveTestModel(columnsOf(iris, MEASUREMENTS), iris.getCol('Species'), iris.name);
  let view: EvaluationView | undefined;
  try {
    const model = await forgeDb.models.get(id);
    const loaded = await loadModel(id);
    view = await EvaluationView.open(model, loaded, iris, iris.getCol('Species'), exactMapping(loaded.features, iris));
    const titles = rowTitles(view)[0];
    const roc = rocTitleOf(view);
    const pr = prTitleOf(view);
    expectRows(titles, ['Predicted Species vs Actual', roc, pr]);
    expect(rowTitles(view)[1].includes('Performance'), true, 'Performance is not in row 2');
    const metrics = gridWith(view.resultsPane.root, 'Metric').dataFrame;
    const aucRow = metricRow(gridWith(view.resultsPane.root, 'Metric'), 'AUC-ROC');
    expect(roc, `ROC Curve (AUC = ${metrics.get('Evaluation', aucRow).toFixed(3)})`);
    const apRow = metricRow(gridWith(view.resultsPane.root, 'Metric'), 'Average precision');
    expect(pr, `PR Curve (AP = ${metrics.get('Evaluation', apRow).toFixed(3)})`);
  } finally {
    view?.close();
    await deleteModel(id);
  }
}

/** A viewer that throws while it is built, then the Wrong predictions grid: each time the frames go back and no result
 * is shown. */
async function expectFailedBuildReleased(): Promise<void> {
  // One iteration of depth 1 gets rows wrong.
  const unbalanced = await unbalancedIris();
  const viewerFailure = 'forge-test-viewer-failure';
  const gridFailure = 'forge-test-grid-failure';
  const {confusionMatrix, grid} = DG.Viewer;
  let view: TrainView | undefined;
  try {
    DG.Viewer.confusionMatrix = () => {
      throw new ForgeError(viewerFailure);
    };
    const [opened, built] = await framesSharing(unbalanced.columns.toList(), async () => {
      const trainView = await openTrainView(unbalanced);
      await awaitCheck(() => balloonWith(viewerFailure) !== null, 'The failing viewer was not reported', WAIT_MS);
      return trainView;
    });
    view = opened;
    expect(view.resultsPane.frames.length, 0, 'A result was shown without its viewers');
    expectReleased(built);

    DG.Viewer.confusionMatrix = confusionMatrix;
    // Only the frame of Wrong predictions, every column of the table and the prediction, gets a grid with both.
    DG.Viewer.grid = (df, options) => {
      if (df.col('Species') !== null && df.col('Predicted Species') !== null)
        throw new ForgeError(gridFailure);
      return grid(df, options);
    };
    const trainView = view;
    const [, rebuilt] = await framesSharing(unbalanced.columns.toList(), async () => {
      setHyperparameter(trainView, 'iterations', 1);
      setHyperparameter(trainView, 'maxDepth', 1);
      await awaitTrained(trainView);
      await awaitCheck(() => balloonWith(gridFailure) !== null, 'The failing grid was not reported', WAIT_MS);
    });
    expect(view.resultsPane.frames.length, 0, 'A result was shown without Wrong predictions');
    expectReleased(rebuilt);
  } finally {
    DG.Viewer.confusionMatrix = confusionMatrix;
    DG.Viewer.grid = grid;
    if (view !== undefined)
      await closeTrainView(view, [unbalanced]);
  }
}

category('UI: Results', () => {
  test('resultsLayout: the template rows, the absent panels dropped, the method\'s rows', async () => {
    const all = new Set<PanelKey>(['predictedVsActual', 'rocCurve', 'prCurve', 'perClass', 'metrics', 'confusions',
      'insights', 'wrongPredictions', 'features', 'correlation', 'classBalance']);
    const binary = resultsLayout('binary', all, 0);
    expectRows(binary.map((r) => r.panels), [['predictedVsActual', 'rocCurve', 'prCurve'],
      ['insights', 'confusions', 'metrics'], ['features', 'correlation', 'classBalance'], ['wrongPredictions']]);
    expectArray(binary.map((r) => r.height), [350, 350, 450, 300]);
    expectRows(binary.map((r) => r.weights), [[1, 1, 1], [1, 1.5, 1.5], [1, 1, 1], [1]]);
    const multiclass = resultsLayout('multiclass', new Set<PanelKey>(['predictedVsActual', 'features', 'correlation',
      'classBalance']), 0);
    expectRows(multiclass.map((r) => r.panels), [['predictedVsActual'], ['insights', 'metrics', 'classBalance'],
      ['features', 'correlation']]);
    expectRows(multiclass.map((r) => r.weights), [[1], [1, 1.5, 1.5], [1, 1]]);
    const regression = resultsLayout('regression', new Set<PanelKey>(['predictedVsActual', 'residuals',
      'distributions', 'statistics']), 4);
    expectRows(regression.map((r) => r.panels), [['predictedVsActual', 'residuals', 'distributions'],
      ['insights', 'metrics', 'statistics'], ['visualize:0', 'visualize:1', 'visualize:2'], ['visualize:3']]);
    expectArray(regression.map((r) => r.height), [350, 350, 350, 350]);
    expectRows(regression.map((r) => r.weights), [[1, 1, 1], [1, 1.5, 1.5], [1, 1, 1], [1]]);
    expectRows(resultsLayout('binary', new Set<PanelKey>(['metrics']), 0).map((r) => r.panels),
      [['insights', 'metrics']]);
  });

  test('resultsLayout: Insights & Tips is always in row 2, Class balance in its template row', async () => {
    const withoutWrong = new Set<PanelKey>(['predictedVsActual', 'confusions', 'insights', 'features', 'correlation',
      'classBalance']);
    expectRows(resultsLayout('multiclass', withoutWrong, 0).map((r) => r.panels), [['predictedVsActual', 'confusions'],
      ['insights', 'metrics', 'classBalance'], ['features', 'correlation']]);
    expectRows(resultsLayout('binary', withoutWrong, 0).map((r) => r.panels), [['predictedVsActual'],
      ['insights', 'confusions', 'metrics'], ['features', 'correlation', 'classBalance']]);
    const full = new Set<PanelKey>([...withoutWrong, 'wrongPredictions', 'perClass']);
    expectRows(resultsLayout('multiclass', full, 0).map((r) => r.panels), [['predictedVsActual', 'confusions',
      'perClass'], ['insights', 'metrics', 'classBalance'], ['features', 'correlation'], ['wrongPredictions']]);
  });

  test('a hovered row is linked through the kept rows', async () => {
    const iris = await openIris();
    iris.getCol('Sepal.Width').set(3, null);
    const view = await openTrainView(iris);
    try {
      view.methodInput.value = 'SVM';
      await awaitTrained(view);
      await awaitCheck(() => view.resultsPane.frames[0]?.rowCount === 149, 'The results show every row', WAIT_MS);
      const frame = view.resultsPane.frames[0];
      frame.mouseOverRowIdx = 5;
      // The kept rows are 0, 1, 2, 4, 5, 6, ...: the 6th is row 6.
      await awaitCheck(() => iris.mouseOverRowIdx === 6, `The table hovers row ${iris.mouseOverRowIdx}`, WAIT_MS);
      iris.mouseOverRowIdx = 3;
      await awaitCheck(() => frame.mouseOverRowIdx === -1, `The frame hovers row ${frame.mouseOverRowIdx}`, WAIT_MS);
      expect(iris.mouseOverRowIdx, 3, 'The skipped row came back as none');
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: TIMEOUT});

  test('a retrain keeps the grid, the rows and the viewers; a Features change rebuilds the viewers; a hint ' +
    'releases the frames and keeps the grid', async () => {
    const iris = await openIris();
    const view = await openTrainView(iris);
    try {
      const grid = gridWith(view.resultsPane.root, 'Metric').root;
      const rows = rowsOf(view);
      const titles = rowTitles(view);
      for (const row of rows) {
        expect(row.style.height !== '' && slotsOf(row).every((s) => s.style.flexGrow !== ''), true,
          `A row or a slot without a size: ${row.style.height}`);
      }
      expect(view.resultsPane.root.closest('.forge-train-results')?.querySelector(':scope > .ui-box') ?? null, null,
        'A box is a direct child of the Results panel');
      const viewerRoots = () => ['Predicted Species vs Actual', 'Features'].map((t) => panelViewer(view, t).root);
      const roots = viewerRoots();
      const [frame] = view.resultsPane.frames;
      setHyperparameter(view, 'maxDepth', 1);
      await awaitTrained(view);
      expect(gridWith(view.resultsPane.root, 'Metric').root === grid, true, 'A retrain built a new metrics grid');
      const after = rowsOf(view);
      const retrained = rowTitles(view);
      for (let i = 0; i < rows.length; i++) {
        const isSame = JSON.stringify(titles[i]) === JSON.stringify(retrained[i]);
        expect(after[i] === rows[i], isSame, `Row ${i + 1}: ${JSON.stringify(titles[i])}`);
      }
      expect(after[0] === rows[0], true, 'Row 1 was rebuilt');
      expect(viewerRoots().every((root, i) => root === roots[i]), true, 'A retrain rebuilt the viewers');
      expect(view.resultsPane.frames[0] === frame && frame.columns.length > 0, true, 'A retrain replaced the frame');
      const training = view.lastTraining;
      if (training === undefined)
        throw new Error('No training');
      expectArray(frame.getCol('Predicted Species').toList().map(String),
        trainPredictionOf(training.result).toList().map(String));

      const shown = [...view.resultsPane.frames];
      view.featuresInput.value = view.featuresInput.value.slice(0, 3);
      await awaitTrained(view);
      await awaitCheck(() => view.resultsPane.frames[0] !== frame, 'A Features change kept the frame', WAIT_MS);
      expect(viewerRoots().every((root, i) => root !== roots[i]), true, 'A Features change kept the viewers');
      expectArray(shown.map((f) => f.columns.length), shown.map(() => 0));

      const current = [...view.resultsPane.frames];
      view.resultsPane.showHint('');
      expect(rowsOf(view).length, 0, 'The rows stay with a hint');
      expect(view.resultsPane.frames.length, 0, 'Frames stay with a hint');
      expectArray(current.map((f) => f.columns.length), current.map(() => 0));
      expect(resultsText(view).includes('Performance'), true, 'No Performance header with a hint');
      setHyperparameter(view, 'maxDepth', 6);
      await awaitTrained(view);
      expect(gridWith(view.resultsPane.root, 'Metric').root === grid, true, 'The hint dropped the metrics grid');
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: TIMEOUT});

  test('switching to another view and back keeps the rows and their sizes', async () => {
    const iris = await openIris();
    const view = await openTrainView(iris);
    let other: DG.ViewBase | undefined;
    try {
      await awaitCheck(() => rowsOf(view).every((row) => slotsOf(row).every((s) => s.offsetWidth > 0)),
        'The rows are not laid out', WAIT_MS);
      const heights = rowsOf(view).map((row) => row.offsetHeight);
      other = grok.shell.addView(DG.View.create());
      await delay(1000);
      grok.shell.v = view;
      await delay(1000);
      const rows = rowsOf(view);
      expectArray(rows.map((row) => row.offsetHeight), heights);
      expect(rows.every((row) => slotsOf(row).length > 0 && slotsOf(row).every((s) => s.offsetWidth > 0 &&
        s.offsetHeight > 0)), true, 'A panel collapsed');
    } finally {
      other?.close();
      await closeTrainView(view, [iris]);
    }
  }, {timeout: TIMEOUT});

  test('unbalanced iris: the class shares, Class balance\'s description, Wrong predictions under its title with row ' +
    'numbers, its frame given back when it goes', async () => {
    const unbalanced = await unbalancedIris();
    const view = await openTrainView(unbalanced);
    let shown: DG.DataFrame[] = [];
    try {
      const balance = panelViewer(view, 'Class balance');
      expect(balance.props.description, 'The classes are unbalanced.');
      expect(balance.props.descriptionVisibilityMode, 'Always');
      const line = insightItem(view, 'The classes are unbalanced.');
      expect(line.textContent, 'The classes are unbalanced.');
      // 41.67% each for setosa and versicolor take the two tenths left over: the shares sum to 100%. The pointer is the
      // tooltip's last line, in grey.
      expect(await enterTooltip(line), 'setosa: 41.7%versicolor: 41.7%virginica: 16.6%See Class balance.');
      const rect = line.getBoundingClientRect();
      line.dispatchEvent(new MouseEvent('mouseenter', {clientX: rect.left + 2, clientY: rect.top + 2}));
      await delay(600);
      const pointer = ui.tooltip.root.querySelector('.forge-tooltip-minor');
      ui.tooltip.hide();
      expect(pointer?.textContent, 'See Class balance.');

      setHyperparameter(view, 'iterations', 1);
      setHyperparameter(view, 'maxDepth', 1);
      await awaitTrained(view);
      await awaitCheck(() => rowTitles(view).flat().includes('Wrong predictions'), 'No Wrong predictions', WAIT_MS);
      expectRows(rowTitles(view).at(-1)!, ['Wrong predictions']);
      const wrong = gridWith(view.resultsPane.root, 'Predicted Species');
      expect(wrong.props.showRowHeader, true, 'Wrong predictions shows no row numbers');
      const title = panelOf(view, 'Wrong predictions').firstElementChild?.getBoundingClientRect();
      const canvas = wrong.canvas.getBoundingClientRect();
      expect(title !== undefined && canvas.top >= title.bottom - 0.5 && canvas.height > 0, true,
        `The grid from ${canvas.top}, the title to ${title?.bottom}`);
      expectUnderTitles(view);

      shown = [...view.resultsPane.frames];
      expect(shown.length, 2, 'Wrong predictions has no frame of its own');
      const mismatches = shown[1];
      setHyperparameter(view, 'iterations', 100);
      setHyperparameter(view, 'maxDepth', 6);
      await awaitTrained(view);
      await awaitCheck(() => !rowTitles(view).flat().includes('Wrong predictions'), 'Wrong predictions stays',
        WAIT_MS);
      expect(view.resultsPane.frames[0] === shown[0], true, 'The retrain replaced the viewers\' frame');
      expect(mismatches.columns.length, 0, 'The frame of Wrong predictions was not given back');
    } finally {
      await closeTrainView(view, [unbalanced]);
    }
    expectReleased(shown);
  }, {timeout: TIMEOUT});

  test('a viewer or the Wrong predictions grid that fails to build gives its frame back', expectFailedBuildReleased,
    {timeout: TIMEOUT});

  test('a retrain of a regression redraws the regression line', async () => {
    const iris = await openIris();
    const view = await openTrainView(iris);
    try {
      view.targetInput.value = iris.getCol('Petal.Width');
      view.featuresInput.value = columnsOf(iris, ['Sepal.Length', 'Sepal.Width', 'Petal.Length']);
      await awaitTrained(view);
      view.methodInput.value = 'XGBoost';
      await awaitCheck(() => view.lastTraining?.engine.name === 'XGBoost' && !view.isTraining, 'No XGBoost regression',
        TIMEOUT);
      const plot = panelViewer(view, 'Predicted Petal.Width vs Actual');
      if (!(plot instanceof DG.ScatterPlotViewer))
        throw new Error('Predicted Petal.Width vs Actual is not a scatter plot');
      const [frame] = view.resultsPane.frames;
      // The plot recomputes its line when the state of its inputs changes, the filter's version among them
      // (regression_line.dart getRegressionState), not the columns' values: the line cannot be read, its trigger can.
      let filterChanges = 0;
      const sub = frame.onFilterChanged.subscribe(() => filterChanges++);
      try {
        setHyperparameter(view, 'iterations', 1);
        await awaitTrained(view);
        await delay(500);
      } finally {
        sub.unsubscribe();
      }
      expect(view.resultsPane.frames[0] === frame && plot.root.isConnected, true, 'The retrain rebuilt the plot');
      expect(filterChanges > 0, true, 'The values changed in place without a filter change');
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: TIMEOUT});

  test('ROC Curve and PR Curve: lines over their own frames of Crux\'s curves, the AUC and AP of Performance; a ' +
    'threshold change moves only their marks, a retrain gives new frames', expectRocCurve, {timeout: 120000});

  test('ROC Curve and PR Curve from the class probabilities of a two-class classifier, Per-class for more classes; ' +
    'new frames on a retrain', expectProbabilityRocCurve, {timeout: 180000});

  test('ROC Curve stays when Predict probability goes on and off; the mark follows it, the grid is kept',
    expectProbabilityToggle, {timeout: 180000});

  test('Run Evaluation of a two-class model shows ROC Curve from the class probabilities', expectEvaluationRocCurve,
    {timeout: TIMEOUT});

  test('MetricsTable: the Evaluation column and bullet, a finding in the failure color\'s tooltip', async () => {
    const evaluated = new MetricsTable({metrics: {train: {accuracy: 0.9, f1: 0.9}, validation: {accuracy: 0.85,
      f1: 0.8}}, rowCount: 150, skippedRows: 0, evaluation: {values: {accuracy: 0.9, f1: 0.9}, tableName: 't',
      rowCount: 150, skippedRows: 0}});
    expectArray(gridWith(evaluated.root, 'Metric').dataFrame.columns.names(),
      ['Metric', 'Train', 'Validation', 'Evaluation']);
    expect(Array.from(evaluated.root.querySelectorAll('ul > li'), (li) => li.textContent)
      .includes('Evaluation: t (150 rows)'), true, evaluated.root.textContent ?? '');
    const low = new MetricsTable({metrics: {train: {mse: 1, rmse: 1, mae: 1, r2: 0.1}, validation: {mse: 1, rmse: 1,
      mae: 1, r2: 0.12}}, rowCount: 150, skippedRows: 0});
    const grid = gridWith(low.root, 'Metric');
    expect(gridTooltip(grid, grid.cell('Train', metricRow(grid, 'R2'))), 'R2 is below 0.3.');
    const plain = gridWith(evaluated.root, 'Metric');
    expect(gridTooltip(plain, plain.cell('Train', metricRow(plain, 'Accuracy'))), '', 'A plain value has a tooltip');

    const evaluationOnly = gridWith(new MetricsTable({metrics: {train: {accuracy: 0.9}, validation: {accuracy: 0.85}},
      rowCount: 150, skippedRows: 0, evaluation: {values: {accuracy: 0.9, auc: 0.8}, tableName: 't', rowCount: 150,
        skippedRows: 0}}).root, 'Metric');
    const metrics = evaluationOnly.dataFrame;
    expectArray(metrics.getCol('Metric').toList(), ['Accuracy', 'AUC-ROC']);
    const aucRow = metricRow(evaluationOnly, 'AUC-ROC');
    expect(metrics.getCol('Train').isNone(aucRow) && metrics.getCol('Validation').isNone(aucRow), true,
      'AUC-ROC has Train or Validation values');
    expectFloat(metrics.get('Evaluation', aucRow), 0.8, 1e-6);
  });

  test('EvaluationView: the Evaluation column equals Train on the training data; no column, no application',
    async () => {
      const iris = await openIris();
      const {id} = await saveIrisModel(iris);
      let view: EvaluationView | undefined;
      try {
        const model = await forgeDb.models.get(id);
        const loaded = await loadModel(id);
        const columnCount = iris.columns.length;
        const [opened, frames] = await framesSharing(iris.columns.toList(), () => EvaluationView.open(model, loaded,
          iris, iris.getCol('Species'), exactMapping(loaded.features, iris)));
        view = opened;
        expect(view.name, `Evaluation: ${model.name}`);
        expect(resultsState(view), 'false/null');
        const metrics = gridWith(view.resultsPane.root, 'Metric').dataFrame;
        expectArray(metrics.columns.names(), ['Metric', 'Train', 'Validation', 'Evaluation']);
        for (let i = 0; i < metrics.rowCount; i++)
          expectFloat(metrics.get('Evaluation', i), metrics.get('Train', i), 1e-6);
        const text = view.resultsPane.root.textContent ?? '';
        expect(text.includes(`Evaluation: ${iris.name} (150 rows)`), true, text);
        const titles = rowTitles(view);
        expectRows(titles[0], ['Predicted Species vs Actual', 'Confusions', 'Per-class']);
        expect(titles[1]?.includes('Performance'), true, JSON.stringify(titles));
        expect(titles[1][0], 'Insights & Tips', 'The evaluation hides Insights & Tips');
        // This model gets the overfitting tip, so the metric line stays out.
        expectArray(insightLines(view), ['The model may be overfitting.']);
        expect(iris.columns.length, columnCount, 'A column was added');
        expect(await forgeDb.applications.query().where('model_id', '=', id).count(), 0);
        view.close();
        view = undefined;
        expectReleased(frames);
      } finally {
        view?.close();
        await deleteModel(id);
      }
    }, {timeout: TIMEOUT});

  test('Run Evaluation opens a reference model\'s data without a dialog and outside the workspace', async () => {
    const iris = await openIrisFromFile();
    const ref = datasetRefOf(iris);
    let id: string | undefined;
    let view: DG.ViewBase | undefined;
    try {
      if (ref === null)
        throw new Error('iris has no origin');
      id = (await saveTestModel(columnsOf(iris, MEASUREMENTS), iris.getCol('Species'), iris.name,
        {storage_mode: 'reference', dataset_ref: ref})).id;
      grok.shell.closeTable(iris);
      const tableCount = grok.shell.tables.length;
      const dialogCount = DG.Dialog.getOpenDialogs().length;
      const model = await forgeDb.models.get(id);
      await runEvaluation(model);
      view = Array.from(grok.shell.views).find((v) => v.name === `Evaluation: ${model.name}`);
      expect(view instanceof EvaluationView, true, 'No evaluation view');
      expect(DG.Dialog.getOpenDialogs().length, dialogCount, 'A dialog opened');
      expect(grok.shell.tables.length, tableCount, 'The data joined the workspace');
    } finally {
      view?.close();
      if (isOpen(iris))
        grok.shell.closeTable(iris);
      if (id !== undefined)
        await deleteModel(id);
    }
  }, {timeout: TIMEOUT});

  test('Run Evaluation of a model without data asks for an open table, and warns without one', async () => {
    const iris = await openIris();
    const {id} = await saveIrisModel(iris);
    const unfit = await saveSkipUniqueModel();
    const tableView = grok.shell.addTableView(iris);
    let dialog: DG.Dialog | undefined;
    try {
      const shown = new Set(DG.Dialog.getOpenDialogs().map((d) => d.root));
      await runEvaluation(await forgeDb.models.get(id));
      dialog = DG.Dialog.getOpenDialogs().find((d) => !shown.has(d.root));
      expect(dialog?.title, 'Evaluate model');
      expect(dialog?.input('Table').value?.name, iris.name);
      expect(dialog?.input('Target').value?.name, 'Species');
      dialog?.close();
      await runEvaluation(await forgeDb.models.get(unfit.id));
      const warning = await readUntil(() => balloonWith('No open table has the columns the model needs. Open its ' +
        'training data or a similar table.'), (b) => b !== null);
      expect(warning !== null, true, 'No warning without a fitting table');
    } finally {
      dialog?.close();
      tableView.close();
      if (isOpen(iris))
        grok.shell.closeTable(iris);
      await deleteModel(id);
      await deleteModel(unfit.id);
    }
  }, {timeout: TIMEOUT});

  test('Performance has Run Evaluation for a saved model only', async () => {
    const iris = await openIris();
    const {id} = await saveIrisModel(iris);
    const manual = await insertModelRow({name: `forge-test-model-${Date.now()}-manual`});
    const accordions: DG.Accordion[] = [];
    const buttonsOf = async (modelId: string) => {
      const model = await forgeDb.models.get(modelId);
      const accordion = modelAccordion(ForgeModelHandler.rowOf(model), model);
      accordions.push(accordion);
      document.body.append(accordion.root);
      const pane = accordion.getPane('Performance');
      pane.expanded = true;
      return Array.from(pane.root.querySelectorAll('button'));
    };
    const texts = (buttons: HTMLElement[]) => buttons.map((b) => (b.textContent ?? '').trim());
    try {
      const buttons = await buttonsOf(id);
      expectArray(texts(buttons), ['Run Evaluation']);
      expect(await enterTooltip(buttons[0]),
        'Apply the model to its training data, or to an open table, and show how well it predicts.');
      expectArray(texts(await buttonsOf(manual)), []);
    } finally {
      for (const accordion of accordions)
        accordion.root.remove();
      await deleteModel(id);
      await forgeDb.models.delete(manual);
    }
  }, {timeout: TIMEOUT});
});

/** The lines of Insights & Tips of [view] with the metric's value and percentage masked. */
const metricLines = (view: TrainView) =>
  insightLines(view).map((l) => l.replace(/ -?\d\.\d{3} /, ' <v> ').replace(/\d+%/, '<p>%'));
const BALANCED_ACCURACY_LINE =
  'Balanced accuracy <v> on the validation rows: on average <p>% of each class is predicted correctly.';

// One category, split for its length.
category('UI: Results', () => {
  test('a metric value\'s tooltip: findings only; a tiny value in exponent form', async () => {
    expectArray([0.000032, -0.00004, 0, 0.0004, 0.0012].map(tinyValueText), ['3.2e-05', '-4.0e-05', null, '4.0e-04',
      null]);
    const finding: MetricFinding = {id: 'r2', kind: 'threshold', splits: ['train'], text: 'R2 is below 0.3.',
      detail: ''};
    expect(metricTooltip(0.5, []) === null, true, 'A plain value has a tooltip');
    expect(metricTooltip(0, [finding]), 'R2 is below 0.3.');
    expect(metricTooltip(0.000032, []), '3.2e-05.');
    expect(metricTooltip(undefined, [finding]) === null, true, 'No value has a tooltip');
    const grid = gridWith(new MetricsTable({metrics: {train: {r2: 0.9}, validation: {r2: 0.8}}, rowCount: 150,
      skippedRows: 0}).root, 'Metric');
    expect(headerTooltip(grid, 'Train').startsWith('Value on the training rows'), true, headerTooltip(grid, 'Train'));
    expect(headerTooltip(grid, 'Validation').startsWith('Value on rows the model did not see'), true,
      headerTooltip(grid, 'Validation'));
  });

  test('the layout stays while the threshold empties and fills Insights & Tips; one line reads the metric when ' +
    'nothing else does', async () => {
    const demog = await openDemog();
    const view = await TrainView.create(demog);
    grok.shell.addView(view);
    const {predictProbability, threshold, setThresholdButton} = view.outputInputs;
    const options = () => view.lastTraining?.result.options;
    let iris: DG.DataFrame | undefined;
    let irisView: TrainView | undefined;
    try {
      view.targetInput.value = demog.getCol('SEX');
      view.featuresInput.value = columnsOf(demog, ['AGE', 'HEIGHT', 'WEIGHT']);
      predictProbability.value = true;
      await awaitTrained(view);
      await awaitCheck(() => options()?.postprocessingInfo.includes('class-probabilities') === true && !view.isTraining,
        'Predict probability does not cut the result', WAIT_MS);
      const titles = rowTitles(view);
      const rows = rowsOf(view);
      expect(titles[1][0], 'Insights & Tips');
      expect(titles[2].includes('Class balance'), true, JSON.stringify(titles));

      threshold.value = 0.03;
      await awaitCheck(() => options()?.binaryClassificationThreshold === 0.03, 'The threshold does not re-cut',
        WAIT_MS);
      expectRows(rowTitles(view), titles);
      expect(rowsOf(view).every((r, i) => r === rows[i]), true, 'A re-cut rebuilt a row');
      expect(insightLines(view).some((l) => l.startsWith('The model finds ')), true, insightLines(view).join(' | '));
      setThresholdButton.click();
      await awaitCheck(() => options()?.thresholdMetric === 'f1', 'Set does not optimize for F1', WAIT_MS);
      expectRows(rowTitles(view), titles);
      expect(isShown(panelOf(view, 'Insights & Tips')), true, 'Insights & Tips is hidden');
      expectArray(metricLines(view), [BALANCED_ACCURACY_LINE]);

      iris = await twoSpeciesIris();
      irisView = await openTrainView(iris);
      // XGBoost gets the overfitting tip on these rows; LDA gets none.
      irisView.methodInput.value = 'LDA';
      await awaitCheck(() => irisView?.lastTraining?.engine.name === 'LDA' && !irisView.isTraining, 'No LDA training',
        TIMEOUT);
      expectArray(metricLines(irisView), [BALANCED_ACCURACY_LINE]);
      const tooltip = await enterTooltip(insightItem(irisView, 'Balanced accuracy'));
      expect(tooltip.includes('See Performance.'), true, tooltip);
      expectRows(rowTitles(irisView)[1], ['Insights & Tips', 'Confusions', 'Performance']);
    } finally {
      await closeTrainView(view, [demog]);
      if (irisView !== undefined && iris !== undefined)
        await closeTrainView(irisView, [iris]);
    }
  }, {timeout: TIMEOUT});
});

category('UI: Save model', () => {
  test('Save model offers Data storage: None and Copy without an origin', async () => {
    let saved: ModelStorageMode | undefined;
    const dialog = saveModelDialog('Iris model', async (info, choice) => {
      saved = choice.mode;
    }, {ref: null, rowCount: 150}).show();
    try {
      const storage = dialog.input('Data storage');
      expectArray(optionsOfRadio(storage), ['None', 'Copy']);
      expect(storage.value, 'None');
      const text = () => dialog.root.textContent ?? '';
      expect(text().includes('Only a summary of the data is saved.'), true, text());
      storage.value = 'Copy';
      await awaitCheck(() => text().includes('The training columns (150 rows) will be uploaded to the server.'),
        'No Copy line', WAIT_MS);
      dialog.getButton('OK').click();
      await awaitCheck(() => saved !== undefined, 'OK did not save', WAIT_MS);
      expect(saved, 'copy');
    } finally {
      dialog.close();
    }
  });

  test('Save model presets Reference for a table opened from a file', async () => {
    const iris = await openIrisFromFile();
    const ref = datasetRefOf(iris);
    grok.shell.closeTable(iris);
    const dialog = saveModelDialog('Iris model', async () => {}, {ref, rowCount: 150}).show();
    try {
      const storage = dialog.input('Data storage');
      expectArray(optionsOfRadio(storage), ['Reference', 'None', 'Copy']);
      expect(storage.value, 'Reference');
      const text = dialog.root.textContent ?? '';
      expect(text.includes(`A link to ${IRIS} is saved; the data stays where it is.`), true, text);
      expect(text.includes('training rows') || text.includes('stores them on the server'), false, text);
    } finally {
      dialog.close();
    }
  });

  test('Save in copy and reference modes; the balloon links to Apply and the catalog', async () => {
    const iris = await openIrisFromFile();
    iris.name = `forge-test-iris-${Date.now()}`;
    const view = await openTrainView(iris);
    const stamp = Date.now();
    const models = () => forgeDb.models.query().where('dataset_name', '=', iris.name);
    let catalog: DG.ViewBase | undefined;
    let apply: DG.Dialog | undefined;
    try {
      const save = saveDialogOf(view);
      expect(save?.input('Data storage').value, 'Reference');
      save?.close();
      const refName = `forge-test-model-${stamp}-ref`;
      const origin = datasetRefOf(iris);
      if (origin === null)
        throw new Error('iris has no origin');
      await view.saveModelAs({name: refName, description: '', tags: []}, {mode: 'reference', ref: origin});
      const balloon = await readUntil(() => balloonWith(`Model "${refName}" saved.`), (b) => b !== null);
      const links = Array.from(balloon?.querySelectorAll('a') ?? [], (a) => a.textContent);
      expectArray(links, ['Apply...', 'Show in the catalog']);
      balloon?.querySelectorAll('a')[0].click();
      await awaitCheck(() => DG.Dialog.getOpenDialogs().some((d) => d.title === 'Apply predictive model'),
        'Apply... opens no dialog', WAIT_MS);
      apply = DG.Dialog.getOpenDialogs().find((d) => d.title === 'Apply predictive model');
      expect(apply?.input('Model').value, refName);
      expect(apply?.input('Table').value?.name, iris.name);
      apply?.close();
      balloon?.querySelectorAll('a')[1].click();
      await awaitCheck(() => grok.shell.v instanceof ForgeApp, 'Show in the catalog opens no catalog', WAIT_MS);
      catalog = grok.shell.v;
      const [ref] = await models();
      expect(ref.storage_mode, 'reference');
      expect(storedDatasetRef(ref.dataset_ref)?.path, IRIS);

      setHyperparameter(view, 'iterations', 25);
      await awaitTrained(view);
      const copyName = `forge-test-model-${stamp}-copy`;
      await view.saveModelAs({name: copyName, description: '', tags: []}, {mode: 'copy'});
      const copy = (await models()).find((m) => m.storage_mode === 'copy');
      const tableId = copy?.dataset_table_id;
      if (copy === undefined || !tableId)
        throw new Error('No copy-mode model, or it has no table');
      const info: DG.TableInfo | undefined = await grok.dapi.tables.find(tableId);
      expect(info?.friendlyName, `${copyName} (training data)`);
      const accordion = modelAccordion(ForgeModelHandler.rowOf(copy), copy);
      await expectPaneText(accordion, 'Details', ['Data storage', 'Copy', 'Data copy', '(training data)']);
      const refAccordion = modelAccordion(ForgeModelHandler.rowOf(ref), ref);
      await expectPaneText(refAccordion, 'Details', ['Data storage', 'Reference', 'Data source', IRIS]);
      await deleteModel(copy.id);
      const left: DG.TableInfo | undefined = await grok.dapi.tables.find(tableId);
      expect(left === undefined || left === null, true, 'The copy outlived its model');
    } finally {
      apply?.close();
      catalog?.close();
      for (const model of await models())
        await deleteModel(model.id);
      await closeTrainView(view, [iris]);
      grok.shell.closeTable(iris);
    }
  }, {timeout: TIMEOUT});
});

category('UI: Apply dialog', () => {
  test('Apply dialog opens prefilled', async () => {
    const iris = await openIris();
    const {id, name} = await saveIrisModel(iris);
    let opened: DG.Dialog | undefined;
    try {
      const dialog = (await applyModelDialog({table: iris, modelId: id})).show();
      opened = dialog;
      const ok = dialog.getButton('OK');
      const text = dialog.root.textContent ?? '';
      for (const caption of ['Table', 'Model', 'Batch size', ...MEASUREMENTS])
        expect(text.includes(caption), true, `${caption} is missing`);
      expect(isShown(dialog.input('Missing values').root), false, 'Missing values is shown without gaps');
      expect(dialog.input('Model').value, name);
      expect(dialog.input('Sepal.Length').value?.name, 'Sepal.Length');
      expect(dialog.input('Batch size').value, 10000);
      await awaitCheck(() => !isDisabled(ok), 'OK stays disabled', WAIT_MS);
      const columns = groupOf(dialog.root, 'Columns');
      expect(columns.summary.textContent, '4 of 4 matched');
      expect(columns.summary.classList.contains('forge-group-invalid'), false, 'The summary is red');
      expect(columns.summary.getAttribute('aria-invalid'), 'false');
      expect(isShown(columns.body), false, 'Columns is expanded with every feature matched');
      expect(columns.body.contains(dialog.input('Sepal.Length').root), true, 'The rows are not under Columns');
      const moreOptions = groupOf(dialog.root, 'More options');
      expect(isShown(moreOptions.body), false, 'More options is expanded');
      expect(moreOptions.body.contains(dialog.input('Batch size').root), true, 'Batch size is not under More options');

      dialog.input('Sepal.Width').value = iris.getCol('Sepal.Length');
      await awaitCheck(() => isDisabled(ok) && (dialog.input('Sepal.Width').validity ?? '').includes('also used'),
        'A column used twice is not marked', WAIT_MS);
      expect(columns.summary.textContent, '3 of 4 matched');
      expect(columns.summary.classList.contains('forge-group-invalid'), true, 'The summary is not red');
      expect(columns.summary.getAttribute('aria-invalid'), 'true');
      expect(isShown(columns.body), true, 'Columns stays collapsed with a problem');
      expect(columns.summary.parentElement?.getAttribute('aria-expanded'), 'true');
      expect(isShown(moreOptions.body), false, 'More options opened');
      expectAlignedCaptions(dialog.root, ['Columns', 'More options']);
    } finally {
      opened?.close();
      await deleteModel(id);
    }
  }, {timeout: TIMEOUT});

  test('Apply dialog does not ask for the columns Skip unique categories left out', async () => {
    const {id} = await saveSkipUniqueModel();
    const table = DG.DataFrame.fromColumns([DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'x', [1, 2, 3, 4])]);
    let opened: DG.Dialog | undefined;
    try {
      const dialog = (await applyModelDialog({table, modelId: id})).show();
      opened = dialog;
      await awaitCheck(() => !isDisabled(dialog.getButton('OK')), 'OK stays disabled without subject', WAIT_MS);
      expect(groupOf(dialog.root, 'Columns').summary.textContent, '1 of 1 matched');
      expect((dialog.root.textContent ?? '').includes('subject'), false, 'The dialog asks for subject');
    } finally {
      opened?.close();
      await deleteModel(id);
    }
  }, {timeout: TIMEOUT});

  test('Apply dialog closes on OK and applies in the task bar', async () => {
    const iris = await openIris();
    const {id} = await saveIrisModel(iris);
    let opened: DG.Dialog | undefined;
    let applying = false;
    try {
      const dialog = (await applyModelDialog({table: iris, modelId: id})).show();
      opened = dialog;
      const ok = dialog.getButton('OK');
      await awaitCheck(() => !isDisabled(ok), 'OK stays disabled', WAIT_MS);
      applying = true;
      ok.click();
      await expectApplied(dialog, iris, id);
    } finally {
      opened?.close();
      if (applying)
        await applicationSettled(id);
      await deleteModel(id);
    }
  }, {timeout: TIMEOUT});

  test('Apply dialog of a model saved with Predict probability adds its probability column', async () => {
    const iris = await twoSpeciesIris();
    const {id} = await saveTestModel(columnsOf(iris, MEASUREMENTS), iris.getCol('Species'), iris.name, {},
      PREDICT_PROBABILITY);
    let opened: DG.Dialog | undefined;
    let applying = false;
    try {
      const dialog = (await applyModelDialog({table: iris, modelId: id})).show();
      opened = dialog;
      const ok = dialog.getButton('OK');
      await awaitCheck(() => !isDisabled(ok), 'OK stays disabled', WAIT_MS);
      applying = true;
      ok.click();
      await expectApplied(dialog, iris, id);
      await awaitCheck(() => iris.col('P(Species = versicolor)') !== null, 'No probability column', WAIT_MS);
      await awaitCheck(() => balloonWith('Added the columns "Species (predicted)", "P(Species = versicolor)"') !== null,
        'No balloon naming both columns', WAIT_MS);
    } finally {
      opened?.close();
      if (applying)
        await applicationSettled(id);
      await deleteModel(id);
    }
  }, {timeout: TIMEOUT});

  test('Apply dialog opens the Columns rows scrolled to the top', async () => {
    const rows = 40;
    const features = Array.from({length: 30}, (_, k) => DG.Column.fromList(DG.COLUMN_TYPE.FLOAT,
      `feature ${String(k + 1).padStart(2, '0')}`, valuesOf(rows, (i) => (i * (k + 3)) % 17)));
    const y = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'y', valuesOf(rows, (i) => i % 5));
    const {id} = await saveTestModel(features, y, `forge-test-wide-${Date.now()}`);
    let opened: DG.Dialog | undefined;
    try {
      const dialog = (await applyModelDialog({table: await openIris(), modelId: id})).show();
      opened = dialog;
      const block = dialog.root.querySelector('.forge-apply-rows');
      if (!(block instanceof HTMLElement))
        throw new Error('No rows block');
      // The dialog moves the focus to an input on the next tick after show().
      await delay(500);
      expect(isShown(groupOf(dialog.root, 'Columns').body), true, 'Columns is collapsed with 30 empty rows');
      expect(getComputedStyle(block).overflowY, 'auto');
      expect(block.scrollHeight > block.clientHeight, true, 'The rows do not overflow; the check would prove nothing');
      expect(block.scrollTop, 0);
    } finally {
      opened?.close();
      await deleteModel(id);
    }
  }, {timeout: TIMEOUT});

  test('Apply dialog shows Missing values and the impute inputs for a Skip rows XGBoost model', async () => {
    const iris = await openIris();
    const {id} = await saveIrisModel(iris);
    iris.getCol('Sepal.Width').set(3, null);
    let opened: DG.Dialog | undefined;
    try {
      const dialog = (await applyModelDialog({table: iris, modelId: id})).show();
      opened = dialog;
      const missingValues = dialog.input('Missing values');
      expect(isShown(missingValues.root), true, 'Missing values is hidden with a gap');
      expect(missingValues.inputType, 'Radio');
      expect(missingValues.value, 'Skip rows');
      expectOptionsInRow(missingValues);
      expect(isShown(dialog.input('Neighbors').root), false, 'Neighbors is shown with Skip rows');
      await awaitCheck(() => !isDisabled(dialog.getButton('OK')), 'OK is disabled with Skip rows', WAIT_MS);

      missingValues.value = 'Impute';
      await awaitCheck(() => isShown(dialog.input('Neighbors').root) && isShown(dialog.input('Distance').root),
        'The impute inputs do not appear', WAIT_MS);
      expect(dialog.input('Neighbors').value, 4);
      expect(dialog.input('Distance').value, 'Euclidean');
      missingValues.value = 'Skip rows';
      await awaitCheck(() => !isShown(dialog.input('Neighbors').root), 'The impute inputs stay', WAIT_MS);
    } finally {
      opened?.close();
      await deleteModel(id);
    }
  }, {timeout: TIMEOUT});

  test('Apply dialog hides Missing values for a model trained with the gaps kept', async () => {
    const iris = await openIris();
    const {id} = await saveKeptGapsModel(iris);
    iris.getCol('Sepal.Width').set(3, null);
    let opened: DG.Dialog | undefined;
    try {
      const dialog = (await applyModelDialog({table: iris, modelId: id})).show();
      opened = dialog;
      expect(isShown(dialog.input('Missing values').root), false, 'Missing values is shown for XGBoost');
      await awaitCheck(() => !isDisabled(dialog.getButton('OK')), 'OK is disabled', WAIT_MS);
    } finally {
      opened?.close();
      await deleteModel(id);
    }
  }, {timeout: TIMEOUT});

  test('Apply dialog: rows without a close column stay empty, a refused model shows why', async () => {
    const iris = await openIris();
    const {id} = await saveIrisModel(iris);
    const other = DG.DataFrame.fromColumns(['AGE', 'HEIGHT', 'WEIGHT'].map((name) =>
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, name, [30, 170, 70])));
    const name = `forge-test-model-${Date.now()}-manual`;
    let manualId: string | undefined;
    let opened: DG.Dialog | undefined;
    try {
      manualId = (await forgeDb.models.insert({...XGBOOST_FIELDS, name, task: 'classification',
        target_name: 'Species', storage_mode: 'none'}))[0].id;
      const empty = (await applyModelDialog({table: other, modelId: id})).show();
      opened = empty;
      const ok = empty.getButton('OK');
      for (const feature of MEASUREMENTS)
        expect(empty.input(feature).value === null, true, `${feature} is prefilled`);
      await awaitCheck(() => isDisabled(ok), 'OK is enabled with empty rows', WAIT_MS);
      const columns = groupOf(empty.root, 'Columns');
      expect(columns.summary.textContent, '0 of 4 matched');
      expect(columns.summary.classList.contains('forge-group-invalid'), true, 'The summary is not red');
      expect(isShown(columns.body), true, 'Columns is collapsed with empty rows');

      empty.input('Model').value = name;
      await awaitCheck(() => (empty.root.textContent ?? '').includes(`The model '${name}' has no model file`),
        'The refusal is not shown', WAIT_MS);
      expect(isDisabled(ok), true, 'OK is enabled for a refused model');
    } finally {
      opened?.close();
      if (manualId !== undefined)
        await forgeDb.models.delete(manualId);
      await deleteModel(id);
    }
  }, {timeout: TIMEOUT});

  test('modelLabels tells apart models with the same name', async () => {
    const t0 = dayjs('2026-01-02T10:00:05');
    const rows = [
      {id: 'a', name: 'N', created_on: t0},
      {id: 'b', name: 'N', created_on: t0.add(1, 'minute')},
      {id: 'c', name: 'N', created_on: t0.add(25, 'second')},
      {id: 'd', name: 'N', created_on: t0.add(25, 'second')},
      {id: 'e', name: 'M', created_on: t0},
    ];
    const labels = [...modelLabels(rows)];
    const labelOf = (id: string) => labels.find(([, row]) => row.id === id)?.[0];
    expect(labels.length, rows.length);
    expect(labelOf('e'), 'M');
    expect(labelOf('b'), 'N (2026-01-02 10:01)');
    expect(labelOf('a'), 'N (2026-01-02 10:00:05)');
    expect(labelOf('c'), 'N (2026-01-02 10:00:30)');
    expect(labelOf('d'), 'N (2026-01-02 10:00:30) #2');
  });

  test('Apply dialog lists every model of a repeated name and opens the chosen one', async () => {
    const name = `forge-test-model-${Date.now()}-same`;
    const ids: string[] = [];
    let opened: DG.Dialog | undefined;
    try {
      for (let i = 0; i < 3; i++) {
        ids.push((await forgeDb.models.insert({...XGBOOST_FIELDS, name, task: 'classification',
          target_name: 'Species', storage_mode: 'none'}))[0].id);
      }
      const labels = [...modelLabels(await forgeDb.models.query().where('name', '=', name))];
      const dialog = (await applyModelDialog({table: await openIris(), modelId: ids[1]})).show();
      opened = dialog;
      const listed = optionsOf(dialog, 'Model').filter((label) => label.startsWith(name));
      expectArray(listed.sort(), labels.map(([label]) => label).sort());
      expect(listed.length, 3);
      expect(dialog.input('Model').value, labels.find(([, row]) => row.id === ids[1])?.[0]);
    } finally {
      opened?.close();
      for (const id of ids)
        await forgeDb.models.delete(id);
    }
  }, {timeout: TIMEOUT});

  test('Apply dialog re-orders the models and re-fills the rows on a Table change', async () => {
    const rows = 40;
    const alpha = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'alpha', valuesOf(rows, (i) => i % 7));
    const beta = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'beta', valuesOf(rows, (i) => i / 3));
    const gamma = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'gamma', valuesOf(rows, (i) => i % 7 + i / 3));
    const other = DG.DataFrame.fromColumns([alpha, beta, gamma]);
    other.name = `forge-test-other-${Date.now()}`;
    const tables = [grok.shell.addTable(await openIris()), grok.shell.addTable(other)];
    const iris = tables[0];
    const ids: string[] = [];
    let opened: DG.Dialog | undefined;
    try {
      const otherModel = await saveTestModel([alpha, beta], gamma, other.name);
      ids.push(otherModel.id);
      const irisModel = await saveIrisModel(iris);
      ids.push(irisModel.id);
      const dialog = (await applyModelDialog({table: iris, modelId: irisModel.id})).show();
      opened = dialog;
      const columns = groupOf(dialog.root, 'Columns');
      const isFirst = (first: string, second: string) =>
        optionsOf(dialog, 'Model').indexOf(first) < optionsOf(dialog, 'Model').indexOf(second);
      expect(isFirst(irisModel.name, otherModel.name), true, optionsOf(dialog, 'Model').join(', '));

      dialog.input('Table').value = other;
      await awaitCheck(() => columns.summary.textContent === '0 of 4 matched', 'The rows were not re-filled', WAIT_MS);
      expect(dialog.input('Model').value, irisModel.name);
      expect(isFirst(otherModel.name, irisModel.name), true, optionsOf(dialog, 'Model').join(', '));
      for (const feature of MEASUREMENTS)
        expect(dialog.input(feature).value === null, true, `${feature} is prefilled from the other table`);

      dialog.input('Table').value = iris;
      await awaitCheck(() => columns.summary.textContent === '4 of 4 matched', 'The iris rows were not re-filled',
        WAIT_MS);
      expect(isFirst(irisModel.name, otherModel.name), true, optionsOf(dialog, 'Model').join(', '));
      expect(dialog.input('Sepal.Length').value?.name, 'Sepal.Length');
    } finally {
      opened?.close();
      for (const id of ids)
        await deleteModel(id);
      for (const table of tables)
        grok.shell.closeTable(table);
    }
  }, {timeout: TIMEOUT});

  test('Apply dialog marks a mapped column removed while it is open, and Enter waits for a fix', async () => {
    const iris = await openIris();
    const {id} = await saveIrisModel(iris);
    let opened: DG.Dialog | undefined;
    let applying = false;
    try {
      const dialog = (await applyModelDialog({table: iris, modelId: id})).show();
      opened = dialog;
      const ok = dialog.getButton('OK');
      await awaitCheck(() => !isDisabled(ok), 'OK stays disabled', WAIT_MS);
      const sepalWidth = iris.getCol('Sepal.Width');
      iris.columns.remove('Sepal.Width');
      await awaitCheck(() => isDisabled(ok) &&
        (dialog.input('Sepal.Width').validity ?? '').includes('is no longer in the table'),
      'The removed column is not marked', WAIT_MS);
      const columns = groupOf(dialog.root, 'Columns');
      expect(columns.summary.textContent, '3 of 4 matched');
      expect(isShown(columns.body), true, 'Columns stays collapsed with a problem');
      pressEnter(dialog, 'Model');
      await delay(500);
      expect(dialog.root.isConnected, true, 'Enter closed the dialog with a problem');
      expect(iris.col('Species (predicted)') === null, true, 'Enter applied the model with a problem');

      iris.columns.add(sepalWidth);
      await awaitCheck(() => !isDisabled(ok), 'OK stays disabled with the column back', WAIT_MS);
      applying = true;
      pressEnter(dialog, 'Model');
      await expectApplied(dialog, iris, id);
    } finally {
      opened?.close();
      if (applying)
        await applicationSettled(id);
      await deleteModel(id);
    }
  }, {timeout: TIMEOUT});
});

let catalogModel: {id: string; iris: DG.DataFrame} | undefined;
const sharedModel = () => savedFixture(catalogModel);

/** Types [text] into the Tags box inside [root] and presses Enter, as a user adds a tag; waits for the chip. */
async function typeTag(root: HTMLElement, text: string): Promise<void> {
  const box = root.querySelector('input.d4-tags-selector-input');
  if (!(box instanceof HTMLInputElement))
    throw new Error('No Tags text box');
  typeInto(box, text);
  for (const type of ['keydown', 'keypress', 'keyup'])
    box.dispatchEvent(new KeyboardEvent(type, {key: 'Enter', keyCode: 13, bubbles: true}));
  const chips = () => Array.from(root.querySelectorAll('.d4-tags-selector-tags-container > *'), (c) => c.textContent);
  await awaitCheck(() => chips().includes(text), `No chip '${text}': ${chips()}`, WAIT_MS);
}

/** Opens the accordion pane [name] and waits until its text has every one of [texts]; fails with the text it has. */
async function expectPaneText(accordion: DG.Accordion, name: string, texts: string[]): Promise<void> {
  const pane = accordion.getPane(name);
  pane.expanded = true;
  const hasAll = (text: string) => texts.every((t) => text.includes(t));
  const text = await readUntil(() => pane.root.textContent ?? '', hasAll);
  expect(hasAll(text), true, `${name}: ${text}`);
}

category('UI: Catalog', () => {
  before(async () => {
    useTestGalleryStorage();
    const iris = await openIris();
    catalogModel = {id: (await saveIrisModel(iris)).id, iris};
  });

  after(async () => {
    clearTestGalleryStorage();
    if (catalogModel !== undefined)
      await deleteModel(catalogModel.id);
  });

  test('the autostart registered the model handler', async () => {
    expect(DG.ObjectHandler.list().some((h) => h.name === 'Forge model handler'), true);
    const row = ForgeModelHandler.rowOf(await forgeDb.models.get(sharedModel().id));
    expect(DG.ObjectHandler.forEntity(row)?.name, 'Forge model handler');
  });

  test('the catalog toolbar has Refresh alone, a search without a magnifier, Applicable to under it', async () => {
    const view = await ForgeApp.create();
    grok.shell.addView(view);
    try {
      const header = view.root.querySelector('.grok-gallery-search-bar .forge-gallery-actions');
      const icons = Array.from(header?.querySelectorAll<HTMLElement>('.grok-icon') ?? []);
      expect(icons.length, 1, 'The toolbar has other icons than Refresh');
      const searchBar = view.modelsGallery.search.root.closest('.d4-search-bar');
      expect(searchBar !== null, true, 'The search box is not in the platform\'s search bar');
      const magnifier = view.modelsGallery.search.root.querySelector<HTMLElement>('i.ui-input-icon-left');
      expect(magnifier === null || getComputedStyle(magnifier).display === 'none', true,
        'The search box shows a magnifier');
      expect(header?.parentElement?.nextElementSibling, view.applicableToInput.root,
        'Applicable to is not on the line under the toolbar');
      const lefts = () => [view.applicableToInput.captionLabel, header]
        .map((e) => e?.getBoundingClientRect().left ?? -100);
      await waitUntil(() => Math.abs(lefts()[0] - lefts()[1]) < 2,
        `Applicable to and the icons start at ${lefts()}`, WAIT_MS);
      const colorOf = (e: Element | null | undefined) => e instanceof Element ? getComputedStyle(e).color : 'none';
      expect(colorOf(icons[0]), colorOf(view.applicableToInput.root.querySelector('.ui-input-options > i')),
        'The Refresh icon is not the blue of the input\'s icons');
      expect(view.applicableToInput.inputType, DG.InputType.Table);
      expect(view.applicableToInput.value, null);
      expect(view.models.col('tags') !== null && view.models.col('features') !== null, true,
        view.models.columns.names().join(', '));
    } finally {
      view.close();
    }
  });

  test('a chosen row is the current object; Compare follows the selection', async () => {
    const {id} = sharedModel();
    const otherId = await insertModelRow({name: `forge-test-model-${Date.now()}-compare`});
    const view = await ForgeApp.create();
    grok.shell.addView(view);
    try {
      const ids = view.models.getCol('id').toList();
      view.models.currentRowIdx = ids.indexOf(id);
      await awaitCheck(() => shows(id),
        'The chosen model is not the current object', WAIT_MS);
      const current: unknown = grok.shell.o;
      expect(current instanceof DG.DomainRow && current.typeName === MODEL_TYPE, true, 'Not a model');

      view.models.selection.set(ids.indexOf(id), true);
      view.models.selection.set(ids.indexOf(otherId), true);
      await awaitCheck(() => comparedIds() === [id, otherId].sort().join(),
        `The context panel does not compare the two models: ${grok.shell.o}`, WAIT_MS);
      expect(DG.ObjectHandler.forEntity(grok.shell.o)?.name, 'Forge model comparison handler');
      await view.compareSelected();
      const compared = grok.shell.tv;
      try {
        expect(compared.name, 'Compare models');
        expect(compared.dataFrame.rowCount, 2);
        const types = Array.from(compared.viewers, (v) => v.type);
        // The viewer reports its registered name or its class name, from run to run.
        expect(types.some((t) => t === 'Forms' || t === 'FormsViewer'), hasFormsViewer(),
          `The Compare view's viewers: ${types}`);
      } finally {
        compared.close();
      }
      view.models.selection.set(ids.indexOf(otherId), false);
      await awaitCheck(() => shows(id),
        'One selected row does not show its model again', WAIT_MS);

      view.models.currentRowIdx = -1;
      view.models.selection.set(ids.indexOf(otherId), true);
      await awaitCheck(() => comparedIds() === [id, otherId].sort().join(), 'No comparison without a current row',
        WAIT_MS);
      view.models.selection.set(ids.indexOf(id), false);
      await awaitCheck(() => shows(otherId),
        'Without a current row, the one selected row does not show its model', WAIT_MS);
      view.models.selection.set(ids.indexOf(id), true);
      await awaitCheck(() => comparedIds() !== '', 'The two are not compared again', WAIT_MS);
      view.models.selection.setAll(false);
      await awaitCheck(() => (grok.shell.o ?? null) === null, 'With nothing chosen the panel keeps the comparison',
        WAIT_MS);
    } finally {
      view.close();
      await forgeDb.models.delete(otherId);
    }
  }, {timeout: TIMEOUT});

  test('a refresh keeps the current row and the selection', async () => {
    const {id} = sharedModel();
    const otherId = await insertModelRow({name: `forge-test-model-${Date.now()}-refresh`});
    const view = await ForgeApp.create();
    grok.shell.addView(view);
    try {
      const indexOf = (modelId: string) => view.models.getCol('id').toList().indexOf(modelId);
      const selectedIds = () => Array.from(view.models.selection.getSelectedIndexes(),
        (i) => view.models.get('id', i)).sort();
      const reload = async (description: string, load: number) => {
        const frame = view.models.dart;
        await updateModelInfo(otherId, undefined, {description});
        expect(catalogState(view.root), `${load - 1}/true`, 'The catalog is not busy before the debounced reload');
        await awaitCheck(() => view.models.dart !== frame, 'The catalog does not reload', WAIT_MS);
        expect(catalogState(view.root), `${load}/false`);
        await awaitCheck(() => view.models.currentRowIdx === indexOf(id), 'The current row is lost', WAIT_MS);
      };
      view.models.currentRowIdx = indexOf(id);
      await awaitCheck(() => shows(id),
        'The chosen model is not the current object', WAIT_MS);
      const shown: unknown = grok.shell.o;
      const shownDart: unknown = shown instanceof DG.DomainRow ? shown.dart : null;
      await reload('refreshed', 1);
      const now: unknown = grok.shell.o;
      expect(now instanceof DG.DomainRow && now.dart === shownDart, true,
        'The context panel got a new object for the same model');

      view.models.selection.set(indexOf(id), true);
      view.models.selection.set(indexOf(otherId), true);
      await awaitCheck(() => comparedIds() !== '', 'Two selected rows are not compared', WAIT_MS);
      await reload('refreshed again', 2);
      expectArray(selectedIds(), [id, otherId].sort());
      await awaitCheck(() => comparedIds() === [id, otherId].sort().join(), 'The comparison is lost', WAIT_MS);
      const raced = await refreshThenChange(view.root, () => view.models.dart);
      expect(`${raced.later} ${raced.isRebound}`, `${raced.first} false`, 'A Refresh in flight published first');
    } finally {
      view.close();
      await forgeDb.models.delete(otherId);
    }
  });

  test('a deleted model leaves the context panel', async () => {
    const {id} = sharedModel();
    const stamp = Date.now();
    const shownId = await insertModelRow({name: `forge-test-model-${stamp}-deleted`});
    const comparedId = await insertModelRow({name: `forge-test-model-${stamp}-compared`});
    const view = await ForgeApp.create();
    grok.shell.addView(view);
    const indexOf = (modelId: string) => view.models.getCol('id').toList().indexOf(modelId);
    // The test bundle's delete fires its own modelsChanged, which the follower _initForge started does not hear.
    const follow = ForgeModelHandler.followChanges();
    try {
      view.models.currentRowIdx = indexOf(shownId);
      await awaitCheck(() => shows(shownId), 'The model is not the current object', WAIT_MS);
      const menuDelete = CATALOG_ACTIONS.find((a) => a.name === 'Delete model');
      await menuDelete?.run(ForgeModelHandler.rowOf(await forgeDb.models.get(shownId)), null);
      const confirm = DG.Dialog.getOpenDialogs().find((d) => d.title === 'Delete model');
      if (confirm === undefined)
        throw new Error('No Delete model dialog');
      confirm.getButton('OK').click();
      await awaitCheck(() => indexOf(shownId) < 0, 'The deleted model stays in the catalog', WAIT_MS);
      await awaitCheck(() => !shows(shownId), 'The context panel keeps the deleted model', WAIT_MS);

      view.models.selection.set(indexOf(id), true);
      view.models.selection.set(indexOf(comparedId), true);
      await awaitCheck(() => comparedIds() === [id, comparedId].sort().join(), 'The two are not compared', WAIT_MS);
      await menuDelete?.run(ForgeModelHandler.rowOf(await forgeDb.models.get(comparedId)), null);
      DG.Dialog.getOpenDialogs().find((d) => d.title === 'Delete model')?.getButton('OK').click();
      await awaitCheck(() => shows(id), 'A comparison that lost a model does not show the model left', WAIT_MS);
    } finally {
      follow.unsubscribe();
      view.close();
      const left: {id: string}[] = await forgeDb.models.query().where('id', '=', [shownId, comparedId]).select('name')
        .top(2);
      for (const row of left)
        await forgeDb.models.delete(row.id);
    }
  }, {timeout: TIMEOUT});

  test('Applicable to keeps only the models a chosen table fits', async () => {
    const {id} = sharedModel();
    const tables = [grok.shell.addTable(await openIris()),
      grok.shell.addTable(await grok.data.files.openTable('System:DemoFiles/cars.csv'))];
    const view = await ForgeApp.create();
    grok.shell.addView(view);
    const isShown = () => view.models.filter.get(view.models.getCol('id').toList().indexOf(id));
    try {
      view.applicableToInput.value = tables[1];
      await awaitCheck(() => !isShown(), 'The iris model stays with cars', WAIT_MS);
      view.applicableToInput.value = tables[0];
      await awaitCheck(() => isShown(), 'The iris model is hidden with iris', WAIT_MS);
      view.applicableToInput.value = null;
      await awaitCheck(() => isShown() && view.models.filter.trueCount === view.models.rowCount,
        'An empty Applicable to does not show every model', WAIT_MS);
      view.applicableToInput.value = tables[1];
      await awaitCheck(() => !isShown(), 'The iris model stays with cars', WAIT_MS);
      grok.shell.closeTable(tables[1]);
      await awaitCheck(() => view.applicableToInput.value === null && isShown(), 'A closed table stays chosen',
        WAIT_MS);
    } finally {
      view.close();
      for (const table of tables.filter(isOpen))
        grok.shell.closeTable(table);
    }
  });

  test('Apply... opens on the preferred, the current or the first table the model fits', async () => {
    const {id} = sharedModel();
    const iris = await openIris();
    const cars = await grok.data.files.openTable('System:DemoFiles/cars.csv');
    const noFeatures = await insertModelRow({name: `forge-test-model-${Date.now()}-preset`});
    const tableOf = async (options: ApplyDialogOptions): Promise<DG.DataFrame> => {
      const dialog = (await applyModelDialog(options)).show();
      const table: DG.DataFrame | null = dialog.input('Table').value;
      dialog.close();
      if (table === null)
        throw new Error('No table');
      return table;
    };
    const views = [grok.shell.addTableView(iris), grok.shell.addTableView(cars)];
    try {
      expect((await tableOf({table: cars, modelId: id, preferredTable: iris})).name, iris.name);
      const fitting = await tableOf({table: cars, modelId: id, preferredTable: null});
      expect(fitting.name !== cars.name && MEASUREMENTS.every((m) => fitting.col(m) !== null), true,
        `With cars current the dialog opens on ${fitting.name}`);
      expect((await tableOf({table: cars, modelId: noFeatures})).name, cars.name);
      grok.shell.v = views[0];
      expect((await tableOf({table: cars, modelId: id, preferredTable: null})).name, iris.name);
      expect((await tableOf({table: cars})).name, cars.name);
    } finally {
      for (const view of views)
        view.close();
      await forgeDb.models.delete(noFeatures);
    }
  });

  test('the accordion title\'s name is bound: the platform\'s name attribute and the model\'s tooltip', async () => {
    const model = await forgeDb.models.get(sharedModel().id);
    const accordion = modelAccordion(ForgeModelHandler.rowOf(model), model);
    document.body.append(accordion.root);
    try {
      const name = titleName(accordion.root, model.name);
      if (name === undefined)
        throw new Error('The title has no name label');
      expect(name.getAttribute('name')?.startsWith('label-'), true, `name="${name.getAttribute('name')}"`);
      const tooltip = await enterTooltip(name);
      expect(tooltip.includes('Predicts Species with XGBoost from Sepal.Length, Sepal.Width'), true, tooltip);
    } finally {
      ui.tooltip.hide();
      accordion.root.remove();
    }
  });

  test('the model accordion has the five panes and reads the model', async () => {
    const {id, iris} = sharedModel();
    const model = await forgeDb.models.get(id);
    const row = ForgeModelHandler.rowOf(model);
    const properties = new ForgeModelHandler().renderProperties(row);
    await awaitCheck(() => [model.name, 'Details', 'Performance', 'Activity', 'Sharing', 'History']
      .every((name) => (properties.textContent ?? '').includes(name)), 'The context panel lacks its title or a pane',
    WAIT_MS);
    expect(properties.querySelector('.svg-model') !== null, true, 'The title has no model icon');

    // Attached at the context panel's width, so the grids are measured as shown.
    const host = ui.div([]);
    host.style.width = '320px';
    document.body.append(host);
    const measured: string[] = [];
    try {
      const accordion = modelAccordion(row, model);
      host.append(accordion.root);
      const title = accordion.root.firstElementChild;
      expect(title?.classList.contains('d4-accordion-title'), true, `The first element is ${title?.className}`);
      expect(title?.querySelector('.d4-star') !== null && title?.textContent?.includes(model.name), true,
        'The title has no star or no name');
      expectArray(accordion.panes.map((p) => p.name), ['Details', 'Performance', 'Activity', 'Sharing', 'History']);
      await expectPaneText(accordion, 'Details', ['Author', 'Created', 'Updated', 'Table', `${iris.name} (150 rows)`,
        'Data storage', 'None', 'Last run', 'Never', 'Applications', 'Features', 'Sepal.Length', 'Target', 'Species',
        'Method', 'XGBoost', 'Task', 'Tags']);
      expect(tagsAutofill(accordion.getPane('Details').root), 'off');
      await expectPaneText(accordion, 'Performance', [`Split seed: ${model.seed}`]);
      const performance = accordion.getPane('Performance').root;
      const metrics = gridWith(performance, 'Metric');
      expectArray(metrics.dataFrame.columns.names(), ['Metric', 'Train', 'Validation']);
      expect(metrics.dataFrame.getCol('Metric').toList().includes('Accuracy') && !metrics.props.allowEdit, true,
        'No Accuracy, or the metrics are editable');
      expect(headerTooltip(metrics, 'Validation').startsWith('Value on rows the model did not see'), true,
        'No Validation header tooltip');
      expect(gridTooltip(metrics, metrics.cell('Metric', 0)), METRIC_DESCRIPTIONS.accuracy);
      expectArray(Array.from(performance.querySelectorAll('ul > li'), (li) => li.textContent?.trim()),
        ['Validation: 5-fold cross-validation on 150 rows, stratified', `Split seed: ${model.seed}`]);
      expect(performance.querySelector('.forge-seed')?.nextElementSibling?.classList.contains('fa-copy'), true,
        'The copy icon is not right after the seed');
      measured.push(expectRowsFit('Performance', metrics));
      await expectPaneText(accordion, 'Activity', ['Not applied yet.']);
      await expectPaneText(accordion, 'Sharing', ['Not shared yet.', 'Share...']);
      expect(accordion.getPane('Sharing').root.textContent?.includes('ask an administrator'), false,
        'The author is told to ask for the Share permission');
      await expectPaneText(accordion, 'History', ['insert', 'promote']);

      const loaded = await loadModel(id);
      await applyAndRecord({model: loaded, table: iris, mapping: exactMapping(loaded.features, iris),
        batchSize: DEFAULT_BATCH_SIZE, missingValues: {mode: 'skip'}}, 'ui');
      await expectPaneText(accordion, 'Activity', ['1 application']);
      const activity = gridWith(accordion.getPane('Activity').root, 'Status');
      expectArray(activity.dataFrame.columns.names(), ['When', 'Who', 'Table', 'Rows', 'Prediction column', 'Status',
        'Source', 'Duration (ms)']);
      const first = (column: string) => activity.dataFrame.get(column, 0);
      expectArray([activity.dataFrame.rowCount, first('Who'), first('Status'), first('Prediction column')],
        [1, DG.User.current().login, 'completed', 'Species (predicted)']);
      expect(headerTooltip(activity, 'Status'), 'Outcome of the application.');
      expect(gridTooltip(activity, activity.cell('Status', 0)), 'Completed: the column was added.');
      expect(gridTooltip(activity, activity.cell('Source', 0)), 'The Apply dialog or the catalog.');
      measured.push(expectRowsFit('Activity', activity));
      await expectGridsFillPanes(host, {Performance: metrics, Activity: activity});
      await expectPaneText(accordion, 'Details', ['Applications1', `Last run${dayjs().format('YYYY-MM-DD')}`]);
    } finally {
      host.remove();
    }
    return measured.join('; ');
  }, {timeout: TIMEOUT});

  test('the Sharing pane is read again in place', async () => {
    const row = ForgeModelHandler.rowOf(await forgeDb.models.get(sharedModel().id));
    const host = ui.div([ui.loader()]);
    const entity = grok.dapi.getEntities([row.id]).then((found) => found[0] ?? null);
    await refreshSharing(host, row, entity);
    await refreshSharing(host, row, entity);
    expect(host.children.length, 1, 'The pane keeps an old reading');
    expect(host.textContent?.includes('Not shared yet.') && host.textContent.includes('Share...'), true,
      host.textContent ?? '');
  });

  test('Edit model is prefilled and OK writes the name, description and tags', async () => {
    const name = `forge-test-model-${Date.now()}`;
    const id = await insertModelRow({name, description: 'd', tags: 'a'});
    // The test bundle's writes fire its own modelsChanged, which the follower _initForge started does not hear.
    const follow = ForgeModelHandler.followChanges();
    let opened: DG.Dialog | undefined;
    try {
      const model = await forgeDb.models.get(id);
      grok.shell.o = ForgeModelHandler.rowOf(model);
      const dialog = editModelDialog(model).show();
      opened = dialog;
      expect(dialog.input('Name').value, name);
      expect(dialog.input('Description').value, 'd');
      expectArray(tagsOfInput(dialog.input('Tags')), ['a']);
      expect(tagsAutofill(dialog.root), 'off');
      const ok = dialog.getButton('OK');
      dialog.input('Name').value = ' ';
      await awaitCheck(() => isDisabled(ok), 'OK is enabled without a name', WAIT_MS);
      dialog.input('Name').value = `${name}-edited`;
      dialog.input('Description').value = 'edited';
      dialog.input('Tags').value = ['a', 'b'];
      await awaitCheck(() => !isDisabled(ok), 'OK stays disabled with a name', WAIT_MS);
      ok.click();
      await awaitCheck(() => !dialog.root.isConnected, 'The dialog stays open', WAIT_MS);
      const row = await readUntil(() => forgeDb.models.get(id), (r) => r.name === `${name}-edited`);
      expect(row.name, `${name}-edited`);
      expect(row.description, 'edited');
      expect(row.tags, 'a, b');
      const shownName = () => {
        const shown: unknown = grok.shell.o;
        return shown instanceof DG.DomainRow ? `${shown.id === id} ${shown.displayName}` : `${shown}`;
      };
      expect(await readUntil(shownName, (shown) => shown === `true ${name}-edited`), `true ${name}-edited`,
        'The context panel does not show the edited model');
    } finally {
      follow.unsubscribe();
      opened?.close();
      await forgeDb.models.delete(id);
    }
  });

  test('a typed tag is saved', async () => {
    let savedTags: string[] | undefined;
    const save = saveModelDialog('forge-test-typed', async ({tags}) => {
      savedTags = tags;
    }).show();
    try {
      await typeTag(save.input('Tags').root, 'demo');
      expect(save.root.isConnected, true, `Enter in the Tags box submitted Save model with tags ${savedTags}`);
      const value: unknown = save.input('Tags').value;
      expectArray(tagsOfInput(save.input('Tags')), ['demo']);
      expect(Array.isArray(value), true, `The Tags value is ${typeof value}: ${value}`);
      save.getButton('OK').click();
      await awaitCheck(() => savedTags !== undefined, 'OK did not save', WAIT_MS);
      expectArray(savedTags ?? [], ['demo']);
    } finally {
      save.close();
    }

    const id = await insertModelRow({name: `forge-test-model-${Date.now()}-typed`});
    let opened: DG.Dialog | undefined;
    try {
      const edit = editModelDialog(await forgeDb.models.get(id)).show();
      opened = edit;
      await typeTag(edit.input('Tags').root, 'demo');
      expect(edit.root.isConnected, true, 'Enter in the Tags box submitted Edit model');
      expectArray(tagsOfInput(edit.input('Tags')), ['demo']);
      edit.getButton('OK').click();
      const storedTags = (tags: string) => readUntil(() => forgeDb.models.get(id), (row) => row.tags === tags);
      const edited = await storedTags('demo');
      expect(edited.tags, 'demo');

      // Two chips in a row in Details: the second change comes while the first one is being written.
      const accordion = modelAccordion(ForgeModelHandler.rowOf(edited), edited);
      document.body.append(accordion.root);
      const dialogs = new Set(DG.Dialog.getOpenDialogs().map((d) => d.root));
      const newDialogs = () => DG.Dialog.getOpenDialogs().filter((d) => !dialogs.has(d.root));
      try {
        await expectPaneText(accordion, 'Details', ['Tags']);
        const details = accordion.getPane('Details').root;
        await typeTag(details, 'x');
        await typeTag(details, 'y');
        const written = await storedTags('demo, x, y');
        expectArray([written.tags, written.version, newDialogs().length], ['demo, x, y', edited.version + 2, 0]);
        await recordApplication({model_id: id, row_count: 1, source: 'api', status: 'completed'});
        await expectPaneText(accordion, 'Details', ['Applications1']);
        await typeTag(details, 'z');
        expectArray([(await storedTags('demo, x, y, z')).version, newDialogs().length], [edited.version + 3, 0]);
      } finally {
        for (const dialog of newDialogs())
          dialog.close();
        accordion.root.remove();
      }
    } finally {
      opened?.close();
      await deleteModel(id);
    }
  });

  test('Save model has Name, Description and Tags', async () => {
    const dialog = saveModelDialog('Iris model', async () => {}).show();
    try {
      for (const caption of ['Name', 'Description', 'Tags'])
        expect(dialog.input(caption).caption, caption);
      expectArray(tagsOfInput(dialog.input('Tags')), []);
      const box = dialog.input('Tags').root.querySelector('input.d4-tags-selector-input');
      expect(box instanceof HTMLInputElement && box.placeholder, 'Type a tag and press Enter');
      expect(tagsAutofill(dialog.root), 'off');
    } finally {
      dialog.close();
    }
  });
});

// One category, split for its length.
category('UI: Catalog', () => {
  test('the model commands', async () => {
    expectArray(MODEL_ACTIONS.map((a) => a.name), ['Apply...', 'Download']);
    expectArray(CATALOG_ACTIONS.map((a) => a.name),
      ['Apply...', 'Edit model...', 'Share...', 'Download', 'Delete model']);
    const download = MODEL_ACTIONS[1];
    const saved = ForgeModelHandler.rowOf(await forgeDb.models.get(sharedModel().id));
    expect(download.isApplicable?.(saved), true);
    const manual = ForgeModelHandler.rowOf({id: DG.Utils.uuid4(), name: 'forge-test-manual'});
    expect(download.isApplicable?.(manual), false);

    await grok.functions.call('Forge:_initForge');
    const actions = ui.contextActions(saved);
    document.body.append(actions);
    const labels = () => Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label'),
      (label) => label.textContent?.trim() ?? '');
    const count = (name: string) => labels().filter((label) => label === name).length;
    try {
      actions.click();
      await awaitCheck(() => labels().includes('Apply...'), `The row's menu has no Apply...: ${labels()}`, WAIT_MS);
      expectArray(['Apply...', 'Download', 'Edit model...', 'Delete model'].map(count), [1, 1, 0, 0]);
    } finally {
      for (const popup of Array.from(document.querySelectorAll('.d4-menu-popup')))
        popup.remove();
      actions.remove();
    }
  });

  test('openShareDialog opens the platform\'s Share dialog of the model and returns', async () => {
    const name = `forge-test-model-${Date.now()}-share`;
    const id = await insertModelRow({name});
    const shareDialog = () => DG.Dialog.getOpenDialogs().find((d) => d.title.includes(name));
    try {
      await openShareDialog(id);
      await awaitCheck(() => shareDialog() !== undefined, 'No Share dialog of the model', WAIT_MS);
      expect(shareDialog()?.title.startsWith('Share'), true, shareDialog()?.title);
    } finally {
      shareDialog()?.root.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', keyCode: 27, bubbles: true}));
      const isClosed = await awaitCheck(() => shareDialog() === undefined, '', WAIT_MS).then(() => true, () => false);
      shareDialog()?.close();
      await forgeDb.models.delete(id);
      expect(isClosed, true, 'Escape does not close the Share dialog');
    }
  });

  test('Forge:shareModel refuses a missing model without a dialog', async () => {
    const dialogCount = DG.Dialog.getOpenDialogs().length;
    const error = await grok.functions.call('Forge:shareModel', {model: 'forge-test-none'})
      .then(() => null, (e: unknown) => e);
    expect(errorMessage(error).includes('No Forge model \'forge-test-none\' is available to you.'), true,
      String(error));
    expect(DG.Dialog.getOpenDialogs().length, dialogCount, 'A dialog opened');
  });
});

let paneModel: {id: string; name: string; iris: DG.DataFrame} | undefined;
const paneFixture = () => savedFixture(paneModel);

/** The line of the table pane [root] that names [modelName]. */
function paneLine(root: HTMLElement, modelName: string): HTMLElement | undefined {
  return Array.from(root.querySelectorAll<HTMLElement>('.forge-model-line'))
    .find((line) => (line.textContent ?? '').includes(modelName));
}

/** The model name and the run icon of [modelName]'s line in the table pane [root]; fails without them. */
function paneParts(root: HTMLElement, modelName: string): {link: HTMLElement; run: HTMLElement} {
  const line = paneLine(root, modelName);
  const link = line?.querySelector<HTMLElement>('.d4-link-label');
  const run = runIconOf(line);
  if (link === undefined || link === null || run === undefined || run === null)
    throw new Error('No model name or run icon');
  return {link, run};
}

/** The name label of an accordion title under [root] that reads [modelName]. */
function titleName(root: ParentNode, modelName: string): HTMLElement | undefined {
  return Array.from(root.querySelectorAll<HTMLElement>('.d4-accordion-title label'))
    .find((label) => label.textContent === modelName);
}

const trainButtonOf = (root: HTMLElement) =>
  Array.from(root.querySelectorAll<HTMLElement>('button')).find((b) => b.textContent?.trim() === 'Train');
const runIconOf = (line: HTMLElement | undefined) => line?.querySelector<HTMLElement>('.grok-icon.svg-model-run');

category('UI: Table panel', () => {
  before(async () => {
    useTestTuneStorage();
    const iris = await openIris();
    const {id, name} = await saveIrisModel(iris);
    paneModel = {id, name, iris};
  });

  after(async () => {
    clearTestTuneStorage();
    if (paneModel === undefined)
      return;
    if (isOpen(paneModel.iris))
      grok.shell.closeTable(paneModel.iris);
    await deleteModel(paneModel.id);
  });

  test('tableModelsPanel is a panel of a table', async () => {
    const panel = DG.Func.find({package: 'Forge', name: 'tableModelsPanel'})[0];
    expect(panel?.friendlyName, PANEL_TABLE_MODELS);
    expect(panel?.options['role'], 'panel');
    expect(panel?.options['condition'] === undefined, true, `A condition: ${panel?.options['condition']}`);
    expectArray(panel?.inputs.map((p) => p.propertyType) ?? [], [DG.TYPE.DATA_FRAME]);
  });

  test('the pane lists the fitting models as the built-in pane does, then Train', async () => {
    const {name, iris} = paneFixture();
    const pane = tableModelsPane(iris);
    await awaitCheck(() => paneLine(pane, name) !== undefined && trainButtonOf(pane) !== undefined,
      `The pane does not list the model: ${pane.textContent}`, WAIT_MS);
    const line = paneLine(pane, name);
    expect(runIconOf(line)?.classList.contains('d4-entity-markup-icon'), true, 'The line has no run icon');
    expect(line?.querySelector('.svg-model') === null, true, 'The line has a model icon');
    expect(line?.querySelector('.d4-link-label')?.textContent, name);
    expect(line?.lastElementChild?.classList.contains('grok-icon'), true, 'The line ends without the menu arrow');
    expect(trainButtonOf(pane)?.classList.contains('ui-btn-raised'), true, 'Train is not a big button');
    expect(pane.textContent?.includes('No model fits this table.'), false);

    const xy = DG.DataFrame.fromColumns([DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'x', [1, 2]),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'y', [3, 4])]);
    const empty = tableModelsPane(xy);
    await awaitCheck(() => trainButtonOf(empty) !== undefined, 'No Train for a table no model fits', WAIT_MS);
    expect(empty.querySelector('.forge-model-line') === null, true, `A model fits x, y: ${empty.textContent}`);
    expect(empty.textContent?.includes('No model fits this table.'), true, empty.textContent ?? '');
  });

  test('the run icon opens Apply predictive model with the model and the table', async () => {
    const {name, iris} = paneFixture();
    const tableView = grok.shell.addTableView(iris);
    const shown = new Set(DG.Dialog.getOpenDialogs().map((d) => d.root));
    const opened = () => DG.Dialog.getOpenDialogs().find((d) => !shown.has(d.root));
    try {
      const pane = tableModelsPane(iris);
      await awaitCheck(() => paneLine(pane, name) !== undefined, 'The pane does not list the model', WAIT_MS);
      runIconOf(paneLine(pane, name))?.click();
      await awaitCheck(() => opened() !== undefined, 'The run icon opens no dialog', WAIT_MS);
      const dialog = opened();
      expect(dialog?.title, 'Apply predictive model');
      expect(dialog?.input('Model').value, name);
      expect(dialog?.input('Table').value?.name, iris.name);
    } finally {
      opened()?.close();
      tableView.close();
    }
  });

  test('in the page, the platform\'s row CSS applies: the run icon is 12 px and centred on the name', async () => {
    const {name, iris} = paneFixture();
    const pane = tableModelsPane(iris);
    document.body.append(pane);
    try {
      await awaitCheck(() => paneLine(pane, name) !== undefined, 'The pane does not list the model', WAIT_MS);
      const {link, run} = paneParts(pane, name);
      expect(run.offsetHeight, 12, 'The run icon is not the built-in pane\'s 12 px');
      const middle = (e: Element) => {
        const r = e.getBoundingClientRect();
        return (r.top + r.bottom) / 2;
      };
      expect(Math.abs(middle(run) - middle(link)) <= 1, true,
        `The run icon's centre ${middle(run)} is off the name's ${middle(link)}`);
    } finally {
      pane.remove();
    }
  });

  test('the model name has the model\'s tooltip and the model\'s menu', async () => {
    const {name, iris} = paneFixture();
    const pane = tableModelsPane(iris);
    document.body.append(pane);
    const labels = () => Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label'),
      (label) => label.textContent?.trim() ?? '');
    try {
      await awaitCheck(() => paneLine(pane, name) !== undefined, 'The pane does not list the model', WAIT_MS);
      const {link, run} = paneParts(pane, name);
      expect(await enterTooltip(run), 'Run prediction');
      ui.tooltip.hide();
      const r = link.getBoundingClientRect();
      link.dispatchEvent(new MouseEvent('mouseenter', {clientX: r.left + 2, clientY: r.top + 2}));
      await awaitCheck(() => [name, 'Predicts ', 'Sepal.Length', 'Created ', ' by ']
        .every((t) => (ui.tooltip.root.textContent ?? '').includes(t)) &&
        !(ui.tooltip.root.textContent ?? '').includes('Trained on'),
      `The tooltip is not the model's: ${ui.tooltip.root.textContent}`, WAIT_MS);
      ui.tooltip.hide();

      link.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, cancelable: true, clientX: r.left + 2,
        clientY: r.top + 2}));
      await awaitCheck(() => labels().includes('Apply...') && labels().includes('Download'),
        `The name's menu lacks the model commands: ${labels()}`, WAIT_MS);
    } finally {
      ui.tooltip.hide();
      for (const popup of Array.from(document.querySelectorAll('.d4-menu-popup')))
        popup.remove();
      pane.remove();
    }
  });

  test('the panel function\'s pane shows the model\'s tooltip on a name, with several tables open', async () => {
    const {name, iris} = paneFixture();
    const tables = [iris, await openDemog(), await openDemoTable(WINE, 'wine'),
      await openDemoTable('System:DemoFiles/cars.csv', 'cars')];
    const views = tables.map((t) => grok.shell.addTableView(t));
    const widget: unknown = await DG.Func.find({package: 'Forge', name: 'tableModelsPanel'})[0].apply({table: iris});
    if (!(widget instanceof DG.Widget))
      throw new Error('The panel function gives no widget');
    const pane = widget.root;
    document.body.append(pane);
    try {
      await awaitCheck(() => paneLine(pane, name) !== undefined, 'The pane does not list the model', WAIT_MS);
      const link = paneLine(pane, name)?.querySelector<HTMLElement>('.d4-link-label');
      if (link === undefined || link === null)
        throw new Error('No model name');
      ui.tooltip.hide();
      const r = link.getBoundingClientRect();
      link.dispatchEvent(new MouseEvent('mouseenter', {clientX: r.left + 3, clientY: r.top + 3}));
      await awaitCheck(() => ui.tooltip.isVisible && [name, 'Predicts ', 'Created ']
        .every((t) => (ui.tooltip.root.textContent ?? '').includes(t)),
      `No model tooltip on the name: ${ui.tooltip.root.textContent}`, WAIT_MS);
    } finally {
      ui.tooltip.hide();
      pane.remove();
      for (const view of views)
        view.close();
      for (const t of tables.slice(1))
        grok.shell.closeTable(t);
    }
  }, {timeout: TIMEOUT});

  test('Train opens the Train view on the table', async () => {
    const {iris} = paneFixture();
    const pane = tableModelsPane(iris);
    await awaitCheck(() => trainButtonOf(pane) !== undefined, 'The pane has no Train', WAIT_MS);
    const views = new Set(Array.from(grok.shell.views));
    trainButtonOf(pane)?.click();
    const trainView = () => Array.from(grok.shell.views).find((v) => v instanceof TrainView && !views.has(v));
    await awaitCheck(() => trainView() !== undefined, 'Train opens no Train view', WAIT_MS);
    const view = trainView();
    if (!(view instanceof TrainView))
      throw new Error('No Train view');
    try {
      expect(view.tableInput.value?.name, iris.name);
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: TIMEOUT});

  test('the catalog has a description for the AI assistant', async () => {
    const view = await ForgeApp.create();
    try {
      expect(view.aiDescription?.includes('Forge'), true, `aiDescription: ${view.aiDescription}`);
    } finally {
      view.close();
    }
  });
});

/** The **Runs** grid of [view]; fails without one. */
function runsGrid(view: TrainView): DG.Grid {
  const grid = view.runsPanel.grid;
  if (grid === null)
    throw new Error('No Runs grid');
  return grid;
}

const runsColumn = (view: TrainView, name: string) => runsGrid(view).dataFrame.getCol(name).toList();

/** Clicks the **#** cell of [row] in the Runs [grid], as the user does. */
function clickRunsRow(grid: DG.Grid, row: number): void {
  const cell = grid.cell('#', row).bounds;
  const canvas = grid.overlay.getBoundingClientRect();
  const at = {bubbles: true, cancelable: true, button: 0, clientX: canvas.left + cell.x + cell.width / 2,
    clientY: canvas.top + cell.y + cell.height / 2};
  for (const type of ['mousedown', 'mouseup', 'click'])
    grid.overlay.dispatchEvent(new MouseEvent(type, at));
}

category('UI: Runs', () => {
  test('every recorded training is a row, newest first and current; a restored run trains again', async () => {
    const iris = await openIris();
    const view = await openTrainView(iris);
    try {
      expect(view.runs.length, 1);
      expect(isShown(view.runsPanel.root), true, 'The Runs panel is hidden');
      const grid = runsGrid(view);
      expect(grid.dataFrame.rowCount, 1);
      expectArray([grid.dataFrame.get('#', 0), grid.dataFrame.get('Method', 0), grid.dataFrame.get('Status', 0)],
        [1, 'XGBoost', 'completed']);
      expect(typeof grid.dataFrame.get('Accuracy (validation)', 0), 'number');
      expect(grid.dataFrame.columns.names().includes('AUC-ROC (validation)'), true,
        grid.dataFrame.columns.names().join(', '));
      expect(grid.dataFrame.columns.names().includes('MCC (validation)'), true,
        grid.dataFrame.columns.names().join(', '));
      expect(grid.dataFrame.currentRowIdx, 0);
      expect(grid.props.showCurrentRowIndicator, true);
      expect(gridTooltip(grid, grid.cell('Method', 0)).includes('iterations: 20'), true,
        gridTooltip(grid, grid.cell('Method', 0)));
      expect(gridTooltip(grid, grid.cell('Features', 0)), MEASUREMENTS.join('\n'));
      const accuracyTooltip = gridTooltip(grid, grid.cell('Accuracy (validation)', 0));
      expect(/^Train (\d\.\d{3}|-?\d\.\de-\d\d)$/.test(accuracyTooltip), true, accuracyTooltip);

      const firstShown = shownRun(view);
      expect(firstShown !== null && firstShown === view.lastTraining?.runId, true, `Shown run ${firstShown}`);
      setHyperparameter(view, 'iterations', 30);
      await awaitTrained(view);
      expect(shownRun(view) !== firstShown && shownRun(view) === view.lastTraining?.runId, true, 'The shown run');
      expectArray(runsColumn(view, '#'), [2, 1]);
      expect(runsGrid(view).dataFrame.currentRowIdx, 0);

      view.restoreRun(view.runs[1]);
      expect(view.hyperparameterInputs.get('iterations')?.value, 20);
      await awaitTrained(view);
      expect(view.runs.length, 3);
      expectArray(runsColumn(view, '#'), [3, 2, 1]);
      expect(view.hyperparameterInputs.get('iterations')?.value, 20);
      expect(view.runs[0].hyperparameters['iterations'], 20);
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: TIMEOUT});

  test('a click on a row whose inputs differ sets them again and trains; otherwise and on a move nothing', async () => {
    const iris = await openIris();
    const view = await openTrainView(iris);
    try {
      setHyperparameter(view, 'iterations', 30);
      await awaitTrained(view);
      runsGrid(view).dataFrame.currentRowIdx = 1;
      await delay(CHECK_DELAY_MS * 3);
      expect(view.hyperparameterInputs.get('iterations')?.value === 30 && view.runs.length === 2, true,
        'Moving the current row restored a run');

      clickRunsRow(runsGrid(view), 1);
      await awaitCheck(() => view.hyperparameterInputs.get('iterations')?.value === 20,
        'A click on a row does not set its inputs', WAIT_MS);
      await awaitTrained(view);
      expectArray(runsColumn(view, '#'), [3, 2, 1]);
      expect(runsGrid(view).dataFrame.currentRowIdx, 0);

      clickRunsRow(runsGrid(view), 0);
      await delay(CHECK_DELAY_MS * 3);
      expect(!view.isTraining && view.runs.length === 3, true, 'A click on a row with the form\'s inputs trained');

      view.featuresInput.value = [];
      await delay(CHECK_DELAY_MS * 3);
      expect(view.runs.length, 3, 'An empty Features trained');
      clickRunsRow(runsGrid(view), 0);
      await awaitCheck(() => view.featuresInput.value.length === 4, 'A click does not set the features again', WAIT_MS);
      await awaitTrained(view);
      await awaitCheck(() => view.runs.length === 4, 'A click on a row with other inputs does not train', TIMEOUT);
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: TIMEOUT});

  test('a click on the current run after a threshold move restores the threshold without training', async () => {
    const iris = await twoSpeciesIris();
    const view = await openTrainView(iris);
    const {predictProbability, threshold} = view.outputInputs;
    try {
      predictProbability.value = true;
      await delay(CHECK_DELAY_MS * 3);
      await awaitTrained(view);
      const [count, shown] = [view.runs.length, shownRun(view)];
      threshold.value = 0.7;
      await delay(CHECK_DELAY_MS * 3);
      expect(shown !== null && shownRun(view) === shown, true, `A re-cut shows ${shownRun(view)}, not ${shown}`);
      clickRunsRow(runsGrid(view), 0);
      await delay(CHECK_DELAY_MS * 3);
      expect(predictProbability.value, view.runs[0].steps.predictProbability);
      expect(`${view.runs[0].steps.predictProbability} ${threshold.value}`, 'false 0.7',
        'Unchecking by a run without Predict probability did not keep the hidden threshold');
      expect(!view.isTraining && view.runs.length === count && shownRun(view) === shown, true, 'The click trained');
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: TIMEOUT});

  test('a run whose columns the table lost sets nothing', async () => {
    const iris = await openIris();
    const view = await openTrainView(iris);
    try {
      const first = view.runs[0];
      iris.columns.remove('Petal.Width');
      await delay(CHECK_DELAY_MS * 2);
      await awaitCheck(() => !view.isTraining, 'The view keeps training', TIMEOUT);
      const count = view.runs.length;
      const features = view.featuresInput.value.map((c) => c.name);
      view.restoreRun(first);
      const balloon = await readUntil(() => balloonWith('The table no longer has the columns of this run'),
        (b) => b !== null);
      expect(balloon?.textContent?.includes('The table no longer has the columns of this run: Petal.Width.'), true,
        balloon?.textContent ?? '');
      await delay(CHECK_DELAY_MS * 3);
      expect(view.isTraining, false, 'A run with a missing column trains');
      expect(view.runs.length, count);
      expectArray(view.featuresInput.value.map((c) => c.name), features);
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: TIMEOUT});

  test('a Table change empties Runs and hides it', async () => {
    const iris = await openIris();
    const demog = grok.shell.addTable(await openDemog());
    const view = await openTrainView(iris);
    try {
      view.tableInput.value = demog;
      expect(view.runs.length, 0);
      expect(isShown(view.runsPanel.root), false, 'Runs stays shown after a Table change');
      await awaitTrained(view);
      expectArray(view.runs.map((run) => `${run.number} ${run.targetName}`), [`1 ${view.targetInput.value?.name}`]);
    } finally {
      await closeTrainView(view, [iris, demog]);
      grok.shell.closeTable(demog);
    }
  }, {timeout: TIMEOUT});
});

const SEARCH_TIMEOUT = 180000;

/** Switches [view] to Auto-tune with Search options expanded. */
function autoTune(view: TrainView): TuneSettings {
  const tune = view.tuneSettings;
  tune.setMode(AUTO_TUNE);
  tune.setExpanded(true);
  return tune;
}

const statusOf = (view: TrainView) => view.tuneSettings.statusText.textContent;

/** The row of the setting [name] in Search options; fails without one. */
function rowOf(tune: TuneSettings, name: string): TuneRow {
  const row = tune.row(name);
  if (row === undefined)
    throw new Error(`No row for ${name}`);
  return row;
}

/** Limit Trials, [n] of them. */
function stopAfterTrials(limit: LimitInputs, n: number): void {
  limit.kind.value = TRIALS;
  limit.trials.value = n;
}

/** Types [text] into [field], as the user does. */
function typeInto(field: HTMLInputElement, text: string): void {
  field.value = text;
  field.dispatchEvent(new Event('input', {bubbles: true}));
}

/** Picks [value] in [select], as the user does. */
function pick(select: HTMLSelectElement, value: string): void {
  select.value = value;
  select.dispatchEvent(new Event('change', {bubbles: true}));
}

/** The box of the text inside [element], as a DOM Range over its text nodes gives it (not the boxes around them). */
function textRect(element: HTMLElement): DOMRect {
  const walker = document.createTreeWalker(element, NodeFilter.SHOW_TEXT);
  const texts: Text[] = [];
  for (let node = walker.nextNode(); node !== null; node = walker.nextNode()) {
    if ((node.textContent ?? '').trim() !== '')
      texts.push(node as Text);
  }
  const range = document.createRange();
  range.setStart(texts[0], 0);
  range.setEnd(texts[texts.length - 1], texts[texts.length - 1].length);
  return range.getBoundingClientRect();
}

/** The font, color and alignment of [e]'s text. */
function textStyleOf(e: Element): string {
  const s = getComputedStyle(e);
  return [s.font, s.color, s.textAlign].join(' ');
}

/** The label [text] of a Search options [grid]. */
const gridLabelOf = (grid: HTMLElement, text: string) =>
  Array.from(grid.querySelectorAll<HTMLElement>('.forge-tune-label')).find((l) => l.textContent === text)!;

/** What the user sees of a Search options block's status line, toggle, Metric, Limit and Max time rows (and Max trials
 * once switched) and its [button]'s bar; the rows' heights and the gaps between them. Leaves the limit at Time. */
function searchBlockLines(block: {statusText: HTMLElement; toggle: HTMLElement; grid: HTMLElement; limit: LimitInputs},
  button: HTMLElement): string[] {
  const at = (n: number) => Math.round(n);
  block.limit.kind.value = TIME;
  const [status, toggleText] = [block.statusText, block.toggle.querySelector('span')!].map(textRect);
  const chevron = Array.from(block.toggle.querySelectorAll<HTMLElement>('.grok-icon')).find(isShown)!;
  const [metric, stop, time] = ['Metric', 'Limit', 'Max time'].map((text) => gridLabelOf(block.grid, text));
  const {kind, minutes, trials} = block.limit;
  const [radios, field] = [kind.input, minutes.input].map((e) => e.getBoundingClientRect());
  const lines = [
    `status ${textStyleOf(block.statusText)} at ${at(status.left)}`,
    `toggle ${textStyleOf(block.toggle)}, chevron at ${at(chevron.getBoundingClientRect().left)}, ` +
    `${at(toggleText.top - status.bottom)} px under the status`,
    ...[metric, stop, time].map((label) => `${label.textContent} ${textStyleOf(label)} ending at ` +
      `${at(label.getBoundingClientRect().right)}, ${at(label.getBoundingClientRect().height)} high`),
    `Limit ${at(radios.top - stop.getBoundingClientRect().top)} below its label, ${at(radios.left)}-` +
    `${at(radios.right)}, ${at(time.getBoundingClientRect().top - stop.getBoundingClientRect().bottom)} px ` +
    'above Max time',
    `Max time ${textStyleOf(minutes.input)} at ${at(field.left)}-${at(field.right)}, ${at(field.height)} high`,
    `bar ${button.className} ${textStyleOf(button)} ending at ${at(textRect(button).right)}`];
  kind.value = TRIALS;
  const trialsField = trials.input.getBoundingClientRect();
  lines.push(`Max trials ${textStyleOf(trials.input)} at ${at(trialsField.left)}-${at(trialsField.right)}, ` +
    `${at(trialsField.height)} high`);
  kind.value = TIME;
  return lines;
}

/** The option [caption] of the radio [input]. */
const radioOption = (input: DG.InputBase, caption: string) =>
  Array.from(input.root.querySelectorAll<HTMLElement>('.ui-radio-button')).find((o) => o.textContent === caption);

/** Starts a search of 1000 trials with [start] and cancels it from its own task-bar entry (`<action> (<p>%)`, p the
 * finished share) once at least 1% finished, so some completed; checks the best so far trains. */
async function cancelSearchOnceSettled(view: TrainView, action: string,
  start: (maxTrials: number) => Promise<void>): Promise<void> {
  // Only the search's own entry counts, as an earlier one may still be in the task bar.
  const entries = () => Array.from(document.querySelectorAll('.d4-task-bar-entry'));
  const earlier = new Set(entries());
  const searching = start(1000);
  const pattern = new RegExp(`${action} \\((\\d+)%\\)`);
  const finished = (e: Element) => Number(pattern.exec(e.textContent ?? '')?.[1] ?? 0);
  const entry = () => entries().find((e) => !earlier.has(e) && finished(e) >= 1);
  await awaitCheck(() => entry() !== undefined, 'The search launches no trial after its first ones', TIMEOUT);
  cancelTaskBarEntry(entry());
  await awaitCheck(() => !view.isTraining, 'The cancelled search does not end', TIMEOUT);
  await searching;
  expect(view.lastTraining?.result.search?.stopped, 'cancelled');
  expect(balloonWith('Search cancelled') !== null, true, 'No balloon on a cancel');
}

/** The **Search** grid of [view]; fails without one. */
function searchGrid(view: TrainView): DG.Grid {
  const grid = view.searchPanel.grid;
  if (grid === null)
    throw new Error('No Search grid');
  return grid;
}

/** The bars of [chart], top to bottom, are [grid]'s rows: the top one the chosen method's. */
function expectChartOrder(chart: DG.Viewer, grid: DG.Grid): void {
  const {hitAreas} = chart.getWidgetStatus();
  const bars = Object.keys(hitAreas).filter((k) => k.startsWith('bar ')).sort((a, b) => hitAreas[a].y - hitAreas[b].y);
  expectArray(bars, grid.dataFrame.getCol('Method').toList().slice(0, bars.length).map((m) => `bar ${m}`));
  expect(grid.dataFrame.get('Chosen', 0), true, 'The top bar is not the chosen method\'s');
}

/** [grid] ranks [leaders] (of a Find best model by accuracy) as the chart: by the score less its SE as stored, best
 * first, a tie by the method's name; the chosen one first and current; the tooltips name the rule. */
function expectSearchOrder(grid: DG.Grid, leaders: Leader[]): void {
  const frame = grid.dataFrame;
  const shown = (l: Leader) => Math.fround(l.score - l.se);
  const sorted = [...leaders].sort((a, b) => shown(b) - shown(a) ||
    (a.method < b.method ? -1 : a.method > b.method ? 1 : 0));
  expectArray(frame.getCol('Method').toList().slice(0, sorted.length), sorted.map((l) => l.method));
  expectArray(frame.getCol('#').toList(), Array.from({length: frame.rowCount}, (_, i) => i + 1));
  const cautious = frame.getCol('Accuracy − SE');
  expectArray(sorted.map((l, i) => Math.abs(cautious.get(i) - (l.score - l.se)) < 1e-6), sorted.map(() => true));
  expect(`${sorted[0].chosen} ${frame.getCol('Chosen').get(0)} ${frame.currentRowIdx}`, 'true true 0',
    'The chosen one is not first');
  expectArray([headerTooltip(grid, '#'), headerTooltip(grid, 'Accuracy − SE'), headerTooltip(grid, 'Chosen'),
    gridTooltip(grid, grid.cell('Chosen', 0))], ['Rank by Accuracy − SE.', 'Accuracy minus its standard error ' +
    'over the folds: a cautious estimate of Accuracy. The method with the highest value is chosen.',
  'The method with the highest Accuracy − SE.', 'Chosen: the highest Accuracy − SE.']);
}

const searchRecord = (holdOut?: {rows: number; fraction: number}): SearchRecord => ({mode: 'tune',
  metric: 'accuracy', budgetSeconds: 60, seed: 1, trials: 12, completed: 10, pruned: 1, failed: 1, elapsedMs: 1000,
  stopped: 'budget', ...(holdOut === undefined ? {} : {holdOut}), leaders: [{method: 'XGBoost',
    hyperparameters: {iterations: 20}, score: 0.8, se: 0.01, foldScores: [0.8, 0.8, 0.8, 0.8, 0.8], trials: 12,
    chosen: true}]});

/** The current model card of [view]; fails without one. */
function currentCard(view: TrainView): HTMLElement {
  const card = view.modelCards.root.querySelector<HTMLElement>('.forge-card-current');
  if (card === null)
    throw new Error('No current model card');
  return card;
}

/** The pinned model cards of [view], top to bottom. */
const pinnedCards = (view: TrainView) =>
  Array.from(view.modelCards.root.querySelectorAll<HTMLElement>('.forge-card:not(.forge-card-current)'));
const cardNumbers = (cards: HTMLElement[]) => cards.map((card) => card.dataset.forgeRun ?? '');
/** The text of [card]'s element of the class [cls]; '' without one. */
const cardText = (card: HTMLElement, cls: string) => card.querySelector(`.${cls}`)?.textContent ?? '';
const isBest = (card: HTMLElement) => isShown(card.querySelector<HTMLElement>('.forge-card-best')!);
const pinIcon = (card: HTMLElement) => card.querySelector<HTMLElement>('.forge-card-pin')!;
const unpinIcon = (card: HTMLElement) => card.querySelector<HTMLElement>('.forge-card-unpin')!;
/** The icons [card] shows, as "+" and "trash". */
const visibleIcons = (card: HTMLElement) => Array.from(card.querySelectorAll<HTMLElement>('.forge-card-icon'))
  .filter((icon) => isShown(icon) && icon.style.visibility !== 'hidden')
  .map((icon) => icon.classList.contains('forge-card-pin') ? '+' : 'trash');
/** Each row of [table] as its cells' texts. */
const tableRows = (table: HTMLElement) =>
  Array.from(table.querySelectorAll('tr'), (row) => Array.from(row.cells, (cell) => cell.textContent ?? ''));
/** The validation value of [metric] of [run]. */
const validationOf = (run: SessionRun, metric: MetricId) => run.metrics?.validation[metric];

/** The numbers of [runs] (newest first) in the order of their pinned cards by [metric]: the best validation value as
 * shown (3 decimals) first, a run without one last, ties newest first. */
function sortedNumbers(runs: SessionRun[], metric: MetricId): string[] {
  const rank = (run: SessionRun) => {
    const value = validationOf(run, metric);
    return value === undefined ? -Number.MAX_VALUE : directionOf(metric) * roundTo3(value);
  };
  return [...runs].sort((a, b) => rank(b) - rank(a)).map((run) => `${run.number}`);
}

category('UI: Search', () => {
  test('Tune keeps Predict probability; the Metric lists Log loss and Brier score last', async () => {
    const iris = await twoSpeciesIris();
    const view = await openTrainView(iris);
    try {
      view.outputInputs.predictProbability.value = true;
      autoTune(view);
      await awaitCheck(() => !isDisabled(view.tuneButton), 'Tune stays disabled', WAIT_MS);
      expect(await enterTooltip(view.tuneButton), 'Find the best settings');
      const metrics = Array.from(view.tuneSettings.metric.options, (o) => o.textContent ?? '');
      expectArray(metrics.slice(-2), ['Log loss', 'Brier score']);
    } finally {
      await closeTrainView(view, [iris]);
      clearTestTuneStorage();
    }
  }, {timeout: TIMEOUT});

  test('Tune optimizes the winner\'s threshold for the chosen metric, unless it was edited by hand', async () => {
    const iris = await twoSpeciesIris();
    const view = await openTrainView(iris);
    const {predictProbability, threshold, thresholdMetric} = view.outputInputs;
    const options = () => view.lastTraining?.result.options;
    try {
      predictProbability.value = true;
      await awaitCheck(() => options()?.thresholdMetric === 'f1', 'Checking does not optimize', WAIT_MS);
      thresholdMetric.value = 'MCC';
      await view.startSearch('tune', {metric: 'accuracy', maxTrials: 2});
      await awaitTrained(view);
      const best = bestThresholdOf(view.lastTraining!.result, 'mcc');
      expect(view.lastTraining?.request.search !== undefined, true, 'The winner is not the shown training');
      expect(options()?.postprocessingInfo.includes('class-probabilities'), true,
        'The re-cut after Tune dropped Predict probability');
      expect(view.runs[0].steps.predictProbability, true, 'The winner\'s run row lost Predict probability');
      expect(options()?.thresholdMetric, 'mcc');
      expect(options()?.binaryClassificationThreshold, best);
      expect(threshold.value, best, 'The winner\'s threshold is not optimized');

      threshold.value = 0.65;
      await awaitCheck(() => options()?.binaryClassificationThreshold === 0.65, 'The hand edit does not re-cut',
        WAIT_MS);
      await view.startSearch('tune', {metric: 'accuracy', maxTrials: 2});
      await awaitTrained(view);
      expect(threshold.value, 0.65, 'A search replaced the hand edit');
      expect(options()?.binaryClassificationThreshold, 0.65);
      expect(options()?.thresholdMetric === undefined, true, 'A search optimized a hand edit');
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: SEARCH_TIMEOUT});

  test('Tune follows the selection; a method without settings has no Settings row', async () => {
    const iris = await openIris();
    const view = await openTrainView(iris);
    try {
      expect(view.tuneSettings.isAuto, false, 'Auto-tune is the default mode');
      expect(view.tuneButton.offsetParent === null, true, 'Tune is shown in Manual');
      autoTune(view);
      expect(view.tuneButton.textContent, 'Tune');
      expect(isDisabled(view.tuneButton), false, 'Tune is disabled on a valid selection');
      grok.userSettings.add(TuneSettings.storage, 'Naive Bayes', JSON.stringify({mode: AUTO_TUNE}));
      view.methodInput.value = 'Naive Bayes';
      await awaitTrained(view);
      const bayes = view.tuneSettings;
      expect(bayes.isAuto, false, 'Naive Bayes, which has no settings, is in Auto-tune');
      const shown = [bayes.mode.root, bayes.resetButton, bayes.statusText, view.tuneButton]
        .filter((e) => e.offsetParent !== null);
      expect(shown.length, 0, 'Naive Bayes, which has no settings, shows Settings, its reset, its status or Tune');
      expect(bayes.root.offsetHeight, 0, 'Naive Bayes leaves a gap where Settings was');
      view.methodInput.value = 'XGBoost';
      await awaitTrained(view);
      expect(view.tuneSettings.isAuto, true, 'XGBoost forgot Auto-tune');
      expect(view.tuneSettings.mode.root.offsetParent !== null, true, 'XGBoost does not show Settings again');
      grok.userSettings.delete(TuneSettings.storage, 'Naive Bayes');
      expect(isDisabled(view.tuneButton), false, 'Tune stays disabled for XGBoost');
      view.featuresInput.value = [];
      await awaitCheck(() => resultsState(view) === 'false/true', 'An empty Features is not a problem', WAIT_MS);
      expect(isDisabled(view.tuneButton), true, 'Tune is enabled without features');
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: TIMEOUT});

  test('each Search options row shows its setting\'s Manual tooltip, and a re-layout leaves none behind', async () => {
    const iris = await openIris();
    const view = await openTrainView(iris);
    try {
      const tune = autoTune(view);
      tune.togglePin('maxDepth');
      const caption = (name: string) => view.hyperparameterInputs.get(name)!.caption;
      const description = (name: string) => rowOf(tune, name).property.description ?? '';
      const wrong: string[] = [];
      for (const name of view.hyperparameterInputs.keys()) {
        const row = rowOf(tune, name);
        const controls = [row.label, ...row.from === undefined || row.isPinned ? [row.value] : [row.from, row.to!]];
        // A disabled control gets no mouse events, so its tooltip cannot show (its label's can).
        for (const control of controls.filter((c) => !(c as HTMLInputElement).disabled)) {
          const text = await enterTooltip(control as HTMLElement);
          if (text !== description(name))
            wrong.push(`${name} ${control.className || control.tagName}: "${text}"`);
        }
        const pin = await enterTooltip(row.pin);
        if (!pin.includes(`${caption(name)}:`))
          wrong.push(`${name} pin: "${pin}"`);
      }
      expectArray(wrong, []);
      expect(description('iterations') !== '', true, 'XGBoost Iterations has no description');
      const label = (text: string) => Array.from(tune.grid.querySelectorAll<HTMLElement>('.forge-tune-label'))
        .find((l) => l.textContent === text)!;
      const settingsTooltip = 'How the settings are chosen: by hand or by a search on this data';
      expectArray([await enterTooltip(tune.mode.captionLabel), await enterTooltip(radioOption(tune.mode, MANUAL)!),
        await enterTooltip(radioOption(tune.mode, AUTO_TUNE)!)], [settingsTooltip, settingsTooltip, settingsTooltip]);
      const {limit} = tune;
      limit.kind.value = TIME;
      const tooltips = [await enterTooltip(label('Metric')), await enterTooltip(tune.metric),
        await enterTooltip(label('Limit')), await enterTooltip(radioOption(limit.kind, TIME)!),
        await enterTooltip(radioOption(limit.kind, TRIALS)!), await enterTooltip(label('Max time')),
        await enterTooltip(limit.minutes.input), await enterTooltip(tune.toggle)];
      limit.kind.value = TRIALS;
      tooltips.push(await enterTooltip(label('Max trials')), await enterTooltip(limit.trials.input));
      const limitTooltip = 'What ends the search: a time or a trial limit';
      expectArray(tooltips, ['Metric the search optimizes', 'Metric the search optimizes', limitTooltip,
        limitTooltip, limitTooltip, 'Minutes the search may run', 'Minutes the search may run',
        'Hide search options', 'Maximum number of trials, pruned and failed ones included',
        'Maximum number of trials, pruned and failed ones included']);

      const [rate, depth] = [rowOf(tune, 'eta').label, rowOf(tune, 'maxDepth').label];
      rate.dispatchEvent(new MouseEvent('mouseenter'));
      await delay(400);
      rate.dispatchEvent(new MouseEvent('mouseleave'));
      depth.dispatchEvent(new MouseEvent('mouseenter'));
      await delay(400);
      expect(ui.tooltip.root.textContent, description('maxDepth'));
      depth.dispatchEvent(new MouseEvent('mouseleave'));
      rate.dispatchEvent(new MouseEvent('mouseenter'));
      await delay(400);
      tune.togglePin('eta');
      expect(ui.tooltip.isVisible, false, 'A re-layout leaves the tooltip of a hovered row');
    } finally {
      // Closing writes the pending save first.
      await closeTrainView(view, [iris]);
      grok.userSettings.delete(TuneSettings.storage, 'XGBoost');
    }
  }, {timeout: 120000});

  test('Search options are remembered per user and method; a stale item falls back to its default', async () => {
    const iris = await openIris();
    let view = await openTrainView(iris);
    try {
      const tune = autoTune(view);
      tune.togglePin('maxDepth');
      typeInto(rowOf(tune, 'iterations').from!, '15');
      tune.limit.minutes.value = 5;
      stopAfterTrials(tune.limit, 20);
      pick(tune.metric, 'f1');
      tune.flushSave();
      const stored = TuneSettings.storedRecord('XGBoost') as {[key: string]: any};
      expect(JSON.stringify([stored.mode, stored.expanded, stored.pinned, stored.ranges.iterations, stored.metric,
        stored.limit, stored.minutes, stored.trials]), JSON.stringify([AUTO_TUNE, true, ['maxDepth'],
        {low: 15, high: 100}, 'f1', TRIALS, 5, 20]));
      await closeTrainView(view, [iris]);

      view = await openTrainView(iris);
      const again = view.tuneSettings;
      expectArray([again.isAuto, again.isExpanded, rowOf(again, 'maxDepth').pin.getAttribute('aria-pressed'),
        rowOf(again, 'iterations').from!.value, again.limit.kind.value, again.limit.minutes.value,
        again.limit.trials.value, again.metric.value], [true, true, 'true', '15', TRIALS, 5, 20, 'f1']);
      await closeTrainView(view, [iris]);

      // The old time select's seconds: whole minutes, rounded up.
      grok.userSettings.add(TuneSettings.storage, 'XGBoost', JSON.stringify({mode: AUTO_TUNE,
        pinned: ['speed', 'eta'], ranges: {iterations: {low: 5, high: 500}, eta: {low: 0.1, high: 0.2}}, metric: 'r2',
        trials: 0, seconds: 90}));
      view = await openTrainView(iris);
      const stale = view.tuneSettings;
      expectArray([stale.isAuto, stale.isExpanded, rowOf(stale, 'eta').pin.getAttribute('aria-pressed'),
        rowOf(stale, 'iterations').from!.value, rowOf(stale, 'eta').from!.value, stale.limit.kind.value,
        stale.limit.minutes.value, stale.limit.trials.value, stale.metric.value],
      [true, false, 'true', '10', '0.1', TIME, 2, 50, 'balancedAccuracy']);
      const engine = view.lastTraining!.engine;
      const loaded = (value: unknown) => {
        grok.userSettings.add(TuneSettings.storage, 'XGBoost', JSON.stringify(value));
        const {limit, minutes, trials} = TuneSettings.load(engine, stale.space, hyperparametersOf(engine));
        return `${limit} ${minutes} ${trials}`;
      };
      expectArray([loaded({seconds: 30}), loaded({seconds: 900, trials: 70}), loaded({seconds: 7200}),
        loaded({limit: TRIALS, minutes: 61, trials: 10001}), loaded({limit: TRIALS, minutes: 3, trials: 9}),
        loaded({limit: TIME, minutes: 4}), loaded({mode: AUTO_TUNE}), loaded({limit: 'Steps'})],
      [`${TIME} 1 50`, `${TIME} 15 70`, `${TIME} 60 50`, `${TRIALS} 1 50`, `${TRIALS} 3 9`, `${TIME} 4 50`,
        `${TRIALS} 1 50`, `${TRIALS} 1 50`]);
      grok.userSettings.add(TuneSettings.storage, 'XGBoost', '{not json');
      expect(JSON.stringify(TuneSettings.load(engine, stale.space, hyperparametersOf(engine))),
        JSON.stringify({mode: MANUAL, expanded: false, pinned: [], ranges: {}, limit: TRIALS, minutes: 1, trials: 50}));
    } finally {
      await closeTrainView(view, [iris]);
      grok.userSettings.delete(TuneSettings.storage, 'XGBoost');
    }
  }, {timeout: 120000});

  test('a click on a Runs or Search row switches to Manual', async () => {
    const iris = await openIris();
    const view = await openTrainView(iris);
    try {
      setHyperparameter(view, 'iterations', 30);
      await awaitTrained(view);
      autoTune(view);
      await delay(CHECK_DELAY_MS * 3);
      clickRunsRow(runsGrid(view), 1);
      await awaitCheck(() => !view.tuneSettings.isAuto, 'A Runs row click stays in Auto-tune', WAIT_MS);
      await awaitTrained(view);
      autoTune(view);
      view.restoreLeader({method: 'XGBoost', hyperparameters: {iterations: 11}});
      expect(view.tuneSettings.isAuto, false, 'A Search row click stays in Auto-tune');
      await awaitTrained(view);
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: TIMEOUT});

  test('Tune trains the chosen settings in Manual as Auto-tuned, without the Search panel of one method', async () => {
    const iris = await openIris();
    const demog = grok.shell.addTable(await openDemog());
    const view = await openTrainView(iris);
    try {
      expect(isShown(view.searchPanel.root), false, 'Search is shown before a search');
      const tune = autoTune(view);
      const [cardRight, tuneRight] = [currentCard(view), view.tuneButton].map((e) => e.getBoundingClientRect().right);
      expect(tuneRight > 0 && cardRight <= tuneRight + 1, true, `The card ends at ${cardRight}, Tune at ${tuneRight}`);
      tune.togglePin('maxDepth');
      stopAfterTrials(tune.limit, 6);
      pick(tune.metric, 'accuracy');
      await view.startSearch('tune', tune.choice());
      expect(view.isTraining, false);
      const search = view.lastTraining?.result.search;
      expect(search?.mode, 'tune');
      expect(`${view.tuneSettings.isAuto} ${statusOf(view)}`, 'false Auto-tuned');
      expect(view.runs[0].settingsStatus, 'Auto-tuned');
      const card = currentCard(view);
      expect(cardText(card, 'forge-card-kind'), 'Auto-tuned settings');
      expect(card.querySelector('.d4-tab-host'), null, 'The card has tabs');
      expect(card.querySelector('.forge-tune-status-text'), null, 'The card shows a status');
      expectArray(search?.space?.pinned ?? [], ['maxDepth']);
      expect(search?.space?.searched.some((d) => d.name === 'maxDepth'), false, 'A pinned setting was searched');
      expect(isShown(view.searchPanel.root) || view.searchPanel.grid !== null, false, 'Search shows one method');
      const leader = search!.leaders[0];
      expect(`${search!.leaders.length} ${leader.method} ${leader.chosen}`, '1 XGBoost true');
      expectFloat(view.lastTraining!.result.metrics.validation.accuracy!, leader.score, 1e-9);
      for (const [name, value] of Object.entries(leader.hyperparameters)) {
        const shown: unknown = view.hyperparameterInputs.get(name)?.value;
        if (typeof value === 'number')
          expectFloat(Number(shown), value, 1e-6);
        else
          expect(shown, value, name);
      }
      expectArray(gridWith(view.resultsPane.root, 'Metric').dataFrame.columns.names(),
        ['Metric', 'Train', 'Validation']);
      expect(resultsText(view).includes('Search: Tune by Accuracy'), true, resultsText(view));
      expect(resultsText(view).includes('Validation after a search is optimistic'), true, resultsText(view));
      expect(shownRun(view) !== null && shownRun(view) === view.lastTraining?.runId, true, 'The shown run');
      const runTooltip = gridTooltip(runsGrid(view), runsGrid(view).cell('Method', 0));
      expect(runTooltip.startsWith(`Found by Tune: ${search!.trials} trials`), true, runTooltip);

      const tuned = view.hyperparameterInputs.get('iterations')!;
      const value = tuned.value;
      tuned.value = Number(value) === 50 ? 51 : 50;
      expect(statusOf(view), 'Customized');
      tuned.value = value;
      expect(statusOf(view), 'Auto-tuned');
      await awaitTrained(view);

      view.tableInput.value = demog;
      expect(isShown(view.searchPanel.root), false, 'Search stays shown after a Table change');
    } finally {
      await closeTrainView(view, [iris, demog]);
      grok.shell.closeTable(demog);
    }
  }, {timeout: SEARCH_TIMEOUT});

  test('Tune again continues the search: Tune says so, the trials and the limit add up; a change starts afresh',
    async () => {
      const iris = await openIris();
      const view = await openTrainView(iris);
      try {
        let tune = autoTune(view);
        stopAfterTrials(tune.limit, 3);
        await awaitCheck(() => !isDisabled(view.tuneButton), 'Tune stays disabled', WAIT_MS);
        expect(await enterTooltip(view.tuneButton), 'Find the best settings');
        await view.startSearch('tune', tune.choice());
        const first = view.lastTraining!.result.search!;
        tune = autoTune(view);
        await awaitCheck(() => !isDisabled(view.tuneButton), 'Tune stays disabled after a tune', WAIT_MS);
        expect(await enterTooltip(view.tuneButton), `Continue the search: ${first.trials} trials so far`);
        await view.startSearch('tune', tune.choice());
        const second = view.lastTraining!.result.search!;
        expect(`${second.trials > first.trials} ${second.trials <= first.trials + 3} ${second.maxTrials}`,
          'true true 6', `${first.trials} trials, then ${second.trials}`);

        tune = autoTune(view);
        const metric = tune.metric.value === 'f1' ? 'mcc' : 'f1';
        pick(tune.metric, metric);
        await awaitCheck(() => !isDisabled(view.tuneButton), 'Tune stays disabled with another metric', WAIT_MS);
        expect(await enterTooltip(view.tuneButton), 'Find the best settings');
        await view.startSearch('tune', tune.choice());
        const fresh = view.lastTraining!.result.search!;
        expect(`${fresh.trials <= 3} ${fresh.maxTrials} ${fresh.metric}`, `true 3 ${metric}`);
      } finally {
        await closeTrainView(view, [iris]);
        clearTestTuneStorage();
      }
    }, {timeout: SEARCH_TIMEOUT});

  test('a search after the threshold moved starts afresh, and neither Tune nor the task bar promises to continue',
    async () => {
      const iris = await twoSpeciesIris();
      const view = await openTrainView(iris);
      const {predictProbability, threshold} = view.outputInputs;
      const texts = new Set<string>();
      const observer = new MutationObserver(() => {
        for (const entry of Array.from(document.querySelectorAll('.d4-task-bar-entry')))
          texts.add(entry.textContent ?? '');
      });
      try {
        predictProbability.value = true;
        await awaitCheck(() => view.lastTraining?.result.options.thresholdMetric === 'f1', 'Checking does not optimize',
          WAIT_MS);
        const tune = autoTune(view);
        stopAfterTrials(tune.limit, 2);
        await view.startSearch('tune', tune.choice());
        await awaitTrained(view);
        const moved = threshold.value === 0.3 ? 0.35 : 0.3;
        threshold.value = moved;
        await awaitCheck(() => view.lastTraining?.result.options.binaryClassificationThreshold === moved,
          'The threshold does not re-cut', WAIT_MS);
        const again = autoTune(view);
        await awaitCheck(() => !isDisabled(view.tuneButton), 'Tune stays disabled', WAIT_MS);
        expect(await enterTooltip(view.tuneButton), 'Find the best settings');
        observer.observe(document.body, {subtree: true, childList: true, characterData: true});
        await view.startSearch('tune', again.choice());
        const fresh = view.lastTraining!.result.search!;
        expect(`${fresh.trials <= 2} ${fresh.maxTrials}`, 'true 2', 'The search was continued');
        expect([...texts].some((t) => t.startsWith('Continuing')), false, [...texts].join('; '));
      } finally {
        observer.disconnect();
        await closeTrainView(view, [iris]);
        clearTestTuneStorage();
      }
    }, {timeout: SEARCH_TIMEOUT});

  test('the status reads Auto-tuned only on the table, target and features of the last tune', async () => {
    const iris = await openIris();
    const view = await openTrainView(iris);
    try {
      await view.startSearch('tune', {metric: 'accuracy', maxTrials: 2});
      const tune = view.tuneSettings;
      expect(statusOf(view), 'Auto-tuned');
      const species = view.targetInput.value!;
      expect(species.name, 'Species');
      view.targetInput.value = iris.getCol('Petal.Width');
      expect(statusOf(view), 'Customized');
      expect(await enterTooltip(tune.statusText), `Tuned for Species on ${iris.name}; tune again for this data`);
      view.targetInput.value = species;
      expect(statusOf(view), 'Auto-tuned');
      expect(await enterTooltip(tune.statusText), 'Found by the search');
      const features = view.featuresInput.value;
      view.featuresInput.value = [...features].reverse();
      expect(statusOf(view), 'Auto-tuned', 'The order of the features counts');
      view.featuresInput.value = features.slice(1);
      expect(statusOf(view), 'Customized');
      expect(view.tuneSettings === tune, true, 'The method changed');
      await awaitTrained(view);
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: SEARCH_TIMEOUT});

  test('a hold-out search names the score column <Metric> (hold-out)', async () => {
    const iris = await openIris();
    const view = await openTrainView(iris);
    try {
      view.validationInputs.scheme.value = 'Hold-out';
      await awaitTrained(view);
      await chooseAuto(view);
      await view.startSearch('best', {metric: 'accuracy', maxTrials: 4, methods: ['LDA', 'Decision Tree']});
      const grid = searchGrid(view);
      expectArray(grid.dataFrame.columns.names(), ['#', 'Method', 'Accuracy (hold-out)', 'SE', 'Accuracy − SE',
        'Trials', 'Chosen']);
      expect(headerTooltip(grid, 'Accuracy (hold-out)').includes('On the hold out'), true,
        headerTooltip(grid, 'Accuracy (hold-out)'));
      expect(gridTooltip(grid, grid.cell('Accuracy (hold-out)', 0)).startsWith('Folds: '), true,
        gridTooltip(grid, grid.cell('Accuracy (hold-out)', 0)));
      expectArray([headerTooltip(grid, 'SE').startsWith('The standard error'), headerTooltip(grid, 'Trials')],
        [true, 'Every trial of the method, pruned and failed ones included.']);
      const run = view.runs[0];
      expect(ModelCards.kindLine(run), 'Auto-tuned model');
      if (visibleHyperparametersOf(run.engine).length > 0)
        expect(cardText(currentCard(view), 'forge-card-kind'), 'Auto-tuned model');
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: SEARCH_TIMEOUT});

  test('cancelling the search trains the best so far', async () => {
    const iris = await openIris();
    const view = await openTrainView(iris);
    try {
      await cancelSearchOnceSettled(view, 'Tuning XGBoost',
        (maxTrials) => view.startSearch('tune', {metric: 'accuracy', maxTrials}));
      expect(view.lastTraining?.result.search?.leaders.length, 1);
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: SEARCH_TIMEOUT});

  test('a search disables the form and ignores Runs and Search rows; each input\'s state comes back after it',
    async () => {
      const iris = await openIris();
      const view = await openTrainView(iris);
      try {
        setHyperparameter(view, 'iterations', 30);
        view.validationInputs.shuffle.value = false;
        await awaitTrained(view);
        const states = () => formStates(view);
        const before = states();
        expect(before.includes('Split seed false'), true, before.join(', '));
        const count = view.runs.length;
        const searching = view.startSearch('tune', {metric: 'accuracy', maxTrials: 6});
        expectArray(enabledInputs(view), []);
        // Tune is gated only in Auto-tune: in Manual it is hidden.
        expect(isDisabled(view.tuneButton) || view.tuneButton.offsetParent === null, true,
          'Tune is shown and enabled during the search');
        view.restoreRun(view.runs[1]);
        clickRunsRow(runsGrid(view), 1);
        view.restoreLeader({method: 'XGBoost', hyperparameters: {iterations: 11}});
        expectArray(enabledInputs(view), []);
        await searching;
        expect(view.lastTraining?.result.search?.mode, 'tune', 'The search was superseded');
        expect(view.runs.length, count + 1);
        // The winner's settings have their own rules.
        const rules = enabledRules(EngineRegistry.discover().find((e) => e.name === 'XGBoost')!) ?? {};
        const values = view.lastTraining!.result.hyperparameters;
        const settings = [...view.hyperparameterInputs].map(([name, input]) => `${input.caption} ` +
          `${isEnabled(rules, name, values, 'classification')}`);
        expectArray(states(), [...before.slice(0, before.length - settings.length), ...settings]);
      } finally {
        await closeTrainView(view, [iris]);
      }
    }, {timeout: SEARCH_TIMEOUT});

  test('a model with test and search shows the Test column, the bullets and Details', async () => {
    const iris = await openIris();
    const {id} = await saveIrisModel(iris, {metrics: {train: {accuracy: 0.9}, validation: {accuracy: 0.8},
      test: {accuracy: 0.85}}, search: searchRecord({rows: 1170, fraction: 0.2})});
    try {
      const model = await forgeDb.models.get(id);
      const summary = storedSummary(model);
      expect(summary?.search?.holdOut?.rows, 1170);
      const table = new MetricsTable(summary!);
      expectArray(gridWith(table.root, 'Metric').dataFrame.columns.names(), ['Metric', 'Train', 'Validation', 'Test']);
      const grid = gridWith(table.root, 'Metric');
      expect(headerTooltip(grid, 'Test').startsWith('Value on the rows held out before the search'), true,
        headerTooltip(grid, 'Test'));
      const lines = Array.from(table.root.querySelectorAll('ul > li'), (li) => li.textContent ?? '');
      expect(lines.includes('Search: Tune by Accuracy, limit 1 min: 12 trials of 1 method, chosen by the ' +
        'one-standard-error rule'), true, lines.join(' | '));
      const halving = new MetricsTable({...summary!, search: {...searchRecord(), mode: 'best', budgetSeconds: undefined,
        maxTrials: 100, strategy: 'halving'}});
      expect((halving.root.textContent ?? '').includes('Search: Find best model by Accuracy, limit 100 trials, Focus ' +
        'on leaders: 1 method, chosen by the highest Accuracy − SE'), true, halving.root.textContent ?? '');
      const continued = new MetricsTable({...summary!, search: {...searchRecord(), maxTrials: 50}});
      expect((continued.root.textContent ?? '').includes('Search: Tune by Accuracy, limit 1 min + 50 trials: 12 ' +
        'trials of 1 method'), true, continued.root.textContent ?? '');
      expect(lines.includes('Test: 1170 rows held out before the search'), true, lines.join(' | '));
      const optimistic = new MetricsTable({...summary!, search: searchRecord()});
      expect((optimistic.root.textContent ?? '').includes('Validation after a search is optimistic'), true,
        optimistic.root.textContent ?? '');

      const accordion = modelAccordion(ForgeModelHandler.rowOf(model), model);
      await expectPaneText(accordion, 'Details', ['Search', 'Tune (Accuracy, limit 1 min, 12 trials)']);
    } finally {
      await deleteModel(id);
    }
  }, {timeout: TIMEOUT});
});

category('UI: Cards', () => {
  test('bestMarks: 3 decimals, ties, all when all tie, not one card; kindLine, settingsOf, settingsBlock', async () => {
    expectArray(ModelCards.bestMarks([0.9801, 0.9804, 0.97], 1), [true, true, false]);
    expectArray(ModelCards.bestMarks([0.03, 0.02, 0.0204], -1), [false, true, true]);
    expectArray(ModelCards.bestMarks([0.5], 1), [false]);
    expectArray(ModelCards.bestMarks([0.5, 0.5], 1), [true, true]);
    expectArray(ModelCards.bestMarks([0.5001, 0.4999], 1), [true, true]);
    expectArray(ModelCards.bestMarks([undefined, undefined], 1), [false, false]);
    expectArray(ModelCards.bestMarks([0.5, undefined], 1), [true, false]);

    const engine = xgboost();
    const defaults = defaultHyperparameters(engine);
    const candidate = (method: string): CandidateRecord => ({method, trials: 4, outcome: 'complete'});
    const find = (methods: string[]): SearchRecord =>
      ({...searchRecord(), mode: 'best', methods: methods.map(candidate)});
    expectArray([ModelCards.kindLine({settingsStatus: 'Default'}), ModelCards.kindLine({settingsStatus: 'Customized'}),
      ModelCards.kindLine({settingsStatus: 'Auto-tuned', search: searchRecord()}),
      ModelCards.kindLine({settingsStatus: 'Auto-tuned', search: find(['XGBoost'])}),
      ModelCards.kindLine({settingsStatus: 'Auto-tuned', search: find(['XGBoost', 'SVM'])})],
    ['Default settings', 'Customized settings', 'Auto-tuned settings', 'Auto-tuned settings', 'Auto-tuned model']);

    const properties = visibleHyperparametersOf(engine);
    const entries = (values: Record<string, unknown>) => properties
      .map((p) => `${DG.StringUtils.camelCaseToSentence(p.caption || p.name)} ${String(values[p.name])}`);
    expect(properties.some((p) => p.name === 'iterations'), true, 'XGBoost has no iterations');
    const shown = (values: typeof defaults) => ModelCards.settingsOf({engine, hyperparameters: values})
      .map((s) => `${s.caption} ${s.value}`);
    expectArray(shown(defaults), entries(defaults));
    expectArray(shown({...defaults, iterations: 30}), entries({...defaults, iterations: 30}));
    expectArray(ModelCards.settingsOf({engine, hyperparameters: defaults}).map((s) => s.description),
      properties.map((p) => p.description ?? ''));

    expectArray(ModelCards.settingsBlock('Default settings', []), []);
    const settings = (n: number) => Array.from({length: n}, (_, i) =>
      ({caption: `Setting ${i}`, value: `${i}`, description: i === 0 ? '' : `What ${i} does`}));
    /** The text of [block]'s kind line, then its one table's rows, each row's cells joined by spaces. */
    const blockOf = (block: HTMLElement[]) => [block[0].textContent ?? '',
      ...tableRows(block[1].querySelector('table')!).map((r) => r.join(' '))];
    expect(ModelCards.settingsBlock('Auto-tuned settings', settings(5))[1].querySelectorAll('table').length, 1);
    expectArray(blockOf(ModelCards.settingsBlock('Auto-tuned model', settings(1))),
      ['Auto-tuned model', 'Setting 0 0']);
    expectArray(blockOf(ModelCards.settingsBlock('Customized settings', settings(3))),
      ['Customized settings', 'Setting 0 0 Setting 1 1', 'Setting 2 2']);
    expectArray(blockOf(ModelCards.settingsBlock('Auto-tuned model', settings(4))),
      ['Auto-tuned model', 'Setting 0 0 Setting 1 1', 'Setting 2 2 Setting 3 3']);
    expectArray(blockOf(ModelCards.settingsBlock('Auto-tuned settings', settings(5))),
      ['Auto-tuned settings', 'Setting 0 0 Setting 1 1', 'Setting 2 2 Setting 3 3', 'Setting 4 4']);

    const block = ModelCards.settingsBlock('Default settings', settings(2));
    const host = ui.divV(block);
    document.body.append(host);
    try {
      const cell = (k: number) => host.querySelectorAll('td')[k].firstElementChild as HTMLElement;
      expectArray([await enterTooltip(cell(2)), await enterTooltip(cell(3)), await enterTooltip(cell(0)),
        await enterTooltip(block[0])], ['What 1 does', 'What 1 does', '', '']);
    } finally {
      host.remove();
    }
  });

  test('the current card in Trained models: the method, the settings status and line, Compare by\'s validation ' +
    'value and the status; a Table change empties it', async () => {
    const iris = await openIris();
    const demog = grok.shell.addTable(await openDemog());
    const view = await openTrainView(iris);
    const cards = view.modelCards;
    try {
      const group = cards.group.root;
      expect(group.previousElementSibling === view.trainButton.closest('.forge-train-row'), true,
        'Trained models is not right after the Train row');
      expect(group.querySelector('.forge-group-header label')?.textContent, 'Trained models');
      expect(cards.current === view.runs[0], true, 'The current card is not the newest run');
      const first = currentCard(view);
      expect(isShown(group) && isShown(first), true, 'The current card is hidden');
      const run = view.runs[0];
      const metrics = run.metrics!;
      const ids = metricIdsOf({metrics});
      expect(cards.compareBy.value, 'Balanced accuracy');
      expectArray(cards.compareBy.items, ids.map((id) => METRIC_LABELS[id]));
      expectArray([first.dataset.forgeRun ?? '', cardText(first, 'forge-card-title')], ['1', 'XGBoost']);
      expect(run.settingsStatus, 'Default');
      /** [card]'s kind line, then its settings as `<name> <value>`, row by row. */
      const shownSettings = (card: HTMLElement) => [cardText(card, 'forge-card-settings'),
        ...tableRows(card.querySelector<HTMLElement>('.forge-card-settings-table')!)
          .flatMap((c) => c.length > 2 ? [`${c[0]} ${c[1]}`, `${c[2]} ${c[3]}`] : [`${c[0]} ${c[1]}`])];
      expectArray(shownSettings(first),
        ['Default settings', ...ModelCards.settingsOf(run).map((s) => `${s.caption} ${s.value}`)]);
      const [train, validation] = [metrics.train.balancedAccuracy!, metrics.validation.balancedAccuracy!];
      expectArray([cardText(first, 'forge-card-metric'), cardText(first, 'forge-card-value')],
        ['Balanced accuracy', validation.toFixed(3)]);
      const {isOverfitting} = overfittingOf('balancedAccuracy', train, validation);
      const statusPhrase = isOverfitting ? 'Overfitting' : 'No overfitting';
      const status = first.querySelector<HTMLElement>('.forge-card-status')!;
      expect(`${status.textContent} ${status.classList.contains(isOverfitting ? 'forge-card-warn' : 'forge-card-ok')}`,
        `${isOverfitting ? '!' : ''}${statusPhrase} true`);
      expect(status.querySelector('.forge-card-mark') !== null, isOverfitting, 'The Overfitting mark is not a "!"');
      const rect = (cls: string) => first.querySelector(`.${cls}`)!.getBoundingClientRect();
      const [title, label, icon, statusRect] = ['forge-card-title', 'forge-card-metric', 'forge-card-pin',
        'forge-card-status'].map(rect);
      expect(Math.abs(label.left - title.left) <= 1 && Math.abs(statusRect.right - icon.right) <= 1, true,
        `metric ${label.left} vs title ${title.left}, status ${statusRect.right} vs icon ${icon.right}`);
      const settingsRight = Math.max(...Array.from(first.querySelectorAll('.forge-card-settings, ' +
        '.forge-card-settings-table table'), (e) => e.getBoundingClientRect().right));
      expect(settingsRight <= icon.left + 1, true, `The settings end at ${settingsRight}, the icon at ${icon.left}`);
      const [kindLeft, tableLeft] = ['forge-card-kind', 'forge-card-settings-table'].map((c) => rect(c).left);
      expect(Math.abs(kindLeft - tableLeft) <= 1 && kindLeft > title.left + 5, true,
        `The kind starts at ${kindLeft}, the table at ${tableLeft}, the title at ${title.left}`);
      const labelTip = await enterTooltip(first.querySelector<HTMLElement>('.forge-card-metric')!);
      const valueTip = await enterTooltip(first.querySelector<HTMLElement>('.forge-card-value')!);
      expect(`${labelTip === METRIC_DESCRIPTIONS.balancedAccuracy} ${valueTip.startsWith('MetricTrainValidation')}`,
        'true true', `${labelTip} | ${valueTip}`);
      expect(isBest(first), false, 'Best on the only card');
      expectArray(visibleIcons(first), ['+']);
      expect(pinIcon(first).getAttribute('aria-label'), 'Pin to compare with other models');

      const rows = tableRows(ModelCards.metricsTooltip(run, 'balancedAccuracy'));
      expectArray(rows[0], ['Metric', 'Train', 'Validation']);
      expectArray(rows.slice(1).map((row) => row[0]), ids.map((id) => METRIC_LABELS[id]));
      const selected =
        ModelCards.metricsTooltip(run, 'balancedAccuracy').querySelectorAll('tr.forge-card-tip-selected');
      expectArray(Array.from(selected, (row) => row.firstElementChild?.textContent ?? ''), ['Balanced accuracy']);
      const withTest = {...run, metrics: {...metrics, test: metrics.validation}};
      expectArray(tableRows(ModelCards.metricsTooltip(withTest, 'balancedAccuracy'))[0],
        ['Metric', 'Train', 'Validation']);
      const statusTip = ModelCards.statusTooltip(run, 'balancedAccuracy')!;
      const tipLines = Array.from(statusTip.children, (child) => child.textContent ?? '');
      expectArray([tipLines.length, tipLines[0], tipLines[2]], [3, 'Balanced accuracy:',
        isOverfitting ? 'May perform worse on new data.' : 'Should perform about as well on new data.']);
      const tipRows = tableRows(statusTip.querySelector('table')!);
      expectArray(tipRows.map((row) => row[0]), ['Train', 'Validation', 'Gap']);
      const [shownTrain, shownValidation, shownGap] = tipRows.map((row) => Number(row[1]));
      expect(shownGap.toFixed(3), (shownTrain - shownValidation).toFixed(3), 'The Gap is not the shown difference');
      expect(`${shownTrain.toFixed(3)} ${shownValidation.toFixed(3)}`,
        `${train.toFixed(3)} ${validation.toFixed(3)}`);
      /** The status tooltip of [run] with [id]'s train and validation values replaced: its table rows and sentence. */
      const tipOf = (id: MetricId, t: number, v: number) => {
        const tip = ModelCards.statusTooltip({...run, metrics: {...metrics, train: {...metrics.train, [id]: t},
          validation: {...metrics.validation, [id]: v}}}, id)!;
        return [...tableRows(tip.querySelector('table')!).map((row) => row.join(' ')),
          tip.lastElementChild?.textContent];
      };
      expectArray(tipOf('balancedAccuracy', 0.9734, 0.9466),
        ['Train 0.973', 'Validation 0.947', 'Gap 0.026', 'Should perform about as well on new data.']);
      expectArray(tipOf('logLoss', 0.1, 0.151),
        ['Train 0.100', 'Validation 0.151', 'Ratio 1.51', 'May perform worse on new data.']);
      expectArray(tipOf('logLoss', 0, 0.2),
        ['Train 0.000', 'Validation 0.200', 'Ratio -', 'May perform worse on new data.']);
      const form = view.root.querySelector<HTMLElement>('.forge-train-form')!;
      expect(form.scrollWidth <= form.clientWidth, true, `${form.scrollWidth} > ${form.clientWidth}`);

      setHyperparameter(view, 'iterations', 30);
      await awaitTrained(view);
      expect(cards.current === view.runs[0] && view.runs.length === 2, true, 'The current card is not the new run');
      expect(first.isConnected, false, 'The previous current card stays');
      const second = currentCard(view);
      expect(second.dataset.forgeRun, '2');
      expect(cards.root.querySelectorAll('.forge-card').length, 1);
      const iterations30 = `${view.hyperparameterInputs.get('iterations')!.caption} 30`;
      const secondSettings = shownSettings(second);
      expectArray(secondSettings,
        ['Customized settings', ...ModelCards.settingsOf(view.runs[0]).map((s) => `${s.caption} ${s.value}`)]);
      expect(secondSettings.some((s) => s.endsWith(iterations30)), true, secondSettings.join(' | '));

      view.tableInput.value = demog;
      expect(cards.current === undefined && cards.pinned.length === 0, true, 'A Table change keeps the cards');
      expect(isShown(cards.group.root), false, 'Trained models stays shown after a Table change');
    } finally {
      await closeTrainView(view, [iris, demog]);
      grok.shell.closeTable(demog);
    }
  }, {timeout: TIMEOUT});

  /** A click while the form is locked is restoreRun's own guard (UI: Train, UI: Search): not repeated here. */
  test('pin and unpin: a pinned card follows Compare by; the trash unpins; a click on a pinned card restores it',
    async () => {
      const iris = await openIris();
      const view = await openTrainView(iris);
      const cards = view.modelCards;
      try {
        const first = currentCard(view);
        pinIcon(first).click();
        expectArray(cards.pinned.map((run) => run.number), [1]);
        const [pinned] = pinnedCards(view);
        expect(first.nextElementSibling === pinned, true, 'The pinned card is not under the current one');
        expect(pinned.dataset.forgeRun, '1');
        expectArray(visibleIcons(pinned), ['trash']);
        expect(unpinIcon(pinned).getAttribute('aria-label'), 'Remove from the comparison');
        expectArray(visibleIcons(first), []);
        expect(first.querySelector('.forge-card-unpin'), null, 'The current card has a trash');

        unpinIcon(pinned).click();
        expect(cards.pinned.length === 0 && !pinned.isConnected, true, 'The pinned card\'s trash does not unpin');
        expectArray(visibleIcons(first), ['+']);

        pinIcon(first).click();
        setHyperparameter(view, 'iterations', 30);
        await awaitTrained(view);
        const second = currentCard(view);
        expectArray(visibleIcons(second), ['+']);
        pinIcon(second).click();
        expectArray(cardNumbers(pinnedCards(view)), sortedNumbers(cards.pinned, 'balancedAccuracy'));
        expect(second.nextElementSibling === pinnedCards(view)[0], true, 'A pinned card is above the current one');
        const pinnedOf = (number: string) => pinnedCards(view).find((card) => card.dataset.forgeRun === number)!;

        pinnedOf('1').querySelector<HTMLElement>('.forge-card-header')!.click();
        await awaitCheck(() => view.hyperparameterInputs.get('iterations')?.value === 20,
          'A click on a pinned card does not set its inputs', WAIT_MS);
        await awaitTrained(view);
        expect(view.runs.length, 3);
        const third = currentCard(view);
        expect(third.dataset.forgeRun, '3');
        expectArray(visibleIcons(third), []);

        const other = pinnedOf('2');
        unpinIcon(other).click();
        await delay(CHECK_DELAY_MS * 3);
        expectArray(cardNumbers(pinnedCards(view)), ['1']);
        const iterations = view.hyperparameterInputs.get('iterations')?.value;
        expect(!other.isConnected && !view.isTraining && view.runs.length === 3 && iterations === 20, true,
          'The trash of a pinned card of other inputs does not remove it only');
        expectArray(visibleIcons(third), []);

        unpinIcon(pinnedCards(view)[0]).click();
        expect(cards.pinned.length, 0);
        expectArray(visibleIcons(third), ['+']);
      } finally {
        await closeTrainView(view, [iris]);
      }
    }, {timeout: TIMEOUT});

  test('Compare by sorts the pinned cards and moves Best; the pick stays after a training', async () => {
    const iris = await openIris();
    const view = await openTrainView(iris);
    const cards = view.modelCards;
    /** The pinned cards are in [metric]'s order, every shown card has its value and Best is where bestMarks puts it. */
    const expectCompared = (metric: MetricId) => {
      expectArray(cardNumbers(pinnedCards(view)), sortedNumbers(cards.pinned, metric));
      const shown = [currentCard(view), ...pinnedCards(view)];
      const runs = shown.map((card) => view.runs.find((run) => `${run.number}` === card.dataset.forgeRun)!);
      expectArray(shown.map((card) => `${cardText(card, 'forge-card-metric')} ${cardText(card, 'forge-card-value')}`),
        runs.map((run) => `${METRIC_LABELS[metric]} ${validationOf(run, metric)?.toFixed(3) ?? '-'}`));
      expectArray(shown.map(isBest), ModelCards.bestMarks(runs.map((run) => validationOf(run, metric)),
        directionOf(metric)));
    };
    try {
      pinIcon(currentCard(view)).click();
      setHyperparameter(view, 'iterations', 1);
      await awaitTrained(view);
      pinIcon(currentCard(view)).click();
      setHyperparameter(view, 'iterations', 30);
      await awaitTrained(view);
      expect(`${view.runs.length} ${cards.pinned.length}`, '3 2');
      expect(cards.metricId, 'balancedAccuracy');
      expectCompared('balancedAccuracy');

      cards.compareBy.value = 'Log loss';
      await waitUntil(() => cards.metricId === 'logLoss', `Compare by is ${cards.metricId}`, WAIT_MS);
      expectCompared('logLoss');

      setHyperparameter(view, 'iterations', 40);
      await awaitTrained(view);
      expect(`${view.runs.length} ${cards.compareBy.value} ${cards.metricId}`, '4 Log loss logLoss',
        'The pick did not stay after a training');
      expectCompared('logLoss');
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: TIMEOUT});

  test('cards hide while the target, the validation or the seed differ and come back at once', async () => {
    const iris = await openIris();
    const view = await openTrainView(iris);
    const {seed, scheme} = view.validationInputs;
    try {
      const first = currentCard(view);
      pinIcon(first).click();
      const [pinned] = pinnedCards(view);
      const species = view.targetInput.value!;
      view.featuresInput.value = columnsOf(iris, ['Sepal.Length', 'Sepal.Width', 'Petal.Width']);
      expect(isShown(first) && isShown(pinned), true, 'A Features change hides a card');
      view.targetInput.value = iris.getCol('Petal.Length');
      expect(isShown(first) || isShown(pinned), false, 'A card of another target is shown');
      expect(isShown(view.modelCards.group.root), false, 'Trained models is shown without a card of this data');
      await awaitTrained(view);
      const regression = currentCard(view);
      expect(isShown(regression) && !isShown(pinned), true, 'The regression card is hidden or the pinned one shown');
      expect(view.modelCards.metricId, 'r2');
      expectArray(view.modelCards.compareBy.items.slice(0, 4), ['MSE', 'RMSE', 'MAE', 'R2']);

      view.targetInput.value = species;
      expect(isShown(pinned) && !isShown(regression), true, 'The pinned card is not back at once');
      const shownSeed = seed.value;
      seed.value = (shownSeed ?? 0) + 1;
      expect(isShown(pinned), false, 'Shown on another seed');
      seed.value = shownSeed;
      expect(isShown(pinned), true, 'Not back on its seed');
      const shownScheme = scheme.value;
      scheme.value = 'Hold-out';
      expect(isShown(pinned), false, 'Shown on another scheme');
      scheme.value = shownScheme;
      expect(isShown(pinned), true, 'Not back on its scheme');
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: TIMEOUT});
});

/** Chooses Auto in [view]'s Method and waits until the check has listed the methods. */
async function chooseAuto(view: TrainView): Promise<FindSettings> {
  view.methodInput.value = 'Auto';
  await awaitCheck(() => view.findSettings.checkedMethods.length > 0 && resultsState(view) === 'false/false',
    'Auto lists no checked method', TIMEOUT);
  return view.findSettings;
}

/** The Find best model settings [view] remembers, its pending save written first. */
function storedFind(view: TrainView): {[key: string]: any} {
  view.findSettings.flushSave();
  return TuneSettings.storedRecord(FIND_KEY);
}

/** The items of the Candidates of [find]: each listed method's name and check. */
const findItems = (find: FindSettings) =>
  Array.from(find.candidates.root.querySelectorAll<HTMLElement>('.ui-input-multi-choice-checks > div'),
    (item) => ({name: item.textContent ?? '', check: item.querySelector('input')!}));
const listedMethods = (find: FindSettings) => findItems(find).map((i) => i.name);

function findItem(find: FindSettings, name: string): {check: HTMLInputElement} {
  const found = findItems(find).find((i) => i.name === name);
  if (found === undefined)
    throw new Error(`No method ${name} in Candidates`);
  return found;
}

/** Toggles the method [name] of the Candidates of [find], as a click does. */
const clickMethod = (find: FindSettings, name: string) => findItem(find, name).check.click();

category('UI: Find best', () => {
  test('Auto is the first Method item; choosing it shows the Find block, keeps the model shown and runs nothing',
    async () => {
      const iris = await openIris();
      const view = await openTrainView(iris);
      try {
        expect(view.methodInput.items[0], 'Auto');
        const trained = `${view.runs.length} ${shownRun(view)} ${view.resultsPane.frames.length}`;
        const training = view.lastTraining;
        const find = await chooseAuto(view);
        expect(view.isAutoMethod, true);
        expect(`${isShown(view.tuneSettings.root)} ${isShown(find.root)} ${isTrainShown(view)}`, 'false true false',
          'Tune, Find, Train shown');
        expect([...view.hyperparameterInputs.values()].every((i) => i.root.offsetParent === null), true,
          'A setting of the method is shown under Auto');
        expect(await enterTooltip(view.methodInput.captionLabel),
          'Auto: Find tries every checked method and trains the best one');
        await delay(CHECK_DELAY_MS * 3);
        expect(`${view.runs.length} ${shownRun(view)} ${view.resultsPane.frames.length}`, trained,
          'Auto trained or hid the model\'s results');
        expect(view.lastTraining === training && !view.isTraining, true, 'Auto dropped or retrained the model');
        expect(resultsText(view).includes('Click Find'), false, resultsText(view));
        view.featuresInput.value = view.featuresInput.value.slice(1);
        await awaitCheck(() => resultsText(view).includes('Click Find'), 'No Click Find without a model', WAIT_MS);
        expect(shownRun(view), null, 'A model is shown after a Features change');
        expect(find.statusText.textContent, 'Not searched');
        expect(await enterTooltip(find.statusText), 'Click Find to search');
        const classifiers = COST_RANK.filter((name) => CLASSIFIERS.includes(name));
        expectArray(listedMethods(find), classifiers);
        expectArray(find.checkedMethods, classifiers.filter((name) => !['RealMLP', 'SVM'].includes(name)));
        find.limit.kind.value = TIME;
        find.limit.minutes.value = 5;
        expect(findItem(find, 'RealMLP').check.checked, true, 'RealMLP unchecked at 5 minutes');
        find.limit.minutes.value = 1;
        expect(findItem(find, 'RealMLP').check.checked, false, 'RealMLP checked at 1 minute');
        const runs = view.runs.length;
        clickMethod(find, 'XGBoost');
        expect(storedFind(view).methods?.XGBoost, false, 'The toggle is not remembered');
        await delay(CHECK_DELAY_MS * 3);
        expect(`${view.runs.length} ${view.isTraining}`, `${runs} false`, 'Auto trained');
        for (const name of find.checkedMethods)
          clickMethod(find, name);
        expect(await hoverText(view.findButton), 'Check at least one method');
        view.methodInput.value = 'XGBoost';
        expect(`${isShown(view.tuneSettings.root)} ${isShown(find.root)}`, 'true false', 'Tune, Find shown');
        await awaitTrained(view);
      } finally {
        await closeTrainView(view, [iris]);
        clearTestTuneStorage();
      }
    }, {timeout: TIMEOUT});

  test('Find remembers Strategy and the limit; an old stored time becomes whole minutes', async () => {
    const iris = await openIris();
    grok.userSettings.add(TuneSettings.storage, FIND_KEY, JSON.stringify({seconds: 30, trials: 70}));
    let view = await openTrainView(iris);
    try {
      let find = await chooseAuto(view);
      expect(`${find.strategy.value} ${find.limit.kind.value} ${find.limit.minutes.value} ${find.limit.trials.value}`,
        `Equal ${TIME} 1 70`);
      find.strategy.value = 'Leaders';
      find.limit.minutes.value = 7;
      stopAfterTrials(find.limit, 40);
      const stored = storedFind(view);
      expect(`${stored.strategy} ${stored.limit} ${stored.minutes} ${stored.trials} ${stored.seconds}`,
        `halving ${TRIALS} 7 40 undefined`);
      await closeTrainView(view, [iris]);
      view = await openTrainView(iris);
      find = await chooseAuto(view);
      expect(`${find.strategy.value} ${find.limit.kind.value} ${find.limit.minutes.value} ${find.limit.trials.value}`,
        `Leaders ${TRIALS} 7 40`);
      grok.userSettings.add(TuneSettings.storage, FIND_KEY, JSON.stringify({strategy: 'fast', seconds: 3600.5}));
      const {strategy, limit, minutes, trials} = FindSettings.load();
      expect(`${strategy} ${limit} ${minutes} ${trials}`, `equal ${TIME} 60 100`);
    } finally {
      await closeTrainView(view, [iris]);
      clearTestTuneStorage();
    }
  }, {timeout: TIMEOUT});

  test('Search: the cautious scores as a bar chart left of the grid over its frame, described, no title, line, ' +
    'coloring or selectors; splitters between the blocks; a halving round in the tooltips; a bar click trains its ' +
    'method',
  async () => {
    const iris = await openIris();
    const view = await openTrainView(iris);
    try {
      await chooseAuto(view);
      await view.startSearch('best', {metric: 'accuracy', maxTrials: 16, strategy: 'halving',
        methods: ['Naive Bayes', 'Decision Tree', 'k-Nearest Neighbors', 'XGBoost', 'Extra Trees', 'Random Forest']});
      const search = view.lastTraining!.result.search!;
      expect(`${search.strategy} ${search.maxTrials}`, 'halving 16');
      const naiveBayes = searchGrid(view).dataFrame.getCol('Method').toList().indexOf('Naive Bayes');
      expect(gridTooltip(searchGrid(view), searchGrid(view).cell('Trials', naiveBayes)), 'No settings to tune');
      expect(resultsText(view).includes('Search: Find best model by Accuracy, limit 16 trials, Focus on leaders: ' +
        `${search.leaders.length} method`), true, resultsText(view));
      const chart = view.searchPanel.chart!;
      const grid = searchGrid(view);
      await delay(500);
      const [chartBox, gridBox] = [chart.root, grid.root].map((e) => e.getBoundingClientRect());
      const share = chartBox.width / (chartBox.width + gridBox.width);
      expect(chartBox.right <= gridBox.left + 1 && Math.abs(share - 2 / 3) < 0.03 &&
        Math.abs(chartBox.top - gridBox.top) <= 1 && Math.abs(chartBox.height - gridBox.height) <= 1, true,
      `The chart at ${chartBox.left}-${chartBox.right} (${share} of the width), ${chartBox.height} px high at ` +
        `${chartBox.top}; the grid at ${gridBox.left}, ${gridBox.height} px at ${gridBox.top}`);
      const splitOf = (root: HTMLElement) => root.parentElement?.parentElement;
      expect(splitOf(chart.root) === splitOf(grid.root) && splitOf(chart.root)?.classList.contains('ui-split-h'), true,
        'The chart and the grid share no splitter');
      // Wrong predictions, Search and Runs in one vertical splitter.
      const split = view.searchPanel.root.closest('.ui-split-v');
      const wrong = Array.from(view.resultsPane.root.querySelectorAll<HTMLElement>('.forge-results-title'))
        .find((t) => t.textContent === 'Wrong predictions');
      expect(split !== null && split.contains(view.runsPanel.root) && wrong !== undefined && split.contains(wrong),
        true, 'Wrong predictions, Search and Runs are not in one splitter');
      const {props} = chart;
      expectArray([props.showValueSelector, props.showCategorySelector, props.showStackSelector, props.legendVisibility,
        props.showTitle, props.description, props.descriptionVisibilityMode, props.stackColumnName ?? '',
        props.colorColumnName ?? '', props.showSelectedRows, props.barSortType, props.barSortOrder,
        props.valueColumnName],
      [false, false, false, 'Never', false, 'Accuracy minus its standard error: a cautious estimate', 'Always', '', '',
        false, 'by value', 'desc', 'Accuracy − SE']);
      const title = view.searchPanel.root.querySelector<HTMLElement>('.forge-results-title span')!;
      expect(await enterTooltip(title), 'Each method\'s best settings, ranked by Accuracy − SE.');
      const frame = chart.dataFrame;
      expect(frame === grid.dataFrame, true, 'The chart and the grid show different frames');
      expect(frame.meta.formulaLines.items.length, 0);
      const leaders = search.leaders;
      expectSearchOrder(grid, leaders);
      const bars = Object.keys(chart.getWidgetStatus().hitAreas).filter((k) => k.startsWith('bar '));
      expectArray(bars.sort(), leaders.map((l) => `bar ${l.method}`).sort());
      expectChartOrder(chart, grid);
      const eliminated = (search.methods ?? []).filter((m) => m.eliminatedRound !== undefined);
      expect(eliminated.length > 0, true, 'No method was eliminated by the halving');
      const row = grid.dataFrame.getCol('Method').toList().indexOf(eliminated[0].method);
      const round = `Eliminated in round ${eliminated[0].eliminatedRound! + 1}`;
      expect(gridTooltip(grid, grid.cell('Trials', row)), round);
      const methodTooltip = gridTooltip(grid, grid.cell('Method', row));
      expect(methodTooltip.includes(round), true, methodTooltip);

      const other = search.leaders.find((l) => !l.chosen)!;
      const area = chart.getWidgetStatus().hitAreas[`bar ${other.method}`];
      expect(area !== undefined, true, `No bar of ${other.method}`);
      const canvas = chart.getInfo()['canvas'] as HTMLCanvasElement;
      const rect = canvas.getBoundingClientRect();
      const at = {clientX: rect.left + area.x + area.width / 2, clientY: rect.top + area.y + area.height / 2,
        bubbles: true};
      canvas.dispatchEvent(new MouseEvent('mousemove', at));
      canvas.dispatchEvent(new MouseEvent('click', at));
      await awaitCheck(() => view.methodInput.value === other.method, 'A bar click does not set its method', WAIT_MS);
      await awaitTrained(view);
      expect(view.lastTraining?.engine.name, other.method);
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: SEARCH_TIMEOUT});

  test('Find with a cancel trains the best so far and the form is locked meanwhile', async () => {
    const iris = await openIris();
    const view = await openTrainView(iris);
    try {
      const find = await chooseAuto(view);
      await cancelSearchOnceSettled(view, 'Finding best model', (maxTrials) => {
        const searching = view.startSearch('best', {metric: 'accuracy', maxTrials,
          methods: ['XGBoost', 'Decision Tree', 'Softmax']});
        expectArray(enabledInputs(view), []);
        expect(find.inputs.every((i) => !i.enabled) && find.metric.disabled &&
          findItems(find).every((i) => i.check.disabled), true, 'The Find block is not locked');
        return searching;
      });
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: SEARCH_TIMEOUT});

  test('the Find block looks like the Auto-tune block', async () => {
    const iris = await openIris();
    const view = await openTrainView(iris);
    try {
      // The platform sets the form's label width a while after an input shows (the block's Time among them).
      const labelWidths = () => Array.from(view.root.querySelectorAll<HTMLElement>('label.ui-input-label'))
        .filter((l) => l.offsetParent !== null).map((l) => l.style.minWidth).join(',');
      const settle = async () => {
        for (let last = '', i = 0; i < 20 && labelWidths() !== last; i++) {
          last = labelWidths();
          await delay(600);
        }
      };
      autoTune(view);
      await awaitCheck(() => !isDisabled(view.tuneButton), 'Tune stays disabled', WAIT_MS);
      await settle();
      const tune = searchBlockLines(view.tuneSettings, view.tuneButton);
      const find = await chooseAuto(view);
      await awaitCheck(() => !isDisabled(view.findButton), 'Find stays disabled', WAIT_MS);
      expect(find.isExpanded, false, 'Search options is expanded by default');
      find.setExpanded(true);
      expect(find.statusText.textContent, 'Not searched');
      await settle();
      expectArray(searchBlockLines(find, view.findButton), tune);
      const methodLabel = view.methodInput.captionLabel.getBoundingClientRect();
      const status = textRect(find.statusText);
      expect(status.top >= methodLabel.bottom - 1 &&
        (find.statusText.parentElement!.compareDocumentPosition(find.toggle) & Node.DOCUMENT_POSITION_FOLLOWING) !== 0,
      true, 'The status is not between Method and Search options');
      expectArray(Array.from(find.grid.querySelectorAll('.forge-tune-label'), (l) => l.textContent ?? ''),
        ['Metric', 'Strategy', 'Limit', 'Max time', 'Max trials', 'Candidates']);
      const settings = [...view.hyperparameterInputs.values()].map((i) => i.caption);
      expectArray(Array.from(view.tuneSettings.grid.querySelectorAll('.forge-tune-label'), (l) => l.textContent ?? ''),
        ['Metric', 'Limit', 'Max time', 'Max trials', ...settings]);
      expect((find.grid.compareDocumentPosition(view.findButton) & Node.DOCUMENT_POSITION_FOLLOWING) !== 0, true,
        'Find is not after Search options');
      expect(view.findButton.parentElement!.classList.contains('forge-tune-bar'), true, 'Find is not in a Tune bar');
      const [candidates, metric] = [find.candidatesLabel, gridLabelOf(find.grid, 'Metric')]
        .map((e) => e.getBoundingClientRect());
      // Labels beside the inputs: the first check level with Candidates; above them: under it.
      const isCondensed = find.root.closest('.ui-form-condensed') !== null;
      const firstCheck = findItems(find)[0].check.getBoundingClientRect();
      const middle = (r: DOMRect) => r.top + r.height / 2;
      const isPlaced = isCondensed ? firstCheck.top >= candidates.bottom - 1 :
        Math.abs(middle(firstCheck) - middle(candidates)) < 4;
      expect(Math.abs(candidates.right - metric.right) <= 1 && isPlaced &&
        textStyleOf(find.candidatesLabel) === textStyleOf(gridLabelOf(find.grid, 'Metric')), true,
      `Condensed ${isCondensed}: Candidates ends at ${candidates.right}, Metric at ${metric.right}; the first check ` +
        `at ${firstCheck.top}-${firstCheck.bottom}, the label at ${candidates.top}-${candidates.bottom}`);
    } finally {
      await closeTrainView(view, [iris]);
      clearTestTuneStorage();
    }
  }, {timeout: TIMEOUT});

  test('Predict probability under Auto follows the checked methods; SVM is unchecked by default at a quick limit on ' +
    'a classifier',
  async () => {
    const iris = await twoSpeciesIris();
    const view = await openTrainView(iris);
    const {predictProbability} = view.outputInputs;
    try {
      const find = await chooseAuto(view);
      expect(`${findItem(find, 'SVM').check.checked} ${findItem(find, 'XGBoost').check.checked}`, 'false true',
        'SVM, XGBoost checked');
      expect(predictProbability.enabled, true, 'Predict probability is disabled under Auto');
      predictProbability.value = true;
      find.limit.kind.value = TIME;
      find.limit.minutes.value = 5;
      expect(findItem(find, 'SVM').check.checked, true, 'SVM unchecked at 5 minutes');
      find.limit.minutes.value = 1;
      expect(findItem(find, 'SVM').check.checked, false, 'SVM checked at 1 minute');
      pick(find.metric, 'auc');
      predictProbability.value = false;
      expect(listedMethods(find).includes('SVM'), true, 'SVM is not listed under AUC-ROC without Predict probability');
      const tooltip = await enterTooltip(find.candidatesLabel);
      expect(tooltip.includes('Not listed: SVM'), false, tooltip);
      for (const name of find.checkedMethods)
        clickMethod(find, name);
      expect(predictProbability.enabled, true, 'Unchecking every method disabled Predict probability');
      expect(await hoverText(view.findButton), 'Check at least one method');
    } finally {
      await closeTrainView(view, [iris]);
      clearTestTuneStorage();
    }
  }, {timeout: TIMEOUT});

  test('Find leaves SVM unchecked by default above 1000 rows, at any limit and task; the Candidates tooltip says why',
    async () => {
      const table = twoClassTable(1001);
      const view = await openTrainView(table);
      try {
        const find = await chooseAuto(view);
        find.limit.kind.value = TIME;
        find.limit.minutes.value = 5;
        expect(findItem(find, 'SVM').check.checked, false, 'SVM checked on 1001 rows');
        find.setRowCount(1000);
        expect(findItem(find, 'SVM').check.checked, true, 'SVM unchecked on 1000 rows');
        find.setRowCount(1001);
        expect(findItem(find, 'SVM').check.checked, false, 'SVM checked on 1001 rows after 1000');
        const tooltip = await enterTooltip(find.candidatesLabel);
        expect(tooltip.includes('SVM: unchecked by default on more than 1000 rows'), true, tooltip);
        view.targetInput.value = table.getCol('a');
        view.featuresInput.value = [table.getCol('b')];
        await waitUntil(() => listedMethods(find).includes('Linear Regression'),
          () => `Candidates not rebuilt for a regression: ${listedMethods(find).join(', ')}`, TIMEOUT);
        expect(findItem(find, 'SVM').check.checked, false, 'SVM checked on a regression of 1001 rows');
        clickMethod(find, 'SVM');
        expect(findItem(find, 'SVM').check.checked, true, 'A click does not check SVM');
      } finally {
        await closeTrainView(view, [table]);
        clearTestTuneStorage();
      }
    }, {timeout: TIMEOUT});
});

/** The two models of `UI: Gallery`: one saved with a file and a tag, one bare row; [name] starts both names. */
let galleryModels: {saved: string; bare: string; name: string; tag: string} | undefined;

const menuLabels = () => Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label'),
  (label) => label.textContent?.trim() ?? '');

function closeMenus(): void {
  for (const popup of Array.from(document.querySelectorAll('.d4-menu-popup')))
    popup.remove();
}

const clickItem = (item: HTMLElement | undefined, init: MouseEventInit = {}) =>
  item?.dispatchEvent(new MouseEvent('click', {bubbles: true, cancelable: true, ...init}));

/** Right-clicks [item] and returns the labels of the menu it shows; the menu is closed. */
async function itemMenu(item: HTMLElement | undefined): Promise<string[]> {
  closeMenus();
  item?.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, cancelable: true, clientX: 20, clientY: 20}));
  await waitUntil(() => menuLabels().length > 0, 'The item shows no menu', WAIT_MS);
  const labels = menuLabels();
  closeMenus();
  return labels;
}

/** A gallery over [frame] (a `name` column) with plain cards and the hint `.forge-test-empty`. */
function testGallery(frame: DG.DataFrame): Gallery {
  const nameOf = (i: number) => ui.divText(frame.getCol('name').get(i));
  return new Gallery({key: 'test', noun: 'items', placeholder: 'Search', grid: DG.Viewer.grid(frame),
    searchColumns: ['name'], sortFields: [{caption: 'Name', column: 'name'}],
    defaultSort: {column: 'name', desc: false}, card: nameOf, listItem: nameOf, tooltip: () => null, menu: () => {},
    emptyHint: ui.divText('Nothing yet', 'forge-test-empty')});
}

category('UI: Gallery', () => {
  before(async () => {
    useTestTuneStorage();
    useTestGalleryStorage();
    const stamp = Date.now();
    const name = `forge-test-gallery-${stamp}`;
    const tag = `forge-test-tag-${stamp}`;
    const saved = (await saveIrisModel(await openIris(), {name: `${name}-saved`, tags: tag})).id;
    galleryModels = {saved, bare: await insertModelRow({name: `${name}-bare`}), name, tag};
  });

  after(async () => {
    clearTestTuneStorage();
    clearTestGalleryStorage();
    if (galleryModels !== undefined) {
      await deleteModel(galleryModels.saved);
      await forgeDb.models.delete(galleryModels.bare);
    }
  });

  test('the gallery shows cards, models newest first', async () => {
    const view = await ForgeApp.create();
    grok.shell.addView(view);
    try {
      expect(view.modelsGallery.mode, 'cards');
      expect(JSON.stringify(view.modelsGallery.sort), '{"column":"created_on","desc":true}');
      const created = view.models.getCol('created_on');
      const times = view.modelsGallery.shownRows().map((i) => dayjs(created.get(i)).valueOf());
      expect(times.every((t, k) => k === 0 || times[k - 1] >= t), true, 'The cards are not newest first');
      expect(view.modelsGallery.renderCount, 1, 'A new view renders its cards more than once');
    } finally {
      view.close();
      useTestGalleryStorage();
    }
  });

  test('the gallery remembers its mode and sort', async () => {
    const view = await ForgeApp.create();
    try {
      view.modelsGallery.setMode('list');
      view.modelsGallery.setSort('name', false);
    } finally {
      view.close();
    }
    const again = await ForgeApp.create();
    try {
      expect(`${again.modelsGallery.mode} ${JSON.stringify(again.modelsGallery.sort)}`,
        'list {"column":"name","desc":false}');
    } finally {
      again.close();
      useTestGalleryStorage();
    }
  });

  test('the search filters by name and #tag and the count follows it', async () => {
    const {saved, bare, name, tag} = savedFixture(galleryModels);
    const view = await ForgeApp.create();
    grok.shell.addView(view);
    const gallery = view.modelsGallery;
    const count = () => gallery.root.querySelector('.grok-items-view-counts')?.textContent ?? '';
    const shownIds = () => gallery.shownRows().map((i) => view.models.get('id', i));
    const total = view.models.rowCount;
    try {
      expect(count(), `${total}`);
      const bareRow = view.models.getCol('id').toList().indexOf(bare);
      const bareItem = gallery.itemOf(bareRow);
      // A search value change filters and renders at once.
      gallery.search.value = `${name}-BARE`;
      expect(shownIds().join(), bare, 'The search');
      expect(count(), `1 / ${total}`);
      expect(gallery.itemOf(bareRow) === bareItem, true, 'A search built the item again');
      gallery.search.value = `#${tag}`;
      expect(shownIds().join(), saved, 'The tag search');
      gallery.search.value = 'forge-test-no-such-model';
      expect(gallery.root.querySelector('.grok-items-view-empty') !== null, true, 'No empty-set icon without a match');
      gallery.search.value = '';
      expect(count(), `${total}`);
    } finally {
      view.close();
    }
  });

  test('a click shows the model, Ctrl+click compares, the menu has the catalog commands', async () => {
    const {saved, bare} = savedFixture(galleryModels);
    const view = await ForgeApp.create();
    grok.shell.addView(view);
    const gallery = view.modelsGallery;
    const ids = view.models.getCol('id').toList();
    const [a, b] = [ids.indexOf(saved), ids.indexOf(bare)];
    try {
      const marked = () => [a, b].filter((i) => gallery.itemOf(i)?.classList.contains('d4-current'));
      clickItem(gallery.itemOf(a));
      await waitUntil(() => shows(saved), 'A clicked card is not the current object', WAIT_MS);
      expectArray(marked(), [a]);
      clickItem(gallery.itemOf(b), {ctrlKey: true});
      await waitUntil(() => comparedIds() === [saved, bare].sort().join(), 'Ctrl+click compares nothing', WAIT_MS);
      expectArray(await itemMenu(gallery.itemOf(a)),
        ['Apply...', 'Edit model...', 'Share...', 'Download', 'Delete model', 'Compare']);
      clickItem(gallery.itemOf(a), {ctrlKey: true});
      expectArray(marked(), [b]);
      expect(view.models.currentRowIdx, b, 'The card left picked is not the current row');
      await waitUntil(() => shows(bare), 'The card left picked is not the current object', WAIT_MS);
      clickItem(gallery.itemOf(b));
      expectArray(await itemMenu(gallery.itemOf(b)), ['Apply...', 'Edit model...', 'Share...', 'Delete model']);
    } finally {
      closeMenus();
      view.close();
    }
  }, {timeout: TIMEOUT});

  test('the list shows each model as a link that makes it current', async () => {
    const {saved, name} = savedFixture(galleryModels);
    const view = await ForgeApp.create();
    grok.shell.addView(view);
    const gallery = view.modelsGallery;
    try {
      gallery.setMode('list');
      const i = view.models.getCol('id').toList().indexOf(saved);
      const line = gallery.itemOf(i);
      expect(line?.classList.contains('d4-link-label'), true, 'The line is not the platform\'s link');
      expect(`${line?.textContent} ${line?.children.length}`, `${name}-saved 0`, 'The line has more than the name');
      clickItem(line);
      await waitUntil(() => shows(saved), 'A clicked link is not the current object', WAIT_MS);
      expect(line?.classList.contains('d4-current'), true, 'The current line is not marked');
    } finally {
      view.close();
      useTestGalleryStorage();
    }
  });

  test('the empty states: the hint of an empty frame, the empty-set icon, at most 500 items', async () => {
    const empty = testGallery(DG.DataFrame.fromColumns([DG.Column.fromStrings('name', [])]));
    expect(empty.root.querySelector('.forge-test-empty') !== null, true, 'An empty frame shows no hint');
    empty.detach();
    const many = testGallery(DG.DataFrame.fromColumns([DG.Column.fromStrings('name',
      Array.from({length: MAX_ITEMS + 1}, (_, i) => `item ${i}`))]));
    expect(many.shownRows().length, MAX_ITEMS);
    expect((many.root.textContent ?? '').includes(`Showing ${MAX_ITEMS} of ${MAX_ITEMS + 1}.`), true,
      'No line about the items not shown');
    many.detach();
    const view = await ForgeApp.create();
    try {
      expect(view.emptyHint.textContent?.includes('No models yet.'), true, view.emptyHint.textContent ?? '');
      expect(view.emptyHint.querySelector('a')?.textContent, 'Train a model...');
    } finally {
      view.close();
    }
  });

  test('launchTraining opens the Train view with the method; a method that cannot learn says why', async () => {
    const iris = await openIris();
    const tableView = grok.shell.addTableView(iris);
    const trainViews = (): TrainView[] => {
      const views: unknown[] = Array.from(grok.shell.views);
      return views.filter((v): v is TrainView => v instanceof TrainView);
    };
    try {
      // launchTraining resolves once the view is added.
      await launchTraining('XGBoost');
      expect(trainViews().length, 1, 'launchTraining opens no Train view');
      const [xgboost] = trainViews();
      expect(xgboost.methodInput.value, 'XGBoost');
      await closeTrainView(xgboost, [iris]);

      grok.shell.v = tableView;
      await launchTraining('Linear Regression');
      expect(trainViews().length, 1, 'launchTraining opens no Train view');
      const [linear] = trainViews();
      try {
        const mismatch = 'Linear Regression cannot learn from this selection';
        await waitUntil(() => balloonWith(mismatch) !== null, 'No balloon for the mismatch', WAIT_MS);
        const balloon = balloonWith(mismatch)?.textContent ?? '';
        expect(balloon.includes('is chosen for now.'), true, balloon);
        await waitUntil(() => linear.methodInput.value !== 'Linear Regression' && linear.methodInput.value !== null,
          'Method keeps Linear Regression', WAIT_MS);
      } finally {
        await closeTrainView(linear, [iris]);
      }
    } finally {
      tableView.close();
    }
  }, {timeout: TIMEOUT});
});

// After the categories, so their order stays: the Train views of these keep their tune settings in the tests' storage
// ('UI: Table panel' uses it in its own hooks).
for (const name of ['UI', 'UI: Train', 'UI: Preparation', 'UI: Results', 'UI: Save model', 'UI: Runs', 'UI: Search',
  'UI: Cards', 'UI: Find best']) {
  category(name, () => {
    before(async () => {
      useTestTuneStorage();
      useTestGalleryStorage();
    });
    after(async () => {
      clearTestTuneStorage();
      clearTestGalleryStorage();
    });
  });
}

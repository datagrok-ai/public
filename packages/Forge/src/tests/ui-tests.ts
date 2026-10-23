/* eslint-disable max-lines */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {after, awaitCheck, before, category, delay, expect, expectArray, expectFloat, test}
  from '@datagrok-libraries/test/src/test';
import {recordApplication} from '../apply/application-store';
import {applyAndRecord, DEFAULT_BATCH_SIZE, loadModel} from '../apply/apply-model';
import {exactMapping} from '../apply/column-matching';
import {updateModelInfo} from '../catalog/model-edit';
import {APP_NAME, MENU_PATH, MODEL_TYPE, PANEL_TABLE_MODELS} from '../constants';
import {EngineRegistry} from '../engines/engine-registry';
import {errorMessage, ForgeError} from '../forge-error';
import {forgeDb, ModelStorageMode} from '../generated/db';
import {METRIC_DESCRIPTIONS} from '../metrics/metrics';
import {preparationOptionsOf} from '../preparation/preparation-options';
import {datasetRefOf, storedDatasetRef} from '../storage/dataset-ref';
import {deleteModel, modelsChanged} from '../storage/model-store';
import {metricsRecordOf, trainPredictionOf} from '../training/train-model';
import {applyModelDialog, ApplyDialogOptions, isOpen, modelLabels} from '../ui/apply-model-dialog';
import {ButtonGate} from '../ui/button-gate';
import {gridTooltip} from '../ui/data-grid';
import {editModelDialog} from '../ui/edit-model-dialog';
import {EvaluationView, runEvaluation} from '../ui/evaluation-view';
import {ForgeApp} from '../ui/forge-app';
import {CATALOG_ACTIONS, MODEL_ACTIONS, openShareDialog} from '../ui/model-actions';
import {hasFormsViewer, ModelComparison} from '../ui/model-comparison';
import {ForgeModelHandler} from '../ui/model-handler';
import {MetricsTable} from '../ui/metrics-table';
import {modelAccordion, refreshSharing} from '../ui/model-panes';
import {PanelKey, resultsLayout} from '../ui/results-layout';
import {saveModelDialog} from '../ui/save-model-dialog';
import {tableModelsPane} from '../ui/table-models-panel';
import {tagsOfInput} from '../ui/tags-input';
import {CHECK_DELAY_MS, SHOWN_RUN, TrainView} from '../ui/train-view';
import {balloonWith, cancelOnceShown, catalogState, CLASSIFIERS, columnsOf, deleteRunsOf, demoDemog, expectNoFrames,
  expectReleased, expectGridsFillPanes, framesSharing, insertModelRow, IRIS, MEASUREMENTS, openDemoTable, openIris,
  openIrisFromFile, readUntil, refreshThenChange, REGRESSORS, runsOf, savedFixture, saveIrisModel, saveKeptGapsModel,
  saveSkipUniqueModel, saveTestModel, taskBarEntry, twoSpeciesIris, unbalancedIris, valuesOf, WINE,
  XGBOOST_FIELDS} from './test-data';

const WAIT_MS = 5000;
const TIMEOUT = 60000;
/** The methods Method offers on iris, sorted: for Species, and for Petal.Length. */
const IRIS_CLASSIFIERS = [...CLASSIFIERS].sort();
const IRIS_REGRESSORS = [...REGRESSORS].sort();

/** Opens the Train view on [table] and waits for the training it starts by itself. */
async function openTrainView(table: DG.DataFrame): Promise<TrainView> {
  const view = await TrainView.create(table);
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

/** The tooltip the mouse entering the enabled [element] shows; '' if none. */
async function enterTooltip(element: HTMLElement): Promise<string> {
  ui.tooltip.hide();
  const rect = element.getBoundingClientRect();
  element.dispatchEvent(new MouseEvent('mouseenter', {clientX: rect.left + 2, clientY: rect.top + 2}));
  // The platform shows a tooltip 250 ms after the mouse enters.
  await delay(600);
  const text = ui.tooltip.isVisible ? ui.tooltip.root.textContent ?? '' : '';
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

category('UI', () => {
  test('Forge app opens', async () => {
    const view = await ForgeApp.create();
    expect(view.models.currentRowIdx, -1);
    expect(`${view.root.classList.contains('forge-app')} ${catalogState(view.root)}`, 'true 0/false', 'Class, load');
    expect(isDisabled(view.deleteIcon), true, 'Delete is enabled without a current row');
    expect(isDisabled(view.applyIcon), true, 'Apply is enabled without a current row');
    grok.shell.addView(view);
    let measured = '';
    try {
      expect(view.name, APP_NAME);
      for (const caption of ['Methods', 'Models'])
        expect((view.root.textContent ?? '').includes(caption), true, `${caption} is missing`);
      const methods = view.methodsGrid;
      expectArray(methods.dataFrame.columns.names(), ['Method', 'Package', 'Method type', 'Roles', 'Hyperparameters']);
      expect(methods.dataFrame.rowCount, EngineRegistry.discover().length);
      const xgboost = methods.dataFrame.getCol('Method').toList().indexOf('XGBoost');
      expect(xgboost >= 0 && methods.props.allowEdit === false, true, 'XGBoost is not listed, or Methods is editable');
      expect(headerTooltip(methods, 'Roles'), 'What the method can do.');
      expect(gridTooltip(methods, methods.cell('Method type', xgboost)), 'A package function.');
      const roles = gridTooltip(methods, methods.cell('Roles', xgboost));
      expect(roles.includes('train: Trains a model'), true, `The Roles tooltip: ${roles}`);
      measured = expectRowsFit('Methods', methods);
      expect(view.root.contains(view.deleteIcon), true, 'The delete icon is missing');
      expect(view.root.contains(view.applyIcon), true, 'The apply icon is missing');
      const header = view.root.querySelector('.forge-pane-header');
      const grid = header?.parentElement?.querySelector('.d4-grid');
      if (!(grid instanceof HTMLElement) || !(header instanceof HTMLElement))
        throw new Error('No catalog grid or header');
      await awaitCheck(() => Math.abs(grid.getBoundingClientRect().width - header.getBoundingClientRect().width) < 1,
        'The catalog grid does not take the width of its pane', WAIT_MS);
      expect(view.models.col('id') !== null, true, 'The catalog has no id column');
      const captions = {name: 'Name', engine_name: 'Method', task: 'Task', target_name: 'Target',
        storage_mode: 'Data storage', row_count: 'Training rows', tags: 'Tags', created_on: 'Created'};
      for (const [column, caption] of Object.entries(captions))
        expect(view.models.col(column)?.meta.friendlyName, caption, column);
    } finally {
      view.close();
    }
    return measured;
  });

  test('the catalog refreshes itself and the delete icon follows the current row', async () => {
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
      view.models.currentRowIdx = -1;
      await awaitCheck(() => isDisabled(view.deleteIcon) && isDisabled(view.applyIcon),
        'Delete or Apply is enabled without a current row', WAIT_MS);
      view.models.currentRowIdx = rowOf();
      await awaitCheck(() => !isDisabled(view.deleteIcon) && !isDisabled(view.applyIcon),
        'Delete or Apply is disabled with a current row', WAIT_MS);

      await deleteModel(modelId);
      modelId = undefined;
      await awaitCheck(() => rowOf() < 0, 'The deleted model stays', WAIT_MS);
      expect(isDisabled(view.deleteIcon), view.models.currentRowIdx < 0, 'Delete does not follow the current row');
    } finally {
      if (modelId !== undefined)
        await deleteModel(modelId);
      view.close();
    }
  }, {timeout: 30000});

  test('app function is registered', async () => {
    const apps = DG.Func.find({package: 'Forge', tags: [DG.FUNC_TYPES.APP]});
    expect(apps.length, 1);
    expect(apps[0].friendlyName, APP_NAME);
  });

  test('menu function is registered', async () => {
    const topMenu = DG.Func.find({package: 'Forge', name: 'forgeModels'})[0]?.topMenu;
    expect(topMenu?.startsWith(MENU_PATH), true, `Unexpected top menu: ${topMenu}`);
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
      expect(forms[0].contains(view.trainButton), true, 'Train is not in the form');
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
      for (const caption of ['Accuracy', 'F1'])
        expect(metrics.includes(caption), true, `${caption} is missing: ${metrics}`);
      expect(metrics.includes('Sensitivity'), false, 'Sensitivity on three classes');
      expect(text.includes('Seed: '), true, 'No seed line');
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
      expectArray([...view.methodInput.items].sort(), IRIS_REGRESSORS);
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
      expectArray([...view.methodInput.items].sort(), IRIS_CLASSIFIERS);
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
      const balloon = 'Softmax cannot be used with this selection; Linear Regression is chosen.';
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
      expect(isShown(data.body) && isShown(preparation.body) && isShown(method.body), true, 'A group starts collapsed');
      expect(data.body.contains(view.featuresInput.root), true, 'Features is not under Data');
      expect(preparation.body.contains(view.missingValuesInputs.choice.root), true,
        'Missing values is not under Preparation');
      const groupRoot = view.preparationInputs.group.root;
      expect(!isShown(groupRoot) && groupRoot.offsetHeight === 0, true, 'Preparation is shown for iris without gaps');
      expectArray(view.groups.map((g) => Array.from(g.root.parentElement?.children ?? []).indexOf(g.root)), [0, 1, 2]);
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

  test('train function is registered', async () => {
    const topMenu = DG.Func.find({package: 'Forge', name: 'forgeTrain'})[0]?.topMenu;
    expect(topMenu?.startsWith(MENU_PATH), true, `Unexpected top menu: ${topMenu}`);
  });

  test('apply functions are registered', async () => {
    const topMenu = DG.Func.find({package: 'Forge', name: 'forgeApply'})[0]?.topMenu;
    expect(topMenu?.startsWith(MENU_PATH), true, `Unexpected top menu: ${topMenu}`);
    const inputs = DG.Func.find({package: 'Forge', name: 'applyModel'})[0]?.inputs.map((p) => p.name) ?? [];
    expectArray(inputs, ['model', 'table', 'columnNamesMap', 'showProgress']);
  });

  test('Train view shows Missing values for a gapped feature and gives the columns back', async () => {
    const iris = await openIris();
    iris.getCol('Sepal.Width').set(3, null);
    const runs = () => runsOf([iris]);
    const view = await openTrainView(iris);
    try {
      view.methodInput.value = 'SVM';
      await awaitTrained(view);
      const choice = view.missingValuesInputs.choice;
      expect(isShown(choice.root), true, 'Missing values is hidden with a gap');
      expect(choice.inputType, 'Radio');
      expect(choice.value, 'Skip rows');
      expectOptionsInRow(choice);
      // The skipped row makes every training column a copy: only the Results pane's frames hold a table column.
      const [, frames] = await framesSharing(iris.columns.toList(), () => view.train());
      expectNoFrames(frames.filter((f) => !view.resultsPane.frames.includes(f)));
      expect(Array.from(view.root.querySelectorAll('ul > li'), (li) => li.textContent)
        .includes('Rows: 149 used, 1 skipped (missing values)'), true, 'The Results rows item is missing');
      expect((await runs()).find((r) => r.id === view.lastTraining?.runId)?.row_count, 149);

      const imputeLabel = Array.from(choice.root.querySelectorAll('.ui-radio-button label'))
        .find((label) => label.textContent === 'Impute');
      if (!(imputeLabel instanceof HTMLElement))
        throw new Error('No Impute option');
      imputeLabel.click();
      expect(choice.value, 'Impute');
      const neighbors = view.missingValuesInputs.inputs.find((input) => input.caption === 'Neighbors');
      if (neighbors === undefined)
        throw new Error('No Neighbors input');
      await awaitCheck(() => isShown(neighbors.root), 'Neighbors does not appear', WAIT_MS);
      neighbors.value = 0;
      await awaitCheck(() => isDisabled(view.trainButton) && !view.isTraining, 'Train is enabled with Neighbors 0',
        TIMEOUT);
      expect(view.lastTraining === undefined, true, 'Neighbors 0 trained');
      await awaitCheck(() => resultsText(view).includes('Fix the settings.'),
        'The Results hint is not Fix the settings. for Neighbors 0', WAIT_MS);
      expect(isTrainShown(view), false, 'Train is shown for SVM with Neighbors 0');
      neighbors.value = 3;
      await awaitTrained(view);
      expect(isTrainShown(view), false, 'Train is shown for SVM with Neighbors 3');

      view.featuresInput.value = columnsOf(iris, ['Sepal.Length', 'Petal.Length', 'Petal.Width']);
      await awaitCheck(() => !isShown(choice.root), 'Missing values stays without gaps', WAIT_MS);
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: TIMEOUT});
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
      expect(isShown(view.preparationInputs.group.root), false, 'Preparation is shown for XGBoost');
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
    const rows = 10001;
    const a = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'a', valuesOf(rows, (i) => i % 2 + (i % 7) / 100));
    const b = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'b', valuesOf(rows, (i) => (i % 11) / 10));
    const label = DG.Column.fromStrings('label', Array.from({length: rows}, (_, i) => i % 2 === 0 ? 'even' : 'odd'));
    const table = DG.DataFrame.fromColumns([a, b, label]);
    table.name = `forge-test-slow-${Date.now()}`;
    const view = await openTrainView(table);
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

  test('a training completed after a cancelled one clears the invalid mark', async () => {
    const demog = demoDemog(50000);
    const view = await TrainView.create(demog);
    view.targetInput.value = demog.getCol('sex');
    view.methodInput.value = 'SVM';
    grok.shell.addView(view);
    try {
      await awaitCheck(() => isTrainShown(view) && !isDisabled(view.trainButton), 'Train is not offered for SVM',
        WAIT_MS);
      const cancelled = view.train();
      await cancelOnceShown('Training SVM model', WAIT_MS);
      await cancelled;
      expect(`${view.runs[0]?.status} ${resultsState(view)}`, 'cancelled false/true',
        'The cancel came after the training ended');
      await view.train();
      await awaitCheck(() => resultsState(view).startsWith('false/'), 'The results stay busy', WAIT_MS);
      expect(`${view.runs[0]?.status} ${resultsState(view)}`, 'completed false/false');
      expect(shownRun(view), view.lastTraining?.runId ?? null, 'The results show another run than the last training');
    } finally {
      await closeTrainView(view, [demog]);
    }
  }, {timeout: TIMEOUT});

  test('Preparation follows demog: Missing values, One-hot encoding, Skip unique categories', async () => {
    const demog = await openDemog();
    const view = await TrainView.create(demog);
    grok.shell.addView(view);
    const preparation = view.preparationInputs;
    const {oneHot, skipUniqueCategories, predictProbability, cutoff} = preparation;
    try {
      view.methodInput.value = 'SVM';
      await awaitTrained(view);
      expect(isShown(preparation.group.root) && isShown(view.missingValuesInputs.choice.root), true,
        'Preparation does not show Missing values for HEIGHT and WEIGHT');
      expect([oneHot, skipUniqueCategories, predictProbability, cutoff].some((input) => isShown(input.root)), false,
        'A step is shown for numerical features and a many-class target');
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

  test('Predict probability trains on the classes; the cutoff re-cuts without a run', async () => {
    const demog = await openDemog();
    const view = await TrainView.create(demog);
    grok.shell.addView(view);
    const {predictProbability, cutoff} = view.preparationInputs;
    const models = () => forgeDb.models.query().where('dataset_name', '=', demog.name);
    try {
      view.targetInput.value = demog.getCol('SEX');
      view.featuresInput.value = columnsOf(demog, ['AGE', 'HEIGHT', 'WEIGHT']);
      await awaitCheck(() => isShown(predictProbability.root), 'Predict probability is not offered for SEX', WAIT_MS);
      expect(predictProbability.value || isShown(cutoff.root), false, 'Predict probability starts checked');
      await awaitTrained(view);
      expect(view.lastTraining?.result.task, 'classification');
      const classificationGrid = gridWith(view.root, 'Metric').root;

      predictProbability.value = true;
      await awaitCheck(() => isShown(cutoff.root), 'No Positive class cutoff', WAIT_MS);
      expect(cutoff.value, 0.5);
      await awaitTrained(view);
      const trained = view.lastTraining?.result;
      expect(trained?.task, 'regression');
      expect(view.methodInput.items.includes('Softmax'), false, 'Softmax is listed for the probability');
      const auc = trained?.metrics.validation.auc;
      expect(auc !== undefined && auc > 0.5 && auc <= 1, true, `Validation AUC-ROC ${auc}`);
      expect(gridWith(view.root, 'Metric').dataFrame.getCol('Metric').toList().includes('AUC-ROC'), true,
        'No AUC-ROC row in Results');
      const probabilityGrid = gridWith(view.root, 'Metric');
      expect(probabilityGrid.root === classificationGrid, false, 'A task change kept the classification grid');
      const runCount = (await runsOf([demog])).length;

      cutoff.value = 0.9;
      await awaitCheck(() => view.lastTraining?.result.options.binaryClassificationThreshold === 0.9,
        'The cutoff does not re-cut the result', WAIT_MS);
      await delay(CHECK_DELAY_MS * 3);
      expect(view.isTraining, false, 'The cutoff retrains');
      expect((await runsOf([demog])).length, runCount, 'The cutoff recorded a run');
      const recut = view.lastTraining?.result.metrics.train;
      expect(recut?.auc, trained?.metrics.train.auc, 'AUC-ROC changed with the cutoff');
      expect(recut?.sensitivity !== trained?.metrics.train.sensitivity ||
        recut?.specificity !== trained?.metrics.train.specificity, true, 'The cutoff left the shares as they were');
      const recutGrid = gridWith(view.root, 'Metric');
      expect(recutGrid.root === probabilityGrid.root, true, 'The cutoff built a new metrics grid');
      const sensitivityRow = recutGrid.dataFrame.getCol('Metric').toList().indexOf('Sensitivity');
      // The column stores single-precision floats.
      expect(Math.abs((recutGrid.dataFrame.getCol('Train').get(sensitivityRow) ?? NaN) - (recut?.sensitivity ?? NaN)) <
        1e-6, true, 'The grid shows the old sensitivity');
      expect(gridTooltip(recutGrid, recutGrid.cell('Train', sensitivityRow)).startsWith(`${recut?.sensitivity}.`), true,
        'The cell tooltip shows the old sensitivity');

      cutoff.value = 1.5;
      await awaitCheck(() => cutoff.validity !== null, 'No mark on a cutoff of 1.5', WAIT_MS);
      await awaitCheck(() => resultsText(view).includes('Fix the settings.'), 'No hint for a cutoff of 1.5', WAIT_MS);
      expect(isDisabled(view.saveButton), true, 'Save stays enabled with a cutoff of 1.5');
      expect(view.lastTraining?.result.options.binaryClassificationThreshold, 0.9, 'A cutoff of 1.5 re-cut the result');
      cutoff.value = 0.9;
      await awaitCheck(() => cutoff.validity === null && !resultsText(view).includes('Fix the settings.'),
        'A valid cutoff leaves the hint', WAIT_MS);
      expect(isDisabled(view.saveButton), false, 'Save stays disabled with a valid cutoff');

      const name = `forge-test-model-${Date.now()}`;
      await view.saveModelAs({name, description: '', tags: []});
      const [saved] = await models();
      const options = preparationOptionsOf(saved?.options);
      expect(options.binaryClassificationThreshold, 0.9);
      expect(metricsRecordOf(saved?.metrics)?.train.sensitivity, recut?.sensitivity);

      // A change while the cutoff is invalid leaves nothing to re-cut; fixing the cutoff then trains.
      cutoff.value = 1.5;
      await awaitCheck(() => resultsText(view).includes('Fix the settings.'), 'No hint for a cutoff of 1.5', WAIT_MS);
      view.featuresInput.value = columnsOf(demog, ['AGE', 'HEIGHT']);
      await delay(CHECK_DELAY_MS * 5);
      await awaitCheck(() => !view.isTraining && resultsText(view).includes('Fix the settings.'),
        'A change with an invalid cutoff does not show the hint', WAIT_MS);
      expect(view.lastTraining === undefined, true, 'A change with an invalid cutoff trained');
      cutoff.value = 0.5;
      await awaitTrained(view);
      expect(view.lastTraining?.result.options.binaryClassificationThreshold, 0.5);
      expect(resultsText(view).includes('Fix the settings.'), false, 'Results still shows the hint');
      expect(gridWith(view.root, 'Metric').dataFrame.getCol('Metric').toList().includes('AUC-ROC'), true,
        'Results shows no metrics');
    } finally {
      for (const model of await models())
        await deleteModel(model.id);
      await closeTrainView(view, [demog]);
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

/** The text of the Insights & Tips block of [view]; '' while it is hidden. */
function insightsText(view: TrainView): string {
  const block = view.resultsPane.root.querySelector('.forge-insights');
  return block instanceof HTMLElement && isShown(block) ? block.textContent ?? '' : '';
}

/** The text of each line of Insights & Tips of [view], its pointer line included. */
const insightLines = (view: TrainView) =>
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

/** Checks that every panel title of [view] answers the mouse with a tooltip. */
async function expectTitleTooltips(view: TrainView): Promise<void> {
  for (const span of view.resultsPane.root.querySelectorAll<HTMLElement>('.forge-results-title > span')) {
    expect(getComputedStyle(span).pointerEvents, 'auto', `${span.textContent}: the title lets the mouse through`);
    expect((await enterTooltip(span)).length > 10, true, `${span.textContent} has no tooltip`);
  }
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

/** Predict probability on two iris species: ROC Curve over its own frame, kept by a cutoff change and a retrain. */
async function expectRocCurve(): Promise<void> {
  const iris = await twoSpeciesIris();
  const view = await openTrainView(iris);
  const {predictProbability, cutoff} = view.preparationInputs;
  try {
    predictProbability.value = true;
    const rocTitle = () => rowTitles(view).flat().find((t) => t.startsWith('ROC Curve')) ?? '';
    await awaitCheck(() => rocTitle() !== '' && view.lastTraining?.result.scores !== undefined && !view.isTraining,
      'No ROC Curve', TIMEOUT);
    const aucTitle = () => `ROC Curve (AUC = ${view.lastTraining?.result.metrics.train.auc?.toFixed(3)})`;
    expect(rocTitle(), aucTitle());
    const roc = panelViewer(view, rocTitle());
    const frame = roc.dataFrame;
    expect(roc.type, DG.VIEWER.SCATTER_PLOT);
    const [fprName, tprName] = ['False Positive Rate', 'True Positive Rate'];
    expectArray([roc.props.xColumnName, roc.props.yColumnName, roc.props.linesOrderColumnName],
      [fprName, tprName, fprName]);
    expectArray(frame.columns.names(), ['Threshold', fprName, tprName]);
    expect(frame.rowCount, 101);
    expect(view.resultsPane.frames.some((f) => f.dart === frame.dart), false, 'ROC Curve is on the preview frame');
    expect(frame.meta.formulaLines.items.some((line) => line.title === 'Non-informative baseline'), true,
      'No baseline');
    expect(frame.currentRowIdx, 50, 'The cutoff 0.5 is not marked');
    const content = (title: string) => panelOf(view, title).lastElementChild;
    const confusions = content('Confusions');
    let curveChanges = 0;
    const sub = frame.onValuesChanged.subscribe(() => curveChanges++);
    try {
      cutoff.value = 0.9;
      await awaitCheck(() => view.lastTraining?.result.options.binaryClassificationThreshold === 0.9,
        'The cutoff does not re-cut the result', WAIT_MS);
      await delay(CHECK_DELAY_MS * 3);
    } finally {
      sub.unsubscribe();
    }
    const recut = panelViewer(view, rocTitle());
    expect(recut.root === roc.root && roc.root.isConnected && recut.dataFrame.dart === frame.dart, true,
      'The cutoff rebuilt ROC Curve');
    expect(curveChanges, 0, 'The cutoff changed the curve');
    expect(frame.currentRowIdx, 90, 'The mark did not move to the cutoff 0.9');
    expect(content('Confusions') === confusions, true, 'The cutoff rebuilt Confusions');

    const fpr = frame.getCol(fprName).toList();
    const tpr = frame.getCol(tprName).toList();
    const before = view.lastTraining;
    if (before?.engine.name === 'XGBoost')
      setHyperparameter(view, 'iterations', 1);
    else
      view.methodInput.value = 'XGBoost';
    await awaitCheck(() => view.lastTraining !== before && view.lastTraining?.engine.name === 'XGBoost' &&
      !view.isTraining, 'No XGBoost retrain', TIMEOUT);
    await delay(CHECK_DELAY_MS * 3);
    const retrained = panelViewer(view, rocTitle());
    expect(retrained.root === roc.root && retrained.dataFrame.dart === frame.dart && frame.rowCount === 101, true,
      'The retrain rebuilt ROC Curve or its frame');
    expect(rocTitle(), aucTitle());
    expect(frame.getCol(fprName).toList().some((v, i) => v !== fpr[i]) ||
      frame.getCol(tprName).toList().some((v, i) => v !== tpr[i]), true, 'The retrain left the curve as it was');
    expect(frame.currentRowIdx, 90, 'The retrain moved the mark');
  } finally {
    await closeTrainView(view, [iris]);
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
    const all = new Set<PanelKey>(['predictedVsActual', 'rocCurve', 'metrics', 'confusions', 'insights',
      'wrongPredictions', 'features', 'correlation', 'classBalance']);
    const binary = resultsLayout('binary', all, 0);
    expectRows(binary.map((r) => r.panels), [['predictedVsActual', 'rocCurve', 'metrics'],
      ['insights', 'confusions', 'wrongPredictions'], ['features', 'correlation', 'classBalance']]);
    expectArray(binary.map((r) => r.height), [350, 300, 450]);
    expectRows(binary.map((r) => r.weights), [[2, 2, 1], [1, 1.5, 2.5], [1, 1, 1]]);
    const multiclass = resultsLayout('multiclass', new Set<PanelKey>(['predictedVsActual', 'features', 'correlation',
      'classBalance']), 0);
    expectRows(multiclass.map((r) => r.panels), [['predictedVsActual'], ['metrics', 'classBalance'],
      ['features', 'correlation']]);
    expectRows(multiclass.map((r) => r.weights), [[1], [1.5, 1.5], [1, 1]]);
    const regression = resultsLayout('regression', new Set<PanelKey>(['predictedVsActual', 'residuals',
      'distributions', 'statistics']), 4);
    expectRows(regression.map((r) => r.panels), [['predictedVsActual', 'residuals', 'distributions'],
      ['metrics', 'statistics'], ['visualize:0', 'visualize:1', 'visualize:2'], ['visualize:3']]);
    expectArray(regression.map((r) => r.height), [350, 300, 350, 350]);
    expectRows(regression.map((r) => r.weights), [[1, 1, 1], [1.5, 1.5], [1, 1, 1], [1]]);
    expectRows(resultsLayout('binary', new Set<PanelKey>(['metrics']), 0).map((r) => r.panels), [['metrics']]);
  });

  test('resultsLayout: a row 2 of two panels takes Class balance from row 3', async () => {
    const withoutWrong = new Set<PanelKey>(['predictedVsActual', 'confusions', 'insights', 'features', 'correlation',
      'classBalance']);
    expectRows(resultsLayout('multiclass', withoutWrong, 0).map((r) => r.panels), [['predictedVsActual', 'confusions'],
      ['insights', 'metrics', 'classBalance'], ['features', 'correlation']]);
    expectRows(resultsLayout('binary', withoutWrong, 0).map((r) => r.panels), [['predictedVsActual', 'metrics'],
      ['insights', 'confusions', 'classBalance'], ['features', 'correlation']]);
    const withoutInsights = new Set<PanelKey>(['predictedVsActual', 'confusions', 'wrongPredictions', 'classBalance']);
    expectRows(resultsLayout('binary', withoutInsights, 0).map((r) => r.panels), [['predictedVsActual', 'metrics'],
      ['confusions', 'wrongPredictions', 'classBalance']]);
    const full = new Set<PanelKey>([...withoutWrong, 'wrongPredictions']);
    expectRows(resultsLayout('multiclass', full, 0).map((r) => r.panels), [['predictedVsActual', 'confusions'],
      ['insights', 'metrics', 'wrongPredictions'], ['features', 'correlation', 'classBalance']]);
  });

  test('iris: Insights & Tips, the classifier rows, Class balance, the frames given back on close', async () => {
    const iris = await openIris();
    const [view, frames] = await framesSharing(iris.columns.toList(), () => openTrainView(iris));
    try {
      expect(insightsText(view).includes('Insights & Tips'), true, 'No Insights & Tips');
      expect(insightsText(view).includes('Features are highly correlated.'), true, insightsText(view));
      const training = view.lastTraining;
      if (training === undefined)
        throw new Error('No training');
      const labels = trainPredictionOf(training.result);
      const target = training.request.target;
      const hasWrong = Array.from({length: target.length}, (_, i) => i)
        .some((i) => String(target.get(i)) !== String(labels.get(i)));
      // Without wrong predictions row 2 has two panels and takes Class balance.
      expectRows(rowTitles(view), [['Predicted Species vs Actual', 'Confusions'],
        ['Insights & Tips', 'Performance', hasWrong ? 'Wrong predictions' : 'Class balance'],
        ['Features', 'Features correlation', ...(hasWrong ? ['Class balance'] : [])]]);
      if (hasWrong)
        expect(gridWith(view.resultsPane.root, 'Predicted Species').root.classList.contains('d4-grid'), true);

      const chart = panelViewer(view, 'Class balance');
      expect(chart.type, DG.VIEWER.BAR_CHART);
      expect(chart.dataFrame.dart === view.resultsPane.frames[0].dart, true, 'Class balance has a frame of its own');
      expect(chart.props.splitColumnName, target.name);
      expect(chart.props.stackColumnName, target.name);
      expect(chart.props.valueAggrType, DG.AGG.TOTAL_COUNT);
      expect(chart.props.showStackSelector || chart.props.showCategorySelector || chart.props.showValueSelector, false,
        'Class balance shows a column selector');
      expect(panelViewer(view, 'Features correlation').props.description ?? '', 'Features are highly correlated.');
      // The panel shows the short line only; the pairs and the pointer are in its tooltip.
      const correlated = insightItem(view, 'Features are highly correlated.');
      expect(correlated.textContent, 'Features are highly correlated.');
      const correlatedTooltip = await enterTooltip(correlated);
      expect(correlatedTooltip.includes('Petal.Length and Petal.Width: 0.96') &&
        correlatedTooltip.endsWith('See Features correlation.'), true, correlatedTooltip);
      expect(insightsText(view).includes('See '), false, `A pointer is on the panel: ${insightsText(view)}`);
      await expectTitleTooltips(view);
      expectUnderTitles(view);

      // The panels start in the ratio of their weights and every viewer fills its panel.
      const [first, second, third] = rowsOf(view);
      await awaitCheck(() => slotsOf(first).every((s) => s.offsetWidth > 0), 'Row 1 is not laid out', WAIT_MS);
      const ratio = (row: HTMLElement, i: number) => slotsOf(row)[i].offsetWidth / slotsOf(row)[0].offsetWidth;
      expectFloat(ratio(first, 1), 1, 0.05);
      expectFloat(ratio(second, 1), 1.5, 0.05);
      expectFloat(ratio(third, 1), 1, 0.05);
      for (const panel of view.resultsPane.root.querySelectorAll<HTMLElement>(
        '.forge-results-panel:not(.forge-results-metrics):not(.forge-insights)')) {
        const content = panel.lastElementChild;
        expect(content instanceof HTMLElement && Math.abs(content.offsetWidth - panel.clientWidth) <= 1 &&
          content.offsetHeight > 0, true, `${titleOf(panel)}: ${content instanceof HTMLElement ?
          content.offsetWidth : 0} of ${panel.clientWidth} px`);
      }
      expect(panelOf(view, 'Performance').firstElementChild?.classList.contains('d4-dialog-header'), true,
        'Performance has its own header');
      const shown = view.resultsPane.frames;
      expect(shown.length > 0 && shown.every((f) => frames.includes(f) && f.columns.length > 0), true,
        'The previews hold no frame of the table\'s columns');
    } finally {
      await closeTrainView(view, [iris]);
    }
    expectReleased(frames);
  }, {timeout: TIMEOUT});

  test('winequality: the regression rows, Statistics and the method\'s own views', async () => {
    const wine = await openDemoTable(WINE, 'wine');
    const view = await openTrainView(wine);
    try {
      view.targetInput.value = wine.getCol('quality');
      await awaitCheck(() => view.featuresInput.validity !== null, 'Features is not marked', WAIT_MS);
      view.featuresInput.value = view.featuresInput.value.filter((c) => c.name !== 'quality');
      await awaitTrained(view);
      expect(view.lastTraining?.engine.name, 'PLS Regression');
      await awaitCheck(() => rowTitles(view).length > 3, `No PLS Regression row: ${JSON.stringify(rowTitles(view))}`,
        TIMEOUT);
      const titles = rowTitles(view);
      const second = titles[1].filter((t) => t !== 'Insights & Tips');
      expectRows([titles[0], second, titles[2]], [['Predicted quality vs Actual', 'Residuals',
        'Distributions of predicted and actual quality'], ['Performance', 'Statistics'],
      ['Features', 'Features correlation']]);
      expect(titles[1].indexOf('Insights & Tips') <= 0, true, 'Insights & Tips is not first in row 2');
      expect(titles.slice(3).every((row) => row.length > 0 && row.length <= 3 &&
        row.every((t) => t === 'PLS Regression')), true, JSON.stringify(titles));
      const statistics = gridWith(view.resultsPane.root, 'Statistic');
      expect(statistics.dataFrame.rowCount, 11);
      expect(headerTooltip(statistics, 'Predicted'), 'The statistic of the model\'s predictions on the same rows.');
      // One line for every statistic far from the target's, pointing to Statistics and shown as its description.
      const deviation = 'The predictions\' statistics deviate significantly from the actual ones.';
      const lines = insightLines(view);
      expect(lines.filter((l) => l.includes('deviate significantly')).length <= 1 &&
        !lines.some((l) => l.includes('is far from the actual one')), true, lines.join(' | '));
      const deviating = lines.filter((l) => l.startsWith(deviation));
      expectArray(deviating, deviating.length === 0 ? [] : [deviation]);
      expect(statistics.props.description ?? '', deviating.length === 0 ? '' : deviation);
      if (deviating.length > 0) {
        const tooltip = await enterTooltip(insightItem(view, deviation));
        expect(tooltip.endsWith('See Statistics.'), true, tooltip);
      }
      await expectTitleTooltips(view);
      expectUnderTitles(view);
      const r2 = view.lastTraining?.result.metrics.train.r2 ?? NaN;
      const metrics = gridWith(view.resultsPane.root, 'Metric');
      const tooltip = gridTooltip(metrics, metrics.cell('Train', metricRow(metrics, 'R2')));
      expect(tooltip.endsWith('R2 is below 0.3.'), r2 < 0.3, `Train R2 ${r2}: ${tooltip}`);
    } finally {
      await closeTrainView(view, [wine]);
    }
  }, {timeout: 120000});

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

  test('ROC Curve: a line over its own frame of 101 thresholds, the AUC of Performance; a cutoff moves only its ' +
    'mark, a retrain refills the frame', expectRocCurve, {timeout: 120000});

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
    const tooltip = gridTooltip(grid, grid.cell('Train', metricRow(grid, 'R2')));
    expect(tooltip.endsWith('R2 is below 0.3.'), true, tooltip);
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
        expectRows(titles[0], ['Predicted Species vs Actual', 'Confusions']);
        expect(titles[1]?.includes('Performance'), true, JSON.stringify(titles));
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
  box.value = text;
  box.dispatchEvent(new Event('input', {bubbles: true}));
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
    const iris = await openIris();
    catalogModel = {id: (await saveIrisModel(iris)).id, iris};
  });

  after(async () => {
    if (catalogModel !== undefined)
      await deleteModel(catalogModel.id);
  });

  test('the autostart registered the model handler', async () => {
    expect(DG.ObjectHandler.list().some((h) => h.name === 'Forge model handler'), true);
    const row = ForgeModelHandler.rowOf(await forgeDb.models.get(sharedModel().id));
    expect(DG.ObjectHandler.forEntity(row)?.name, 'Forge model handler');
  });

  test('the catalog header has Applicable to and a disabled Compare icon', async () => {
    const view = await ForgeApp.create();
    grok.shell.addView(view);
    try {
      const header = view.root.querySelector('.forge-pane-header');
      expect(header?.contains(view.compareIcon) && header.contains(view.applyIcon), true, 'An icon is not in header');
      expect(header?.nextElementSibling, view.applicableToInput.root, 'Applicable to is not on the line under Models');
      const colorOf = (e: Element | null | undefined) => e instanceof Element ? getComputedStyle(e).color : 'none';
      const grey = ui.div([]);
      grey.style.color = 'var(--grey-3)';
      view.root.append(grey);
      const refresh = Array.from(header?.querySelectorAll('.grok-icon') ?? [])
        .find((i) => i instanceof HTMLElement && !isDisabled(i));
      expect(colorOf(refresh), colorOf(view.applicableToInput.root.querySelector('.ui-input-options > i')),
        'An enabled header icon is not the blue of the input\'s icons');
      expect(colorOf(view.compareIcon), colorOf(grey), 'A disabled header icon is not the platform\'s grey');
      grey.remove();
      expect(view.applicableToInput.inputType, DG.InputType.Table);
      expect(view.applicableToInput.value, null);
      expect(isDisabled(view.compareIcon), true, 'Compare is enabled without a selection');
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
      expect(isDisabled(view.compareIcon), true, 'Compare is enabled with one selected row');
      view.models.selection.set(ids.indexOf(otherId), true);
      await awaitCheck(() => !isDisabled(view.compareIcon), 'Compare stays disabled with two selected rows', WAIT_MS);
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
      expect(isDisabled(view.applyIcon) || isDisabled(view.compareIcon), false, 'An icon greyed out');
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
      view.deleteIcon.click();
      const confirm = DG.Dialog.getOpenDialogs().find((d) => d.title === 'Delete model');
      if (confirm === undefined)
        throw new Error('No Delete model dialog');
      confirm.getButton('OK').click();
      await awaitCheck(() => indexOf(shownId) < 0, 'The deleted model stays in the catalog', WAIT_MS);
      await awaitCheck(() => !shows(shownId), 'The context panel keeps the deleted model', WAIT_MS);

      view.models.selection.set(indexOf(id), true);
      view.models.selection.set(indexOf(comparedId), true);
      await awaitCheck(() => comparedIds() === [id, comparedId].sort().join(), 'The two are not compared', WAIT_MS);
      const menuDelete = CATALOG_ACTIONS.find((a) => a.name === 'Delete model');
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
      expect(tooltip.includes('Predicts Species from 4 features with XGBoost'), true, tooltip);
    } finally {
      ui.tooltip.hide();
      accordion.root.remove();
    }
  });

  test('dragging the panel title\'s name starts the drag of the model row, also for a row without a file', async () => {
    const manualId = await insertModelRow({name: `forge-test-model-${Date.now()}-manual`});
    try {
      for (const id of [sharedModel().id, manualId]) {
        const model = await forgeDb.models.get(id);
        grok.shell.setCurrentObject(ForgeModelHandler.rowOf({id: model.id, name: model.name}), false, true);
        await awaitCheck(() => titleName(document, model.name) !== undefined,
          'The context panel does not show the model', WAIT_MS);
        const name = titleName(document, model.name);
        if (name === undefined)
          throw new Error('No name in the panel title');
        expect(getComputedStyle(name).userSelect, 'none');
        const before = new Set(Array.from(document.querySelectorAll('.d4-balloon')));
        const r = name.getBoundingClientRect();
        const at = (type: string, dx: number) => new MouseEvent(type, {bubbles: true, cancelable: true, button: 0,
          buttons: 1, clientX: r.left + 2 + dx, clientY: r.top + 2});
        name.dispatchEvent(at('mousedown', 0));
        for (const dx of [3, 10, 20, 40])
          document.dispatchEvent(at('mousemove', dx));
        try {
          expect(document.body.classList.contains(`d4-drag-${MODEL_TYPE}`), true, document.body.className);
        } finally {
          document.dispatchEvent(at('mouseup', 40));
        }
        await delay(300);
        const errors = Array.from(document.querySelectorAll('.d4-balloon'))
          .filter((b) => !before.has(b) && (b.textContent ?? '').includes('Error'));
        expect(errors.length, 0, errors.map((b) => b.textContent).join(' | '));
      }
    } finally {
      grok.shell.setCurrentObject(null, false, true);
      await forgeDb.models.delete(manualId);
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
      await expectPaneText(accordion, 'Performance', [`Seed: ${model.seed}`]);
      const performance = accordion.getPane('Performance').root;
      const metrics = gridWith(performance, 'Metric');
      expectArray(metrics.dataFrame.columns.names(), ['Metric', 'Train', 'Validation']);
      expect(metrics.dataFrame.getCol('Metric').toList().includes('Accuracy') && !metrics.props.allowEdit, true,
        'No Accuracy, or the metrics are editable');
      expect(headerTooltip(metrics, 'Validation').startsWith('Value on rows the model did not see'), true,
        'No Validation header tooltip');
      expect(gridTooltip(metrics, metrics.cell('Metric', 0)), METRIC_DESCRIPTIONS.accuracy);
      expectArray(Array.from(performance.querySelectorAll('ul > li'), (li) => li.textContent?.trim()),
        ['Validation: 5-fold cross-validation on 150 rows', `Seed: ${model.seed}`]);
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
    expectArray(CATALOG_ACTIONS.map((a) => a.name), ['Apply...', 'Edit model...', 'Download', 'Delete model']);
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

  test('the context panel follows Forge:updateModel of its model only', async () => {
    const stamp = Date.now();
    const name = `forge-test-model-${stamp}-followed`;
    const id = await insertModelRow({name, description: 'd', tags: 'a'});
    const otherId = await insertModelRow({name: `forge-test-model-${stamp}-other`});
    const renamed = `${name}-renamed`;
    const details = () => titleName(document, renamed)?.closest('.d4-accordion')?.textContent ?? '';
    try {
      grok.shell.setCurrentObject(ForgeModelHandler.rowOf(await forgeDb.models.get(id)), false, true);
      await awaitCheck(() => titleName(document, name) !== undefined, 'The context panel does not show the model',
        WAIT_MS);
      const title = titleName(document, name);
      await grok.functions.call('Forge:updateModel', {model: otherId, description: 'other'});
      await delay(1500);
      expect(title?.isConnected, true, 'Another model\'s change rebuilt the context panel');
      await grok.functions.call('Forge:updateModel', {model: id, name: renamed, description: 'forge-test-edited',
        tags: 'b'});
      await awaitCheck(() => details().includes('forge-test-edited'), `The Details do not follow: ${details()}`,
        WAIT_MS);
      const shown: unknown = grok.shell.o;
      expect(shown instanceof DG.DomainRow && shown.id === id && shown.values.tags === 'b', true,
        `The context panel shows ${shown}`);
    } finally {
      grok.shell.setCurrentObject(null, false, true);
      await forgeDb.models.delete(id);
      await forgeDb.models.delete(otherId);
    }
  }, {timeout: TIMEOUT});

  test('the context panel leaves the model Forge:deleteModel deletes', async () => {
    const name = `forge-test-model-${Date.now()}-gone`;
    const id = await insertModelRow({name});
    try {
      grok.shell.setCurrentObject(ForgeModelHandler.rowOf(await forgeDb.models.get(id)), false, true);
      await awaitCheck(() => titleName(document, name) !== undefined, 'The context panel does not show the model',
        WAIT_MS);
      await grok.functions.call('Forge:deleteModel', {model: id});
      await awaitCheck(() => !shows(id) && titleName(document, name) === undefined,
        'The context panel keeps the deleted model', WAIT_MS);
    } finally {
      if (await forgeDb.models.query().where('id', '=', id).select('name').first() !== null)
        await forgeDb.models.delete(id);
    }
  }, {timeout: TIMEOUT});

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
    const iris = await openIris();
    const {id, name} = await saveIrisModel(iris);
    paneModel = {id, name, iris};
  });

  after(async () => {
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

  test('the model name has the built-in tooltip and the model\'s menu', async () => {
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
      await awaitCheck(() => ['Created by', 'Trained on', 'Inputs', 'Outputs', 'Method', 'Sepal.Length']
        .every((t) => (ui.tooltip.root.textContent ?? '').includes(t)),
      `The tooltip is not the built-in one: ${ui.tooltip.root.textContent}`, WAIT_MS);
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

  test('the panel function\'s pane shows the built-in tooltip on a name, with several tables open', async () => {
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
      await awaitCheck(() => ui.tooltip.isVisible && ['Created by', 'Trained on', 'Inputs', 'Method']
        .every((t) => (ui.tooltip.root.textContent ?? '').includes(t)),
      `No built-in tooltip on the name: ${ui.tooltip.root.textContent}`, WAIT_MS);
    } finally {
      ui.tooltip.hide();
      pane.remove();
      for (const view of views)
        view.close();
      for (const t of tables.slice(1))
        grok.shell.closeTable(t);
    }
  }, {timeout: TIMEOUT});

  test('a click on a name right after the pane is expanded shows the model', async () => {
    const {name, iris} = paneFixture();
    const tableView = grok.shell.addTableView(iris);
    const isPanelShown = grok.shell.windows.showContextPanel;
    grok.shell.windows.showContextPanel = true;
    const headers = () =>
      Array.from(document.querySelectorAll<HTMLElement>('.grok-prop-panel .d4-accordion-pane-header'));
    const header = () => headers().find((h) => h.textContent?.trim() === PANEL_TABLE_MODELS);
    const panel = () => document.querySelector<HTMLElement>('.grok-prop-panel') ?? ui.div();
    try {
      // The view just added may set its own current object after this one: the table is set until it stays.
      await awaitCheck(() => {
        if (!(grok.shell.o instanceof DG.DataFrame && grok.shell.o.dart === iris.dart))
          grok.shell.setCurrentObject(iris, false, true);
        return header() !== undefined;
      }, 'No Forge models pane', WAIT_MS).catch(() => {
        throw new Error(`No Forge models pane for ${String(grok.shell.o)}: ${headers().map((h) => h.textContent)}`);
      });
      if (header()?.getAttribute('aria-expanded') !== 'true')
        header()?.click();
      await awaitCheck(() => paneLine(panel(), name) !== undefined, 'The pane does not list the model', WAIT_MS);
      // Expanding any accordion pane makes the platform drop unforced current-object changes for 2 s.
      const other = headers().find((h) => h !== header());
      if (other === undefined)
        throw new Error('The context panel has no other pane');
      if (other.getAttribute('aria-expanded') === 'true')
        other.click();
      other.click();
      paneParts(panel(), name).link.click();
      await awaitCheck(() => grok.shell.o instanceof DG.DomainRow && grok.shell.o.values?.name === name,
        `The click left ${String(grok.shell.o)} as the current object`, WAIT_MS);
    } finally {
      grok.shell.windows.showContextPanel = isPanelShown;
      tableView.close();
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
      expect(grid.dataFrame.currentRowIdx, 0);
      expect(grid.props.showCurrentRowIndicator, true);
      expect(gridTooltip(grid, grid.cell('Method', 0)).includes('iterations: 20'), true,
        gridTooltip(grid, grid.cell('Method', 0)));
      expect(gridTooltip(grid, grid.cell('Features', 0)), MEASUREMENTS.join('\n'));
      expect(gridTooltip(grid, grid.cell('Accuracy (validation)', 0)).startsWith('Train '), true);

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

  test('a cancelled training is a cancelled row; a click on it trains again', async () => {
    const demog = demoDemog(SLOW_ROWS);
    const view = await TrainView.create(demog);
    grok.shell.addView(view);
    try {
      trainSlowly(view, demog);
      await cancelOnceShown('Training XGBoost model', TIMEOUT);
      await awaitCheck(() => !view.isTraining, 'The cancelled training does not end', TIMEOUT);
      expect(view.runs[0]?.status, 'cancelled');
      const grid = runsGrid(view);
      expect(grid.dataFrame.get('Status', 0), 'cancelled');
      expect(gridTooltip(grid, grid.cell('Status', 0)), 'Stopped from the task bar.');
      clickRunsRow(grid, 0);
      await awaitCheck(() => view.isTraining, 'A click on the cancelled run does not train', WAIT_MS);
    } finally {
      await closeTrainView(view, [demog]);
    }
  }, {timeout: TIMEOUT});

  test('a click on the current run after a cutoff move restores the cutoff without training', async () => {
    const iris = await twoSpeciesIris();
    const view = await openTrainView(iris);
    const {predictProbability, cutoff} = view.preparationInputs;
    try {
      predictProbability.value = true;
      await awaitTrained(view);
      const [count, shown] = [view.runs.length, shownRun(view)];
      cutoff.value = 0.7;
      await delay(CHECK_DELAY_MS * 3);
      expect(shown !== null && shownRun(view) === shown, true, `A re-cut shows ${shownRun(view)}, not ${shown}`);
      clickRunsRow(runsGrid(view), 0);
      await delay(CHECK_DELAY_MS * 3);
      expect(cutoff.value, view.runs[0].steps.cutoff);
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

  test('a restored run sets Missing values and One-hot encoding again', async () => {
    const demog = await openDemog();
    const view = await TrainView.create(demog);
    grok.shell.addView(view);
    const {oneHot} = view.preparationInputs;
    const choice = view.missingValuesInputs.choice;
    try {
      view.methodInput.value = 'SVM';
      await awaitTrained(view);
      view.targetInput.value = demog.getCol('AGE');
      view.featuresInput.value = columnsOf(demog, ['HEIGHT', 'WEIGHT', 'SEX', 'RACE']);
      await awaitCheck(() => isShown(oneHot.root) && oneHot.value, 'One-hot encoding is not checked', WAIT_MS);
      choice.value = 'Impute';
      await awaitTrained(view);
      const imputed = view.runs[0];
      expect(imputed.missingValues.mode, 'impute');
      expect(imputed.steps.oneHot, true);
      expect(runsGrid(view).dataFrame.get('Preparation', 0).startsWith('impute ('), true,
        runsGrid(view).dataFrame.get('Preparation', 0));

      choice.value = 'Skip rows';
      oneHot.value = false;
      view.restoreRun(imputed);
      expect(choice.value, 'Impute');
      expect(oneHot.value, true);
      await awaitTrained(view);
      expect(view.runs[0].missingValues.mode, 'impute');
    } finally {
      await closeTrainView(view, [demog]);
    }
  }, {timeout: TIMEOUT});
});

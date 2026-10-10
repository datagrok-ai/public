import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {after, awaitCheck, before, category, delay, expect, expectArray, test}
  from '@datagrok-libraries/test/src/test';
import {applyAndRecord, DEFAULT_BATCH_SIZE, loadModel} from '../apply/apply-model';
import {exactMapping} from '../apply/column-matching';
import {updateModelInfo} from '../catalog/model-edit';
import {APP_NAME, MENU_PATH, MODEL_TYPE} from '../constants';
import {EngineRegistry} from '../engines/engine-registry';
import {forgeDb, ModelStorageMode} from '../generated/db';
import {METRIC_DESCRIPTIONS} from '../metrics/metrics';
import {preparationOptionsOf} from '../preparation/preparation-options';
import {datasetRefOf, storedDatasetRef} from '../storage/dataset-ref';
import {deleteModel, modelsChanged} from '../storage/model-store';
import {metricsRecordOf} from '../training/train-model';
import {applyModelDialog, ApplyDialogOptions, isOpen, modelLabels} from '../ui/apply-model-dialog';
import {ButtonGate} from '../ui/button-gate';
import {gridTooltip} from '../ui/data-grid';
import {editModelDialog} from '../ui/edit-model-dialog';
import {ForgeApp} from '../ui/forge-app';
import {CATALOG_ACTIONS, MODEL_ACTIONS} from '../ui/model-actions';
import {hasFormsViewer, ModelComparison} from '../ui/model-comparison';
import {ForgeModelHandler} from '../ui/model-handler';
import {modelAccordion, refreshSharing} from '../ui/model-panes';
import {saveModelDialog} from '../ui/save-model-dialog';
import {tagsOfInput} from '../ui/tags-input';
import {CHECK_DELAY_MS, TrainView} from '../ui/train-view';
import {columnsOf, expectNoFrames, framesSharing, insertModelRow, IRIS, MEASUREMENTS, openDemoTable, openIris,
  openIrisFromFile, savedFixture, saveIrisModel, saveSkipUniqueModel, saveTestModel, valuesOf, XGBOOST_FIELDS}
  from './test-data';

const WAIT_MS = 5000;
const TIMEOUT = 60000;

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

const runsOf = async (tables: DG.DataFrame[]) =>
  await forgeDb.trainingRuns.query().where('dataset_name', '=', tables.map((t) => t.name));

/** Closes [view], waits until no training of it is left, then deletes the runs of [tables]. */
async function closeTrainView(view: TrainView, tables: DG.DataFrame[]): Promise<void> {
  view.close();
  await awaitCheck(() => !view.isTraining, 'A training outlives the view', TIMEOUT);
  for (const run of await runsOf(tables))
    await forgeDb.trainingRuns.delete(run.id);
}

/** The open balloon that contains [text], or null. */
function balloonWith(text: string): HTMLElement | null {
  return Array.from(document.querySelectorAll<HTMLElement>('.d4-balloon'))
    .find((b) => (b.textContent ?? '').includes(text)) ?? null;
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

/** The task-bar entry whose text includes [text], or undefined. */
const taskBarEntry = (text: string) => Array.from(document.querySelectorAll('.d4-task-bar-entry'))
  .find((e) => (e.textContent ?? '').includes(text));

/** Clicks the cancel icon of the task-bar [entry], as the user does. */
function cancelTaskBarEntry(entry: Element | undefined): void {
  const cancel = entry?.querySelector('.fa-times');
  if (!(cancel instanceof HTMLElement) || !cancel.isConnected)
    throw new Error('No task-bar entry to cancel');
  cancel.click();
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

/** Reads [read] every 100 ms until [done] holds of what it read or WAIT_MS pass; the last reading. */
async function readUntil<T>(read: () => Promise<T> | T, done: (value: T) => boolean): Promise<T> {
  let value = await read();
  for (let i = 0; i < WAIT_MS / 100 && !done(value); i++) {
    await delay(100);
    value = await read();
  }
  return value;
}

/** Waits until [dialog] is closed, [table] has the Species prediction and the application of [modelId] is recorded. */
async function expectApplied(dialog: DG.Dialog, table: DG.DataFrame, modelId: string): Promise<void> {
  await awaitCheck(() => !dialog.root.isConnected, 'The dialog stays open', WAIT_MS);
  await awaitCheck(() => table.col('Species (predicted)') !== null, 'No prediction column', WAIT_MS);
  expect(await readUntil(() => forgeDb.applications.query().where('model_id', '=', modelId).count(),
    (count) => count > 0), 1);
}

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
    expect(isDisabled(view.deleteIcon), true, 'Delete is enabled without a current row');
    expect(isDisabled(view.applyIcon), true, 'Apply is enabled without a current row');
    grok.shell.addView(view);
    let measured = '';
    try {
      expect(view.name, APP_NAME);
      const text = view.root.textContent ?? '';
      for (const caption of ['Methods', 'Models'])
        expect(text.includes(caption), true, `${caption} is missing`);
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
      for (const caption of ['Table', 'Target', 'Features', 'Method', 'Results'])
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
      expectArray([...view.methodInput.items].sort(), ['Linear Regression', 'PLS Regression', 'SVM', 'XGBoost']);
      expect(view.lastTraining?.engine.name, 'Linear Regression');
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: TIMEOUT});

  test('a change during a training supersedes it without a run row', async () => {
    const wine = await grok.data.files.openTable('System:DemoFiles/winequality.csv');
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
    const demog = grok.data.demo.demog(SLOW_ROWS);
    demog.name = `forge-test-demog-${Date.now()}`;
    const view = await TrainView.create(demog);
    grok.shell.addView(view);
    try {
      trainSlowly(view, demog);
      for (let i = 0; i < TIMEOUT / 10 && !view.isTraining; i++)
        await delay(10);
      expect(view.isTraining && view.methodInput.value === 'XGBoost', true, 'XGBoost does not train by itself');
      await awaitCheck(() => taskBarEntry('Training XGBoost model') !== undefined, 'No task-bar progress', TIMEOUT);
      // The entry's text becomes the fold progress after the first fit.
      const entry = taskBarEntry('Training XGBoost model');
      await delay(100);
      cancelTaskBarEntry(entry);
      await awaitCheck(() => !view.isTraining, 'The cancelled training does not end', TIMEOUT);
      const recorded = await runsOf([demog]);
      expectArray(recorded.map((r) => `${r.status} ${r.engine_name}`), ['cancelled XGBoost']);
      expect(view.lastTraining === undefined && isDisabled(view.saveButton), true, 'A cancelled training is saveable');

      view.methodInput.value = 'SVM';
      const hint = `SVM takes a while on ${SLOW_ROWS} rows, so it does not retrain on every change. Click Train.`;
      await awaitCheck(() => (view.root.textContent ?? '').includes(hint), 'No hint for a slow method', WAIT_MS);
      await delay(500);
      expect(view.isTraining, false, `SVM trains by itself on ${SLOW_ROWS} rows`);
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
      expectArray([...view.methodInput.items].sort(), ['SVM', 'Softmax', 'XGBoost']);
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
      const choice = view.missingValuesInputs.choice;
      expect(isShown(choice.root), true, 'Missing values is hidden with a gap');
      expect(choice.inputType, 'Radio');
      expect(choice.value, 'Skip rows');
      expectOptionsInRow(choice);
      // The skipped row makes every training column a copy: no frame holds a column of the table.
      const [, frames] = await framesSharing(iris.columns.toList(), () => view.train());
      expectNoFrames(frames);
      expect(Array.from(view.root.querySelectorAll('ul > li'), (li) => li.textContent)
        .includes('Rows: 149 used, 1 skipped (missing values)'), true, 'The Results rows item is missing');
      expect((await runs())[0]?.row_count, 149);

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
      expect(isTrainShown(view), false, 'Train is shown for XGBoost with Neighbors 0');
      neighbors.value = 3;
      await awaitTrained(view);
      expect(isTrainShown(view), false, 'Train is shown for XGBoost with Neighbors 3');

      view.featuresInput.value = columnsOf(iris, ['Sepal.Length', 'Petal.Length', 'Petal.Width']);
      await awaitCheck(() => !isShown(choice.root), 'Missing values stays without gaps', WAIT_MS);
    } finally {
      await closeTrainView(view, [iris]);
    }
  }, {timeout: TIMEOUT});
});

category('UI: Train', () => {
  test('a training superseded while its task-bar progress is shown writes no run', async () => {
    const demog = grok.data.demo.demog(SLOW_ROWS);
    demog.name = `forge-test-demog-${Date.now()}`;
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
    } finally {
      await closeTrainView(view, [demog]);
    }
  }, {timeout: TIMEOUT});

  test('an invalid hyperparameter blocks training and saving', async () => {
    const iris = await openIris();
    const view = await openTrainView(iris);
    try {
      setHyperparameter(view, 'iterations', 500);
      await awaitCheck(() => view.hyperparameterInputs.get('iterations')?.validity !== null,
        'Iterations is not marked', WAIT_MS);
      await delay(CHECK_DELAY_MS * 2);
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

  test('Preparation follows demog: Missing values, One-hot encoding, Skip unique categories', async () => {
    const demog = await openDemog();
    const view = await TrainView.create(demog);
    grok.shell.addView(view);
    const preparation = view.preparationInputs;
    const {oneHot, skipUniqueCategories, predictProbability, cutoff} = preparation;
    try {
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
      expect(isShown(columns.body), true, 'Columns stays collapsed with a problem');
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
    try {
      const dialog = (await applyModelDialog({table: iris, modelId: id})).show();
      opened = dialog;
      const ok = dialog.getButton('OK');
      await awaitCheck(() => !isDisabled(ok), 'OK stays disabled', WAIT_MS);
      ok.click();
      await expectApplied(dialog, iris, id);
    } finally {
      opened?.close();
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

  test('Apply dialog shows Missing values and the impute inputs', async () => {
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
      pressEnter(dialog, 'Model');
      await expectApplied(dialog, iris, id);
    } finally {
      opened?.close();
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
      const reload = async (description: string) => {
        const frame = view.models.dart;
        await updateModelInfo(otherId, undefined, {description});
        await awaitCheck(() => view.models.dart !== frame, 'The catalog does not reload', WAIT_MS);
        await awaitCheck(() => view.models.currentRowIdx === indexOf(id), 'The current row is lost', WAIT_MS);
      };
      view.models.currentRowIdx = indexOf(id);
      await awaitCheck(() => shows(id),
        'The chosen model is not the current object', WAIT_MS);
      const shown: unknown = grok.shell.o;
      const shownDart: unknown = shown instanceof DG.DomainRow ? shown.dart : null;
      await reload('refreshed');
      const now: unknown = grok.shell.o;
      expect(now instanceof DG.DomainRow && now.dart === shownDart, true,
        'The context panel got a new object for the same model');

      view.models.selection.set(indexOf(id), true);
      view.models.selection.set(indexOf(otherId), true);
      await awaitCheck(() => comparedIds() !== '', 'Two selected rows are not compared', WAIT_MS);
      await reload('refreshed again');
      expectArray(selectedIds(), [id, otherId].sort());
      expect(isDisabled(view.applyIcon) || isDisabled(view.compareIcon), false, 'An icon greyed out');
      await awaitCheck(() => comparedIds() === [id, otherId].sort().join(), 'The comparison is lost', WAIT_MS);
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
      let accordion = modelAccordion(row, model);
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
      accordion.root.remove();
      accordion = modelAccordion(row, model);
      host.append(accordion.root);
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
      await expectPaneText(accordion, 'Details', ['Last run', dayjs().format('YYYY-MM-DD')]);
      expect((accordion.getPane('Details').root.textContent ?? '').includes('Never'), false, 'Last run is Never');
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
      } finally {
        for (const dialog of newDialogs())
          dialog.close();
        accordion.root.remove();
      }
    } finally {
      opened?.close();
      await forgeDb.models.delete(id);
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
});

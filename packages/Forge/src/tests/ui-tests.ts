import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {awaitCheck, category, delay, expect, expectArray, test} from '@datagrok-libraries/test/src/test';
import {APP_NAME, MENU_PATH} from '../constants';
import {forgeDb} from '../generated/db';
import {deleteModel, modelsChanged} from '../storage/model-store';
import {applyModelDialog, modelLabels} from '../ui/apply-model-dialog';
import {ForgeApp} from '../ui/forge-app';
import {saveModelDialog} from '../ui/save-model-dialog';
import {TrainView} from '../ui/train-view';
import {columnsOf, expectReleased, framesSharing, MEASUREMENTS, openIris, saveIrisModel, saveTestModel, valuesOf,
  XGBOOST_FIELDS} from './test-data';

const WAIT_MS = 5000;
const TIMEOUT = 60000;
const ALIGN_PX = 4;

async function openTrainView(iris: DG.DataFrame): Promise<TrainView> {
  const view = await TrainView.create(iris);
  grok.shell.addView(view);
  await awaitCheck(() => !isDisabled(view.trainButton), 'Train stays disabled', WAIT_MS);
  return view;
}

const isDisabled = (button: HTMLElement) => button.classList.contains('d4-disabled');
const isShown = (element: HTMLElement) => element.style.display !== 'none';

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

/** Checks that the two options of the shown radio [input] sit on one line. */
function expectOptionsInRow(input: DG.InputBase): void {
  const options = Array.from(input.root.querySelectorAll<HTMLElement>('.ui-radio-button'));
  expect(options.length === 2 && options.every((o) => o.offsetHeight > 0 && o.offsetTop === options[0].offsetTop),
    true, `Option positions ${options.map((o) => o.offsetTop).join(', ')}`);
}

/** Waits until [dialog] is closed, [table] has the Species prediction and the application of [modelId] is recorded. */
async function expectApplied(dialog: DG.Dialog, table: DG.DataFrame, modelId: string): Promise<void> {
  const applications = () => forgeDb.applications.query().where('model_id', '=', modelId).count();
  await awaitCheck(() => !dialog.root.isConnected, 'The dialog stays open', WAIT_MS);
  await awaitCheck(() => table.col('Species (predicted)') !== null, 'No prediction column', WAIT_MS);
  for (let i = 0; i < 50 && await applications() === 0; i++)
    await delay(100);
  expect(await applications(), 1);
}

/** The texts of the options of the choice input captioned [caption] in [dialog], in list order. */
const optionsOf = (dialog: DG.Dialog, caption: string) =>
  Array.from(dialog.input(caption).root.querySelectorAll('option')).map((o) => o.textContent ?? '');

const pressEnter = (dialog: DG.Dialog, caption: string) => dialog.input(caption).input.dispatchEvent(
  new KeyboardEvent('keydown', {key: 'Enter', keyCode: 13, bubbles: true}));

category('UI', () => {
  test('Forge app opens', async () => {
    const view = await ForgeApp.create();
    expect(view.models.currentRowIdx, -1);
    expect(isDisabled(view.deleteIcon), true, 'Delete is enabled without a current row');
    expect(isDisabled(view.applyIcon), true, 'Apply is enabled without a current row');
    grok.shell.addView(view);
    try {
      expect(view.name, APP_NAME);
      const text = view.root.textContent ?? '';
      for (const caption of ['Methods', 'Method type', 'XGBoost', 'Models'])
        expect(text.includes(caption), true, `${caption} is missing`);
      expect(view.root.contains(view.deleteIcon), true, 'The delete icon is missing');
      expect(view.root.contains(view.applyIcon), true, 'The apply icon is missing');
      const grid = view.root.querySelector('.d4-grid');
      const header = view.root.querySelector('.forge-pane-header');
      if (!(grid instanceof HTMLElement) || !(header instanceof HTMLElement))
        throw new Error('No catalog grid or header');
      await awaitCheck(() => Math.abs(grid.getBoundingClientRect().width - header.getBoundingClientRect().width) < 1,
        'The catalog grid does not take the width of its pane', WAIT_MS);
      expect(view.models.col('id') !== null, true, 'The catalog has no id column');
      const captions = {name: 'Name', engine_name: 'Method', task: 'Task', target_name: 'Target',
        storage_mode: 'Data storage', row_count: 'Training rows', created_on: 'Created'};
      for (const [column, caption] of Object.entries(captions))
        expect(view.models.col(column)?.meta.friendlyName, caption, column);
    } finally {
      view.close();
    }
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

  test('Train view opens', async () => {
    const view = await TrainView.create(await openIris());
    grok.shell.addView(view);
    try {
      expect(view.name, 'Predictive model');
      expect(view.getIcon().classList.contains('svg-model'), true, 'The tab icon is not the model icon');
      const text = view.root.textContent ?? '';
      for (const caption of ['Table', 'Target', 'Features', 'Method', 'Train', 'Results'])
        expect(text.includes(caption), true, `${caption} is missing`);
      expect(text.includes('Engine'), false, 'The word Engine is shown');
      expect(view.getRibbonPanels().flat().some((e) => e.contains(view.saveButton)), true, 'Save is not in the ribbon');
      expect(isDisabled(view.saveButton), true, 'Save is enabled before training');
      expect(view.targetInput.value?.name, 'Species');
      expectArray(view.featuresInput.value.map((c) => c.name), MEASUREMENTS);
    } finally {
      view.close();
    }
  });

  test('Train view trains iris and saves the model once', async () => {
    const iris = await openIris();
    const runs = () => forgeDb.trainingRuns.query().where('dataset_name', '=', iris.name);
    const models = () => forgeDb.models.query().where('dataset_name', '=', iris.name);
    const view = await openTrainView(iris);
    try {
      const training = view.train();
      expect(isDisabled(view.trainButton), true, 'Train is enabled while training');
      expect(isDisabled(view.saveButton), true, 'Save is enabled while training');
      await Promise.all([training, view.train()]);
      expect((await runs()).length, 1, 'A second Train started another training');
      expect(isDisabled(view.trainButton), false, 'Train stays disabled after training');
      expect(isDisabled(view.saveButton), false, 'Save stays disabled after training');
      expect(view.lastTraining !== undefined, true, 'No training result');
      const text = view.root.textContent ?? '';
      for (const caption of ['Accuracy', 'F1'])
        expect(text.includes(caption), true, `${caption} is missing`);
      expect(text.includes('Sensitivity'), false, 'Sensitivity on three classes');
      expect(text.includes('Save...'), false, 'Save... is still in Results');

      const name = `forge-test-model-${Date.now()}`;
      await Promise.all([view.saveModelAs(name, ''), view.saveModelAs(`${name}-again`, '')]);
      expect(isDisabled(view.saveButton), true, 'Save stays enabled after saving');
      await view.saveModelAs(`${name}-later`, '');
      const saved = await models();
      expect(saved.length, 1, 'The same training was saved more than once');
      expect(saved[0].name, name);
      expect((await runs())[0].model_id, saved[0].id);
    } finally {
      for (const model of await models())
        await deleteModel(model.id);
      for (const run of await runs())
        await forgeDb.trainingRuns.delete(run.id);
      view.close();
    }
  }, {timeout: 60000});

  test('a failing training records a failed run', async () => {
    const iris = await openIris();
    const runs = () => forgeDb.trainingRuns.query().where('dataset_name', '=', iris.name);
    const view = await openTrainView(iris);
    try {
      const iterations = view.hyperparameterInputs.get('iterations');
      if (iterations === undefined)
        throw new Error('No iterations input');
      iterations.value = 0;
      await view.train();
      const recorded = await runs();
      expect(recorded.length, 1);
      expect(recorded[0].status, 'failed');
      expect((recorded[0].error ?? '') !== '', true, 'The failed run has no error');
      expect(view.lastTraining === undefined, true, 'A failed training left a result');
      expect(isDisabled(view.saveButton), true, 'Save is enabled after a failed training');
    } finally {
      for (const run of await runs())
        await forgeDb.trainingRuns.delete(run.id);
      view.close();
    }
  }, {timeout: 60000});

  test('a bad selection disables Train and marks the input', async () => {
    const iris = await openIris();
    const view = await TrainView.create(iris);
    grok.shell.addView(view);
    const columns = (names: string[]) => names.map((name) => iris.getCol(name));
    try {
      view.targetInput.value = iris.getCol('Petal.Length');
      view.featuresInput.value = columns(['Species', 'Sepal.Length', 'Sepal.Width', 'Petal.Width']);
      await awaitCheck(() => view.featuresInput.validity !== null, 'Features is not marked', WAIT_MS);
      expect(view.featuresInput.validity?.includes('Species'), true, `${view.featuresInput.validity}`);
      expect(isDisabled(view.trainButton), true, 'Train is enabled');
      expect(view.targetInput.validity === null, true, `Target is marked: ${view.targetInput.validity}`);

      view.featuresInput.value = columns(['Sepal.Length', 'Sepal.Width', 'Petal.Width']);
      await awaitCheck(() => view.featuresInput.validity === null && !isDisabled(view.trainButton),
        'Features stays marked or Train stays disabled', WAIT_MS);
    } finally {
      view.close();
    }
  }, {timeout: 30000});

  test('Train view groups its inputs and opens a group with a problem', async () => {
    const view = await openTrainView(await openIris());
    try {
      const data = groupOf(view.root, 'Data');
      const method = groupOf(view.root, 'Method');
      expect(isShown(data.body) && isShown(method.body), true, 'A group starts collapsed');
      expect(data.body.contains(view.featuresInput.root), true, 'Features is not under Data');
      expect(data.body.contains(view.missingValuesInputs.choice.root), true, 'Missing values is not under Data');
      const iterations = view.hyperparameterInputs.get('iterations');
      expect(iterations !== undefined && method.body.contains(iterations.root), true, 'Iterations is not under Method');
      expect(view.groups.every((g) => !g.root.contains(view.trainButton)), true, 'Train is inside a group');
      const trainRight = view.trainButton.getBoundingClientRect().right;
      const inputsRight = view.featuresInput.root.querySelector('.ui-input-editor')?.getBoundingClientRect().right;
      expect(inputsRight !== undefined && Math.abs(trainRight - inputsRight) < ALIGN_PX, true,
        `Train ends at ${trainRight}, the inputs at ${inputsRight}`);

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
      view.close();
    }
  }, {timeout: 30000});

  test('Train view drops the old form\'s subscriptions on a Table change', async () => {
    const tables = [grok.shell.addTable(await openIris()), grok.shell.addTable(await openIris())];
    const view = await TrainView.create(tables[0]);
    grok.shell.addView(view);
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
      view.close();
      for (const table of tables)
        grok.shell.closeTable(table);
    }
  }, {timeout: 30000});

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
    const runs = () => forgeDb.trainingRuns.query().where('dataset_name', '=', iris.name);
    const view = await openTrainView(iris);
    try {
      const choice = view.missingValuesInputs.choice;
      expect(isShown(choice.root), true, 'Missing values is hidden with a gap');
      expect(choice.inputType, 'Radio');
      expect(choice.value, 'Skip rows');
      expectOptionsInRow(choice);
      const [, frames] = await framesSharing(iris.columns.toList(), () => view.train());
      expectReleased(frames);
      expect((view.root.textContent ?? '').includes('Rows: 149 used, 1 skipped (missing values).'), true,
        'The Results rows line is missing');
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
      await awaitCheck(() => isDisabled(view.trainButton), 'Train is enabled with Neighbors 0', WAIT_MS);
      neighbors.value = 3;
      await awaitCheck(() => !isDisabled(view.trainButton), 'Train stays disabled with Neighbors 3', WAIT_MS);

      view.featuresInput.value = columnsOf(iris, ['Sepal.Length', 'Petal.Length', 'Petal.Width']);
      await awaitCheck(() => !isShown(choice.root), 'Missing values stays without gaps', WAIT_MS);
    } finally {
      for (const run of await runs())
        await forgeDb.trainingRuns.delete(run.id);
      view.close();
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

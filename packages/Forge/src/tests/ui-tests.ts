import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {awaitCheck, category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';
import {APP_NAME, MENU_PATH} from '../constants';
import {forgeDb} from '../generated/db';
import {deleteModel, modelsChanged} from '../storage/model-store';
import {ForgeApp} from '../ui/forge-app';
import {TrainView} from '../ui/train-view';
import {MEASUREMENTS, openIris, XGBOOST_FIELDS} from './test-data';

const WAIT_MS = 5000;

async function openTrainView(iris: DG.DataFrame): Promise<TrainView> {
  const view = await TrainView.create(iris);
  grok.shell.addView(view);
  await awaitCheck(() => !isDisabled(view.trainButton), 'Train stays disabled', WAIT_MS);
  return view;
}

const isDisabled = (button: HTMLElement) => button.classList.contains('d4-disabled');

category('UI', () => {
  test('Forge app opens', async () => {
    const view = await ForgeApp.create();
    expect(view.models.currentRowIdx, -1);
    expect(isDisabled(view.deleteIcon), true, 'Delete is enabled without a current row');
    grok.shell.addView(view);
    try {
      expect(view.name, APP_NAME);
      const text = view.root.textContent ?? '';
      for (const caption of ['Methods', 'Method type', 'XGBoost', 'Models'])
        expect(text.includes(caption), true, `${caption} is missing`);
      expect(view.root.contains(view.deleteIcon), true, 'The delete icon is missing');
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
      await awaitCheck(() => isDisabled(view.deleteIcon), 'Delete is enabled without a current row', WAIT_MS);
      view.models.currentRowIdx = rowOf();
      await awaitCheck(() => !isDisabled(view.deleteIcon), 'Delete is disabled with a current row', WAIT_MS);

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

  test('train function is registered', async () => {
    const topMenu = DG.Func.find({package: 'Forge', name: 'forgeTrain'})[0]?.topMenu;
    expect(topMenu?.startsWith(MENU_PATH), true, `Unexpected top menu: ${topMenu}`);
  });
});

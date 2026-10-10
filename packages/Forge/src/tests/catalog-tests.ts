import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {after, awaitCheck, before, category, expect, expectArray, expectExceptionAsync, expectFloat, test}
  from '@datagrok-libraries/test/src/test';
import {applyAndRecord, ApplyRequest, DEFAULT_BATCH_SIZE, loadModel} from '../apply/apply-model';
import {exactMapping} from '../apply/column-matching';
import {applicableTables} from '../catalog/applicable-tables';
import {compareFormFields, compareModels, CompareModelRow} from '../catalog/compare-models';
import {modelActivity} from '../catalog/model-activity';
import {updateModelInfo} from '../catalog/model-edit';
import {readModelFile} from '../catalog/model-file';
import {PANEL_PREDICTED_BY, PREDICTION_TAG} from '../constants';
import {ForgeError} from '../forge-error';
import {forgeDb} from '../generated/db';
import {tagsOf} from '../storage/model-fields';
import {deleteModel, modelsChanged} from '../storage/model-store';
import {comparisonForms, hasFormsViewer, ModelComparison, NO_FORMS} from '../ui/model-comparison';
import {ForgeModelHandler, forgeModelHandler} from '../ui/model-handler';
import {isForgePrediction, predictedByPane} from '../ui/prediction-column-panel';
import {trainModel} from '../training/train-model';
import {columnsOf, insertModelRow, MEASUREMENTS, openIris, PREDICT_PROBABILITY, requestOf, savedFixture,
  saveIrisModel, twoSpeciesIris} from './test-data';

const TIMEOUT = 90000;
const WAIT_MS = 10000;
const IRIS_ROW = {features: {columns: MEASUREMENTS.map((name) => ({name, type: DG.COLUMN_TYPE.FLOAT}))}};

interface Fixture { id: string; blob: Uint8Array; iris: DG.DataFrame }

let fixture: Fixture | undefined;
const shared = () => savedFixture(fixture);

/** A regressor and a classifier as a comparison reads them. */
function compareRows(): CompareModelRow[] {
  const created = dayjs('2026-10-06T10:00:00Z');
  return [{id: 'r', name: 'forge-test-regressor', description: 'Petal length', created_on: created,
    engine_name: 'XGBoost', task: 'regression', target_name: 'Petal.Length', row_count: 150,
    metrics: {train: {mse: 0.1, rmse: 0.3, mae: 0.2, r2: 0.9}, validation: {mse: 0.2, rmse: 0.4, r2: 0.8}}},
  {id: 'c', name: 'forge-test-classifier', created_on: created, engine_name: 'XGBoost',
    task: 'classification', target_name: 'Species', row_count: 140,
    metrics: {train: {accuracy: 1, f1: 1}, validation: {accuracy: 0.95, f1: 0.94}, positiveClass: 'setosa'}}];
}

category('Catalog', () => {
  before(async () => {
    const iris = await openIris();
    const {id, result} = await saveIrisModel(iris);
    fixture = {id, blob: result.blob, iris};
  });

  after(async () => {
    if (fixture !== undefined)
      await deleteModel(fixture.id);
  });

  test('compareModels lists the models and their metrics side by side', async () => {
    const [regressor, classifier] = compareRows();
    const df = compareModels([regressor, classifier]);
    expect(df.name, 'Compare models');
    expect(df.rowCount, 2);
    expectArray(df.columns.names(), ['Name', 'Description', 'Method', 'Task', 'Target', 'Training rows', 'Created',
      'MSE (train)', 'MSE (validation)', 'RMSE (train)', 'RMSE (validation)', 'MAE (train)', 'MAE (validation)',
      'R2 (train)', 'R2 (validation)', 'Accuracy (train)', 'Accuracy (validation)', 'F1 (train)',
      'F1 (validation)']);
    expectArray(df.columns.toList().slice(5).map((c) => c.type), [DG.COLUMN_TYPE.INT, DG.COLUMN_TYPE.DATE_TIME,
      ...Array<string>(12).fill(DG.COLUMN_TYPE.FLOAT)]);
    expectArray(df.getCol('Name').toList(), ['forge-test-regressor', 'forge-test-classifier']);
    expect(df.get('Training rows', 1), 140);
    expectFloat(df.get('MSE (train)', 0), 0.1, 1e-6);
    expect(df.getCol('MSE (train)').isNone(1), true, 'The classifier has an MSE');
    expect(df.getCol('MAE (validation)').isNone(0), true, 'The regressor has a validation MAE');
    expectFloat(df.get('Accuracy (validation)', 1), 0.95, 1e-6);
    expect(df.getCol('Accuracy (validation)').isNone(0), true, 'The regressor has an accuracy');
    expectArray(compareFormFields(df), ['Name', 'Method', 'Task', 'Target', 'Training rows', 'Created',
      'MSE (validation)', 'RMSE (validation)', 'MAE (validation)', 'R2 (validation)', 'Accuracy (validation)',
      'F1 (validation)']);

    await expectExceptionAsync(async () => {
      compareModels([regressor]);
    }, (e) => e instanceof ForgeError && e.message === 'Select at least two models to compare.');
  });

  test('compareModels lists AUC-ROC (train) and (validation) for a predict-probability model', async () => {
    const iris = await twoSpeciesIris();
    const result = await trainModel(await requestOf(columnsOf(iris, MEASUREMENTS), iris.getCol('Species'),
      PREDICT_PROBABILITY));
    const [regressor] = compareRows();
    const probability: CompareModelRow = {id: 'p', name: 'forge-test-probability', created_on: dayjs(),
      engine_name: 'XGBoost', task: result.task, target_name: 'Species', row_count: result.rowCount,
      metrics: result.metrics};
    const df = compareModels([regressor, probability]);
    const names = df.columns.names();
    expect(names.includes('AUC-ROC (train)') && names.includes('AUC-ROC (validation)'), true, names.join(', '));
    expectFloat(df.get('AUC-ROC (validation)', 1), result.metrics.validation.auc ?? NaN, 1e-6);
    expect(df.getCol('AUC-ROC (train)').isNone(0), true, 'The regressor has an AUC-ROC');
  }, {timeout: TIMEOUT});

  test('the comparison handler claims a comparison and shows it as forms', async () => {
    const comparison = new ModelComparison(compareRows());
    expect(DG.ObjectHandler.forEntity(comparison)?.name, 'Forge model comparison handler');
    const root = comparisonForms(comparison.rows);
    const title = root.firstElementChild;
    expect(title?.classList.contains('d4-accordion-title') && title.textContent, 'Compare 2 models');
    expect(title?.querySelector('.svg-model') !== null, true, 'The title has no model icon');
    // At the context panel's width: the viewer lays its forms out once attached.
    const host = ui.div([root]);
    host.style.width = '320px';
    document.body.append(host);
    try {
      const labels = () => Array.from(root.querySelectorAll('.d4-multi-form-column-name'), (e) => e.textContent);
      if (!hasFormsViewer())
        expect(root.textContent?.includes(NO_FORMS), true, root.textContent ?? '');
      else {
        await awaitCheck(() => labels().includes('Name') && labels().includes('F1 (validation)'),
          `No forms with Name and F1 (validation): ${root.outerHTML.slice(0, 500)}`, WAIT_MS);
        const widthOf = (e: Element | null) => Math.round(e?.getBoundingClientRect().width ?? -1);
        const widths = [root.children[1], root.querySelector('.d4-multi-form')].map(widthOf);
        expectArray(widths, [widthOf(host), widthOf(host)]);
      }
    } finally {
      host.remove();
    }
  });

  test('applicableTables keeps the tables with a close column for every feature', async () => {
    const cars = await grok.data.files.openTable('System:DemoFiles/cars.csv');
    const {iris} = shared();
    expectArray(applicableTables(IRIS_ROW, [cars, iris]).map((t) => t.name), [iris.name]);
    expect(applicableTables(IRIS_ROW, [cars]).length, 0);
    expect(applicableTables({}, [iris]).length, 0);
  });

  test('applicableTables does not require a feature that Skip unique categories left out', async () => {
    const features = {columns: [{name: 'x', type: DG.COLUMN_TYPE.FLOAT}, {name: 'y', type: DG.COLUMN_TYPE.FLOAT},
      {name: 'subject', type: DG.COLUMN_TYPE.STRING}]};
    const table = DG.DataFrame.fromColumns([DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'x', [1, 2]),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'y', [3, 4])]);
    expect(applicableTables({features, options: {skippedColumns: ['subject']}}, [table]).length, 1);
    expect(applicableTables({features}, [table]).length, 0);
  });

  test('modelActivity lists the applications newest first and the last run of any status', async () => {
    const {id, iris} = shared();
    const empty = await modelActivity(id);
    expect(empty.count, 0);
    expect(empty.applications.length, 0);
    expect(empty.lastRun === undefined, true, 'A model never applied has a last run');

    const model = await loadModel(id);
    const request: ApplyRequest = {model, table: iris, mapping: exactMapping(model.features, iris),
      batchSize: DEFAULT_BATCH_SIZE, missingValues: {mode: 'skip'}};
    await applyAndRecord(request, 'api');
    await expectExceptionAsync(async () => {
      await applyAndRecord(request, 'ui', {canceled: true, update: () => {}});
    }, (e) => e instanceof ForgeError && e.message === 'Application was cancelled.');

    const activity = await modelActivity(id);
    expect(activity.count, 2);
    expectArray(activity.applications.map((a) => a.status), ['cancelled', 'completed']);
    expect(activity.lastRun?.status, 'cancelled');
    expect(activity.lastRun?.when.isSame(activity.applications[0].created_on), true, 'Last run is not the newest');
  }, {timeout: TIMEOUT});

  test('readModelFile returns the saved blob and refuses other files', async () => {
    const {id, blob} = shared();
    const row = await forgeDb.models.get(id);
    const file = await readModelFile(row);
    expect(file.name, `${row.name}.bin`);
    expectArray(Array.from(file.bytes), Array.from(blob));
    expect((await readModelFile({name: 'forge-test: <a/b>', blob: row.blob})).name, 'forge-test_ _a_b_.bin');

    const refusal = (name: string) => (e: unknown) => e instanceof ForgeError &&
      e.message === `The model '${name}' has no model file in Forge's storage, so there is nothing to download.`;
    await expectExceptionAsync(() => readModelFile({name: 'forge-test-elsewhere',
      blob: 'file://System:DemoFiles/iris.csv'}).then(() => {}), refusal('forge-test-elsewhere'));
    await expectExceptionAsync(() => readModelFile({name: 'forge-test-no-file'}).then(() => {}),
      refusal('forge-test-no-file'));
  }, {timeout: TIMEOUT});

  test('the model handler claims model rows and renders the card and the tooltip', async () => {
    const model = await forgeDb.models.get(shared().id);
    const row = ForgeModelHandler.rowOf(model);
    expect(forgeModelHandler.isApplicable(row), true);
    expect(forgeModelHandler.isApplicable('forge-test-string'), false);
    expect(forgeModelHandler.getCaption(row), model.name);
    const card = forgeModelHandler.renderCard(row).textContent ?? '';
    for (const text of ['Predict Species', 'by Sepal.Length, Sepal.Width, Petal.Length, Petal.Width', 'using XGBoost',
      `Created on ${model.created_on.format('YYYY-MM-DD')}`])
      expect(card.includes(text), true, `The card misses '${text}': ${card}`);
    const tooltip = forgeModelHandler.renderTooltip(row).textContent ?? '';
    for (const text of [model.name, 'XGBoost', 'classification', 'Species', 'Training rows', '150'])
      expect(tooltip.includes(text), true, `The tooltip misses '${text}': ${tooltip}`);
  });

  test('isPredictionColumn and the Predicted by panel', async () => {
    const col = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'p', ['a']);
    expect(isForgePrediction(col), false);
    expect(await grok.functions.call('Forge:isPredictionColumn', {col}), false);
    col.setTag(PREDICTION_TAG, shared().id);
    expect(isForgePrediction(col), true);
    expect(await grok.functions.call('Forge:isPredictionColumn', {col}), true);
    // The server names a function by its code name and keeps the annotation's name as the friendly name.
    const panel = DG.Func.find({package: 'Forge', name: 'predictedByPanel'})[0];
    expect(panel?.friendlyName, PANEL_PREDICTED_BY);
    expect(panel?.options['condition'], 'Forge:isPredictionColumn(col)');
    expect(panel?.options['role'], 'panel');

    const pane = predictedByPane(col);
    await awaitCheck(() => (pane.textContent ?? '').includes('Predict Species'), 'The card is not shown', WAIT_MS);
    col.setTag(PREDICTION_TAG, DG.Utils.uuid4());
    const gone = predictedByPane(col);
    const noModel = 'The model that predicted this column is no longer available to you.';
    await awaitCheck(() => (gone.textContent ?? '') === noModel, 'A missing model is not explained', WAIT_MS);
  });

  test('updateModelInfo writes name and description and refuses a stale version', async () => {
    const stamp = Date.now();
    const id = await insertModelRow({name: `forge-test-model-${stamp}`});
    let changes = 0;
    const sub = modelsChanged.subscribe(() => changes++);
    try {
      const {version} = await forgeDb.models.get(id);
      const edited = await updateModelInfo(id, version, {name: `forge-test-edited-${stamp}`, description: 'Edited'});
      let row = await forgeDb.models.get(id);
      expect(row.version, edited);
      expect(row.name, `forge-test-edited-${stamp}`);
      expect(row.description, 'Edited');
      expect(changes, 1);

      await updateModelInfo(id, edited, {description: 'Edited again'});
      row = await forgeDb.models.get(id);
      expect(row.name, `forge-test-edited-${stamp}`);
      expect(row.description, 'Edited again');

      await expectExceptionAsync(() => updateModelInfo(id, edited, {description: 'Stale'}).then(() => {}),
        (e) => e instanceof DG.DomainVersionConflictError);
      expect((await forgeDb.models.get(id)).description, 'Edited again');
      await updateModelInfo(id, undefined, {description: ''});
      expect((await forgeDb.models.get(id)).description ?? '', '');
    } finally {
      sub.unsubscribe();
      await forgeDb.models.delete(id);
    }
  });

  test('updateModelInfo writes the tags as their text, trimmed, without empty ones and repeats', async () => {
    const id = await insertModelRow({name: `forge-test-model-${Date.now()}`});
    try {
      const {version} = await forgeDb.models.get(id);
      const edited = await updateModelInfo(id, version, {tags: [' a ', 'b', '', 'a']});
      const {tags} = await forgeDb.models.get(id);
      expect(tags, 'a, b');
      expectArray(tagsOf(tags), ['a', 'b']);
      await updateModelInfo(id, edited, {tags: []});
      const cleared = (await forgeDb.models.get(id)).tags;
      expect(cleared == null, true, `Cleared tags are stored as '${cleared}'`);
    } finally {
      await forgeDb.models.delete(id);
    }
  });
});

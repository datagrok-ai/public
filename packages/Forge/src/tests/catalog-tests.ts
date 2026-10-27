import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {after, awaitCheck, before, category, expect, expectArray, expectExceptionAsync, expectFloat, test}
  from '@datagrok-libraries/test/src/test';
import {applyAndRecord, ApplyRequest, DEFAULT_BATCH_SIZE, loadModel} from '../apply/apply-model';
import {exactMapping} from '../apply/column-matching';
import {applicableTables, modelsFor} from '../catalog/applicable-tables';
import {compareFormFields, compareModels, CompareModelRow} from '../catalog/compare-models';
import {runGetModel} from '../catalog/get-model';
import {ListModelsInput, runListModels} from '../catalog/list-models';
import {modelActivity} from '../catalog/model-activity';
import {runUpdateModel, updateModelInfo} from '../catalog/model-edit';
import {readModelFile} from '../catalog/model-file';
import {PANEL_PREDICTED_BY, PREDICTION_TAG} from '../constants';
import {forgeDb} from '../generated/db';
import {metricFindings} from '../metrics/metrics';
import {NO_NAME, tagsOf, tagsText} from '../storage/model-fields';
import {deleteModel, modelsChanged} from '../storage/model-store';
import {ForgeApp} from '../ui/forge-app';
import {comparisonForms, hasFormsViewer, ModelComparison, NO_FORMS} from '../ui/model-comparison';
import {ForgeModelHandler, forgeModelHandler} from '../ui/model-handler';
import {isForgePrediction, predictedByPane} from '../ui/prediction-column-panel';
import {runTrainModel} from '../training/train-api';
import {metricsRecordOf, trainModel} from '../training/train-model';
import {columnsOf, deleteModelsNamed, deleteRunsOf, expectMetrics, failsWith, insertModelRow, MEASUREMENTS,
  openDemoTable, openIris, PREDICT_PROBABILITY, requestOf, savedFixture, saveIrisModel, twoSpeciesIris, WINE}
  from './test-data';

const TIMEOUT = 90000;
const WAIT_MS = 10000;
const IRIS_ROW = {features: {columns: MEASUREMENTS.map((name) => ({name, type: DG.COLUMN_TYPE.FLOAT}))}};

const API_TAGS = ['forge-test', 'api'];

interface Fixture { id: string; blob: Uint8Array; iris: DG.DataFrame }
interface QualityFixture { id: string; name: string; table: DG.DataFrame }

let fixture: Fixture | undefined;
let quality: QualityFixture | undefined;
let qualityName: string | undefined;
let qualityTable: DG.DataFrame | undefined;
const shared = () => savedFixture(fixture);
const sharedQuality = () => savedFixture(quality);

async function listedIds(input: ListModelsInput): Promise<string[]> {
  return (await runListModels(input)).getCol('Id').toList();
}

/** Whether the listed [ids] hold the iris fixture and the quality fixture. */
function holds(ids: string[]): [boolean, boolean] {
  return [ids.includes(shared().id), ids.includes(sharedQuality().id)];
}

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
    const {id, result} = await saveIrisModel(iris, {tags: tagsText(API_TAGS)});
    fixture = {id, blob: result.blob, iris};
    const table = await openDemoTable(WINE, 'cat-quality');
    qualityTable = table;
    const name = `forge-test-cat-quality-${Date.now()}`;
    qualityName = name;
    const {modelId} = await runTrainModel({table, target: 'quality', name, showProgress: false});
    if (modelId !== null)
      quality = {id: modelId, name, table};
  });

  after(async () => {
    if (fixture !== undefined)
      await deleteModel(fixture.id);
    if (qualityName !== undefined)
      await deleteModelsNamed(qualityName);
    if (qualityTable !== undefined)
      await deleteRunsOf([qualityTable]);
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
    }, failsWith('Select at least two models to compare.'));
  });

  test('compareModels lists AUC-ROC for a predict-probability model and a plain classifier', async () => {
    const iris = await twoSpeciesIris();
    const result = await trainModel(await requestOf(columnsOf(iris, MEASUREMENTS), iris.getCol('Species'),
      PREDICT_PROBABILITY));
    const [regressor] = compareRows();
    const probability: CompareModelRow = {id: 'p', name: 'forge-test-probability', created_on: dayjs(),
      engine_name: 'XGBoost', task: result.task, target_name: 'Species', row_count: result.rowCount,
      metrics: result.metrics};
    const threeSpecies = await openIris();
    const classifier = await trainModel(await requestOf(columnsOf(threeSpecies, MEASUREMENTS),
      threeSpecies.getCol('Species')));
    const plain: CompareModelRow = {...probability, id: 'c', name: 'forge-test-classifier', task: classifier.task,
      row_count: classifier.rowCount, metrics: classifier.metrics};
    const df = compareModels([regressor, probability, plain]);
    const names = df.columns.names();
    expect(names.includes('AUC-ROC (train)') && names.includes('AUC-ROC (validation)'), true, names.join(', '));
    expectFloat(df.get('AUC-ROC (validation)', 1), result.metrics.validation.auc ?? NaN, 1e-6);
    expectFloat(df.get('AUC-ROC (validation)', 2), classifier.metrics.validation.auc ?? NaN, 1e-6);
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

  test('modelsFor keeps the rows whose models fit the table, in their order', async () => {
    const {iris} = shared();
    const xy = {columns: [{name: 'x', type: DG.COLUMN_TYPE.FLOAT}, {name: 'y', type: DG.COLUMN_TYPE.FLOAT}]};
    const rows = [{id: 'xy', features: xy}, {id: 'iris', ...IRIS_ROW}, {id: 'iris-2', ...IRIS_ROW}];
    expectArray(modelsFor(rows, iris).map((row) => row.id), ['iris', 'iris-2']);
    expect(modelsFor([], iris).length, 0);
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
    }, failsWith('Application was cancelled.'));

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

    const refusal = (name: string) =>
      failsWith(`The model '${name}' has no model file in Forge's storage, so there is nothing to download.`);
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
    const card = forgeModelHandler.renderCard(row);
    expect(card.children.length, 5);
    for (const text of ['Predict Species', 'by Sepal.Length, Sepal.Width, Petal.Length, Petal.Width', 'using XGBoost',
      `Created on ${model.created_on.format('YYYY-MM-DD')}`])
      expect((card.textContent ?? '').includes(text), true, `The card misses '${text}': ${card.textContent}`);
    const applicable = card.querySelector('.forge-model-card-title > .grok-icon');
    expect(applicable !== null, applicableTables(IRIS_ROW, grok.shell.tables).length > 0,
      'The Applicable to icon does not follow the open tables');
    const tooltip = forgeModelHandler.renderTooltip(row);
    const text = tooltip.textContent ?? '';
    for (const part of [model.name, 'Predicts Species from 4 features with XGBoost · classification · 150 rows',
      'Sepal.Length, Sepal.Width, Petal.Length, Petal.Width',
      `Created ${model.created_on.format('YYYY-MM-DD')}`, 'Data storage: None'])
      expect(text.includes(part), true, `The tooltip misses '${part}': ${text}`);
    expect(tooltip.querySelector('.forge-tooltip-title')?.textContent, model.name);
    await awaitCheck(() => (tooltip.textContent ?? '').includes(`by ${DG.User.current().login}`),
      `The tooltip does not name the author: ${tooltip.textContent}`, WAIT_MS);
  });

  test('a row without its json columns, as the Domains view gives it, gets its features read once', async () => {
    const model = await forgeDb.models.get(shared().id);
    const row = ForgeModelHandler.rowOf({...model, features: undefined, options: undefined});
    expect(row.values.features == null && row.values.options == null, true, 'The row carries its json columns');
    const card = forgeModelHandler.renderCard(row);
    const byLine = () => card.children[2]?.textContent ?? '';
    expect(byLine(), 'by 4 features');
    const tooltip = forgeModelHandler.renderTooltip(row);
    const tooltipText = () => tooltip.textContent ?? '';
    expect(tooltipText().includes('Predicts Species from 4 features with XGBoost'), true, tooltipText());
    const names = 'Sepal.Length, Sepal.Width, Petal.Length, Petal.Width';
    await awaitCheck(() => byLine() === `by ${names}`, `The by line stays '${byLine()}'`, WAIT_MS);
    await awaitCheck(() => tooltipText().includes(names), `The tooltip lists no features: ${tooltipText()}`, WAIT_MS);
    expect(tooltipText().includes('0 features'), false, tooltipText());
    expect(card.children.length, 5);
    const applicable = card.querySelector('.forge-model-card-title > .grok-icon');
    expect(applicable !== null, applicableTables(IRIS_ROW, grok.shell.tables).length > 0,
      'The Applicable to icon does not follow the open tables');
    const again = forgeModelHandler.renderCard(row);
    await awaitCheck(() => again.children[2]?.textContent === `by ${names}`, 'A second card is not filled', WAIT_MS);
  });

  test('a catalog row\'s Name cell has the model\'s tooltip; other cells keep the grid\'s', async () => {
    const {id} = shared();
    const view = await ForgeApp.create();
    grok.shell.addView(view);
    try {
      const grid = Array.from(view.root.querySelectorAll('.d4-grid'), (e) => DG.toJs(DG.Widget.find(e)))
        .find((g): g is DG.Grid => g instanceof DG.Grid && g.dataFrame.col('id') !== null);
      if (grid === undefined)
        throw new Error('No catalog grid');
      grid.sort(['name'], [false]);
      const row = grid.tableRowToGrid(view.models.getCol('id').toList().indexOf(id));
      const name = view.modelTooltip(grid.cell('name', row))?.textContent ?? '';
      expect(name.includes('Predicts Species from 4 features with XGBoost'), true, `Name cell tooltip: ${name}`);
      expect(view.modelTooltip(grid.cell('engine_name', row)) === null, true, 'Another cell has the model tooltip');
    } finally {
      view.close();
    }
  });

  test('the tooltip and the card of a row with only an id and a name do not throw', async () => {
    const row = ForgeModelHandler.rowOf({id: DG.Utils.uuid4(), name: 'forge-test-bare'});
    const tooltip = forgeModelHandler.renderTooltip(row);
    expect(tooltip.querySelector('.forge-tooltip-title')?.textContent, 'forge-test-bare');
    expect((tooltip.textContent ?? '').includes('undefined'), false, tooltip.textContent ?? '');
    expect(forgeModelHandler.renderCard(row).children.length, 5);
  });

  test('model cards have one height whatever the model holds; the tooltip cuts a long feature list', async () => {
    const features = (n: number) => ({columns: Array.from({length: n},
      (_, i) => ({name: `forge-test-feature-${i + 1}`, type: DG.COLUMN_TYPE.FLOAT}))});
    const rows = [{name: 'forge-test-short', features: features(1)},
      {name: `forge-test-${'long-name-'.repeat(8)}`, features: features(12)}]
      .map((values) => ForgeModelHandler.rowOf({id: DG.Utils.uuid4(), target_name: 'y', engine_name: 'XGBoost',
        options: {}, ...values}));
    const cards = rows.map((row) => forgeModelHandler.renderCard(row));
    // As the Domains view shows a card: ui.card's class in the gallery's wrapper (xamgle grok_api.dart:1238).
    for (const card of cards)
      card.classList.add('d4-gallery-card');
    const host = ui.div(cards.map((card) => ui.div([card], 'grok-gallery-grid-item-wrapper')));
    document.body.append(host);
    try {
      const heights = cards.map((card) => card.offsetHeight);
      expect(heights[0] > 0 && heights[0] === heights[1], true, `Card heights ${heights.join(', ')}`);
      for (const card of cards) {
        expect(card.scrollHeight <= card.clientHeight, true,
          `The card's lines do not fit: ${card.scrollHeight} > ${card.clientHeight}`);
        const last = card.lastElementChild;
        expect(last instanceof HTMLElement && last.getBoundingClientRect().bottom <=
          card.getBoundingClientRect().bottom, true, 'The last line is cut');
      }
    } finally {
      host.remove();
    }
    const tooltip = forgeModelHandler.renderTooltip(rows[1]).textContent ?? '';
    expect(/\+\d+ more/.test(tooltip) && !tooltip.includes('forge-test-feature-12'), true, tooltip);
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

// The API functions, on the fixtures of the category's before.
category('Catalog', () => {
  test('listModels lists the models newest first and filters by task, method, name and target', async () => {
    const all = await listedIds({});
    expectArray(holds(all), [true, true]);
    expect(all.indexOf(sharedQuality().id) < all.indexOf(shared().id), true, 'The newer quality model is not first');
    expectArray(holds(await listedIds({task: 'classification'})), [true, false]);
    expect(holds(await listedIds({method: 'XGBoost'}))[0], true, 'XGBoost does not list the iris model');
    expectArray(holds(await listedIds({name: 'CAT-QUALITY'})), [false, true]);
    expectArray(holds(await listedIds({target: 'species'})), [true, false]);
  }, {timeout: TIMEOUT});

  test('listModels filters by tags, every one of them, in any case', async () => {
    expectArray(holds(await listedIds({tags: 'API'})), [true, false]);
    expectArray(holds(await listedIds({tags: 'api, zzz'})), [false, false]);
  });

  test('listModels filters by the table the model fits and limits the rows', async () => {
    expectArray(holds(await listedIds({table: shared().iris})), [true, false]);
    expectArray(holds(await listedIds({table: sharedQuality().table})), [false, true]);
    expect((await runListModels({limit: 1})).rowCount, 1);
  });

  test('listModels refuses an unknown task and a limit below 1', async () => {
    await expectExceptionAsync(() => runListModels({task: 'x'}).then(() => {}),
      failsWith('The task must be regression or classification.'));
    await expectExceptionAsync(() => runListModels({limit: 0}).then(() => {}),
      failsWith('The limit must be a whole number of at least 1.'));
  });

  test('listModels gives the columns of a model and its validation metrics', async () => {
    const {id} = shared();
    const df = await runListModels({tags: 'api'});
    expect(df.name, 'Forge models');
    expectArray(df.columns.names().slice(0, 12), ['Id', 'Name', 'Description', 'Method', 'Task', 'Target',
      'Features', 'Training rows', 'Data storage', 'Tags', 'Author', 'Created']);
    expect(df.columns.names().some((name) => name.endsWith(' (train)')), false, df.columns.names().join(', '));
    const row = df.getCol('Id').toList().indexOf(id);
    const stored = metricsRecordOf((await forgeDb.models.get(id)).metrics);
    expectFloat(df.get('Accuracy (validation)', row), stored?.validation.accuracy ?? NaN, 1e-6);
    expect(df.get('Features', row), MEASUREMENTS.join(', '));
    expect(df.get('Author', row), DG.User.current().login);
    expect(df.get('Tags', row), 'forge-test, api');
    expect(df.get('Method', row), 'XGBoost');
    expect(df.get('Training rows', row), 150);
    expect(df.get('Data storage', row), 'none');
  });

  test('getModel by id and by name gives everything Forge keeps about the model', async () => {
    const {id, iris} = shared();
    const row = await forgeDb.models.get(id);
    const stored = metricsRecordOf(row.metrics);
    if (stored === null)
      throw new Error('The fixture has no metrics');
    for (const idOrName of [id, row.name]) {
      const details = await runGetModel(idOrName, [iris]);
      expect(details.id, id);
      expect(details.version, row.version);
      expect(details.name, row.name);
      expectArray(details.tags, API_TAGS);
      expect(details.method, 'XGBoost');
      expect(details.package, 'Eda');
      expect(details.task, 'classification');
      expect(details.target, 'Species');
      expect(details.targetType, DG.COLUMN_TYPE.STRING);
      expectArray([...details.classes ?? []].sort(), ['setosa', 'versicolor', 'virginica']);
      expectArray(details.features, MEASUREMENTS);
      expectArray(MEASUREMENTS.map((name) => details.featureTypes[name]), MEASUREMENTS.map(() => DG.COLUMN_TYPE.FLOAT));
      expectMetrics(details.metrics?.train ?? {}, stored.train, 1e-9);
      expectMetrics(details.metrics?.validation ?? {}, stored.validation, 1e-9);
      expectArray(details.findings.map((f) => f.text), metricFindings(stored).map((f) => f.text));
      expect(details.storage.mode, 'none');
      expect(details.storage.source === undefined && details.storage.tableId === undefined, true,
        JSON.stringify(details.storage));
      expect(details.dataset.rows, 150);
      expect(details.hasModelFile, true);
      expect(details.keepsTrainingRows, false);
      expect(details.author, DG.User.current().login);
      expect(dayjs(details.created).isSame(row.created_on) && dayjs(details.updated).isSame(row.updated_on), true,
        `${details.created}, ${details.updated}`);
      expect(details.applicableTo.includes(iris.name), true, details.applicableTo.join(', '));
    }
    expect((await runGetModel(id, [])).applicableTo.length, 0);
  });

  test('getModel counts the applications and gives the last run', async () => {
    const {id, table} = sharedQuality();
    const before = await runGetModel(id, []);
    expect(before.lastRun === null, true, JSON.stringify(before.lastRun));
    expect(before.applications, 0);
    const model = await loadModel(id);
    const applied = table.clone();
    await applyAndRecord({model, table: applied, mapping: exactMapping(model.features, applied),
      batchSize: DEFAULT_BATCH_SIZE, missingValues: {mode: 'skip'}}, 'api');
    const after = await runGetModel(id, []);
    expect(after.applications, 1);
    expect(after.lastRun?.status, 'completed');
    expect(typeof after.lastRun?.when === 'string' && dayjs(after.lastRun.when).isValid(), true);
  }, {timeout: TIMEOUT});

  test('getModel of a model without a file, and of a missing model', async () => {
    const id = await insertModelRow({name: `forge-test-model-${Date.now()}`});
    try {
      const details = await runGetModel(id, []);
      expect(details.hasModelFile, false);
      expect(details.metrics === null, true, JSON.stringify(details.metrics));
      expect(details.findings.length, 0);
    } finally {
      await forgeDb.models.delete(id);
    }
    await expectExceptionAsync(() => runGetModel('forge-test-missing', []).then(() => {}),
      failsWith('No Forge model \'forge-test-missing\' is available to you.'));
  });

  test('updateModel renames, describes and tags a model, the rest staying', async () => {
    const stamp = Date.now();
    const id = await insertModelRow({name: `forge-test-model-${stamp}`, description: 'Before'});
    let changes = 0;
    const sub = modelsChanged.subscribe(() => changes++);
    try {
      const {version} = await forgeDb.models.get(id);
      const renamed = await runUpdateModel({model: id, name: `forge-test-renamed-${stamp}`});
      expect(renamed.id, id);
      expect(renamed.name, `forge-test-renamed-${stamp}`);
      expect(renamed.description, 'Before');
      expect(renamed.version > version, true, `${renamed.version} after ${version}`);
      expect((await forgeDb.models.get(id)).name, `forge-test-renamed-${stamp}`);
      expect(changes, 1);

      const described = await runUpdateModel({model: `forge-test-renamed-${stamp}`, description: 'After'});
      expect(described.name, `forge-test-renamed-${stamp}`);
      const row = await forgeDb.models.get(id);
      expect(row.name, `forge-test-renamed-${stamp}`);
      expect(row.description, 'After');

      expectArray((await runUpdateModel({model: id, tags: 'x, y,x'})).tags, ['x', 'y']);
      expect((await forgeDb.models.get(id)).tags, 'x, y');
      expectArray((await runUpdateModel({model: id, tags: ''})).tags, []);
      const cleared = (await forgeDb.models.get(id)).tags;
      expect(cleared == null, true, `Cleared tags are stored as '${cleared}'`);
    } finally {
      sub.unsubscribe();
      await forgeDb.models.delete(id);
    }
  });

  test('updateModel refuses a blank name and a call that changes nothing', async () => {
    const id = await insertModelRow({name: `forge-test-model-${Date.now()}`});
    try {
      await expectExceptionAsync(() => runUpdateModel({model: id, name: ' '}).then(() => {}),
        failsWith(NO_NAME));
      await expectExceptionAsync(() => runUpdateModel({model: id}).then(() => {}),
        failsWith('Give a name, a description or tags to change.'));
    } finally {
      await forgeDb.models.delete(id);
    }
  });

  test('Forge:listModels, getModel and updateModel are registered and callable', async () => {
    const signature = (name: string) => {
      const func = DG.Func.find({package: 'Forge', name})[0];
      return [...(func?.inputs ?? []).map((p) => `${p.name}: ${p.propertyType}`),
        ...(func?.outputs ?? []).map((p) => `-> ${p.propertyType}`)];
    };
    expectArray(signature('listModels'), ['name: string', 'target: string', 'task: string', 'method: string',
      'tags: string', 'table: dataframe', 'limit: int', '-> dataframe']);
    expectArray(signature('getModel'), ['model: string', '-> object']);
    expectArray(signature('updateModel'), ['model: string', 'name: string', 'description: string', 'tags: string',
      '-> object']);

    const {id} = shared();
    const models: unknown = await grok.functions.call('Forge:listModels', {});
    expect(models instanceof DG.DataFrame && models.col('Id') !== null, true, `Result: ${models}`);
    const details: unknown = await grok.functions.call('Forge:getModel', {model: id});
    expect(typeof details === 'object' && details !== null && 'id' in details && details.id === id, true,
      JSON.stringify(details));
    const stamp = Date.now();
    const edited = await insertModelRow({name: `forge-test-model-${stamp}`});
    try {
      const updated: unknown = await grok.functions.call('Forge:updateModel',
        {model: edited, name: `forge-test-called-${stamp}`});
      expect(typeof updated === 'object' && updated !== null && 'name' in updated &&
        updated.name === `forge-test-called-${stamp}`, true, JSON.stringify(updated));
    } finally {
      await forgeDb.models.delete(edited);
    }
  });
});

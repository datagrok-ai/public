import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, expectArray, expectExceptionAsync, expectFloat, test}
  from '@datagrok-libraries/test/src/test';
import schema from '../../databases/forge/schema.json';
import {applyAndRecord, DEFAULT_BATCH_SIZE, loadModel} from '../apply/apply-model';
import {exactMapping} from '../apply/column-matching';
import {defaultHyperparameters} from '../engines/engine';
import {EngineRegistry} from '../engines/engine-registry';
import {ForgeError} from '../forge-error';
import {forgeDb} from '../generated/db';
import {deleteTrainingCopy, trainingCopyName, uploadTrainingCopy} from '../storage/dataset-copy';
import {datasetFingerprint} from '../storage/dataset-fingerprint';
import {datasetRefOf, openDatasetRef, storedDatasetRef} from '../storage/dataset-ref';
import {ModelStorage, modelFieldsOf, normalizedTags, tagsOf, tagsText, trainingRunOf} from '../storage/model-fields';
import {BLOB_ROOT, deleteModel, modelsChanged, ownBlob, runDeleteModel, saveModel, saveTrainedModel}
  from '../storage/model-store';
import {linkTrainingRun, recordTrainingRun} from '../storage/training-run-store';
import {runTrainModel} from '../training/train-api';
import {FeaturesSchema, prepareTraining, TrainingResult, trainModel} from '../training/train-model';
import {columnsOf, deleteModelsNamed, deleteRunsOf, engineByName, expectReleased, failsWith, framesSharing,
  insertModelRow, IRIS, MEASUREMENTS, openDemoTable, openIris, openIrisFromFile, requestOf, selectionOf, XGBOOST_FIELDS}
  from './test-data';

const IRIS_COLUMNS = [...MEASUREMENTS, 'Species'];

/** Trains [engineName] on iris (Species by the measurements) and saves it with [storage]; returns the model id. */
async function saveIrisWith(iris: DG.DataFrame, engineName: string, storage: ModelStorage): Promise<string> {
  const engine = engineByName(EngineRegistry.discover(), engineName);
  const request = await prepareTraining({...selectionOf(columnsOf(iris, MEASUREMENTS), iris.getCol('Species')),
    engine, hyperparameters: defaultHyperparameters(engine)});
  const result = await trainModel(request);
  return await saveModel(modelFieldsOf({name: `forge-test-model-${Date.now()}`, description: '', tags: [],
    engine, datasetName: iris.name, result, fingerprint: datasetFingerprint(request.features, request.target),
    storage}), result.blob);
}

const smallTable = () => DG.DataFrame.fromColumns([DG.Column.fromList(DG.COLUMN_TYPE.INT, 'x', [1, 2])]);

async function tableExists(id: string): Promise<boolean> {
  return (await grok.dapi.tables.find(id)) != null;
}

category('Storage', () => {
  test('forge tables are registered', async () => {
    expect(Array.isArray(await forgeDb.models.query({limit: 1})), true);
    expect(Array.isArray(await forgeDb.trainingRuns.query({limit: 1})), true);
    expect(Array.isArray(await forgeDb.applications.query({limit: 1})), true);
  });

  test('a training run describes its options, features and data summary as a model does', async () => {
    const descriptions = async (table: string) => new Map((await grok.dapi.domains.registry.rowProperties(table))
      .map((p) => [p.name, p.description]));
    const [model, run] = [await descriptions('forge.model'), await descriptions('forge.training_run')];
    for (const column of ['options', 'features', 'hyperparameters', 'splitting', 'metrics', 'dataset_fingerprint']) {
      expect(run.get(column), model.get(column), column);
      expect((run.get(column) ?? '') !== '', true, `${column} has no description`);
    }
    expect(run.get('options')?.startsWith('Preparation replay; keys compatible with the old models'), true);
  });

  test('model, run and application round trip', async () => {
    const [{id: modelId}] = await forgeDb.models.insert({
      ...XGBOOST_FIELDS,
      name: `forge-test-model-${Date.now()}`,
      task: 'classification',
      target_name: 'species',
      features: {columns: [{name: 'sepal len', type: 'double'}]},
      options: {preprocessingInfo: ['one-hot']},
      storage_mode: 'none',
    });
    const applicationCount = () => forgeDb.applications.query().where('model_id', '=', modelId).count();
    let runId: string | undefined;
    try {
      const model = await forgeDb.models.get(modelId);
      expect(model.author_id, (await grok.dapi.users.current()).id);
      expect(model.features?.columns?.[0]?.name, 'sepal len');
      expect(model.options?.preprocessingInfo?.[0], 'one-hot');
      expect(model.storage_mode, 'none');

      [{id: runId}] = await forgeDb.trainingRuns.insert({
        ...XGBOOST_FIELDS,
        model_id: modelId,
        task: 'classification',
        target_name: 'species',
        status: 'completed',
        source: 'ui',
        started_on: new Date().toISOString(),
      });
      await forgeDb.applications.insert({model_id: modelId, row_count: 150, source: 'api', status: 'completed'});
      const [{id: cancelledId}] = await forgeDb.applications.insert({model_id: modelId, row_count: 150,
        skipped_rows: 3, source: 'ui', status: 'cancelled'});
      expect(await applicationCount(), 2);
      const cancelled = await forgeDb.applications.get(cancelledId);
      expect(cancelled.status, 'cancelled');
      expect(cancelled.skipped_rows, 3);
    } finally {
      try {
        if (runId !== undefined)
          await forgeDb.trainingRuns.delete(runId);
      } finally {
        await forgeDb.models.delete(modelId);
      }
    }
    expect(await applicationCount(), 0);
  });

  test('tags are stored as one comma-separated text and read back as a list', async () => {
    expectArray(normalizedTags([' b ', 'a', '', 'b', '  ', 'c']), ['b', 'a', 'c']);
    expect(tagsText([' a ', 'b', '', 'a']), 'a, b');
    expect(tagsText(['', ' ']), '');
    expectArray(tagsOf('a, b'), ['a', 'b']);
    expectArray(tagsOf(' a ,b,, a , '), ['a', 'b']);
    expectArray(tagsOf(''), []);
    expectArray(tagsOf(null), []);
    expectArray(tagsOf(undefined), []);
    expectArray(tagsOf(tagsText(['x, y', 'z'])), ['x', 'y', 'z']);
  });

  test('datasetFingerprint of iris', async () => {
    const iris = await grok.data.files.openTable(IRIS);
    const species = iris.getCol('Species');
    const fingerprint = datasetFingerprint(columnsOf(iris, MEASUREMENTS), species);
    expect(fingerprint.rowCount, 150);
    expect(fingerprint.columnCount, 5);
    expect(/^[0-9a-f]{8}$/.test(fingerprint.hash), true, `Hash ${fingerprint.hash}`);
    expect(fingerprint.columns[4].name, 'Species');
    expect(fingerprint.columns[4].categories?.length, 3);
    expectFloat(fingerprint.columns[0].min ?? NaN, 4.3, 1e-6, 'Sepal.Length min');
    expect(datasetFingerprint(columnsOf(iris, MEASUREMENTS), species).hash, fingerprint.hash);
    const swapped = ['Sepal.Width', 'Sepal.Length', 'Petal.Length', 'Petal.Width'];
    expect(datasetFingerprint(columnsOf(iris, swapped), species).hash !== fingerprint.hash, true,
      'Swapping two features keeps the hash');
  });

  test('saveModel writes the blob and the row in no-data mode', async () => {
    const stamp = Date.now();
    const iris = await openIris();
    const request = await requestOf(iris.clone(null, MEASUREMENTS), iris.getCol('Species'));
    const result = await trainModel(request);
    const id = await saveModel(modelFieldsOf({name: `forge-test-model-${stamp}`, description: '',
      tags: [], engine: request.engine, datasetName: iris.name, result,
      fingerprint: datasetFingerprint(request.features, request.target), storage: {mode: 'none'}}), result.blob);
    let path = '';
    let folder = '';
    try {
      const model = await forgeDb.models.get(id);
      expect(model.storage_mode, 'none');
      expect(model.tags == null, true, `Tags of a model without tags: '${model.tags}'`);
      expect((await grok.dapi.getEntities([id])).length, 1, 'The saved model is not a platform entity');
      expect((await forgeDb.models.get(id, {withAccess: true}))['~can_share'], true, 'Its author cannot share it');
      expect(model.blob?.startsWith(`file://${BLOB_ROOT}/`), true, `Blob ${model.blob}`);
      expect(model.has_training_rows, false);
      expect(typeof model.metrics?.validation?.accuracy, 'number');
      expect(model.features?.columns?.length, 4);
      expect(model.dataset_table_id == null, true, 'dataset_table_id is set');
      expect(model.dataset_ref == null, true, 'dataset_ref is set');
      expect(model.dataset_fingerprint?.rowCount, 150);
      path = (model.blob ?? '').substring('file://'.length);
      expect(await grok.dapi.files.exists(path), true, `${path} does not exist`);
      folder = path.substring(0, path.lastIndexOf('/'));
      expectArray((await grok.dapi.files.list(folder)).map((f) => f.name), ['model.bin']);
      expect((await grok.dapi.tables.filter(`name = "${iris.name}"`).list()).length, 0);
    } finally {
      await deleteModel(id);
    }
    expect(await grok.dapi.files.exists(path), false, `${path} is left`);
    expect(await grok.dapi.files.exists(folder), false, `${folder} is left`);
    expect(await forgeDb.models.query().where('id', '=', id).count(), 0);
  }, {timeout: 60000});

  test('deleteModel never deletes outside the model folders', async () => {
    const stamp = Date.now();
    const decoyFolder = `System:DomainFiles/forge/forge-test-${stamp}`;
    const decoy = `${decoyFolder}/decoy.bin`;
    const blobs = [`file://${decoy}`, `file://${BLOB_ROOT}//model.bin`, `file://${BLOB_ROOT}/../decoy.bin`];
    const modelFolders = async () => (await grok.dapi.files.list(BLOB_ROOT)).map((f) => f.name).sort();
    await grok.dapi.files.write(decoy, new Uint8Array([1, 2, 3]));
    let pendingId: string | undefined;
    try {
      const before = await modelFolders();
      for (const blob of blobs) {
        [{id: pendingId}] = await forgeDb.models.insert({...XGBOOST_FIELDS, name: `forge-test-model-${stamp}`,
          task: 'regression', target_name: 'y', storage_mode: 'none', blob});
        await deleteModel(pendingId);
        pendingId = undefined;
      }
      expect(await grok.dapi.files.exists(decoy), true, 'The decoy file was deleted');
      expectArray(await modelFolders(), before);
    } finally {
      try {
        if (pendingId !== undefined)
          await forgeDb.models.delete(pendingId);
      } finally {
        await grok.dapi.files.delete(decoyFolder);
      }
    }
  });

  test('datasetRefOf reads the origin the platform recorded', async () => {
    const iris = await openIrisFromFile();
    try {
      const ref = datasetRefOf(iris);
      expect(ref?.kind, 'file', `The tags: ${JSON.stringify([...iris.tags.entries()])}`);
      expect(ref?.path, IRIS);
      expect(/^\w+ = OpenFile\("System:DemoFiles\/iris\.csv"\)$/.test(ref?.script ?? ''), true, ref?.script);
    } finally {
      grok.shell.closeTable(iris);
    }
    expect(datasetRefOf(await grok.data.files.openTable(IRIS)) === null, true,
      'grok.data.files.openTable recorded an origin');
    expect(datasetRefOf(smallTable()) === null, true, 'A frame built in code has an origin');

    const local = smallTable();
    local.setTag(DG.Tags.SourceFile, 'iris.csv');
    expect(datasetRefOf(local) === null, true, 'A local file has an origin');
    local.setTag(DG.Tags.SourceFile, IRIS);
    expect(datasetRefOf(local)?.path, IRIS);

    const query = smallTable();
    query.setTag(DG.Tags.CreationScript, 'orders = Samples:Orders(country="USA") //{"timestamp": 1}');
    query.setTag(DG.Tags.DataQueryId, 'f0e1d2c3-0000-4000-8000-000000000000');
    query.setTag(DG.Tags.DataQueryName, 'Orders');
    expect(JSON.stringify(datasetRefOf(query)), JSON.stringify({kind: 'query',
      script: 'orders = Samples:Orders(country="USA")', id: 'f0e1d2c3-0000-4000-8000-000000000000', name: 'Orders'}));
  });

  test('openDatasetRef opens a file reference outside the workspace', async () => {
    const iris = await openIrisFromFile();
    const ref = datasetRefOf(iris);
    grok.shell.closeTable(iris);
    if (ref === null)
      throw new Error('iris has no origin');
    const tables = grok.shell.tables.length;
    const opened = await openDatasetRef(ref);
    expect(opened.rowCount, 150);
    expectArray(IRIS_COLUMNS.map((name) => opened.col(name) !== null), IRIS_COLUMNS.map(() => true));
    expect(grok.shell.tables.length, tables, 'The table was added to the workspace');
    const stored = storedDatasetRef(JSON.parse(JSON.stringify(ref)));
    expect(JSON.stringify(stored), JSON.stringify(ref));
    expect(storedDatasetRef({kind: 'url', script: 'x = F()'}) === null, true, 'An unknown kind is read');
    expect(storedDatasetRef({kind: 'file'}) === null, true, 'A reference without a script is read');
  });

  test('openDatasetRef refuses a script that is not one plain call per line', async () => {
    const scripts = ['data = OpenFile("x"); grok.shell.info("x")', 'data = OpenFile(Evil())', 'grok.shell.info("x")',
      'data = OpenFile("x")\nEvil()'];
    const isRefusal = (e: unknown) => e instanceof ForgeError && e.message.includes('not a script Forge can run');
    for (const script of scripts) {
      await expectExceptionAsync(async () => {
        await openDatasetRef({kind: 'script', script});
      }, isRefusal);
    }
  });

  test('uploadTrainingCopy and deleteTrainingCopy', async () => {
    const iris = await openIris();
    const columns = columnsOf(iris, IRIS_COLUMNS);
    const modelName = `forge-test-copy-${Date.now()}`;
    const [id, frames] = await framesSharing(columns, () => uploadTrainingCopy(columns, modelName));
    let decoyId: string | undefined;
    try {
      expectReleased(frames);
      const info = await grok.dapi.tables.find(id);
      expect(info?.friendlyName, trainingCopyName(modelName));
      const copy = await grok.dapi.tables.getTable(id);
      expect(copy.rowCount, 150);
      expectArray(copy.columns.names(), IRIS_COLUMNS);
      expect(iris.name.startsWith('forge-test-iris-'), true, 'The table was renamed');

      const decoy = smallTable();
      decoy.name = `forge-test-decoy-${Date.now()}`;
      decoyId = await grok.dapi.tables.uploadDataFrame(decoy);
      await deleteTrainingCopy(decoyId);
      expect(await tableExists(decoyId), true, 'A table that is not a training copy was deleted');
    } finally {
      try {
        await deleteTrainingCopy(id);
      } finally {
        if (decoyId !== undefined)
          await grok.dapi.tables.delete(await grok.dapi.tables.find(decoyId));
      }
    }
    expect(await tableExists(id), false, 'The copy is left');
    await deleteTrainingCopy(id);
  }, {timeout: 60000});

  test('a copy-mode model keeps the table id, and deleteModel deletes the table', async () => {
    const iris = await openIris();
    const tableId = await uploadTrainingCopy(columnsOf(iris, IRIS_COLUMNS), `forge-test-model-${Date.now()}`);
    let modelId: string | undefined;
    try {
      modelId = await saveIrisWith(iris, 'XGBoost', {mode: 'copy', tableId});
      const model = await forgeDb.models.get(modelId);
      expect(model.storage_mode, 'copy');
      expect(model.dataset_table_id, tableId);
      expect(model.dataset_ref == null, true, 'dataset_ref is set');
      expect(model.has_training_rows, false);
    } finally {
      if (modelId !== undefined)
        await deleteModel(modelId);
      else
        await deleteTrainingCopy(tableId);
    }
    expect(await tableExists(tableId), false, 'The copy outlived its model');
  }, {timeout: 60000});

  test('a reference-mode SVM model keeps the reference and says it holds training rows', async () => {
    const iris = await openIrisFromFile();
    const ref = datasetRefOf(iris);
    iris.name = `forge-test-iris-${Date.now()}`;
    const saved = ref === null ? Promise.reject(new Error('iris has no origin')) :
      saveIrisWith(iris, 'SVM', {mode: 'reference', ref});
    const id = await saved.finally(() => grok.shell.closeTable(iris));
    try {
      const model = await forgeDb.models.get(id);
      expect(model.storage_mode, 'reference');
      expect(model.engine_name, 'SVM');
      expect(model.has_training_rows, true);
      expect(model.dataset_table_id == null, true, 'dataset_table_id is set');
      expect(JSON.stringify(storedDatasetRef(model.dataset_ref)), JSON.stringify(ref));
    } finally {
      await deleteModel(id);
    }
  }, {timeout: 60000});

  test('recordTrainingRun and linkTrainingRun', async () => {
    const stamp = Date.now();
    const values = Array.from({length: 12}, (_, i) => i + 1);
    const x = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'x', values);
    const y = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'y', values.map((v) => v * 2));
    const request = await requestOf(DG.DataFrame.fromColumns([x]), y);
    const runId = await recordTrainingRun(trainingRunOf({request, datasetName: `forge-test-run-${stamp}`,
      fingerprint: datasetFingerprint(request.features, y), status: 'completed', source: 'ui',
      startedOn: new Date().toISOString(), durationMs: 0}));
    let modelId: string | undefined;
    try {
      const run = await forgeDb.trainingRuns.get(runId);
      expect(run.status, 'completed');
      expect(run.task, 'regression');
      expect(run.features?.columns?.[0]?.name, 'x');
      expect(run.model_id == null, true, `model_id ${run.model_id}`);
      [{id: modelId}] = await forgeDb.models.insert({...XGBOOST_FIELDS, name: `forge-test-model-${stamp}`,
        task: 'regression', target_name: 'y', storage_mode: 'none'});
      await linkTrainingRun(runId, modelId);
      expect((await forgeDb.trainingRuns.get(runId)).model_id, modelId);
    } finally {
      try {
        await forgeDb.trainingRuns.delete(runId);
      } finally {
        if (modelId !== undefined)
          await forgeDb.models.delete(modelId);
      }
    }
  });
});

/** The uploaded training copies of the model named [modelName]. */
async function trainingCopiesOf(modelName: string): Promise<DG.TableInfo[]> {
  const tables: DG.TableInfo[] = await grok.dapi.tables.filter(modelName).list();
  return tables.filter((t) => t.friendlyName === trainingCopyName(modelName));
}

const modelCount = (id: string) => forgeDb.models.query().where('id', '=', id).count();

// Saving through Forge:trainModel and saveTrainedModel, Forge:deleteModel, the schema.
category('Storage', () => {
  test('the schema is the package\'s and a training run records its source', async () => {
    const manifest = await grok.dapi.domains.schema('forge').manifest();
    expect(manifest.version, schema.version);
    const source = (await grok.dapi.domains.registry.rowProperties('forge.training_run'))
      .find((p) => p.name === 'source');
    expect(source?.description, 'Where the training was started: the Train view or Forge:trainModel');
    expectArray(source?.choices ?? [], ['ui', 'api']);
  });

  test('models and training runs are promoted eagerly', async () => {
    const {tables} = await grok.dapi.domains.schema('forge').manifest();
    expect(tables?.model?.promotion, 'eager');
    expect(tables?.training_run?.promotion, 'eager');
  });

  test('Forge:trainModel with a name saves the model without data and links the run', async () => {
    const iris = await openDemoTable(IRIS, 'api');
    const name = `forge-test-api-${Date.now()}`;
    let changes = 0;
    const sub = modelsChanged.subscribe(() => changes++);
    let modelId: string | null = null;
    try {
      const result = await runTrainModel({table: iris, target: 'Species', name: ` ${name} `, description: 'API',
        tags: 'a, b,a', showProgress: false});
      modelId = result.modelId;
      expect(result.name, name);
      expect(changes, 1);
      const model = await forgeDb.models.get(modelId ?? '');
      expect(model.name, name);
      expect(model.storage_mode, 'none');
      expect(model.description, 'API');
      expect(model.tags, 'a, b');
      expect(model.row_count, 150);
      expect(model.dataset_name, iris.name);
      expect(model.dataset_table_id == null && model.dataset_ref == null, true, 'Data was stored');
      expect((await forgeDb.trainingRuns.get(result.runId)).model_id, modelId);
      expect((await trainingCopiesOf(name)).length, 0);
    } finally {
      sub.unsubscribe();
      await deleteModelsNamed(name);
      await deleteRunsOf([iris]);
    }
  }, {timeout: 60000});

  test('Forge:trainModel with storage copy uploads the training columns; Forge:deleteModel deletes them', async () => {
    const iris = await openDemoTable(IRIS, 'api');
    const name = `forge-test-api-copy-${Date.now()}`;
    let modelId: string | null = null;
    try {
      modelId = (await runTrainModel({table: iris, target: 'Species', name, storage: 'copy', showProgress: false}))
        .modelId;
      const model = await forgeDb.models.get(modelId ?? '');
      expect(model.storage_mode, 'copy');
      const tableId = model.dataset_table_id ?? '';
      expectArray((await trainingCopiesOf(name)).map((t) => t.id), [tableId]);
      expectArray((await grok.dapi.tables.getTable(tableId)).columns.names(), IRIS_COLUMNS);
      expect(await runDeleteModel(name), modelId);
      expect((await trainingCopiesOf(name)).length, 0, 'The copy outlived its model');
    } finally {
      await deleteModelsNamed(name);
      await deleteRunsOf([iris]);
    }
  }, {timeout: 60000});

  test('Forge:trainModel with storage reference keeps the table\'s origin', async () => {
    const iris = await openIrisFromFile();
    iris.name = `forge-test-api-${Date.now()}`;
    let modelId: string | null = null;
    try {
      modelId = (await runTrainModel({table: iris, target: 'Species', name: iris.name, storage: 'reference',
        showProgress: false})).modelId;
      const model = await forgeDb.models.get(modelId ?? '');
      expect(model.storage_mode, 'reference');
      expect(storedDatasetRef(model.dataset_ref)?.path, IRIS);
      expect(model.dataset_table_id == null, true, 'dataset_table_id is set');
    } finally {
      grok.shell.closeTable(iris);
      await deleteModelsNamed(iris.name);
      await deleteRunsOf([iris]);
    }
  }, {timeout: 60000});

  test('saveTrainedModel deletes the uploaded copy when the save fails', async () => {
    const iris = await openIris();
    const request = await requestOf(columnsOf(iris, MEASUREMENTS), iris.getCol('Species'));
    const result = await trainModel(request);
    const name = `forge-test-api-failed-${Date.now()}`;
    // A result without a feature list fails the row's fields after the upload, before any model file is written.
    const broken: TrainingResult = {...result, features: {} as FeaturesSchema};
    try {
      await expectExceptionAsync(() => saveTrainedModel({name, description: '', tags: [], engine: request.engine,
        datasetName: iris.name, result: broken, fingerprint: datasetFingerprint(request.features, request.target),
        columns: columnsOf(iris, IRIS_COLUMNS), runId: '', storage: {mode: 'copy'}}).then(() => {}));
      expect((await trainingCopiesOf(name)).length, 0, 'The copy outlived the failed save');
      expect(await forgeDb.models.query().where('name', '=', name).count(), 0);
    } finally {
      await deleteModelsNamed(name);
    }
  }, {timeout: 60000});

  test('Forge:deleteModel by id and by name: the row, the file and the applications go, the run stays', async () => {
    const iris = await openDemoTable(IRIS, 'api');
    const names = ['id', 'name'].map((by) => `forge-test-api-delete-${by}-${Date.now()}`);
    try {
      for (const [i, by] of ['id', 'name'].entries()) {
        const name = names[i];
        const {modelId, runId} = await runTrainModel({table: iris, target: 'Species', name, showProgress: false});
        if (modelId === null)
          throw new Error('The model was not saved');
        const model = await loadModel(modelId);
        const table = iris.clone();
        await applyAndRecord({model, table, mapping: exactMapping(model.features, table),
          batchSize: DEFAULT_BATCH_SIZE, missingValues: {mode: 'skip'}}, 'api');
        const applications = () => forgeDb.applications.query().where('model_id', '=', modelId).count();
        expect(await applications(), 1);
        const folder = ownBlob((await forgeDb.models.get(modelId)).blob)?.folder ?? '';
        expect(await grok.dapi.files.exists(folder), true, folder);
        expect(await runDeleteModel(by === 'id' ? modelId : name), modelId);
        expect(await modelCount(modelId), 0, `The model is left (by ${by})`);
        expect(await grok.dapi.files.exists(folder), false, `${folder} is left`);
        expect(await applications(), 0);
        expect((await forgeDb.trainingRuns.get(runId)).model_id == null, true, 'The run keeps the model');
      }
    } finally {
      await Promise.all(names.map(deleteModelsNamed));
      await deleteRunsOf([iris]);
    }
  }, {timeout: 90000});

  test('Forge:deleteModel refuses a missing model and a name several models have', async () => {
    const stamp = Date.now();
    const name = `forge-test-api-dup-${stamp}`;
    const missing = `forge-test-api-none-${stamp}`;
    const ids = [await insertModelRow({name}), await insertModelRow({name})];
    try {
      await expectExceptionAsync(() => runDeleteModel(name).then(() => {}),
        failsWith(`Several models are named '${name}'. Use the model id.`));
      await expectExceptionAsync(() => runDeleteModel(missing).then(() => {}),
        failsWith(`No Forge model '${missing}' is available to you.`));
      expect(await forgeDb.models.query().where('name', '=', name).count(), 2);
    } finally {
      for (const id of ids)
        await forgeDb.models.delete(id);
    }
  });

  test('Forge:deleteModel through grok.functions.call, and its registration', async () => {
    const name = `forge-test-api-call-${Date.now()}`;
    const id = await insertModelRow({name});
    try {
      expect(await grok.functions.call('Forge:deleteModel', {model: name}), id);
      expect(await modelCount(id), 0);
    } finally {
      await deleteModelsNamed(name);
    }
    const func = DG.Func.find({package: 'Forge', name: 'deleteModel'})[0];
    expectArray(func.inputs.map((p) => `${p.name}: ${p.propertyType}`), ['model: string']);
    expectArray(func.outputs.map((p) => `${p.name}: ${p.propertyType}`), ['id: string']);
  });
});

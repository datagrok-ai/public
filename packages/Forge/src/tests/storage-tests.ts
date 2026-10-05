import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, expectArray, expectFloat, test} from '@datagrok-libraries/test/src/test';
import {forgeDb} from '../generated/db';
import {datasetFingerprint} from '../storage/dataset-fingerprint';
import {modelFieldsOf, trainingRunOf} from '../storage/model-fields';
import {BLOB_ROOT, deleteModel, saveModel} from '../storage/model-store';
import {linkTrainingRun, recordTrainingRun} from '../storage/training-run-store';
import {trainModel} from '../training/train-model';
import {IRIS, MEASUREMENTS, openIris, requestOf, XGBOOST_FIELDS} from './test-data';

category('Storage', () => {
  test('forge tables are registered', async () => {
    expect(Array.isArray(await forgeDb.models.query({limit: 1})), true);
    expect(Array.isArray(await forgeDb.trainingRuns.query({limit: 1})), true);
    expect(Array.isArray(await forgeDb.applications.query({limit: 1})), true);
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

  test('datasetFingerprint of iris', async () => {
    const iris = await grok.data.files.openTable(IRIS);
    const species = iris.getCol('Species');
    const fingerprint = datasetFingerprint(iris.clone(null, MEASUREMENTS), species);
    expect(fingerprint.rowCount, 150);
    expect(fingerprint.columnCount, 5);
    expect(/^[0-9a-f]{8}$/.test(fingerprint.hash), true, `Hash ${fingerprint.hash}`);
    expect(fingerprint.columns[4].name, 'Species');
    expect(fingerprint.columns[4].categories?.length, 3);
    expectFloat(fingerprint.columns[0].min ?? NaN, 4.3, 1e-6, 'Sepal.Length min');
    expect(datasetFingerprint(iris.clone(null, MEASUREMENTS), species).hash, fingerprint.hash);
    const swapped = ['Sepal.Width', 'Sepal.Length', 'Petal.Length', 'Petal.Width'];
    expect(datasetFingerprint(iris.clone(null, swapped), species).hash !== fingerprint.hash, true,
      'Swapping two features keeps the hash');
  });

  test('saveModel writes the blob and the row in no-data mode', async () => {
    const stamp = Date.now();
    const iris = await openIris();
    const request = await requestOf(iris.clone(null, MEASUREMENTS), iris.getCol('Species'));
    const result = await trainModel(request);
    const id = await saveModel(modelFieldsOf({name: `forge-test-model-${stamp}`, description: '',
      engine: request.engine, datasetName: iris.name, result,
      fingerprint: datasetFingerprint(request.features, request.target)}), result.blob);
    let path = '';
    let folder = '';
    try {
      const model = await forgeDb.models.get(id);
      expect(model.storage_mode, 'none');
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

  test('recordTrainingRun and linkTrainingRun', async () => {
    const stamp = Date.now();
    const values = Array.from({length: 12}, (_, i) => i + 1);
    const x = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'x', values);
    const y = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'y', values.map((v) => v * 2));
    const request = await requestOf(DG.DataFrame.fromColumns([x]), y);
    const runId = await recordTrainingRun(trainingRunOf({request, datasetName: `forge-test-run-${stamp}`,
      fingerprint: datasetFingerprint(request.features, y), status: 'completed',
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

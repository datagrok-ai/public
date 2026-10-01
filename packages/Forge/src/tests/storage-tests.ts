import * as grok from 'datagrok-api/grok';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {forgeDb, ModelInsert} from '../generated/db';

category('Storage', () => {
  test('forge tables are registered', async () => {
    expect(Array.isArray(await forgeDb.models.query({limit: 1})), true);
    expect(Array.isArray(await forgeDb.trainingRuns.query({limit: 1})), true);
    expect(Array.isArray(await forgeDb.applications.query({limit: 1})), true);
  });

  test('model, run and application round trip', async () => {
    const engine: Pick<ModelInsert, 'engine_name' | 'engine_namespace' | 'engine_kind'> =
      {engine_name: 'XGBoost', engine_namespace: 'Eda', engine_kind: 'function'};
    const [{id: modelId}] = await forgeDb.models.insert({
      ...engine,
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
        ...engine,
        model_id: modelId,
        task: 'classification',
        target_name: 'species',
        status: 'completed',
        started_on: new Date().toISOString(),
      });
      await forgeDb.applications.insert({model_id: modelId, row_count: 150, source: 'api', status: 'completed'});
      expect(await applicationCount(), 1);
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
});

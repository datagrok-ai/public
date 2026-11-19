import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, expectArray, expectExceptionAsync, expectFloat, test}
  from '@datagrok-libraries/test/src/test';
import {apiMapping, DEFAULT_BATCH_SIZE, featureSchemasOf, loadModel, predictRows, runApplyModel}
  from '../apply/apply-model';
import {runGetModel} from '../catalog/get-model';
import {CHECK_METHOD, checkApply, checkTrain, isCheckApplicable} from '../engines/check-method';
import {featureSemTypesOf, isComplete, isServerEngine} from '../engines/engine';
import {EngineRegistry} from '../engines/engine-registry';
import {errorMessage} from '../forge-error';
import {forgeDb} from '../generated/db';
import {MetricId, MetricValues, regressionMetrics} from '../metrics/metrics';
import {statisticsOf} from '../metrics/statistics';
import {runEvaluateModel} from '../training/evaluate-api';
import {runTrainModel} from '../training/train-api';
import {deleteModelsNamed, deleteRunsOf, engineByName, failsWith, IMPUTE, openSolubility, smilesColumn, smilesSeries,
  valuesOf} from './test-data';

// The temporary check method Forge Check (molecules): every test here goes with it.
const TRAIN_TIMEOUT = 120000;
const APPLY_TIMEOUT = 180000;
const REGRESSION_IDS: MetricId[] = ['mse', 'rmse', 'mae', 'r2'];

function floats(name: string, values: (number | null)[]): DG.Column {
  return DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, name, values);
}

function expectFinite(values: MetricValues, label: string): void {
  for (const id of REGRESSION_IDS)
    expect(Number.isFinite(values[id]), true, `${label} ${id}: ${values[id]}`);
}

category('Engines', () => {
  test('Forge Check: discovered as a browser method that takes molecules', async () => {
    const engine = engineByName(EngineRegistry.discover(), CHECK_METHOD);
    expect(isComplete(engine), true, 'Not complete');
    expect(`${engine.kind} ${engine.namespace}`, 'function Forge');
    expect(isServerEngine(engine), false, 'A server method');
    expectArray(featureSemTypesOf(engine), ['molecule']);
    const smiles = smilesColumn('smiles', ['CO', 'CCO', 'CCCO']);
    const x = floats('x', [1, 2, 3]);
    const text = DG.Column.fromStrings('text', ['a', 'b', 'a']);
    const y = floats('y', [1, 2, 3]);
    expect(isCheckApplicable(DG.DataFrame.fromColumns([smiles, x]), y), true, '[smiles, x]');
    expect(isCheckApplicable(DG.DataFrame.fromColumns([floats('x', [1, 2, 3])]), y), false, '[x]');
    expect(isCheckApplicable(DG.DataFrame.fromColumns([smilesColumn('smiles', ['CO', 'CCO', 'CCCO']), text]), y),
      false, '[smiles, text]');
    expect(isCheckApplicable(DG.DataFrame.fromColumns([smilesColumn('smiles', ['CO', 'CCO', 'CCCO'])]),
      DG.Column.fromStrings('label', ['a', 'b', 'a'])), false, 'A text target');
  });

  test('Forge Check: fits a line and predicts zeros from molecules alone', async () => {
    const rows = 20;
    const smiles = smilesColumn('smiles', smilesSeries(rows));
    const x = floats('x', valuesOf(rows, (i) => i));
    const y = floats('y', valuesOf(rows, (i) => 2 * i + 1));
    const fitted = checkApply(DG.DataFrame.fromColumns([smiles, x]), checkTrain(DG.DataFrame.fromColumns([smiles, x]),
      y)).columns.byIndex(0);
    for (let i = 0; i < rows; i++)
      expectFloat(fitted.getNumber(i), 2 * i + 1, 1e-9, `row ${i}`);
    const molecules = DG.DataFrame.fromColumns([smilesColumn('smiles', smiles.toList())]);
    expectArray(checkApply(molecules, checkTrain(molecules, y)).columns.byIndex(0).toList(),
      Array.from({length: rows}, () => 0));
  });
});

category('Training', () => {
  test('Forge Check: all-zero predictions give finite metrics', async () => {
    const table = await openSolubility();
    try {
      const sol = table.getCol('SOL');
      const zeros = DG.Column.fromFloat64Array('zeros', new Float64Array(sol.length));
      const metrics = regressionMetrics(sol, zeros);
      expectFinite(metrics, 'zeros');
      const values = sol.toList().map(Number);
      const mean = values.reduce((sum, v) => sum + v, 0) / values.length;
      const r2 = 1 - values.reduce((sum, v) => sum + v * v, 0) / values.reduce((sum, v) => sum + (v - mean) ** 2, 0);
      expectFloat(metrics.r2 ?? NaN, r2, 1e-9, 'r2');
      // The skewness and kurtosis of a constant column are undefined.
      for (const {id, actual, predicted} of statisticsOf(sol, zeros)) {
        expect(Number.isFinite(actual), true, `actual ${id}: ${actual}`);
        if (id !== 'skew' && id !== 'kurt')
          expect(predicted, 0, `predicted ${id}`);
      }
      const result = await runTrainModel({table, target: 'SOL', features: ['smiles'], method: CHECK_METHOD,
        showProgress: false});
      expect(result.method, CHECK_METHOD);
      expectFinite(result.metrics.train, 'train');
      expectFinite(result.metrics.validation, 'validation');
    } finally {
      await deleteRunsOf([table]);
    }
  }, {timeout: TRAIN_TIMEOUT});

  test('Forge Check: suggested for molecules and numbers, refused for XGBoost', async () => {
    const table = await openSolubility();
    try {
      const mixed = await runTrainModel({table, target: 'SOL', features: ['smiles', 'MolLogP', 'TPSA'],
        showProgress: false});
      expect(mixed.method, CHECK_METHOD);
      expect((mixed.metrics.train.r2 ?? 0) > 0, true, `R2 ${mixed.metrics.train.r2}`);
      await expectExceptionAsync(() => runTrainModel({table, target: 'SOL', features: ['smiles', 'MolLogP', 'TPSA'],
        method: 'XGBoost', showProgress: false}).then(() => {}),
      failsWith('XGBoost cannot learn from Molecule features. Uncheck: smiles.'));
      const numbers = await runTrainModel({table, target: 'SOL', features: ['MolLogP', 'TPSA'], showProgress: false});
      expect(numbers.method !== CHECK_METHOD, true, numbers.method);
    } finally {
      await deleteRunsOf([table]);
    }
  }, {timeout: TRAIN_TIMEOUT});
});

category('Apply', () => {
  test('Forge Check: a saved molecule model applies and evaluates', async () => {
    const table = await openSolubility();
    const name = `forge-test-mol-${Date.now()}`;
    let modelId: string | null = null;
    try {
      const result = await runTrainModel({table, target: 'SOL', features: ['smiles', 'MolLogP'], name,
        storage: 'none', showProgress: false});
      modelId = result.modelId;
      if (modelId === null)
        throw new Error('The model was not saved');
      expect((await runGetModel(modelId, [])).featureTypes.smiles, 'string');
      const stored = featureSchemasOf((await forgeDb.models.get(modelId)).features) ?? [];
      expect(stored.find((f) => f.name === 'smiles')?.semType, 'Molecule');

      const fresh = await openSolubility();
      await runApplyModel(modelId, fresh, null, false);
      const loaded = await loadModel(modelId);
      const blob = await grok.dapi.files.readAsBytes(loaded.blobPath);
      const expected = checkApply(fresh.clone(null, ['smiles', 'MolLogP']), blob).columns.byIndex(0);
      const predicted = fresh.getCol('SOL (predicted)');
      for (let i = 0; i < fresh.rowCount; i++)
        expectFloat(predicted.getNumber(i), expected.getNumber(i), 1e-9, `row ${i}`);

      const unmarked = DG.DataFrame.fromColumns([DG.Column.fromStrings('smiles', fresh.getCol('smiles').toList()),
        fresh.getCol('MolLogP').clone()]);
      expect(unmarked.getCol('smiles').semType ?? '', '', 'The copy is marked');
      await expectExceptionAsync(() => runApplyModel(modelId ?? '', unmarked, null, false).then(() => {}),
        (e) => errorMessage(e).startsWith('\'smiles\' is not a Molecule column'));

      const gapped = fresh.clone(null, ['smiles', 'MolLogP']);
      gapped.getCol('smiles').semType = DG.SEMTYPE.MOLECULE;
      gapped.getCol('smiles').set(3, null);
      await expectExceptionAsync(() => predictRows({model: loaded, table: gapped, mapping: apiMapping(loaded, gapped,
        null), batchSize: DEFAULT_BATCH_SIZE, missingValues: IMPUTE}).then(() => {}),
      (e) => errorMessage(e).includes('Impute cannot fill the Molecule column \'smiles\'. Choose Skip rows.'));

      const evaluation = await runEvaluateModel({model: modelId, table: fresh, showProgress: false});
      expectFinite(evaluation.metrics.evaluation, 'evaluation');
    } finally {
      await deleteModelsNamed(name);
      await deleteRunsOf([table]);
    }
  }, {timeout: APPLY_TIMEOUT});
});

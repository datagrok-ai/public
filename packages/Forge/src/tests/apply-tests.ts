import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, expect, expectArray, expectExceptionAsync, expectFloat, test}
  from '@datagrok-libraries/test/src/test';
import {applyAndRecord, ApplyRequest, DEFAULT_BATCH_SIZE, LoadedModel, loadModel, predictionName}
  from '../apply/apply-model';
import {compatibility, exactMapping, isSuggested, MAX_NAME_DISTANCE, mappedColumns, mappingProblems, nameDistance,
  suggestMapping} from '../apply/column-matching';
import {featureFrame} from '../apply/feature-frame';
import {PREDICTION_TAG} from '../constants';
import {Engine} from '../engines/engine';
import {apply, LoopProgress} from '../engines/engine-calls';
import {errorMessage, ForgeError} from '../forge-error';
import {forgeDb, ModelInsert} from '../generated/db';
import {metricsOf, regressionMetrics} from '../metrics/metrics';
import {MissingValuesSettings} from '../preparation/missing-values';
import {ONE_HOT} from '../preparation/preparation-options';
import {releaseFrame, sharedFrame} from '../preparation/shared-frame';
import {BLOB_ROOT, deleteModel} from '../storage/model-store';
import {ColumnSchema, prepareTraining, TrainingResult, trainModel} from '../training/train-model';
import {columnsOf, expectMetrics, expectReleased, framesSharing, IMPUTE, MEASUREMENTS, openIris, rawValues, requestOf,
  saveIrisModel, saveTestModel, selectionOf, valuesOf, XGBOOST_FIELDS} from './test-data';

const TIMEOUT = 90000;
const IRIS_FEATURES: ColumnSchema[] = MEASUREMENTS.map((name) => ({name, type: DG.COLUMN_TYPE.FLOAT}));

interface Fixture { id: string; result: TrainingResult; iris: DG.DataFrame; model: LoadedModel; reference: unknown[] }

let fixture: Fixture | undefined;
let fixtureId: string | undefined;

/** The engine's predictions for [columns], as a list. */
async function predictionsOf(engine: Engine, columns: DG.Column[], blob: Uint8Array): Promise<unknown[]> {
  const frame = sharedFrame(columns);
  try {
    return (await apply(engine, frame, blob)).toList();
  } finally {
    releaseFrame(frame);
  }
}

function shared(): Fixture {
  if (fixture === undefined)
    throw new Error('The test model was not saved');
  return fixture;
}

function requestFor(model: LoadedModel, table: DG.DataFrame, options: {mapping?: Map<string, string>;
  batchSize?: number; missingValues?: MissingValuesSettings} = {}): ApplyRequest {
  return {model, table, mapping: options.mapping ?? exactMapping(model.features, table),
    batchSize: options.batchSize ?? DEFAULT_BATCH_SIZE, missingValues: options.missingValues ?? {mode: 'skip'}};
}

function irisCopy(name: string): DG.DataFrame {
  const copy = shared().iris.clone();
  copy.name = `${shared().iris.name}-${name}`;
  return copy;
}

async function applicationsOf(table: DG.DataFrame) {
  return forgeDb.applications.query().where('model_id', '=', shared().id).where('table_name', '=', table.name);
}

const throwsWith = (text: string) => (e: unknown) => errorMessage(e).includes(text);

category('Apply', () => {
  before(async () => {
    const iris = await openIris();
    const {id, result} = await saveIrisModel(iris);
    fixtureId = id;
    const model = await loadModel(id);
    fixture = {id, result, iris, model,
      reference: await predictionsOf(model.engine, columnsOf(iris, MEASUREMENTS), result.blob)};
  });

  after(async () => {
    if (fixtureId !== undefined)
      await deleteModel(fixtureId);
  });

  test('nameDistance splits the calibration pairs at MAX_NAME_DISTANCE', async () => {
    const pairs: [string, string, number][] = [['Sepal.Width', 'sepal width', 0.04],
      ['Sepal.Length', 'SepalLengthCm', 0.05], ['age', 'age_years', 0.16], ['Petal.Width', 'pH', 0.42],
      ['Sepal.Length', 'density', 0.55], ['Petal.Width', 'col 1', 0.57]];
    pairs.forEach(([feature, column, figure], i) => {
      const distance = nameDistance(feature, column);
      expectFloat(distance, figure, 0.021, `${feature} / ${column}`);
      expect(distance <= MAX_NAME_DISTANCE, i < 3, `${feature} / ${column}: ${distance}`);
    });
  });

  test('exactMapping matches names ignoring case and expands one-hot columns', async () => {
    const table = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'sepal.length', [1, 2]),
      DG.Column.fromList(DG.COLUMN_TYPE.INT, 'color=red', [1, 0]),
      DG.Column.fromList(DG.COLUMN_TYPE.INT, 'color=blue', [0, 1]),
    ]);
    const color: ColumnSchema = {name: 'color', type: DG.COLUMN_TYPE.STRING};
    const features: ColumnSchema[] = [{name: 'Sepal.Length', type: DG.COLUMN_TYPE.FLOAT}, color,
      {name: 'z', type: DG.COLUMN_TYPE.FLOAT}];
    const mapping = exactMapping(features, table);
    expect(mapping.get('Sepal.Length'), 'sepal.length');
    expectArray(mappedColumns(color, mapping).map(([, column]) => column), ['color=red', 'color=blue']);
    expect(mapping.has('z'), false);
    const problems = mappingProblems(features, mapping, table, 'XGBoost');
    expectArray(problems.map((p) => p.feature), ['z']);
    expectArray(featureFrame(table, features.slice(0, 2), mapping).columns.names(),
      ['Sepal.Length', 'color=red', 'color=blue']);
  });

  test('suggestMapping takes close names of compatible columns only', async () => {
    const iris = await openIris();
    iris.getCol('Sepal.Length').name = 'sepal length';
    iris.getCol('Sepal.Width').name = 'sepal width';
    const sepal: ColumnSchema[] = IRIS_FEATURES.slice(0, 2);
    const mapping = suggestMapping(sepal, iris);
    expect(mapping.get('Sepal.Length'), 'sepal length');
    expect(mapping.get('Sepal.Width'), 'sepal width');
    expect(isSuggested(sepal, iris), true);

    const petal: ColumnSchema = {name: 'Petal.Width', type: DG.COLUMN_TYPE.FLOAT};
    const ph = DG.DataFrame.fromColumns([DG.Column.fromList(DG.COLUMN_TYPE.INT, 'col 1', [1, 2]),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'pH', [3.1, 3.2])]);
    expect(suggestMapping([petal], ph).size, 0);
    expect(isSuggested([petal], ph), false);

    const single = (col: DG.Column, feature: ColumnSchema) =>
      suggestMapping([feature], DG.DataFrame.fromColumns([col])).get(feature.name) ?? 'unmapped';
    const number: ColumnSchema = {name: 'Value', type: DG.COLUMN_TYPE.FLOAT};
    expect(single(DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'value', ['a', 'b']), number), 'unmapped');
    expect(single(DG.Column.fromList(DG.COLUMN_TYPE.INT, 'value', [1, 2]), number), 'value');
    expect(single(DG.Column.dateTime('value', 2), number), 'unmapped');
    expect(single(DG.Column.fromBigInt64Array('value', new BigInt64Array([1n, 2n])), number), 'unmapped');

    const molecule: ColumnSchema = {name: 'smiles', type: DG.COLUMN_TYPE.STRING, semType: 'Molecule'};
    const smiles = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'smiles', ['C', 'CC']);
    expect(single(smiles, molecule), 'unmapped');
    smiles.semType = 'molecule';
    expect(single(smiles, molecule), 'smiles');
  });

  test('compatibility texts', async () => {
    const x: ColumnSchema = {name: 'x', type: DG.COLUMN_TYPE.FLOAT};
    const message = (feature: ColumnSchema, col: DG.Column) => {
      const fit = compatibility(feature, col, 'XGBoost');
      return fit.kind === 'ok' ? 'ok' : `${fit.kind}: ${fit.message}`;
    };
    expect(message(x, DG.Column.fromBigInt64Array('big', new BigInt64Array([1n]))), 'error: \'big\' holds very ' +
      'large whole numbers, which XGBoost cannot read. Convert the column to a decimal type or choose another column.');
    expect(message(x, DG.Column.dateTime('when', 1)),
      'error: \'when\' holds dates but \'x\' needs numbers. Choose a column with numbers.');
    expect(message(x, DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'name', ['a'])),
      'error: \'name\' holds text but \'x\' needs numbers. Choose a column with numbers.');
    expect(message({name: 'label', type: DG.COLUMN_TYPE.STRING}, DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'v', [1])),
      'error: \'v\' holds numbers but \'label\' needs text. Choose a column with text.');
    const conc: ColumnSchema = {name: 'conc', type: DG.COLUMN_TYPE.FLOAT, semType: 'Concentration'};
    expect(message(conc, DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'conc', [1])),
      'hint: \'conc\' is not marked as Concentration; check that it is the same kind of value.');
    const molecule: ColumnSchema = {name: 'smiles', type: DG.COLUMN_TYPE.STRING, semType: 'Molecule'};
    expect(message(molecule, DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'smiles', ['C'])),
      'error: \'smiles\' is not a Molecule column, which \'smiles\' needs. Choose a Molecule column.');
    expect(message(x, DG.Column.fromList(DG.COLUMN_TYPE.INT, 'n', [1])), 'ok');
  });

  test('mappingProblems texts', async () => {
    const table = DG.DataFrame.fromColumns([DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'a', [1]),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'b', [2])]);
    const features: ColumnSchema[] = ['a', 'b', 'c', 'd'].map((name) => ({name, type: DG.COLUMN_TYPE.FLOAT}));
    const mapping = new Map([['a', 'a'], ['b', 'a'], ['c', 'gone']]);
    const problems = mappingProblems(features, mapping, table, 'XGBoost');
    expectArray(problems.map((p) => `${p.feature}${p.isUnmapped ? '!' : ''}: ${p.message}`), [
      'b: \'a\' is also used for \'a\'. Choose a different column.',
      'c: The column \'gone\' is no longer in the table. Choose another column for \'c\'.',
      'd!: Choose a column for \'d\'.',
    ]);
  });

  test('predictionName adds a number for every repeat, ignoring case', async () => {
    const table = DG.DataFrame.fromColumns([DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Species', ['a'])]);
    expect(predictionName(table, 'Species'), 'Species (predicted)');
    table.columns.addNewString('Species (predicted)');
    expect(predictionName(table, 'Species'), 'Species (predicted 2)');
    table.columns.addNewString('Species (predicted 2)');
    expect(predictionName(table, 'Species'), 'Species (predicted 3)');
    const other = DG.DataFrame.fromColumns([DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'species (PREDICTED)', ['a'])]);
    expect(predictionName(other, 'Species'), 'Species (predicted 2)');
  });

  test('reproduces the training predictions and records the application', async () => {
    const iris = await openIris();
    const {id, result} = await saveIrisModel(iris);
    const applications = () => forgeDb.applications.query().where('model_id', '=', id);
    try {
      const names = iris.columns.names();
      const sepalLength = rawValues(iris.getCol('Sepal.Length'));
      const model = await loadModel(id);
      const frame = featureFrame(iris, model.features, exactMapping(model.features, iris));
      expect(frame.getCol('Sepal.Length').getRawData().buffer === iris.getCol('Sepal.Length').getRawData().buffer,
        true, 'The feature frame copied a column');
      releaseFrame(frame);

      const [{column, skippedRows}, frames] = await framesSharing(iris.columns.toList(),
        () => applyAndRecord(requestFor(model, iris), 'api'));
      expectReleased(frames);
      expect(column.dataFrame?.dart === iris.dart, true, 'The prediction column belongs to another frame');
      expect(column.name, 'Species (predicted)');
      expect(column.type, DG.COLUMN_TYPE.STRING);
      expect(column.getTag(PREDICTION_TAG), id);
      expect(skippedRows, 0);
      expectMetrics(metricsOf('classification', iris.getCol('Species'), column, result.metrics.positiveClass),
        result.metrics.train, 1e-9);
      expectArray(column.toList(), await predictionsOf(model.engine, columnsOf(iris, MEASUREMENTS), result.blob));
      expectArray(iris.columns.names(), [...names, 'Species (predicted)']);
      expect(iris.getCol('col 1').type, DG.COLUMN_TYPE.INT);
      expectArray(rawValues(iris.getCol('Sepal.Length')), sepalLength);

      const [record] = await applications();
      expect(record.status, 'completed');
      expect(record.row_count, 150);
      expect(record.skipped_rows, 0);
      expect(record.column_name, 'Species (predicted)');
      expect(record.source, 'api');
      expect(record.table_name, iris.name);

      expect((await applyAndRecord(requestFor(model, iris), 'api')).column.name, 'Species (predicted 2)');
      expect((await applications()).length, 2);
    } finally {
      await deleteModel(id);
    }
    expect(await applications().count(), 0);
  }, {timeout: TIMEOUT});

  test('reproduces the training predictions of a regressor', async () => {
    const iris = await openIris();
    const {id, result} = await saveTestModel(columnsOf(iris, ['Sepal.Length', 'Sepal.Width', 'Petal.Width']),
      iris.getCol('Petal.Length'), iris.name);
    try {
      const model = await loadModel(id);
      const {column} = await applyAndRecord(requestFor(model, iris), 'api');
      expect(column.name, 'Petal.Length (predicted)');
      expect(column.type, DG.COLUMN_TYPE.FLOAT);
      expectMetrics(regressionMetrics(iris.getCol('Petal.Length'), column), result.metrics.train, 1e-6);

      // Numerical predictions of several batches and around a skipped row are put together from their raw data.
      const gapped = iris.clone();
      gapped.getCol('Sepal.Width').set(3, null);
      const batched = (await applyAndRecord(requestFor(model, gapped, {batchSize: 40}), 'api')).column;
      expect(batched.isNone(3), true, 'The skipped row has a prediction');
      for (let i = 0; i < iris.rowCount; i++) {
        if (i !== 3)
          expectFloat(batched.get(i), column.get(i), 1e-9, `Row ${i}`);
      }
    } finally {
      await deleteModel(id);
    }
  }, {timeout: TIMEOUT});

  test('batches give the same predictions, report progress and stop on cancel', async () => {
    const {model, reference} = shared();
    const percents: number[] = [];
    const progress: LoopProgress = {canceled: false, update: (percent) => percents.push(percent)};
    const {column} = await applyAndRecord(requestFor(model, irisCopy('batches'), {batchSize: 40}), 'api', progress);
    expectArray(column.toList(), reference);
    expect(percents.length, 4);
    expect(percents.every((p, i) => i === 0 || p > percents[i - 1]), true, percents.join(', '));
    expect(percents[percents.length - 1], 100);

    const cancelled = irisCopy('cancel');
    const columnCount = cancelled.columns.length;
    await expectExceptionAsync(async () => {
      await applyAndRecord(requestFor(model, cancelled, {batchSize: 40}), 'api', {canceled: true, update: () => {}});
    }, (e) => e instanceof ForgeError && e.message === 'Application was cancelled.');
    expect(cancelled.columns.length, columnCount);
    const records = await applicationsOf(cancelled);
    expect(records.length, 1);
    expect(records[0].status, 'cancelled');
  }, {timeout: TIMEOUT});

  test('a cancel from the event loop stops an application between batches', async () => {
    const {model, iris} = shared();
    // Iris four times, so 600 one-row batches run well past one pause interval.
    const repeated = (name: string) => DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, name,
      Array.from({length: 4 * iris.rowCount}, (_, i) => iris.get(name, i % iris.rowCount)));
    const table = DG.DataFrame.fromColumns(MEASUREMENTS.map(repeated));
    table.name = `${iris.name}-cancel-timer`;
    const columnCount = table.columns.length;
    let isScheduled = false;
    // The flag flips in a timer after the first batch, as a click on the task bar's cross does: only a batch loop that
    // yields to the event loop sees it before the last batch.
    const progress = {canceled: false, update: () => {
      if (!isScheduled) {
        isScheduled = true;
        setTimeout(() => {
          progress.canceled = true;
        }, 0);
      }
    }};
    await expectExceptionAsync(async () => {
      await applyAndRecord(requestFor(model, table, {batchSize: 1}), 'api', progress);
    }, (e) => e instanceof ForgeError && e.message === 'Application was cancelled.');
    expect(table.columns.length, columnCount);
    const records = await applicationsOf(table);
    expect(records.length, 1);
    expect(records[0].status, 'cancelled');
  }, {timeout: TIMEOUT});

  test('the API maps by exact names, the dialog prefill by close names', async () => {
    const {model, reference} = shared();
    const renamed = irisCopy('renamed');
    renamed.getCol('Sepal.Length').name = 'sepal length';
    await expectExceptionAsync(async () => {
      await applyAndRecord(requestFor(model, renamed), 'api');
    }, (e) => e instanceof ForgeError && e.message.includes('Choose a column for \'Sepal.Length\''));
    expect((await applicationsOf(renamed)).length, 0);

    const mapping = suggestMapping(model.features, renamed);
    const {column} = await applyAndRecord(requestFor(model, renamed, {mapping}), 'api');
    expectArray(column.toList(), reference);
    expect(renamed.col('sepal length') !== null && renamed.col('Sepal.Length') === null, true,
      renamed.columns.names().join(', '));
  }, {timeout: TIMEOUT});

  test('Skip rows leaves the rows with missing values without a prediction', async () => {
    const {model, reference} = shared();
    const table = irisCopy('skip');
    table.getCol('Sepal.Width').set(0, null);
    table.getCol('Sepal.Width').set(5, null);
    const {column, skippedRows} = await applyAndRecord(requestFor(model, table), 'api');
    expect(skippedRows, 2);
    for (let i = 0; i < table.rowCount; i++) {
      if (i === 0 || i === 5)
        expect(column.isNone(i), true, `Row ${i} has a prediction`);
      else
        expect(column.get(i), reference[i], `Row ${i}`);
    }
    const [record] = await applicationsOf(table);
    expect(record.skipped_rows, 2);
  }, {timeout: TIMEOUT});

  test('Impute predicts every row and leaves the table\'s gaps', async () => {
    const {model} = shared();
    const table = irisCopy('impute');
    table.getCol('Sepal.Width').set(0, null);
    table.getCol('Sepal.Width').set(5, null);
    const [{column, skippedRows}, frames] = await framesSharing(table.columns.toList(),
      () => applyAndRecord(requestFor(model, table, {missingValues: IMPUTE}), 'api'));
    expect(skippedRows, 0);
    expect(column.stats.missingValueCount, 0);
    expect(table.getCol('Sepal.Width').isNone(0) && table.getCol('Sepal.Width').isNone(5), true,
      'The table\'s missing values were filled');
    expectReleased(frames);
  }, {timeout: TIMEOUT});

  test('Impute leaves a row without feature values unpredicted and counts it', async () => {
    const {model, reference} = shared();
    const table = irisCopy('impute-empty-row');
    const emptyRow = 7;
    for (const name of MEASUREMENTS)
      table.getCol(name).set(emptyRow, null);
    table.getCol('Sepal.Width').set(0, null);
    const [{column, skippedRows}, frames] = await framesSharing(table.columns.toList(),
      () => applyAndRecord(requestFor(model, table, {missingValues: IMPUTE}), 'api'));
    expect(skippedRows, 1);
    expect(column.isNone(emptyRow), true, 'The row without feature values has a prediction');
    expect(column.isNone(0), false, 'The imputed row has no prediction');
    for (let i = 1; i < table.rowCount; i++) {
      if (i !== emptyRow)
        expect(column.get(i), reference[i], `Row ${i}`);
    }
    expectReleased(frames);
    expect((await applicationsOf(table))[0]?.skipped_rows, 1);
  }, {timeout: TIMEOUT});

  test('two applications at once on one table get their own names and release their frames', async () => {
    const {model, reference} = shared();
    const table = irisCopy('concurrent');
    const [results, frames] = await framesSharing(table.columns.toList(), () => Promise.all([
      applyAndRecord(requestFor(model, table), 'api'), applyAndRecord(requestFor(model, table), 'api')]));
    expectArray(results.map((r) => r.column.name).sort(), ['Species (predicted)', 'Species (predicted 2)'].sort());
    for (const {column} of results)
      expectArray(column.toList(), reference);
    expectReleased(frames);
    expect((await applicationsOf(table)).length, 2);
  }, {timeout: TIMEOUT});

  test('a one-hot model is applied end to end and releases its frames', async () => {
    const rows = 40;
    const x = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'x', valuesOf(rows, (i) => i / 4));
    const color = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'color',
      Array.from({length: rows}, (_, i) => ['red', 'green', 'blue'][i % 3]));
    const y = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'y', valuesOf(rows, (i) => i / 4 + 5 * (i % 3)));
    // The columns the replay builds from 'color', so the method learns what it will be given.
    const encoded = color.categories.map((category) => DG.Column.fromInt32Array(`color=${category}`,
      Int32Array.from({length: rows}, (_, i) => color.get(i) === category ? 1 : 0)));
    const features = [{name: 'x', type: DG.COLUMN_TYPE.FLOAT}, {name: 'color', type: DG.COLUMN_TYPE.STRING}];
    const {id, result} = await saveTestModel([x, ...encoded], y, `forge-test-one-hot-${Date.now()}`,
      {features: {columns: features}, options: {preprocessingInfo: [ONE_HOT], postprocessingInfo: []}});
    try {
      const model = await loadModel(id);
      const expected = await predictionsOf(model.engine, [x, ...encoded], result.blob);
      const table = DG.DataFrame.fromColumns([x, color]);
      const [{column}, frames] = await framesSharing(table.columns.toList(),
        () => applyAndRecord(requestFor(model, table), 'api'));
      expect(column.name, 'y (predicted)');
      expectArray(column.toList(), expected);
      expectReleased(frames);
      expectArray(table.columns.names(), ['x', 'color', 'y (predicted)']);
      expect(color.dataFrame?.dart === table.dart, true, 'The text column left the table');
    } finally {
      await deleteModel(id);
    }
  }, {timeout: TIMEOUT});

  test('a table whose every row has a missing value gets no prediction', async () => {
    const {model} = shared();
    const table = irisCopy('all-missing');
    const sepalWidth = table.getCol('Sepal.Width');
    for (let i = 0; i < table.rowCount; i++)
      sepalWidth.set(i, null, false);
    const columnCount = table.columns.length;
    await expectExceptionAsync(async () => {
      await applyAndRecord(requestFor(model, table), 'api');
    }, (e) => e instanceof ForgeError && e.message === 'Every row has a missing value in the columns the model ' +
      'needs, so nothing can be predicted. Fill the missing values or choose Impute.');
    expect(table.columns.length, columnCount);
  }, {timeout: TIMEOUT});

  test('a refused mapping is not recorded, a failing method call is', async () => {
    const {model} = shared();
    const table = irisCopy('refused');
    const mapping = exactMapping(model.features, table);
    mapping.set('Sepal.Length', 'Species');
    await expectExceptionAsync(async () => {
      await applyAndRecord(requestFor(model, table, {mapping}), 'api');
    }, throwsWith('\'Species\' holds text but \'Sepal.Length\' needs numbers.'));
    expect((await applicationsOf(table)).length, 0);

    const folder = `${BLOB_ROOT}/${DG.Utils.uuid4()}`;
    await grok.dapi.files.write(`${folder}/model.bin`, new Uint8Array([1, 2, 3]));
    let id: string | undefined;
    try {
      [{id}] = await forgeDb.models.insert({...XGBOOST_FIELDS, name: `forge-test-model-${Date.now()}`,
        task: 'classification', target_name: 'Species', target: {name: 'Species', type: DG.COLUMN_TYPE.STRING},
        features: {columns: IRIS_FEATURES}, storage_mode: 'none', blob: `file://${folder}/model.bin`});
      const broken = await loadModel(id);
      await expectExceptionAsync(async () => {
        await applyAndRecord(requestFor(broken, table), 'api');
      });
      const records = await forgeDb.applications.query().where('model_id', '=', id);
      expect(records.length, 1);
      expect(records[0].status, 'failed');
      expect((records[0].error ?? '') !== '', true, 'The failed application has no error');
    } finally {
      if (id !== undefined)
        await deleteModel(id);
      else
        await grok.dapi.files.delete(folder);
    }
  }, {timeout: TIMEOUT});

  test('loadModel refuses models it cannot apply', async () => {
    const stamp = Date.now();
    const ownBlob = () => `file://${BLOB_ROOT}/${DG.Utils.uuid4()}/model.bin`;
    const features = {columns: IRIS_FEATURES};
    const cases: [Partial<ModelInsert>, string][] = [
      [{features}, `The model 'forge-test-model-${stamp}' has no model file, so it cannot be applied. ` +
        'It was not trained and saved in Forge.'],
      [{features, blob: 'file://System:DemoFiles/iris.csv'},
        `The model 'forge-test-model-${stamp}' points to a file outside Forge's storage and cannot be applied.`],
      [{blob: ownBlob()}, `The model 'forge-test-model-${stamp}' has no feature list.`],
      [{features, blob: ownBlob(), engine_name: 'forge-test-engine'},
        'The method \'forge-test-engine\' is not installed. Install the Eda package.'],
    ];
    for (const [fields, text] of cases) {
      const [{id}] = await forgeDb.models.insert({...XGBOOST_FIELDS, name: `forge-test-model-${stamp}`,
        task: 'classification', target_name: 'Species', storage_mode: 'none', ...fields});
      try {
        await expectExceptionAsync(() => loadModel(id).then(() => {}),
          (e) => e instanceof ForgeError && e.message === text);
      } finally {
        await forgeDb.models.delete(id);
      }
    }
  });

  test('loadModel finds a model by id or by a unique name', async () => {
    const stamp = Date.now();
    const name = `forge-test-model-${stamp}`;
    const insert = async (modelName: string) => (await forgeDb.models.insert({...XGBOOST_FIELDS, name: modelName,
      task: 'regression', target_name: 'y', storage_mode: 'none'}))[0].id;
    // The rows have no model file, so a row that is found is refused with its own name.
    const isFound = (e: unknown) => e instanceof ForgeError && e.message.startsWith(`The model '${name}' has no model`);
    const ids: string[] = [];
    try {
      ids.push(await insert(name));
      ids.push(await insert(`forge-test-dup-${stamp}`), await insert(`forge-test-dup-${stamp}`));
      await expectExceptionAsync(() => loadModel(ids[0]).then(() => {}), isFound);
      await expectExceptionAsync(() => loadModel(name).then(() => {}), isFound);
      await expectExceptionAsync(() => loadModel(`forge-test-none-${stamp}`).then(() => {}),
        (e) => e instanceof ForgeError &&
          e.message === `No Forge model 'forge-test-none-${stamp}' is available to you.`);
      await expectExceptionAsync(() => loadModel(`forge-test-dup-${stamp}`).then(() => {}),
        (e) => e instanceof ForgeError &&
          e.message === `Several models are named 'forge-test-dup-${stamp}'. Use the model id.`);
    } finally {
      for (const id of ids)
        await forgeDb.models.delete(id);
    }
  });

  test('Forge:applyModel applies by exact names and the given pairs', async () => {
    const {id} = shared();
    const call = (table: DG.DataFrame, columnNamesMap: {[feature: string]: string}) =>
      grok.functions.call('Forge:applyModel', {model: id, table, columnNamesMap, showProgress: false});

    const exact = irisCopy('api');
    await call(exact, {});
    expect(exact.col('Species (predicted)') !== null, true, exact.columns.names().join(', '));

    const renamed = irisCopy('api-renamed');
    renamed.getCol('Sepal.Length').name = 'sepal length';
    await call(renamed, {'Sepal.Length': 'sepal length'});
    expect(renamed.col('Species (predicted)') !== null, true, renamed.columns.names().join(', '));

    const unmapped = irisCopy('api-unmapped');
    unmapped.getCol('Sepal.Length').name = 'sepal length';
    await expectExceptionAsync(() => call(unmapped, {}),
      throwsWith('The table has no column \'Sepal.Length\' the model needs. Map it in columnNamesMap.'));

    const dated = irisCopy('api-dated');
    dated.columns.add(DG.Column.dateTime('when', dated.rowCount));
    await expectExceptionAsync(() => call(dated, {'Sepal.Length': 'when'}),
      throwsWith('\'when\' holds dates but \'Sepal.Length\' needs numbers. Choose a column with numbers.'));

    for (const table of [exact, renamed])
      expectArray((await applicationsOf(table)).map((r) => r.source), ['api']);
    for (const table of [unmapped, dated])
      expect((await applicationsOf(table)).length, 0);
  }, {timeout: TIMEOUT});

  test('an integer feature with missing values never reaches the method', async () => {
    const rows = 40;
    const nullRows = [5, 17];
    const k = DG.Column.fromList(DG.COLUMN_TYPE.INT, 'k', valuesOf(rows, (i) => i % 7, nullRows));
    const x = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'x', valuesOf(rows, (i) => i / 2));
    const y = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'y', valuesOf(rows, (i) => 2 * (i % 7) + i / 2 + (i % 3) / 10));
    const table = DG.DataFrame.fromColumns([k, x, y]);
    const original = table.columns.toList().map(rawValues);

    const request = await prepareTraining(selectionOf(columnsOf(table, ['k', 'x']), y));
    expect(request.options.missingValues?.skippedRows, 2);
    const metrics = (await trainModel(request)).metrics;
    releaseFrame(request.features);

    const kept = DG.BitSet.create(rows, (i) => !nullRows.includes(i));
    const byHand = await requestOf([k.clone(kept), x.clone(kept)], y.clone(kept));
    const expected = (await trainModel(byHand)).metrics;
    expectMetrics(metrics.train, expected.train, 1e-9);
    expectMetrics(metrics.validation, expected.validation, 1e-9);

    expectArray(table.columns.names(), ['k', 'x', 'y']);
    table.columns.toList().forEach((col, i) => expectArray(rawValues(col), original[i]));
    expect(k.isNone(5) && k.isNone(17), true, 'The gaps in k were filled');
  }, {timeout: TIMEOUT});
});

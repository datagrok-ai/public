import * as DG from 'datagrok-api/dg';
import {category, expect, expectArray, expectExceptionAsync, test} from '@datagrok-libraries/test/src/test';
import {defaultValuesOf} from '../engines/engine';
import {ForgeError} from '../forge-error';
import {imputeFunction, imputeSettingsOf, prepareMissingValues} from '../preparation/missing-values';
import {replayPostprocessing, replayPreprocessing} from '../preparation/pipeline';
import {preparationOptionsOf} from '../preparation/preparation-options';
import {columnNamed, expectReleased, framesSharing, IMPUTE, names, NO_OPTIONS, valuesOf} from './test-data';

const gapCount = (columns: DG.Column[]) => columns.reduce((sum, c) => sum + c.stats.missingValueCount, 0);

category('Preparation', () => {
  test('one-hot replaces the text and boolean columns', async () => {
    const columns = [
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'x', [0.5, 1.5, 2.5]),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'c', ['a', 'b', 'a']),
      DG.Column.fromList(DG.COLUMN_TYPE.BOOL, 'f', [true, false, true]),
    ];
    const encoded = replayPreprocessing(columns, {...NO_OPTIONS, preprocessingInfo: ['one-hot']});
    expectArray(names(encoded), ['x', 'c=a', 'c=b', 'f=false', 'f=true']);
    expectArray(columnNamed(encoded, 'c=a').toList(), [1, 0, 1]);
    expectArray(columnNamed(encoded, 'c=b').toList(), [0, 1, 0]);
    expectArray(columnNamed(encoded, 'f=false').toList(), [0, 1, 0]);
    expectArray(columnNamed(encoded, 'f=true').toList(), [1, 0, 1]);
    expectArray(names(columns), ['x', 'c', 'f']);
  });

  test('skip-unique-categories removes the all-unique text columns', async () => {
    const columns = [
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'id', ['a', 'b', 'c']),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'group', ['x', 'x', 'y']),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'v', [1, 2, 3]),
    ];
    const kept = replayPreprocessing(columns, {...NO_OPTIONS, preprocessingInfo: ['skip-unique-categories']});
    expectArray(names(kept), ['group', 'v']);
    // A record of the skipped names replaces the rule: a repeating id is dropped, a unique group is kept.
    const repeated = [
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'id', ['a', 'a', 'b']),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'group', ['x', 'y', 'z']),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'v', [1, 2, 3]),
    ];
    const recorded = {...NO_OPTIONS, preprocessingInfo: ['skip-unique-categories'], skippedColumns: ['id']};
    expectArray(names(replayPreprocessing(repeated, recorded)), ['group', 'v']);
    expectArray(names(replayPreprocessing(repeated, {...recorded, skippedColumns: undefined})), ['id', 'v']);
  });

  test('binary-classification maps scores to the two classes', async () => {
    const scores = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'y', [0.2, 0.5, 0.9]);
    const classes = replayPostprocessing(scores, {...NO_OPTIONS, postprocessingInfo: ['binary-classification'],
      positiveClass: 'yes', negativeClass: 'no', binaryClassificationThreshold: 0.5, targetType: 'string'});
    expectArray(classes.toList(), ['no', 'yes', 'yes']);
  });

  test('empty step lists leave the data untouched', async () => {
    const columns = [DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'c', ['a', 'b'])];
    const prediction = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'y', [0.2, 0.9]);
    expect(replayPreprocessing(columns, NO_OPTIONS) === columns, true, 'New columns without steps');
    expect(replayPostprocessing(prediction, NO_OPTIONS) === prediction, true, 'A new column without steps');
  });

  test('missing-value ids are skipped, unknown ids refused', async () => {
    const columns = [DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'c', ['a', 'b'])];
    const recorded = {...NO_OPTIONS, preprocessingInfo: ['impute-missing', 'ignore-missing']};
    expect(replayPreprocessing(columns, recorded) === columns, true, 'The recorded ids were replayed');
    await expectExceptionAsync(async () => {
      replayPreprocessing(columns, {...NO_OPTIONS, preprocessingInfo: ['one-hot-x']});
    }, (e) => e instanceof ForgeError &&
      e.message === 'The model uses the preparation step \'one-hot-x\', which Forge cannot replay yet.');
  });

  test('preparationOptionsOf tolerates missing and malformed options', async () => {
    for (const value of [{}, null, undefined, 'x', [1]])
      expect(JSON.stringify(preparationOptionsOf(value)), JSON.stringify(NO_OPTIONS), JSON.stringify(value));
    const options = preparationOptionsOf({preprocessingInfo: ['one-hot', 3], positiveClass: 'a', targetType: 7,
      binaryClassificationThreshold: 0.4, missingValues: {mode: 'impute', neighbors: 4, skippedRows: 2},
      oneHotCategories: {c: ['x', 1, 'y'], d: 'z'}, skippedColumns: ['id', 2]});
    expect(JSON.stringify(options), JSON.stringify({preprocessingInfo: ['one-hot'], postprocessingInfo: [],
      missingValues: {mode: 'impute', skippedRows: 2, neighbors: 4},
      oneHotCategories: {c: ['x', 'y'], d: []}, skippedColumns: ['id'],
      positiveClass: 'a', binaryClassificationThreshold: 0.4}));
  });

  test('Skip rows drops the rows with a missing feature or target value, without touching the input', async () => {
    const rows = 12;
    const a = DG.Column.fromList(DG.COLUMN_TYPE.INT, 'a', valuesOf(rows, (i) => i, [3]));
    const b = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'b', valuesOf(rows, (i) => i / 2, [7]));
    const y = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'y', valuesOf(rows, (i) => i % 3, [9]).map((v) =>
      v === null ? null : `class ${v}`));
    const prepared = await prepareMissingValues([a, b], y, {mode: 'skip'});
    expect(prepared.keptRows?.trueCount, 9);
    expect(prepared.skippedRows, 3);
    expectArray(prepared.features.map((c) => c.length), [9, 9]);
    expect(prepared.target?.length, 9);
    expect(gapCount(prepared.features), 0);
    expect(a.isNone(3) && b.isNone(7) && y.isNone(9), true, 'The input columns changed');
    expect(a.length, rows);

    const full = [DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'x', valuesOf(rows, (i) => i))];
    const untouched = await prepareMissingValues(full, undefined, {mode: 'skip'});
    expect(untouched.keptRows === null, true, 'Rows kept without gaps');
    expect(untouched.features === full, true, 'The columns were copied without gaps');
  });

  test('Skip rows without a target counts the feature gaps only', async () => {
    const rows = 12;
    const columns = [
      DG.Column.fromList(DG.COLUMN_TYPE.INT, 'a', valuesOf(rows, (i) => i, [3])),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'b', valuesOf(rows, (i) => i / 2, [3, 7])),
    ];
    const prepared = await prepareMissingValues(columns, undefined, {mode: 'skip'});
    expect(prepared.skippedRows, 2);
    expect(prepared.target === undefined, true, 'A target appeared');
    const skipped = Array.from({length: rows}, (_, i) => i).filter((i) => !(prepared.keptRows?.get(i) ?? true));
    expectArray(skipped, [3, 7]);
  });

  test('a copied text column keeps only the categories its rows have', async () => {
    const rows = 12;
    const c = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'c',
      Array.from({length: rows}, (_, i) => i === 3 ? null : i % 3 === 0 ? 'x' : 'y'));
    const columns = [
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'a', valuesOf(rows, (i) => i)),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'b', valuesOf(rows, (i) => 2 * i)),
      c,
    ];
    const skipped = await prepareMissingValues(columns, undefined, {mode: 'skip'});
    expectArray(columnNamed(skipped.features, 'c').categories, ['x', 'y']);
    const imputed = await prepareMissingValues(columns, undefined, IMPUTE);
    expect(gapCount(imputed.features), 0);
    const imputedCategories = columnNamed(imputed.features, 'c').categories;
    expect(imputedCategories.includes(''), false, imputedCategories.join(', '));
    expect(c.categories.includes(''), true, 'The input column lost its empty category');
  }, {timeout: 30000});

  test('Impute fills the gaps in copies of the gapped columns and skips the rows it cannot fill', async () => {
    const func = imputeFunction();
    if (func === undefined)
      throw new Error('Eda:knnImpute is not available');
    expect(JSON.stringify(defaultValuesOf(imputeSettingsOf(func))),
      JSON.stringify({neighbors: 4, distance: 'Euclidean'}));

    const rows = 30;
    const columnsOf = (allMissing: number[]) => [
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'a', valuesOf(rows, (i) => i, [3, 10, ...allMissing])),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'b', valuesOf(rows, (i) => 2 * i, [7, ...allMissing])),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'c', valuesOf(rows, (i) => rows - i, allMissing)),
    ];
    const [a, b, c] = columnsOf([]);
    const prepared = await prepareMissingValues([a, b, c], undefined, IMPUTE);
    expect(gapCount(prepared.features), 0);
    expectArray(prepared.imputedColumns, ['a', 'b']);
    expect(prepared.skippedRows, 0);
    expect(prepared.keptRows === null, true, 'Rows were skipped');
    expect(a.isNone(3) && a.isNone(10) && b.isNone(7), true, 'The input columns were imputed');
    expect(columnNamed(prepared.features, 'c') === c, true, 'The column without gaps was copied');

    // A yes/no column without gaps stays shared in the imputation frame, which is given back.
    const gapped = columnsOf([12]);
    gapped.push(DG.Column.fromList(DG.COLUMN_TYPE.BOOL, 'flag', Array.from({length: rows}, (_, i) => i % 2 === 1)));
    const [withEmptyRow, frames] = await framesSharing(gapped,
      () => prepareMissingValues(gapped, undefined, IMPUTE));
    expect(withEmptyRow.failedRows, 1);
    expect(withEmptyRow.skippedRows, 1);
    expectArray(withEmptyRow.features.map((col) => col.length), [rows - 1, rows - 1, rows - 1, rows - 1]);
    expect(withEmptyRow.keptRows?.get(12), false);
    expect(gapCount(withEmptyRow.features), 0);
    expectReleased(frames);
  }, {timeout: 30000});
});

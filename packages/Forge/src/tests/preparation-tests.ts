import * as DG from 'datagrok-api/dg';
import {category, expect, expectArray, expectExceptionAsync, test} from '@datagrok-libraries/test/src/test';
import {defaultValuesOf} from '../engines/engine';
import {ForgeError} from '../forge-error';
import {imputeFunction, imputeSettingsOf, prepareMissingValues} from '../preparation/missing-values';
import {PreparationOptions, preparationOptionsOf} from '../preparation/preparation-options';
import {replayPostprocessing, replayPreprocessing} from '../preparation/preparation-steps';
import {releaseFrame} from '../preparation/shared-frame';
import {expectReleased, framesSharing, IMPUTE, valuesOf} from './test-data';

const NO_STEPS: PreparationOptions = {preprocessingInfo: [], postprocessingInfo: []};

const gapCount = (frame: DG.DataFrame) => frame.columns.toList().reduce((sum, c) => sum + c.stats.missingValueCount, 0);

category('Preparation', () => {
  test('one-hot replaces the text and boolean columns', async () => {
    const frame = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'x', [0.5, 1.5, 2.5]),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'c', ['a', 'b', 'a']),
      DG.Column.fromList(DG.COLUMN_TYPE.BOOL, 'f', [true, false, true]),
    ]);
    const encoded = replayPreprocessing(frame, {...NO_STEPS, preprocessingInfo: ['one-hot']});
    expectArray(encoded.columns.names(), ['x', 'c=a', 'c=b', 'f=false', 'f=true']);
    expectArray(encoded.getCol('c=a').toList(), [1, 0, 1]);
    expectArray(encoded.getCol('c=b').toList(), [0, 1, 0]);
    expectArray(encoded.getCol('f=false').toList(), [0, 1, 0]);
    expectArray(encoded.getCol('f=true').toList(), [1, 0, 1]);
    expectArray(frame.columns.names(), ['x', 'c', 'f']);
  });

  test('skip-unique-categories removes the all-unique text columns', async () => {
    const frame = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'id', ['a', 'b', 'c']),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'group', ['x', 'x', 'y']),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'v', [1, 2, 3]),
    ]);
    const kept = replayPreprocessing(frame, {...NO_STEPS, preprocessingInfo: ['skip-unique-categories']});
    expectArray(kept.columns.names(), ['group', 'v']);
  });

  test('binary-classification maps scores to the two classes', async () => {
    const scores = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'y', [0.2, 0.5, 0.9]);
    const classes = replayPostprocessing(scores, {...NO_STEPS, postprocessingInfo: ['binary-classification'],
      positiveClass: 'yes', negativeClass: 'no', binaryClassificationThreshold: 0.5, targetType: 'string'});
    expectArray(classes.toList(), ['no', 'yes', 'yes']);
  });

  test('empty step lists leave the data untouched', async () => {
    const frame = DG.DataFrame.fromColumns([DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'c', ['a', 'b'])]);
    const prediction = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'y', [0.2, 0.9]);
    expect(replayPreprocessing(frame, NO_STEPS) === frame, true, 'A new frame without steps');
    expect(replayPostprocessing(prediction, NO_STEPS) === prediction, true, 'A new column without steps');
  });

  test('missing-value ids are skipped, unknown ids refused', async () => {
    const frame = DG.DataFrame.fromColumns([DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'c', ['a', 'b'])]);
    const recorded = {...NO_STEPS, preprocessingInfo: ['impute-missing', 'ignore-missing']};
    expect(replayPreprocessing(frame, recorded) === frame, true, 'The recorded ids were replayed');
    await expectExceptionAsync(async () => {
      replayPreprocessing(frame, {...NO_STEPS, preprocessingInfo: ['one-hot-x']});
    }, (e) => e instanceof ForgeError &&
      e.message === 'The model uses the preparation step \'one-hot-x\', which Forge cannot replay yet.');
  });

  test('preparationOptionsOf tolerates missing and malformed options', async () => {
    for (const value of [{}, null, undefined, 'x', [1]])
      expect(JSON.stringify(preparationOptionsOf(value)), JSON.stringify(NO_STEPS), JSON.stringify(value));
    const options = preparationOptionsOf({preprocessingInfo: ['one-hot', 3], positiveClass: 'a', targetType: 7,
      binaryClassificationThreshold: 0.4, missingValues: {mode: 'impute', neighbors: 4, skippedRows: 2}});
    expect(JSON.stringify(options), JSON.stringify({preprocessingInfo: ['one-hot'], postprocessingInfo: [],
      missingValues: {mode: 'impute', skippedRows: 2, neighbors: 4}, positiveClass: 'a',
      binaryClassificationThreshold: 0.4}));
  });

  test('Skip rows drops the rows with a missing feature or target value, without touching the input', async () => {
    const rows = 12;
    const a = DG.Column.fromList(DG.COLUMN_TYPE.INT, 'a', valuesOf(rows, (i) => i, [3]));
    const b = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'b', valuesOf(rows, (i) => i / 2, [7]));
    const y = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'y', valuesOf(rows, (i) => i % 3, [9]).map((v) =>
      v === null ? null : `class ${v}`));
    const frame = DG.DataFrame.fromColumns([a, b]);
    const prepared = await prepareMissingValues(frame, y, {mode: 'skip'});
    expect(prepared.keptRows?.trueCount, 9);
    expect(prepared.skippedRows, 3);
    expect(prepared.features.rowCount, 9);
    expect(prepared.target?.length, 9);
    expect(gapCount(prepared.features), 0);
    expect(a.isNone(3) && b.isNone(7) && y.isNone(9), true, 'The input columns changed');
    expect(frame.rowCount, rows);

    const full = DG.DataFrame.fromColumns([DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'x', valuesOf(rows, (i) => i))]);
    const untouched = await prepareMissingValues(full, undefined, {mode: 'skip'});
    expect(untouched.keptRows === null, true, 'Rows kept without gaps');
    expect(untouched.features === full, true, 'The frame was copied without gaps');
  });

  test('Skip rows without a target counts the feature gaps only', async () => {
    const rows = 12;
    const frame = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.INT, 'a', valuesOf(rows, (i) => i, [3])),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'b', valuesOf(rows, (i) => i / 2, [3, 7])),
    ]);
    const prepared = await prepareMissingValues(frame, undefined, {mode: 'skip'});
    expect(prepared.skippedRows, 2);
    expect(prepared.target === undefined, true, 'A target appeared');
    const skipped = Array.from({length: rows}, (_, i) => i).filter((i) => !(prepared.keptRows?.get(i) ?? true));
    expectArray(skipped, [3, 7]);
  });

  test('a copied text column keeps only the categories its rows have', async () => {
    const rows = 12;
    const c = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'c',
      Array.from({length: rows}, (_, i) => i === 3 ? null : i % 3 === 0 ? 'x' : 'y'));
    const frame = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'a', valuesOf(rows, (i) => i)),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'b', valuesOf(rows, (i) => 2 * i)),
      c,
    ]);
    const skipped = await prepareMissingValues(frame, undefined, {mode: 'skip'});
    expectArray(skipped.features.getCol('c').categories, ['x', 'y']);
    const imputed = await prepareMissingValues(frame, undefined, IMPUTE);
    expect(gapCount(imputed.features), 0);
    expect(imputed.features.getCol('c').categories.includes(''), false,
      imputed.features.getCol('c').categories.join(', '));
    releaseFrame(imputed.features);
    expect(c.categories.includes(''), true, 'The input column lost its empty category');
  }, {timeout: 30000});

  test('Impute fills the gaps in copies of the gapped columns and skips the rows it cannot fill', async () => {
    const func = imputeFunction();
    if (func === undefined)
      throw new Error('Eda:knnImpute is not available');
    expect(JSON.stringify(defaultValuesOf(imputeSettingsOf(func))),
      JSON.stringify({neighbors: 4, distance: 'Euclidean'}));

    const rows = 30;
    const frameOf = (allMissing: number[]) => DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'a', valuesOf(rows, (i) => i, [3, 10, ...allMissing])),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'b', valuesOf(rows, (i) => 2 * i, [7, ...allMissing])),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'c', valuesOf(rows, (i) => rows - i, allMissing)),
    ]);
    const frame = frameOf([]);
    const prepared = await prepareMissingValues(frame, undefined, IMPUTE);
    expect(gapCount(prepared.features), 0);
    expectArray(prepared.imputedColumns, ['a', 'b']);
    expect(prepared.skippedRows, 0);
    expect(prepared.keptRows === null, true, 'Rows were skipped');
    expect(frame.getCol('a').isNone(3) && frame.getCol('a').isNone(10) && frame.getCol('b').isNone(7), true,
      'The input columns were imputed');
    expect(prepared.features.getCol('c').getRawData().buffer === frame.getCol('c').getRawData().buffer, true,
      'The column without gaps was copied');
    releaseFrame(prepared.features);

    // A yes/no column without gaps stays shared in the imputation frame, which is given back when rows are skipped.
    const gapped = frameOf([12]);
    const flags = Array.from({length: rows}, (_, i) => i % 2 === 1);
    gapped.columns.add(DG.Column.fromList(DG.COLUMN_TYPE.BOOL, 'flag', flags));
    const [withEmptyRow, frames] = await framesSharing(gapped.columns.toList(),
      () => prepareMissingValues(gapped, undefined, IMPUTE));
    expect(withEmptyRow.failedRows, 1);
    expect(withEmptyRow.skippedRows, 1);
    expect(withEmptyRow.features.rowCount, rows - 1);
    expect(withEmptyRow.keptRows?.get(12), false);
    expect(gapCount(withEmptyRow.features), 0);
    expectReleased(frames);
  }, {timeout: 30000});
});

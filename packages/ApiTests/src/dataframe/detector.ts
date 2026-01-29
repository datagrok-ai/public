import * as DG from 'datagrok-api/dg';
import {category, test} from '@datagrok-libraries/test/src/test';

category('Detector', () => {
  test('SingleEmptyCategory', async () => {
    const df: DG.DataFrame = DG.DataFrame.fromObjects([
      {col1: null},
      {col1: null},
      {col1: null},
      {col1: null},
      {col1: null},
    ])!;
    const col: DG.Column = df.col('col1')!;
    const res: boolean = DG.Detector.sampleCategories(col, (v) => false, 1);
    if (res)
      throw new Error('DG.Detector.sampleCategories() always returns true on single empty category.');
  });
});

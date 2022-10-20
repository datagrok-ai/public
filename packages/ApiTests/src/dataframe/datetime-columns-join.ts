import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/utils/src/test';


category('DataFrame', () => {
  test('Join DF with DF (with datetime columns) should be without errors', async () => {
    const df1 = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.TYPE.INT, 'A', [1, 2, 3]),
      DG.Column.fromList(DG.TYPE.INT, 'B', [4, 5, 6]),
      DG.Column.fromList(DG.TYPE.INT, 'C', [7, 8, 9]),
    ]);

    const df2 = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.TYPE.INT, 'AA', [1, 2, 3]),
      DG.Column.fromList(DG.TYPE.INT, 'BB', [4, 5, 6]),
      DG.Column.fromList(DG.TYPE.INT, 'CC', [7, 8, 9]),
    ]);

    df1.columns.addNew('D', DG.COLUMN_TYPE.DATE_TIME);
    df1.set('D', 0, ''); // <- will be a problem on join
    df1.set('D', 1, '2017-10-18T00:00:00.000Z');
    df1.set('D', 2, '2017-11-18T00:00:00.000Z');
    
    const result = DG.DataFrame.fromColumns([
      df2.columns.byName('AA'),
      df2.columns.byName('BB'),
      df2.columns.byName('CC'),
    ]).join(df1, ['AA'], ['A'], ['AA', 'BB', 'CC'], ['B', 'C', 'D'], 'inner', false);

    expect(result != null, true);
  });
});

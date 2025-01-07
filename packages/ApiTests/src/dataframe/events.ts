import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {hashDataFrame} from '@datagrok-libraries/utils/src/dataframe-utils';
import {category, expect, test, expectArray} from '@datagrok-libraries/utils/src/test';
import dayjs from 'dayjs';

category('DataFrame', () => {
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.TYPE.FLOAT, 'x', [1, 2, 3]),
    DG.Column.fromList(DG.TYPE.FLOAT, 'y', [4, 5, 6]),
    DG.Column.fromList(DG.TYPE.FLOAT, 'z', [7, 8, 9]),
  ]);

  test('raise an event on set', async () => {
    let dataChanged = false;
    df.onDataChanged.subscribe(() => {
      dataChanged = true;
    });
    df.set('x', 1, 0);
    expect(dataChanged, true);
  });

  test('raise an event on setValues', async () => {
    let dataChanged = false;
    df.onDataChanged.subscribe(() => {
      dataChanged = true;
    });
    df.rows.setValues(1, [0, 0, 0]);
    expect(dataChanged, true);
  });

  test('hash', async () => {
    const df = grok.data.demo.demog(10);
    expect(hashDataFrame(df).length, 32);

    const df1 = DG.DataFrame.fromCsv(`a,b\n1,0\n2,0\n3,0`);
    const df2 = DG.DataFrame.fromCsv(`a,b\n2,0\n1,0\n3,0`);
    expectArray(hashDataFrame(df1), hashDataFrame(df2));

    const df3 = DG.DataFrame.fromCsv(`a,b\n"abc",0\n"dce",0\n"xyz",0`);
    const df4 = DG.DataFrame.fromCsv(`a,b\n"dce",0\n"abc",0\n"xyz",0`);
    expectArray(hashDataFrame(df3), hashDataFrame(df4));
  });

  test('datetime column', async () => {
    const t = grok.data.testData('demog');
    const c = t.columns.byName('started');
    c.set(1, dayjs.utc('2022-01-01'));
    expect(c.get(1).valueOf(), 1640995200000);
    c.set(1, null);
    expect(c.get(1), null);
  });
}, {owner: 'aparamonov@datagrok.ai'});

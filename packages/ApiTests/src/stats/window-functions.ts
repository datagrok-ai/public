import * as DG from 'datagrok-api/dg';
import {category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';

category('Stats: window functions', () => {
  const t = DG.DataFrame.fromCsv('region,date,amount\na,3,10\nb,1,20\na,2,\nb,4,40');
  const amount = t.getCol('amount');
  const list = (c: DG.Column) => Array.from({length: c.length}, (_, i) => c.isNone(i) ? null : c.get(i));

  test('cumSum', async () => {
    expectArray(list(amount.stats.cumSum()), [10, 30, null, 70]);
    expect(amount.stats.cumSum().type, DG.COLUMN_TYPE.FLOAT);
    expect(amount.stats.cumSum().length, 4);
  });

  test('cumSum: by', async () => {
    expectArray(list(amount.stats.cumSum({by: ['region']})), [10, 20, null, 60]);
  });

  test('cumSum: orderBy and ascending', async () => {
    expectArray(list(amount.stats.cumSum({orderBy: ['date']})), [30, 20, null, 70]);
    expectArray(list(amount.stats.cumSum({orderBy: [t.getCol('date')], ascending: [false]})), [50, 70, null, 40]);
  });

  test('cumSum: explicit order', async () => {
    expectArray(list(amount.stats.cumSum({order: Int32Array.from([3, 2, 1, 0])})), [70, 60, null, 40]);
    expectArray(list(amount.stats.cumSum({order: [3, 2, 1, 0]})), [70, 60, null, 40]);
  });

  test('cumSum: mask', async () => {
    const mask = DG.BitSet.create(4, (i) => i === 0 || i === 3);
    expectArray(list(DG.Stats.fromColumn(amount, mask).cumSum()), [10, null, null, 50]);
  });

  test('movingAvg', async () => {
    expectArray(list(amount.stats.movingAvg(2)), [10, 15, 20, 40]);
    expectArray(list(amount.stats.movingAvg(2, {minPeriods: 2})), [null, 15, null, null]);
    expectArray(list(amount.stats.movingAvg(3, {by: ['region']})), [10, 20, 10, 30]);
  });
}, {owner: 'askalkin@datagrok.ai'});

import * as DG from 'datagrok-api/dg';
// import * as grok from 'datagrok-api/grok';
import {category, expect, test} from '@datagrok-libraries/utils/src/test';

category('Stats: fromColumn', () => {
  const t = DG.DataFrame.create(3);
  t.columns.add(DG.Column.fromInt32Array('number', Int32Array.from([12, 10, 15])));
  t.columns.add(DG.Column.fromInt32Array('population', Int32Array.from([3, 4, 5])));

  const stats = DG.Stats.fromColumn(t.getCol('population'));

  test('avg', async () => {
    expect(stats.avg, 4);
  });

  test('kurt', async () => {
    expect(stats.kurt, -1.5);
  });

  test('max', async () => {
    expect(stats.max, 5);
  });

  test('med', async () => {
    expect(stats.med, 4);
  });

  test('min', async () => {
    expect(stats.min, 3);
  });

  test('missingValueCount', async () => {
    expect(stats.missingValueCount, 0);
  });

  test('q1', async () => {
    expect(stats.q1, 3);
  });

  test('q2', async () => {
    expect(stats.q2, 4);
  });

  test('q3', async () => {
    expect(stats.q3, 5);
  });

  test('skew', async () => {
    expect(stats.skew, 0);
  });

  test('stdev', async () => {
    expect(stats.stdev, 1);
  });

  test('sum', async () => {
    expect(stats.sum, 12);
  });

  test('totalCount', async () => {
    expect(stats.totalCount, 3);
  });

  test('uniqueCount', async () => {
    expect(stats.uniqueCount, 3);
  });

  test('valueCount', async () => {
    expect(stats.valueCount, 3);
  });

  test('variance', async () => {
    expect(stats.variance, 1);
  });

  test('corr', async () => {
    expect(stats.corr(t.getCol('number')), 0.5960395606792696);
  });

  test('spearmanCorr', async () => {
    expect(stats.spearmanCorr(t.getCol('number')), 0.4999999999999999);
  });

  test('toString()', async () => {
    const t = DG.DataFrame.create(1);
    t.columns.add(DG.Column.fromInt32Array('number', Int32Array.from([12])));
    t.columns.add(DG.Column.fromInt32Array('population', Int32Array.from([3])));
    const stats = DG.Stats.histogramsByCategories(t.getCol('number'), t.getCol('population'));
    expect(stats.toString(), '1,0,0,0,0,0,0,0,0,0');
  });
}, {owner: 'dkovalyov@datagrok.ai'});

category('Stats: fromValues', () => {
  const floatArr = new Float32Array([1.5, 2.3, 3.7, 4.8, 5.23]);
  const intArr = new Int8Array([1, 2, 3, 4, 5]);
  const jsArr = [1, 2, 3, 4, 5];
  
  const floatStats = DG.Stats.fromValues(floatArr);
  const intStats = DG.Stats.fromValues(intArr);
  const jsStats = DG.Stats.fromValues(jsArr);

  test('avg float', async () => {
    expect(floatStats.avg, 3.50600004196167);
  });
  
  test('avg int', async () => {
    expect(intStats.avg, 3);
  });
  
  test('avg js', async () => {
    expect(jsStats.avg, 3);
  });
}, {owner: 'dkovalyov@datagrok.ai'});

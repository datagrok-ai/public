import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {category, expect, test} from '@datagrok-libraries/utils/src/test';

category('Stats', () => {
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
});

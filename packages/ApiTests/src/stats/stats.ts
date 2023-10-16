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
});

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

  test('kurt float', async () => {
    expect(floatStats.kurt, -1.5502044653377132);
  });

  test('kurt int', async () => {
    expect(intStats.kurt, -1.3);
  });

  test('kurt js', async () => {
    expect(jsStats.kurt, -1.3);
  });

  test('max float', async () => {
    expect(floatStats.max, 5.230000019073486);
  });

  test('max int', async () => {
    expect(intStats.max, 5);
  });

  test('max js', async () => {
    expect(jsStats.max, 5);
  });

  test('med float', async () => {
    expect(floatStats.med, 3.700000047683716);
  });

  test('med int', async () => {
    expect(intStats.med, 3);
  });

  test('med js', async () => {
    expect(jsStats.med, 3);
  });

  test('min float', async () => {
    expect(floatStats.min, 1.5);
  });

  test('min int', async () => {
    expect(intStats.min, 1);
  });

  test('min js', async () => {
    expect(jsStats.min, 1);
  });

  test('missingValueCount float', async () => {
    expect(floatStats.missingValueCount, 0);
  });

  test('missingValueCount int', async () => {
    expect(intStats.missingValueCount, 0);
  });

  test('missingValueCount js', async () => {
    expect(jsStats.missingValueCount, 0);
  });

  test('q1 float', async () => {
    expect(floatStats.q1, 2.299999952316284);
  });

  test('q1 int', async () => {
    expect(intStats.q1, 2);
  });

  test('q1 js', async () => {
    expect(jsStats.q1, 2);
  });

  test('q2 float', async () => {
    expect(floatStats.q2, 3.700000047683716);
  });

  test('q2 int', async () => {
    expect(intStats.q2, 3);
  });

  test('q2 js', async () => {
    expect(jsStats.q2, 3);
  });

  test('q3 float', async () => {
    expect(floatStats.q3, 4.800000190734863);
  });

  test('q3 int', async () => {
    expect(intStats.q3, 4);
  });

  test('q3 js', async () => {
    expect(jsStats.q3, 4);
  });

  test('skew float', async () => {
    expect(floatStats.skew, -0.17449530458852078);
  });

  test('skew int', async () => {
    expect(intStats.skew, 0);
  });

  test('skew js', async () => {
    expect(jsStats.skew, 0);
  });

  test('stdev float', async () => {
    expect(floatStats.stdev, 1.5939511200866263);
  });

  test('stdev int', async () => {
    expect(intStats.stdev, 1.5811388300841898);
  });

  test('stdev js', async () => {
    expect(jsStats.stdev, 1.5811388300841898);
  });

  test('sum float', async () => {
    expect(floatStats.sum, 17.53000020980835);
  });

  test('sum int', async () => {
    expect(intStats.sum, 15);
  });

  test('sum js', async () => {
    expect(jsStats.sum, 15);
  });

  test('totalCount float', async () => {
    expect(floatStats.totalCount, 5);
  });

  test('totalCount int', async () => {
    expect(intStats.totalCount, 5);
  });

  test('totalCount js', async () => {
    expect(jsStats.totalCount, 5);
  });

  test('uniqueCount float', async () => {
    expect(floatStats.uniqueCount, 5);
  });

  test('uniqueCount int', async () => {
    expect(intStats.uniqueCount, 5);
  });

  test('uniqueCount js', async () => {
    expect(jsStats.uniqueCount, 5);
  });

  test('valueCount float', async () => {
    expect(floatStats.valueCount, 5);
  });

  test('valueCount int', async () => {
    expect(intStats.valueCount, 5);
  });

  test('valueCount js', async () => {
    expect(jsStats.valueCount, 5);
  });

  test('variance float', async () => {
    expect(floatStats.variance, 2.540680173225411);
  });

  test('variance int', async () => {
    expect(intStats.variance, 2.5);
  });

  test('variance js', async () => {
    expect(jsStats.variance, 2.5);
  });
});

/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';

import {category, test, expect, expectFloat, awaitCheck} from '@datagrok-libraries/test/src/test';
import {FitConstants} from '@datagrok-libraries/statistics/src/fit/const';

const CONCENTRATIONS = [1e-9, 3e-9, 1e-8, 3e-8, 1e-7, 3e-7, 1e-6, 3e-6, 1e-5, 3e-5, 1e-4];

/** Sigmoid dose-response curve JSON with its inflection point at 10^-`logIC50`. */
function curveJson(logIC50: number): string {
  return JSON.stringify({
    chartOptions: {logX: true},
    series: [{fitFunction: 'sigmoid', name: 'series', points: CONCENTRATIONS.map((x) => ({
      x: x, y: 5 + 95 / (1 + Math.pow(10, Math.log10(x) - logIC50)),
    }))}],
  });
}

function curveTable(name: string, logIC50s: number[]): DG.DataFrame {
  const col = DG.Column.fromStrings('curve', logIC50s.map(curveJson));
  col.semType = FitConstants.FIT_SEM_TYPE;
  const df = DG.DataFrame.fromColumns([col]);
  df.name = name;
  return df;
}

/** Runs a curve statistic function the way the property panel's "+" button does. */
async function addStatisticColumn(df: DG.DataFrame, funcName: string,
  params: {[key: string]: string | number}): Promise<void> {
  await DG.Func.find({name: funcName})[0]
    .prepare({table: df, curveColumn: df.col('curve')!, ...params})
    .call(false, undefined, {processed: false});
}

category('calculated columns', () => {
  test('curveStatistic adds a calculated column', async () => {
    const df = curveTable('calcColAdd', [-6.5, -6.5]);
    await addStatisticColumn(df, 'curveStatistic', {propName: 'ic50', seriesNumber: 0});

    const col = df.col('curve 1 ic50');
    expect(col !== null);
    expect(col!.type, DG.COLUMN_TYPE.FLOAT);
    // the statistic must come back in data space, not as the log-space fitted parameter
    expectFloat(Math.log10(col!.get(0)), -6.5, 0.15);

    // the platform owns the formula, and binds the curve column as a column rather than a per-row value
    expect(col!.getTag(DG.TAGS.FORMULA).includes('Curves:curveStatistic'));
    expect(col!.getTag(DG.TAGS.FORMULA).includes('${curve}'));
    expect(col!.getTag('.%formula-column-is-vector-func'), 'true');
  });

  test('recalculates when the source curve changes', async () => {
    const df = curveTable('calcColRecalc', [-6.5, -6.5]);
    await addStatisticColumn(df, 'curveStatistic', {propName: 'ic50', seriesNumber: 0});
    const before = df.col('curve 1 ic50')!.get(1);

    // regression guard: an unstable result column name, or reading the value through
    // column.dataFrame, both leave this silently stale instead of failing
    df.col('curve')!.set(1, curveJson(-5.0), true);
    await awaitCheck(() => df.col('curve 1 ic50')!.get(1) !== before,
      'statistic column did not recalculate after the source curve changed', 5000);

    expectFloat(Math.log10(df.col('curve 1 ic50')!.get(1)), -5.0, 0.3);
    expect(df.col('curve 1 ic50')!.get(0), before, 'untouched row must not be recalculated');
  });

  test('legacy statistic names still resolve', async () => {
    const df = curveTable('calcColLegacy', [-6.5]);
    await addStatisticColumn(df, 'curveStatistic', {propName: 'interceptX', seriesNumber: 0});

    // interceptX is the pre-refactor name for the sigmoid inflection point
    const col = df.col('curve 1 interceptX');
    expect(col !== null);
    expectFloat(Math.log10(col!.get(0)), -6.5, 0.15);
  });

  test('curveAggrStatistic aggregates across series', async () => {
    const df = curveTable('calcColAggr', [-6.5, -5.5]);
    await addStatisticColumn(df, 'curveAggrStatistic', {propName: 'ic50', aggrType: 'med'});

    const col = df.col('curve med ic50');
    expect(col !== null);
    // single-series cells, so the median across series is that series' own IC50
    expectFloat(Math.log10(col!.get(0)), -6.5, 0.15);
    expectFloat(Math.log10(col!.get(1)), -5.5, 0.15);
  });
});

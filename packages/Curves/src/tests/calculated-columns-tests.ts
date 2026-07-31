/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';

import {category, test, expect, expectFloat, awaitCheck} from '@datagrok-libraries/test/src/test';
import {FitConstants} from '@datagrok-libraries/statistics/src/fit/const';
import {IFitChartData} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {getStatistic} from '@datagrok-libraries/statistics/src/fit/fit-engine';
import {getChartDataAggrStats, aggregatedStatisticsProperties, calculateSeriesFit} from '../fit/fit-statistics';
import {setOutlier} from '../fit/fit-renderer';
import {getOrCreateParsedChartData} from '../fit/fit-chart-data';

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

/** One cell holding several sigmoid series, each with its own inflection point. */
function multiSeriesCurveJson(logIC50s: number[]): string {
  return JSON.stringify({
    chartOptions: {logX: true},
    series: logIC50s.map((logIC50, i) => ({
      fitFunction: 'sigmoid', name: `series ${i}`, points: CONCENTRATIONS.map((x) => ({
        x: x, y: 5 + 95 / (1 + Math.pow(10, Math.log10(x) - logIC50)),
      })),
    })),
  });
}

/** Same curve, but with logX left to the column level rather than baked into the cell. */
function curveJsonNoChartOptions(logIC50: number): string {
  const parsed = JSON.parse(curveJson(logIC50));
  delete parsed.chartOptions;
  return JSON.stringify(parsed);
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

  test('averaging IC50 across series is a geometric mean', async () => {
    const data = JSON.parse(multiSeriesCurveJson([-7, -5])) as IFitChartData;

    // aggregated in log space, so 1e-7 and 1e-5 give 1e-6; per-series conversion would give 5.05e-6
    expectFloat(Math.log10(getChartDataAggrStats(data, 'avg').ic50!), -6, 0.1);
  });

  test('mixed fit functions only aggregate statistics common to all series', async () => {
    const data = JSON.parse(multiSeriesCurveJson([-7, -5])) as IFitChartData;
    data.series![1].fitFunction = 'linear';

    // ic50 belongs to the sigmoid only; rSquared is produced by both
    const names = aggregatedStatisticsProperties(data).map((p) => p.name);
    expect(names.includes('ic50'), false, 'ic50 is not viable across a sigmoid and a linear fit');
    expect(names.includes('rSquared'), true);
    expect(getChartDataAggrStats(data, 'avg').ic50 === undefined, true, 'ic50 must not be aggregated over a subset');
  });

  test('repeated calculation does not re-log stored parameters', async () => {
    // an IC50 above 1 (nM/uM-scale data) stays positive after one log, so an in-place conversion
    // logs it again on the next call: 100 -> 2 -> 0.301
    const series = JSON.parse(curveJson(-6.5)).series[0];
    series.parameters = [100, 1, 100, 5];
    const logX = {logX: true, logY: false};

    const first = getStatistic(calculateSeriesFit(series, 0, logX, undefined, false), 'ic50');
    const second = getStatistic(calculateSeriesFit(series, 0, logX, undefined, false), 'ic50');
    expect(first, second, 'stored parameters drifted between calls');
    expectFloat(first!, 100, 1);
  });

  test('fitting a series does not write fit-space parameters back onto it', async () => {
    // no stored parameters, so the first call fits. If the fitted (log-space) parameters are written
    // back, the second call reads them as a concentration and logs them again
    const data = JSON.parse(multiSeriesCurveJson([2])) as IFitChartData;
    const series = data.series![0];
    const logX = {logX: true, logY: false};

    const first = getStatistic(calculateSeriesFit(series, 0, logX, undefined, false), 'ic50');
    expect(series.parameters === undefined, true, 'the series must keep its data-space contract');
    const second = getStatistic(calculateSeriesFit(series, 0, logX, undefined, false), 'ic50');
    expect(first, second, 'the reported IC50 drifted between calls');
  });

  test('recalculation applies column-level chart options', async () => {
    const col = DG.Column.fromStrings('curve', [curveJsonNoChartOptions(-6.5), curveJsonNoChartOptions(-6.5)]);
    col.semType = FitConstants.FIT_SEM_TYPE;
    col.setTag(FitConstants.TAG_FIT, JSON.stringify({chartOptions: {logX: true}}));
    const df = DG.DataFrame.fromColumns([col]);
    df.name = 'calcColColumnOptions';
    await addStatisticColumn(df, 'curveStatistic', {propName: 'ic50', seriesNumber: 0});
    const before = df.col('curve 1 ic50')!.get(1);

    df.col('curve')!.set(1, curveJsonNoChartOptions(-5.0), true);
    await awaitCheck(() => df.col('curve 1 ic50')!.get(1) !== before,
      'statistic column did not recalculate', 5000);

    // the detached recalculation column has no dataframe, so it must read logX off its own tags
    expectFloat(Math.log10(df.col('curve 1 ic50')!.get(1)), -5.0, 0.3);
  });

  test('log-space statistics are not offered for aggregations we do not convert', async () => {
    const data = JSON.parse(multiSeriesCurveJson([-7, -5])) as IFitChartData;

    // ic50 lives in log space under logX, and pIC50 is only derived during the conversion. Under an
    // aggregation we do not convert they would be a log-space number and a null, so both are dropped -
    // this also guards the regression where unlogging every aggregation turned a count of 2 into 100
    const counted = aggregatedStatisticsProperties(data, 'count').map((p) => p.name);
    expect(counted.includes('ic50'), false, 'ic50 would be a log-space number under count');
    expect(counted.includes('pIC50'), false, 'pIC50 is only derived in data space');
    expect(counted.includes('rSquared'), true, 'space-independent statistics stay');

    // they come back for an aggregation that is converted
    const averaged = aggregatedStatisticsProperties(data, 'avg').map((p) => p.name);
    expect(averaged.includes('ic50'), true);
    expect(getChartDataAggrStats(data, 'count').ic50 === undefined, true);
  });

  test('outlier toggle refreshes statistic columns named by the new API', async () => {
    const df = curveTable('calcColOutlierToggle', [-6.5]);
    const stat = DG.Column.float('ic50 col', df.rowCount);
    // tags the legacy "+" writes; ic50 has no slot in the seven legacy keys
    stat.tags['.sourceColumn'] = 'curve';
    stat.tags['.seriesNumber'] = '0';
    stat.tags['.statistics'] = 'ic50';
    stat.init(() => 1e-6);
    df.columns.add(stat);

    const gridCell = DG.Viewer.grid(df).cell('curve', 0);
    const data = getOrCreateParsedChartData(gridCell.cell, true);
    setOutlier(gridCell, data.series![0].points[0], 0, 0, data);

    const value = df.col('ic50 col')!.get(0);
    expect(value != null && !isNaN(value), true, 'outlier toggle blanked the ic50 column');
    expectFloat(Math.log10(value), -6.5, 0.3);
  });
});

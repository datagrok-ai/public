/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';

import {category, test, expect, expectArray, expectFloat, awaitCheck} from '@datagrok-libraries/test/src/test';
import {FitConstants} from '@datagrok-libraries/statistics/src/fit/const';
import {IFitChartData} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {getStatistic} from '@datagrok-libraries/statistics/src/fit/fit-engine';
import {getChartDataAggrStats, aggregatedStatisticsProperties, calculateSeriesFit, curveStatisticAt} from '../fit/fit-statistics';
import {setOutlier} from '../fit/fit-renderer';
import {getOrCreateParsedChartData} from '../fit/fit-chart-data';
import {CONCENTRATIONS, curveJson, multiSeriesCurveJson} from './curve-data';
import {addStatisticColumn as addStatisticColumnFromPanel} from '../fit/fit-grid-cell-handler';

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

/** Grid column names in grid order; index 0 is the row header, which has no name. */
function gridColumnNames(grid: DG.Grid): string[] {
  const names: string[] = [];
  for (let i = 0; i < grid.columns.length; i++) {
    const name = grid.columns.byIndex(i)?.name;
    if (name)
      names.push(name);
  }
  return names;
}

category('calculated columns', () => {
  test('the statistic column is added right after the curve column', async () => {
    // `join(table)` appends, so the position the legacy addStatisticsColumn inserted at is restored
    const df = curveTable('calcColPosition', [-6.5]);
    df.columns.addNewString('trailing').set(0, 'x');
    const gridCell = DG.Viewer.grid(df).cell('curve', 0);

    await addStatisticColumnFromPanel(gridCell, 'curveStatistic', {propName: 'ic50', seriesNumber: 0});

    expectArray(df.columns.names(), ['curve', 'curve 1 ic50', 'trailing']);
    // the grid only takes a position from the dataframe when it first builds a GridColumn
    expectArray(gridColumnNames(gridCell.grid), ['curve', 'curve 1 ic50', 'trailing']);
  });

  test('placing the statistic column keeps a manually reordered grid', async () => {
    // a dragged grid is no longer in dataframe order, and must not be snapped back to it
    const df = curveTable('calcColGridOrder', [-6.5]);
    df.columns.addNewString('a').set(0, 'x');
    df.columns.addNewString('b').set(0, 'y');
    const grid = DG.Viewer.grid(df);
    grid.columns.setOrder(['b', 'curve', 'a']);

    await addStatisticColumnFromPanel(grid.cell('curve', 0), 'curveStatistic', {propName: 'ic50', seriesNumber: 0});

    expectArray(gridColumnNames(grid), ['b', 'curve', 'curve 1 ic50', 'a']);
    expectArray(df.columns.names(), ['curve', 'curve 1 ic50', 'a', 'b']);
  });

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

    // an unstable result name, or reading through column.dataFrame, leaves this silently stale
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
    // an IC50 above 1 stays positive after one log, so an in-place conversion drifts: 100 -> 2 -> 0.301
    const series = JSON.parse(curveJson(-6.5)).series[0];
    series.parameters = [100, 1, 100, 5];
    const logX = {logX: true, logY: false};

    const first = getStatistic(calculateSeriesFit(series, 0, logX, undefined, false), 'ic50');
    const second = getStatistic(calculateSeriesFit(series, 0, logX, undefined, false), 'ic50');
    expect(first, second, 'stored parameters drifted between calls');
    expectFloat(first!, 100, 1);
  });

  test('fitting a series does not write fit-space parameters back onto it', async () => {
    // the first call fits; writing the log-space parameters back would make the second log them again
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

    // pIC50 is only derived during the conversion, so it is the one that cannot be produced
    const counted = aggregatedStatisticsProperties(data, 'count').map((p) => p.name);
    expect(counted.includes('pIC50'), false, 'pIC50 is only derived in data space');
    expect(counted.includes('rSquared'), true, 'space-independent statistics stay');

    // ic50 stays available and is reported in fit space, which is what a stdev of it should be, and
    // what a saved project recorded before this refactor - it must not come back as a null column
    expect(counted.includes('ic50'), true, 'ic50 must stay available for aggregations we do not convert');
    expectFloat(getChartDataAggrStats(data, 'avg').ic50!, 1e-6, 1e-7);
    expectFloat(getChartDataAggrStats(data, 'count').ic50!, 2, 0.001);
    expectFloat(getChartDataAggrStats(data, 'count').interceptX!, 2, 0.001);
  });

  test('stored parameters round-trip through fit space on both axes', async () => {
    // seriesInFitSpace and toDataSpace must be inverses, or a top of 100 comes back as 1e100
    for (const logOptions of [{logX: false, logY: false}, {logX: true, logY: false},
      {logX: false, logY: true}, {logX: true, logY: true}]) {
      const series = JSON.parse(curveJson(-6.5)).series[0];
      series.parameters = [100, 1, 1e-6, 5];
      const fit = calculateSeriesFit(series, 0, logOptions, undefined, false);

      expectFloat(getStatistic(fit, 'top')!, 100, 1e-6, `top under ${JSON.stringify(logOptions)}`);
      expectFloat(getStatistic(fit, 'bottom')!, 5, 1e-6, `bottom under ${JSON.stringify(logOptions)}`);
      expectFloat(getStatistic(fit, 'ic50')!, 1e-6, 1e-12, `ic50 under ${JSON.stringify(logOptions)}`);
    }
  });

  test('a detached column strips the stray pipe like the initial parse does', async () => {
    // recalculation hands the function a detached column, whose parse must sanitize the value too
    const col = DG.Column.fromStrings('curve', [curveJson(-6.5).replace('{', '{|')]);
    col.semType = FitConstants.FIT_SEM_TYPE;
    expect(col.dataFrame === null || col.dataFrame === undefined, true, 'the column must be detached');

    const value = curveStatisticAt(col, 0, 'ic50', 0);
    expect(value !== null, true, 'a pipe in the cell value blanked the row');
    expectFloat(Math.log10(value!), -6.5, 0.3);
  });

  test('a legacy statistic reads the same per series and aggregated', async () => {
    // `top` on a linear fit resolves through the positional fallback, per series and aggregated alike
    const data = JSON.parse(multiSeriesCurveJson([-6.5])) as IFitChartData;
    data.series![0].fitFunction = 'linear';
    const logX = {logX: true, logY: false};

    const perSeries = getStatistic(calculateSeriesFit(data.series![0], 0, logX, undefined, false), 'top');
    const aggregated = getChartDataAggrStats(data, 'med').top;
    expect(perSeries !== undefined, true, 'the per-series path lost the fallback');
    expect(aggregated !== undefined, true, 'the aggregated path drops a statistic the per-series one reports');
  });

  test('a legacy name and its canonical field agree after aggregation', async () => {
    // interceptX aliases onto ic50, so collecting it before the conversion leaves it in fit space
    const data = JSON.parse(multiSeriesCurveJson([-7, -5])) as IFitChartData;
    const stats = getChartDataAggrStats(data, 'avg');

    expect(stats.interceptX, stats.ic50, 'interceptX must be reported in the same space as ic50');
    expectFloat(Math.log10(stats.interceptX!), -6, 0.1);
  });

  test('mergeSeries does not change a statistic', async () => {
    // it merges on a copy for the plot only, so it cannot change an extracted statistic
    const plain = JSON.parse(multiSeriesCurveJson([-7, -5])) as IFitChartData;
    const merged = JSON.parse(multiSeriesCurveJson([-7, -5])) as IFitChartData;
    merged.chartOptions!.mergeSeries = true;
    const logX = {logX: true, logY: false};

    for (const i of [0, 1]) {
      const a = getStatistic(calculateSeriesFit(plain.series![i], i, logX, undefined, false), 'ic50');
      const b = getStatistic(calculateSeriesFit(merged.series![i], i, logX, undefined, false), 'ic50');
      expect(a, b, `series ${i} changed when mergeSeries was set`);
    }
    expect(getChartDataAggrStats(merged, 'med').ic50, getChartDataAggrStats(plain, 'med').ic50);
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

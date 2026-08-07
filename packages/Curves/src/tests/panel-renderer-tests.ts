/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {category, test, expect, expectArray, expectFloat, awaitCheck, delay} from '@datagrok-libraries/test/src/test';
import {FitConstants} from '@datagrok-libraries/statistics/src/fit/const';
import {IFitChartData} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {FitChartCellRenderer} from '../fit/fit-renderer';
import {FitGridCellHandler} from '../fit/fit-grid-cell-handler';
import {normalizeStatisticNames, chartPropertiesFor, changeCurvesOptions} from '../fit/fit-options';
import {sigmoidPoints} from './curve-data';
import {getOrCreateParsedChartData, getColumnChartOptions} from '../fit/fit-chart-data';

/** Curve carrying stored parameters whose inflection point is above 1, where converting in place
 * would log it again on the next pass. */
function curveWithParameters(): IFitChartData {
  return {
    chartOptions: {logX: true},
    series: [{fitFunction: 'sigmoid', name: 'series', parameters: [100, 1, 100, 5], points: sigmoidPoints(-6.5)}],
  } as IFitChartData;
}

function curveTable(name: string, chartOptions: boolean, showStatistics?: string[]): DG.DataFrame {
  const cell = JSON.stringify({
    ...(chartOptions ? {chartOptions: {logX: true, ...(showStatistics ? {showStatistics} : {})}} : {}),
    series: [{fitFunction: 'sigmoid', name: 'series', points: sigmoidPoints(-6.5)}],
  });
  const col = DG.Column.fromStrings('curve', [cell, cell]);
  col.semType = FitConstants.FIT_SEM_TYPE;
  const df = DG.DataFrame.fromColumns([col]);
  df.name = name;
  return df;
}

async function addStatisticColumn(df: DG.DataFrame, params: {[key: string]: string | number}): Promise<void> {
  await DG.Func.find({package: 'Curves', name: 'curveStatistic'})[0]
    .prepare({table: df, curveColumn: df.col('curve')!, ...params})
    .call(false, undefined, {processed: false});
}

category('panel and renderer', () => {
  test('every registered aggregation input offers only convertible aggregations', async () => {
    // enumerated from the registry: a hand-written list of the ones we fixed missed a third
    const CONVERTIBLE = ['min', 'max', 'avg', 'med', 'q1', 'q2', 'q3'];
    const checked: string[] = [];
    for (const func of DG.Func.find({package: 'Curves'})) {
      for (const input of func.inputs.filter((p: DG.Property) => /^(aggrType|aggregation)$/.test(p.name))) {
        checked.push(func.name);
        expect(input.choices !== undefined && input.choices !== null, true,
          `${func.name} takes a free-text aggregation`);
        for (const choice of input.choices!)
          expect(CONVERTIBLE.includes(choice), true, `${func.name} offers ${choice}`);
      }
    }
    expect(checked.length >= 3, true, `expected every aggregation function to be checked, saw ${checked}`);
  });

  test('rendering does not mutate the cached series parameters', async () => {
    const chartData = curveWithParameters();
    const before = [...chartData.series![0].parameters!];

    const canvas = ui.canvas(200, 120);
    const g = canvas.getContext('2d')!;
    const renderer = new FitChartCellRenderer();
    const bounds = FitChartCellRenderer.inflateScreenBounds(new DG.Rect(0, 0, 200, 120));
    // the cached chart data is shared with the statistics, so a repaint must not convert in place
    renderer.renderCurves(g, bounds, chartData);
    renderer.renderCurves(g, bounds, chartData);

    expectArray(chartData.series![0].parameters!, before);
  });

  test('firing values-changed on a curve column recalculates its statistic columns', async () => {
    // logX lives only on the column, so flipping it changes the space the curve is fitted in
    const df = curveTable('panelColumnOption', false);
    const col = df.col('curve')!;
    col.setTag(FitConstants.TAG_FIT, JSON.stringify({chartOptions: {logX: true}}));
    await addStatisticColumn(df, {propName: 'ic50', seriesNumber: 0});
    const before = df.col('curve 1 ic50')!.get(0);
    expectFloat(Math.log10(before), -6.5, 0.3);

    col.setTag(FitConstants.TAG_FIT, JSON.stringify({chartOptions: {logX: false}}));
    // pins the platform contract the panel relies on, not the call site in changeCurvesOptions
    col.fireValuesChanged();
    await awaitCheck(() => df.col('curve 1 ic50')!.get(0) !== before,
      'statistic column did not recalculate after the column-level option changed', 5000);
  });

  test('clearing a cell-level override does not recalculate on its own', async () => {
    // the rewrite that drops cell claims must not read as a data change; a recalculation bumps version
    const cell = JSON.stringify({
      chartOptions: {logX: true},
      series: [{fitFunction: 'sigmoid', name: 'series', pointColor: '#ff0000', points: sigmoidPoints(-6.5)}],
    });
    const col = DG.Column.fromStrings('curve', [cell, cell]);
    col.semType = FitConstants.FIT_SEM_TYPE;
    const df = DG.DataFrame.fromColumns([col]);
    df.name = 'panelNoRecalcOnCosmetic';
    await addStatisticColumn(df, {propName: 'ic50', seriesNumber: 0});

    const stat = df.col('curve 1 ic50')!;
    const before = stat.version;
    for (let j = 0; j < col.length; j++) {
      const data = JSON.parse(col.get(j));
      delete data.series[0].pointColor;
      col.set(j, JSON.stringify(data), false);
    }
    await delay(500);
    expect(stat.version, before, 'a cosmetic rewrite recalculated the statistic column');
  });

  test('a column-level option overrides a cell-level one, without recalculating', async () => {
    // both halves in one test: notifying recalculates, not notifying used to leave the cache stale
    const cell = JSON.stringify({
      chartOptions: {logX: true},
      series: [{fitFunction: 'sigmoid', name: 'series', pointColor: '#ff0000', points: sigmoidPoints(-6.5)}],
    });
    const col = DG.Column.fromStrings('curve', [cell, cell]);
    col.semType = FitConstants.FIT_SEM_TYPE;
    const df = DG.DataFrame.fromColumns([col]);
    df.name = 'panelColumnOverridesCell';
    await addStatisticColumn(df, {propName: 'ic50', seriesNumber: 0});

    const stat = df.col('curve 1 ic50')!;
    const version = stat.version;
    const gridCell = DG.Viewer.grid(df).cell('curve', 0);
    const input = {property: {name: 'pointColor'}, value: '#00ff00'} as unknown as DG.InputBase;

    changeCurvesOptions(gridCell, input, 'seriesOptions', 'Column');
    await delay(500);

    expect(getOrCreateParsedChartData(df.cell(0, 'curve')).series![0].pointColor === '#ff0000', false,
      'the cell-level value still wins after the column-level change');
    expect(stat.version, version, 'a cosmetic option change recalculated the statistic column');
  });

  test('clearing cell overrides marks the table changed without recalculating', async () => {
    // silent so a cosmetic option does not refit, but the table still has to read as modified
    const cell = JSON.stringify({
      chartOptions: {logX: true},
      series: [{fitFunction: 'sigmoid', name: 'series', showCurveConfidenceInterval: true, points: sigmoidPoints(-6.5)}],
    });
    const col = DG.Column.fromStrings('curve', [cell, cell]);
    col.semType = FitConstants.FIT_SEM_TYPE;
    const df = DG.DataFrame.fromColumns([col]);
    df.name = 'panelMarksTableChanged';
    await addStatisticColumn(df, {propName: 'ic50', seriesNumber: 0});

    const stat = df.col('curve 1 ic50')!;
    const version = stat.version;
    const gridCell = DG.Viewer.grid(df).cell('curve', 0);
    let fired = 0;
    const sub = df.onValuesChanged.subscribe(() => fired++);

    changeCurvesOptions(gridCell,
      {property: {name: 'showCurveConfidenceInterval'}, value: false} as unknown as DG.InputBase,
      'seriesOptions', 'Column');
    await delay(500);
    sub.unsubscribe();

    // the value the data declares is left alone - it just stops winning
    expect(JSON.parse(col.get(0)).series[0].showCurveConfidenceInterval, true, 'the cell data was mutated');
    expect(getOrCreateParsedChartData(df.cell(0, 'curve')).series![0].showCurveConfidenceInterval, false,
      'the column-level value did not outrank the value declared in the cell');
    expect(fired > 0, true, 'nothing marked the table changed, so the option would not be saved');
    expect(stat.version, version, 'a cosmetic option change recalculated the statistic column');
  });

  test('a column-level change overrides a value set at cell level first', async () => {
    // detectSettings recorded "no cell claims this" once, so a later cell write was invisible to it
    const df = curveTable('panelCellThenColumn', true);
    const gridCell = DG.Viewer.grid(df).cell('curve', 0);
    const prop = (value: any) => ({property: {name: 'connectDots'}, value}) as unknown as DG.InputBase;

    changeCurvesOptions(gridCell, prop(true), 'seriesOptions', 'Cell');
    expect(JSON.parse(df.col('curve')!.get(0)).series[0].connectDots, true, 'cell-level write failed');
    expect(JSON.parse(df.col('curve')!.get(0)).explicit.seriesOptions.includes('connectDots'), true,
      'a cell-level change did not claim the property');

    changeCurvesOptions(gridCell, prop(false), 'seriesOptions', 'Column');
    expect(JSON.parse(df.col('curve')!.get(0)).explicit.seriesOptions.includes('connectDots'), false,
      'the cell-level claim survived a column-level change');
    expect(getOrCreateParsedChartData(df.cell(0, 'curve')).series![0].connectDots, false,
      'the column-level value did not reach the cell');
  });

  test('a column-level option outranks a value the fresh data declares', async () => {
    // a layout applied to fresh data: the column tag comes back, the cells still declare the property
    const cell = JSON.stringify({
      chartOptions: {logX: true},
      series: [{fitFunction: 'sigmoid', name: 'series', showCurveConfidenceInterval: true, points: sigmoidPoints(-6.5)}],
    });
    const col = DG.Column.fromStrings('curve', [cell, cell]);
    col.semType = FitConstants.FIT_SEM_TYPE;
    const df = DG.DataFrame.fromColumns([col]);
    df.name = 'panelFreshDataClaim';
    // the column tag as a layout would restore it: the value plus the record that it was set here
    col.setTag(FitConstants.TAG_FIT, JSON.stringify({
      seriesOptions: {showCurveConfidenceInterval: false},
      explicit: {seriesOptions: ['showCurveConfidenceInterval']},
    }));

    expect(getOrCreateParsedChartData(df.cell(0, 'curve')).series![0].showCurveConfidenceInterval, false,
      'the restored column-level option lost to the value the cell declares');
    // and without the record it stays advisory, so a series keeps whatever the data gave it
    col.setTag(FitConstants.TAG_FIT, JSON.stringify({seriesOptions: {showCurveConfidenceInterval: false}}));
    expect(getOrCreateParsedChartData(df.cell(1, 'curve')).series![0].showCurveConfidenceInterval, true,
      'an option nobody set explicitly overrode the data');
  });

  test('options stored under the legacy tag are read and migrate onto the layout-carried one', async () => {
    // only '.%' tags reach a layout, so options in '.fit' never travelled to a fresh table
    expect(FitConstants.TAG_FIT.startsWith('.%'), true, 'the fit tag is no longer carried by layouts');

    const df = curveTable('panelLegacyTag', false);
    const col = df.col('curve')!;
    col.setTag(FitConstants.TAG_FIT_LEGACY, JSON.stringify({chartOptions: {logX: true}}));

    expect(getColumnChartOptions(col).chartOptions!.logX, true, 'options under the legacy tag were ignored');
    expect(JSON.parse(col.getTag(FitConstants.TAG_FIT)).chartOptions.logX, true,
      'reading did not migrate the options onto the layout-carried tag');
  });

  test('property panel renders for a saved legacy statistic', async () => {
    // no test covered this panel before, and it is the surface this ticket changed most
    const df = curveTable('panelBinding', true, ['interceptX']);
    const gridCell = DG.Viewer.grid(df).cell('curve', 0);

    // the panes build their contents lazily, so this asserts the panel is constructed without
    // throwing - which nothing covered before, on the surface this ticket changed most
    const root = new FitGridCellHandler().renderProperties(gridCell);
    expect(root instanceof HTMLElement, true, 'the property panel failed to render');
  });

  test('statistic names stored under a legacy name map onto the current ones', async () => {
    const chartData = JSON.parse(JSON.stringify({
      chartOptions: {logX: true},
      series: [{fitFunction: 'sigmoid', name: 's', points: sigmoidPoints(-6.5)}],
    })) as IFitChartData;

    // a saved project ticks its checkbox only if the stored name resolves to the current one
    expectArray(normalizeStatisticNames(chartData, ['interceptX']), ['ic50']);
    expectArray(normalizeStatisticNames(chartData, ['auc', 'rSquared']), ['auc', 'rSquared']);
    // the same statistic under both names collapses to one entry
    expectArray(normalizeStatisticNames(chartData, ['interceptX', 'ic50']), ['ic50']);

    // and the choices the panel offers come from the cell's own fit function
    const names = chartPropertiesFor(chartData).find((p) => p.name === 'showStatistics')!.choices!;
    expect(names.includes('ic50'), true);
    expect(names.includes('ec50'), false, 'a sigmoid cell must not offer the 4PL-regression name');
  });
});

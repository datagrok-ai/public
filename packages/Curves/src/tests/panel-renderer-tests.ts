/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {category, test, expect, expectArray, expectFloat, awaitCheck, delay} from '@datagrok-libraries/test/src/test';
import {FitConstants} from '@datagrok-libraries/statistics/src/fit/const';
import {IFitChartData} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {FitChartCellRenderer} from '../fit/fit-renderer';
import {stdevWhisker} from '../fit/render-utils';
import {FitGridCellHandler} from '../fit/fit-grid-cell-handler';
import {normalizeStatisticNames, chartPropertiesFor, changeCurvesOptions, seriesPropertiesFor} from '../fit/fit-options';
import {sigmoidPoints, ciCurveJson, labelledCurveJson, curveJson, renderedTexts} from './curve-data';
import {multiSeriesCurveJson} from './curve-data';
import {getOrCreateParsedChartData, getColumnChartOptions} from '../fit/fit-chart-data';
import {getSeriesFitFunction} from '@datagrok-libraries/statistics/src/fit/fit-data';
import {fitFunctions, fitFunctionDescriptions} from '@datagrok-libraries/statistics/src/fit/fit-engine';

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
    const col = DG.Column.fromStrings('curve', [ciCurveJson(), ciCurveJson()]);
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
    const col = DG.Column.fromStrings('curve', [ciCurveJson(), ciCurveJson()]);
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

  test('a column-level option survives a layout applied to a fresh table', async () => {
    // the round trip this work exists for: the layout carries the column tag, the fresh table brings
    // the source's own cells back, and the option has to win anyway
    const cells = () => [ciCurveJson(), ciCurveJson()];
    const table = (name: string) => {
      const col = DG.Column.fromStrings('curve', cells());
      col.semType = FitConstants.FIT_SEM_TYPE;
      const df = DG.DataFrame.fromColumns([col]);
      df.name = name;
      return df;
    };

    const source = table('layoutSource');
    const sourceView = grok.shell.addTableView(source);
    const fresh = table('layoutTarget');
    const freshView = grok.shell.addTableView(fresh);
    try {
      expect(getOrCreateParsedChartData(fresh.cell(0, 'curve')).series![0].showCurveConfidenceInterval, true,
        'the fresh table should start with the value its data declares');

      changeCurvesOptions(sourceView.grid.cell('curve', 0),
        {property: {name: 'showCurveConfidenceInterval'}, value: false} as unknown as DG.InputBase,
        'seriesOptions', 'Column');
      changeCurvesOptions(sourceView.grid.cell('curve', 0),
        {property: {name: 'statisticsMode'}, value: 'both'} as unknown as DG.InputBase,
        'chartOptions', 'Column');
      const layout = sourceView.saveLayout();

      freshView.loadLayout(layout);
      await delay(500);

      expect(getOrCreateParsedChartData(fresh.cell(0, 'curve')).series![0].showCurveConfidenceInterval, false,
        'the column-level option did not survive the layout');
      expect(getOrCreateParsedChartData(fresh.cell(0, 'curve')).chartOptions?.statisticsMode, 'both',
        'a chart option did not survive the layout');
    } finally {
      sourceView.close();
      freshView.close();
    }
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

  test('the panel offers the label names the cell carries', async () => {
    const chartData = JSON.parse(labelledCurveJson()) as IFitChartData;
    const choices = chartPropertiesFor(chartData).find((p) => p.name === 'showLabels')!.choices!;

    // plot-level first, then whatever the curves add - a name is offered once however many carry it
    expectArray(choices, ['Z prime', 'compound']);
  });

  test('a plot-level label is drawn once, a per-curve one per curve', async () => {
    const drawn = renderedTexts(labelledCurveJson(), 'panelLabels');

    // the plate statistic describes the cell, so it appears once however many curves there are
    expect(drawn.filter((t) => t.startsWith('Z prime')).length, 1, 'the plot-level label was not drawn once');
    expect(drawn.filter((t) => t.startsWith('compound')).length, 2, 'each curve should carry its own label');
    expect(drawn.includes('compound: GRK-1'), true);
    expect(drawn.includes('compound: GRK-2'), true);
  });

  test('only the labels named in showLabels are drawn', async () => {
    const chartData = JSON.parse(labelledCurveJson()) as IFitChartData;
    chartData.chartOptions!.showLabels = ['compound'];
    const drawn = renderedTexts(JSON.stringify(chartData), 'panelLabelsFiltered');

    expect(drawn.some((t) => t.startsWith('Z prime')), false, 'a label not named in showLabels was drawn');
    expect(drawn.filter((t) => t.startsWith('compound')).length, 2);
  });

  test('an error bar starts at the marker edge and hides when the marker covers it', async () => {
    // screen y grows downwards, so the top of the bar is the smaller number
    const marker4 = stdevWhisker(100, 70, 130, 2)!;
    expect(marker4 !== null, true, 'a deviation larger than the marker should be drawn');
    expectArray(marker4.top, [70, 98]);
    expectArray(marker4.bottom, [102, 130]);

    // a point can set its own size, and a large marker used to swallow the bar whole
    const bigMarker = stdevWhisker(100, 70, 130, 40);
    expect(bigMarker === null, true, 'a deviation smaller than the marker must not be drawn inside it');

    // the boundary: exactly at the marker edge is still nothing to show
    expect(stdevWhisker(100, 80, 120, 20) === null, true);
    expect(stdevWhisker(100, 79, 121, 20) !== null, true);
  });

  test('one list of statistics, and the mode decides where it is drawn', async () => {
    const data = JSON.parse(multiSeriesCurveJson([-7, -5])) as IFitChartData;
    data.chartOptions!.showStatistics = ['ic50'];
    data.chartOptions!.statisticsMode = 'both';
    data.chartOptions!.aggrType = 'med';
    const both = renderedTexts(JSON.stringify(data), 'panelAggrBoth');

    // named the way the extracted column is, so the plot and the column agree
    const line = both.find((t) => t.startsWith('med ic50'));
    expect(line !== undefined, true, 'the aggregated statistic was not drawn');
    expect(both.filter((t) => t.startsWith('med ic50')).length, 1, 'it describes the cell, so it is drawn once');
    expect(both.filter((t) => t.startsWith('IC50')).length, 2, 'each curve should still show its own');
    // aggregated in fit space, so two decades apart give the middle one, not the arithmetic mean
    expect(line!.includes('1.0') || line!.includes('1e-6') || line!.includes('9.9'), true, `unexpected value: ${line}`);

    data.chartOptions!.statisticsMode = 'aggregated';
    const aggregated = renderedTexts(JSON.stringify(data), 'panelAggrOnly');
    expect(aggregated.filter((t) => t.startsWith('med ic50')).length, 1);
    expect(aggregated.some((t) => t.startsWith('IC50')), false, 'the per-series values should give way to the summary');

    data.chartOptions!.statisticsMode = 'series';
    const series = renderedTexts(JSON.stringify(data), 'panelAggrSeries');
    expect(series.some((t) => t.startsWith('med ')), false, 'nothing should be summarised in series mode');
    expect(series.filter((t) => t.startsWith('IC50')).length, 2);
  });

  test('a single-series cell has nothing to aggregate', async () => {
    const data = JSON.parse(multiSeriesCurveJson([-6.5])) as IFitChartData;
    data.chartOptions!.showStatistics = ['ic50'];
    data.chartOptions!.statisticsMode = 'aggregated';
    const drawn = renderedTexts(JSON.stringify(data), 'panelAggrSingle');

    expect(drawn.some((t) => t.startsWith('med ')), false,
      'with one curve the aggregate is that curve, so it should not be drawn');
    // and the mode falls back to per series rather than drawing nothing at all
    expect(drawn.some((t) => t.startsWith('IC50')), true, 'the single curve should still show its own');
  });

  test('the panel offers droplines only for a fit that has asymptotes', async () => {
    const named = (json: string): string[] =>
      seriesPropertiesFor(JSON.parse(json) as IFitChartData, null).map((p) => p.name);

    expect(named(curveJson(-6.5)).includes('droplines'), true, 'a sigmoid levels off at both ends');

    const linear = JSON.parse(curveJson(-6.5)) as IFitChartData;
    linear.series![0].fitFunction = 'linear';
    expect(named(JSON.stringify(linear)).includes('droplines'), false,
      'a line has no asymptotes, so an ICxx would never be drawn');
    // the rest of the series options are untouched
    expect(named(JSON.stringify(linear)).includes('markerType'), true);
    // the render stamps this one on every paint, so editing it never took
    expect(named(curveJson(-6.5)).includes('columnName'), false, 'the column a series came from is not a setting');
  });

  test('a custom fit function is offered once, everywhere, and never an empty one', async () => {
    const custom = {name: 'TestCustomFit', function: '(params, x) => params[0] * x',
      getInitialParameters: '(x, y) => new Float32Array([1])', parameterNames: ['Slope']};
    const withCustom = JSON.parse(curveJson(-6.5)) as IFitChartData;
    withCustom.series![0].fitFunction = custom as any;
    try {
      // rendering the cell is what registers it, and the panel is opened after that
      getSeriesFitFunction(withCustom.series![0]);
      const choices = seriesPropertiesFor(withCustom, custom as any).find((p) => p.name === 'fitFunction')!.choices!;
      expect(choices.filter((c) => c === custom.name).length, 1, 'a registered function is named once');
      expect(choices.includes(''), false, 'a series is fitted with something, so none is not a choice');

      // registered once, it can be picked for any curve, not only the one that brought it
      const plain = seriesPropertiesFor(JSON.parse(curveJson(-6.5)) as IFitChartData, null)
        .find((p) => p.name === 'fitFunction')!.choices!;
      expect(plain.includes(custom.name), true, 'and offered wherever a fit function is picked');
      expect(plain.includes(''), false, 'here too');
    } finally {
      delete fitFunctions[custom.name];
    }
  });

  test('picking a custom fit function stores the notation, not a name nothing can resolve', async () => {
    const custom = {name: 'TestStoredFit', function: '(params, x) => params[0] * x',
      getInitialParameters: '(x, y) => new Float32Array([1])', parameterNames: ['Slope']};
    try {
      // rendering the column that carries it is what registers it, and its notation with it
      const carrier = JSON.parse(curveJson(-6.5)) as IFitChartData;
      carrier.series![0].fitFunction = custom as any;
      getSeriesFitFunction(carrier.series![0]);

      // a curve that has never seen that function, the way another column would be
      const cell = curveJson(-6.5);
      const col = DG.Column.fromStrings('curve', [cell, cell]);
      col.semType = FitConstants.FIT_SEM_TYPE;
      const df = DG.DataFrame.fromColumns([col]);
      df.name = 'customFitNotationStored';
      const gridCell = DG.Viewer.grid(df).cell('curve', 0);

      changeCurvesOptions(gridCell, {property: {name: 'fitFunction'}, value: custom.name} as unknown as DG.InputBase,
        'seriesOptions', 'Cell');
      await delay(300);
      const stored = JSON.parse(df.col('curve')!.get(0)!).series[0].fitFunction;
      expect(typeof stored, 'object', 'the name alone would not resolve in a session that never loaded it');
      expect(stored.name, custom.name);
      expect(!!stored.function && !!stored.getInitialParameters, true, 'the notation travels with the curve');

      changeCurvesOptions(gridCell, {property: {name: 'fitFunction'}, value: 'linear'} as unknown as DG.InputBase,
        'seriesOptions', 'Cell');
      await delay(300);
      expect(JSON.parse(df.col('curve')!.get(0)!).series[0].fitFunction, 'linear',
        'a built-in one is always there, so it stays a name');
    } finally {
      delete fitFunctions[custom.name];
      delete fitFunctionDescriptions[custom.name];
    }
  });

  test('statistics are written in the plot font, not the one the caller left', async () => {
    const data = JSON.parse(curveJson(-6.5)) as IFitChartData;
    data.chartOptions!.showStatistics = ['auc', 'rSquared'];
    const col = DG.Column.fromStrings('curve', [JSON.stringify(data)]);
    col.semType = FitConstants.FIT_SEM_TYPE;
    const df = DG.DataFrame.fromColumns([col]);
    df.name = 'statisticsFont';

    // the grid hands the renderer whatever font it last drew with, which is what used to change
    // the statistics when the filter opened and change them back when the panel was scrolled
    const fontsWhenCallerLeft = (callerFont: string): string[] => {
      const fonts: string[] = [];
      const g = ui.canvas(400, 300).getContext('2d')!;
      const fillText = g.fillText.bind(g);
      g.fillText = ((text: string, x: number, y: number) => {
        if (text.includes(':')) fonts.push(g.font);
        fillText(text, x, y);
      }) as any;
      g.font = callerFont;
      new FitChartCellRenderer().renderCurves(g, new DG.Rect(0, 0, 400, 300),
        getOrCreateParsedChartData(df.cell(0, 'curve')));
      return [...new Set(fonts)];
    };

    const seeded = fontsWhenCallerLeft('20px Comic Sans MS');
    expect(seeded.length, 1, 'every statistic should be written alike');
    expect(seeded[0].includes('Comic Sans'), false, 'and not in whatever the caller was using');
    expect(fontsWhenCallerLeft('11px Roboto, "Roboto Local"')[0], seeded[0],
      'the same font however the canvas reached the renderer');
  });

  test('several droplines are named, and one off the plot is not drawn', async () => {
    // the curve levels off at 5 and 100, so IC50 sits at the inflection and IC90 a decade later
    const data = JSON.parse(curveJson(-6.5)) as IFitChartData;
    data.series![0].droplines = ['IC50', 'IC90'];
    const both = renderedTexts(JSON.stringify(data), 'droplinesBoth');
    expect(both.filter((t) => t === 'IC50').length, 1, 'each drawn dropline should be named');
    expect(both.filter((t) => t === 'IC90').length, 1);

    // a single line has nothing to be confused with, so it stays unlabelled as before
    data.series![0].droplines = ['IC50'];
    expect(renderedTexts(JSON.stringify(data), 'droplinesOne').includes('IC50'), false,
      'one dropline needs no caption');

    // a grid cell has no room for captions, and crammed ones read worse than none
    data.series![0].droplines = ['IC50', 'IC90'];
    expect(renderedTexts(JSON.stringify(data), 'droplinesSmall', 200, 150).includes('IC50'), false,
      'captions should give way when the plot is small');

    // IC99 of a curve centred at 1e-5 lands at 1e-3, past the highest concentration tested
    const late = JSON.parse(curveJson(-5)) as IFitChartData;
    late.series![0].droplines = ['IC50', 'IC99'];
    const drawn = renderedTexts(JSON.stringify(late), 'droplinesOffPlot');
    expect(drawn.includes('IC50'), true);
    expect(drawn.includes('IC99'), false, 'a dropline off the plot should not be drawn at all');
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

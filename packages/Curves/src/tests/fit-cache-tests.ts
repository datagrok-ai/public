import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {FitConstants} from '@datagrok-libraries/statistics/src/fit/const';
import {IFitChartData} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {getSeriesFitFunction} from '@datagrok-libraries/statistics/src/fit/fit-data';
import {FitChartCellRenderer} from '../fit/fit-renderer';
import {chartDataId, getOrCreateCachedFitCurve, getOrCreateParsedChartData, viewerFits} from '../fit/fit-chart-data';

import {curveJson, multiSeriesCurveJson} from './curve-data';

/** A one-curve table, and the cell the caches are keyed on. */
function curveCell(name: string, logIC50: number = -6.5): DG.Cell {
  const col = DG.Column.fromStrings('curve', [curveJson(logIC50)]);
  col.semType = FitConstants.FIT_SEM_TYPE;
  const df = DG.DataFrame.fromColumns([col]);
  df.name = name;
  return df.cell(0, 'curve');
}

category('fit cache', () => {
  test('a curve fitted in log space is not handed to a linear plot', async () => {
    // the grid, the property panel and a trellis can each show one cell in a different space at once
    const cell = curveCell('fitCacheLogSpace');
    const series = getOrCreateParsedChartData(cell).series![0];
    const fitFunc = getSeriesFitFunction(series);

    const inLog = getOrCreateCachedFitCurve(series, 0, fitFunc, {logX: true, logY: false}, cell);
    const linear = getOrCreateCachedFitCurve(series, 0, fitFunc, {logX: false, logY: false}, cell);
    expect(Array.from(inLog.parameters).toString() !== Array.from(linear.parameters).toString(), true,
      `the two spaces should fit differently, both gave ${Array.from(inLog.parameters)}`);

    const again = getOrCreateCachedFitCurve(series, 0, fitFunc, {logX: true, logY: false}, cell);
    expect(again === inLog, true, 'asking twice in the same space should come from the cache');
  });

  test('a viewer fit is reused, and never confused with another cell', async () => {
    // both at row 0, which with an unsaved table is most of what the grid's key looks at
    const one = curveCell('fitCacheRowsOne', -7);
    const other = curveCell('fitCacheRowsOther', -5);
    const oneData = getOrCreateParsedChartData(one);
    const otherData = getOrCreateParsedChartData(other);
    expect(chartDataId(oneData) !== chartDataId(otherData), true, 'two cells should not share an identity');

    const logOptions = {logX: true, logY: false};
    const fitFunc = getSeriesFitFunction(oneData.series![0]);
    const first = getOrCreateCachedFitCurve(oneData.series![0], 0, fitFunc, logOptions, undefined, false,
      `${chartDataId(oneData)} || idx: 0`);
    const second = getOrCreateCachedFitCurve(otherData.series![0], 0, fitFunc, logOptions, undefined, false,
      `${chartDataId(otherData)} || idx: 0`);
    expect(first === second, false, 'two different curves must not answer for each other');
    expect(Array.from(first.parameters).toString() !== Array.from(second.parameters).toString(), true,
      'and they fit differently, so a collision would be visible');

    const repeat = getOrCreateCachedFitCurve(oneData.series![0], 0, fitFunc, logOptions, undefined, false,
      `${chartDataId(oneData)} || idx: 0`);
    expect(repeat === first, true, 'the same curve should be fitted once');
  });

  test('the viewer cache is kept apart from the grid one, and agrees with it', async () => {
    const cell = curveCell('fitCacheApart');
    const data = getOrCreateParsedChartData(cell);
    const series = data.series![0];
    const fitFunc = getSeriesFitFunction(series);
    const logOptions = {logX: true, logY: false};

    const forGrid = getOrCreateCachedFitCurve(series, 0, fitFunc, logOptions, cell);
    const forViewer = getOrCreateCachedFitCurve(series, 0, fitFunc, logOptions, cell, true,
      `${chartDataId(data)} || idx: 0`);
    expect(forGrid === forViewer, false, 'the two caches must not share storage');
    expect(Array.from(forGrid.parameters).toString(), Array.from(forViewer.parameters).toString(),
      'but keeping them apart must not change the answer');
  });

  test('a render that names its curves fits them once', async () => {
    const data = JSON.parse(multiSeriesCurveJson([-7, -5])) as IFitChartData;
    const g = ui.canvas(300, 200).getContext('2d')!;
    const bounds = new DG.Rect(0, 0, 300, 200);
    const identities = ['probeCell || idx: 0', 'probeCell || idx: 1'];
    new FitChartCellRenderer().renderCurves(g, bounds, data, undefined, identities);

    const keys = (viewerFits as unknown as {K: string[]}).K.filter(Boolean);
    expect(identities.every((id) => keys.some((k) => k.startsWith(id))), true,
      `the render should have cached a fit per identity, holds ${keys.length}`);

    let refitted = false;
    const options = data.chartOptions;
    viewerFits.getOrCreate(`${identities[0]} || logX: ${options?.logX} || logY: ${options?.logY}`, () => {
      refitted = true;
      return null as any;
    });
    expect(refitted, false, 'a second render of the same curve should not fit it again');
  });

  test('a pipe in a title survives, and one that breaks the value is still stripped', async () => {
    const titled = JSON.parse(curveJson(-6.5)) as IFitChartData;
    titled.chartOptions!.title = 'compound 10 | Collapse Assay | run 4';
    // the second value is valid only once the stray pipe is gone
    const col = DG.Column.fromStrings('curve', [JSON.stringify(titled), '{"series":|[{"points":[]}]}']);
    col.semType = FitConstants.FIT_SEM_TYPE;
    const df = DG.DataFrame.fromColumns([col]);
    df.name = 'pipeTitles';

    expect(getOrCreateParsedChartData(df.cell(0, 'curve')).chartOptions!.title,
      'compound 10 | Collapse Assay | run 4', 'a pipe inside a title is not part of the value');
    // the rescue still applies to a value the pipes actually break
    expect(getOrCreateParsedChartData(df.cell(1, 'curve')).series!.length, 1,
      'stripping should still rescue a value that does not parse as it stands');
  });

  test('editing a cell retires the fits made from it', async () => {
    const cell = curveCell('fitCacheEdit', -7);
    const before = getOrCreateParsedChartData(cell);
    const fitFunc = getSeriesFitFunction(before.series![0]);
    const logOptions = {logX: true, logY: false};
    const beforeFit = getOrCreateCachedFitCurve(before.series![0], 0, fitFunc, logOptions, undefined, false,
      `${chartDataId(before)} || idx: 0`);

    cell.column.set(0, curveJson(-5));
    const after = getOrCreateParsedChartData(cell);
    expect(chartDataId(after) !== chartDataId(before), true, 'a re-parsed cell should get a new identity');
    const afterFit = getOrCreateCachedFitCurve(after.series![0], 0, fitFunc, logOptions, undefined, false,
      `${chartDataId(after)} || idx: 0`);
    expect(Array.from(afterFit.parameters).toString() !== Array.from(beforeFit.parameters).toString(), true,
      'so the edited curve is fitted afresh');
  });
});

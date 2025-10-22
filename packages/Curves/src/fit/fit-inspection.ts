/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import {FitSeries, FitFunctionType} from '@datagrok-libraries/statistics/src/fit/new-fit-API';
import {inspectCurve} from './fit-renderer';
import * as grok from 'datagrok-api/grok';

export function randomizeTableId() {
  return `${Math.random()}-${Math.random()}-${Math.random()}-${Math.random()}-${Math.random()}-${Math.random()}-${Math.random()}-${Math.random()}`;
}

export function inspectSeriesByName(records: Record<string, FitSeries>, seriesName: string, fitFunctionName: FitFunctionType): void {
  if (!seriesName || !fitFunctionName)
    return;
  inspectSeries(records[seriesName], fitFunctionName);
}

export function inspectSeries(series: FitSeries, fitFunctionName: FitFunctionType) {
  const seriesName = series?.name ?? 'Series';
  if (!series || !fitFunctionName)
    return;
  const curveCol = DG.Column.string('Curve', 1);
  curveCol.set(0, JSON.stringify({
    chartOptions: {
      logX: true,
      title: seriesName,
    },
    series: [{...series, fit: undefined, fitFunction: fitFunctionName, clickToToggle: true, droplines: ['IC50'], name: seriesName}],
  }), false);
  const df = DG.DataFrame.fromColumns([curveCol]);
  df.name = seriesName;
  df.id = randomizeTableId();
  curveCol.semType = 'fit';
  curveCol.tags[DG.Tags.CellRenderer] = 'fit';
  const grid = DG.Viewer.grid(df);
  const gridCell = grid.cell('Curve', 0);
  grok.shell.windows.showContextPanel = true;
  grok.shell.o = gridCell;
  inspectCurve(gridCell, {width: 480, height: 370}, true);
}

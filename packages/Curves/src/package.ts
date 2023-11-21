/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';

import {FitGridCellHandler, calculateSeriesStats, getChartDataAggrStats} from './fit/fit-grid-cell-handler';
import {FitChartCellRenderer, TAG_FIT_CHART_FORMAT, TAG_FIT_CHART_FORMAT_3DX, getChartData, substituteZeroes} from './fit/fit-renderer';
import {MultiCurveViewer} from './fit/multi-curve-viewer';
import {curveDemo} from './fit/fit-demo';
import {convertXMLToIFitChartData} from './fit/fit-parser';
import {LogOptions} from '@datagrok-libraries/statistics/src/fit/fit-data';


export const _package = new DG.Package();
const SOURCE_COLUMN_TAG = '.sourceColumn';
const SERIES_NUMBER_TAG = '.seriesNumber';
const SERIES_AGGREGATION_TAG = '.seriesAggregation';
const STATISTICS_TAG = '.statistics';


//name: Fit
//tags: cellRenderer
//meta.cellType: fit
//meta.virtual: true
//output: grid_cell_renderer result
export function fitCellRenderer(): FitChartCellRenderer {
  return new FitChartCellRenderer();
}

//name: MultiCurveViewer
//description: A viewer that superimposes multiple in-cell curves on one chart
//tags: viewer
//output: viewer result
export function _FitViewer(): MultiCurveViewer {
  return new MultiCurveViewer();
}

//name: Curve fitting
//description: Curve fitting is the process of constructing a curve, or mathematical function, that has the best fit to a series of data points
//meta.demoPath: Curves | Curve fitting
//test: curveFitDemo() //wait: 2000
export async function curveFitDemo(): Promise<void> {
  await curveDemo();
}

//tags: init
export function _initCurves(): void {
  DG.ObjectHandler.register(new FitGridCellHandler());
}

//name: addStatisticsColumn
//input: dynamic chartColumn
//input: dynamic p
//input: string seriesName
//input: int seriesNumber
export function addStatisticsColumn(chartColumn: DG.GridColumn, p: DG.Property, seriesName: string, seriesNumber: number): void {
  const grid = chartColumn.grid;
  const column = DG.Column.float(`${chartColumn.name}_${seriesName}_${p.name}`, chartColumn.column?.length);
  column.tags[SOURCE_COLUMN_TAG] = chartColumn.name;
  column.tags[SERIES_NUMBER_TAG] = seriesNumber;
  column.tags[STATISTICS_TAG] = p.name;

  column
    .init((i) => {
      const gridCell = DG.GridCell.fromColumnRow(grid, chartColumn.name, grid.tableRowToGrid(i));
      if (gridCell.cell.value === '')
        return null;
      const chartData = gridCell.cell.column.getTag(TAG_FIT_CHART_FORMAT) === TAG_FIT_CHART_FORMAT_3DX ?
        convertXMLToIFitChartData(gridCell.cell.value) : getChartData(gridCell);
      if (chartData.series![seriesNumber] === undefined || chartData.series![seriesNumber].points.every((p) => p.outlier))
        return null;
      if (chartData.chartOptions?.allowXZeroes && chartData.chartOptions?.logX &&
        chartData.series?.some((series) => series.points.some((p) => p.x === 0)))
        substituteZeroes(chartData);
      const chartLogOptions: LogOptions = {logX: chartData.chartOptions?.logX, logY: chartData.chartOptions?.logY};
      const fitResult = calculateSeriesStats(chartData.series![seriesNumber], chartLogOptions);
      return p.get(fitResult);
    });
  grid.dataFrame.columns.insert(column, chartColumn.idx);
}

//name: addAggrStatisticsColumn
//input: dynamic chartColumn
//input: dynamic p
//input: string aggrType
export function addAggrStatisticsColumn(chartColumn: DG.GridColumn, p: DG.Property, aggrType: string): void {
  const grid = chartColumn.grid;
  const column = DG.Column.float(`${chartColumn.name}_${aggrType}_${p.name}`, chartColumn.column?.length);
  column.tags[SOURCE_COLUMN_TAG] = chartColumn.name;
  column.tags[SERIES_AGGREGATION_TAG] = aggrType;
  column.tags[STATISTICS_TAG] = p.name;

  column
    .init((i) => {
      const gridCell = DG.GridCell.fromColumnRow(grid, chartColumn.name, grid.tableRowToGrid(i));
      if (gridCell.cell.value === '')
        return null;
      const chartData = gridCell.cell.column.getTag(TAG_FIT_CHART_FORMAT) === TAG_FIT_CHART_FORMAT_3DX ?
        convertXMLToIFitChartData(gridCell.cell.value) : getChartData(gridCell);
      if (chartData.series?.every((series) => series.points.every((p) => p.outlier)))
        return null;
      if (chartData.chartOptions?.allowXZeroes && chartData.chartOptions?.logX &&
        chartData.series?.some((series) => series.points.some((p) => p.x === 0)))
        substituteZeroes(chartData);
      const fitResult = getChartDataAggrStats(chartData, aggrType);
      return p.get(fitResult);
    });
  grid.dataFrame.columns.insert(column, chartColumn.idx);
}

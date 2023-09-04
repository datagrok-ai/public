import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {
  fitSeries,
  getColumnChartOptions,
  getSeriesStatistics,
  getSeriesFitFunction,
} from '@datagrok-libraries/statistics/src/fit/fit-data';
import {statisticsProperties, fitSeriesProperties, fitChartDataProperties, FIT_CELL_TYPE, TAG_FIT, IFitChartData, IFitSeries, IFitChartOptions} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {TAG_FIT_CHART_FORMAT, TAG_FIT_CHART_FORMAT_3DX, getChartData, mergeProperties} from './fit-renderer';
import {MultiCurveViewer} from './multi-curve-viewer';
import {convertXMLToIFitChartData} from './fit-parser';


const SOURCE_COLUMN_TAG = '.sourceColumn';
const SERIES_NUMBER_TAG = '.seriesNumber';
const STATISTICS_TAG = '.statistics';
const CHART_OPTIONS = 'chartOptions';
const SERIES_OPTIONS = 'seriesOptions';


function addStatisticsColumn(chartColumn: DG.GridColumn, p: DG.Property, seriesNumber: number): void {
  const grid = chartColumn.grid;
  const column = DG.Column.float(p.name, chartColumn.column?.length);
  column.tags[SOURCE_COLUMN_TAG] = chartColumn.name;
  column.tags[SERIES_NUMBER_TAG] = seriesNumber;
  column.tags[STATISTICS_TAG] = p.name;

  column
    .init((i) => {
      const chartData = getChartData(
        DG.GridCell.fromColumnRow(grid, chartColumn.name, grid.tableRowToGrid(i)));
      const fitResult = getSeriesStatistics(chartData.series![0], getSeriesFitFunction(chartData.series![0]));
      return p.get(fitResult);
    });
  grid.dataFrame.columns.add(column);
}

function changeCurvesOptions(gridCell: DG.GridCell, columnChartOptions: IFitChartData,
  inputBase: DG.InputBase, options: string): void {
  gridCell.cell.column.tags[TAG_FIT] = JSON.stringify(columnChartOptions);
  for (let i = 0; i < gridCell.cell.column.length; i++) {
    const value = gridCell.cell.column.get(i);
    if (value === '') continue;
    const chartData: IFitChartData = gridCell.cell.column.getTag(TAG_FIT_CHART_FORMAT) === TAG_FIT_CHART_FORMAT_3DX ?
      convertXMLToIFitChartData(value) : JSON.parse(value ?? '{}') ?? {};

    if (options === CHART_OPTIONS) {
      if (chartData.chartOptions === undefined) continue;
      (chartData.chartOptions[inputBase.property.caption as keyof IFitChartOptions] as any) = inputBase.value;
    }
    else if (options === SERIES_OPTIONS) {
      if (chartData.series === undefined) continue;
      for (let j = 0; j < chartData.series.length; j++)
        (chartData.series[j][inputBase.property.caption as keyof IFitSeries] as any) = inputBase.value;
    }
    gridCell.cell.column.set(i, JSON.stringify(chartData));
  }
  if (gridCell.cell.column.getTag(TAG_FIT_CHART_FORMAT) === TAG_FIT_CHART_FORMAT_3DX)
    gridCell.cell.column.setTag(TAG_FIT_CHART_FORMAT, '');
  gridCell.grid.invalidate();
}


export class FitGridCellHandler extends DG.ObjectHandler {
  get type(): string {
    return 'GridCell';
  }

  isApplicable(x: any): boolean {
    return x instanceof DG.GridCell && x.cellType === FIT_CELL_TYPE;
  }

  renderProperties(gridCell: DG.GridCell, context: any = null): HTMLElement {
    const acc = ui.accordion();
    const chartData = getChartData(gridCell);
    const columnChartOptions = getColumnChartOptions(gridCell.cell.column);
    const seriesOptionsRefresh = {onValueChanged: (inputBase: DG.InputBase) => 
      changeCurvesOptions(gridCell, columnChartOptions, inputBase, SERIES_OPTIONS)};
    const chartOptionsRefresh = {onValueChanged: (inputBase: DG.InputBase) =>
      changeCurvesOptions(gridCell, columnChartOptions, inputBase, CHART_OPTIONS)};

    mergeProperties(fitSeriesProperties, columnChartOptions.seriesOptions, chartData.seriesOptions ? chartData.seriesOptions :
      chartData.series ? chartData.series[0] ?? {} : {});
    mergeProperties(fitChartDataProperties, columnChartOptions.chartOptions,
      chartData.chartOptions ? chartData.chartOptions : {});

    acc.addPane('Options', () => ui.divV([
      ui.h3('Series options'),
      ui.input.form(chartData.seriesOptions ? chartData.seriesOptions : chartData.series![0], fitSeriesProperties, seriesOptionsRefresh),
      ui.h3('Chart options'),
      ui.input.form(chartData.chartOptions, fitChartDataProperties, chartOptionsRefresh),
    ]));

    acc.addPane('Fit', () => {
      const host = ui.divV([]);

      for (let i = 0; i < chartData.series!.length; i++) {
        const series = chartData.series![i];
        const fitFunction = getSeriesFitFunction(chartData.series![i]);
        if (!series.parameters)
          series.parameters = fitSeries(series, fitFunction).parameters;
        const seriesStatistics = getSeriesStatistics(series, fitFunction);
        const color = series.fitLineColor ? DG.Color.fromHtml(series.fitLineColor) ?
          series.fitLineColor : DG.Color.toHtml(DG.Color.getCategoricalColor(i)) : DG.Color.toHtml(DG.Color.getCategoricalColor(i));
        host.appendChild(ui.panel([
          ui.h1(series.name ?? 'series ' + i, {style: {color: color}}),
          ui.input.form(seriesStatistics, statisticsProperties, {
            onCreated: (input) => input.root.appendChild(
              ui.iconFA('plus',
                () => addStatisticsColumn(gridCell.gridColumn, input.property, i),
                `Calculate ${input.property.name} for the whole column`))
          })
        ]));
      }

      return host;
    });

    acc.addPane('Chart', () => MultiCurveViewer.fromChartData(chartData).root);

    return acc.root;
  }
}

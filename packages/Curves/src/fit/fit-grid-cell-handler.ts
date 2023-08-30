import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {
  fitSeries,
  getColumnChartOptions,
  getSeriesStatistics,
  getSeriesFitFunction,
} from '@datagrok-libraries/statistics/src/fit/fit-data';
import {statisticsProperties, fitSeriesProperties, fitChartDataProperties, FIT_CELL_TYPE, TAG_FIT} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {getChartData} from './fit-renderer';
import {MultiCurveViewer} from './multi-curve-viewer';


const SOURCE_COLUMN_TAG = '.sourceColumn';
const SERIES_NUMBER_TAG = '.seriesNumber';
const STATISTICS_TAG = '.statistics';


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
    const refresh = {onValueChanged: (_: any) => {
      gridCell.cell.column.tags[TAG_FIT] = JSON.stringify(columnChartOptions);
      gridCell.grid.invalidate();
    }};

    acc.addPane('Options', () => ui.divV([
      ui.input.form(columnChartOptions.seriesOptions, fitSeriesProperties, refresh),
      ui.input.form(columnChartOptions.chartOptions, fitChartDataProperties, refresh),
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
